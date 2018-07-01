/* BOREXINO Reconstruction program
 *
 * Author: Marcin Misiaszek <marcin.misiaszek@lngs.infn.it>
 * Maintainer: Marcin Misiaszek <marcin.misiaszek@lngs.infn.it>
 *
 *
 * Implementation of bx_energy_reco_lik
 *
 */

#include "bx_energy_reco_lik.hh"
#include "bx_echidna_event.hh"
#include "db_channel.hh"
#include "db_calib.hh"
#include "bx_dbi.hh"
#include "constants.hh"
#include "TMinuit.h"
#include "barn_interface.hh"

#include <math.h> 
#define pi 3.14159

namespace { 
  // MINUIT ROOT workaround
  bx_energy_reco_lik *current_module;
  void fcn (int32_t &npar, double *gin, double &f, double *x, int32_t iflag) 
  { f = current_module->lik_fcn (x); }
};

// ctor
bx_energy_reco_lik::bx_energy_reco_lik (): bx_base_module("bx_energy_reco_lik", bx_base_module::main_loop) {
  require_event_stage (bx_detector::laben, bx_base_event::clustered);
  current_module = this;
}

// module interface
void bx_energy_reco_lik::begin () {

  method = 1;
	
  // Initialization of Minuit
  p_minuit = new TMinuit(1);
  p_minuit->SetFCN(fcn);
  p_minuit->SetPrintLevel(-1);

  // Local cashing of pmts' positions and normal vectors
  f4_pmt_positions = new std::vector<vector3>(constants::laben::channels);
  f4_pmt_normal_vectors = new std::vector<vector3>(constants::laben::channels);
  b_pmt_cone = new std::vector<bool>(constants::laben::channels);
  
  float radius = bx_dbi::get()->get_profile().light_collector_radius();
//  f4_lg_length = bx_dbi::get()->get_profile().light_guide_lenght();
  f4_lg_entry_radius = bx_dbi::get()->get_profile().light_guide_entry_aperture();
  f4_cathode_radius = bx_dbi::get()->get_profile().pmt_chathode_radius();

  for (int32_t i=0; i < constants::laben::channels; ++i)  {
    const db_channel_laben& channel_info = dynamic_cast<const db_channel_laben&>(bx_dbi::get()->get_channel(i+1));
    ((*f4_pmt_positions)[i])[0] = channel_info.pmt_x();
    ((*f4_pmt_positions)[i])[1] = channel_info.pmt_y();
    ((*f4_pmt_positions)[i])[2] = channel_info.pmt_z();
    ((*f4_pmt_normal_vectors)[i])[0] = channel_info.pmt_x()/radius;
    ((*f4_pmt_normal_vectors)[i])[1] = channel_info.pmt_y()/radius;
    ((*f4_pmt_normal_vectors)[i])[2] = channel_info.pmt_z()/radius;
    (*b_pmt_cone)[i] = channel_info.pmt_has_cone(); 
  }
			    
}

bx_echidna_event* bx_energy_reco_lik::doit (bx_echidna_event *ev) {

  double vstart[1] = { 50. };
  double step[1] = { 1. };

  double f8_value, f8_error;

  // Loop on every cluster
  for (int32_t i = 0; i < ev->get_laben().get_nclusters(); i++)  {
  
    // Get position from rec_cluster
    //bx_energy_lik& e_lik = ev->get_laben().get_cluster(i).get_energy_lik();
    const bx_position& rec_p = 	ev->get_laben ().get_rec_cluster (i).get_position();

    // Get starting point and initialize minuit parameters
    p_minuit->DefineParameter (0, "E", vstart[0], step[0], 1.0, 3000.0);

    u1_collected_charge = new std::vector<uint8_t>(constants::laben::channels, 0);
    const bx_laben_cluster& cluster = ev->get_laben().get_cluster(i);

    // Loop on every hit
    for (int32_t i = 0; i < cluster.get_clustered_nhits(); ++i)  {

	//  db_channel for each hit
	const bx_laben_clustered_hit& hit = cluster.get_clustered_hit(i);
	const bx_laben_decoded_hit& d_hit = hit.get_decoded_hit();
	const db_channel* ch_info = d_hit.get_db_channel();
	if ( !ch_info->is_ordinary() ) continue;
    
	// collected charge, used for likelihood function calculation 
	(*u1_collected_charge)[ch_info->get_lg()-1] += 1;
	}
	
    vector3 position(rec_p.get_x(),rec_p.get_y(),rec_p.get_z());
    vector3 distance, distance_normal_vector;
    double f8_omega=0.;
    double f8_attenuation_factor=0.;
    v_omega_attenuation = new std::vector<double>(constants::laben::channels, 888.0);

//__________________________________________________________


    // Loop on every channel to reconstruct event's energy
    for (int32_t k=0; k < constants::laben::channels; ++k)  {
	const db_channel& ch_info = bx_dbi::get()->get_channel(k+1);
	if (!ch_info.is_ordinary()) continue;

if( method == 0)
  {
  distance = (*f4_pmt_positions)[k] - position;
  double path = distance.Mag();	  
	  
  if ( (*b_pmt_cone)[k] )  {
        distance_normal_vector = distance*(1/path);
        double costheta = (*f4_pmt_normal_vectors)[k] * distance_normal_vector;

        // Igor Machulin - Geant4 fit for the cone
        if (costheta < 0.92)  { path += 0.2898; }
        else                  { path += ( 3.134 - costheta*5.621 + costheta*costheta*2.748 ); }

        f8_omega = (pi*f4_lg_entry_radius*f4_lg_entry_radius*pow(costheta,1.7))/(path*path);
        f8_attenuation_factor = 1.;
//      f8_attenuation_factor = exp(path/attenuation_length);  // Scint.
	}
	// no cone
	else {
        distance_normal_vector = distance*(1/path);
        double costheta = (*f4_pmt_normal_vectors)[k] * distance_normal_vector;

        f8_omega = (pi*f4_cathode_radius*f4_cathode_radius*costheta)/(path*path);
        f8_attenuation_factor = 1.;
//      f8_attenuation_factor = exp(path/attenuation_length);  // Scint.
	}
  }

double path=0.;

if( method == 1 )
  {
  double attenuation_length = 8.0;
    
    
	double A;				//A - length of square sides 				
  	if ( (*b_pmt_cone)[k] ) {
	    A = f4_lg_entry_radius * sqrt(pi);
	    distance = 
	      (*f4_pmt_positions)[k] * (1. - 0.3 * ( 1. / ((*f4_pmt_positions)[k]).Mag()) ) - position;
	  }
	  else {
	    A = f4_cathode_radius * sqrt(pi);
	    distance = (*f4_pmt_positions)[k] - position;
	  }
	
	  
	path = distance.Mag();  
	distance_normal_vector = distance*(1/path);
	double costheta = (*f4_pmt_normal_vectors)[k] * distance_normal_vector;
	double sintheta = sqrt( 1. - ( costheta * costheta ) );
	double sectheta = 1 / costheta;

	f8_omega = 
	  2*(-atan((A*sectheta*(-A/2. + path*sintheta))/(path*sqrt(2*pow(A,2) + 
	  4*pow(path,2)*pow(costheta,2) - 4*A*path*sintheta + 4*pow(path,2)*pow(sintheta,2)))) + 
	  atan((A*sectheta*(A + 2*path*sintheta))/(2.*path*sqrt(2*pow(A,2) + 
	  4*pow(path,2)*pow(costheta,2) + 4*A*path*sintheta + 4*pow(path,2)*pow(sintheta,2)))));
	  f8_attenuation_factor = exp( -path/attenuation_length);  // Scint.
	
	//f8_attenuation_factor = 1.;
  }	

  if (path < 4.5 ) continue;

  (*v_omega_attenuation)[k] = f8_omega * f8_attenuation_factor;
    
  }  
    
    
    // Minimization with migrad
    p_minuit->Command ("SET NOW");
    p_minuit->Command ("MIGRAD");
    
    p_minuit->GetParameter (0, f8_value, f8_error);
 
#if 0   
    double cf=0.; //callibration factor
    
      switch (method) {
        case 0: cf = 0.00875274;  break;
	case 1: cf = 0.00481538;  break;
	}
      
//    e_lik.f4_e = f8_value * cf;
#endif

    
    //    std::cout << p.f4_e << std::endl;

    
    delete u1_collected_charge;
    delete v_omega_attenuation;
  }
  
    return ev;
}




// the value to be minimized
double bx_energy_reco_lik::lik_fcn (double *x) {

    f8_result = 0.;

    // Loop on every channel to reconstruct event's energy
    for (int32_t k=0; k < constants::laben::channels; ++k)  {
	const db_channel& ch_info = bx_dbi::get()->get_channel(k+1);
	if (!ch_info.is_ordinary()) continue;	
        if ( (*v_omega_attenuation)[k] == 888.0 ) continue;

	expected_charge = (*v_omega_attenuation)[k] * x[0];
    
	if ((*u1_collected_charge)[k]) f8_prob = ( 1.0 - exp( - expected_charge ));
	    else  
		f8_prob = exp( - expected_charge );


//	double real_q = (*u1_collected_charge)[k];
//	f8_prob = TMath::Poisson(real_q, expected_charge);
		
	
	f8_result -= log( f8_prob );
    }

  return f8_result;
}



void bx_energy_reco_lik::end () {
}


bx_energy_reco_lik::~bx_energy_reco_lik () {    
}


