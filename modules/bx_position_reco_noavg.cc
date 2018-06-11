/*
 * Copy of position_lngs algorithm with noavg correction 
 *
 * $Id: bx_position_reco_noavg.cc,v 1.3 2015/08/17 17:24:22 misiaszek Exp $
 */

#include "bx_position_reco_noavg.hh"
#include "bx_echidna_event.hh"
#include "db_channel.hh"
#include "db_calib.hh"
#include "bx_dbi.hh"
#include "constants.hh"
#include "TMinuit.h"
#include "barn_interface.hh"

#include "bx_pdf_noavg.hh"

#include <math.h> 
#define pi 3.141592653589793238
#define speed_of_light_m_ns 0.3

namespace { 
  // MINUIT ROOT workaround
  bx_position_reco_noavg *current_module_n;
  void fcn_n (int32_t &npar, double *gin, double &f, double *x, int32_t iflag) { f = current_module_n->my_fcn (x); }
};

// ctor
bx_position_reco_noavg::bx_position_reco_noavg() : bx_base_module("bx_position_reco_noavg", bx_base_module::main_loop) {
  require_event_stage (bx_detector::laben, bx_base_event::clustered);
  current_module_n = this; 	// see minuit root workaround
}

void bx_position_reco_noavg::begin () {
//
// Cut on arrival time of photons
//
   f8_tmax = get_parameter ("tmax").get_float ();
//
   f4_minuit_ok=0;
   f4_minuit_blank=0;
   f4_minuit_unreadable=0;
   f4_minuit_unknown=0;
   f4_minuit_abnormal=0;

   f4_iconv_0=0;
   f4_iconv_1=0;
   f4_iconv_2=0;
   f4_iconv_3=0;
//
// To initialize Minuit
  p_minuit = new TMinuit(4);
  p_minuit->SetPrintLevel(-1);
  p_minuit->SetFCN(fcn_n);

  // Local cashing of pmts' positions and normal vectors
  f4_pmt_positions = new std::vector<TVector3>(constants::laben::channels);
  f4_pmt_normal_vectors = new std::vector<TVector3>(constants::laben::channels);
  b_pmt_cone = new std::vector<bool>(constants::laben::channels);
  
  if (get_parameter ("refidx").get_float () > 0) f4_ref_index = get_parameter ("refidx").get_float ();
  else f4_ref_index = bx_dbi::get()->get_calib().get_refraction_index_data();
//  
  float radius = bx_dbi::get()->get_profile().light_collector_radius();
//
  f4_t_0 = 35. + (f4_ref_index*radius/speed_of_light_m_ns);
//  get_message(bx_message::log) << "T0 starting point" << f4_t_0 << dispatch;
  
  for ( int32_t i=0; i<constants::laben::channels; i++ ) {
    const db_channel_laben& channel_info = dynamic_cast<const db_channel_laben&>(bx_dbi::get()->get_channel(i+1));
    ((*f4_pmt_positions)[i])(0) = channel_info.pmt_x();
    ((*f4_pmt_positions)[i])(1) = channel_info.pmt_y();
    ((*f4_pmt_positions)[i])(2) = channel_info.pmt_z();
    float raggio = sqrt (channel_info.pmt_x()*channel_info.pmt_x()+channel_info.pmt_y()*channel_info.pmt_y() + channel_info.pmt_z()*channel_info.pmt_z());    
    ((*f4_pmt_normal_vectors)[i])(0) = channel_info.pmt_x()/raggio;
    ((*f4_pmt_normal_vectors)[i])(1) = channel_info.pmt_y()/raggio;
    ((*f4_pmt_normal_vectors)[i])(2) = channel_info.pmt_z()/raggio;
    (*b_pmt_cone)[i] = channel_info.pmt_has_cone(); 
  }
  get_message(bx_message::debug) << "Local cashing  of pmts' positions done" << dispatch;

  // Create the light guide object and set its base_module pointer
  // coment out the following two lines for speed up (yusuke koshio)
  //light_guide::get();
  //get_message(bx_message::debug) << "Light guide object created" << dispatch;
  
  // Create a new bx_interpolator for the pdf function
  const std::vector<float> pdf_time_v = bx_dbi::get()->get_calib().get_pdf_time_vector();
  const std::vector<float> pdf_data_v = bx_dbi::get()->get_calib().get_pdf_data_vector();
  double pdf_time[pdf_time_v.size()];
  double pdf_data[pdf_time_v.size()];
  
  double pdf_data2[pdf_time_v.size()],pdf_data3[pdf_time_v.size()];
  double pdf_data4[pdf_time_v.size()],pdf_data5[pdf_time_v.size()];
  double pdf_data6[pdf_time_v.size()],pdf_data7[pdf_time_v.size()];
  double pdf_data8[pdf_time_v.size()],pdf_data9[pdf_time_v.size()];
  double pdf_data10[pdf_time_v.size()];

  for ( int32_t i=0; i<(int32_t)pdf_time_v.size(); i++ ) {
    pdf_time[i] = pdf_time_v[i];
    pdf_data[i] = pdf_data_v[i];

    pdf_data2[i] = pdf_n_02[i];
    pdf_data3[i] = pdf_n_03[i];
    pdf_data4[i] = pdf_n_04[i];
    pdf_data5[i] = pdf_n_05[i];
    pdf_data6[i] = pdf_n_06[i];
    pdf_data7[i] = pdf_n_07[i];
    pdf_data8[i] = pdf_n_08[i];
    pdf_data9[i] = pdf_n_09[i];
    pdf_data10[i] = pdf_n_10[i];
  }
  pdf1 = new interpolator(pdf_time_v.size(),pdf_time,pdf_data); 
  pdf2 = new interpolator(pdf_time_v.size(),pdf_time,pdf_data2); 
  pdf3 = new interpolator(pdf_time_v.size(),pdf_time,pdf_data3); 
  pdf4 = new interpolator(pdf_time_v.size(),pdf_time,pdf_data4); 
  pdf5 = new interpolator(pdf_time_v.size(),pdf_time,pdf_data5); 
  pdf6 = new interpolator(pdf_time_v.size(),pdf_time,pdf_data6); 
  pdf7 = new interpolator(pdf_time_v.size(),pdf_time,pdf_data7); 
  pdf8 = new interpolator(pdf_time_v.size(),pdf_time,pdf_data8); 
  pdf9 = new interpolator(pdf_time_v.size(),pdf_time,pdf_data9); 
  pdf10 = new interpolator(pdf_time_v.size(),pdf_time,pdf_data10); 

  get_message(bx_message::debug) << "Pdf interpolator object created" << dispatch;

  get_message(bx_message::debug) << "begin" << dispatch;
}

bx_echidna_event* bx_position_reco_noavg::doit (bx_echidna_event *ev) {

  // Loop on every cluster
  trigger_type = ev->get_trigger ().get_trigger_type ();
  Int_t size_clusters = ev->get_laben ().get_nclusters ();
  Int_t size_mclusters = ev->get_laben ().get_nclusters_muons (); 
  for (int32_t i = 0; i < size_clusters + size_mclusters; i++) {
    bx_laben_cluster& cluster = (i < size_clusters ) ? ev->get_laben().get_cluster(i) : ev->get_laben().get_cluster_muon(i-size_clusters);
    //const bx_baricenter& b = cluster.get_baricenter();
    const bx_position_lngs& b_pos=cluster.get_position_lngs ();
    bx_position_noavg& pos = cluster.get_position_noavg ();

    float f4_t_b;
    float xb=b_pos.get_x();
    float yb=b_pos.get_y();
    float zb=b_pos.get_z();

    float f4_t_r=f4_ref_index/speed_of_light_m_ns*sqrt((xb*xb)+(yb*yb)+(zb*zb));
    f4_t_b=f4_t_0-f4_t_r;

    // Get starting point and initialize minuit parameters
    p_minuit->DefineParameter (0, "x position", xb, 1.e-2, -10., 10.);
    p_minuit->DefineParameter (1, "y position", yb, 1.e-2, -10., 10.);
    p_minuit->DefineParameter (2, "z position", zb, 1.e-2, -10., 10.);
    p_minuit->DefineParameter (3, "time event", f4_t_b, 1., -10., 100.);
    
    // To be used in my_fcn (see minuit root workaround)
//    p_fit_ev = ev;
//    i4_fit_cluster = i;
    p_fit_clu = &cluster;
    i4_n_iterations = 0;
    p_minuit->Command("SET NOW");

//  Minimization with MINIMIZE (same as MIGRAD, but calls SIMPLEX when MIGRAD fails and then calls MIGRAD again)
//  ierr=0 ok; ierr=1 command blank; ierr=2 command line unreadable; ierr=3 unknown command; ierr=4 abnormal termination;
//    
    double position[4], er[4], position_min[3];
    int32_t ierr = 999;
    for (int32_t k=0; k<10;k++){
      ierr= p_minuit->Command("MIN");
      if(ierr==0) {
	for ( int32_t j=0; j<4; j++ ) p_minuit->GetParameter(j, position[j], er[j]);
	position_min[0]=position[0];
	position_min[1]=position[1];
	position_min[2]=position[2];

	p_minuit->Command("IMP");
	for ( int32_t j=0; j<4; j++ ) p_minuit->GetParameter(j, position[j], er[j]);
	double dist_min=sqrt((position[0]-position_min[0])*(position[0]-position_min[0])+(position[1]-position_min[1])*(position[1]-position_min[1])+(position[2]-position_min[2])*(position[2]-position_min[2]));
	if(dist_min < 0.1) break;
      }
    }
    p_minuit->Command ("HESSE");
//

    int32_t nvpar,nparx,iconv;
//  nvpar=number of variable parameters
//  nparx=highest parameter number defined by user
//  iconv indicator of goodness of covariance 
//  iconv=0 not calculated at all
//  iconv=1 approximation only, not accurate
//  iconv=2 full matrix, but forced positive-definite
//  iconv=3 full accurate cov matrix
//
    double amin,edm,errdef;
//  amin= best function value found so far
//  edm= estimated vertical distance remaining to minimum
//  errdef= value of UP defining parameter uncertainties: by default=1
//
    p_minuit->mnstat(amin,edm,errdef,nvpar,nparx,iconv);
//    
    for ( int32_t j=0; j<4; j++ ) {
      p_minuit->GetParameter(j, position[j], er[j]);
//      p_minuit->mnerrs(j,eplus[j],eminus[j],eparab,globcc);
//      error[j] = ( fabs(eplus[j]) + fabs(eminus[j]) )/2.;
    }

    switch ( ierr ) {
      case 0:
//	get_message(bx_message::debug) << "MINUIT: Command executed normally" << dispatch;
        b_iconv = true;
	f4_minuit_ok++;
        break;
      case 1:
//	get_message(bx_message::debug) << "MINUIT: Command is blank, ignored" << dispatch;
        b_iconv = false;
	f4_minuit_blank++;
	break;
      case 2:
//	get_message(bx_message::debug) << "MINUIT: Command line unreadable, ignored" << dispatch;
        b_iconv = false;
        f4_minuit_unreadable++;
	break;
      case 3:
//	get_message(bx_message::debug) << "MINUIT: Unknown command, ignored" << dispatch;
        b_iconv = false;
	f4_minuit_unknown++;
	break;
      case 4: 
//	get_message(bx_message::debug) << "MINUIT: Abnormal termination" << dispatch;
        b_iconv = false;
	f4_minuit_abnormal++;
	break;
    }
    switch ( iconv ) {
      case 0:
//	get_message(bx_message::debug) << "HESSE: no cov matrix calculated " << dispatch;
	i4_hesse = -1;
	f4_iconv_0++;
        break;
      case 1:
//	get_message(bx_message::debug) << "HESSE: cov matrix approximated " << dispatch;
	i4_hesse = 0;
	f4_iconv_1++;
	break;
      case 2:
//	get_message(bx_message::debug) << "HESSE: full cov matrix, but forced pos-def " << dispatch;
	i4_hesse = 0;
        f4_iconv_2++;
	break;
      case 3:
//	get_message(bx_message::debug) << "HESSE: cov matrix ok " << dispatch;
	i4_hesse = 1;
	f4_iconv_3++;
	break;
    }
    double likelihood = my_fcn (position);
    pos.f4_x = position[0];
    pos.f4_y = position[1];
    pos.f4_z = position[2];
    pos.f4_t = position[3];
    pos.f4_dx = er[0];
    pos.f4_dy = er[1];
    pos.f4_dz = er[2];
    pos.f4_dt = er[3];
    pos.f4_user = likelihood;
    pos.i4_matrix = i4_hesse;
    pos.b_converged = b_iconv;

//
  }
//  ev->get_laben ().mark_stage (bx_base_event::reconstructed_noavg);
  
  return ev;
}

void bx_position_reco_noavg::end () {

  // Delete vectors and pointers
  delete f4_pmt_positions;
  delete f4_pmt_normal_vectors;
  delete b_pmt_cone;
  delete pdf1;
  delete pdf2;
  delete pdf3;
  delete pdf4;
  delete pdf5;
  delete pdf6;
  delete pdf7;
  delete pdf8;
  delete pdf9;
  delete pdf10;
  delete p_minuit;  
  get_message(bx_message::log) << "MINUIT: Command executed normally " << f4_minuit_ok <<" times"<< dispatch;
  get_message(bx_message::log) << "MINUIT: Command blank " << f4_minuit_blank <<"times"<< dispatch;
  get_message(bx_message::log) << "MINUIT: Command unreadable " << f4_minuit_unreadable <<" times"<< dispatch;
  get_message(bx_message::log) << "MINUIT: Command unknown " << f4_minuit_unknown <<" times"<< dispatch;
  get_message(bx_message::log) << "MINUIT: Abnormal termination " << f4_minuit_abnormal <<" times"<< dispatch;

  get_message(bx_message::log) << "HESSE: No cov matrix " << f4_iconv_0 <<" times"<< dispatch;
  get_message(bx_message::log) << "HESSE: cov matrix approximated " << f4_iconv_1 <<" times"<< dispatch;
  get_message(bx_message::log) << "HESSE: ful cov matrix, but forced pos-def  " << f4_iconv_2 <<" times"<< dispatch;
  get_message(bx_message::log) << "HESSE: cov matrix ok " << f4_iconv_3 <<" times"<< dispatch;

  get_message(bx_message::debug) << "end" << dispatch;
}

// Driver function for the position minimization: computes the value to be minimized.
double bx_position_reco_noavg::my_fcn (double *x) {

//  const bx_laben_cluster& cluster = p_fit_ev->get_laben().get_cluster(i4_fit_cluster);
//  int32_t n_evento = p_fit_ev->get_event_number();
  TVector3 ev_position(x);
  TVector3 distance,distance_normal_vector;
  double ev_time = -x[3];
  float lg_lenght = bx_dbi::get()->get_profile().light_guide_lenght();
//  float radius = bx_dbi::get()->get_profile().stainless_steel_sphere_radius();
//  double ev_time = (ev_position.Mag() - 2*radius) * (f4_ref_index/speed_of_light_m_ns);
  double f8_result=0.;
  i4_n_iterations++;

      
  // Loop on every hit
  for (int32_t i = 0; i < p_fit_clu->get_clustered_nhits (); i++) {
    // Get clustered hit reference and db_channel pointer
    const bx_laben_clustered_hit& hit = p_fit_clu->get_clustered_hit(i);
    const bx_laben_decoded_hit& dhit = hit.get_decoded_hit();
    const db_channel* ch_info = dhit.get_db_channel();
    if (!ch_info->is_ordinary ()) continue;
    if (hit.get_time() > f8_tmax) continue;
    
    int32_t index = ch_info->get_lg()-1;
    distance = (*f4_pmt_positions)[index] - ev_position;
    double path = distance.Mag(); 
      
    // Add the light guide contribution
    // change for speed up (yusuke koshio)
    if ( (*b_pmt_cone)[index] ) {
      //distance_normal_vector = distance*(1/path);
      //double costheta = (*f4_pmt_normal_vectors)[index] * distance_normal_vector;
      //if(costheta > 1.) {
//       get_message(bx_message::warn) <<" costheta >1 !!! in event number =    "<< n_evento << "    costheta =   "<< costheta << dispatch;
      // costheta=1.;
      //}
      //if(costheta < -1.) {
//       get_message(bx_message::warn) <<" costheta < -1 !!! in event number =    "<< n_evento << "    costheta =   "<< costheta << dispatch;
      // costheta =-1.;
      //}
      //double theta = 180.*acos(costheta)/pi;
      //path += light_guide::get()->get_path(theta);
      path += 0.27; // add 0.04 for reflection effect by cone
    }
    else {
      path += lg_lenght;
    }
    

    // Time of flight correction
    double time = hit.get_time();
    double delay = time - (ev_time + (f4_ref_index*path)/speed_of_light_m_ns);
    float cc = dhit.get_charge_pe();
    if(cc<0.1) 	cc=1;//noavg correction for charge - only difference from lngs code
    double pdf_value;
    if( cc < 1. ) {
      pdf_value = pdf1->get_value(delay);
    }
    else if( cc < 2. ) {
      pdf_value = (2.-cc) * pdf1->get_value(delay) + (cc-1.) * pdf2->get_value(delay) ;
    }
    else if( cc < 3. ) {
      pdf_value = (3.-cc) * pdf2->get_value(delay) + (cc-2.) * pdf3->get_value(delay);
    }
    else if( cc < 4. ) {
      pdf_value = (4.-cc) * pdf3->get_value(delay) + (cc-3.) * pdf4->get_value(delay);
    }
    else if( cc < 5. ) {
      pdf_value = (5.-cc) * pdf4->get_value(delay) + (cc-4.) * pdf5->get_value(delay);
    }
    else if( cc < 6. ) {
      pdf_value = (6.-cc) * pdf5->get_value(delay) + (cc-5.) * pdf6->get_value(delay);
    }
    else if( cc < 7. ) {
      pdf_value = (7.-cc) * pdf6->get_value(delay) + (cc-6.) * pdf7->get_value(delay);
    }
    else if( cc < 8. ) {
      pdf_value = (8.-cc) * pdf7->get_value(delay) + (cc-7.) * pdf8->get_value(delay);
    }
    else if( cc < 9. ) {
      pdf_value = (9.-cc) * pdf8->get_value(delay) + (cc-8.) * pdf9->get_value(delay);
    }
    else if( cc < 10. ) {
      pdf_value = (10.-cc) * pdf9->get_value(delay) + (cc-9.) * pdf10->get_value(delay);
    }
    else {
      pdf_value = pdf10->get_value(delay);
    }

    if(pdf_value > 0.){ 
     f8_result -= log(pdf_value);
    }
    if(pdf_value <= 0.){
     get_message(bx_message::warn) << "bx_position_reco_noavg: in MY_FCN pdf_value < 0!!! " << dispatch;
    }
  }

  return f8_result;
}

/*  
 *  $Log: bx_position_reco_noavg.cc,v $
 *  Revision 1.3  2015/08/17 17:24:22  misiaszek
 *  Id & log keywords added
 *
 */
