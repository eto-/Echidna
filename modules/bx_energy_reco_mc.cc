/* BOREXINO Reconstruction program
 *
 * Author: Marcin Misiaszek <marcin.misiaszek@lngs.infn.it>
 * Maintainer: Marcin Misiaszek <marcin.misiaszek@lngs.infn.it>
 *
 * $Id: bx_energy_reco_mc.cc,v 1.10 2009/10/26 14:41:38 misiaszek Exp $
 *
 * Implementation of bx_energy_reco_mc
 *
 */
#include "bx_energy_reco_mc.hh"
#include "messenger.hh"
//#include "db_channel.hh"
#include "bx_echidna_event.hh"
//#include "barn_interface.hh"
#include "interpolator2d.hh"

// ctor
bx_energy_reco_mc::bx_energy_reco_mc (): bx_base_module("bx_energy_reco_mc", bx_base_module::main_loop) {
  require_event_stage (bx_detector::laben, bx_base_event::clustered);
  require_trigger_type (bx_trigger_event::neutrino);
}


//calibration tables for nhits, charge and npe

const 									//clusters.nhits
 double bx_energy_reco_mc::nhits[5][3] = { { 79.6, 155.4, 300.2},    	//R = 0 , E = 0.2, 0.4, 0.8 MeV	
				  	   { 79.6, 155.5, 300.0},	//R = 1 m
				  	   { 79.2, 154.5, 296.9},	//R = 2 m
				  	   { 77.8, 151.0, 287.3},	//R = 3 m
				  	   { 72.3, 138.4, 258.3} };	//R = 4 m	

const									//clusters.charge
 double bx_energy_reco_mc::charge[5][3] ={ { 89.0, 177.2, 355.0},    	//R = 0 , E = 0.2, 0.4, 0.8 MeV	
				  	   { 89.1, 177.7, 356.5},	//R = 1 m
				  	   { 88.8, 177.5, 356.8},	//R = 2 m
				  	   { 87.8, 175.6, 355.6},	//R = 3 m
					   { 82.8, 166.9, 351.4} };	//R = 4 m	

const 									//clusters.npe
 double bx_energy_reco_mc::npe[5][3] =   { { 95.5, 189.0, 376.6},    	//R = 0 , E = 0.2, 0.4, 0.8 MeV	
				  	   { 95.6, 189.5, 377.5},	//R = 1 m
				  	   { 95.3, 189.2, 377.6},	//R = 2 m
				  	   { 94.1, 186.9, 375.4},	//R = 3 m
				  	   { 88.7, 177.2, 368.8} };	//R = 4 m	


// module interface
void bx_energy_reco_mc::begin () {

  static double r[5] = { 0., 1., 2., 3., 4. };          //radius of event
  static double e[3] = { 0.2, 0.4, 0.8 };               //energy

  const double (*n)[3];
  n = new double[5][3];

//choose calibration table for energy reconstruction  
  method = 2;    
  
  switch (method) {
    case 0: n = nhits;  break;
    case 1: n = charge; break;
    case 2: n = npe;    break;

  }
 			      
//  for( int32_t i = 0; i < 5.q; i++ )
//    for( int32_t j = 0; j < 3; j++ ) std::cout << r[i] << " " << e[j] << " " << n[i][j] << std::endl;

//initialization of interpolator2d object
   const int32_t M = 5, N = 3;
   int32_t i, j;
         			    
   interpolator2d::Vec_DP x1(M),x2(N);
   interpolator2d::Mat_DP y(M,N);
	    
   for (i=0;i<M;i++) x1[i]=r[i];
   for (i=0;i<N;i++) x2[i]=e[i];

   for (i=0;i<M;i++) {
	 for (j=0;j<N;j++) {
			     y[i][j]=n[i][j];
			   };
		     };
		     
   interpol = new interpolator2d(x1, x2, y);					       
			    
}

bx_echidna_event* bx_energy_reco_mc::doit (bx_echidna_event *ev) {
    // Loop on every cluster
    for (int32_t i = 0; i < ev->get_laben ().get_nclusters (); i++) {
        
	// Get rec_cluster reference
	//bx_laben_cluster& c = ev->get_laben ().get_cluster (i);
	const bx_laben_rec_cluster& rec_c = ev->get_laben ().get_rec_cluster (i);

	float n=0.;
	
	  switch (method) {
	    case 0: n = rec_c.get_rec_nhits (); break;
	    case 1: n = rec_c.get_cluster().get_charge(); break;
	    case 2: n = rec_c.get_cluster().get_npe();  break;
	  }
	
	// Get position from rec_cluster 
	//bx_energy_mc& e_mc = c.get_energy_mc ();
	const bx_position& rec_p = rec_c.get_position();
	
	float r = rec_p.get_r();	    	
	float e = m_get_energy(r,n);

        r = e = 0.; //remove warnings

	//std::cout <<  n  << " " <<  r  << " " << e << std::endl;
        //e_mc.f4_e = e;

    };
    return ev;
}




void bx_energy_reco_mc::end () {
    delete interpol;
}



bx_energy_reco_mc::~bx_energy_reco_mc () {    
}



float bx_energy_reco_mc::m_get_energy(float r, float nhits) {
    float e_low = 0.2 - 0.2;
    float e_up  = 0.8 + 0.2;
    float n_low = interpol->get_value(r, e_low);
    float n_up =  interpol->get_value(r, e_up);
    
    //std::cout << "low = " << n_low << " up = " << n_up << std::endl;

    if (nhits < n_low) return 0.0;
    if (nhits > n_up ) return 0.0;
    
    float e_index;
    
    while (e_up - e_low > 0.005) {
      e_index = (e_low + e_up) / 2.;
      if ( interpol->get_value(r,e_index) > nhits ) e_up = e_index;
      else e_low = e_index;
    }
    return e_index;
}

