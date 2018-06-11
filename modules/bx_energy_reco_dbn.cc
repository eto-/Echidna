/* BOREXINO Reconstruction program
 *
 * Authors:  A.Derbin, V.Muratova, O.Smirnov 
 * Maintainer: A.Derbin, V.Muratova,O.Smirnov
 *
 * Implementation of bx_energy_reco_dbn
 *
 */
#include <TH1F.h>
#include "bx_energy_reco_dbn.hh"
#include "bx_echidna_event.hh"
#include "barn_interface.hh"

#define nbins 55
#define NearlyZero 1E-10

// ctor
bx_energy_reco_dbn::bx_energy_reco_dbn (): bx_base_module("bx_energy_reco_dbn", bx_base_module::main_loop) {
  require_event_stage (bx_detector::laben, bx_base_event::clustered);
  require_trigger_type (bx_trigger_event::neutrino);
}


//Tables of smoothing spline coefficients (MC for 400 keV)
//First index: Q==0,Npe==1,Nhits==2
const double bx_energy_reco_dbn::SSpline_Y[3][55]={
          {1.003026e+00,1.003014e+00,1.002975e+00,1.002903e+00,1.002779e+00,
           1.002569e+00,1.002245e+00,1.001878e+00,1.001539e+00,1.001264e+00,
           1.001024e+00,1.000799e+00,1.000557e+00,1.000364e+00,1.000261e+00,
           1.000230e+00,1.000221e+00,1.000176e+00,1.000045e+00,9.997661e-01,
           9.992422e-01,9.984395e-01,9.974210e-01,9.961780e-01,9.947008e-01,
           9.930738e-01,9.914269e-01,9.897485e-01,9.880473e-01,9.862764e-01,
           9.843652e-01,9.821678e-01,9.794756e-01,9.760789e-01,9.717669e-01,
           9.664576e-01,9.601846e-01,9.530402e-01,9.452048e-01,9.369215e-01,
           9.285433e-01,9.204675e-01,9.130162e-01,9.063740e-01,9.006844e-01,
           8.960587e-01,8.925471e-01,8.925471e-01,8.925471e-01,8.925471e-01,
           8.925471e-01,8.925471e-01,8.925471e-01,8.925471e-01,8.925471e-01},
          {1.002446e+00,1.002438e+00,1.002413e+00,1.002363e+00,1.002272e+00,
           1.002108e+00,1.001846e+00,1.001559e+00,1.001313e+00,1.001131e+00,
           1.000988e+00,1.000863e+00,1.000695e+00,1.000521e+00,1.000362e+00,
           1.000209e+00,1.000034e+00,9.998042e-01,9.994892e-01,9.990559e-01,
           9.984154e-01,9.975171e-01,9.964054e-01,9.950698e-01,9.935070e-01,
           9.917887e-01,9.900266e-01,9.881956e-01,9.862940e-01,9.842908e-01,
           9.821329e-01,9.796847e-01,9.767505e-01,9.731188e-01,9.685606e-01,
           9.629732e-01,9.563766e-01,9.488632e-01,9.406114e-01,9.318686e-01,
           9.230143e-01,9.144583e-01,9.065193e-01,8.993839e-01,8.931958e-01,
           8.880756e-01,8.840911e-01,8.840911e-01,8.840911e-01,8.840911e-01,
           8.840911e-01,8.840911e-01,8.840911e-01,8.840911e-01,8.840911e-01},
	  {1.001786e+00,1.001816e+00,1.001902e+00,1.002031e+00,1.002180e+00,
           1.002309e+00,1.002377e+00,1.002413e+00,1.002444e+00,1.002454e+00,
           1.002389e+00,1.002213e+00,1.001877e+00,1.001448e+00,1.000974e+00,
           1.000425e+00,9.997361e-01,9.988768e-01,9.978491e-01,9.966367e-01,
           9.951888e-01,9.934724e-01,9.915416e-01,9.893928e-01,9.869982e-01,
           9.843910e-01,9.816759e-01,9.788440e-01,9.759023e-01,9.728109e-01,
           9.694800e-01,9.657201e-01,9.612960e-01,9.559799e-01,9.495365e-01,
           9.418544e-01,9.329526e-01,9.229323e-01,9.120249e-01,9.005680e-01,
           8.890370e-01,8.779189e-01,8.676252e-01,8.584195e-01,8.505031e-01,
           8.440144e-01,8.390148e-01,8.390148e-01,8.390148e-01,8.390148e-01,
           8.390148e-01,8.390148e-01,8.390148e-01,8.390148e-01,8.390148e-01}};	   

const double bx_energy_reco_dbn::SSpline_M[3][55]={
          {-2.312826e-05,-2.601709e-05,-3.178272e-05,-4.862312e-05,-8.355785e-05,
           -1.392121e-04,-3.761327e-05,3.142430e-05,8.146408e-05,2.386112e-05,
            3.062154e-05,-4.980373e-05,6.069130e-05,1.018337e-04,7.551924e-05,
            2.389218e-05,-3.704569e-05,-8.704677e-05,-1.360411e-04,-2.566778e-04,
           -3.063080e-04,-1.905319e-04,-2.261131e-04,-2.525947e-04,-1.681367e-04,
            2.563220e-05,-5.308577e-05,-2.621006e-06,-7.305098e-05,-1.237362e-04,
           -2.734809e-04,-4.998513e-04,-6.956962e-04,-9.445379e-04,-1.017763e-03,
           -9.681761e-04,-8.913422e-04,-6.955321e-04,-4.722941e-04,-1.030325e-04,
            3.155109e-04,6.551660e-04,8.107346e-04,9.570806e-04,1.075996e-03,
            1.122269e-03,1.120042e-03,0.000000e+00,0.000000e+00,0.000000e+00,
            0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00},
          {-1.463187e-05,-1.732250e-05,-2.269600e-05,-3.793334e-05,-7.050144e-05,
           -1.232141e-04,-1.885734e-05,4.600036e-05,7.757805e-05,3.075631e-05,
            3.659452e-05,-7.389691e-05,3.775415e-06,2.056153e-05,8.821717e-06,
           -2.226500e-05,-5.419853e-05,-8.870704e-05,-1.017936e-04,-2.139957e-04,
           -2.853090e-04,-1.917546e-04,-2.277571e-04,-2.405790e-04,-1.731777e-04,
           -1.147875e-07,-8.893467e-05,-5.757609e-05,-1.041318e-04,-1.360812e-04,
           -2.792733e-04,-4.881596e-04,-6.850166e-04,-9.560188e-04,-1.050466e-03,
           -1.017331e-03,-9.350036e-04,-7.434293e-04,-5.219101e-04,-1.153016e-04,
            3.148557e-04,6.451754e-04,8.066498e-04,9.498391e-04,1.077772e-03,
            1.146469e-03,1.150515e-03,0.000000e+00,0.000000e+00,0.000000e+00,
            0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00},
           {6.008790e-05,5.775239e-05,4.499711e-05,2.287729e-05,-1.653282e-05,
           -7.924088e-05,-2.846352e-05,1.259834e-06,-9.779450e-06,-8.632139e-05,
           -9.608932e-05,-1.961069e-04,-8.252274e-05,-2.796875e-05,-6.988503e-05,
           -1.485196e-04,-1.769230e-04,-1.641907e-04,-1.770763e-04,-2.351622e-04,
           -2.954733e-04,-1.940057e-04,-2.146468e-04,-2.556562e-04,-2.375686e-04,
           -6.913202e-05,-1.341689e-04,-9.443958e-05,-1.467714e-04,-2.171673e-04,
           -4.213696e-04,-6.707696e-04,-8.815548e-04,-1.154788e-03,-1.262707e-03,
           -1.227230e-03,-1.146182e-03,-8.988406e-04,-5.815432e-04,-7.155143e-05,
            4.228215e-04,8.578679e-04,1.092258e-03,1.300397e-03,1.442451e-03,
            1.496336e-03,1.506171e-03,0.000000e+00,0.000000e+00,0.000000e+00,
            0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00}};



// module interface
void bx_energy_reco_dbn::begin () {

  pHisto_nhits = new TH1F("pHisto_nhits","Corrected hits",200,0.,600.);			    
  barn_interface::get ()->store (barn_interface::file, pHisto_nhits, this);
  pHisto_npe  = new TH1F("pHisto_npe","Corrected npe",200,0.,600.);
  barn_interface::get ()->store (barn_interface::file, pHisto_npe, this);		    
  pHisto_pe   = new TH1F("pHisto_pe","Corrected pe",200,0.,600.);
  barn_interface::get ()->store (barn_interface::file, pHisto_pe, this);  			    

}

bx_echidna_event* bx_energy_reco_dbn::doit (bx_echidna_event *ev) {
  
  float nhits,npe,pe,r;
  double corr;  
  const double fRadialBinning = 0.1;				       
	  
  for (int32_t i = 0; i < ev->get_laben ().get_nclusters (); i++) {
        
    //bx_laben_cluster& c = ev->get_laben ().get_cluster (i);
    const bx_laben_rec_cluster& rec_c = ev->get_laben ().get_rec_cluster (i);
    
    nhits = rec_c.get_rec_nhits ();       //Nhits
    pe    = rec_c.get_cluster().get_charge(); //Q
    npe   = rec_c.get_cluster().get_npe();//Npe  
    // Get position 
    const bx_position& recp_dbn = rec_c.get_position();
	
    r = recp_dbn.get_r()/fRadialBinning;	    	
    //get_message(bx_message::debug)<<" r="<<recp_dbn.get_r()<<dispatch;
    
    //bx_energy_reco_dbn& p=ev->get_laben ().get_cluster (i).get_energy_reco_dbn;   
    
    Spline(nbins,r,&SSpline_Y[0][0],&SSpline_M[0][0],&corr);  
    f4_pe_dbn=pe/corr;
    pHisto_pe->Fill(f4_pe_dbn,1.);

    //get_message(bx_message::debug)<<"corr pe="<<corr<<dispatch;

    Spline(nbins,(double)r,&SSpline_Y[1][0],&SSpline_M[1][0],&corr);  
    f4_npe_dbn=npe/corr;   
    pHisto_npe->Fill(f4_npe_dbn,1.);

    //get_message(bx_message::debug)<<"corr npe="<<corr<<dispatch;

    Spline(nbins,(double)r,&SSpline_Y[2][0],&SSpline_M[2][0],&corr);  
    f4_nhits_dbn=nhits/corr;        
    pHisto_nhits->Fill(f4_nhits_dbn,1.);

    //get_message(bx_message::debug)<<"corr nhits="<<corr<<dispatch;
    }
  return ev;
}




void bx_energy_reco_dbn::end () {

}



bx_energy_reco_dbn::~bx_energy_reco_dbn () {    
}



float bx_energy_reco_dbn::Spline(int32_t N, float X , const double *data,const double *m,double *y) {
  int32_t i;
  double dX,dY;
  double a,b,c;

  if ((X>(double)N)) {*y=data[N-1];return (data[N-1]);}
  if ((X<0.)) {*y=data[0];return (data[0]);}
  
  i=(int32_t)floor(X);

  dX=(double)X-(double)i;
  if (dX<=NearlyZero) {
    *y=data[i];
    return(data[i]);
    }

  dY=data[i+1]-data[i];
  a=(m[i+1]-m[i])/6.;
  b=m[i]/2.;
  c=dY-m[i]/2.-(m[i+1]-m[i])/6.;
  *y=a*dX*dX*dX+b*dX*dX+c*dX+data[i];
  return ((float)*y);
  }
