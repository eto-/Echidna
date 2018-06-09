/* BOREXINO Reconstruction program
 *
 * Author: Werner Maneschg <werner.maneschg@mpi-hd.mpg.de>, Michael Wurm 
 * Maintainer: Davide Dangelo 
 *
 * $Id: bx_global_tracker.cc,v 1.28 2014/12/11 21:27:12 wurm Exp $
 *
 * Implemenentation of bx_global_tracker
 *
 * last modification: 2009-06-23
*/


// NOTE: For saving the new output-variables of the calculated axis intercept and slopes in the 2
// planes one has to define a new class similar to "bx_track.hh"; this class belong to the
// "event" (see: offline/Echidna/event$), which is maintained by ddangelo and only changed by
// him.
 
#include "bx_global_tracker.hh"
#include "bx_echidna_event.hh"
#include "db_channel.hh"
#include "constants.hh"
#include "bx_track.hh"


// NOTE: Selection of events which has clustered in the outer AND inner detector
// So, the bx_echidna_event contains automatically only laben_events which are correlated
// with a muon_event passed through the inner and outer detector

bx_global_tracker::bx_global_tracker (): bx_base_module("bx_global_tracker", bx_base_module::main_loop) {
  require_event_stage (bx_detector::laben, bx_base_event::clustered);
  require_event_stage (bx_detector::muon, bx_base_event::clustered);
  require_trigger_type (bx_trigger_event::neutrino);
}


// NOTE: Initialisation of variables which are declared in the hh-file; the variables and methods
// declared in the hh-file are typically interesting for the enduser. The enduser can change the
// default values in the echidna.cfg file.

void bx_global_tracker::begin () {
  i4_global_events=0;
  i4_tracked_events=0;

  f4_max_dev   = get_parameter ("max_dev")  .get_float();
  f4_omega     = get_parameter ("omega")    .get_float();
  f4_chi2max44 = get_parameter ("chi2max44").get_float();
  f4_max_par   = get_parameter ("max_par")  .get_float();
  i4_min_nhits = 2250;
}



bx_echidna_event* bx_global_tracker::doit (bx_echidna_event *ev) {

  i4_global_events++;

  const bx_track_by_points& t_m = ev->get_muon().get_track();
  bx_track_by_points t_l = ev->get_laben().get_track_energy();
  bool fit_ep_tof = bool(ev->get_laben().get_track_tof().get_x1() && ev->get_laben().get_decoded_nhits()>1e3 && ev->get_laben().get_decoded_nhits()<i4_min_nhits);
  if (ev->get_laben().is_tracked_tof() || fit_ep_tof ) t_l= ev->get_laben().get_track_tof(); //prefer tof over energy, use tof also for events between 1000-2250 hits when only the entry point is usable

  //get_message (bx_message::log) << "event " << ev->get_event_number () << dispatch;
  
  // Select the coordinates of the points

  Int_t total_pts = 4;		// Number of all possible reconstructed points in ID and OD (default 4)

  bool m1 = (t_m.get_dx1()!=0);
  bool m2 = (t_m.get_dx2()!=0);
  bool l1 = (t_l.get_dx1()!=0);
  bool l2 = (t_l.get_dx2()!=0);
  if ( fit_ep_tof ) {
    l2 = false; // use xp tof only if at least 2250 hits are registered
    //get_message (bx_message::log) << ev->get_laben().get_decoded_nhits() << " laben hits => laben XP is not included in fit"  << dispatch; 
  }

  // check whether laben n decoded hits and track impact parameters are consistent

  if (t_m.get_impact() > 7.2) { m1=false; m2=false; }
  if (ev->get_laben().get_decoded_nhits()>7500 && t_m.get_impact() > 6.0 && t_l.get_impact() < 6.0 && l1 && l2) { m1=false; m2=false; }
  if (ev->get_laben().get_decoded_nhits()>7500 && t_m.get_impact() < 6.0 && t_l.get_impact() > 6.0 && m1 && m2) { l1=false; l2=false; }





  Int_t all_pts = int(m1)+int(m2)+int(l1)+int(l2);  // number of reconstructed points
    
  //get_message (bx_message::log) << "available points: m1=" << m1 << " m2=" << m2 << " l1=" << l1 << " l2=" << l2 << " -> " << all_pts << dispatch; 
  if (all_pts<2) {
    //get_message (bx_message::log) << "not enough points on track. discard ..." << dispatch;
    return ev;
  }


  float x[total_pts], y[total_pts], z[total_pts];
  float ex[total_pts], ey[total_pts], ez[total_pts];
    
  int ip=0;  
  if (m1) {					// if present, write OD entry point
    x[ip] = t_m.get_x1();
    y[ip] = t_m.get_y1();
    z[ip] = t_m.get_z1();
    ex[ip] = .5;//t_m.get_dx1();
    ey[ip] = .5;//t_m.get_dy1();
    ez[ip] = .5;//t_m.get_dz1();     
    ip++;
  }

  if (l1) {					// if present, write ID entry point
    x[ip] = t_l.get_x1();
    y[ip] = t_l.get_y1();
    z[ip] = t_l.get_z1();
    ex[ip] = .5;//t_l.get_dx1();
    ey[ip] = .5;//t_l.get_dy1();
    ez[ip] = .5;//t_l.get_dz1();     
    ip++;
  }

  if (l2) {					// if present, write ID exit point
    x[ip] = t_l.get_x2();
    y[ip] = t_l.get_y2();
    z[ip] = t_l.get_z2();
    ex[ip] = .5;//t_l.get_dx2();
    ey[ip] = .5;//t_l.get_dy2();
    ez[ip] = .5;//t_l.get_dz2();     
    ip++;
  }

  if (m2) {					// if present, write OD exit point
    x[ip] = t_m.get_x2();
    y[ip] = t_m.get_y2();
    z[ip] = t_m.get_z2();
    ex[ip] = .5;//t_m.get_dx2();
    ey[ip] = .5;//t_m.get_dy2();
    ez[ip] = .5;//t_m.get_dz2();     
  }

  //for (int ip=0; ip<all_pts; ip++) {
    //get_message (bx_message::log) << "P" << ip << ": x=" << x[ip] << "(" << ex[ip] << "), y=" << y[ip] << "(" << ey[ip] << "), z=" << z[ip] << "(" << ez[ip] << ")" << dispatch; //hu? 
  //}

// retrieve track from m_fit_track method

  float *par = m_fit_track(x, y, z, ex, ey, ez, all_pts, all_pts);
  float chi2   = par[0];
  float alpha  = par[1];
  float beta   = par[2];
  float gamma  = par[3];
  float delta  = par[4];
  float ealpha = par[5];
  float ebeta  = par[6];
  float egamma = par[7];
  float edelta = par[8];
  int   ndf    = (all_pts-1)*2;
  float rchi2  = chi2/ndf;

  /*get_message(bx_message::log) << "FIT OF " << all_pts << " POINTS:" << dispatch;
  get_message(bx_message::log) << "alpha = " << alpha << ", ealpha = " << ealpha << dispatch;
  get_message(bx_message::log) << "beta = " << beta << ", ebeta = " << ebeta << dispatch;
  get_message(bx_message::log) << "gamma = " << gamma << ", egamma = " << egamma << dispatch;
  get_message(bx_message::log) << "delta = " << delta << ", edelta = " << edelta << dispatch; 
  get_message(bx_message::log) << "chi2 = " << chi2 << "/" << ndf << " (" << chi2/ndf << ")" << dispatch;
 */
 
// if fit is bad, check if there is a better fit with 3 points
// Declaration of reconstructed points which are finally used in the fit

/*bool p1 = true; // Entrypoint OD
bool p2 = true; // Exitpoint OD
bool p3 = true; // Entrypoint ID
bool p4 = true; // Exitpoint ID

if(all_pts < 4){
  if(m1 = bool(t_m.get_dx1()<0.0000001)) {p1 = false;};
  if(m2 = bool(t_m.get_dx2()<0.0000001)) {p2 = false;};
  if(l1 = bool(t_l.get_dx1()<0.0000001)) {p3 = false;};
  if(l2 = bool(t_l.get_dx2()<0.0000001)) {p4 = false;};
  } 
*/

  float fit_nr = all_pts;
  if (all_pts==total_pts && rchi2>f4_chi2max44) {
    float chi2_min = rchi2;
     for (int i=0; i<all_pts; i++) {
       par = m_fit_track(x, y, z, ex, ey, ez, all_pts, i);
       if (par[0]/(all_pts-2)/2 < chi2_min) {
         chi2_min = par[0];
	 fit_nr   = i;
         chi2   = par[0];
         alpha  = par[1];
         beta   = par[2];
         gamma  = par[3];
         delta  = par[4];
         ealpha = par[5];
         ebeta  = par[6];
         egamma = par[7];
         edelta = par[8];
	 ndf    = (all_pts-2)*2;
	 rchi2  = chi2/ndf;
	 //if (i == 0) {p1 = false;} 
	 //if (i == 1) {p2 = false;}
	 //if (i == 2) {p3 = false;}
	 //if (i == 3) {p4 = false;}
	 /*get_message(bx_message::log) << "BETTER FIT WITHOUT POINT " << i << dispatch;
	 get_message(bx_message::log) << "alpha = " << alpha << ", ealpha = " << ealpha << dispatch;
         get_message(bx_message::log) << "beta = " << beta << ", ebeta = " << ebeta << dispatch;
         get_message(bx_message::log) << "gamma = " << gamma << ", egamma = " << egamma << dispatch;
         get_message(bx_message::log) << "delta = " << delta << ", edelta = " << edelta << dispatch; 
         get_message(bx_message::log) << "chi2 = " << chi2 << "/" << ndf << " (" << chi2/ndf << ")" << dispatch;
*/      }
    }
  }

  if (fit_nr != all_pts) {
    //get_message(bx_message::log) << "Fit of 3 points returned better result. Used instead ... " << dispatch;
    if (fit_nr==0) m1 = false;
    if (fit_nr==1) l1 = false;
    if (fit_nr==2) l2 = false;
    if (fit_nr==3) m2 = false;
    all_pts--;
  }


// Finally used reconstructed points for the fit

  Char_t points = 0x00;
  if(m1)  {points = points | 0x01;} //0001
  if(l1)  {points = points | 0x02;} //0010
  if(l2)  {points = points | 0x04;} //0100
  if(m2)  {points = points | 0x08;} //1000


// check for tracks without converging fit

 if (fabs(alpha)>f4_max_par || fabs(beta)>f4_max_par || fabs(gamma)>f4_max_par || fabs(delta)>f4_max_par) {
   //get_message(bx_message::log) << "Too large parameters. Fit did not converge. Skip the event ..." << dispatch;
   return ev;
 }



// Define the direction of the muon: upward or downward
// use the two most remote track points
  bool muon_downward = false;
  int entry_nr = 0, exit_nr = 3;
  if (!m2) {
    exit_nr--;
    if (!l2) exit_nr--;
  }
  if (!m1) {
    entry_nr++;
    if (!l1) entry_nr++;
  }
  if (exit_nr == entry_nr) return ev;
 
  // if there is doubt if the information of OD (which is mainly used for determining track direction) is correct,
  // rely on ID reco only

  if (ev->get_laben().get_nclusters()) {
    if (  (t_m.get_t1()-ev->get_laben().get_cluster(0).get_start_time()+ev->get_laben().get_trigger_rawt()<-35.) ||
  	  (!t_m.get_dx2()) || (t_m.get_t2()<t_m.get_t1()) || (t_m.get_z1()<-5.) ) {
      entry_nr = 1;
      exit_nr = 2;
    }
  }

  if (fabs(delta)>1.) {    // track is steep, z coordinates are reliable
    if (z[exit_nr] < z[entry_nr]) muon_downward = true;
  }
  else {     // track is gentle, not good for z
    if (fabs(beta)>1.) {     // track is steep in xy plane, y coordinate better than x
      if ( (y[exit_nr]-y[entry_nr])*beta*delta < 0 ) muon_downward = true;
    }
    else {     // x coordinate better
      if ( (x[exit_nr]-x[entry_nr])*delta < 0 ) muon_downward = true;
    }
  }


// Fill variables

  bx_track_fitted& t = ev->get_track_global(); 
  t.f8_alpha         = alpha;
  t.f8_alpha_error   = fabs(ealpha);
  t.f8_beta          = beta;
  t.f8_beta_error    = fabs(ebeta);
  t.f8_gamma         = gamma;
  t.f8_gamma_error   = fabs(egamma);
  t.f8_delta         = delta;
  t.f8_delta_error   = fabs(edelta);
  t.f4_chi2          = chi2/ndf;
  t.u1_points	     = points;
  t.b_downward       = muon_downward;

  float phi = atan2(beta,1.);
  if ( (delta<0 && muon_downward) || (delta>0 && !muon_downward) ) phi += constants::number::pi;
  while ( phi<0 ) phi += 2.*constants::number::pi;
  while ( phi>2.*constants::number::pi ) phi -= 2.*constants::number::pi;

  float theta = acos(delta/sqrt(1.+beta*beta+delta*delta));
  if (theta < 0) theta += constants::number::pi;
  if ( (delta<0 && muon_downward) || (delta>0 && !muon_downward) ) theta = constants::number::pi - theta;

  // if only OD or ID tracking is used for determining track direction, compare with the angles derived from there

  /*if (entry_nr==0 && exit_nr==3) {
    if (fabs(phi-t_m.get_phi()) > 0.5*constants::number::pi && fabs(phi-t_m.get_phi()) < 1.5*constants::number::pi) {
      phi += constants::number::pi;
      if (phi > 2.*constants::number::pi) phi -= 2.*constants::number::pi;
      theta = constants::number::pi-theta;
    }
  }

 if (entry_nr==1 && exit_nr==2) {
    if (fabs(phi-t_l.get_phi()) > 0.5*constants::number::pi && fabs(phi-t_l.get_phi()) < 1.5*constants::number::pi) {
      phi += constants::number::pi;
      if (phi > 2.*constants::number::pi) phi -= 2.*constants::number::pi;
      theta = constants::number::pi-theta;
    }
  }*/

  float x0 = -(alpha*beta+gamma*delta) / (1.+beta*beta+delta*delta);
  float impact = sqrt( x0*x0 + pow(alpha+beta*x0, 2) + pow(gamma+delta*x0, 2) );

  t.f4_phi     = phi;
  t.f4_theta   = theta;
  t.f4_impact  = impact;
  t.f4_dphi    = 0.;
  t.f4_dtheta  = 0.;
  t.f4_dimpact = 0.;
  
  if (fabs(ealpha)<1e-10 && fabs(ebeta)<1e-10 && fabs(egamma)<1e-10 && fabs (edelta<1e-10)) { get_message(bx_message::warn) << "Fit to track did not converge! Skip event ..." << dispatch; return ev; } 

  // write normalized number of decoded hits into track variable
  if ( ev->get_laben().get_nclusters()>0 )
    t.f4_labennormhits = ev->get_laben().get_decoded_nhits() / (ev->get_laben().get_n_live_pmts() - ev->get_laben().get_invalid_pmts()) * 2000.;
           
  t.b_is_valid       = true;
 
  i4_tracked_events++; 
  
  // ONLY FOR TEST: UNCOMMENT AFTER FINAL TESTS
  /*get_message(bx_message::debug) << "FINAL TRACK:" << dispatch;
  get_message(bx_message::debug) << "alpha = " << alpha << ", ealpha = " << ealpha << dispatch;
  get_message(bx_message::debug) << "beta = " << beta << ", ebeta = " << ebeta << dispatch;
  get_message(bx_message::debug) << "gamma = " << gamma << ", egamma = " << egamma << dispatch;
  get_message(bx_message::debug) << "delta = " << delta << ", edelta = " << edelta << dispatch; 
  get_message(bx_message::debug) << "chi2 = " << chi2 << dispatch;
  */
  return ev;  
}


void bx_global_tracker::end () {

  if (i4_global_events > 0) get_message (bx_message::info) << "run summary: " << i4_tracked_events
<< " of " << i4_global_events << " events were reconstructed." << dispatch;
}	



float* bx_global_tracker::m_fit_track(float *px, float *py, float *pz, float *pdx, float *pdy, float *pdz, int npts, int bpt) {

  //get_message(bx_message::debug) << "fitting " << npts << " points, omitting point " << bpt << dispatch;

  float omega = f4_omega/180.*3.14152;

  double x_arr[npts];
  double y_arr[npts];
  double z_arr[npts];
  double dx_arr[npts];
  double dy_arr[npts];
  double dz_arr[npts];
  int jp = 0;
  int gpts = npts;
  if (bpt<npts) gpts--;
  for (int ip=0; ip<npts; ip++) {
    if (bpt==ip) continue;
    x_arr[jp]  = px[ip];
    y_arr[jp]  = py[ip];
    z_arr[jp]  = pz[ip];
    dx_arr[jp] = pdx[ip];
    dy_arr[jp] = pdy[ip];
    dz_arr[jp] = pdz[ip];
    jp++;
  }

  //for (int i=0; i<npts; i++) get_message(bx_message::debug) << i << "th point: "<< x_arr[i] << " : " << y_arr[i] << " : " << z_arr[i] << dispatch;

  bool rotate_y=false, rotate_z=false;
  double c = cos(omega);
  double s = sin(omega);

  // rotate around y-axis
  if (fabs((z_arr[gpts-1]-z_arr[0])/(x_arr[gpts-1]-x_arr[0]))>f4_max_dev) {
    rotate_y = true;
    //get_message(bx_message::debug) << "rotating around y-axis ..." << dispatch;
    for (int i=0; i<gpts; i++) {
      float x_temp = x_arr[i], dx_temp = dx_arr[i];
      x_arr[i] = c*x_temp + s*z_arr[i];
      z_arr[i] = -s*x_temp + c*z_arr[i];
      dx_arr[i] = c*dx_temp + s*dz_arr[i];
      dz_arr[i] = -s*dx_temp + c*dz_arr[i];
      // get_message(bx_message::debug) << x_arr[i] << " : " << y_arr[i] << " : " << z_arr[i] <<dispatch;
    }
  }

    // rotate around z-axis
  if (fabs((y_arr[gpts-1]-y_arr[1])/(x_arr[gpts-1]-x_arr[0]))>f4_max_dev) {
    rotate_z = true;
    //get_message(bx_message::debug) << "rotating around z-axis ..." << dispatch;
    for (int i=0; i<gpts; i++) {
      float x_temp = x_arr[i], dx_temp = dx_arr[i];
      x_arr[i] = c*x_temp + s*y_arr[i];
      y_arr[i] = -s*x_temp + c*y_arr[i];
      dx_arr[i] = sqrt( pow(c*dx_temp,2) + pow(s*dz_arr[i],2) );
      dy_arr[i] = sqrt( pow(-s*dx_temp,2) + pow(c*y_arr[i],2) );
    }
  } 

  TGraphErrors g_xy(gpts, x_arr, y_arr, dx_arr, dy_arr);
  TGraphErrors g_xz(gpts, x_arr, z_arr, dx_arr, dz_arr);

  TF1 fit_xy("fitxy","[0]+[1]*x",-10,10);
  Double_t fit_par_xy = (y_arr[gpts-1]-y_arr[0])/(x_arr[gpts-1]-x_arr[0]);
  //get_message(bx_message::debug) << "xy start: alpha=" << fit_par_xy << ", beta=" << z_arr[0]-fit_par_xy*x_arr[0] << dispatch;
  fit_xy.SetParameter(1,fit_par_xy);  
  fit_xy.SetParameter(0,y_arr[0]-fit_par_xy*x_arr[0]);
  g_xy.Fit("fitxy","Q0");
  double alpha  = fit_xy.GetParameter(0);
  double beta   = fit_xy.GetParameter(1);
  double dalpha = fit_xy.GetParError(0);
  double dbeta  = fit_xy.GetParError(1);
  float chi2  = fit_xy.GetChisquare();
  TF1 fit_xz("fitxz","[0]+[1]*x",-10,10);
  Double_t fit_par_xz = (z_arr[gpts-1]-z_arr[0])/(x_arr[gpts-1]-x_arr[0]);
  //get_message(bx_message::debug) << "xz start: gamma=" << fit_par_xz << ", delta=" << z_arr[0]-fit_par_xz*x_arr[0] << dispatch;
  fit_xz.SetParameter(1,fit_par_xz);
  fit_xz.SetParameter(0,z_arr[0]-fit_par_xz*x_arr[0]);
  g_xz.Fit("fitxz","Q0");
  double gamma = fit_xz.GetParameter(0);
  double delta = fit_xz.GetParameter(1);
  double dgamma = fit_xy.GetParError(0);
  double ddelta = fit_xy.GetParError(1);
  chi2 += fit_xz.GetChisquare();
  //get_message(bx_message::debug) << "final: alpha=" << alpha << ", beta=" << beta << ", gamma=" << gamma << ", delta=" << delta << dispatch;

  if (rotate_z) {
    //get_message(bx_message::debug) << "rotating back z ... " << dispatch;
    double temp_alpha = alpha;
    double temp_beta = beta;
    alpha = c*alpha+(s*alpha*(s+c*beta))/(c-s*beta);
    gamma = gamma+(s*temp_alpha*delta)/(c-s*beta);
    beta  = (s+c*beta)/(c-s*beta);
    delta = delta/(c-s*temp_beta);
    double temp_dalpha = dalpha;
    double temp_dbeta = dbeta;
    dalpha = c*dalpha+(s*dalpha*(s+c*dbeta))/(c-s*dbeta);
    dgamma = dgamma+(s*temp_dalpha*ddelta)/(c-s*dbeta);
    dbeta  = (s+c*dbeta)/(c-s*dbeta);
    ddelta = ddelta/(c-s*temp_dbeta);
  }

  if (rotate_y) {
    //bx_message(bx_message::debug) << "rotating back y ..." << dispatch;
    alpha = alpha+(s*beta*gamma)/(c-s*delta);
    gamma = c*gamma+(s*gamma*(s+c*delta))/(c-s*delta);
    beta  = beta/(c-s*delta);
    delta = (s+c*delta)/(c-s*delta);
    dalpha = dalpha+(s*dbeta*dgamma)/(c-s*ddelta); 
    dgamma = c*dgamma+(s*dgamma*(s+c*ddelta))/(c-s*ddelta);
    dbeta  = dbeta/(c-s*ddelta);
    ddelta = (s+c*ddelta)/(c-s*ddelta);
  }  

  float *parameters = new float[9];
  parameters[0] = chi2;
  parameters[1] = alpha;
  parameters[2] = beta;
  parameters[3] = gamma;
  parameters[4] = delta;
  parameters[5] = dalpha;
  parameters[6] = dbeta;
  parameters[7] = dgamma;
  parameters[8] = ddelta;
  
  return parameters;
}

/*
 * $Log: bx_global_tracker.cc,v $
 * Revision 1.28  2014/12/11 21:27:12  wurm
 * added getters and variables to BxEvent that concern C-11 likelihood tagging
 * * for tracks: new distance-to- and projection-on-track getters, dE/dx calculation
 * * counting of neutron-like clusters in tt128, new flag of neutron clusters
 *
 * Revision 1.27  2011/03/11 11:59:04  wurm
 * added some cross-checks on laben nhits and impact parameter of OD/ID tracking
 *
 * Revision 1.26  2011-03-09 10:29:30  wurm
 * lowered verbosity
 *
 * Revision 1.25  2011-03-08 18:22:05  wurm
 * bug-fix on muon-laben-consistency check
 *
 * Revision 1.24  2010-11-30 10:36:34  wurm
 * small change to phi/theta calculation
 *
 * Revision 1.23  2010-10-05 13:01:56  wurm
 * fixed points variable and improved phi calculation
 *
 * Revision 1.22  2010-08-26 08:25:09  wurm
 * fixed theta variable
 *
 * Revision 1.21  2010-08-05 15:23:08  wurm
 * changed the conditions for laben.is_tracked_tof, adjusted global tracker
 *
 * Revision 1.20  2010-08-04 08:29:35  wurm
 * removed unneccesary variables, set uncertainties for global tracker to 0.5 m
 *
 * Revision 1.19  2010-08-03 15:57:34  wurm
 * 1) introduced theta, phi, impact variables for fitted tracks
 * 2) set fixed uncertainties for global tracking
 * 3) condition for execution of laben tracking changed to MTF OR MCF
 * 4) lowered verbosity of all modules
 *
 * Revision 1.18  2010-07-27 15:00:58  wurm
 * changed laben track from reference to variable
 *
 * Revision 1.17  2010-06-28 10:44:02  wurm
 * new laben tof tracking, modifications to muon tracking (introducing phi,theta,impact to the event) and modification of global tracker to prefer laben tof over laben energy (and conditions)
 *
 * Revision 1.16  2010-05-21 16:09:12  ddangelo
 * patched to comply with new class names
 *
 * Revision 1.15  2009-10-23 10:24:48  wurm
 * mute
 *
 * Revision 1.14  2009-10-23 10:09:13  wurm
 * removed warnings
 *
 * Revision 1.13  2009-10-23 09:29:00  wurm
 * commented debug messages
 *
 * Revision 1.12  2009-09-30 10:02:22  wurm
 * improved downward flag for horizontal tracks
 *
 * Revision 1.11  2009-08-03 10:07:04  wurm
 * fill the is_valid variable
 *
 * Revision 1.10  2009-06-24 16:06:13  maneschg
 * Add bitfield for reconstructed points used for the global fit
 *
 * Revision 1.9  2009-04-21 07:30:15  maneschg
 * Old variable for used_points commented out. New corrected variable not yet operational.
 *
 * Revision 1.8  2008-12-15 14:32:24  wurm
 * lowered verbosity
 *
 * Revision 1.7  2008-12-09 09:10:37  wurm
 * bug fix
 *
 * Revision 1.6  2008-12-04 15:38:25  wurm
 * added large coordinates rotations, still very verbose
 *
 * Revision 1.5  2008-08-10 14:16:57  maneschg
 * Add approach to discard non-well reconstructed exit points
 *
 * Revision 1.4  2008-08-06 09:50:03  maneschg
 * Discrimination of non reconstructed points and adding some new functions
 *
 * Revision 1.3  2008-07-11 15:27:44  maneschg
 * add new module with global fitter for muontrack-reconstruction
 *
 * Revision 1.2  2008-04-29 14:05:54  ddangelo
 * fixed the name
 *
 * Revision 1.1  2008-04-29 13:45:55  ddangelo
 * added empty
 *
 * 
 */

