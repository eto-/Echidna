/* BOREXINO Reconstruction program
 *
 * Author: Michael Wurm <mwurm@ph.tum.de>, Davide Franco <davide.franco@mi.infn.it> 
 * Maintainer: Davide D'Angelo <davide.dangelo@mi.infn.it>
 *
 * $Id: bx_laben_energy_tracker.cc,v 1.7 2014/12/11 21:27:12 wurm Exp $
 *
 * Implemenentation of bx_laben_energy_tracker
 *
*/

#include "bx_laben_energy_tracker.hh"
#include "bx_echidna_event.hh"
//#include "bx_dbi.hh"
//#include "db_profile.hh"
#include "db_channel.hh"
//#include "db_calib.hh"
#include "constants.hh"
//#include "barn_interface.hh"
#include "bx_track.hh"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TF1.h"

static Float_t pi = constants::number::pi;

bx_laben_energy_tracker::bx_laben_energy_tracker (): bx_base_module("bx_laben_energy_tracker", bx_base_module::main_loop) {
  require_event_stage (bx_detector::laben, bx_base_event::clustered);
//  require_event_stage (bx_detector::muon, bx_base_event::clustered);
  require_trigger_type (bx_trigger_event::neutrino);
}

void bx_laben_energy_tracker::begin () {
  i4_laben_events = 0;
  i4_tracked_events = 0;

  i4_enable_histos  			= get_parameter ("enable_histos")				.get_int  ();
  f4_tau            			= get_parameter ("tau")          				.get_float();
  f4_time_limit        			= get_parameter ("time_limit")         				.get_float();  
  f4_long_limit = 40.;
  i4_nsteps = 8;
}


bx_echidna_event* bx_laben_energy_tracker::doit (bx_echidna_event *ev) {
  
  const bx_laben_event& er = ev->get_laben();
  
  if (ev->get_trigger().get_btb_inputs()!=4 && !ev->get_muon().has_cluster()) return ev;

  //get_message(bx_message::debug) << ev->get_event_number() << " >>> ID tracking" << dispatch;
  
  i4_laben_events++;
  
  bx_laben_tracked_event &ew = dynamic_cast<bx_laben_tracked_event&>(ev->get_laben());

  if (!er.get_nclusters()) { // work only if at least 1 cluster is present
    //get_message(bx_message::debug) << "no ID cluster! skip ..." << dispatch;
    return ev;
  }

  Int_t nhits = er.get_cluster(0).get_clustered_nhits();
  std::vector<Float_t> x_v, y_v, z_v, t_v;
  std::vector<Int_t>   s_v;
  Int_t entryhits=0;

  for (Int_t i = 0; i<nhits; i++) {
    const bx_laben_clustered_hit& chit = er.get_cluster(0).get_clustered_hit(i); 
    if ((chit.get_time() > f4_long_limit) || (chit.get_order_in_channel() > 1)) continue;  //AAA 0-based or 1-based ??? to be checked
    const db_channel_laben* c = chit.get_decoded_hit().get_db_channel();
    s_v.push_back(Int_t(chit.get_time()/f4_time_limit));
    x_v.push_back(c->pmt_x());
    y_v.push_back(c->pmt_y());
    z_v.push_back(c->pmt_z());
    t_v.push_back(exp(-chit.get_time()/f4_tau));
    entryhits++;
  }
  Int_t nhits_used = t_v.size();

  Float_t entry_x = 0., entry_y = 0., entry_z = 0., entry_r=0., weight = 0.;
  for(Int_t i = 0; i < nhits_used; i++) {
    if (s_v[i]>1) continue;
    entry_x += x_v[i];//*t_v[i];
    entry_y += y_v[i];//*t_v[i];
    entry_z += z_v[i];//*t_v[i];
    entry_r += sqrt(pow(x_v[i],2)+pow(y_v[i],2)+pow(z_v[i],2));//*t_v[i];
    weight += 1;//t_v[i];
  }

  Float_t entry_length = sqrt( entry_x*entry_x + entry_y*entry_y + entry_z*entry_z ); 

  entry_r /= weight;
  entry_x *= entry_r/entry_length;
  entry_y *= entry_r/entry_length;
  entry_z *= entry_r/entry_length;

  Float_t entry_dx = 0., entry_dy = 0., entry_dz = 0.;
  for(Int_t i = 0; i < nhits_used; i++) {
    entry_dx += pow(entry_x - x_v[i],2);//*t_v[i];    
    entry_dy += pow(entry_y - y_v[i],2);//*t_v[i];    
    entry_dz += pow(entry_z - z_v[i],2);//*t_v[i];    
  }

  entry_dx = sqrt(entry_dx/weight/entryhits+0.3);
  entry_dy = sqrt(entry_dy/weight/entryhits+0.3);
  entry_dz = sqrt(entry_dz/weight/entryhits+0.3);

  //get_message(bx_message::debug) << "ID entrypoint:  x " << entry_x << "+/-" << entry_dx << ", y " << entry_y << "+/-" << entry_dy << ", z " << entry_z << "+/-" << entry_dz << ", t " << er.get_cluster(0).get_start_time() - er.get_trigger_rawt() << dispatch;

  // filling event variables
  bx_track_by_points& t = ew.get_track_energy(); 
  t.f4_t1 = er.get_cluster(0).get_start_time() - er.get_trigger_rawt();
  t.f4_x1 = entry_x;
  t.f4_y1 = entry_y;
  t.f4_z1 = entry_z;
  t.f4_dx1 = entry_dx;
  t.f4_dy1 = entry_dy;
  t.f4_dz1 = entry_dz;
  t.b_downward = true;
  ew.b_is_tracked_energy = false;

  // GET TRACK ORIENTIATION ON SPHERE

  //get_message(bx_message::debug) << "find sphere orientation: steps ..." << dispatch;

  // build the center of hits in time bins of f4_time_limit, referred to as steps

  Float_t step_x[i4_nsteps], step_y[i4_nsteps], step_z[i4_nsteps], step_r[i4_nsteps];
  Float_t step_dx[i4_nsteps], step_dy[i4_nsteps], step_dz[i4_nsteps];
  Int_t step_n[i4_nsteps];
  for (Int_t s = 0; s < i4_nsteps; s++) step_x[s] = step_y[s] = step_z[s] = step_r[s] = step_dx[s] = step_dy[s] = step_dz[s] = step_n[s] = 0;

  Int_t s;
  for (Int_t i = 0; i < nhits_used; i++) {
    s = s_v[i];
    step_x[s] += x_v[i];
    step_y[s] += y_v[i];
    step_z[s] += z_v[i];
    step_n[s]++;
  }

  // normalize steps and compute radius length

  for (Int_t s=0; s<i4_nsteps; s++) {
    if (step_n[s]==0) continue;
    step_x[s] /= step_n[s];
    step_y[s] /= step_n[s];
    step_z[s] /= step_n[s];
    step_r[s] = sqrt( pow(step_x[s],2) + pow(step_y[s],2) + pow(step_z[s],2) );
    //get_message(bx_message::debug) << s << " #; x = " << step_x[s] << ", y = " << step_y[s] << ", z = " << step_z[s] << dispatch;
  }

  // see at which step to stop (the first that is nearer to the entry point than the step before), or the first step without any hits
  
  Int_t stop = i4_nsteps;
  Float_t step_max = 0.;
  Float_t step_width = 0.;

  for (Int_t s=0; s<i4_nsteps; s++) {
    if (step_n[s]==0) {
      stop = s;
      break;
    }
    step_width = pow(step_x[s]-entry_x,2) + pow(step_y[s]-entry_y,2) + pow(step_z[s]-entry_z,2);
    if (step_width > step_max) {
      step_max = step_width;
    }
    else {
      if (s>3) {
        stop = s;
        break;
      }
    }
  }
  
  if (stop<=1) {
    get_message(bx_message::error) << "no hits surrounding the entry point found. exiting." << dispatch;
    return ev;
  }
  //get_message(bx_message::debug) << "last step to be used: " << stop << dispatch;
  
  // finally, compute center of all steps

  Float_t com_x=0, com_y=0, com_z=0, com_r=0;
  {
    for (Int_t s=0; s<stop; s++) {
      com_x += step_x[s];
      com_y += step_y[s];
      com_z += step_z[s];
    }
    com_x /= stop;
    com_y /= stop;
    com_z /= stop;
    com_r = sqrt(com_x*com_x+com_y*com_y+com_z*com_z);
  }

  //get_message(bx_message::debug) << "center: x = " << com_x << ", y = " << com_y << ", z = " << com_z << dispatch;


  // choose the polar coordinates in a way that step_angles are far from phi=0, theta=0,180
  Bool_t turn_X = false, turn_Z = false;
  if (com_z>com_r/sqrt(2.) || com_z<com_r/sqrt(2.)) turn_Z = true;
  if (!turn_Z && com_x>0) turn_X = true;
  if (turn_Z && com_x<0) turn_X = true;

  //get_message(bx_message::debug) << " turn Z? " << turn_Z << ", turn X? " << turn_X << dispatch;
  
  std::vector<Float_t> theta_v, phi_v, r_v;
  theta_v.resize(nhits_used);
  phi_v.resize(nhits_used);
  r_v.resize(nhits_used);
  Float_t com_theta, com_phi;
  Float_t step_theta[i4_nsteps], step_phi[i4_nsteps];
  Float_t step_dtheta[i4_nsteps], step_dphi[i4_nsteps];
  for (Int_t i=0; i<i4_nsteps; i++) { step_dtheta[i]=0; step_dphi[i]=0; }
  Float_t entry_theta, entry_phi;

  if (!turn_Z) {
    if (!turn_X) {
      com_theta   = acos(com_z/com_r);
      com_phi     = atan2(com_y, com_x);
      entry_theta = acos(entry_z/entry_r);
      entry_phi    = atan2(entry_y, entry_x);
      for (Int_t s=0; s<stop; s++) {
        step_theta[s] = acos(step_z[s]/step_r[s]);
	step_phi[s]   = atan2(step_y[s], step_x[s]);
      }
      for (Int_t i=0; i<nhits_used; i++) {
        r_v[i]     = sqrt( x_v[i]*x_v[i] + y_v[i]*y_v[i] + z_v[i]*z_v[i] );
        theta_v[i] = acos(z_v[i]/r_v[i]);
        phi_v[i]   = atan2(y_v[i], x_v[i]);
      }
    }
    else {
      com_theta   = acos(com_z/com_r);
      com_phi     = atan2(-com_y, -com_x);
      entry_theta = acos(entry_z/entry_r);
      entry_phi    = atan2(-entry_y, -entry_x);
      for (Int_t s=0; s<stop; s++) {
        step_theta[s] = acos(step_z[s]/step_r[s]);
        step_phi[s]   = atan2(-step_y[s], -step_x[s]);
      }
      for (Int_t i=0; i<nhits_used; i++) {
        r_v[i]     = sqrt( x_v[i]*x_v[i] + y_v[i]*y_v[i] + z_v[i]*z_v[i] );
        theta_v[i] = acos(z_v[i]/r_v[i]);
        phi_v[i]   = atan2(-y_v[i], -x_v[i]);
      }
    }
  }
  else {
    if (!turn_X) {
      com_theta   = acos(com_y/com_r);
      com_phi     = atan2(-com_z, com_x);
      entry_theta = acos(entry_y/entry_r);
      entry_phi    = atan2(-entry_z, entry_x);
      for (Int_t s=0; s<stop; s++) {
        step_theta[s] = acos(step_y[s]/step_r[s]);
        step_phi[s]   = atan2(-step_z[s], step_x[s]);
      }
      for (Int_t i=0; i<nhits_used; i++) {
        r_v[i]     = sqrt( x_v[i]*x_v[i] + y_v[i]*y_v[i] + z_v[i]*z_v[i] );
        theta_v[i] = acos(y_v[i]/r_v[i]);
        phi_v[i]   = atan2(-z_v[i], x_v[i]);
      }
    }
    else {
      com_theta   = acos(com_y/com_r);
      com_phi     = atan2(com_z, -com_x);
      entry_theta = acos(entry_y/entry_r);
      entry_phi    = atan2(entry_z, -entry_x);
      for (Int_t s=0; s<stop; s++) {
        step_theta[s] = acos(step_y[s]/step_r[s]);
        step_phi[s]   = atan2(step_z[s], -step_x[s]);
      }
      for (Int_t i=0; i<nhits_used; i++) {
        r_v[i]     = sqrt( x_v[i]*x_v[i] + y_v[i]*y_v[i] + z_v[i]*z_v[i] );
        theta_v[i] = acos(y_v[i]/r_v[i]);
        phi_v[i]   = atan2(z_v[i], -x_v[i]);
      }
    }
  }

  //get_message(bx_message::debug) << " Entry  :  theta = " << entry_theta << ", phi = "  << entry_phi << dispatch;
  //get_message(bx_message::debug) << " Center :  theta = " << com_theta << ", phi = "  << com_phi << dispatch;


  // find angular uncertainty of the remaining steps
  
  for (Int_t i=0; i < nhits_used; i++) {
    Int_t s = s_v[i];
    if (s>=stop) continue;
    //get_message(bx_message::error) << theta_v[i] << " " << step_theta[s] << dispatch;
    step_dtheta[s] += (theta_v[i]-step_theta[s])*(theta_v[i]-step_theta[s]);
    step_dphi[s] += (phi_v[i]-step_phi[s])*(phi_v[i]-step_phi[s]);
  } 
  for (Int_t s=0; s<stop; s++) {
    step_dtheta[s] /= step_n[s];
    step_dphi[s]   /= step_n[s];
    //get_message(bx_message::error) << step_n[s] << " " << step_dtheta[s] << " " << step_dphi[s] << dispatch;
  }

  // fit in the theta-phi-plane

  //get_message(bx_message::error) << stop << dispatch;
  //for (Int_t i=0; i<stop; i++) get_message(bx_message::error) << step_theta[i] << " " << step_phi[i] << " " << step_dtheta[i] << " " << step_dphi[i] << dispatch;
  TGraphErrors gtp(stop, step_theta, step_phi, step_dtheta, step_dphi);
  //TGraphErrors gtp(stop, step_theta, step_phi);
  TF1 ftp("ftp", "[0]*x+([1]-[0]*[2])", pi, pi);
  ftp.SetParameter(0, Float_t((com_phi-entry_phi)/(com_theta-entry_theta)));
  ftp.FixParameter(1, Float_t(entry_phi));
  ftp.FixParameter(2, Float_t(entry_theta));
  gtp.Fit(&ftp, "RQ0");
  double m = ftp.GetParameter(0);
  double dm = ftp.GetParError(0);
  //get_message(bx_message::error) << ftp.GetParameter(0) << " " << ftp.GetParError(0) << dispatch;
  double m_theta = 0.05/sqrt(1.+m*m);
  double m_phi   = 0.05*m/sqrt(1.+m*m);
  double dir_theta = entry_theta, dir_phi = entry_phi;
  double dir_dtm = -m/(1.+m*m)*m_theta;
  double dir_dpm = (1./m-m/(1.+m*m))*m_phi;
  if (com_theta > entry_theta) {
    dir_theta += m_theta;
    dir_phi   += m_phi;
  }
  else {
    dir_theta -= m_theta;
    dir_phi   -= m_phi;
    dir_dtm   *= -1.;
    dir_dpm   *= -1.;
  }

  //get_message(bx_message::debug) << " Displ. :  theta = " << dir_theta << ", phi = " << dir_phi << dispatch;

  //get_message(bx_message::debug) << " Fit: m-start = " << (com_phi-entry_phi)/(com_theta-entry_theta) << dispatch;
  //get_message(bx_message::debug) << " Fit: m-stop  = " << m << "+/-" << dm  << dispatch; 

  if (!(m<0 || m>=0)) {
    get_message(bx_message::error) << "Fit did not converge. Exiting." << dispatch;
    return ev;
  }

  // GET IMPACT PARAMETER FROM DECODED HITS AND MEANTIME

  float ndechits_impact[10] = {0,3500,15000,17500,20500,23000,25000,26500,27500,28000};
  float impact_ndechits[10] = {6.85,4.5,3.7,3.5,3,2.5,2,1.5,1,.5};
  float meantime_impact[8]  = {0,70,75,1200,5300,5600,5800,6000};
  float impact_meantime[8]  = {6.85,6.2,5.75,4.5,3.5,3,2.5,0.5};

  int   ndechits = ev->get_laben().get_decoded_nhits();
  float meantime = ev->get_laben().get_cluster(0).get_mean_time();

  int idec = 0;
  while (idec<10) {
    if (ndechits < ndechits_impact[idec]) break;
    idec++;
  }
  float impact_nh  = 1.5;
  float impact_dnh = 2.;
  if (idec<10) {
    impact_nh  = (impact_ndechits[idec]-impact_ndechits[idec-1])/(ndechits_impact[idec]-ndechits_impact[idec-1])*(ndechits-ndechits_impact[idec-1])+impact_ndechits[idec-1];
    impact_dnh = .75;
  }

  //get_message(bx_message::debug) << "Ndechits = " << ndechits << " => impact = " << impact_nh << "+/-" << impact_dnh << dispatch;

  int imt = 0;
  while (imt<8) {
    if (meantime < meantime_impact[imt]) break;
    imt++;
  }
  float impact_mt = 1.5;
  float impact_dmt = 2.;
  if (imt<8) {
    impact_mt  = (impact_meantime[imt]-impact_meantime[imt-1])/(meantime_impact[imt]-meantime_impact[imt-1])*(meantime-meantime_impact[imt-1])+impact_meantime[imt-1];
    impact_dmt = .75;
  }

  //get_message(bx_message::debug) << "meantime = " << meantime << " => impact = " << impact_mt << "+/-" << impact_dmt << dispatch;
  
  double impact_r = ( impact_nh/pow(impact_dnh,2) + impact_mt/pow(impact_dmt,2) ) / ( pow(impact_dnh,-2) + pow(impact_dmt,-2) ); 
  double impact_dr = 1./( pow(impact_dnh,-2) + pow(impact_dmt,-2) );

  //get_message(bx_message::debug) << "=> Impact Parameter = " << impact_r << "+/-" << impact_dr << dispatch;
  
  // GET IMPACT POINT
 
  // get center point in the right distance
  double impact_x = entry_x*impact_r*impact_r/entry_r/entry_r;
  double impact_y = entry_y*impact_r*impact_r/entry_r/entry_r;
  double impact_z = entry_z*impact_r*impact_r/entry_r/entry_r;
  double impact_r1 = sqrt(impact_x*impact_x+impact_y*impact_y+impact_z*impact_z);

  double displ_x, displ_y, displ_z, displ_dx, displ_dy, displ_dz;
  // get small displacement vector from orientation
  if (!turn_Z) {
    if (!turn_X) {
      displ_x = sin(dir_theta)*cos(dir_phi);
      displ_y = sin(dir_theta)*sin(dir_phi);
      displ_z = cos(dir_theta);
      displ_dx = fabs((cos(dir_theta)*dir_dtm*cos(dir_phi)-sin(dir_theta)*sin(dir_phi)*dir_dpm)*dm);
      displ_dy = fabs((cos(dir_theta)*dir_dtm*sin(dir_phi)+sin(dir_theta)*cos(dir_phi)*dir_dpm)*dm);
      displ_dz = fabs(sin(dir_theta)*dir_dtm*dm);
    }
    else {
      displ_x = -sin(dir_theta)*cos(dir_phi);
      displ_y = -sin(dir_theta)*sin(dir_phi);
      displ_z = cos(dir_theta);
      displ_dx = fabs((cos(dir_theta)*dir_dtm*cos(dir_phi)-sin(dir_theta)*sin(dir_phi)*dir_dpm)*dm);
      displ_dy = fabs((cos(dir_theta)*dir_dtm*sin(dir_phi)+sin(dir_theta)*cos(dir_phi)*dir_dpm)*dm);
      displ_dz = fabs(sin(dir_theta)*dir_dtm*dm);
    }
  }
  else {
    if (!turn_X) {
      displ_x = sin(dir_theta)*cos(dir_phi);
      displ_z = -sin(dir_theta)*sin(dir_phi);
      displ_y = cos(dir_theta);
      displ_dx = fabs((cos(dir_theta)*dir_dtm*cos(dir_phi)-sin(dir_theta)*sin(dir_phi)*dir_dpm)*dm);
      displ_dz = fabs((cos(dir_theta)*dir_dtm*sin(dir_phi)+sin(dir_theta)*cos(dir_phi)*dir_dpm)*dm);
      displ_dy = fabs(sin(dir_theta)*dir_dtm*dm);
    }
    else {
      displ_x = -sin(dir_theta)*cos(dir_phi);
      displ_z = sin(dir_theta)*sin(dir_phi);
      displ_y = cos(dir_theta);
      displ_dx = fabs((cos(dir_theta)*dir_dtm*cos(dir_phi)-sin(dir_theta)*sin(dir_phi)*dir_dpm)*dm);
      displ_dz = fabs((cos(dir_theta)*dir_dtm*sin(dir_phi)+sin(dir_theta)*cos(dir_phi)*dir_dpm)*dm);
      displ_dy = fabs(sin(dir_theta)*dir_dtm*dm);
    }
  }

  //get_message(bx_message::debug) << "displ: x = " << displ_x << ", y = " << displ_y << ", z = " << displ_z << dispatch;

  double cos_displ_entry = (displ_x*entry_x+displ_y*entry_y+displ_z*entry_z)/entry_r;

  displ_x *= impact_r1/cos_displ_entry;
  displ_y *= impact_r1/cos_displ_entry;
  displ_z *= impact_r1/cos_displ_entry;

  //get_message(bx_message::debug) << "1st step impact: x = " << impact_x << ", y = " << impact_y << ", z = " << impact_z << dispatch;
  //get_message(bx_message::debug) << "1st step displacement: x = " << displ_x << ", y = " << displ_y << ", z = " << displ_z << dispatch;

  double shift_length = impact_r/entry_r*sqrt(entry_r*entry_r-impact_r*impact_r);
  if (impact_r>entry_r) shift_length = 0;
  double displ_length = sqrt( pow(displ_x-impact_x,2) + pow(displ_y-impact_y,2) + pow(displ_z-impact_z,2) );

  //get_message(bx_message::debug) << "shift length = " << shift_length << ", displ length = " << displ_length << dispatch;

  impact_x += (displ_x-impact_x)*shift_length/displ_length;
  impact_y += (displ_y-impact_y)*shift_length/displ_length;
  impact_z += (displ_z-impact_z)*shift_length/displ_length;

  // approximate the uncertainty
 
//  double f_entry = impact_r1-shift_length;
//  double f_displ = shift_length;
//  double f_impact_e = (2.*impact_r1/impact_r - shift_length/impact_r + impact_r/(impact_r*impact_r-entry_r*entry_r)*shift_length)/entry_r;
//  double f_impact_d = (shift_length/impact_r - impact_r/(impact_r*impact_r-entry_r*entry_r)*shift_length)/impact_r1*cos_displ_entry; 
//  double impact_dx = sqrt( pow(f_entry*sqrt(entry_dx*entry_dx-0.3),2) + pow(f_displ*displ_dx,2) + pow((f_impact_e*entry_x+f_impact_d*displ_x)*impact_dr,2) + 0.4); 
//  double impact_dy = sqrt( pow(f_entry*sqrt(entry_dy*entry_dy-0.3),2) + pow(f_displ*displ_dy,2) + pow((f_impact_e*entry_y+f_impact_d*displ_y)*impact_dr,2) + 0.4);
//  double impact_dz = sqrt( pow(f_entry*sqrt(entry_dz*entry_dz-0.3),2) + pow(f_displ*displ_dz,2) + pow((f_impact_e*entry_z+f_impact_d*displ_z)*impact_dr,2) + 0.4);  
  double impact_dx = sqrt( impact_dr*impact_dr + 0.5 ); 
  double impact_dy = sqrt( impact_dr*impact_dr + 0.5 );
  double impact_dz = sqrt( impact_dr*impact_dr + 0.5 );

  //get_message(bx_message::debug) << "Impact Point: x = " << impact_x << "+/-" << impact_dx << ", y = " << impact_y << "+/-" << impact_dy << ", z = " << impact_z << "+/-" << impact_dz << dispatch;
   float dex = sqrt(pow(entry_x-impact_x,2)+pow(entry_y-impact_y,2)+pow(entry_z-impact_z,2));
   float impact_rad = sqrt(pow(impact_x,2)+pow(impact_y,2)+pow(impact_z,2));

  t.f4_t2 = t.f4_t1+dex/0.3;   //er.get_cluster(0).get_start_time() - er.get_trigger_rawt();  
  t.f4_x2 = impact_x;
  t.f4_y2 = impact_y;
  t.f4_z2 = impact_z;
  t.f4_dx2 = impact_dx;
  t.f4_dy2 = impact_dy;
  t.f4_dz2 = impact_dz;
  t.f4_phi = atan2(entry_y-impact_y, entry_x-impact_x);
  while (t.f4_phi<0) t.f4_phi += 2*constants::number::pi;
  t.f4_theta = acos((entry_z-impact_z)/dex);
  t.f4_impact = impact_rad;
  t.b_downward = (entry_z>impact_z);
  
  // write normalized number of decoded hits into track variable
  if ( ev->get_laben().get_nclusters()>0 )
    t.f4_labennormhits = ev->get_laben().get_decoded_nhits() / (ev->get_laben().get_n_live_pmts() - ev->get_laben().get_invalid_pmts()) * 2000.;
  
  ew.b_is_tracked_energy = true;


  ev->get_laben().mark_stage (bx_base_event::tracked);

  i4_tracked_events++;
  return ev;  
}


void bx_laben_energy_tracker::end () {
  if (i4_laben_events > 0) get_message (bx_message::info) << "run summary: " << i4_tracked_events << " of " << i4_laben_events << " events were reconstructed." << dispatch;
}	

/*
 * $Log: bx_laben_energy_tracker.cc,v $
 * Revision 1.7  2014/12/11 21:27:12  wurm
 * added getters and variables to BxEvent that concern C-11 likelihood tagging
 * * for tracks: new distance-to- and projection-on-track getters, dE/dx calculation
 * * counting of neutron-like clusters in tt128, new flag of neutron clusters
 *
 * Revision 1.6  2010/08/03 15:57:34  wurm
 * 1) introduced theta, phi, impact variables for fitted tracks
 * 2) set fixed uncertainties for global tracking
 * 3) condition for execution of laben tracking changed to MTF OR MCF
 * 4) lowered verbosity of all modules
 *
 * Revision 1.5  2010-07-30 17:06:59  wurm
 * require EP and XP for is_tracked variable to be set true
 *
 * Revision 1.4  2010-07-01 12:06:37  wurm
 * added computation of track angles, impact
 *
 * Revision 1.3  2010-05-26 15:11:14  ddangelo
 * last commits on c12 brought to c13:
 * 1) error meassages restored
 * 2) use of root types
 * 3) bug fix for initializion of vectors in exit point calculation
 *
 * Revision 1.11.2.2  2010-05-19 12:27:49  wurm
 * changed data types for 64b
 *
 * Revision 1.11.2.1  2010-05-19 09:26:01  ddangelo
 * error messages restored
 *
 * Revision 1.11  2009-10-23 09:28:46  wurm
 * commented debug messages
 *
 * Revision 1.10  2009-07-30 11:25:18  wurm
 * removed the removal
 *
 * Revision 1.9  2009-07-30 09:41:17  wurm
 * removed bug for filling is_tracked variable
 *
 * Revision 1.8  2009-01-24 12:03:43  razeto
 * Amin commit: call fit with N to avoid fit drawing
 *
 * Revision 1.7  2008-12-15 08:31:10  wurm
 * mended non-converging fits
 *
 * Revision 1.6  2008-09-02 07:28:11  wurm
 * new uncertainties
 *
 * Revision 1.5  2008-07-21 15:00:52  wurm
 * added impact point
 *
 * Revision 1.3  2008-02-22 12:02:28  ddangelo
 * required muon event
 * required min 1 cluster laben
 * event stage marking
 *
 * Revision 1.2  2008-02-21 19:06:41  ddangelo
 * computes entry point (preliminary) with errors
 * to be tested
 *
 * Revision 1.1  2008-02-21 17:17:08  ddangelo
 * added. empty.
 *
 * 
 */

