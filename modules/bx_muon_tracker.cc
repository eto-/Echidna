/* BOREXINO Reconstruction program
 *
 * Author: Michael Wurm <mwurm@ph.tum.de>, Davide D'Angelo <davide.dangelo@mi.infn.it> 
 * Maintainer: Davide D'Angelo <davide.dangelo@mi.infn.it>
 *
 * $Id: bx_muon_tracker.cc,v 1.29 2014/12/11 21:27:12 wurm Exp $
 *
 * Implemenentation of bx_muon_tracker
 *
*/

#include "bx_muon_tracker.hh"
#include "bx_echidna_event.hh"
#include "bx_laben_event.hh"
#include "bx_dbi.hh"
#include "db_profile.hh"
#include "db_channel.hh"
#include "db_calib.hh"
#include "constants.hh"
#include "barn_interface.hh"
#include "bx_track.hh"

static float pi = constants::number::pi;

bx_muon_tracker::bx_muon_tracker (): bx_base_module("bx_muon_tracker", bx_base_module::main_loop) {
  require_event_stage (bx_detector::muon, bx_base_event::clustered);
  require_trigger_type (bx_trigger_event::neutrino);
  require_trigger_type (bx_trigger_event::muon);
}

void bx_muon_tracker::begin () {
  // statistics ctrs
  i4_muon_events = 0;
  i4_rec_events  = 0;

  i4_enable_histos        = get_parameter ("enable_histos")       .get_int  ();
  f4_tau                  = get_parameter ("tau")                 .get_float();
  f4_hit_charge_threshold = get_parameter ("hit_charge_threshold").get_float();
    
  f4_entry_tau            = get_parameter ("entry_tau")           .get_float();
  f4_entry_mean           = get_parameter ("entry_mean")          .get_float();
  f4_entry_sigma          = get_parameter ("entry_sigma")         .get_float();
  f4_exit_incline         = get_parameter ("exit_incline")        .get_float();
  f4_exit_dt_min          = get_parameter ("exit_dt_min")         .get_float();
  f4_exit_sigma           = get_parameter ("exit_sigma")          .get_float();


  if (i4_enable_histos) {
    tracks_costheta =  new TH1F ("tracks_costheta", "bx_muon_tracker_cos(theta)", 30 , -1.,   1. );
    tracks_costheta -> SetFillColor(4);
    tracks_costheta -> SetXTitle("cos(theta)");

    tracks_phi      =  new TH1F ("tracks_phi"     , "bx_muon_tracker_phi"       , 30 ,  0., 360. );
    tracks_phi	-> SetFillColor(4);
    tracks_phi	-> SetXTitle("phi [degree]");

    tracks_distance = new TH1F ("tracks_distance", "bx_muon_tracker_distance"  , 30 ,  0.,   9. );
    tracks_distance -> SetFillColor(4);
    tracks_distance -> SetXTitle("pedal point distance to center [m]");

    barn_interface::get ()->store (barn_interface::file, tracks_costheta, this);
    barn_interface::get ()->store (barn_interface::file, tracks_phi     , this);
    barn_interface::get ()->store (barn_interface::file, tracks_distance, this);
  }
}

bx_echidna_event* bx_muon_tracker::doit (bx_echidna_event *ev) {

  i4_muon_events++;
  const bx_muon_event& er = ev->get_muon();
  bx_muon_tracked_event &ew = dynamic_cast<bx_muon_tracked_event&>(ev->get_muon());

  //get_message (bx_message::debug) << "event " << ev->get_event_number () << "; " << er.get_nclusters() << " cluster(s)" << dispatch;

  // decide: 0 clusters found: stop reco
  //         1 cluster found:  try to split this cluster in two clusters if it is on sphere
  //	   2 clusters found: fine
  //	  >2 clusters found: select the right ones

  /*for (int32_t ic = 0; ic < er.get_nclusters(); ic++) {
    get_message (bx_message::debug) << "# " << ic << " ID " << er.get_cluster(ic).get_id() << " sss " << er.get_cluster(ic).is_sss() << " t " << er.get_cluster(ic).get_start_time() << " x " << er.get_cluster(ic).get_x() << " y " << er.get_cluster(ic).get_y() << " z " << er.get_cluster(ic).get_z() << " Q " << er.get_cluster(ic).get_charge() << dispatch;
  }*/
  // ----------------- check if coordinates of clusters are sensible

  int32_t noc = ev->get_muon().get_nclusters();
  std::vector<int32_t> good_v; // indicates if the cluster should be considered for tracking
  int32_t nogc = noc; // number of good clusters is number of clusters
  for (int32_t ic = 0; ic < noc ; ic++) {
    if ( er.get_cluster(ic).get_x()>-20 && er.get_cluster(ic).get_x()<20 ) { good_v.push_back(1); }
    else {
      good_v.push_back(0);
      nogc--;
    }
  }

  if (nogc<1) {
    //get_message (bx_message::debug) << "no valid cluster from bx_muon_findcluster passed. skip tracking" << dispatch;
    return ev;
  }


  //
  // ----------------- first step: select the entry cluster
  //

  double starttime_ID = ev->get_muon().get_start_time_sss();
  bool noID = true;
  
  // decide which start time to use; if there is info from ID, it is preferred
  
  if (ev->get_trigger().get_trigger_type()==2) {
    //get_message (bx_message::debug) <<  "trigger type 2 >>> track reconstruction is not reliable! start time of OD-SSS " << starttime_ID << " was chosen" << dispatch;
  }
  else {
    if (ev->get_laben().get_nclusters()==0) {
      //get_message (bx_message::debug) <<  "trigger type 1, but no cluster in ID! >>> start time of OD-SSS " << starttime_ID << " was chosen" << dispatch;
    }
    else {
      starttime_ID = ev->get_laben().get_cluster(0).get_start_time() - ev->get_laben().get_trigger_rawt();
       
      //get_message (bx_message::debug) <<  "trigger type 1 and laben cluster present >>> starttime ID is " << starttime_ID  << dispatch;
      noID = false;
    }
  }

  std::vector<float> t_cl;
 
  //get_message (bx_message::debug) << "start_time_sss is " << er.get_start_time_sss() << dispatch;

  // ugly fix to take into account the shift in relative laben/muon gate positions due to the installation of new muon electronics
  if (ev->get_run_number()<17502) f4_entry_mean = +548.;
  else f4_entry_mean = -345.;

  for (int32_t ic=0; ic<noc; ic++) {
    // fill the arrays
    t_cl.push_back( er.get_cluster(ic).get_start_time()+er.get_start_time_sss() );
    // test if time is aligned to first ID cluster (removed for the moment as the relative time delay between ID and OD is changing considerably over time, a less stringent cut of 1mus is introduced)
    //if (t_cl[ic]-starttime_ID<-200.+f4_entry_mean || t_cl[ic]-starttime_ID>f4_entry_mean+200.) {
    if (fabs(t_cl[ic]-starttime_ID)>1000.) {
      t_cl[ic] = 0;
      nogc--;
    }
  }

 // get_message (bx_message::debug) << "after time check: nogc=" << nogc << dispatch;

  if (nogc<1) return ev;

  std::vector<double> entry_weight_cl;
  if (noID) f4_entry_mean = 0.;
  float starttime_sc = 0;

  // compute the weight of each cluster: charge times a steep time function times a gaussian of the expected value relative to ID

  for (int32_t ic=0; ic<noc; ic++) {
    if (good_v[ic]==1) {
      if (starttime_sc == 0) starttime_sc = t_cl[ic];
      entry_weight_cl.push_back( er.get_cluster(ic).get_charge()*exp(-fabs(t_cl[ic]-starttime_sc)/f4_entry_tau) ); 
      // * exp(-pow(t_cl[ic]-starttime_ID-f4_entry_mean,2)/(2.*pow(f4_entry_sigma,2))) );
      // weighting term for ID-OD time correlation removed because of the variation in ID-OD time alignment
    }
    else entry_weight_cl.push_back(0);
    //get_message (bx_message::debug) << "# " << ic << " ID " << er.get_cluster(ic).get_id() << " t=" << t_cl[ic] << " Q=" << er.get_cluster(ic).get_charge() << " WEn=" << entry_weight_cl[ic] << dispatch;
  }

  int32_t entry_id=0;
  double entry_weight=0;
  
  // find "heaviest" cluster
  for (int32_t ic=0; ic<noc; ic++) {
    if (entry_weight_cl[ic]>entry_weight) {
      entry_weight = entry_weight_cl[ic];
      entry_id = er.get_cluster(ic).get_id();
    }
  }

  if (entry_weight<=0) {
    //get_message (bx_message::debug) << "entry point corrupted: Weight is " << entry_weight << dispatch;
    return ev;
  }

  // entry point was found!!

  //get_message (bx_message::debug) << "entry point: # " << entry_id-1 << " ID " << entry_id  << " WEn=" << entry_weight << dispatch;
 
  coo ep;

  ep.x = er.get_cluster(entry_id-1).get_x();
  ep.y = er.get_cluster(entry_id-1).get_y();
  ep.z = er.get_cluster(entry_id-1).get_z();
  ep.dx = 0;
  ep.dy = 0;
  ep.dz = 0;
  ep.time = t_cl[entry_id-1];

  float weightsum = 0;
  // entry uncertainty
  int32_t entryhits = 0;
  double dt_first = -1e9;
  float maxweight = 0;
  for (int32_t ihit=0; ihit<er.get_clustered_nhits(); ihit++) {
    if ( er.get_clustered_hit(ihit).get_affiliation() != entry_id ) continue;
    if (dt_first == -1e9) dt_first = er.get_clustered_hit(ihit).get_time();
    double weight = er.get_clustered_hit(ihit).get_charge()*exp(-(er.get_clustered_hit(ihit).get_time()-dt_first)/f4_tau);
    ep.dx += pow(er.get_clustered_hit(ihit).get_decoded_hit().get_db_channel()->get_x()-ep.x,2)*weight;
    ep.dy += pow(er.get_clustered_hit(ihit).get_decoded_hit().get_db_channel()->get_y()-ep.y,2)*weight;
    ep.dz += pow(er.get_clustered_hit(ihit).get_decoded_hit().get_db_channel()->get_z()-ep.z,2)*weight;
    weightsum += weight;
    entryhits++;
    if (weight>maxweight) maxweight = weight;
  }
  ep.dx = sqrt(ep.dx/weightsum/(weightsum/maxweight)+0.03);
  ep.dy = sqrt(ep.dy/weightsum/(weightsum/maxweight)+0.03);
  ep.dz = sqrt(ep.dz/weightsum/(weightsum/maxweight)+0.03);

  //get_message (bx_message::debug) << "uncertainties computed" << dispatch;
    
  // fill first part of event

  if (entry_id==0) return ev;

  bx_track_by_points& t = ew.get_track();
  t.f4_t1 = ep.time;
  t.f4_x1 = ep.x;
  t.f4_y1 = ep.y;
  t.f4_z1 = ep.z;
  t.f4_dx1 = ep.dx;
  t.f4_dy1 = ep.dy;
  t.f4_dz1 = ep.dz;
  t.b_downward = true;
  ew.b_is_tracked = false;
  //get_message (bx_message::debug) << "entry point written to track" << dispatch;

  // set goodness of entry point to zero, count if there is still a cluster left
  good_v[entry_id-1]=0;
  nogc--;
  if (nogc<1) {
    //get_message (bx_message::debug) << "no 2nd cluster present; only entry point is written!" << dispatch;
    return ev;
  }


  //
  // --------------------- 2nd step: find the exit point
  //

  //get_message(bx_message::debug) << "Search for exit point:" << dispatch;

  std::vector<double> lambda_v;
  std::vector<double> impact_v;
  for (int32_t ic=0; ic<noc; ic++) {
    if (good_v[ic]==1) {
      lambda_v.push_back( - ( (er.get_cluster(ic).get_x()-ep.x)*ep.x + (er.get_cluster(ic).get_y()-ep.y)*ep.y + (er.get_cluster(ic).get_z()-ep.z)*ep.z ) / ( pow(er.get_cluster(ic).get_x()-ep.x,2) + pow(er.get_cluster(ic).get_y()-ep.y,2) + pow(er.get_cluster(ic).get_z()-ep.z,2) ) );
      impact_v.push_back( sqrt( pow(ep.x+(er.get_cluster(ic).get_x()-ep.x)*lambda_v[ic],2) + pow(ep.y+(er.get_cluster(ic).get_y()-ep.y)*lambda_v[ic],2) + pow(ep.z+(er.get_cluster(ic).get_z()-ep.z)*lambda_v[ic],2)) );
      if (!(lambda_v[ic]>0 && lambda_v[ic]<1)) {
        good_v[ic]=0;             // test if the impact point is between the two chosen points and if impact point is inside sphere
	nogc--;
      }
      else {
        if ( ev->get_trigger().get_trigger_type()==1 && impact_v[ic]>7.1 ) {
 	  good_v[ic]=0;
	  nogc--;
	}
	if ( ev->get_trigger().get_trigger_type()==2 && impact_v[ic]<6.5 ) {
          good_v[ic]=0;
          nogc--;
	}
      }
    }
    else {
      lambda_v.push_back(0);
      impact_v.push_back(0);
    }
    //get_message (bx_message::debug) << "# " << ic << " ID " << ic+1 << " lambda = " << lambda_v[ic] << ", impact = " << impact_v[ic] << " => good = " << good_v[ic] << dispatch;
  }
 
  if (nogc<1) {
    //get_message (bx_message::debug) << "no track crossing the ID found! only entry point was written!" << dispatch;
    return ev;
  }

  // check if time alignment is okay
  std::vector<double> exit_weight_cl;
  
  for (int32_t ic=0; ic<noc; ic++) {
    exit_weight_cl.push_back(0);
    //get_message (bx_message::debug) << "ic=" << ic << "exit_weight_cl=" << exit_weight_cl[ic] << dispatch;
    float optimal_dt = 0, dr = 0;
    //get_message (bx_message::debug) << "ic=" << ic << "optimal_dt=" << optimal_dt << ", dr =" << dr << dispatch;
    if (good_v[ic]==1) {
      //get_message (bx_message::debug) << "entered good==1 loop" << dispatch;
      dr = sqrt( pow(er.get_cluster(ic).get_x()-ep.x,2) + pow(er.get_cluster(ic).get_y()-ep.y,2) + pow(er.get_cluster(ic).get_z()-ep.z,2) );
      optimal_dt = f4_exit_incline*dr; 
      //get_message (bx_message::debug) << "optimal_dt=" << optimal_dt << dispatch;
      if (optimal_dt>f4_exit_dt_min) {
        exit_weight_cl[ic] = er.get_cluster(ic).get_charge() * exp(-pow(t_cl[ic]-ep.time-optimal_dt,2)/(2.*pow(f4_exit_sigma*optimal_dt,2)));
        //get_message (bx_message::debug) << "exit_weight_cl=" << exit_weight_cl[ic] << dispatch;
      }
      else {
        good_v[ic] = 0;
	nogc--;
        //get_message (bx_message::debug) << "exit_weight_cl=" << exit_weight_cl[ic] << dispatch;
      }
    }
    //get_message (bx_message::debug) << "# " << ic << " ID " << ic+1 <<" optimal_dt " << optimal_dt << " -> real_dt " << t_cl[ic]-entry_t << ", good=" << good_v[ic] << " => weight " << exit_weight_cl[ic] << dispatch;   
  }
  
  if (nogc<1) {
    //get_message (bx_message::debug) << "bad time alignment of all exit point candidates! only entry point was written!" << dispatch;
    return ev;
  }

  // find the heaviest exit cluster
  int32_t exit_id=0;
  double exit_weight=0;
  for (int32_t ic=0; ic<noc; ic++) {
    if (exit_weight_cl[ic]>exit_weight) {
      exit_weight = exit_weight_cl[ic];
      exit_id = er.get_cluster(ic).get_id();
    }
  }	

  if (exit_id<=0) {
    //get_message (bx_message::debug) << "exit point corrupted! weight is " << exit_weight << "; only entry point was written!" << dispatch;
    return ev;
  }

  //get_message (bx_message::debug) << "WEIGHT RESULTS" << dispatch;
  /*for (int32_t ic=0; ic<noc; ic++) {
    get_message (bx_message::debug) << "# " << ic << " ID " << ic+1 << " good=" << good_v[ic] << ", EnW=" << entry_weight_cl[ic] << ", ExW=" << exit_weight_cl[ic] << dispatch;
  }*/
 
  coo xp;
  xp.x = er.get_cluster(exit_id-1).get_x();
  xp.y = er.get_cluster(exit_id-1).get_y();
  xp.z = er.get_cluster(exit_id-1).get_z();
  xp.dx = 0;
  xp.dy = 0;
  xp.dz = 0;
  xp.time = t_cl[exit_id-1];

  weightsum = 0;
  maxweight = 0;
  dt_first = -1e9;
  int32_t exithits = 0;
  for (int32_t ihit=0; ihit<er.get_clustered_nhits(); ihit++) {
    if ( er.get_clustered_hit(ihit).get_affiliation() != exit_id ) continue;
    if (dt_first==-1e9) dt_first = er.get_clustered_hit(ihit).get_time();
    float weight = er.get_clustered_hit(ihit).get_charge()*exp(-(er.get_clustered_hit(ihit).get_time()-dt_first)/f4_tau);
    xp.dx += pow(er.get_clustered_hit(ihit).get_decoded_hit().get_db_channel()->get_x()-xp.x,2)*weight;
    xp.dy += pow(er.get_clustered_hit(ihit).get_decoded_hit().get_db_channel()->get_y()-xp.y,2)*weight;
    xp.dz += pow(er.get_clustered_hit(ihit).get_decoded_hit().get_db_channel()->get_z()-xp.z,2)*weight;
    weightsum += weight;
    exithits++;
    if (maxweight<weight) maxweight = weight;
  }
  xp.dx = sqrt(xp.dx/weightsum/(weightsum/maxweight)+0.33);
  xp.dy = sqrt(xp.dy/weightsum/(weightsum/maxweight)+0.33);
  xp.dz = sqrt(xp.dz/weightsum/(weightsum/maxweight)+0.33);

  // construct auxiliary parameters
  ep.rad = m_get_radius(ep);
  float dex = sqrt(pow(ep.x-xp.x,2)+pow(ep.y-xp.y,2)+pow(ep.z-xp.z,2));
  float cosalpha = (-ep.x*(xp.x-ep.x)-ep.y*(xp.y-ep.y)-ep.z*(xp.z-ep.z))/(dex*ep.rad);
  coo pedal;
  pedal.x = ep.x + (xp.x-ep.x)*ep.rad*cosalpha/dex;
  pedal.y = ep.y + (xp.y-ep.y)*ep.rad*cosalpha/dex;
  pedal.z = ep.z + (xp.z-ep.z)*ep.rad*cosalpha/dex;
  pedal.rad = m_get_radius(pedal);

  // filling second part of event variables
  t.f4_t2 = xp.time;
  t.f4_x2 = xp.x;
  t.f4_y2 = xp.y;
  t.f4_z2 = xp.z;
  t.f4_dx2 = xp.dx;
  t.f4_dy2 = xp.dy;
  t.f4_dz2 = xp.dz;
  t.f4_phi = m_limit_angle(atan2(ep.y-xp.y, ep.x-xp.x),0,2*pi);
  t.f4_theta = acos((ep.z-xp.z)/dex);
  t.f4_impact = pedal.rad;
  t.b_downward = (ep.x>xp.x);
  ew.b_is_tracked = true;

  // write normalized number of decoded hits into track variable
  if ( ev->get_laben().get_nclusters()>0 )
    t.f4_labennormhits = ev->get_laben().get_decoded_nhits() / (ev->get_laben().get_n_live_pmts() - ev->get_laben().get_invalid_pmts()) * 2000.;
 
  if ( !(xp.x<20 && xp.x>-20) ) {
    //get_message (bx_message::debug) << "exit point spacial information is corrupted! entry point was written!" << dispatch;  
    return ev;
  }


  //if (!t.is_upward()) get_message (bx_message::debug) << "upgoing muon!" << dispatch;
  //get_message (bx_message::debug) << "clusters ID " << entry_id << " and ID " << exit_id << " were chosen." << dispatch;
  //get_message (bx_message::debug) << "both entry and exit point were written!" << dispatch;
  //get_message (bx_message::debug) << "entry point : t " << entry_t << ", x " << entry_x << "+/-" << entry_dx << ", y " << entry_y << "+/-" << entry_dy << ", z " << entry_z << "+/-" << entry_dz << dispatch;
  //get_message (bx_message::debug) << "exit point  : t "  << exit_t << ", x " << exit_x << "+/-" << exit_dx  << ", y " << exit_y << "+/-" << exit_dy << ", z " << exit_z << "+/-" << exit_dz << dispatch;

  // fill histograms
  if (i4_enable_histos==1) {
    tracks_costheta -> Fill( cos(t.get_theta()) );
    tracks_phi	    -> Fill( t.get_phi()*180./pi );
    tracks_distance -> Fill( t.get_impact() );
  }

  //get_message (bx_message::debug) << "pedal point : t " << t.get_pedal_t() << ", x " << t.get_pedal_x() << ", y " << t.get_pedal_y() << ", z " << t.get_pedal_z() << dispatch;	
  //get_message (bx_message::debug) << "muon track : theta " << t.get_theta()*180./pi << "°, phi " << t.get_phi()*180./pi << "°" << dispatch;

  ev->get_muon().mark_stage (bx_base_event::tracked); 

  i4_rec_events++;
  return ev;  
}


void bx_muon_tracker::end () {
  if (i4_muon_events > 0) get_message (bx_message::info) << "run summary: " << i4_rec_events << " of " << i4_muon_events << " events were reconstructed." << dispatch;
}

coo bx_muon_tracker::m_rotate_x(coo input, float phi) {
	coo output;
	output.x = input.x;
	output.y = input.y*cos(phi)-input.z*sin(phi);
	output.z = input.y*sin(phi)+input.z*cos(phi);
	output.dx = input.dx;
	output.dy = input.dy*cos(phi)-input.dz*sin(phi);
	output.dz = input.dy*sin(phi)+input.dz*cos(phi);	
	output.rad = input.rad;
	output.time = input.time;	
	return output;	
}

coo bx_muon_tracker::m_rotate_y(coo input, float phi) {
	coo output;
	output.x = input.x*cos(phi)-input.z*sin(phi);
	output.y = input.y;
	output.z = input.x*sin(phi)+input.z*cos(phi);			
	output.dx = input.dx*cos(phi)-input.dz*sin(phi);
	output.dy = input.dy;
	output.dz = input.dx*sin(phi)+input.dz*cos(phi);		
	output.rad = input.rad;
	output.time = input.time;	
	return output;	
}

coo bx_muon_tracker::m_rotate_z(coo input, float phi) {
	coo output;
	output.x = input.x*cos(phi)-input.y*sin(phi);
	output.y = input.x*sin(phi)+input.y*cos(phi);
	output.z = input.z;	
	output.dx = input.dx*cos(phi)-input.dy*sin(phi);
	output.dy = input.dx*sin(phi)+input.dy*cos(phi);
	output.dz = input.dz;	
	output.rad = input.rad;
	output.time = input.time;	
	return output;	
}
							   
coo bx_muon_tracker::m_construct_xyz(float rad, float phi, float theta, float dphi, float dtheta) {
	coo output;
	output.x = rad*cos(phi)*sin(theta);
	output.y = rad*sin(phi)*sin(theta);
	output.z = rad*cos(theta);
	output.rad = rad;
	output.dx = rad*sqrt(pow(sin(phi)*sin(theta)*dphi,2)+pow(cos(phi)*cos(theta)*dtheta,2));
	output.dy = rad*sqrt(pow(cos(phi)*sin(theta)*dphi,2)+pow(sin(phi)*cos(theta)*dtheta,2));
	output.dz = rad*sin(theta)*dtheta;
	output.phi = phi;
	output.theta = theta;
	output.dphi = dphi;
	output.dtheta = dtheta;
	output.time = 0;
	return output;
}

coo bx_muon_tracker::m_construct_rpt(float x, float y, float z) {
	coo output;
	output.x = x;
	output.y = y;
	output.z = z;
	output.rad = sqrt(x*x+y*y+z*z);
	output.dx = 0;
	output.dy = 0;
	output.dz = 0;
	output.phi = atan2(y, x);
	output.theta = acos(z/output.rad);
	output.time = 0;
	return output;
}

float bx_muon_tracker::m_limit_angle(float phi, float min_angle, float max_angle) {
	float corrector = max_angle-min_angle;
	while (phi<min_angle) phi += corrector;
	while (phi>max_angle) phi -= corrector;
	return phi;
}

float bx_muon_tracker::m_get_radius(coo input) {
	return sqrt(input.x*input.x+input.y*input.y+input.z*input.z);
}

float bx_muon_tracker::m_get_radius(float x, float y, float z=0) {
	return sqrt(x*x+y*y+z*z);
}


coo bx_muon_tracker::m_normalize_vector(coo input) {
	float factor = 1./input.rad;
	coo output;
	output.x = input.x*factor;
	output.y = input.y*factor;
	output.z = input.z*factor;
	output.rad = input.rad*factor;
	output.dx = input.dx*factor;
	output.dy = input.dy*factor;
	output.dz = input.dz*factor;
	output.theta = input.theta;
	output.phi = input.phi;
	output.time = input.time;
	return output;	
}




/*
 * $Log: bx_muon_tracker.cc,v $
 * Revision 1.29  2014/12/11 21:27:12  wurm
 * added getters and variables to BxEvent that concern C-11 likelihood tagging
 * * for tracks: new distance-to- and projection-on-track getters, dE/dx calculation
 * * counting of neutron-like clusters in tt128, new flag of neutron clusters
 *
 * Revision 1.28  2012/10/25 10:11:55  wurm
 * removed/softened the relative timing weight/cut between OD entry point and 1st ID cluster
 *
 * Revision 1.27  2012-05-04 12:46:53  wurm
 * updated the hardcoded relative delays of ID and OD
 *
 * Revision 1.26  2012-02-16 08:49:58  wurm
 * same
 *
 * Revision 1.25  2012-02-16 08:39:31  wurm
 * ugly fix to take into account the shift in relative laben/muon gate positions due to the installation of new muon electronic
 * s
 *
 * Revision 1.24  2012-01-22 18:41:34  wurm
 * fixed hard-coded timing control for entry point
 *
 * Revision 1.23  2010-07-30 17:07:00  wurm
 * require EP and XP for is_tracked variable to be set true
 *
 * Revision 1.22  2010-06-28 10:44:03  wurm
 * new laben tof tracking, modifications to muon tracking (introducing phi,theta,impact to the event) and modification of global tracker to prefer laben tof over laben energy (and conditions)
 *
 * Revision 1.21  2009-10-26 11:19:37  ddangelo
 * trgtype 'muon_tct' (or 'muon_totc') renamed as 'neutron' after real use
 * trgtype 'muon_mtb' renamed as 'muon' as duplicity with the above one disappeared and 'mtb' is associated to btb_inputs flag
 *
 * Revision 1.20  2009-10-23 09:28:28  wurm
 * commented debug messages
 *
 * Revision 1.19  2009-07-30 11:25:03  wurm
 * removed the removal
 *
 * Revision 1.18  2009-07-30 09:40:59  wurm
 * removed bug for filling is_tracked variable
 *
 * Revision 1.17  2008-09-02 07:27:50  wurm
 * new uncertainties
 *
 * Revision 1.16  2008-08-05 15:31:19  ddangelo
 * tmp commit after conflict, to be inspected
 *
 * Revision 1.15  2008-07-21 15:00:14  wurm
 * added check of time dependence of clusters
 *
 * Revision 1.13  2008-02-21 15:59:18  ddangelo
 * adapted to new track indirection level
 * now marking tracked stage
 *
 * Revision 1.12  2008-02-02 15:31:35  wurm
 *
 *
 * added one-hit-cluster
 *
 * Revision 1.11  2008-01-30 17:55:44  wurm
 *
 * added statistical errors to entry and exit point, deactivitated merging. tracker is now checking for time-order of clusters
 *
 * Revision 1.10  2008-01-06 13:12:10  ddangelo
 * lowered verbosity
 *
 * Revision 1.9  2007-11-28 19:29:08  wurm
 *
 * vetoed too early clusters on ground
 *
 * Revision 1.8  2007-11-28 16:03:47  wurm
 *
 *
 * lowered verbosity
 *
 * Revision 1.7  2007-11-27 19:01:09  wurm
 *
 *
 * removed some bugs
 *
 * Revision 1.6  2007-11-26 18:49:59  ddangelo
 * cleaned up and commented
 *
 * Revision 1.5  2007-11-26 14:06:27  ddangelo
 * completely redone
 *
 * Revision 1.4  2007-11-15 12:48:12  wurm
 * lowered verbosity
 *
 * Revision 1.3  2007-11-14 19:26:55  ddangelo
 * quiter
 *
 * Revision 1.2  2007-11-14 19:03:46  ddangelo
 * added michi's version.
 * writing to event
 * much more job to come
 *
 * Revision 1.1  2007-05-03 15:48:30  ddangelo
 * just added. empty skeleton
 *
 */

