/* BOREXINO Reconstruction program
 *
 * Author: Davide D'Angelo <davide.dngelo@mi.infn.it>, Michael Wurm <mwurm@ph.tum.de>
 * Maintainer: Davide D'Angelo <davide.dangelo@mi.infn.it>
 *
 * $Id: bx_muon_findcluster.cc,v 1.35 2009/11/08 15:35:06 wurm Exp $
 *
 * Implemenentation of bx_muon_findcluster
 *
*/

#include <algorithm>
#include "bx_muon_findcluster.hh"
#include "bx_echidna_event.hh"
#include "bx_dbi.hh"
#include "db_profile.hh"
#include "db_channel.hh"
#include "db_calib.hh"
#include "constants.hh"
#include "barn_interface.hh"

bx_muon_findcluster::bx_muon_findcluster (): bx_base_module("bx_muon_findcluster", bx_base_module::main_loop) {
  require_event_stage (bx_detector::muon, bx_base_event::decoded);
  require_trigger_type (bx_trigger_event::neutrino);
  require_trigger_type (bx_trigger_event::muon);
}
      
void bx_muon_findcluster::begin () {
  const db_profile& profile_info = bx_dbi::get ()->get_profile ();
  float time_buffer_ns = profile_info.muon_tdc_range() * constants::muon::tdc::ns_per_clock;

  i4_clustered_events = 0; // just a statistics counter
  i4_start_threshold  = get_parameter ("start_threshold").get_int ();
  i4_bin_width        = get_parameter ("bin_width")      .get_int ();
  i4_n_bins           = int32_t(time_buffer_ns/i4_bin_width);
  i4_count_threshold  = get_parameter ("count_threshold").get_int ();
  i4_ripple_count     = get_parameter ("ripple_count")   .get_int ();
  i4_enable_histos    = get_parameter ("enable_histos")  .get_int(); 
 
  f4_mincharge            = get_parameter ("mincharge")           .get_float();
  f4_maxdt                = get_parameter ("maxdt")               .get_float();
  f4_maxdr_sss            = get_parameter ("maxdr_sss")           .get_float();
  f4_maxdr_floor          = get_parameter ("maxdr_floor")         .get_float();
  f4_tau                  = get_parameter ("tau")                 .get_float();
  i4_max_hits_sss         = get_parameter ("max_hits_sss")        .get_int();
  i4_max_hits_floor       = get_parameter ("max_hits_floor")      .get_int();
  f4_hit_charge_threshold = get_parameter ("hit_charge_threshold").get_float();
 
  i4_split_minimum_hits = get_parameter ("split_minimum_hits").get_int();
  f4_wimp_fraction      = get_parameter ("wimp_fraction")     .get_float();
  f4_split_separation   = get_parameter ("split_separation")  .get_float();
  f4_split_window       = get_parameter ("split_window")      .get_float();
  f4_ds_splitting       = get_parameter ("ds_splitting")      .get_float();
  f4_dt_splitting       = get_parameter ("dt_splitting")      .get_float();

  if (i4_enable_histos) {
    clusters_nhits       = new TH1F ("clusters_nhits"      , "bx_muon_findcluster_time cluster nhits"      , 500 , 0, 500 );
    clusters_npmts       = new TH1F ("clusters_npmts"      , "bx_muon_findcluster_time cluster npmts"      , 208 , 0, 208 );
    clusters_start_times = new TH1F ("clusters_start_times", "bx_muon_findcluster_time cluster start times", 
    	4000, -time_buffer_ns, 0);
    clusters_charge      = new TH1F ("clusters_charge"     , "bx_muon_findcluster_time cluster charge"     , 250 , 0, 500 );

    barn_interface::get ()->store (barn_interface::file, clusters_nhits      , this);
    barn_interface::get ()->store (barn_interface::file, clusters_npmts      , this);
    barn_interface::get ()->store (barn_interface::file, clusters_start_times, this);
    barn_interface::get ()->store (barn_interface::file, clusters_charge     , this);
  }
  
  // Allocate space for internal storage
  binned_times_sss  .reserve (i4_n_bins); 
  binned_times_floor.reserve (i4_n_bins); 
  
  // Calculate the ripple threshold as number of bins for which P(>=5) < 0.1%
  double dark_hits_bin = 0.1;   //FIXME: to be taken from data base
  double p_higher;
  for (i4_ripple_bins = 1; i4_ripple_bins < 10; i4_ripple_bins++) {
    double mu = 2 * dark_hits_bin * i4_ripple_bins;
    double p0 = ::exp (-1 * mu); 
    double pi = p0;
    p_higher = 1 - p0;
    for (int32_t i = 1; i <= i4_ripple_count; i++) { 
      pi = pi * mu / i; 
      p_higher -= pi; 
    }
    if (p_higher > 1E-4) break;
  }
  if (i4_ripple_bins <= 2 && p_higher > 1E-2) 
    get_message (bx_message::error) << "too few ripple_bins " << i4_ripple_bins 
	    << " cause of too high dark rate " << dark_hits_bin << dispatch;
  else 
    get_message (bx_message::log) << "ripple window set to " << i4_ripple_bins 
	    << "bins with count " << i4_ripple_count << "hits and fake probability " << p_higher << dispatch; 
}

bool HitSortPredicate(const bx_muon_clustered_hit& a, const bx_muon_clustered_hit& b) {
  return a.get_time() < b.get_time();
}


bx_echidna_event* bx_muon_findcluster::doit (bx_echidna_event *ev) {

  const db_profile& profile_info = bx_dbi::get ()->get_profile ();
  float time_buffer_ns = profile_info.muon_tdc_range() * constants::muon::tdc::ns_per_clock;

  const bx_muon_event& er = ev->get_muon();
  bx_muon_clustered_event &ew = dynamic_cast<bx_muon_clustered_event&>(ev->get_muon());

  ew.b_has_cluster_sss   = false;
  ew.b_has_cluster_floor = false;

  //get_message (bx_message::debug) << "event " << ev->get_event_number () << dispatch;

  // Fill the 16ns binned time histogram (starting at the gate start)
  binned_times_sss.clear ();
  binned_times_floor.clear ();
  std::fill_n (binned_times_sss  .begin (), i4_n_bins, 0);
  std::fill_n (binned_times_floor.begin (), i4_n_bins, 0);
  for (int32_t i = 0; i < er.get_decoded_nhits (); i++) {
    const bx_muon_decoded_hit &decoded_hit = er.get_decoded_hit (i);
    const db_channel_muon *dbch = decoded_hit.get_db_channel ();
    if (!dbch->is_ordinary ()) continue;
    float hit_time = er.get_decoded_hit (i).get_time () + time_buffer_ns;
    int32_t bin = int32_t(hit_time / i4_bin_width);
    if (bin < 0 && bin > i4_n_bins) { 
      get_message (bx_message::warn) << "hit outside boundaries of binned distribution " << hit_time << std::endl; 
      continue;
    }
    if (dbch->is_sss   ()) binned_times_sss  [bin]++; 
    if (dbch->is_floor ()) binned_times_floor[bin]++;
  }

  // Start looping on bins (sss)
  bool has_cluster_sss = false;
  int32_t start_bin_sss=0, end_bin_sss=0; 
  int32_t n_hits_sss=0;
  for (int32_t bin = 0; bin < i4_n_bins-3; bin++) {

    if (!binned_times_sss[bin]) continue; // skip empty bins

    // Sum 3 consecutive bins 
    int32_t sum_of_3_bins = 0;
    for (int32_t j = 0; j < 3; j++) 
      sum_of_3_bins += binned_times_sss[bin + j]; // bin + j > i4_time_bins is fine since there is 20 spare bins

    // Compare the sum with low threshold
    if (sum_of_3_bins < i4_start_threshold) continue;

    // If here a cluster is found, define the start
    start_bin_sss = bin-1;
 
    // Search the end of cluster
    for (; bin < i4_n_bins-i4_ripple_bins; bin ++) {
      int32_t number_hits = 0;
      for (int32_t j = 0; j < i4_ripple_bins; j++) 
	number_hits += binned_times_sss[bin + j];  // bin + j > i4_time_bins is fine since there is 20 spare bins
      if (number_hits <= i4_ripple_count) {  // END of cluster found (ripple reagion)
	n_hits_sss += number_hits;
	bin += i4_ripple_bins;
	break;
      }
      n_hits_sss += binned_times_sss[bin];
    }
    end_bin_sss = bin;

    has_cluster_sss = true;

//    get_message (bx_message::info) << "event " << ev->get_event_number ()
//	    << " found upper candidate cluster in [" << start_bin_sss * i4_bin_width 
//	    << ", " << end_bin_sss * i4_bin_width << "] ns and " << n_hits_sss << " hits." << dispatch;
  } // end looping on bins (sss)

  // Start looping on bins
  bool has_cluster_floor = false;
  int32_t start_bin_floor=0, end_bin_floor=0; 
  int32_t n_hits_floor=0;
  for (int32_t bin = 0; bin < i4_n_bins-3; bin++) {

    if (!binned_times_floor[bin]) continue; // skip empty bins

    // Sum 3 consecutive bins 
    int32_t sum_of_3_bins = 0;
    for (int32_t j = 0; j < 3; j++) 
      sum_of_3_bins += binned_times_floor[bin + j]; // bin + j > i4_time_bins is fine since there is 20 spare bins

    // Compare the sum with low threshold
    if (sum_of_3_bins < i4_start_threshold) continue;

    // If here a cluster is found, define the start
    start_bin_floor = bin-1;
 
    // Search the end of cluster
    for (; bin < i4_n_bins-i4_ripple_bins; bin ++) {
      int32_t number_hits = 0;
      for (int32_t j = 0; j < i4_ripple_bins; j++) 
	number_hits += binned_times_floor[bin + j];  // bin + j > i4_time_bins is fine since there is 20 spare bins
      if (number_hits <= i4_ripple_count) {  // END of cluster found (ripple reagion)
	n_hits_floor += number_hits;
	bin += i4_ripple_bins;
	break;
      }
      n_hits_floor += binned_times_floor[bin];
    }
	
    end_bin_floor = bin;

    has_cluster_floor=true;
//    get_message (bx_message::info) << "event " << ev->get_event_number ()
//	    << " found floor candidate cluster in [" << start_bin_floor * i4_bin_width 
//	    << ", " << end_bin_floor * i4_bin_width << "] ns and " << n_hits_floor << " hits." << dispatch;
  } // end of loop on bins (floor)

  if (n_hits_sss+n_hits_floor < i4_count_threshold) return ev; // IGNORE this cluster if integral is too low     
//    get_message (bx_message::info) << "event " << ev->get_event_number () << " cluster selected, sss " 
//       << has_cluster_sss << " floor " << has_cluster_floor << dispatch; 

  // set same start and stop times for clustering search
  if (start_bin_sss > start_bin_floor) start_bin_sss = start_bin_floor;
  if (start_bin_sss < start_bin_floor) start_bin_floor = start_bin_sss;
  if (  end_bin_sss <   end_bin_floor)   end_bin_sss =   end_bin_floor;
  if (  end_bin_sss >   end_bin_floor)   end_bin_floor =   end_bin_sss;


  int32_t nhits_sss = 0, nhits_floor = 0, npmts = 0; //, nhits = 0 
  float start_time_sss   = time_buffer_ns;
  float start_time_floor = time_buffer_ns;
  float charge_sss = 0., charge_floor = 0.;
  chits_sss.clear(); chits_floor.clear();

  //get_message (bx_message::debug) << "sss_bin " << start_bin_sss << ", floor_bin " << start_bin_floor << dispatch;

  // loop to find start time of the sphere and of the floor
  for (int32_t hit=0; hit < er.get_decoded_nhits (); hit++) {
    const bx_muon_decoded_hit &decoded_hit = er.get_decoded_hit (hit);
    if (decoded_hit.get_db_channel ()->is_ordinary ()) {
      float hit_time   = er.get_decoded_hit (hit).get_time ()  + time_buffer_ns;
      float hit_charge = er.get_decoded_hit (hit).get_charge ();
      if (hit_charge < f4_hit_charge_threshold) continue;
      int32_t bin = int32_t(hit_time / i4_bin_width);
      if (decoded_hit.get_db_channel ()->is_sss () && bin >= start_bin_sss && bin < end_bin_sss) {
	if (hit_time < start_time_sss) start_time_sss = hit_time; 	
      }
      if (decoded_hit.get_db_channel ()->is_floor () && bin >= start_bin_floor && bin < end_bin_floor) {
	if (hit_time < start_time_floor) start_time_floor = hit_time; 	
      }
      //get_message (bx_message::debug) << "hit " << hit << ":\t t=" << hit_time << ", tb=" << bin << ", q=" << hit_charge << ", sss?" << decoded_hit.get_db_channel()->is_sss() << dispatch;
    } 
  }

  // if there is no sss or floor hit, set start time to the other one
  if (start_time_sss == time_buffer_ns) start_time_sss = start_time_floor;
  if (start_time_floor == time_buffer_ns) start_time_floor = start_time_sss;
  
  // loop over hits: fill charges, nhits and internal hit vectors
  for (int32_t hit = 0; hit < er.get_decoded_nhits (); hit++) {
    const bx_muon_decoded_hit &decoded_hit = er.get_decoded_hit (hit);
    if (decoded_hit.get_db_channel ()->is_ordinary ()) {
      float hit_time   = er.get_decoded_hit (hit).get_time ()  + time_buffer_ns;
      float hit_charge = er.get_decoded_hit (hit).get_charge ();
			if (hit_charge < f4_hit_charge_threshold) continue;
      int32_t bin = int32_t(hit_time / i4_bin_width);
      if (decoded_hit.get_db_channel ()->is_sss () && bin >= start_bin_sss && bin < end_bin_sss) {
        charge_sss += hit_charge;
        nhits_sss++;
	ew.nhits_per_channel[decoded_hit.get_raw_hit().get_muon_channel()]++;
	chits_sss.push_back (bx_muon_clustered_hit(hit_time-start_time_sss, hit_charge, &decoded_hit));
      }
      if (decoded_hit.get_db_channel ()->is_floor () && bin >= start_bin_floor && bin < end_bin_floor) {
        charge_floor += hit_charge;
        nhits_floor++;
	ew.nhits_per_channel[decoded_hit.get_raw_hit().get_muon_channel()]++;
	chits_floor.push_back (bx_muon_clustered_hit(hit_time-start_time_sss, hit_charge, &decoded_hit));
      }
    }
  }

//  for (uint32_t i=0; i<chits_floor.size(); i++) {
//    get_message (bx_message::debug) << "event " << ev->get_event_number () << ": floor hit #" << i << " t " <<  chits_floor[i].get_time() << dispatch;
//  }

  if (start_time_floor < start_time_sss) 
    get_message (bx_message::log) << "event " << ev->get_event_number () << " cluster floor " 
       << start_time_floor << " before cluster sss " << start_time_sss << dispatch;

  //get_message (bx_message::debug) << "event " << ev->get_event_number () << ": clusters present." << dispatch;
  //get_message (bx_message::debug) << "Up    ST " << start_time_sss-time_buffer_ns   << "ns, " << nhits_sss   << " hits, " << charge_sss   << "pe" << dispatch;
  //get_message (bx_message::debug) << "Floor ST " << start_time_floor-time_buffer_ns << "ns, " << nhits_floor << " hits, " << charge_floor << "pe" << dispatch;








  // Sort hit vectors
  std::sort(chits_sss  .begin(), chits_sss  .end(), HitSortPredicate);
  std::sort(chits_floor.begin(), chits_floor.end(), HitSortPredicate);

  // First clustering, affiliate each hit to a a cluster id.
  int32_t n_clusters_sss   = m_affiliate(0             , chits_sss,   f4_maxdt, f4_maxdr_sss  );
  int32_t n_clusters_floor = m_affiliate(n_clusters_sss, chits_floor, f4_maxdt, f4_maxdr_floor);
  int32_t n_clusters = n_clusters_sss + n_clusters_floor;
  
  //get_message (bx_message::debug) << "event " << ev->get_event_number () << " affiliation done, found " << n_clusters_sss << "+" << n_clusters_floor << " clusters" << dispatch;

  // Pre-compute charge of clusters
  std::vector<float> charge_v(n_clusters);
  for (uint32_t i = 0; i < chits_sss.size(); i++) {
    if (chits_sss[i].get_affiliation() > 0) 
      charge_v[chits_sss[i].get_affiliation()-1] += chits_sss[i].get_charge();
  }
  for (uint32_t i=0; i<chits_floor.size(); i++) {
    if (chits_floor[i].get_affiliation() > 0) 
      charge_v[chits_floor[i].get_affiliation()-1] += chits_floor[i].get_charge();
  }

  //for (int32_t i = 0; i < n_clusters; i++) get_message (bx_message::debug) << "#" << i << " : Q " << charge_v[i] << dispatch;

  // identify the one on which to attempt plitting.
  int32_t id_to_be_split = 0;
  float max_charge = 0;
  for (int32_t i = 0; i < n_clusters; i++) {
    if (charge_v[i] < f4_mincharge) { // cluster too small, removing....
      for ( uint32_t hit = 0 ; hit < chits_sss.size(); hit++) {	
        if (chits_sss[hit].get_affiliation() == i+1) chits_sss[hit].i4_affiliation = 0;
        if (chits_sss[hit].get_affiliation() > i+1) chits_sss[hit].i4_affiliation --;
      }
      for ( uint32_t hit = 0; hit < chits_floor.size(); hit++) {
        if (chits_floor[hit].get_affiliation() == i+1) chits_floor[hit].i4_affiliation = 0;
        if (chits_floor[hit].get_affiliation() > i+1) chits_floor[hit].i4_affiliation --;
      }
      if (i<n_clusters_sss) n_clusters_sss--;
      else n_clusters_floor--;
      n_clusters--;
      charge_v.erase(charge_v.begin()+i);
      i--;
    }
    else if (i<n_clusters_sss && charge_v[i] > max_charge) {
      max_charge = charge_v[i];
      id_to_be_split = i+1;
    }
  }

  //for (int32_t i = 0; i < n_clusters_sss; i++) get_message (bx_message::debug) << "sss#" << i << " : Q " << charge_v[i] << dispatch;
  //for (int32_t i = n_clusters_sss; i < n_clusters; i++) get_message (bx_message::debug) << "flo#" << i << " : Q " << charge_v[i] << dispatch;

  // try to split, can do it or not
  bool split = false;
  if (id_to_be_split!=0) split = m_split(id_to_be_split);

  //get_message (bx_message::debug) << "event " << ev->get_event_number () << " split #" << id_to_be_split << ": " << split << dispatch;
  if (split) {
//    get_message (bx_message::debug) << "n_hits_floor " << n_hits_floor << dispatch;
    n_clusters_sss++;
//    get_message (bx_message::debug) << "n_clusters_sss " << n_clusters_sss << dispatch;
    n_clusters++;
//     get_message (bx_message::debug) << "n_clusters " << n_clusters << dispatch;
    for (uint32_t i = 0; i < chits_floor.size(); i++) {
//      get_message (bx_message::debug) << i << "th hit on floor: affiliation is " << chits_floor[i].get_affiliation() << dispatch;
      if (chits_floor[i].get_affiliation()!=0) chits_floor[i].i4_affiliation++;
    }
  }
  //for (int32_t i = 0; i < n_clusters_sss; i++) get_message (bx_message::debug) << "after splitting: sss#" << i << " : Q " << charge_v[i] << dispatch;
  //for (int32_t i = n_clusters_sss; i < n_clusters; i++) get_message (bx_message::debug) << "flo#" << i << " : Q " << charge_v[i] << dispatch;

 //get_message (bx_message::debug) << "sss hits:" << dispatch;
 //for ( uint32_t hit = 0 ; hit < chits_sss.size(); hit++)
   //get_message (bx_message::debug) << "hit " << hit << " aff " << chits_sss[hit].get_affiliation() << dispatch;
 //get_message (bx_message::debug) << "floor hits:" << dispatch;
 //for ( uint32_t hit = 0 ; hit < chits_floor.size(); hit++)
   //get_message (bx_message::debug) << "hit " << hit << " aff " << chits_floor[hit].get_affiliation() << dispatch;


  // SSS clusters creation
  for ( int32_t clu = 0; clu < n_clusters_sss; clu++) {  
    float charge = 0.;
    double x = 0., y = 0., z = 0.;
    float radius_sss = chits_sss[0].get_decoded_hit().get_db_channel()->get_radius();
    float start_time = 1e6;
    int32_t hit_ctr = 0;
    for ( uint32_t hit = 0 ; hit < chits_sss.size(); hit++) {
      const db_channel_muon* dbc = chits_sss[hit].get_decoded_hit().get_db_channel();
      if (chits_sss[hit].get_affiliation() == clu+1) {
        if (start_time == 1e6) start_time = chits_sss[hit].get_time();
        //get_message (bx_message::debug) << "hit " << hit << " aff " << chits_sss[hit].get_affiliation() << " Q " << chits_sss[hit].get_charge() << " x " << dbc->get_x() << " t " << chits_sss[hit].get_time() << dispatch;
        if (clu != 0 && hit_ctr++ > i4_max_hits_sss) break;
        charge += chits_sss[hit].get_charge();
	x += dbc->get_x()*chits_sss[hit].get_charge()*exp(-(chits_sss[hit].get_time()-start_time)/f4_tau);
        y += dbc->get_y()*chits_sss[hit].get_charge()*exp(-(chits_sss[hit].get_time()-start_time)/f4_tau);
        z += dbc->get_z()*chits_sss[hit].get_charge()*exp(-(chits_sss[hit].get_time()-start_time)/f4_tau);
      }
    }
    float r = ::sqrt(x*x+y*y+z*z);
    x*= radius_sss/r; 
    y*= radius_sss/r; 
    z*= radius_sss/r; 
    //get_message (bx_message::debug) << "# " << clu << " ID " << clu+1 << " x " << x << " y " << y << " z " << z << " Q " << charge << " t " << start_time << dispatch;
    
    ew.clusters.push_back(bx_muon_cluster(clu+1, (float) x, (float) y, (float) z, charge, start_time)); 
  }

   //get_message (bx_message::debug) << n_clusters_sss << " sss cluster(s) done ..." << dispatch;

  // Floor clusters creation
  for (int32_t clu = n_clusters_sss; clu < n_clusters; clu++) {  
    float charge = 0., weight = 0.;
    double x = 0., y = 0., z = 0.;
    float start_time = 1e6; 
    int32_t hit_ctr = 0;
    for ( uint32_t hit = 0 ; hit < chits_floor.size(); hit++) {
      const db_channel_muon* dbc = chits_floor[hit].get_decoded_hit().get_db_channel();
      if (chits_floor[hit].get_affiliation() == clu+1) {
        //get_message (bx_message::debug) << "hit " << hit << " aff " << chits_floor[hit].get_affiliation() << " Q " << chits_floor[hit].get_charge() << " x " << dbc->get_x() << " t " << chits_floor[hit].get_time() << dispatch;
        if (start_time == 1e6) start_time = chits_floor[hit].get_time ();
	if (hit_ctr++ > i4_max_hits_floor) break;
        charge += chits_floor[hit].get_charge();
        x += dbc->get_x()*chits_floor[hit].get_charge()*exp(-(chits_floor[hit].get_time()-start_time)/f4_tau);
        y += dbc->get_y()*chits_floor[hit].get_charge()*exp(-(chits_floor[hit].get_time()-start_time)/f4_tau);
        z += dbc->get_z()*chits_floor[hit].get_charge()*exp(-(chits_floor[hit].get_time()-start_time)/f4_tau);
        weight += chits_floor[hit].get_charge()*exp(-(chits_floor[hit].get_time()-start_time)/f4_tau);
      }
    }
    x/= weight;
    y/= weight;
    z/= weight;
    //get_message (bx_message::debug) << "# " << clu << " ID " << clu+1 << " x " << x << " y " << y << " z " << z << " Q " << charge << " t " << start_time << dispatch;
    
    ew.clusters.push_back(bx_muon_cluster(clu+1, (float) x, (float) y, (float) z, charge, start_time)); 
  }
  //get_message (bx_message::debug) << "event " << ev->get_event_number () << " clusters created" << dispatch;

  // Copy hits into event
  ew.clustered_hits.insert(ew.clustered_hits.end(), chits_sss  .begin(), chits_sss  .end());
  ew.clustered_hits.insert(ew.clustered_hits.end(), chits_floor.begin(), chits_floor.end());

  // Fill event level variables
  //nhits = nhits_sss + nhits_floor;
  for (int32_t i = 0; i < constants::muon::channels; i++) { if (ew.nhits_per_channel[i]>0) npmts++; }
  ew.i4_npmts            = npmts;
  ew.f4_start_time_sss   = start_time_sss   - time_buffer_ns;
  ew.f4_start_time_floor = start_time_floor - time_buffer_ns;
  ew.i4_nhits_sss        = nhits_sss;
  ew.i4_nhits_floor      = nhits_floor;
  ew.f4_charge_sss       = charge_sss;
  ew.f4_charge_floor     = charge_floor;
  ew.b_has_cluster_sss   = has_cluster_sss; // mark the event as clustered
  ew.b_has_cluster_floor = has_cluster_floor; // mark the event as clustered

  //get_message (bx_message::debug) << "event " << ev->get_event_number () << " filled" << dispatch;

  // Fill histograms
//  if (i4_enable_histos) {
//    clusters_nhits      ->Fill(nhits);
//    clusters_npmts      ->Fill(npmts);
//    clusters_start_times->Fill(start_time-time_buffer_ns()); //AAA to be fixed;
//    clusters_charge     ->Fill(charge);
//  }
  i4_clustered_events++; // just statistics

  ev->get_muon().mark_stage (bx_base_event::clustered); 
  return ev;  
}

// assign an affiliation number to hits in the vector. This actually defines clusters.
// Space and time correlation is used.
int32_t bx_muon_findcluster::m_affiliate (int32_t offset, std::vector<bx_muon_clustered_hit>& v, float maxdt, float maxdr) {
  int32_t noclusters = offset;
  for (uint32_t i = 0; i < v.size(); i++) {
    for (uint32_t j = i+1; j < v.size(); j++) {
      const db_channel_muon* dbc_i = v[i].get_decoded_hit().get_db_channel();
      const db_channel_muon* dbc_j = v[j].get_decoded_hit().get_db_channel();
      int32_t *aff_i = &(v[i].i4_affiliation);
      int32_t *aff_j = &(v[j].i4_affiliation);
      if (*aff_i !=0 && *aff_i==*aff_j) continue;
      float dt = v[i].get_time() - v[j].get_time();
      if (dt < 0) dt *= -1;
      if (dt > maxdt ) break;  // hits too far in time. Can break because vector is ordered.
      float dx = dbc_i->get_x() - dbc_j->get_x();
      float dy = dbc_i->get_y() - dbc_j->get_y();
      float dz = dbc_i->get_z() - dbc_j->get_z();
      float dr = ::sqrt( dx*dx + dy*dy + dz*dz );
      if (dr > maxdr) continue; // hits too far in space
      // hit pair is related
      if (*aff_i==0 && *aff_j==0) { // both were not belonging to any cluster, start a new one
	noclusters++;
	*aff_i = noclusters;
	*aff_j = noclusters;
      }
      if (*aff_i==0 && *aff_j!=0) *aff_i=*aff_j; // affiliate i-th hit to j-th cluster
      if (*aff_i!=0 && *aff_j==0) *aff_j=*aff_i; // affiliate j-th hit to i-th cluster
      if (*aff_i!=0 && *aff_j!=0 && *aff_i!=*aff_j) { // hit were both affiliated to different clusters, j-th cluster to be removed.
        int32_t needless_aff = *aff_j; // the one to be removed
        //get_message (bx_message::debug) << " m_affiliate: size " << v.size() << " ;extra cluster found, removing " << needless_aff << dispatch;
	for (uint32_t h=0; h<v.size(); h++) {
          if (v[h].get_affiliation() == needless_aff) v[h].i4_affiliation = *aff_i;
	  if (v[h].get_affiliation() >  needless_aff) v[h].i4_affiliation--; // shift down affiliation for following clusters.
	} // end of loop on h
	noclusters--;
      } // end of extra cluster removal
    } // end of loop on j
  } //end of loop on i
  return noclusters-offset;
}

// tries to split the clusters recieved (the highest-charge cluster on sss which is not the first).
// If maxima are too close in time or space it fails.
/*bool bx_muon_findcluster::m_split(int32_t id_to_be_split) {
  int32_t noofhits = 0;
  for (uint32_t i=0; i<chits_sss.size(); i++) {
    if (chits_sss[i].get_affiliation() == id_to_be_split) noofhits++;
  }
  if (!noofhits) return false;

  // hits with low charge are neglected (wimps) for f4_wimp_fraction.
  int32_t wimps = 0;
  for (int32_t i=0; i<noofhits*f4_wimp_fraction; i++) {
    float mincharge = 1e4;
    int32_t minhit = 0;
    for (uint32_t j=0; j<chits_sss.size(); j++) { 
      if ((chits_sss[j].get_affiliation() == id_to_be_split) && chits_sss[j].get_charge() < mincharge) {
	mincharge = chits_sss[j].get_charge();
	minhit = j;
      } 
    }
    chits_sss[minhit].i4_affiliation = -1;
    wimps++;
  }
  get_message (bx_message::debug) << wimps << " of " << noofhits << " hits are wimps." << dispatch;

  // identify start & stop time for search for maxima
  float starttime=1e4, stoptime=0;
  for (uint32_t ihit=0; ihit<chits_sss.size(); ihit++) {
    if (chits_sss[ihit].get_affiliation() != id_to_be_split) continue;
    starttime=chits_sss[ihit].get_time();
    break;
  }
  for (uint32_t ihit=chits_sss.size()-1; ihit>0; ihit--) {
    if (chits_sss[ihit].get_affiliation() != id_to_be_split) continue;
    stoptime=chits_sss[ihit].get_time();
    break;
  }

  //look for maxima
  int32_t   mh_1 = -1 , mh_2 = -1;
  float max1 = 0., max2 = 0.;
  int32_t goodhits = 0;
  for (uint32_t i=0; i<chits_sss.size(); i++) {
    if (chits_sss[i].get_affiliation() != id_to_be_split) continue;
    goodhits ++;
    if (goodhits<(noofhits-wimps)*f4_split_separation && chits_sss[i].get_time()-starttime<f4_split_window && chits_sss[i].get_charge()>max1) {
      max1 = chits_sss[i].get_charge();
      mh_1 = i;
    }
    if (goodhits>(noofhits-wimps)*f4_split_separation && stoptime-chits_sss[i].get_time()<f4_split_window && chits_sss[i].get_charge()>max2) {
      max2 = chits_sss[i].get_charge();
      mh_2 = i;
    }
  }
  if (mh_1==-1 || mh_2==-1) {
    for (uint32_t ihit=0; ihit<chits_sss.size(); ihit++) {
      if (chits_sss[ihit].get_affiliation() == -1) chits_sss[ihit].i4_affiliation = id_to_be_split;
    }
    return false;
  }
}*/

	// alternative splitting method
bool bx_muon_findcluster::m_split(int32_t id_to_be_split) {
	//get_message (bx_message::debug) << " splitting started ..." << dispatch;
  int32_t noofhits = 0;
  for (uint32_t i=0; i<chits_sss.size(); i++) {
    if (chits_sss[i].get_affiliation() == id_to_be_split) noofhits++;
  }
  if (noofhits<i4_split_minimum_hits) return false;

  // hits with low charge are neglected (wimps) for f4_wimp_fraction.
  int32_t wimps = 0;
  for (int32_t i=0; i<noofhits*f4_wimp_fraction; i++) {
    float mincharge = 1e4;
    int32_t minhit = 0;
    for (uint32_t j=0; j<chits_sss.size(); j++) { 
      if ((chits_sss[j].get_affiliation() == id_to_be_split) && chits_sss[j].get_charge() < mincharge) {
	mincharge = chits_sss[j].get_charge();
	minhit = j;
      } 
    }
    chits_sss[minhit].i4_affiliation = -1;
    wimps++;
  }
  //get_message (bx_message::debug) << wimps << " of " << noofhits << " hits are wimps." << dispatch;

  // use first hit as "entry point". look if there are no-wimp hits more than minimum distance away

  float entry_t = 0, entry_x = 0, entry_y = 0, entry_z = 0;

  for (uint32_t ihit=0; ihit<chits_sss.size(); ihit++) {
    if (chits_sss[ihit].get_affiliation() != id_to_be_split) continue;
      entry_x = chits_sss[ihit].get_decoded_hit().get_db_channel()->get_x();
      entry_y = chits_sss[ihit].get_decoded_hit().get_db_channel()->get_y();
      entry_z = chits_sss[ihit].get_decoded_hit().get_db_channel()->get_z();
      entry_t = chits_sss[ihit].get_time();
      break;
    }

  // search for center of charge for exit point

  float exit_x = 0, exit_y = 0, exit_z = 0, exit_q = 0, exit_t = 0;
  float sum_x = 0, sum_y = 0, sum_z = 0, sum_q = 0, distance = 0;

  for (uint32_t ihit=0; ihit<chits_sss.size(); ihit++) {
    if (chits_sss[ihit].get_affiliation() != id_to_be_split) continue;
    exit_x = chits_sss[ihit].get_decoded_hit().get_db_channel()->get_x();
    exit_y = chits_sss[ihit].get_decoded_hit().get_db_channel()->get_y();
    exit_z = chits_sss[ihit].get_decoded_hit().get_db_channel()->get_z();
    distance = sqrt( pow(exit_x-entry_x,2) + pow(exit_y-entry_y,2) + pow(exit_z-entry_z,2) );
    if (distance < f4_ds_splitting) continue;
    if (chits_sss[ihit].get_time() < f4_dt_splitting) continue;
    if (exit_t == 0) exit_t = chits_sss[ihit].get_time();
    exit_q = chits_sss[ihit].get_charge();
    sum_x += exit_q*exit_x;
    sum_y += exit_q*exit_y;
    sum_z += exit_q*exit_z;
    sum_q += exit_q;
  }

  // if sum charge of the exit point is 0, return hits to original affiliation, stop splitting

  if (sum_q==0) {
    for (uint32_t ihit=0; ihit<chits_sss.size(); ihit++) {
    if (chits_sss[ihit].get_affiliation() == -1) chits_sss[ihit].i4_affiliation = id_to_be_split;
    }
    //get_message (bx_message::debug) << "nothing to be split" << dispatch;
    return false;
  }

  exit_x = sum_x/sum_q;
  exit_y = sum_y/sum_q;
  exit_z = sum_z/sum_q;

  //get_message (bx_message::debug) << "entry t " << entry_t << " : exit_t " << exit_t << " : id_to_be_split " << id_to_be_split << dispatch;

  int32_t count_c1_hits=0, count_c2_hits=0;

  // Splitting confirmed, reassign affiliation.
  for (uint32_t i=0; i<chits_sss.size(); i++) {
    if ((chits_sss[i].get_affiliation() < id_to_be_split) && (chits_sss[i].get_affiliation() >= 0)) continue;
    else if (chits_sss[i].get_affiliation() > id_to_be_split) chits_sss[i].i4_affiliation++;
    else {
      const db_channel_muon* dbc = chits_sss[i].get_decoded_hit().get_db_channel();
      float d1 = pow(dbc->get_x()-entry_x, 2) +  pow(dbc->get_y()-entry_y, 2) + pow(dbc->get_z()-entry_z, 2);
      float d2 = pow(dbc->get_x()-exit_x, 2) +  pow(dbc->get_y()-exit_y, 2) + pow(dbc->get_z()-exit_z, 2);
	    //get_message (bx_message::debug) << "hit t " << chits_sss[i].get_time() << dispatch;
      if (d2<d1 &&  ( chits_sss[i].get_time() > (exit_t-entry_t)*0.8 ) && chits_sss[i].get_affiliation()>0) {
        chits_sss[i].i4_affiliation = id_to_be_split+1;
	count_c2_hits++;
      }
      else if (d2>d1 && chits_sss[i].get_time() < (exit_t-entry_t)*1.1 ) {
        chits_sss[i].i4_affiliation = id_to_be_split;
	count_c1_hits++;
      } 
      else chits_sss[i].i4_affiliation = -1;
    }
  }
 
  // if one of the clusters has less than two hits, restore the old clusters 

  if (count_c1_hits<2 || count_c2_hits<2) {
    for (uint32_t i=0; i<chits_sss.size(); i++) {
      if (chits_sss[i].get_affiliation() > id_to_be_split) chits_sss[i].i4_affiliation--;
      if (chits_sss[i].get_affiliation()==-1) chits_sss[i].i4_affiliation = id_to_be_split;
    }
    return false;
  }
  else {
    for (uint32_t i=0; i<chits_sss.size(); i++) {
      if (chits_sss[i].get_affiliation()==-1) chits_sss[i].i4_affiliation = 0;
    }
  }
   

  //for (uint32_t i=0; i<chits_sss.size(); i++) get_message (bx_message::debug) << i << "th hit affiliation: " << chits_sss[i].get_affiliation() << dispatch;
  //get_message (bx_message::debug) << "splitting converged!" << dispatch;
  return true;
}

void bx_muon_findcluster::end () {
  get_message (bx_message::info) << "run summary: found " << i4_clustered_events << " clustered events" << dispatch;
}

/*
 * $Log: bx_muon_findcluster.cc,v $
 * Revision 1.35  2009/11/08 15:35:06  wurm
 * fixed bug in cluster splitting that caused empty clusters
 *
 * Revision 1.34  2009-10-26 11:19:37  ddangelo
 * trgtype 'muon_tct' (or 'muon_totc') renamed as 'neutron' after real use
 * trgtype 'muon_mtb' renamed as 'muon' as duplicity with the above one disappeared and 'mtb' is associated to btb_inputs flag
 *
 * Revision 1.33  2009-10-23 09:27:49  wurm
 * commented debug messages
 *
 * Revision 1.32  2008-09-30 10:10:20  wurm
 * set start and stop time bins for sss and floor to the same value
 *
 * Revision 1.31  2008-09-30 07:46:21  wurm
 * cured large exponents in weight functions
 *
 * Revision 1.30  2008-07-22 07:29:23  wurm
 * removed a bug from floor clusters time alignmnent
 *
 * Revision 1.29.4.1  2008-04-21 13:41:19  wurm
 * removed bug in splitting
 *
 * Revision 1.29  2008-02-26 11:35:51  ddangelo
 * still quiter
 *
 * Revision 1.28  2008-02-22 12:02:40  ddangelo
 * quiter
 *
 * Revision 1.27  2008-02-04 10:43:10  wurm
 * removed a bug from splitting
 *
 * Revision 1.26  2008-02-02 15:29:03  wurm
 *
 *
 * read hit_charge_threshold from echidna.cfg
 *
 * Revision 1.25  2008-02-02 15:23:20  wurm
 *
 *
 * introduced hit charge threshold for clustering
 *
 * Revision 1.24  2008-01-06 13:12:10  ddangelo
 * lowered verbosity
 *
 * Revision 1.23  2007-12-20 18:40:40  ddangelo
 * fixed a bad bug.
 * writing n_hits individual for sss and floor
 *
 * Revision 1.22  2007-12-07 17:33:59  ddangelo
 * compliant with new dynamically loaded tdc clock range
 *
 * Revision 1.21  2007-11-28 19:28:42  wurm
 *
 * removed bug from affiliation
 *
 * Revision 1.20  2007-11-28 18:07:31  wurm
 *
 * removed bug in splitting
 *
 * Revision 1.19  2007-11-28 16:03:19  wurm
 *
 *
 * lowered verbosity
 *
 * Revision 1.18  2007-11-27 19:01:51  wurm
 *
 *
 * removed some bugs-
 *
 * Revision 1.17  2007-11-26 18:50:26  ddangelo
 * debugging, still problems
 *
 * Revision 1.16  2007-11-26 16:16:07  ddangelo
 * cleaned up
 *
 * Revision 1.15  2007-11-26 14:06:26  ddangelo
 * completely redone
 *
 * Revision 1.14  2007-11-14 19:26:55  ddangelo
 * quiter
 *
 * Revision 1.13  2007-11-14 19:00:14  ddangelo
 * writing to event
 *
 * Revision 1.12  2007-11-14 15:46:33  ddangelo
 * new clustering
 * detector segmented in 2/3 parts.
 * individual hits array and cluster definition
 * event writing to be completed
 *
 * Revision 1.11  2007-08-30 15:14:16  saggese
 * Added muon events to clustering (non manteiner commit)
 *
 * Revision 1.10  2007-05-30 16:04:10  ddangelo
 * complaiant to new has_cluster definition
 *
 * Revision 1.9  2007-05-25 17:01:16  ddangelo
 * quiter
 *
 * Revision 1.8  2007-05-25 15:09:56  ddangelo
 * minor things
 *
 * Revision 1.7  2007-03-28 17:50:16  ddangelo
 * code restyled and cleaned up. no longer helper functions.
 * n_hits_per_channel is now filled. fired channels removed.
 * binned array resized correctly. no extra bins.
 * constants used instead of wired numbers.
 *
 * Revision 1.6  2007-03-27 15:21:20  ddangelo
 * cluster start time savbed as negative with respect to trg
 *
 * Revision 1.5  2007-03-22 16:13:02  ddangelo
 * starting to do something...
 *
 * 
 * Revision 1.4  2007-03-10 09:46:05  pallas
 * Removing some warnings (variable init in all cases)
 *
 * Revision 1.3  2007-02-22 19:59:55  ddangelo
 * first run tests, still to go
 *
 * Revision 1.2  2007/02/21 18:50:09  ddangelo
 * first development. still junk
 *
 * Revision 1.1  2007/02/21 15:48:59  ddangelo
 * added. blank.
 *
 */

