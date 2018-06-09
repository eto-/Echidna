/* BOREXINO Reconstruction program
 *
 * Author: Daniela Manuzio <dmanuzio@ge.infn.it>
 * Maintainer: Alessandro Razeto <Alessandro.Razeto@ge.infn.it>
 *
 * $Id: bx_laben_findcluster.cc,v 1.134 2015/07/20 15:23:37 ilia.drachnev Exp $
 *
 * Implemenentation of bx_laben_findcluster
 * Module that finds the candidate fragments in the laben hit 
 * collection. No real splitting is done in this module.
 * Input:  the list of time ordered laben hits
 * Output: a list of bx_cluster objects
 *
 * Algorithm:
 * Candidate: 1. number of hits in 1 bin > (dark_thresh) 
 *            2. number of hits in 3 consecutive bins > (low_thresh)
 * First:     1. number of hits in 3 consecutive bins > (trg_thresh)
 *            2. right delay (trg_delay) with trigger reference time            
 * Following: 1. time distance from previous cluster > (time_separation)
 * 
*/

#include <cmath>
#include <iostream>
#include <vector>
#include <algorithm>
#include "bx_laben_findcluster.hh"
#include "bx_echidna_event.hh"
#include "bx_dbi.hh"
#include "db_run.hh"
#include "db_channel.hh"
#include "constants.hh"
#include "barn_interface.hh"
#include "mach4/Parameters.h"
#include "TMath.h"
#include "TF1.h"


bx_laben_findcluster::bx_laben_findcluster (): bx_base_module("bx_laben_findcluster", bx_base_module::main_loop) {
  require_event_stage (bx_detector::laben, bx_base_event::decoded);
  require_trigger_type (bx_trigger_event::neutrino);
  require_trigger_type (bx_trigger_event::muon);
  require_trigger_type (bx_trigger_event::neutron);
  require_trigger_type (bx_trigger_event::random);
}

void bx_laben_findcluster::begin () {
  i4_he_threshold = get_parameter ("he_threshold").get_int ();
  i4_start_threshold = get_parameter ("start_threshold").get_int ();
  i4_neutron_start_threshold = get_parameter ("neutron_start_threshold").get_int ();
  i4_count_threshold = get_parameter ("count_threshold").get_int ();
  i4_ripple_bins = get_parameter ("ripple_bins").get_int ();
  i4_he_ripple_bins = get_parameter ("he_ripple_bins").get_int ();
  i4_tail_bins = get_parameter ("tail_bins").get_int ();
  i4_long_bins = get_parameter ("long_bins").get_int ();
  b_force_one_large_cluster = get_parameter ("force_one_large_cluster").get_bool ();
  b_strict_gate = get_parameter ("strict_gate").get_bool ();
  b_mach4_cluster = get_parameter ("mach4_cluster").get_bool ();
  
    // Define some parameters
  i4_gate_start = int(bx_dbi::get ()->get_run ().get_laben_gate_start ());
  i4_gate_end = int(bx_dbi::get ()->get_run ().get_laben_gate_width () + i4_gate_start + 500);
  i4_cluster_offset = int(bx_dbi::get ()->get_run ().get_laben_cluster_offset ());
  const float pre_width = bx_dbi::get ()->get_run ().get_laben_gate_width();
  f4_trigger_start_ = bx_dbi::get ()->get_run ().get_laben_gate_start();
  //!(int) rounds down, so we miss the last bit if pre_width is not an integer multiple of _length_
  i4_NFrames_win1_=(int)(pre_width/i4_dt1_len_);
  i4_NFrames_win2_=(int)(pre_width/i4_dt2_len_);
  get_message (bx_message::info)<<"pre_width = "<<pre_width<<dispatch;
  get_message (bx_message::info)<<"f4_trigger_start_ = "<<f4_trigger_start_<<dispatch;
  get_message (bx_message::info)<<"i4_NFrames_win1_ = "<<i4_NFrames_win1_<<dispatch;
  get_message (bx_message::info)<<"i4_NFrames_win2_ = "<<i4_NFrames_win2_<<dispatch;

    // Allocate space for internal storage
  i4_time_bins = (1 << 16) * 50 * 2 / 16;
  binned_times = new short int[i4_time_bins + 20]; // 20 bins are added to have some extra space (set to 0)
  				 		   // on which to loop for (without any need for boundary check)
  fired_channel = new int[constants::laben::channels];
  prev_gps_times[0] = prev_gps_times[1] = 0;
  
  //double dark_rate = bx_dbi::get ()->get_run ()->get_laben_mean_dark_noise () * 2000 * 16e-9
  f_dark_rate = 700; //in Hz, average per PMT

    // Evaluate ripple count (based on dark rate and ripple_bins)
  i4_ripple_count = evaluate_ripple_count (f_dark_rate, i4_ripple_bins, i4_count_threshold, "normal");
  i4_he_ripple_count = evaluate_ripple_count (f_dark_rate, i4_he_ripple_bins, i4_count_threshold, "he");
}

bx_echidna_event* bx_laben_findcluster::doit (bx_echidna_event *ev) {
  bx_laben_event& er = ev->get_laben ();
//noavg code inserted: db_run instance
  run_info = &(bx_dbi::get ()->get_run ());
/*******************************************************************/  

  findcluster (ev);
  findcluster_time (ev);
  
  /*
  for (int i = 0; i < ev->get_laben ().get_nclusters (); i++) {
    double t1 = ev->get_laben ().get_trigger_rawt ();
    double t2 = ev->get_laben ().get_cluster (i).get_start_time ();
    if (::fabs(t1 - t2) > 30000) 
      get_message (bx_message::error) << "internal error, trigger_time - cluster_time(" << i << ") too large (" << t1 - t2 << ") for event " << ev->get_event_number () <<  dispatch;
  }*/
    
  er.mark_stage (bx_base_event::clustered); 
  return ev;  
}

void bx_laben_findcluster::end () {
  delete [] fired_channel;
}

void bx_laben_findcluster::findcluster (bx_echidna_event *ev) {
  bx_laben_event& er = ev->get_laben ();
  bx_laben_clustered_event &ew = dynamic_cast<bx_laben_clustered_event&>(er);

  //random offset for filling array variables for tt1
  int random_offset1 = 0;
  int random_offset2 = 0;
  
  if(ev->get_trigger ().is_neutrino()){
    random_offset1 = rand () % i4_dt1_len_;
    random_offset2 = rand () % i4_dt2_len_;
  }
 
  // Clear internal storage
  clusters.clear ();
  ::memset (binned_times, 0, (i4_time_bins + 20) * sizeof (short int)); // memset is much faster than std::fill_n
  //std::fill_n (binned_times, i4_time_bins + 20, 0); // always the same 20

  //! Prepare the array variables to be filled in the next loop
  std::vector<int> empty_vec(constants::laben::channels,0);
  for(int ii = 0; ii<i4_NFrames_win1_; ++ii) {
    (ew.v_npmts_win1).push_back(0);
    (ew.v_charge_win1).push_back(0);
  }
  for(int ii = 0; ii<i4_NFrames_win2_; ++ii) {
    (ew.v_npmts_win2).push_back(0);
    (ew.v_charge_win2).push_back(0);
  }

  std::vector<std::vector<int> > fc_win1(i4_NFrames_win1_,empty_vec);
  std::vector<std::vector<int> > fc_win2(i4_NFrames_win2_,empty_vec);

    // FILL the 16ns binned time histogram and the array variables
  for (int i = 0; i < er.get_decoded_nhits (); i++) {

    const bx_laben_decoded_hit &decoded_hit = er.get_decoded_hit (i);
    const db_channel_laben* db = decoded_hit.get_db_channel ();
    const int lg = db->get_lg ();

    if (!(decoded_hit.get_flag () & ~bx_laben_decoded_hit::out_of_gate) && decoded_hit.get_db_channel ()->is_ordinary ()) {
      //! Accept only ordinary channels and hits that either have flag=0 or have flag=out_of_gate
      float hit_time = decoded_hit.get_raw_time ();
      int bin = int(hit_time / 16.);
      if (bin >= 0 && bin < i4_time_bins && binned_times[bin] < 0x7FF0) binned_times[bin]++;

      //! Fill array variables; hit time is double, but I don't want to modify old code
      const double d_hit_time = decoded_hit.get_raw_time ();
      if(!(ev->get_trigger ().is_random() || ev->get_trigger ().is_neutrino ()))
	continue;
      const double ref_time = er.get_trigger_rawt ()+f4_trigger_start_;
      const double subt_time = d_hit_time - ref_time;
      if(subt_time<0)
      {
	get_message (bx_message::debug)<<"UNDERFLOW d_hit_time = "<<(int)d_hit_time<<", trigger time = "<<(int)er.get_trigger_rawt ()<<",f4_trigger_start_ = "<<f4_trigger_start_<<dispatch;
	continue;
      }
      
      const double subt_time_random1 = subt_time +  random_offset1;
      const double subt_time_random2 = subt_time +  random_offset2;


      const int frame_win1 = get_gate_index_(i4_dt1_len_,subt_time_random1);
      const int frame_win2 = get_gate_index_(i4_dt2_len_,subt_time_random2);

      if(frame_win1>i4_NFrames_win1_ || frame_win2>i4_NFrames_win2_)
      {
	get_message (bx_message::debug)<<"OVERFLOW frame_win1="<<frame_win1<<", frame_win2="<<frame_win2<<", and there are "<<i4_NFrames_win1_<<" frames in win1 and "<<i4_NFrames_win2_<<" in win2."<<dispatch;
	continue;
      }
      if(frame_win1!=i4_NFrames_win1_)
      {
	if((fc_win1.at(frame_win1)).at(lg-1)==0)
	{
	  ++(ew.v_npmts_win1).at(frame_win1);
          (ew.v_charge_win1).at(frame_win1) += decoded_hit.get_charge_pe();
	  (fc_win1.at(frame_win1)).at(lg-1)=1;
	}
      }
      if(frame_win2!=i4_NFrames_win2_)
      {
	if((fc_win2.at(frame_win2)).at(lg-1)==0)
	{
	  ++(ew.v_npmts_win2).at(frame_win2);
          (ew.v_charge_win2).at(frame_win2) += decoded_hit.get_charge_pe();
	  (fc_win2.at(frame_win2)).at(lg-1)=1;
	}
      }
    }
  } //End loop over decoded hits

    // Define the stating point
  long unsigned gps_times[2];
  ev->get_trigger ().get_gps_time (gps_times[0], gps_times[1]);
  long int gps_dt_s = gps_times[0] - prev_gps_times[0];
  long int gps_dt_ns = gps_times[1] - prev_gps_times[1];
  double dt_gps = double(gps_dt_s) * 1e9 + double(gps_dt_ns);
  ew.f4_window_limit = dt_gps;
  if (dt_gps < 0) get_message (bx_message::error) << "internal error negative gps_dt (" << dt_gps << ") for event " << ev->get_event_number () <<  dispatch;
  double start_time = er.get_trigger_rawt () - 20000; // default start time for the gate
  if (ev->get_trigger ().is_neutron () && i4_gate_start < -8e3) { // Only for new runs with large neutron gate
    start_time = er.get_trigger_rawt () - 1.7e6;
  } else {
    if (b_strict_gate) start_time = er.get_trigger_rawt () + i4_gate_start;
    //else if (dt_gps < 100 && dt_gps > 0) start_time = er.get_trigger_rawt () - dt_gps;
  }

  int start_bin = int(start_time / 16.);
  if (start_bin < 0 || start_bin > i4_time_bins) {
    get_message (bx_message::error) << "internal error start_time (" << start_time << ") out of range for event " << ev->get_event_number () <<  dispatch; 
    return;
  }
  int bin_end = int((er.get_trigger_rawt () + i4_gate_end) / 16.);
  if (bin_end > i4_time_bins) {
    get_message (bx_message::warn) << "internal error end_time (" << bin_end << ") out of range for event " << ev->get_event_number () <<  dispatch; 
    bin_end = i4_time_bins;
  }
  //Now call clusterize ()
  if (ev->get_trigger ().is_neutron ()) clusterize_neutrons (start_bin, bin_end, ev);
  else{
		if( b_mach4_cluster ) clusterize_mach4 (start_bin, bin_end, ev);
		else  clusterize (start_bin, bin_end, ev);

		//bool_muon_clustered and i4_nhits_background are not filled in the main-clustering-subroutines
		//Set to zero here
		int n_main_clusters = clusters.size ();		//these are the main clusters for non-muons
		for (int cl = 0; cl < n_main_clusters; cl++){
			clusters[cl].bool_muon_clustered = 0;		//these clusters were not identified by second muon-clustering (neutron-finding-algorithm), but identified by main-clustering
			clusters[cl].f4_n_hits_background = 0;		//background for main-clustering clusters is calculated in routine "bx_laben_findcluster::findcluster_time"
		}
	}


  //For muons run additional neutron_clustering algorithm
  if (ev->get_trigger ().get_btb_inputs () == 4  && ev->get_trigger ().get_trgtype () == 1) clusterize_neutrons_in_muongate (start_bin, bin_end, ev);

  //At this point, we are done with clustering
  int n_clusters = clusters.size ();		//now includes also possible additional muonclusters
  int num_neutron_clusters_above_cycle13_tresh = 0;	//for recoverability cycle13 within cycle14 of tt128 clusters, count the number of tt128 clusters above the old cycle13 threshold
  int num_clusters_above_neutron_tresh = 0;	// counts the number of clusters above neutron detection threshold (both in tt1 && tt128 events)
  if (n_clusters == 0) ; //get_message (bx_message::info) << "no candidate cluster found for event " << ev->get_event_number () << dispatch;
  else {
      // FOUND some clusters, handle them
      
      // First do some checks, fix some pathological conditions and log
    do { // a do {} while (0) is there only to logically bind toghether a block of statement
        // 1) last cluster not closed
      if (clusters[n_clusters - 1].i4_long_end_bin >= bin_end) get_message (bx_message::info) << "end of cluster " << n_clusters - 1 << " not found for event " << ev->get_event_number () << dispatch;
        // 2) no cluster which is in the right position for the trigger
	// ...
        // 3) too many clusters
      /*if (n_clusters > 3) {
        get_message(bx_message::log) << "found " << n_clusters << " candidate clusters for event " << ev->get_event_number () << " stripping to 3 " << dispatch;
          // simply skip cluster after 3 (in future keep the most significative ones)
        clusters.erase (clusters.begin () + 3, clusters.end ());
        n_clusters = 3;
      } //else get_message (bx_message::debug) << "found " << n_clusters << " candidate clusters for event " << ev->get_event_number () << dispatch;
      */
    } while (0); 

      // Then fill the bx_laben_clustered_hit array with the hits' raw times
      // (still within n_clusters!=0)
    float empty_board_factor = 1. - ev->get_laben().get_empty_boards() / 280.;
    float pmt_factor = (ev->get_laben ().get_n_live_pmts() - ev->get_laben ().get_invalid_pmts() ) / 2000.;
    for (int cl = 0; cl < n_clusters; cl++) {	//loops over all clusters (if muon-event : over main-clusters and muon-clusters identified by neutron algorithm)
      if (clusters[cl].i4_short_end_bin <= 0) clusters[cl].i4_short_end_bin = clusters[cl].i4_long_end_bin;
   
      bx_laben_cluster::bx_laben_cluster_vector *cluster_vector;
      if (!clusters[cl].bool_muon_clustered) cluster_vector = &ew.clusters;
      else cluster_vector = &ew.clusters_muons;  			//cluster in the muon event was clustered by neutron_finding_algorithm

      cluster_vector->push_back (bx_laben_cluster ());
      bx_laben_cluster &c = (*cluster_vector)[cluster_vector->size () - 1];

      //Calculate the num_neutron_clusters_above_cycle13_tresh for tt128
      if (ev->get_trigger ().is_neutron () && clusters[cl].i4_n_hits >= 40. * empty_board_factor * pmt_factor && clusters[cl].i4_n_hits >= 5.) num_neutron_clusters_above_cycle13_tresh++;

      c.f4_clustered_nhits_bkg = clusters[cl].f4_n_hits_background;

      for (int hit = 0; hit < er.get_decoded_nhits (); hit++) {
	const bx_laben_decoded_hit &decoded_hit = er.get_decoded_hit (hit);
	if (!(decoded_hit.get_flag () & ~bx_laben_decoded_hit::out_of_gate) && decoded_hit.get_db_channel ()->is_ordinary ()) {
	  float hit_time = decoded_hit.get_raw_time () / 16;
	  if (hit_time >= clusters[cl].i4_start_bin && hit_time <= clusters[cl].i4_long_end_bin) {
	    c.clustered_hits.push_back (bx_laben_clustered_hit(decoded_hit, hit));
	    bx_laben_clustered_hit &clustered_hit = c.clustered_hits[c.clustered_hits.size () - 1];
	    clustered_hit.f8_time = decoded_hit.get_raw_time();
	    if (hit_time < clusters[cl].i4_short_end_bin) clustered_hit.b_short_cluster = true;
	    else clustered_hit.b_short_cluster = false;
	  }
	}
      }

      if (!c.clustered_hits.size ()) {
	get_message (bx_message::error) << c.clustered_hits.size () << " hits are copied in cluster " << cl + 1 << " for event " << ev->get_event_number () << " while rough integral is " << clusters[cl].i4_n_hits << dispatch;
	cluster_vector->pop_back ();
      }
    }
   
    // for tt128, test whether clusters fulfill the minimum hit vs. empty board condition for neutrons
    if (ev->get_trigger().is_neutron ()) {
      bx_laben_cluster::bx_laben_cluster_vector *cluster_vector;
      cluster_vector = &ew.clusters;
      for (unsigned int cl=0; cl < (*cluster_vector).size(); cl++) {
        if ((*cluster_vector)[cl].get_clustered_nhits() > pmt_factor * (385 - 2*ev->get_laben ().get_empty_boards() ) ) {
          num_clusters_above_neutron_tresh++;
          (*cluster_vector)[cl].b_is_neutron = 1;
        }
      }
    }
  } 
  //Now fill the n_clusters_threshold-variable with meanungful values; for non tt128 events, it is identical to the cluster-size of laben.clusters
  if (ev->get_trigger ().is_neutron ()) ew.i4_nclusters_old = num_neutron_clusters_above_cycle13_tresh;
  else ew.i4_nclusters_old = ew.clusters.size();
  // fill the n_clusters_neutron_variable
  ew.i4_nclusters_neutron = num_clusters_above_neutron_tresh; 

  ev->get_trigger ().get_gps_time (prev_gps_times[0], prev_gps_times[1]);
}

void bx_laben_findcluster::clusterize (int start_bin, int bin_end, const bx_echidna_event *ev) {
    // Search for the cluster candidates 
  int start_threshold = i4_start_threshold;
  if (ev->get_trigger ().is_neutron ()) start_threshold = i4_neutron_start_threshold;
  for (int bin = start_bin; bin < bin_end; bin++) {
      // Avoid looping on void bins
    if (binned_times[bin]) {
        // Then sums the count of 3 consecutive bins and comare it with the low threshold
      int number_hits = 0;
      for (int j = 0; j < 3; j++) number_hits += binned_times[bin + j]; // bin + j > i4_time_bins is fine since there is 20 spare bins
      if (number_hits > start_threshold) {
	  // If here a cluster is found, define the start
	cluster_data data;
	data.i4_n_hits = 0;
        data.i4_start_bin = (bin > 0) ? bin - 1 : 0; // go one bin back for better finding first hit
        data.b_is_neutron = 0;
	
	  // IF force_one_large_cluster force only one LARGE cluster op to end of gate
	if (b_force_one_large_cluster) {
	  data.i4_long_end_bin = bin_end - 1;
	  for (data.i4_n_hits = 0; bin < bin_end; bin++) data.i4_n_hits += binned_times[bin];
	  clusters.push_back (data);
	  break;
	}

	  // Select operating mode: High Energy or normal
	int draft_energy = 0;
	for (int j = 0; j < 10; j++) draft_energy += binned_times[bin + j];  // bin + j > i4_time_bins is fine since there is 20 spare bins
	int ripple_bins, ripple_count;
	if (draft_energy > i4_he_threshold) {
	  ripple_bins = i4_he_ripple_bins;
	  ripple_count = i4_he_ripple_count;
	} else {
	  ripple_bins = i4_ripple_bins;
	  ripple_count = i4_ripple_count;
	}

	  // Search the end of cluster
	for (; bin < bin_end; bin ++) {
          number_hits = 0;
	  for (int j = 0; j < ripple_bins; j++) number_hits += binned_times[bin + j];  // bin + j > i4_time_bins is fine since there is 20 spare bins
	  if (number_hits <= ripple_count) {  // END of cluster found (ripple reagion)
	    data.i4_n_hits += number_hits;
//	    bin += ripple_bins; // will be included later by tail_bins
	    break;
	  }
	  data.i4_n_hits += binned_times[bin];
	}
	
	if (data.i4_n_hits < i4_count_threshold) continue; // IGNORE this cluster if integral is too low

	  // Add the tail contribution which define the very last bin
	int tail_bins = i4_tail_bins;
	int eoc_bin = bin + tail_bins;
	int eoc_long_bin = data.i4_start_bin + i4_long_bins;
	if (eoc_long_bin < eoc_bin) eoc_long_bin = eoc_bin;
	if (eoc_long_bin > bin_end) eoc_long_bin = bin_end; // no more than the gate
	
	  // Look for other clusters in this tail in which case close the cluster prematurelly
	for (; bin < eoc_long_bin; bin ++) { 	
	  number_hits = 0;
          for (int j = 0; j < 3; j++) number_hits += binned_times[bin + j]; // bin + j > i4_time_bins is fine since there is 20 spare bins
          if (number_hits > i4_start_threshold) { bin--; break; }
	  data.i4_n_hits += binned_times[bin];
	}

	  // Fix the end (integral is calulated in the previous loops)
	data.i4_long_end_bin = bin;
	data.i4_short_end_bin = (eoc_bin > bin) ? bin : eoc_bin;

//	get_message (bx_message::debug) << "found candidate cluster at " << data.i4_start_bin * 16 - f4_gate_start << ":" << data.i4_long_end_bin * 16 - f4_gate_start << "ns for event " << ev->get_event_number () << dispatch;

	
	  // Add this to the list of clusters
	clusters.push_back (data);
      }
        // else continue with the next bins
    }
  }
}

void bx_laben_findcluster::clusterize_mach4 (int start_bin__, int bin_end_, const bx_echidna_event *ev) {

	size_t start_bin_ = start_bin__;
	size_t bin_end = bin_end_;

	//assumes binsize == 16
	
	std::vector<uint16_t> ps; //Vector of hits - empty for now
	size_t psvlen = i4_time_bins + 20;  //use length of binned_times.
	
	ps.clear();
	ps.resize(psvlen, 0);
	
	//instead of filling hit times, just duplicate binned_times.
	
	for( size_t i = 0; i < psvlen; i++) ps[i] = binned_times[i];
	
	size_t start_gate_width = cluster_start_bins; //Number of bins to be used to determine start - set in Parameters.h
	size_t end_gate_width = cluster_end_bins; 
	
	//Calculate threshold  for start of cluster
	std::vector<double_t> start_threshold (psvlen, cluster_start_threshold); 
	//Threshold initialized to value in Parameters.h
	
	//Calculate threshold  for end of cluster
	//double dark_noise = 0.0134068; //hard code the value for run 10321 from Mach4		
	double dark_noise = f_dark_rate * 2000 * 16e-9;
	std::vector<double_t> end_threshold (psvlen, dark_noise + cluster_end_sigma*sqrt(dark_noise) ); 
	//Threshold initialized to values in Parameters.h
	
	size_t tail_min_width = cluster_tail_bins; //Minimum number of bins to add to end of cluster - set in Parameters.h
	size_t cluster_min_hits = cluster_integral; // Minimum hits required in each cluster
	size_t log_en = cluster_log_en; //Hit value above which a longer tail is added
	
	enum cluster_status_t { FILLING_START_GATE, SEARCHING_START,
	FILLING_END_GATE, SEARCHING_END }; //State machine - different stages of process
	
	cluster_status_t status = FILLING_START_GATE;
	size_t   bin, start_bin = 0, end_bin = 0, last_end_bin = 0, tail_bin;
	size_t   posn = 0;
	size_t   cluster_hits = 0, gate_hits = 0 /*Sum of hits in each set of start_gate_width bins*/;
	size_t   nclusters = 0;
	bool     store_cluster = false;
	
	//Loop over bins of histogram to identify clusters
	for (bin = start_bin_; bin < bin_end; bin++) {
		
		//#ifdef DEBUG_SPLITTING
		/*std::cout << "die" << std::endl;
		
		std::cout << "Bin " << std::setw(3) << bin << ", hits = "
		<< std::setw(4) << ps[bin] << ": ";
		
		std::cout << "here" << std::endl;*/
		//#endif
		
		switch (status) {
				
			case FILLING_START_GATE:
				gate_hits += ps[bin];	//Fill hits up to start_gate_width
				
				//#ifdef DEBUG_SPLITTING
				/*std::cout << "St Gate  ";
				std::cout << "Gate hits = " << std::setw(4) << gate_hits << "  ";
				std::cout << "Posn = " << std::setw(2) << posn << std::endl;*/
				//#endif
				
				posn++;
				
				if (posn == start_gate_width) { //If we have reached the end, start check for cluster start 
					posn = 0;
					if (gate_hits >= start_gate_width*start_threshold[(int)(bin-start_gate_width*0.5)]) { //Looks like the start of a cluster
						cluster_hits = gate_hits;
						gate_hits = 0;
						start_bin = bin - start_gate_width + 1;
						status = FILLING_END_GATE;
					}
					else
						status = SEARCHING_START;
				}
				break;
				
				case SEARCHING_START:
				//std::cout << "one" << std::endl;
				gate_hits += ps[bin] - ps[bin - start_gate_width];
				//std::cout << "two" << std::endl;
				
				if (gate_hits >= start_gate_width*start_threshold[(int)(bin-start_gate_width*0.5)]) { //Looks like the start of a cluster
					cluster_hits = gate_hits;
					gate_hits = 0;
					start_bin = bin - start_gate_width + 1;
					status = FILLING_END_GATE;
				}
				
				//std::cout << "three" << std::endl;
				
				//#ifdef DEBUG_SPLITTING
				/*std::cout << "Find St  ";
				if (status == FILLING_END_GATE)
					std::cout << "Gate hits = " << std::setw(4) << cluster_hits
					<< std::endl;
				else
					std::cout << "Gate hits = " << std::setw(4) << gate_hits
					<< std::endl;
				
				std::cout << "four" << std::endl;*/
				
				//#endif
				break;
				
				case FILLING_END_GATE:
				cluster_hits += ps[bin];
				gate_hits += ps[bin];	//Fill hits up to start_gate_width
				
				//#ifdef DEBUG_SPLITTING
				/*std::cout << "End Gate ";
				std::cout << "Gate hits = " << std::setw(4) << gate_hits << "  ";
				std::cout << "Posn = " << std::setw(2) << posn << "  ";
				std::cout << "Cluster hits = " << std::setw(4) << cluster_hits
				<< std::endl;*/
				//#endif
				posn++;
				if (posn == end_gate_width) { //If we have reached the end, start check for cluster end 
					posn = 0;
					status = SEARCHING_END;
				}
				break;
				
				case SEARCHING_END:
				gate_hits += ps[bin] - ps[bin - end_gate_width];
				if (gate_hits >=  end_gate_width*end_threshold[(int)(bin-end_gate_width*0.5)]) {
					cluster_hits += ps[bin];
				}
				else { //Looks like the end of the cluster
					end_bin = bin;
					gate_hits = 0;
					store_cluster = true;
					posn = 0;
					status = FILLING_START_GATE;
				}
				//#ifdef DEBUG_SPLITTING
				/*std::cout << "Find End ";
				std::cout << "Gate hits = " << std::setw(4) << gate_hits << "  ";
				std::cout << "           ";
				std::cout << "Cluster hits = " << std::setw(4) << cluster_hits
				<< std::endl;*/
				//#endif*/
				break;
		}
		
		//std::cout << "passed loop" << std::endl;
		
		// deal with the case where we get to the end of laben data while filling
		// a cluster:
		if (bin == bin_end - 1 &&
			(status == FILLING_END_GATE || status == SEARCHING_END)) {
			end_bin = bin_end; //psvlen;
			store_cluster = true;
		}
		
		//std::cout << "passed end case" << std::endl;
		
		if (store_cluster) {
			
			//std::cout << "start store cluster" << std::endl;
			
			if (cluster_hits >= cluster_min_hits) {// only actually store the cluster if it has enough hits
				if (end_bin >= bin_end)
					tail_bin = end_bin = bin_end;
				else {// include some "tail" period after the end gate threshold is reached
					double tail_width = tail_min_width;
					
					if (cluster_hits > log_en) // big event, increase tail length
						tail_width *= (1.0 + std::log(cluster_hits / (double)log_en));
					
					tail_bin = std::min(psvlen, end_bin + (int)(tail_width + 1)); //Ensure that tail isn't outside the window
					/*Note: We have only moved the end bin - we have not yet filled the hits between the tail_bin and end_bin.
					 That will only be done once the cluster is finalized - after the start of the next cluster is identified*/
				}
				
				// Generate the nth cluster if it doesn't already exist
				if (nclusters >= clusters.size() ) 
				{
					cluster_data data;
					clusters.push_back(data);
				}
				
				clusters[ nclusters ].i4_start_bin = start_bin;
				clusters[ nclusters ].i4_long_end_bin = tail_bin;
				clusters[ nclusters ].i4_n_hits = cluster_hits;
				
				// Add in hits in the tail of the *previous* cluster, including
				// truncating it if it overlaps the start of the *current* cluster
				if (nclusters) {
					
					if( clusters[nclusters-1].i4_long_end_bin > clusters[nclusters].i4_start_bin )
						clusters[nclusters-1].i4_long_end_bin = clusters[nclusters].i4_start_bin;
					
					/*Fill in the remaining hits between the tail_bin and the end_bin of the previous cluster*/
					for (size_t prev_bin = last_end_bin;	             //last_end_bin is the end_bin of the previous cluster		
						 prev_bin < size_t(clusters[nclusters-1].i4_long_end_bin); prev_bin++)
					{
						clusters[nclusters-1].i4_n_hits += ps[prev_bin];
					}
					
				}
				
				last_end_bin = end_bin;
				nclusters++;
				
			}
			
			// reset cluster-related variables whether or not cluster was actually
			// stored
		    //std::cout << "done store cluster" << std::endl;
			
			store_cluster = false;
			cluster_hits = 0;
		}
	}
	
	if (nclusters == 0) {
		//#ifdef DEBUG_SPLITTING
	//	Message(DEBUG) << "No valid laben cluster found in event "
	//	<< ev->EventNumber() << std::endl;
		//#endif
		
		return;
	}
	
	// Add in hits in tail of final cluster
	for (size_t prev_bin = last_end_bin;
		 prev_bin < size_t(clusters[nclusters-1].i4_long_end_bin); prev_bin++)
	{
		clusters[nclusters-1].i4_n_hits += ps[prev_bin];
	}
			
	
}


inline bool bigger (float a, float b, float min) { 
  return a > b + min;
}

inline bool smaller (float a, float b, float min) { 
  return a < b - min;
}

inline bool bigger_error (float a, float b) { 
  float err = ::sqrt (a + b);
  if (err < 2) err = 2;
  return bigger (a, b, err);
}

inline bool smaller_error (float a, float b) { 
  float err = ::sqrt (a + b);
  if (err < 2) err = 2;
  return smaller (a, b, err);
}

void bx_laben_findcluster::clusterize_neutrons (int start_bin, int bin_end, const bx_echidna_event *ev) {

    const bool mute = true;

    int num_neutron_clusters = 0;
    bool double_peak = false;
    float tailline_value = 0;

    float single_bin_crit_prob = 0.1;
    float tripple_bin_crit_prob = 1e-8;
    float minimum_mean_population = 0.01;

    float data_highest_bin_population = 0;
    bool data_start = false;
    int data_first_rising_bin = 0;
    int data_start_threshold = 1;
    for (int bin = start_bin; bin < bin_end; bin++) {
    	if (binned_times[bin] >= data_start_threshold && data_start == false){data_first_rising_bin = bin;data_start = true;}  //find where the data starts
	if (binned_times[bin] > data_highest_bin_population) data_highest_bin_population = binned_times[bin];
    }
    start_bin = data_first_rising_bin;

    TF1 poisson_func("poisson_func","TMath::Poisson(x,[0])",0,data_highest_bin_population);//intialize
    poisson_func.SetParameters(1,1);//intialize


    //Protection against misclustering of the gate-start
    //peak search starts after the highest bin which is within the first 160ns after gate-start
    int data_start_max_bin=0;
    int data_start_max_bin_pos=0;
    for (int bin = start_bin; bin < start_bin+10; bin++) {
    	if (binned_times[bin] > data_start_max_bin) {data_start_max_bin = binned_times[bin];data_start_max_bin_pos=bin;}
    }
    start_bin = data_start_max_bin_pos;

    if (!mute) get_message (bx_message::debug) << "Search for peaks starts at " << start_bin * 16 << "ns with entry " << binned_times[start_bin] << " in event " << ev->get_event_number () << dispatch;

    std::vector<float> last_cluster_starttimes_vector;
    std::vector<float> last_cluster_endtimes_vector;
    std::vector<float> noise_distribution_vector_bin;
    std::vector<float> noise_distribution_vector_bincontent;
    std::vector<float> rebinned_vector;
    last_cluster_starttimes_vector.clear();last_cluster_endtimes_vector.clear();
    noise_distribution_vector_bin.clear(); noise_distribution_vector_bincontent.clear();
    rebinned_vector.clear();


    for (int bin = start_bin; bin < bin_end; bin++) {
	    cluster_data data;

	    float baseline_value;
	    float probability=1;
	    float probability_first_bin=1;
	    float probability_second_bin=1;
            float bin_probability; 
	    float upper_integration_limit = data_highest_bin_population; //initialize
	    int integral;
	    float m=0;//initialize value
	    float b=minimum_mean_population;//initialize value

	    if (double_peak == false){

		    // First reject low integral regions (less than 5 hits)
		    integral = 0;
		    for (int j = 0; j < 3; j++) integral += binned_times[bin  + j];
		    if (integral < 5) continue;


		   //Fill the noise_distribution-vectors ignoring regions which have been identified as a cluster
		   noise_distribution_vector_bin.clear();
		   noise_distribution_vector_bincontent.clear();
		   int rejection_zone=5; //don't consider close 5bins = 80ns
		   int u = int(last_cluster_starttimes_vector.size()-1);
		   int nonzerocounter=0;
		   float range=-200;//intialize
 		   for (int j = 0; j >= start_bin-bin; j--) {
			if (j>-rejection_zone && j < rejection_zone) continue; //don't consider close 5bins = 80ns
			if (u >= 0) {
				if (bin+j >= last_cluster_starttimes_vector[u] && bin+j <= last_cluster_endtimes_vector[u]) continue;
				if (bin+j < last_cluster_starttimes_vector[u]) {u--;j++;continue;} //go j++ in case there is a double-cluster
			}
			noise_distribution_vector_bin.push_back(j);
			noise_distribution_vector_bincontent.push_back(binned_times[bin+j]);
			if (binned_times[bin+j] != 0) nonzerocounter++;
			range = j;
			if (nonzerocounter > 200) break;	//take 200 bins which are not empty to have enough statistics
		   }

		  //Fit a line to the noise_distribution using analytical calculation for linear regression
		  //the line is defined as : y = m*x +b
		  float mean_x=0;
		  float mean_y=0;
		  float mean_population;
		  float weight;
		  float weight_sum=0;
		  unsigned int n_bins = noise_distribution_vector_bin.size();

		  if (n_bins >= 2){//Fitting a line only possible for at least 2 points
		     for (unsigned int j = 0; j < n_bins;j++){
			weight = 1;//all weights set to 1; gives importance to strong statistical fluctuations
			mean_x += noise_distribution_vector_bin[j] * weight;
			mean_y += noise_distribution_vector_bincontent[j] * weight;
			weight_sum += weight;
		     }
		     mean_x /= weight_sum;
		     mean_y /= weight_sum;

		     float numerator=0;
		     float denominator=0;		   
		     for (unsigned int j = 0;j<n_bins;j++){
			weight = 1;//all weights set to 1; gives importance to strong statistical fluctuations
			numerator += (noise_distribution_vector_bin[j] - mean_x) * (noise_distribution_vector_bincontent[j] - mean_y) * weight;
			denominator += (noise_distribution_vector_bin[j] - mean_x) * (noise_distribution_vector_bin[j] - mean_x) * weight;
		     }

		     if (denominator != 0) m = numerator/denominator;
		     else m = 0;
		     b = mean_y - m * mean_x;

		     //the regression line is evaluated at x=j=0 which reduced the equation to y = b
		     if (fabs(b - mean_y) < 0.5) m = 0;//if b close to the mean_y -i.e. the inclination of the line is small- set inclination to zero
		     if (m <= 0) mean_population = mean_y; //for falling or constant noise_distribution always use the mean_value (more stable against fluctuations)
		     else mean_population = b;		//for rising noise always use extrapolated b, because mean_y would be too small
		     b = mean_y - m * mean_x;		//m can be set to zero, so b has to be recalculated
		  }
		  else if (n_bins == 1) mean_population = noise_distribution_vector_bincontent[0]; //if there is only one bin, use this bin as estimate of the noise
		  else {//this corresponds to n_bins==0
		   	if (bin - start_bin < rejection_zone) mean_population = data_start_max_bin; //at the gate start the first bins within the rejection zone are canceled. To get an estimate of the noise nevertheless, use data_start_max_bin
			else mean_population = minimum_mean_population;//in all other cases, use the minimum mean_population
		  }
		  if (mean_population < minimum_mean_population) mean_population = minimum_mean_population;
 

  		  //Now Calculate probability
		  poisson_func.SetParameter(0,mean_population);

		  probability=1;
		  rebinned_vector.clear();
		  for (int j=0; j<3; j++) {
	     		//take integral over the poisson-statistic
			//only bins contribute which are higher than the mean_population
			if (binned_times[bin+j] <= int(mean_population+0.5)) bin_probability=1;
			else {
				upper_integration_limit = 10 * binned_times[bin+j];
				if (upper_integration_limit < 20) upper_integration_limit = 20;
				bin_probability =  poisson_func.Integral(binned_times[bin+j],upper_integration_limit);
			}
			if (j==0) probability_first_bin = bin_probability;
			if (j==1) probability_second_bin = bin_probability;

			probability *= bin_probability;
			rebinned_vector.push_back(binned_times[bin+j]);		//remember the bins
		  }
		  if (probability_first_bin > single_bin_crit_prob) continue;	//no cluster if first bin is not unlikely enough

		  //Now check if there is a hit separation because of binning
		  if (probability_second_bin < probability_first_bin){	//if the bins are rising
			float bin_diff = rebinned_vector[0] - mean_population; //put the amount of the first bin which is above average on top of the second
			if (bin_diff < 0) bin_diff = 0;
			rebinned_vector[1] += bin_diff;
			rebinned_vector.erase(rebinned_vector.begin());	//bin_diff has to be subtracted from the first bin; would leave the first bin as an average bin; average bins don't contribute => remove it
			rebinned_vector.push_back(binned_times[bin+3]);	//another bin is added to have again 3 bins inside the vector

			//Calculate the probabilities anew
			probability = 1;
			for (int j=0;j<3;j++){
				if (rebinned_vector[j] <= int(mean_population+0.5)) bin_probability=1;
				else {
					upper_integration_limit = 10 * rebinned_vector[j];
					if (upper_integration_limit < 20) upper_integration_limit = 20;
					bin_probability =  poisson_func.Integral(rebinned_vector[j],upper_integration_limit); //higher integral upper limit because 2 bins were summed
				}
				if (j==0) probability_first_bin = bin_probability;
				if (j==1) probability_second_bin = bin_probability;

				probability *= bin_probability;
			}
		 }
	

		 //If the 3 bins in a row are statistically unlikely then a cluster is detected
		 if (!( probability < tripple_bin_crit_prob )) continue;
		 baseline_value = mean_population;

	         if (!mute) get_message (bx_message::debug) << bin * 16 << " ns; Probability " << probability << dispatch;
		 if (!mute) get_message (bx_message::debug) << bin * 16 << " ns; first bin probability " << probability_first_bin << dispatch;
		 if (!mute) get_message (bx_message::debug) << bin * 16 << " ns; second bin probability " << probability_second_bin << dispatch;
	         if (!mute) get_message (bx_message::debug) << bin * 16 << " ns; mean " << mean_population << dispatch;
	         if (!mute) get_message (bx_message::debug) << bin * 16 << " ns; m " << m << dispatch;
 		 if (!mute) get_message (bx_message::debug) << bin * 16 << " ns; range " << range << dispatch;
		 if (!mute) for (int j=-100;j<100;j++){ get_message (bx_message::debug) << bin * 16 << " ns; j is " << j << " ; in ns " << (bin+j) * 16 << " : " << binned_times[bin+j] << dispatch;}
 		if (!mute) for (int j=0;j<3;j++){ get_message (bx_message::debug) << bin * 16 << " ns; vec j is " << j << " : " << rebinned_vector[j] << dispatch;}

   }
   else {baseline_value = b = tailline_value; m=0; double_peak = false;} //this line comes into play for double_peak search; the tailline is from the previos cluster; for calculation of the background-hits inside the cluster the baseline is assumed to be constant


    //Now the cluster-start has been identified
    data.i4_n_hits = 0;
    data.i4_start_bin = bin;

    if (!mute) get_message (bx_message::debug) << bin * 16 << " ns; rising in bin " << bin << " with " << binned_times[bin] << " hits." << dispatch;
    if (!mute) get_message (bx_message::debug) << "baseline " << baseline_value << dispatch;

    // Finally search for falling edge
    double peak_maximum = binned_times[data.i4_start_bin];
    bool falling_edge = false;
    int stop = bin + 10;
    double falling_edge_threshold;
    if (stop > bin_end) stop = bin_end;
    for (; bin < stop; bin++) {
      if (binned_times[bin] > peak_maximum) peak_maximum = binned_times[bin];
      falling_edge_threshold = 0.5 * (peak_maximum + baseline_value);	//falling edge defined as 50% of peak maximum with respect to initial baseline
      if (!bigger(binned_times[bin], falling_edge_threshold, 0.) && !bigger(binned_times[bin+1], falling_edge_threshold, 0.))
      { falling_edge = true; break;} //falling edge detected if 2 consecutive bins are below falling edge threshold
      data.i4_n_hits += binned_times[bin];
    }

    if (!falling_edge) { if (!mute) get_message (bx_message::debug) << bin * 16 << " ns; falling edge too late" << dispatch; bin = data.i4_start_bin; continue;};
    if (!mute) get_message (bx_message::debug) << bin * 16 << " ns; falling in bin " << bin << dispatch;
    int falling_bin = bin;


    // Peak found: now identify the end and also check if another peak is rising
    bool end = false;
    stop = bin + 25;
    if (stop > bin_end) stop = bin_end;
    for (; bin < stop; bin++) {

     //Check if another peak is rising; analog to peak finding routine from before
     float mean_x = 0;
     float mean_y = 0;
     float weight;
     float weight_sum = 0;
//     float numerator=0;
//     float denominator=0;
//     float m_tail,b_tail;
     bool cluster_rejected;
      if (bin - falling_bin >= 4){	//bin-distance necessary for tailline search

        tailline_value = 0;

        mean_x=0;mean_y=0;weight_sum=0;
        for (int j = -4; j < 0; j++){
	  weight = 1;//all weights set to 1; gives importance to strong statistical fluctuations
          mean_x += j * weight;
          mean_y += binned_times[bin + j] * weight;
          weight_sum += weight;
        }
	mean_x /= weight_sum;
	mean_y /= weight_sum;

/*
        numerator=0;denominator=0;
        for (int j = -5;j<0;j++){
	  weight = 1;//all weights set to 1; gives importance to strong statistical fluctuations
	  numerator += (j - mean_x) * (binned_times[bin + j] - mean_y) * weight;
	  denominator += (j - mean_x) * (j - mean_x) * weight;
	}
	if (denominator != 0) m_tail = numerator/denominator;
	else m_tail = 0;
        b_tail = mean_y - m_tail * mean_x;


        if (fabs(b_tail - mean_y) < 0.5) m_tail = 0;//if b_tail close to the mean_y -i.e. the inclination of the line is small- set inclination to zero
        if (m_tail <= 0) tailline_value = mean_y; //for falling or constant noise_ditstribution always use the mean_value (more stable against fluctuations)
        else tailline_value = b_tail; 
	b_tail = mean_y - m_tail * mean_x;       //m_tail can be set to zero, so b_tail has to be recalculated
*/

	//For the double_peak search, use only the mean for tailline determination (more stable for very close peaks)
	tailline_value = mean_y;

	//Now Calculate probability
	poisson_func.SetParameter(0,tailline_value);

	probability = 1;
	rebinned_vector.clear();
        for (int j=0; j<3; j++) {
	  if (binned_times[bin+j] <= int(tailline_value+0.5)) bin_probability=1;
	  else {
                upper_integration_limit = 10 * binned_times[bin+j];
		if (upper_integration_limit < 20) upper_integration_limit = 20;
	  	bin_probability =  poisson_func.Integral(binned_times[bin+j],upper_integration_limit);
	  }
	  if (j==0) probability_first_bin = bin_probability;
	  if (j==1) probability_second_bin = bin_probability;
	
	  probability *= bin_probability;
	  rebinned_vector.push_back(binned_times[bin+j]);		//remember the bins
	}
	cluster_rejected = false;
	if (probability_first_bin > single_bin_crit_prob) cluster_rejected = true;	//no cluster if first bin is not unlikely enough or only one bin is above average

	//Now check if there is a hit separation because of binning
	if (probability_second_bin < probability_first_bin){
		float bin_diff = rebinned_vector[0] - tailline_value;
		if (bin_diff < 0) bin_diff = 0;
		rebinned_vector[1] += bin_diff;
		rebinned_vector.erase(rebinned_vector.begin());	
		rebinned_vector.push_back(binned_times[bin+3]);

		//Calculate the probabilities anew
		probability = 1;
		for (int j=0;j<3;j++){
			if (rebinned_vector[j] <= int(tailline_value+0.5)) bin_probability=1;
			else {
		                upper_integration_limit = 10 * rebinned_vector[j];
		                if (upper_integration_limit < 20) upper_integration_limit = 20;
				bin_probability =  poisson_func.Integral(rebinned_vector[j],upper_integration_limit);
			     }
			if (j==0) probability_first_bin = bin_probability;
			if (j==1) probability_second_bin = bin_probability;

			probability *= bin_probability;
		}
	}


	if (probability < tripple_bin_crit_prob && cluster_rejected == false) {
          if (!mute) get_message (bx_message::debug) << bin * 16 << " ns; double peak detected in bin " << bin << " with " << binned_times[bin] << " hits." << dispatch;
          if (!mute) get_message (bx_message::debug) << "tailline " << tailline_value << dispatch;

	  double_peak = true;
          data.i4_long_end_bin = bin;
	  bin--;
          break;
        }
      }

      //Identify endline and only accept it if it's not noisy
      integral = 0;
      for (int j = 0; j < 10; j++) integral += binned_times[bin + j];
      float end_line = integral / 10.;

      if (!bigger_error (end_line, baseline_value)) {
	if (!mute) get_message (bx_message::debug) << bin * 16 << " ns; end found with end_line " << end_line << " vs base_line " << baseline_value << dispatch;
	data.i4_n_hits += integral;
	end = true;
	data.i4_long_end_bin = bin = bin + 10;
	bin--;
	break;
      }

      data.i4_n_hits += binned_times[bin];
    }

    //Calculate the Number of Backgroundhits during the clusterlength
    float cluster_length = float(data.i4_long_end_bin - data.i4_start_bin);
    data.f4_n_hits_background = 0.5 * m * cluster_length * cluster_length + b * cluster_length;
    if (data.f4_n_hits_background < 0) data.f4_n_hits_background = 0;


    if (data.i4_n_hits < 5.) {
      if (!mute) get_message (bx_message::debug) << "low integral neutron cluster rejected (" << data.i4_n_hits << " hits) for event " << ev->get_trigger ().get_evid () << dispatch;
      bin = data.i4_start_bin;
      continue; // IGNORE this cluster if integral is too low
    }

    if (!end && !double_peak) {
      if (!mute) get_message (bx_message::debug) << "endline not found; cluster rejected." << dispatch;
      bin = data.i4_start_bin;
      continue; 
    }


    //Remember detected cluster
    last_cluster_starttimes_vector.push_back(data.i4_start_bin);
    last_cluster_endtimes_vector.push_back(data.i4_long_end_bin);
    num_neutron_clusters++;


    // Add this to the list of clusters
    if (!mute) get_message (bx_message::debug) << "neutron cluster found at " << data.i4_start_bin * 16 << " ns; length : " << cluster_length* 16 << " ns with " << data.i4_n_hits << " hits and " << data.f4_n_hits_background << " background hits for evnum " <<  ev->get_event_number () << dispatch;

    //Triggertype128 was clustered
    data.bool_muon_clustered = 0;

    clusters.push_back (data);
  }

   if (num_neutron_clusters != 0 && !mute) get_message (bx_message::debug) << num_neutron_clusters << " neutron clusters found in event " << ev->get_event_number () << dispatch;

}




void bx_laben_findcluster::clusterize_neutrons_in_muongate (int start_bin, int bin_end, const bx_echidna_event *ev) {
     //clusterize_neutrons can't handle the ionic afterpulses in the muon gate because the rising flanks in the time profile are on too small time scales. The poisson-statistics assumption breaks down.
     //clusterize_neutrons_in_muongate fixes this problem by leaving the assumption of poisson-statistics and calculating the mean and rms of the baseline. In addition, it's less sensitive than clusterize_neutrons.
     //both conditions avoid the identification of noise clusters


     int num_neutron_clusters = 0;
     bool double_peak = false;
     float empty_board_factor = 1. - ev->get_laben().get_empty_boards() / 280.;
     float pmt_factor = (ev->get_laben ().get_n_live_pmts() - ev->get_laben ().get_invalid_pmts() ) / 2000.;

     const bool mute = true;

     if (!mute) get_message (bx_message::debug) << "number of empty boards : " << ev->get_laben().get_empty_boards() << " in event " << ev->get_event_number () << dispatch;


    // Use derivative algorithm
    if (start_bin < 11) start_bin = 11;  // Protect for baseline loop which search on bin - 11, bin
    for (int bin = start_bin; bin < bin_end; bin++) {
    cluster_data data;


    // First reject low integral regions
    int integral = 0;
    for (int j = 0; j < 5; j++) integral += binned_times[bin  + j];
    if (integral < 20. * empty_board_factor * pmt_factor) continue;
    if (!mute) {
      bx_message &msg = get_message (bx_message::debug);
      msg << bin * 16 << " ns integral found with " << integral << " ";
      for (int j = 0; j < 25; j++) msg << binned_times[bin  + j] << " ";
      msg << " evnum " <<  ev->get_event_number () << dispatch;
    }


    int baseline_range;
    if (!double_peak) baseline_range = -11;
    else baseline_range = -6;

    // Then search for baseline
    bool baseline = true;
    float baseline_value = 0;
    for (int j = baseline_range; j < -1; j++) { 
      baseline_value += binned_times[bin + j];
      if (bigger_error (binned_times[bin + j], binned_times[bin + j - 1])) baseline = false;
    }
    baseline_value /=  -baseline_range - 1.;  //just mean value; (-baseline_range - 1) is 10 for baseline_range = -11 and 5 for -6

    float baseline_rms = 0;
    for (int j = baseline_range; j < -1; j++) { 
      baseline_rms += (binned_times[bin + j] - baseline_value) * (binned_times[bin + j] - baseline_value);
    }
    baseline_rms = sqrt(baseline_rms / ( -baseline_range - 1. ) );


    // Then search for rising edge
    data.i4_n_hits = 0;     

    if (!bigger (binned_times[bin], baseline_value, baseline_value + 3 * baseline_rms +  10. * empty_board_factor * pmt_factor)) continue;
    data.i4_start_bin = bin;
    data.i4_n_hits += binned_times[bin];
    if (!mute) get_message (bx_message::debug) << bin * 16 << " ns; rising in bin " << bin << " with " << data.i4_n_hits << " hits." << dispatch;
    bin += 1;

    if (!mute) get_message (bx_message::debug) << "baseline " << baseline_value << " ; baseline_rms " << baseline_rms << dispatch;


    // Finally search for falling edge
    double peak_maximum = binned_times[data.i4_start_bin];
    bool falling_edge = false;
    int stop = bin + 10;
    if (stop > bin_end) stop = bin_end;
    for (; bin < stop; bin++) {
      if (binned_times[bin] > peak_maximum) peak_maximum = binned_times[bin];
      double falling_edge_threshold = 0.5 * (peak_maximum + baseline_value);	//falling edge defined as 50% of peak maximum with respect to initial baseline
      if (!bigger(binned_times[bin], falling_edge_threshold, 0.) && !bigger(binned_times[bin+1], falling_edge_threshold, 0.))
      { falling_edge = true; break;} //falling edge detected if 2 consecutive bins are below falling edge threshold
      data.i4_n_hits += binned_times[bin];
    }

    if (!falling_edge) { if (!mute) get_message (bx_message::debug) << bin * 16 << " ns; falling edge too late" << dispatch; bin = data.i4_start_bin; continue;};
    if (!mute) get_message (bx_message::debug) << bin * 16 << " ns; falling in bin " << bin << dispatch;
    int falling_bin = bin;


     // Peak found: now identify the end and also check if another peak is rising
    bool end = false;
    double_peak = false;
    stop = bin + 25;
    if (stop > bin_end) stop = bin_end;
    for (; bin < stop; bin++) {

      //Check if another peak is rising; analog to peak finding routine from before
      if (bin - falling_bin >= 6){	//bin-distance necessary for tailline search

        float tailline_value = 0;
        for (int j = -5; j < -1; j++) {
          tailline_value += binned_times[bin + j];
        }
        tailline_value /= 4.;

        float tailline_rms = 0;
        for (int j = -5; j < -1; j++) {
          tailline_rms += (binned_times[bin + j] - tailline_value) * (binned_times[bin + j] - tailline_value);
        }
        tailline_rms = sqrt(tailline_rms / 4.);

        if (bigger (binned_times[bin], tailline_value, tailline_value + 3 * tailline_rms +  10. * empty_board_factor * pmt_factor)){
          if (!mute) get_message (bx_message::debug) << bin * 16 << " ns; double peak detected in bin " << bin << " with " << binned_times[bin] << " hits." << dispatch;
          if (!mute) get_message (bx_message::debug) << "tailline " << tailline_value << " ; tailline_rms " << tailline_rms << dispatch;
	  double_peak = true;
          data.i4_long_end_bin = bin;
	  bin--;
          break;
        }
      }


      //Identify endline 
      integral = 0;
      for (int j = 0; j < 10; j++) integral += binned_times[bin + j];
      float end_line = integral / 10.;

      if (!bigger_error (end_line, baseline_value)) {
	if (!mute) get_message (bx_message::debug) << bin * 16 << " ns; end found with end_line " << end_line << " vs base_line " << baseline_value << dispatch;
	data.i4_n_hits += integral;
	end = true;
	data.i4_long_end_bin = bin = bin + 10;
	bin--;
	break;
      }

      data.i4_n_hits += binned_times[bin];
    }

    //Calculate the Number of Backgroundhits during the clusterlength
    data.f4_n_hits_background = float(data.i4_long_end_bin - data.i4_start_bin) * baseline_value;


    if (data.i4_n_hits < 40. * empty_board_factor * pmt_factor || data.i4_n_hits < 5.) {
      if (!mute) get_message (bx_message::debug) << "low integral neutron cluster rejected (" << data.i4_n_hits << " hits) for event " << ev->get_trigger ().get_evid () << dispatch;
      continue; // IGNORE this cluster if integral is too low
    }

    if (!end && !double_peak) {
      if (!mute) get_message (bx_message::debug) << "endline not found; cluster rejected." << dispatch;
      bin = falling_bin;
      continue; 
    }

    num_neutron_clusters++;


      // Add this to the list of clusters
    if (!mute) get_message (bx_message::debug) << "neutron cluster found at " << data.i4_start_bin * 16 << " ns; length : " << (data.i4_long_end_bin - data.i4_start_bin) * 16 << " ns with " << data.i4_n_hits << " hits and " << data.f4_n_hits_background << " background hits for evnum " <<  ev->get_event_number () << dispatch;

    //The muon gate was clustered with the algorithm
    data.bool_muon_clustered = 1;
    clusters.push_back (data);
  }

   if (num_neutron_clusters != 0 && mute == false) get_message (bx_message:: debug) << num_neutron_clusters << " neutron clusters found in event " << ev->get_event_number () << dispatch;

}
//! Get the index within the gate, assuming a given bin size for the time
int bx_laben_findcluster::get_gate_index_ (int gate_length,Float_t dT) {
  int N=(int)floor(dT/gate_length);
  return (N); 
}


/*findcluster_time fills i4_npmts and other variables. 
  It finds the start time of each cluster*/
void bx_laben_findcluster::findcluster_time (bx_echidna_event *ev) {
  bx_laben_event& er = ev->get_laben ();
  bx_laben_clustered_event &ew = dynamic_cast<bx_laben_clustered_event&>(ev->get_laben ());
  for (int cl = 0; cl < er.get_nclusters () + er.get_nclusters_muons(); cl++) {
      // Keep the cluster reference
    bx_laben_cluster::bx_laben_cluster_vector *cluster_vector;
    int index = cl;
    if (cl < er.get_nclusters ()) cluster_vector = &ew.clusters;
    else {
      cluster_vector = &ew.clusters_muons;
      index = cl - er.get_nclusters ();
    }
    bx_laben_cluster &cluster_ref = (*cluster_vector)[index];

      // Calculate the starting time separation
    float time_separation = 3;
    int intervals = 5;
      // special handling for low energy and high energy
    if (cluster_ref.get_clustered_nhits () > 500) {
      time_separation = 2;
    } else if (cluster_ref.get_clustered_nhits () < 40) {
      time_separation = 6;
      intervals = 3;
    } else if (cluster_ref.get_clustered_nhits () < 70) {
      time_separation = 5;
      intervals = 3;
    } else if (cluster_ref.get_clustered_nhits () < 150) {
      time_separation = 4;
      intervals = 4;
    }

      // Serch the first it. Pay attention find_start_hit is a recursive function which even modify the second argument
    int start_hit = find_start_hit (cluster_ref, time_separation, intervals);
    if (start_hit < 0) start_hit = find_start_hit (cluster_ref, time_separation, --intervals);
    if (start_hit < 0) {
      start_hit = 0;
      cluster_ref.u1_flag |= bx_laben_cluster::broad;
      get_message (bx_message::log) << "start time of cluster " << cl + 1 << " not found for event " << ev->get_trigger ().get_evid () 
	<< " ( " << cluster_ref.get_clustered_hit (0).get_time () - er.get_trigger_rawt () << " ) " << dispatch;
    } 
    //else 
//      get_message (bx_message::debug) << "start time of cluster " << cluster + 1 << " found for event " << ev->get_trigger ().get_evid ()
//        << " at hit n " << start_hit << " with time_separation of " << time_separation << "ns and intervals " << intervals << dispatch;
   
   //To make it identical to Mach4 clustering:
   if( b_mach4_cluster ) start_hit = 0;

   //for low-nhits tt128 and tt1 neutrons use really the first hit
   if ((ev->get_trigger ().is_neutron () || cl >= er.get_nclusters ()) && cluster_ref.get_clustered_nhits() < 200.) start_hit = 0;

      // Setting of the cluster start time and the rough time
    cluster_ref.f8_start_time = cluster_ref.get_clustered_hit(start_hit).get_time();
//    cluster_ref.f8_rough_time = cluster_ref.get_clustered_hit(0).get_time();
    double relative_time = cluster_ref.f8_start_time - er.get_trigger_rawt ();
    if (relative_time < i4_gate_start || relative_time > i4_gate_end) cluster_ref.u1_flag |= bx_laben_cluster::out_of_gate;
    if (::fabs (relative_time - i4_cluster_offset) < 250) cluster_ref.u1_flag |= bx_laben_cluster::trigger;

      // Some photons might be discarded from the hit list (rely on the strict time ordering).
    if (start_hit) cluster_ref.clustered_hits.erase (cluster_ref.clustered_hits.begin(), cluster_ref.clustered_hits.begin() + start_hit);

      // Setting of the photons' times (relative to the time of the first photon)
    std::vector<bx_laben_clustered_hit>::iterator it = cluster_ref.clustered_hits.begin ();
    for (; it != cluster_ref.clustered_hits.end (); it++) {
        double dt = it->get_time () - cluster_ref.get_start_time ();
	it->f8_time = dt;
      }

      // Warning message
    if (cluster_ref.get_clustered_nhits () == 0) {
      get_message (bx_message::info) << "Cluster " << cl + 1 << " of event " << ev->get_trigger ().get_evid () << " is empty!" << dispatch;
      continue;
    }

      // Setting of the other variables of the cluster
    std::fill_n (fired_channel, constants::laben::channels, 0);


//*****************************************************noavg code***********************************************************
      double dark_rate = 0;
//***********************************************************************************************************************
    double sum_square = 0;
    double sum_square_short = 0;
    for (it = cluster_ref.clustered_hits.begin(); it != cluster_ref.clustered_hits.end (); it++) {
      const bx_laben_decoded_hit& dhit = it->get_decoded_hit ();
      const db_channel_laben* db = dhit.get_db_channel ();
      int lg = db->get_lg ();
      bool cone = db->pmt_has_cone();
      sum_square += it->f8_time * it->f8_time;
      cluster_ref.f4_mean_time += it->f8_time;
      cluster_ref.f4_charge += dhit.get_charge_pe ();
      cluster_ref.f4_charge_mean += dhit.get_charge_mean_pe ();
      cluster_ref.i4_npe += dhit.get_charge_npe ();
//*****************************************************navg code***********************************************************
      double npe = dhit.get_charge_pe ();
      if(npe<0.1) npe=1; if(npe>3) npe=0;//noavg correction
       cluster_ref.f4_charge_noavg+=npe;
       dark_rate += run_info->get_laben_dark_noise(lg);//cumulative for DH probability in 1/10^3s
//****************************************************************************************************************************/
      if (cone) {
	cluster_ref.i4_clustered_nhits_conc++;
        cluster_ref.f4_charge_conc += dhit.get_charge_pe  ();
        cluster_ref.i4_npe_conc    += dhit.get_charge_npe ();
      }	
      if (it->f8_time < 400) {
	cluster_ref.f4_charge_400 += dhit.get_charge_pe();
	cluster_ref.i4_clustered_nhits_400 ++;
	if (dhit.get_charge_pe () >= 0.2) {
	  cluster_ref.i4_clustered_nhits_thresh ++;
	  cluster_ref.f4_charge_thresh += dhit.get_charge_pe ();
	}
      }
      if (!fired_channel[lg-1]) {
	cluster_ref.i4_npmts ++;
	cluster_ref.f4_charge_npmts += dhit.get_charge_pe ();
	if (cone) cluster_ref.i4_npmts_conc ++;
	if (it->b_short_cluster) cluster_ref.i4_npmts_short ++;
	if (it->f8_time < 400) {
	  cluster_ref.i4_npmts_400 ++;
	  if (dhit.get_charge_pe () >= 0.2) {
	    cluster_ref.i4_npmts_thresh ++;
	    if (dhit.get_charge_pe () < 5) cluster_ref.f4_charge_clean += dhit.get_charge_pe ();
	  }
	}
      }
      it->u1_order_in_channel = ++fired_channel[lg-1];
      if (it->b_short_cluster) {
	cluster_ref.i4_clustered_nhits_short++;
	cluster_ref.f4_mean_time_short += it->f8_time;
	cluster_ref.f4_duration_short = (it->f8_time > cluster_ref.f4_duration_short) ? it->f8_time : cluster_ref.f4_duration_short;
	sum_square_short += it->f8_time * it->f8_time;
	cluster_ref.f4_charge_short += dhit.get_charge_pe ();
//*****************************************************noavg code***********************************************************
            double npe = dhit.get_charge_pe ();
           if(npe<0.1) npe=1; if(npe>3) npe=0;//noavg correction
           cluster_ref.f4_charge_noavg_short+=npe;
//****************************************************************************************************************************/
      }
    } //End loop over clustered hits

    cluster_ref.f4_mean_time /= cluster_ref.get_clustered_nhits ();
    cluster_ref.f4_rms_time = sqrt (sum_square / cluster_ref.get_clustered_nhits () - cluster_ref.f4_mean_time * cluster_ref.f4_mean_time); 
    if (cluster_ref.i4_clustered_nhits_short) {
      cluster_ref.f4_mean_time_short /= cluster_ref.i4_clustered_nhits_short;
      cluster_ref.f4_rms_time_short = sqrt (sum_square_short / cluster_ref.i4_clustered_nhits_short - cluster_ref.f4_mean_time_short * cluster_ref.f4_mean_time_short);
    }

    //Loop over decoded hits to find fixed-time-window npmts variables
    if(cl<0) continue;   //Do it for all clusters
    std::fill_n(fired_channel,constants::laben::channels,0);
    /*Since decoded hits are time-ordered, 
      we count the earliest hit that takes place within each PMT, 
      after the start of the 1st cluster. 
      So we need only 1 fired_channel array for 2 fixed-time npmts variables.*/
    for ( int hit = 0; hit < er.get_decoded_nhits (); hit++)
    { //Iterate over all decoded hits in the bx_laben_event
      const bx_laben_decoded_hit &decoded_hit = er.get_decoded_hit (hit);
      const db_channel_laben* db = decoded_hit.get_db_channel ();
      /*Following line skips invalid triggers; copied L209 in negated version*/
      if (decoded_hit.get_flag () & ~bx_laben_decoded_hit::out_of_gate || 
	  !(db->is_ordinary ()))
	continue;
      const double hit_time = decoded_hit.get_raw_time ();
      const double cl_start_t = cluster_ref.get_start_time ();
      const double rel_time = hit_time - cl_start_t;
      const int lg = db->get_lg ();
      if(rel_time >= 0 && !fired_channel[lg-1]) 
      {
	fired_channel[lg-1]=1;
	if(rel_time < i4_dt1_len_) {
	  ++cluster_ref.i4_npmts_dt1;
	  cluster_ref.f4_charge_dt1 += decoded_hit.get_charge_pe ();
//*****************************************************noavg code***********************************************************
           double npe = decoded_hit.get_charge_pe ();
           if(npe<0.1) npe=1; if(npe>3) npe=0;//noavg correction
           cluster_ref.f4_charge_noavg_dt1+=npe;
//***********************************************************************************************************************
	  }
	if(rel_time < i4_dt2_len_) {
	  ++cluster_ref.i4_npmts_dt2;
	  cluster_ref.f4_charge_dt2 += decoded_hit.get_charge_pe ();
//*****************************************************noavg code***********************************************************
           double npe = decoded_hit.get_charge_pe ();
           if(npe<0.1) npe=1; if(npe>3) npe=0;//noavg correction
           cluster_ref.f4_charge_noavg_dt2+=npe;
//***********************************************************************************************************************
          }
      }
    } //End loop over decoded hits
//*****************************************************noavg code***********************************************************
    cluster_ref.f4_charge_noavg_dt1-=dark_rate*i4_dt1_len_*1E-6;//dark rate compensation for math. expectation in window
    cluster_ref.f4_charge_noavg_dt1-=dark_rate*i4_dt2_len_*1E-6;//dark rate compensation for math. expectation in window
    cluster_ref.f4_charge_noavg_short-=dark_rate*cluster_ref.get_duration_short()*1E-6;//dark rate compensation for math. expectation in window
    cluster_ref.f4_charge_noavg-=dark_rate*cluster_ref.get_duration()*1E-6;//dark rate compensation for math. expectation in window
//*****************************************************************************************************************************
  } //End loop over clusters

  if (!(ev->get_trigger ().is_neutron ()) && er.get_nclusters() != 0 ){
  //tt128 neutrons are the only events for which the clustering algorithm (neutron clustering) fills laben.clusters.nhits_bkg
  //for all other events the background is determined as follows:

        int cluster_num = 0;
        double bkg_time_window;   //in ns
        if (ev->get_run_number () <= 12421) bkg_time_window = 350;
	else bkg_time_window = 50;
	double effective_bkg_time_window = bkg_time_window;
        double end_time_of_last_cluster = -1e10;  //because for the first cluster there is no previous cluster, this number has to be reasonable negative so that the whole gate start is included

        bx_laben_cluster::bx_laben_cluster_vector *cluster_vector;
        cluster_vector = &ew.clusters;
        bx_laben_cluster &cluster_ref = (*cluster_vector)[cluster_num];

	// Number of good decoded hits
	double hits_before_cluster = 0;
	for (int i = 0; i < er.get_decoded_nhits (); i++) {
		const bx_laben_decoded_hit &decoded_hit = er.get_decoded_hit (i);
		
		if ((decoded_hit.get_flag () & ~bx_laben_decoded_hit::out_of_gate) && decoded_hit.get_db_channel ()->is_ordinary ()) continue;

		//select decoded hits within the bkg_time_window before the cluster start
		if (decoded_hit.get_raw_time () < cluster_ref.get_start_time () && decoded_hit.get_raw_time () > cluster_ref.get_start_time () - effective_bkg_time_window ) hits_before_cluster++;

                if (decoded_hit.get_raw_time () >= cluster_ref.get_start_time () ) {  //decoded hits arrive at the cluster start

			//now the number of background hits in the long cluster are estimated
			cluster_ref.f4_clustered_nhits_bkg = float(hits_before_cluster / effective_bkg_time_window * cluster_ref.get_duration());

			if (cluster_num == er.get_nclusters () - 1) break;               //last cluster checked
			end_time_of_last_cluster = cluster_ref.get_start_time() + cluster_ref.get_duration();
		
			cluster_num++;  //otherwise take next cluster
			bx_laben_cluster::bx_laben_cluster_vector *cluster_vector;
			cluster_vector = &ew.clusters;
			bx_laben_cluster &cluster_ref = (*cluster_vector)[cluster_num];
		
			hits_before_cluster = 0;
		
			if (cluster_ref.get_start_time () - end_time_of_last_cluster < bkg_time_window) effective_bkg_time_window = cluster_ref.get_start_time () - end_time_of_last_cluster;
			else effective_bkg_time_window = bkg_time_window;
                }
	}
  }

}

int bx_laben_findcluster::find_start_hit (const bx_laben_cluster& cluster_ref, float& time_separation, int intervals) {
  int start_hit = -1;
  
  for (int hit = 0; hit < cluster_ref.get_clustered_nhits () - intervals - 1; hit++) {
      // Only check for the first 64ns (rely on the strict time ordering of the hits)
    if (cluster_ref.get_clustered_hit (hit).get_time () - cluster_ref.get_clustered_hit (0).get_time () > double (64)) break;
    
    bool good = true;
    for (int i = 0; i < intervals; i++) {
      double dt = cluster_ref.get_clustered_hit (hit + i + 1).get_time () - cluster_ref.get_clustered_hit (hit + i).get_time (); 
      if (dt > time_separation) {
	good = false;
	break;
      }
    }
    if (good) {
      start_hit = hit;
      break;
    }
  }

  if (start_hit == -1) {
    if (time_separation > 16) return -1;
    else return find_start_hit (cluster_ref, ++time_separation, intervals);
  } else return start_hit;
}

int bx_laben_findcluster::evaluate_ripple_count (float dark_rate, int ripple_bins, int count_threshold, const char* message) {
    // Calculate the ripple threshold as number of hits in ripple_bins for which P(>=ripple_count) < 0.1% for dark rate
  double mu = dark_rate * 2000 * 16e-9 * ripple_bins;
  double p0 = ::exp (-1 * mu); 
  double pi = p0;
  double p_higher = 1 - p0;
  int ripple_count;
  for (ripple_count = 1; ripple_count < count_threshold; ripple_count++) { // never bigger than threshold
    pi = pi * mu / ripple_count; 
    p_higher -= pi; 
    if (p_higher < 1E-3) break;
  }
  if (p_higher > 1E-3) get_message (bx_message::error) << "too few ripple_bins " << ripple_bins << " cause of too high dark rate " << dark_rate << " for " << message << dispatch;
  else get_message (bx_message::info) << "ripple window set to " << ripple_bins << " bins with count " << ripple_count << "hits and fake probability " << p_higher << " for " << message << dispatch; 
  return ripple_count;
}

/*
 * $Log: bx_laben_findcluster.cc,v $
 * Revision 1.134  2015/07/20 15:23:37  ilia.drachnev
 * noavg charges added
 *
 * Revision 1.133  2014/12/15 12:08:24  dekor
 * Added charge_win1 charge_win2 vars
 *
 * Revision 1.132  2014/12/11 21:27:12  wurm
 * added getters and variables to BxEvent that concern C-11 likelihood tagging
 * * for tracks: new distance-to- and projection-on-track getters, dE/dx calculation
 * * counting of neutron-like clusters in tt128, new flag of neutron clusters
 *
 * Revision 1.131  2014/12/05 15:40:38  smirnov
 * dhit -> decoded_hit
 *
 * Revision 1.130  2014/12/05 12:07:47  smirnov
 *
 * VS: ----------------------------------------------------------------------
 * 05/12/2014 by Oleg S.
 * added filling of charge_dt1 and charge_dt2 variables
 *
 * Revision 1.129  2013/06/05 15:15:06  mosteiro
 * fill npmts_winX variables for neutrinos; shift window by random number
 *
 * The npmts_win1 and npmts_win2 variables were previously for tt64
 * Now we fill them also for neutrinos
 * In case of neutrinos, we shift the start
 * by a random number in the range of 0 to dt1/dt2, respectively
 *
 * Revision 1.128  2013-01-25 16:28:06  mosteiro
 * channel count starts at 1, not 0
 *
 * Revision 1.127  2013-01-25 10:56:23  mosteiro
 * warnings in array variable creation changed to debug messages
 *
 * Revision 1.126  2013-01-22 09:22:17  mosteiro
 * dt1,2 vars calculated for each cluster
 *
 * Revision 1.125  2012-11-25 22:11:03  mosteiro
 * changed two error messages to warn status
 *
 * Revision 1.124  2012-11-09 05:16:08  mosteiro
 * Fill array variables
 *
 * Revision 1.123  2012-11-01 02:26:46  mosteiro
 * Fill i4_npmts_dt1 and i4_npmts_dt2
 *
 * Revision 1.122  2012-05-15 14:04:05  davini
 * bugfix in charge_thresh computation, thanks to Livia
 *
 * Revision 1.121  2012-02-15 16:27:14  davini
 * bugfix in charge_clean computation, thanks to Oleg
 *
 * Revision 1.120  2011-04-15 11:49:07  meindl
 * Modified first hit for neutrons (tt128 & tt1) : first hit set to zero only for low nhits clusters (< 200)
 *
 * Revision 1.119  2011-04-15 11:45:21  meindl
 * Modified first hit for neutrons (tt128 & tt1) : first hit set to zero only for low nhits clusters (< 200)
 *
 * Revision 1.118  2011-04-15 08:45:40  meindl
 * Correct normalization of hits for new variable "NClustersOld"
 *
 * Revision 1.117  2011-04-14 17:12:53  meindl
 * - Filling of new variable "n_clusters_old"
 * - start_hits for tt128 and tt1 neutrons is set to zero
 *
 * Revision 1.116  2011-04-01 12:23:30  meindl
 * New function "clusterize_neutrons_in_muongate" :
 * Uses cycle13 neutron clustering for neutron detection inside the muon gate
 *
 * Revision 1.115  2011-03-24 15:09:06  razeto
 * Added charge_clean
 *
 * Revision 1.114  2011-03-24 10:24:54  meindl
 * Modified criteria in data_start and noise determination for neutron clustering in muon events
 *
 * Revision 1.113  2011-03-23 15:40:50  meindl
 * - Muting debug messages
 * - Removing unused code
 *
 * Revision 1.112  2011-03-14 18:02:17  meindl
 * Parameter tuning.
 *
 * Revision 1.111  2011-03-07 13:33:17  meindl
 * Tailline modification
 *
 * Revision 1.110  2011-03-07 13:26:15  meindl
 * Increase of sensitivity by removing bin-probability requirements
 *
 * Revision 1.109  2011-03-07 08:30:09  meindl
 * - bug fix of incorrect vector access
 * - type declaration float->unsigned int
 *
 * Revision 1.104  2011-03-02 11:20:25  meindl
 * Modification of poisson scan interval.
 *
 * Revision 1.103  2011-03-02 10:45:05  meindl
 * Scan of poisson-statistic enlarged.
 *
 * Revision 1.102  2011-03-02 07:42:34  meindl
 * Fixed warning messages.
 *
 * Revision 1.101  2011-03-01 15:36:34  davini
 * npmts_tresh, nhits_thresh, charge_thresh evalued only in 400 ns since cluster start
 *
 * Revision 1.100  2011-02-28 16:36:23  meindl
 * Linear integration for background estimation.
 *
 * Revision 1.99  2011-02-28 14:38:38  meindl
 * New neutron algorithm implemented.
 *
 * Revision 1.98  2011-02-28 13:02:29  meindl
 * *** empty log message ***
 *
 * Revision 1.97  2011-02-18 16:05:01  davini
 * added npmts_400 nhits_400 charge_400 on Oleg's request; added charge_npmts based on Alessandro, Stefano and Livia idea;
 *
 * Revision 1.96  2010-07-03 09:47:20  meindl
 * Introduced protection against zero cluster events in the background determination.
 *
 * Revision 1.95  2010-07-02 13:15:11  meindl
 * Modification of background_time_window.
 *
 * Revision 1.94  2010-07-01 16:23:04  meindl
 * Alternative calculation of the nhits_bkg for tt1.
 *
 * Revision 1.93  2010-07-01 09:45:31  meindl
 * Errata:
 * The neutron clustering takes 160ns before the cluster start for background determination (80ns if there is a cluster pile-up).
 * CVS ----------------------------------------------------------------------
 *
 * Revision 1.92  2010-07-01 09:34:18  meindl
 * Added estimation of background hits in tt1 short clusters:
 * Variable "laben.clusters.nhits_bkg" is now filled for tt1 short clusters by scaling the number of decoded hits in the whole gate to the (short) cluster-duration.
 *
 * Note: for tt128 clusters, this variable is filled by the neutron clustering algorithm by scaling the background 30ns before cluster start to the whole cluster duration.
 *
 * Revision 1.91  2009-11-18 12:18:50  razeto
 * Added some _short variables
 *
 * Revision 1.90  2009-10-26 11:19:37  ddangelo
 * trgtype 'muon_tct' (or 'muon_totc') renamed as 'neutron' after real use
 * trgtype 'muon_mtb' renamed as 'muon' as duplicity with the above one disappeared and 'mtb' is associated to btb_inputs flag
 *
 * Revision 1.89  2009-10-26 09:52:35  meindl
 * Changed trgtype128 condition to is_muon_tct condition.
 *
 * Revision 1.88  2009-10-26 09:08:45  meindl
 * Fixed bug which set the background variable for triggertype128 events to zero.
 *
 * Revision 1.87  2009-10-24 09:14:49  meindl
 * Change in background determination (one bin more).
 *
 * Revision 1.86  2009-10-23 16:49:15  meindl
 * Now muons are also clustered with the neutron-finding-algorithm.
 * Results are saved in the new cluster-vector "laben.clusters_muons".
 * Also some bug-fixes concerning the background-hits-variable.
 *
 * Revision 1.85  2009-10-23 11:28:40  meindl
 * Muted debug-messages
 *
 * Revision 1.84  2009-10-22 19:04:59  meindl
 * Removed exception of neutron clustering for triggertype128 events with n_decoded_hits<2000.
 * Added background determination for neutron clustering (new variable is zero for main clustering)
 *
 * Revision 1.83  2009-09-17 15:55:41  razeto
 * Compute and assign short values
 *
 * Revision 1.82  2009-09-16 13:20:06  razeto
 * Added short/long clustering marking, waiting for event variable commit
 *
 * Revision 1.81  2009-09-16 11:16:10  razeto
 * Long cluster (aka m4) is now a parameter (defaulting at 1.5us)
 *
 * Revision 1.80  2009-09-10 11:59:51  razeto
 * Added tresh variable calculation (required by oleg)
 *
 * Revision 1.79  2009-08-04 16:01:19  alvaro
 * Updated dark_rate value and fixed unsigned-signed comparison warning
 *
 * Revision 1.77  2009-07-23 10:56:28  razeto
 * Better 500ns definition
 *
 * Revision 1.76  2009-07-17 15:40:20  ddangelo
 * writing rough time commented out
 *
 * Revision 1.75  2009-07-07 08:27:29  razeto
 * Added minimum cluster length at 500ns (M4 merging)
 *
 * Revision 1.74  2008-12-15 12:12:44  razeto
 * Added rms_time to cluster (to improve PID)
 *
 * Revision 1.73  2008-12-10 11:40:24  razeto
 * Added mean charge and peak from tt1
 *
 * Revision 1.72  2008-10-20 14:02:26  meindl
 * Restored version 1.66; added threshold of 5hits for neutron clustering
 *
 * Revision 1.62.2.1  2008-09-30 16:43:06  meindl
 * commented modfied version of the tail contribution ; added alternative double peak identification
 *
 * Revision 1.65  2008-08-27 08:55:02  meindl
 * Modification of falling-edge condition.
 *
 * Revision 1.64  2008-08-26 13:26:10  meindl
 * New neutron clustering algorithm implemented.
 *
 * Revision 1.63  2008-08-26 12:57:24  meindl
 * *** empty log message ***
 *
 * Revision 1.62  2008-02-27 14:49:24  razeto
 * Lower threshold for rising edge (+ some small upgrades)
 *
 * Revision 1.61  2008-02-26 17:53:12  razeto
 * Hmmm positive means > -2?
 *
 * Revision 1.60  2008-02-26 17:26:27  razeto
 * Require only positive defivative for rising edge
 *
 * Revision 1.59  2008-02-26 16:57:06  razeto
 * Stop falling search after 300ns
 *
 * Revision 1.58  2008-02-26 16:42:05  razeto
 * Lowered thresholds for cluster integral
 *
 * Revision 1.57  2008-02-26 16:14:29  razeto
 * Rising edge with only 3 points. (better printout)
 *
 * Revision 1.56  2008-02-26 12:13:25  razeto
 * Minor fixes, baseline rejection still disabled
 *
 * Revision 1.55  2008-02-26 10:21:49  razeto
 * Upgraded to better working
 *
 * Revision 1.54  2008-02-21 18:22:47  razeto
 * Added a preliminary differential clustering for neutrons
 *
 * Revision 1.53  2008-02-15 10:45:58  razeto
 * clusterize internal method created, for more modularity. Added clusterize_neutrons (empty)
 *
 * Revision 1.52  2007-12-21 19:14:59  razeto
 * Fixed long gate with neutrons
 *
 * Revision 1.51  2007-11-15 21:09:45  razeto
 * Fixed bug in flag using real signs for gate parameters
 *
 * Revision 1.50  2007-11-12 15:20:29  razeto
 * Filling new flag variable
 *
 * Revision 1.49  2007-11-12 12:33:41  razeto
 * Ripple is not part of the cluster, but part of the tail (for better cluster lenght)
 *
 * Revision 1.48  2007-11-12 12:08:18  razeto
 * moved dark rate calculation in evaluate_ripple_count
 *
 * Revision 1.47  2007-11-12 09:54:51  razeto
 * Bug fix: n_hits not initalized caused unstable results with optimization.
 * Removed a warning.
 * Added an overflow check.
 *
 * Revision 1.46  2007-11-09 22:56:17  razeto
 * Allow a delayed end of gate (500ns)
 *
 * Revision 1.45  2007-11-09 18:58:40  razeto
 * Fixed a bug for strict_mode (not default) and now use the gate size from db
 *
 * Revision 1.44  2007-11-05 23:46:32  razeto
 * Removed a wrong comment
 *
 * Revision 1.43  2007-11-05 23:42:28  razeto
 * Code more stable
 *
 * Revision 1.42  2007-10-31 17:12:39  razeto
 * Code debugged for normal events; still in beta phase for neutrons
 *
 * Revision 1.41  2007-10-30 18:06:47  razeto
 * Experimental code: work in gateless mode for events close to the previous.
 *
 * Revision 1.40  2007-10-30 15:55:15  ddangelo
 * just variable name updated
 *
 * Revision 1.39  2007-07-05 11:01:34  razeto
 * Doing clustering on type 2/128 events too
 *
 * Revision 1.38  2007-05-05 12:42:04  razeto
 * Better hangling of high and low energy events in clustering
 *
 * Revision 1.37  2007-05-03 12:05:37  razeto
 * Added broad flag and rough time
 *
 * Revision 1.36  2007-04-27 13:53:17  razeto
 * Fill cluster rough time
 *
 * Revision 1.35  2007-04-15 16:53:31  razeto
 * Added high energy mode to handle 140ns dead time
 *
 * Revision 1.34  2007-04-15 14:19:10  razeto
 * New end of cluster
 *
 * Revision 1.33  2007-03-27 15:20:51  ddangelo
 * variables *_conc filled
 *
 * Revision 1.32  2007-03-24 15:57:11  ddangelo
 * debugging
 * varaibles npmt and npmt_conc are now filled correctly
 *
 * Revision 1.31  2006-11-03 13:31:55  razeto
 * Fixed 2 bugs in clustering
 *
 * Revision 1.30  2006/11/03 13:00:47  razeto
 * Added a quality check
 *
 * Revision 1.29  2006/10/23 13:06:27  razeto
 * Reduced the clustering verbosity, and added clustering to random events
 *
 * Revision 1.28  2006/08/21 11:19:02  razeto
 * Updated to new barn_interface
 *
 * Revision 1.27  2006/06/29 15:03:25  razeto
 * Added indexes to lower level hits
 *
 * Revision 1.26  2006/05/08 17:31:31  razeto
 * Added db_channel patch for fadc (sent to the mailing list)
 *
 * Revision 1.25  2006/01/02 21:23:46  razeto
 * Changed test target to file target for db_barn
 *
 * Revision 1.24  2005/11/18 16:46:06  razeto
 * Fill pe and npe variables
 *
 * Revision 1.23  2005/09/28 13:30:36  razeto
 * Introduced a new cut (more relaxed) to select the first hit. To be tested deeply
 *
 * Revision 1.22  2005/09/19 16:33:39  razeto
 * Removed accumulate since is buggy, fixed a bug
 *
 * Revision 1.21  2005/09/19 11:41:29  razeto
 * Added an option to teste the splitting
 *
 * Revision 1.20  2005/09/15 14:52:48  razeto
 * Fixed a bug in second cluster finding
 *
 * Revision 1.19  2005/07/26 10:39:13  razeto
 * debugging the algorithm
 *
 * Revision 1.18  2005/07/15 13:51:10  razeto
 * Updated to use a new algorithm
 *
 * Revision 1.17  2005/07/13 13:02:25  razeto
 * Fixed a typo
 *
 * Revision 1.16  2005/07/13 12:32:08  razeto
 * Merged the the 2 laben clustering modules: now find_cluster_time is
 * a subroutine of the bx_laben_findcluster module
 *
 * Revision 1.15  2005/05/29 17:24:59  razeto
 * Commented a log, fixed a check constant
 *
 * Revision 1.14  2005/05/29 17:08:01  razeto
 * Added a print and fix a small log bug
 *
 * Revision 1.13  2005/05/29 16:48:59  razeto
 * Changed algorithm
 *
 * Revision 1.12  2005/04/29 14:53:43  razeto
 * Moved order_in_channel to bx_laben_findcluster_time
 *
 * Revision 1.11  2005/03/18 11:09:53  razeto
 * Added filling of order_in_channel, fixed indentation and some small upgrades
 *
 * Revision 1.10  2005/03/01 15:18:18  razeto
 * Merged with cycle_2
 *
 * Revision 1.9.2.1  2004/12/13 12:44:17  razeto
 * Commented an unused variable
 *
 * Revision 1.9  2004/11/29 13:21:23  razeto
 * Added Mantainer field
 *
 * Revision 1.8  2004/09/28 12:54:49  dmanuzio
 * Removed ifdef for root barn; added a new parameter to enable the histos in
 * the bx_laben_findcluster_time module and added a check on empty channels
 * in the bx_laben_findcluster module.
 *
 * Revision 1.7  2004/09/22 14:20:47  dmanuzio
 * Changed bx_reco_event to bx_echidna_event
 *
 * Revision 1.6  2004/09/22 12:45:40  dmanuzio
 * Update to follow Alessandro's rfc with the new bx_detector class
 *
 * Revision 1.5  2004/07/15 13:13:47  dmanuzio
 * Removed the bug. Test done!
 *
 * Revision 1.4  2004/07/15 13:12:51  dmanuzio
 * Linking error put on purpose to test EchidnaTest script
 *
 * Revision 1.3  2004/07/15 13:04:13  dmanuzio
 * Test done
 *
 * Revision 1.2  2004/07/15 13:03:29  dmanuzio
 * Test of a new utility
 *
 * Revision 1.1  2004/07/09 13:59:08  dmanuzio
 * Added 2 modules for Laben clustering
 *
 *
 */
