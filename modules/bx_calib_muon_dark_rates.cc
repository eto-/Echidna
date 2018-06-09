/* BOREINO Reconstruction program
 *
 * Author: Davide D'Angelo <davide.dangelo@mi.infn.it>
 * Maintainer: Davide D'Angelo <davide.dangelo@mi.infn.it>
 *             inspired by bx_calib_dark_rates

 * $Id: bx_calib_muon_dark_rates.cc,v 1.4 2007/05/31 14:35:06 ddangelo Exp $
 * 
 * Implementation of bx_calib_muon_dark_rates
 *
 */

#include <cmath>
//#include <iostream>
#include <vector>
#include <algorithm>

#include "bx_calib_muon_dark_rates.hh"
#include "bx_echidna_event.hh"
#include "bx_dbi.hh"
#include "db_profile.hh"
#include "db_channel.hh"
#include "barn_interface.hh"
#include "TF1.h"
#include "TH1F.h"
#include "db_run.hh"

// ctor
bx_calib_muon_dark_rates::bx_calib_muon_dark_rates () : bx_base_module("bx_calib_muon_dark_rates", bx_base_module::main_loop) {
  require_event_stage (bx_detector::muon, bx_base_event::decoded);
  require_trigger_type (bx_trigger_event::random);

  translation_map ["no_pmt"] = no_pmt;
  translation_map ["good"] = good; 
  translation_map ["dead"] = dead; 
  translation_map ["hot"] = hot;
  translation_map ["new_no_pmt"] = new_no_pmt;
  translation_map ["new_good"] = new_good; 
  translation_map ["new_dead"] = new_dead;
  translation_map ["new_hot"] = new_hot;
}
      
// begin
void bx_calib_muon_dark_rates::begin () {
  int nch = constants::muon::channels;

  // reset ctrs
  i4_trg_ctr = 0; 
  v_nhits.resize (nch,0);

  // get user parameters
  f4_dark_rate_thresh_high = get_parameter ("dark_rate_thresh_high").get_float ();
  f4_dark_rate_thresh_low = get_parameter ("dark_rate_thresh_low").get_float ();
  f4_gate_start = get_parameter ("gate_start").get_int();
  f4_gate_end   = get_parameter ("gate_end").get_int();
  f4_gate_width = f4_gate_end-f4_gate_start;

  // Histograms for barn 
  ok_dark_rates = new TH1F ("OD_ok_dark_rates","OD Dark Rates of good PMTs", 10 * (int) (1.1 * f4_dark_rate_thresh_high)  , 0,  1.1 * f4_dark_rate_thresh_high);
  ok_dark_rates->SetXTitle("kHz");
  ok_dark_rates->SetYTitle("channels/100Hz");
  barn_interface::get ()->store (barn_interface::file, ok_dark_rates, this);
  
  dark_rate_vs_channel = new TH1F ("OD_dark_rate_vs_channel","OD Dark Rates vs channel", nch, 1, nch + 1);
  dark_rate_vs_channel->SetXTitle("muon channel number");
  dark_rate_vs_channel->SetYTitle("dark rate (kHz)"); 
  barn_interface::get ()->store (barn_interface::file, dark_rate_vs_channel, this);
   
  // histograms for detector_monitor
  //    random_hits_map_cumulative = new TH2F ("random_hits_map_cumulative","z-axis: number of random hits, uploaded every 200 random triggers", 
  // 					  72, 1, 73, 43, -21.5, 21.5);
  //    random_hits_map_cumulative->SetYTitle("ring number, + 21 north pole, -21 south pole");
  //    random_hits_map_cumulative->SetXTitle("position in ring");
     
  //    random_hits_map = new TH2F ("random_hits_map","z-axis: number of random hits, last 1000 random triggers", 72, 1, 73, 43, -21.5, 21.5);
  //    random_hits_map->SetYTitle("ring number, + 21 north pole, -21 south pole");
  //    random_hits_map->SetXTitle("position in ring");
   
}


// doit
bx_echidna_event* bx_calib_muon_dark_rates::doit (bx_echidna_event *ev) {
  
  i4_trg_ctr++;

  const bx_muon_event& er = ev->get_muon ();
  
  //loop on decoded hits
  for(int i = 0; i < er.get_decoded_nhits (); i++) {

    // skip non-ordinary channels
    if(!er.get_decoded_hit (i).get_db_channel ()->is_ordinary () ) continue;

    int mch = er.get_decoded_hit (i).get_raw_hit().get_muon_channel ();
    float time = er.get_decoded_hit (i).get_time ();
                  
    //fill maps
    //    int hole_id = er.get_decoded_hit (i).get_db_channel ()->pmt_hole_id ();
    //    int ring = hole_id / 100;
    //    int column = abs(hole_id % 100);
        
    //    random_hits_map_cumulative->Fill(column, ring); 	
    //    random_hits_map-> Fill(column, ring); 

    //fill vector with hits
    if ( time > f4_gate_start && time < f4_gate_end ) 
      v_nhits[mch]++; 

    //    if (!(count_random_triggers % 200))   barn_interface::get ()->network_send (random_hits_map_cumulative, this);
    //    if (!(count_random_triggers % 1000)) {
    //      barn_interface::get ()->network_send (random_hits_map, this);
    //      random_hits_map->Reset ();

  }
  return ev;     
}


// end
void bx_calib_muon_dark_rates::end () {

  // effective time for dark noise measurements in ns
  float exposure_time = i4_trg_ctr * f4_gate_width;
  get_message(bx_message::log) << "Number of random triggers: " << i4_trg_ctr << "; exposure time: " << exposure_time / 1e9 << " s." << dispatch;  
  if (i4_trg_ctr < get_parameter ("min_statistics").get_int ()) {
    get_message(bx_message::warn) << "Not enough statistics to calculate dark rates." << dispatch ; 
    return;
  }

  // channel counters
  int n_ordinary_lg = constants::muon::channels; // to be decremented later
  int dead_ctr = 0; 
  int hot_ctr = 0; 

  // variables for calculating mean dark rate of the detector
  // _all: sum/mean over all pmts in ordinary lg; _ok sum/mean excluding hot/dead 
  double sum_rate_all = 0, sum_rate_ok = 0;      
  float mean_rate_all = 0, mean_rate_ok = 0;
  float mean_sigma_ok = 0; 

  std::vector<double> dark_rates (constants::muon::channels);
  std::vector<double> dark_sigmas (constants::muon::channels);  
  std::vector<pmt_status_t> pmt_status(constants::muon::channels);
  
  int nhits_all = 0;

  // loop on mch
  for (int mch = 0; mch < constants::muon::channels; mch++) {
    int lg = constants::muon::channel_offset + mch + 1;
    if ( (bx_dbi::get ()->get_channel (lg).is_ordinary ()) ) {

      // we are optimistic, pmt is good until someone says the opposite
      pmt_status[mch] = good;
	  
      // rate in kHz
      dark_rates[mch] = v_nhits[mch] / (exposure_time * 1e-9) / 1e3; // exposure time is in ns   
      dark_sigmas[mch] = v_nhits[mch] ? dark_rates[mch] / sqrt(v_nhits[mch]) : 0. ;
      dark_rate_vs_channel->SetBinContent(mch, dark_rates[mch]);
      nhits_all += v_nhits[mch];
      
      // high rate, mark as hot 
      if(dark_rates[mch] > f4_dark_rate_thresh_high){
	hot_ctr ++;           
	get_message(bx_message::log) << "PMT in mch " << mch 
				     << " HOT, dark rate " << dark_rates[mch]  
				     << " +- " << dark_sigmas[mch] << " kHz, " 
				     << v_nhits[mch] << " hits" << dispatch;
	//set pmt_status
	pmt_status[mch] = hot;
      }

      // low rate, mark as dead 
      if (dark_rates[mch] < f4_dark_rate_thresh_low){
	//is it 2sigma below the threshold?
	double  n_hits_expected =  exposure_time / 1e9 * (f4_dark_rate_thresh_low * 1000);
	if(v_nhits[mch] <  (n_hits_expected - 6 * sqrt(n_hits_expected))){
	  dead_ctr ++;           
	  get_message(bx_message::log) << "PMT in mch " << mch 
				       << " DEAD, dark rate " << dark_rates[mch] 
				       << " +- " << dark_sigmas[mch] << " kHz, " 
				       << v_nhits[mch] << " hits" << dispatch;
	  //set pmt_status
	  pmt_status[mch] = dead;
	}
      }
	
      //to sum rates of PMTs 
      //all ordinary channels
      sum_rate_all += dark_rates[mch];
      //if PMT is ok
      if(pmt_status[mch] == good){
	ok_dark_rates->Fill(dark_rates[mch]); 
	sum_rate_ok += dark_rates[mch];
      }
    }
    else { // non-ordinary lg
      n_ordinary_lg --;
      dark_rates[mch] = -1.;      
      dark_sigmas[mch] = -1.;
      pmt_status[mch] = no_pmt;
    } 

  } // end of loop on channels

  int n_ok_pmts = n_ordinary_lg - hot_ctr - dead_ctr;

  //barn_interface::get ()->network_send (dark_rate_vs_channel, this);
                                   
  //mean rate in kHz
  mean_rate_all =  sum_rate_all / n_ordinary_lg;
   
  //mean total rate for ok PMTs
  mean_rate_ok  =  ok_dark_rates->GetMean ();
  mean_sigma_ok =  ok_dark_rates->GetRMS ();
 
  get_message(bx_message::log) << "Detector sum dark rate: " << sum_rate_all << " kHz." << dispatch;
  get_message(bx_message::log) << "Detector mean dark rate: " <<  mean_rate_all << " kHz (all ordinary). " 
			       << mean_rate_ok << " +- " << mean_sigma_ok << " kHz (good only). " << dispatch;
   
  get_message(bx_message::log) << "Number of PMTs in ordinary lg: " << n_ordinary_lg << dispatch;
  get_message(bx_message::log) << "Number of hot PMTs: " << hot_ctr << dispatch;    
  get_message(bx_message::log) << "Number of dead PMTs: " << dead_ctr << dispatch;    
  get_message(bx_message::log) << "Number of ok PMTs: " << n_ok_pmts << dispatch;

  //writing to visitors and checking changes in dark rates with respect to the latest calibrated run
  db_run& run_info = bx_dbi::get ()->get_run ();
   
  //warning if mean dark rate of the whole detector has changed
  float previous_mean_dark_noise = run_info.get_muon_mean_dark_noise ();
  float previous_mean_dark_sigma = run_info.get_muon_mean_dark_sigma ();
  if ( ::fabs(previous_mean_dark_noise - mean_rate_ok) > 3 * ((double) mean_sigma_ok + (double) previous_mean_dark_sigma) ) 
    get_message(bx_message::warn) << "Mean dark rate has changed from " <<  previous_mean_dark_noise 
				  << " +- " <<  previous_mean_dark_sigma << " kHz to " 
				  << mean_rate_ok << " +- " <<  mean_sigma_ok << " kHz " << dispatch;

  //set visitors
  run_info.set_muon_mean_dark_noise (mean_rate_ok, this);
  run_info.set_muon_mean_dark_sigma (mean_sigma_ok, this);
  run_info.set_muon_dead (dead_ctr, this);
  run_info.set_muon_hot (hot_ctr, this);           

  // checking individual channels
  for(int mch = 0; mch < constants::muon::channels; mch++) {
    int lg = mch + constants::muon::channel_offset + 1;

    // get values from the previous calibration run     
    float previous_muon_dark_noise = run_info.get_muon_dark_noise (lg) ;
    float previous_muon_dark_sigma = run_info.get_muon_dark_sigma (lg) ;
	
    //compare to present values
    if( ::fabs(previous_muon_dark_noise - dark_rates[mch]) > ::fabs( 3 * ((double) previous_muon_dark_sigma + (double) dark_sigmas [mch])) ) 
      get_message(bx_message::warn) << "Dark rate for channel " << mch 
				    << " has changed from " <<  previous_muon_dark_noise << " +- " <<  previous_muon_dark_sigma 
				    << " kHz to " << dark_rates[mch] << " +- " <<  dark_sigmas[mch] << " kHz " << dispatch;
      
    //set visitors for individual rates
    run_info.set_muon_dark_noise  (lg, dark_rates[mch], this);
    run_info.set_muon_dark_sigma  (lg, dark_sigmas[mch], this);
     
    //compare current and previous status  and set visitors
    const std::string& previous_muon_pmt_status_descr = run_info.get_muon_pmt_status (lg);
    pmt_status_t previous_muon_pmt_status;
    if (!translation_map.check (previous_muon_pmt_status_descr)) {
      get_message(bx_message::error) << "unknown pmt status" << previous_muon_pmt_status_descr << " for ch " << mch << dispatch;
      previous_muon_pmt_status = good;
    } 
    else 
      previous_muon_pmt_status = pmt_status_t(int(translation_map[previous_muon_pmt_status_descr]) % 10);
    
    if(previous_muon_pmt_status != pmt_status[mch]){
      get_message(bx_message::warn) << "Dark rate status changed for mch " << mch 
				    << " from " <<  previous_muon_pmt_status 
				    << " to " << pmt_status[mch] << dispatch;
      pmt_status [mch] = pmt_status_t(int(pmt_status [mch]) + 10);
    }
     
    run_info.set_muon_pmt_status  (lg, translation_map.rfind (pmt_status[mch])->first, this);

  } // end of loop on channels

    // write to the DB
  if(get_parameter ("db_write").get_bool ()) {
    run_info.write_muon_dark_rates (true, this);
    get_message(bx_message::log) << "Writing to DB" << dispatch; 
  }
}
/*
 * $ Log: $ 
 *
 */
 
 



