/* BOREINO Reconstruction program
 *
 * Author: Daniela Manuzio <dmanuzio@ge.infn.it> and Livia Ludhova<livia.ludhova@mi.infn.it> 
 * Maintainer: Livia Ludhova<livia.ludhova@mi.infn.it>
 *
 * 
 * Implementation of bx_calib_laben_dark_rates
 */

#include <cmath>
#include <iostream>
#include <vector>
#include <algorithm>
#include "bx_calib_laben_dark_rates.hh"
#include "bx_echidna_event.hh"
#include "bx_dbi.hh"
#include "db_profile.hh"
#include "db_channel.hh"
#include "barn_interface.hh"
#include "TF1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "db_run.hh"

bx_calib_laben_dark_rates::bx_calib_laben_dark_rates () : bx_base_module("bx_calib_laben_dark_rates", bx_base_module::main_loop) {
  require_event_stage (bx_detector::laben, bx_base_event::decoded);

  translation_map ["no_pmt"] = no_pmt;
  translation_map ["good"] = good; 
  translation_map ["dead"] = dead; 
  translation_map ["hot"] = hot;
  translation_map ["new_no_pmt"] = new_no_pmt;
  translation_map ["new_good"] = new_good; 
  translation_map ["new_dead"] = new_dead;
  translation_map ["new_hot"] = new_hot;
}
      
  //BEGIN
void bx_calib_laben_dark_rates::begin () {


  index_trg_type =  get_parameter ("index_trg_type").get_int ();
  min_dark_time = get_parameter ("min_dark_time").get_float ();
  dark_rate_thresh_high = get_parameter ("dark_rate_thresh_high").get_float ();
  dark_rate_thresh_low = get_parameter ("dark_rate_thresh_low").get_float ();

  count_random_triggers = 0; 
  count_pulser_triggers = 0; 
  count_laser_triggers = 0; 
  count_neutrino_triggers = 0; 

  width_neutrino = 1300;  
  db_run& run_info = bx_dbi::get ()->get_run ();
  gate_width   = run_info.get_laben_gate_width ();	     
  trigger_offset  = -1 * run_info.get_laben_gate_start ();
  time_laser   = trigger_offset - ( -1 * run_info.get_laben_laser_offset ());
  time_pulser  = trigger_offset - ( -1 * run_info.get_laben_pulser_offset ());
  time_neutrino = trigger_offset - (-1 * run_info.get_laben_cluster_offset ());
 
  //time_pulser = 4900.;
  //time_laser =  2800.; 
  //time_neutrino = 980.;
  //trigger_offset = 6660.; 
  //gate_width = 7050.;



  get_message(bx_message::log) << "Trigger parameters (ns): trigger offset " << trigger_offset << " gate width " << gate_width << " time pulser " << time_pulser <<  " time laser " << time_laser << " time_neutrino " << time_neutrino << dispatch; 
  
  //number of time bins for histos, 50 ns binning
  time_after_gate = 2000.; 
  time_before_gate = 2000.; 
  time_bins = (int32_t) (gate_width + time_before_gate + time_after_gate) / 50; 
  

      // Histograms
  int32_t nch = constants::laben::channels;

  ID_ok_dark_rates = new TH1F ("ID_ok_dark_rates","of good inner det. PMTs", 20 * (int32_t) (1.1 * dark_rate_thresh_high)  , 0,  1.1 * dark_rate_thresh_high);
  ID_ok_dark_rates->SetXTitle("dark rate of good ID PMTs [kHZ], 50 Hz bin");

  ID_all_dark_rates = new TH1F ("ID_all_dark_rates","of all alive PMTs", 20 * 100 , 0, 100);
  ID_all_dark_rates->SetXTitle("dark rate of all alive ID PMTs [kHZ], 50 Hz bin");
  
  ID_dark_rate_vs_channel = new TH2F ("ID_dark_rate_vs_channel","z-axes: dark rate in kHz (non ordinary lg z = -1)", nch, 1, nch + 1, 5, 1, 6);
  ID_dark_rate_vs_channel->SetXTitle("logical channel number");
  ID_dark_rate_vs_channel->SetYTitle("1-random, 2-pulser, 3-laser, 4-neutrino, 5-sum");
 
  barn_interface::get ()->store (barn_interface::file, ID_dark_rate_vs_channel, this);
  barn_interface::get ()->store (barn_interface::file, ID_ok_dark_rates, this);
  barn_interface::get ()->store (barn_interface::file, ID_all_dark_rates, this);
   
  //histograms for monitoring
   random_hits_map_cumulative = new TH2F ("random_hits_map_cumulative","z-axis: number of random hits, uploaded every 200 random triggers", 72, 1, 73, 43, -21.5, 21.5);
   random_hits_map_cumulative->SetYTitle("ring number, + 21 north pole, -21 south pole");
   random_hits_map_cumulative->SetXTitle("position in ring");
     

   random_hits_map = new TH2F ("random_hits_map","z-axis: number of random hits, last 1000 random triggers", 72, 1, 73, 43, -21.5, 21.5);
   random_hits_map->SetYTitle("ring number, + 21 north pole, -21 south pole");
   random_hits_map->SetXTitle("position in ring");
   

   nhits_random.resize (nch,0);
   nhits_pulser.resize (nch,0);
   nhits_laser.resize (nch,0);
   nhits_neutrino.resize (nch,0);
}


  //DOIT
bx_echidna_event* bx_calib_laben_dark_rates::doit (bx_echidna_event *ev) {
  
  //  //from random triggers
  if(ev->get_trigger().is_random()) {    
    count_random_triggers++;

    const bx_laben_event& er = ev->get_laben ();
  
      //loop in decoded hits
    for(int32_t i = 0; i < er.get_decoded_nhits (); i++) {
      int32_t lg = er.get_decoded_hit (i).get_raw_hit ().get_logical_channel ();
      float time = er.get_decoded_hit (i).get_raw_time () - er.get_trigger_rawt () + trigger_offset;
                  
      //fill maps
      if(er.get_decoded_hit (i).get_db_channel ()->is_ordinary () ) {
	int32_t hole_id = er.get_decoded_hit (i).get_db_channel ()->pmt_hole_id ();
	int32_t ring = hole_id / 100;
	int32_t azimut = abs(hole_id % 100);
        
        random_hits_map_cumulative->Fill(azimut, ring); 	
	random_hits_map-> Fill(azimut, ring); 
      }
      

      //fill vector with hits
      if ((time > 200.) && (time < (gate_width - 200.)) && !(er.get_decoded_hit (i).is_out_of_gate ()) ) nhits_random[lg - 1] += 1; 
      if ((time > 200.) && (time < (gate_width - 200.)) && (er.get_decoded_hit (i).is_out_of_gate ()) ) {
	get_message(bx_message::warn) << "WEIRD: hit which has time as in gate is marked as out_of_gate " << dispatch ;
      }	
    }
    if (!(count_random_triggers % 200))   barn_interface::get ()->network_send (random_hits_map_cumulative, this);
    if (!(count_random_triggers % 1000)) {
      barn_interface::get ()->network_send (random_hits_map, this);
      random_hits_map->Reset ();
    }
  }
  
  
  //from pulser
  if(ev->get_trigger().is_pulser()) {    
    count_pulser_triggers++;
    const bx_laben_event& er = ev->get_laben ();
    //loop in decoded hits
    for(int32_t i = 0; i < er.get_decoded_nhits (); i++) {
      int32_t lg = er.get_decoded_hit (i).get_raw_hit ().get_logical_channel ();
      float time = er.get_decoded_hit (i).get_raw_time () - er.get_trigger_rawt () + trigger_offset;
      if((time > 200.) &&  (time < (time_pulser - 200.)) && !(er.get_decoded_hit (i).is_out_of_gate ())) nhits_pulser[lg - 1] += 1;
      if((time > 200.) &&  (time < (time_pulser - 200.)) && (er.get_decoded_hit (i).is_out_of_gate ())){
	get_message(bx_message::warn) << "WIERD: hit which has time as in gate is marked as out_of_gate " << dispatch ;
	}
      if( (time > (time_pulser - 200.)) ) break;  
    }
  }
  
  //from timing laser
  if(ev->get_trigger().is_laser394()) {    
    count_laser_triggers++;
    const bx_laben_event& er = ev->get_laben ();
    //loop in decoded hits
    for(int32_t i = 0; i < er.get_decoded_nhits (); i++) {
      int32_t lg = er.get_decoded_hit (i).get_raw_hit ().get_logical_channel ();
      float time = er.get_decoded_hit (i).get_raw_time () - er.get_trigger_rawt () + trigger_offset;
      if((time > 200.) &&  (time < (time_laser - 200.)) && !(er.get_decoded_hit (i).is_out_of_gate ())) nhits_laser[lg - 1] += 1;
      if((time > 200.) &&  (time < (time_laser - 200.)) && (er.get_decoded_hit (i).is_out_of_gate ())){
	get_message(bx_message::warn) << "WIERD: hit which has time as in gate is marked as out_of_gate " << dispatch ;
	}
      if( (time > (time_laser - 200.)) ) break;
    }
  }

      
  // from 1 cluster events, hits after the cluster
  
  if(ev->get_trigger ().is_neutrino ()) { 
    const bx_laben_event& er = ev->get_laben ();
    //only 1-cluster events
    if( (er.get_nclusters ()) == 1.){ 
      int32_t N_clustered_hits =  er.get_cluster (0).get_clustered_nhits ();
        //time of the the last clustered hit
      double time_last = er.get_cluster (0).get_clustered_hit (N_clustered_hits - 1).get_decoded_hit ().get_raw_time () - er.get_trigger_rawt () + trigger_offset;

        //if cluster finishes inside time_neutrino + width_neutrino and there is space after the ned of cluster and end of the gate for some hits
      if(((time_neutrino + width_neutrino) < (gate_width - 200)) && (time_last < (time_neutrino + width_neutrino)) ) {
	count_neutrino_triggers++;
  	  //loop in decoded hits before cluster
	for(int32_t i = N_clustered_hits; i < er.get_decoded_nhits (); i++) {
	  float time = er.get_decoded_hit (i).get_raw_time () - er.get_trigger_rawt () + trigger_offset;
	  int32_t lg = er.get_decoded_hit (i).get_raw_hit ().get_logical_channel ();
	  //hit after the cluster 
          if((time > (time_neutrino + width_neutrino)) && (time < (gate_width - 200))  &&  !(er.get_decoded_hit (i).is_out_of_gate ())) nhits_neutrino[lg - 1] += 1;
	  if((time > 200.) &&  (time < (time_neutrino - 100.)) && (er.get_decoded_hit (i).is_out_of_gate ())){
	    get_message(bx_message::warn) << "WIERD: hit which has time as in gate is marked as out_of_gate " << dispatch ;
	    }
	}
      }
    }
  }

   return ev;     
}


  //END
void bx_calib_laben_dark_rates::end () {

  int32_t index_random = (int32_t) index_trg_type / 1000;
  int32_t index_pulser = (int32_t) (index_trg_type % 1000)/100;
  int32_t index_laser = (int32_t) (index_trg_type % 100)/ 10;
  int32_t index_neutrino = (int32_t) (index_trg_type % 10);


  int32_t n_ordinary_lg = 0; //number of ordinary PMTs     
  int32_t n_ok_pmt = 0; //number of ok PMTs     
  int32_t count_not_ordinary_pmts = 0;
  int32_t count_cone_pmts_dead = 0; 
  int32_t count_cone_pmts_hot = 0; 
  int32_t count_not_cone_pmts_dead = 0; 
  int32_t count_not_cone_pmts_hot = 0; 

    //variables for calculating mean dark rate of the detector
  double sum_rates_all_pmts = 0;
  double sum_rates_ok_pmts = 0;
      
  float mean_rate, mean_ok_rate; // in kHz for the whole detector/per PMT, with/without bad pmts
  float error_mean_rate = 0;
  float error_mean_ok_rate = 0; 

    //effective times for dark hits measurements in ns
  float tot_random_time = count_random_triggers * (gate_width - 400.);
  float tot_pulser_time = count_pulser_triggers * (time_pulser - 400.);
  float tot_laser_time = count_laser_triggers * (time_laser - 400.);
  float tot_neutrino_time = count_neutrino_triggers * (gate_width - time_neutrino - width_neutrino - 200. );
  float tot_dark_time =  index_random * tot_random_time + index_pulser * tot_pulser_time + index_laser * tot_laser_time +  index_neutrino * tot_neutrino_time;

  get_message(bx_message::log) << "Number of random, pulser, laser, 1-cluster neutrino triggers " << count_random_triggers << " , " << count_pulser_triggers << " , " << count_laser_triggers << " , " <<  count_neutrino_triggers << dispatch; 
  get_message(bx_message::log) << "Random, pulser, laser, neutrino - dark times " << tot_random_time / 1e9 << " , " << tot_pulser_time  / 1e9 << " , " << tot_laser_time / 1e9 << " , " <<  tot_neutrino_time / 1e9 << " s " << dispatch; 
  get_message(bx_message::log) << "Total dark time " << tot_dark_time / 1e9 << " s " << dispatch;
 
  std::vector<double> dark_rates_total(constants::laben::channels);
  std::vector<double> error_dark_rates_total(constants::laben::channels);  
  std::vector<pmt_status>  alive(constants::laben::channels);
  
    //for non ordinary channels set dark_rate to -1, alive-pmt_status set to 0
  for(int32_t ch = 1; ch < (constants::laben::channels + 1); ch++) {  
    if(!(bx_dbi::get ()->get_channel (ch).is_ordinary ()) ) {
      count_not_ordinary_pmts ++;
      for(int32_t i = 1; i < 6; i++) ID_dark_rate_vs_channel->SetBinContent (ch, i, -1);
      dark_rates_total[ch - 1] = -1.;      
      error_dark_rates_total[ch - 1] = -1.;
      alive[ch - 1] = no_pmt;
    } //for all ordinary channels to set alive to 1 (ok); later for bad ones will be set to 2 or 3
    else  alive[ch - 1] = good;
  }
  
    //number of ordinary PMTs 
  n_ordinary_lg = constants::laben::channels - count_not_ordinary_pmts;
  n_ok_pmt = n_ordinary_lg; //hot and dead ones will be subtracted later
	
  if (tot_dark_time < min_dark_time * 1e9) get_message(bx_message::warn) << "Not enough statistics to calculate dark rates" << dispatch ; 
  else{ // enough stat
    int32_t nhits_random_all = 0;
    int32_t nhits_pulser_all = 0;
    int32_t nhits_laser_all = 0;
    int32_t nhits_neutrino_all = 0;
      
    //loop on lg
    for(int32_t ch = 1; ch < (constants::laben::channels + 1); ch++) {
      if((bx_dbi::get ()->get_channel (ch).is_ordinary ()) ) {
	int32_t lg_is_ok = 1;
	int32_t lg_is_dead = 0;
	  
	//random triggers (rate in kHz)
	double dark_rate_random = nhits_random[ch - 1] / (tot_random_time * 1e-9) / 1e3;   
	//double error_dark_rate_random = dark_rate_random / sqrt(nhits_random[ch - 1]);
	ID_dark_rate_vs_channel->SetBinContent(ch, 1, dark_rate_random);
	nhits_random_all += nhits_random[ch - 1];

          //pulser triggers (rate in kHz)
	double dark_rate_pulser = nhits_pulser[ch - 1] / (tot_pulser_time * 1e-9) / 1e3;   
	//double error_dark_rate_pulser = dark_rate_pulser / sqrt(nhits_pulser[ch - 1]);
	ID_dark_rate_vs_channel->SetBinContent(ch, 2, dark_rate_pulser);
	nhits_pulser_all += nhits_pulser[ch - 1];
	
          //laser triggers (rate in kHz)
	double dark_rate_laser = nhits_laser[ch - 1] / (tot_laser_time * 1e-9) / 1e3;   
	//double error_dark_rate_laser = dark_rate_laser / sqrt(nhits_laser[ch - 1]);
	ID_dark_rate_vs_channel->SetBinContent(ch, 3, dark_rate_laser);
	nhits_laser_all += nhits_laser[ch - 1];
	
          //neutrino 1-cluster  triggers (rate in kHz)
	double dark_rate_neutrino = nhits_neutrino[ch - 1] / (tot_neutrino_time * 1e-9) / 1e3;   
	//double error_dark_rate_neutrino = dark_rate_neutrino / sqrt(nhits_neutrino[ch - 1]);
	ID_dark_rate_vs_channel->SetBinContent(ch, 4, dark_rate_neutrino);
	nhits_neutrino_all += nhits_neutrino[ch - 1];

 	//all  types of triggers together (rate in kHz)
        double error_dark_rate_total;

	int32_t nhits_total = (index_random * nhits_random[ch - 1] + index_pulser * nhits_pulser[ch - 1] + index_laser * nhits_laser[ch - 1] + index_neutrino * nhits_neutrino[ch - 1]);   
	double dark_rate_total = nhits_total / (tot_dark_time * 1e-9) / 1e3;   
	if(dark_rate_total != 0)  {
	  error_dark_rate_total = dark_rate_total / sqrt(nhits_total); 
	} else {
	  error_dark_rate_total = -10 ;
	}

	ID_dark_rate_vs_channel->SetBinContent(ch, 5, dark_rate_total);
	error_dark_rates_total[ch - 1] = error_dark_rate_total;
	dark_rates_total[ch - 1] = dark_rate_total;
	 
	const db_channel_laben &ch_info = dynamic_cast<const db_channel_laben& >(bx_dbi::get ()->get_channel (ch));
	
	  //is the rate too high? a hot PMT? 
	if(dark_rate_total > dark_rate_thresh_high){
	  if ( ch_info.pmt_has_cone () ) count_cone_pmts_hot ++;           
	  if ( !(ch_info.pmt_has_cone ())) count_not_cone_pmts_hot ++;            
	  get_message(bx_message::log) << "PMT in logical channel " << ch << " HOT, dark rate " << dark_rate_total  << " +- " << error_dark_rate_total << " kHz, " << nhits_total << " hits" << dispatch;
	    //set alive
	  alive[ch - 1] = hot;
	  lg_is_ok = 0;
	  n_ok_pmt --; 
	}

	  //is the rate too low?  a dead PMT? 
        if (dark_rate_total < dark_rate_thresh_low){
	    //is it 2sigma below the threshold?
	  double  Nhits_expected =  tot_dark_time / 1E9 * (dark_rate_thresh_low * 1000);
	  if(nhits_total <  (Nhits_expected - 6 * sqrt(Nhits_expected))){
	    if ( ch_info.pmt_has_cone () ) count_cone_pmts_dead ++;           
	    if ( !(ch_info.pmt_has_cone ())) count_not_cone_pmts_dead ++;            
	    get_message(bx_message::log) << "PMT in logical channel " << ch << " DEAD, dark rate " << dark_rate_total << " +- " << error_dark_rate_total << " kHz, " << nhits_total << " hits" << dispatch;
	      //set alive
	    alive[ch - 1] = dead;
	    lg_is_ok = 0;
	    lg_is_dead = 1;
	    n_ok_pmt --; 
	  }
	}
	

	//to sum rates of PMTs 

	  //all ordinary channels
	sum_rates_all_pmts += dark_rate_total;
	if(!lg_is_dead) ID_all_dark_rates->Fill(dark_rate_total); 

	 //if PMT is ok
	if(lg_is_ok){
	  ID_ok_dark_rates->Fill(dark_rate_total); 
	  sum_rates_ok_pmts += dark_rate_total;
	}
      }
    }

    barn_interface::get ()->network_send (ID_dark_rate_vs_channel, this);
                                   
      //mean rate in kHz
    mean_rate =  sum_rates_all_pmts / n_ordinary_lg;
   
      //mean total rate for ok PMTs
    mean_ok_rate =  ID_ok_dark_rates->GetMean ();
    error_mean_ok_rate =  ID_ok_dark_rates->GetRMS ();

 
      //mean rates in kHz, using all ordinary lg for different trigger types
    double rate_diff_limit =  get_parameter ("rate_diff").get_float (); 
    if(tot_random_time && index_random){ //random
      double mean_rate_random = nhits_random_all / (tot_random_time * 1e-9 * n_ordinary_lg * 1e3);   
      //double error_mean_rate_random =  sqrt(nhits_random_all) / (tot_random_time * 1e-9 * n_ordinary_lg * 1e3);
      get_message(bx_message::log) << "Mean_rate_random (kHz) "   <<   mean_rate_random   << dispatch;
      double diff = std::fabs( 2 * (mean_rate - mean_rate_random)/(mean_rate + mean_rate_random));
      if(diff > rate_diff_limit) get_message(bx_message::warn) << "Non-compatibility mean random dark rates: " << mean_rate_random << " and total:" << mean_rate << " +- (stat) " << error_mean_rate << dispatch;
    }

    if(tot_pulser_time && index_pulser){ //pulser
      double mean_rate_pulser = nhits_pulser_all / (tot_pulser_time * 1e-9 * n_ordinary_lg * 1e3);   
      //double error_mean_rate_pulser =  sqrt(nhits_pulser_all) / (tot_pulser_time * 1e-9 * n_ordinary_lg * 1e3);   
      get_message(bx_message::log) << "Mean_rate_pulser (kHz) "   <<   mean_rate_pulser << dispatch;
      double diff = std::fabs( 2 * (mean_rate - mean_rate_pulser)/(mean_rate + mean_rate_pulser));
      if(diff > rate_diff_limit)  get_message(bx_message::warn) << "Non-compatibility mean pulser dark rates: " << mean_rate_pulser << " and total:" << mean_rate << " +- (stat) " << error_mean_rate << dispatch;
    }

    if(tot_laser_time && index_laser){ //timing laser
      double mean_rate_laser = nhits_laser_all / (tot_laser_time * 1e-9 * n_ordinary_lg * 1e3);   
      // double error_mean_rate_laser = sqrt(nhits_laser_all) / (tot_laser_time * 1e-9 * n_ordinary_lg * 1e3);   
      get_message(bx_message::log) << "Mean_rate_laser (kHz) "    <<   mean_rate_laser  << dispatch;
      double diff = std::fabs( 2 * (mean_rate - mean_rate_laser)/(mean_rate + mean_rate_laser));
      if(diff > rate_diff_limit) get_message(bx_message::warn) << "Non-compatibility mean laser dark rates: " << mean_rate_laser << " and total:" << mean_rate << " +- (stat) " << error_mean_rate << dispatch;    
    }

    if(tot_neutrino_time && index_neutrino){ //neutrino
      double mean_rate_neutrino = nhits_neutrino_all / (tot_neutrino_time * 1e-9 * n_ordinary_lg * 1e3);   
      //double error_mean_rate_neutrino = sqrt(nhits_neutrino_all) / (tot_neutrino_time * 1e-9 * n_ordinary_lg *1e3);   
      get_message(bx_message::log) << "Mean_rate_neutrino (kHz) " <<   mean_rate_neutrino << dispatch;
      double diff = std::fabs( 2 * (mean_rate - mean_rate_neutrino)/(mean_rate + mean_rate_neutrino));
      if(diff > rate_diff_limit) get_message(bx_message::warn) << "Non-compatibility mean neutrino dark rates: " << mean_rate_neutrino << " and total: " << mean_rate << " +- (stat) " << error_mean_rate << dispatch;
    }
  
    
    get_message(bx_message::log) << "Detector dark rate per PMT: " << mean_rate << " kHz (considering all ordinary lg) and " << mean_ok_rate << " +- " << error_mean_ok_rate << " KHz without bad PMTS, using trigger types based on index " << index_trg_type << dispatch;
   
    get_message(bx_message::log) << "Number of PMTs in ordinary lg: " << n_ordinary_lg << dispatch;
    get_message(bx_message::log) << "Number of ok PMTs in ordinary lg: " << n_ok_pmt << dispatch;
    get_message(bx_message::log) << "Number of PMTs with cone in ordinary lg identified as hot: " << count_cone_pmts_hot << dispatch;    
    get_message(bx_message::log) << "Number of PMTs without cone in ordinary lg identified as hot: " << count_not_cone_pmts_hot << dispatch;    
    get_message(bx_message::log) << "Number of PMTs with cone in ordinary lg identified as dead: " << count_cone_pmts_dead << dispatch;      
    get_message(bx_message::log) << "Number of PMTs without cone in ordinary lg identified as dead: " << count_not_cone_pmts_dead << dispatch;    

    //writing to visitors and checking changes in dark rates with respect to the latest calibrated run
    db_run& run_info = bx_dbi::get ()->get_run ();
    
    //warning if mean dark rate of the whole detector has changed
    float previous_mean_dark_noise = run_info.get_laben_mean_dark_noise ();
    float previous_mean_dark_sigma = run_info.get_laben_mean_dark_sigma ();
    if ( ::fabs(previous_mean_dark_noise - mean_ok_rate) > 3 * ((double) error_mean_ok_rate + (double) previous_mean_dark_sigma) ) 
      get_message(bx_message::warn) << "Mean dark rate has changed from " <<  previous_mean_dark_noise << " +- " <<  previous_mean_dark_sigma << " kHz to " << mean_ok_rate << " +- " <<  error_mean_ok_rate << " kHz " << dispatch;

    if(get_parameter ("db_write").get_bool ()){ 
      //    if(!run_info.is_dark_rates_present () && get_parameter ("db_write").get_bool ()){ 
                  
      //set visitors
      run_info.set_laben_mean_dark_noise (mean_ok_rate, this);
      run_info.set_laben_mean_dark_sigma (error_mean_ok_rate, this);
      run_info.set_laben_dead_cone (count_cone_pmts_dead, this);
      run_info.set_laben_dead_no_cone (count_not_cone_pmts_dead, this);
      run_info.set_laben_hot_cone (count_cone_pmts_hot, this);           
      run_info.set_laben_hot_no_cone (count_not_cone_pmts_hot, this);
      
      // checking individual channels
      // get values from the previous calibration run
      for(int32_t ch = 1; ch < (constants::laben::channels + 1); ch++) {
	float previous_laben_dark_noise = run_info.get_laben_dark_noise (ch) ;
	float previous_laben_dark_sigma = run_info.get_laben_dark_sigma (ch) ;
	
        //compare to present values
	if( ::fabs(previous_laben_dark_noise - dark_rates_total [ch - 1]) > ::fabs( 3 * ((double) previous_mean_dark_sigma + (double) error_dark_rates_total [ch - 1])) ) 
	  get_message(bx_message::warn) << "Dark rate for channel " << ch << " has changed from " <<  previous_laben_dark_noise << " +- " <<  previous_laben_dark_sigma << " kHz to " << dark_rates_total[ch - 1] << " +- " <<  error_dark_rates_total[ch - 1] << " kHz " << dispatch;
	
        //set visitors for individula rates
	run_info.set_laben_dark_noise  (ch, dark_rates_total [ch - 1], this);
	run_info.set_laben_dark_sigma  (ch, error_dark_rates_total [ch - 1], this);
	//compare current and previous status  and set visitors
	const std::string& previous_laben_alive_descr = run_info.get_laben_pmt_status (ch);
	pmt_status previous_laben_alive;
	if (!translation_map.check (previous_laben_alive_descr)) {
	  get_message(bx_message::error) << "unknown pmt status" << previous_laben_alive_descr << " for ch " << ch << dispatch;
	  previous_laben_alive = good;
	} else previous_laben_alive = pmt_status(int32_t(translation_map[previous_laben_alive_descr]) % 10);
	
	if(previous_laben_alive != alive[ch - 1]){
	  get_message(bx_message::warn) << "Dark rate status changed for channel " << ch << " from " <<  previous_laben_alive << " to " << alive[ch - 1] << dispatch;
	  alive [ch - 1] = pmt_status(int32_t(alive [ch - 1]) + 10);
	}
	
	
	run_info.set_laben_pmt_status  (ch, translation_map.rfind (alive[ch - 1])->first, this);
      }
    
        //write to the DB
      run_info.write_laben_dark_rates (true, this);
      get_message(bx_message::log) << "Writing to DB" << dispatch; 
    } 
  }
}


/*
  // Fit function. It accepts a 1D histo (the hit time distribution)
  // and fits it in the allowed interval of the gate with a function y = ax + b.
  // The parameter a should be compatible with 0 (flat distribution),  
  //  while the constant b is related with the dark rate.

float bx_calib_laben_dark_rates::fit_rate (TH1D *histo, int32_t npmts) {	
  
  float rate = 0.;
  float slope = 0.;
  float error_slope = 0.;
  float constant = 0.;
  float error_constant =0;

  TF1 *fit = new TF1("fit","[1]*x+[0]", 200, gate_width - 200);
  fit->SetParameter(1,0.001);
  histo->Fit ("fit","NQR");
  slope = fit->GetParameter (1);
  error_slope = fit->GetParError (1);
  constant = fit->GetParameter (0); 
  error_constant = fit->GetParError (0); 
  
    //warning for non flat distribution
  if((fabs(slope) - 2 * error_slope) > 0.01 ){
    get_message(bx_message::warn) << "Non flat distrib. dark noise hits in gate" << dispatch;
    get_message(bx_message::warn) << "Slope: " << slope << " +-" << error_slope <<  dispatch;
    get_message(bx_message::warn) << "Constant: " << constant << " +- " << error_constant << dispatch;
  }
  else{
      //warning if the constant has too big error more than parammeter from echidna.cfg
    if( (fabs(error_constant / constant)) > 0.1 )
      get_message(bx_message::warn) << "big error of the mean dark rate: " << fabs(error_constant / constant) * 100 << " %" << dispatch;
      // number of bins in the fitted part of the histo
    int32_t nbins = ((int32_t) gate_width - 400) / 50; 
    rate = (constant * (float) nbins) / (count_random_triggers * (gate_width - 400.) * 1e-9 ) / npmts * 1e-3; 
  }
  return rate;
}

float bx_calib_laben_dark_rates::fit_error_rate (TH1D *histo, float rate) {	
  
  float error_rate = 0.;
  float constant = 0.;
  float error_constant =0;

  TF1 *fit = new TF1("fit","[1]*x+[0]", 200, gate_width - 200);
  fit->SetParameter(1,0.001);
  histo->Fit ("fit","NQR");
  constant = fit->GetParameter (0); 
  error_constant = fit->GetParError (0); 
     
  error_rate = rate / constant * error_constant ;
  return error_rate;
}
*/

 
 
 
 



