/* BOREXINO Reconstruction program
 *
 * Author: Livia Ludhova <Livia.Ludhova@mi.infn.it>
 * 	   
 * Maintainer: Livia Ludhova <Livia.Ludhova@mi.infn.it>
 *
 *
 * Implementation of bx_calib_laben_retrigger
 */

#include <cmath>
#include <iostream>
#include <vector>
#include <algorithm>
#include "bx_calib_laben_retrigger.hh"
#include "bx_echidna_event.hh"
#include "barn_interface.hh"
#include "db_run.hh"
#include "db_channel.hh"
#include "db_profile.hh"
#include "bx_dbi.hh"
#include "TH1F.h"
#include "TH2F.h"

// ctor
bx_calib_laben_retrigger::bx_calib_laben_retrigger (): bx_base_module("bx_calib_laben_retrigger", bx_base_module::main_loop) {
  require_event_stage (bx_detector::laben, bx_base_event::clustered);
  require_trigger_type (bx_trigger_event::neutrino);
}


// BEGIN
void bx_calib_laben_retrigger::begin () {
  
  n_sigma = get_parameter ("n_sigma").get_float ();
  factor = get_parameter ("factor").get_float ();

    //histograms
  time_vs_lg = new TH2F("time_vs_lg","time_vs_lg_cl1Event",constants::laben::channels,1.,constants::laben::channels + 1.,1000,0,1000); 
  lg_junk = new TH1F("lg_junk","lg_cl1Event_junk",constants::laben::channels,1.,constants::laben::channels + 1.); 
  lg_retrig = new TH1F("lg_retrig","lg_cl1vent_retrig",constants::laben::channels,1.,constants::laben::channels + 1.); 
    
  barn_interface::get ()->store (barn_interface::file, time_vs_lg, this);
  barn_interface::get ()->store (barn_interface::junk, lg_junk, this);    
  barn_interface::get ()->store (barn_interface::file, lg_retrig, this);
}

//DOIT
bx_echidna_event* bx_calib_laben_retrigger::doit (bx_echidna_event *ev) {
  const bx_laben_event& er = ev->get_laben ();
  
  if( (er.get_nclusters ()) == 1.){ 
    double N_hits = er.get_cluster (0).get_clustered_nhits ();
      
    //loop in hits
    for (int32_t hit = 0 ;hit < N_hits; hit++) {    
      double lg = er.get_cluster(0).get_clustered_hit(hit).get_decoded_hit().get_raw_hit ().get_logical_channel ();
      double time = (er.get_cluster (0).get_clustered_hit(hit).get_time ());
      time_vs_lg->Fill (lg,time);
      lg_junk->Fill (lg);
    }    
    
      //find retriggering lg
    for (int32_t ch = 1 ;ch < (constants::laben::channels + 1); ch++) {    
      double hit_or_no = lg_junk->GetBinContent (ch);
      if (hit_or_no > 1) lg_retrig->Fill (ch);
    }
    
    lg_junk->Reset ();
  }  
  return ev;  
}

//END
void bx_calib_laben_retrigger::end () {

  //find mean number of hits per lg in lg-ditrsibution histo
  double mean_retrig = 0;
  int32_t n_work_lg = 0; 
  
  for (int32_t ch = 1 ;ch < (constants::laben::channels + 1); ch++) {    
    double nhit_retrig_one_ch = lg_retrig->GetBinContent (ch);
    if(bx_dbi::get ()->get_channel (ch).is_ordinary ()) {
      if (nhit_retrig_one_ch) {
	n_work_lg ++;
	mean_retrig += nhit_retrig_one_ch;
      }
    }
  }
  
  mean_retrig /= n_work_lg;
  get_message(bx_message::log) << "mean retriggered nhits/lg " << mean_retrig << dispatch;
 
    //total_ratio = fraction of the number of hits in the problematic retrigger region considering all lg  
    // in fact region 110-140 ns is taken, since it is not affected by retriggering and it gives the upper limit 
    // since Nhit_vs_time deacreases 
  TH1D *time_all_channels = time_vs_lg->ProjectionY ("pro",1,constants::laben::channels);
  double n_bad_hits =  time_all_channels->Integral (110,140);
  double n_good_hits = time_all_channels->Integral (0,30);
  double total_ratio = 0;
  double error_total_ratio = 0;
  if(n_bad_hits == 0) get_message(bx_message::warn) << "Not enough statistics" << dispatch;
  else if (n_good_hits && n_bad_hits > 0) {
    total_ratio = n_bad_hits / n_good_hits; 
    error_total_ratio = sqrt ( pow((total_ratio / sqrt(n_bad_hits)), 2) + pow ((total_ratio/sqrt(n_good_hits)),2) ); 
    get_message(bx_message::log) << "total_ratio" << total_ratio << " +- " << error_total_ratio << dispatch;  
  }

  //now to chcek each ordinary lg
  for (int32_t ch = 1 ;ch < (constants::laben::channels + 1); ch++) {    
    if(bx_dbi::get ()->get_channel (ch).is_ordinary ()) {
      index = 0;   
      
        //count the ratio for each lg, look if it is compatible with the total_ratio
      TH1D *time_one_channel = time_vs_lg->ProjectionY ("pro", ch, ch);
      double n_bad_hits_one_ch =  time_one_channel->Integral (150,210);
      double n_good_hits_one_ch = time_one_channel->Integral (0,30);
      double ratio = -10;
      double error_ratio = -10;
      if(n_good_hits_one_ch && n_bad_hits_one_ch) {
	ratio = (n_bad_hits_one_ch / 2) / n_good_hits_one_ch;
	error_ratio = sqrt ( pow((ratio / sqrt(n_bad_hits_one_ch)), 2) + pow ((ratio/sqrt(n_good_hits_one_ch)),2) ); 
	if (ratio > (total_ratio + n_sigma * error_ratio)){
	  index = index + 1;
	}
      }
	
        //now find lg which are not compatible with the mean_retrig
      double nhit_retrig_one_ch = lg_retrig->GetBinContent (ch);   	   
      if (nhit_retrig_one_ch  > (factor * mean_retrig)){
	index = index + 10;
      }
     
   

       if(index > 0)
	 get_message(bx_message::warn) << "lg " << ch << " index " << index << " ratio " << ratio << " +- " << error_ratio << " ratio deviation (nSig) " << (ratio - total_ratio)/error_ratio << " nhit_retrig/mean " << nhit_retrig_one_ch/mean_retrig << dispatch;	  
    }
  }
}
     
