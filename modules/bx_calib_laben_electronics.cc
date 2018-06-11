/* BOREXINO Reconstruction program
 *
 * Author: Daniela Manuzio <dmanuzio@ge.infn.it>, Maria Elena Monzani <monzani@mi.infn.it>,  Livia Ludhova<livia.ludhova@mi.infn.it>
 * Maintainer: Livia Ludhova<livia.ludhova@mi.infn.it>
 * 
 *
 * $Id: bx_calib_laben_electronics.cc,v 1.46 2013/05/24 13:00:42 ludhova Exp $
 *
 * 
 * Implemenentation of bx_calib_laben_electronics
 * Module that studies the status of the inner detector
 * electronics, looking for bad, problematic or dead channels.
 * Input: the decoded hits
 * Output: some histos and a list of bad channels
 * 
*/
#include <cmath>
#include <iostream>
#include <vector>
#include <algorithm>
#include "bx_calib_laben_electronics.hh"
#include "bx_echidna_event.hh"
#include "bx_dbi.hh"
#include "db_profile.hh"
#include "db_channel.hh"
#include "db_run.hh"
#include "barn_interface.hh"
#include "TH1F.h"
#include "TF1.h"

bx_calib_laben_electronics::bx_calib_laben_electronics() : bx_base_module("bx_calib_laben_electronics", bx_base_module::main_loop) {
  require_event_stage (bx_detector::laben, bx_base_event::decoded);
  require_trigger_type (bx_trigger_event::neutrino);
  require_trigger_type (bx_trigger_event::pulser);
  require_trigger_type (bx_trigger_event::laser394);


  ADC_translation_map ["many_zero"] = many_zero;
  ADC_translation_map ["many_FF"]   = many_FF;
  ADC_translation_map ["too_spread"] = too_spread;
  ADC_translation_map ["shifted_from_mean"] = shifted_from_mean;
  ADC_translation_map ["bad_rms"] = bad_rms;
  ADC_translation_map ["very_small_rms"] = very_small_rms;
  ADC_translation_map ["many_negative_values"] = many_negative_values;
  ADC_translation_map ["low_gain"] = low_gain;
  ADC_translation_map ["high_gain"] = high_gain;

  timing_translation_map ["laser_ref_correl"] = laser_ref_correl;       
  timing_translation_map ["trigger_ref_correl_in_pulser"] = trigger_ref_correl_in_pulser  ;     
  timing_translation_map ["trigger_ref_correl_in_laser"] = trigger_ref_correl_in_laser   ;     
  timing_translation_map ["trigger_ref_correl_in_neutrino"] = trigger_ref_correl_in_neutrino;     
  timing_translation_map ["end_of_gate_correl_in_pulser"] = end_of_gate_correl_in_pulser  ;    
  timing_translation_map ["end_of_gate_correl_in_laser"] = end_of_gate_correl_in_laser   ;    
  timing_translation_map ["end_of_gate_correl_in_neutrino"] = end_of_gate_correl_in_neutrino;    
  timing_translation_map ["bad_timing_shape_in_laser"] = bad_timing_shape_in_laser;  

  multiplicity_translation_map["dead_in_pulser"] = dead_in_pulser;
  multiplicity_translation_map["dead_in_laser"] = dead_in_laser;
  multiplicity_translation_map["dead_in_neutrino"] = dead_in_neutrino;
  multiplicity_translation_map["dead_in_raw"] = dead_in_raw;
  multiplicity_translation_map["low_eff_in_pulser"] = low_eff_in_pulser;	      
  multiplicity_translation_map["low_eff_in_laser"] =  low_eff_in_laser;	     
  multiplicity_translation_map["low_eff_in_neutrino"] =  low_eff_in_neutrino;	     
  multiplicity_translation_map["low_eff_in_raw"] =  low_eff_in_raw;	     
  multiplicity_translation_map["hot_in_pulser"] = hot_in_pulser;	     
  multiplicity_translation_map["hot_in_laser"] = hot_in_laser;	     
  multiplicity_translation_map["hot_in_neutrino"] = hot_in_neutrino;	    
  multiplicity_translation_map["hot_in_raw"] = hot_in_raw;	     
  multiplicity_translation_map["retriggering_in_pulser"] = retriggering_in_pulser;  
  multiplicity_translation_map["retriggering_in_laser"] =  retriggering_in_laser;   
  multiplicity_translation_map["retriggering_in_neutrino"] = retriggering_in_neutrino;
  multiplicity_translation_map["fifo_empty"] = fifo_empty;
  multiplicity_translation_map["fifo_full"] = fifo_full;
  multiplicity_translation_map["loosing_raw_hits_in_pulser"]   = loosing_raw_hits_in_pulser;
  multiplicity_translation_map["loosing_raw_hits_in_laser"]    = loosing_raw_hits_in_laser;
  multiplicity_translation_map["loosing_raw_hits_in_neutrino"] = loosing_raw_hits_in_neutrino;

  trg_names[pulser] = "pulser";
  trg_names[laser]  = "laser";
  trg_names[neutrino] = "neutrino";
    
}
      
//BEGIN
void bx_calib_laben_electronics::begin () {

  // Configuration parameters
  f_dead  = get_parameter ("dead_channel_thresh").get_float ();
  f_low_eff[2] = f_low_eff[1] = get_parameter ("low_eff_thresh_l_n").get_float ();
  f_hot[1] = f_hot[2] = get_parameter ("hot_thresh_l_n").get_float ();
  f_low_eff[0] = get_parameter ("low_eff_thresh_p").get_float ();
  f_hot[0] = get_parameter ("hot_thresh_p").get_float ();
  
  f_retriggering = get_parameter ("retriggering_thresh").get_float ();

  f_zero   = get_parameter ("zero_thresh").get_float ();
  f_0xFF   = get_parameter ("0xFF_thresh").get_float ();
  f_too_spread = get_parameter ("too_spread_thresh").get_float ();
  f_mean_offset   = get_parameter ("mean_offset").get_float ();
  f_rms_ADC   = get_parameter ("rms_ADC").get_float ();
  f_rms_charge   = get_parameter ("rms_charge").get_float ();
  f_negative_charge   = get_parameter ("negative_charge_thresh").get_float ();
    

  fifo_full_vs_lg = new TH1F("fifo_full_vs_lg","raw hits",constants::laben::channels, 1., constants::laben::channels + 1.);  
  fifo_empty_vs_lg = new TH1F("fifo_empty_vs_lg","raw hits",constants::laben::channels, 1., constants::laben::channels + 1.);  

  base_vs_lg = new TH2F ("base_vs_lg","base map, pulser triggers, 1st hit", constants::laben::channels, 1., constants::laben::channels + 1., 256, 0., 256.);
  peak_vs_lg = new TH2F ("peak_vs_lg","peak map, pulser triggers, 1st hit", constants::laben::channels, 1., constants::laben::channels + 1., 256, 0., 256.);
  charge_vs_lg = new TH2F ("charge_vs_lg","charge map, pulser triggers, 1st hit", constants::laben::channels, 1., constants::laben::channels + 1., 512, -256., 256.);

  charge_tt1_vs_lg = new TH2F ("charge_tt1_vs_lg","charge map, neutrino triggers, n_decoded_hits < 100", constants::laben::channels, 1., constants::laben::channels + 1., 512, -256., 256.);
 
  pulser_charge_vs_lg = new TH2F ("pulser_charge_vs_lg","raw charg, all hits, pulser trigger",  constants::laben::channels, 1., constants::laben::channels + 1., 512, -256., 256.);
  laser_charge_vs_lg = new TH2F ("laser_charge_vs_lg","raw charge, all hits, laser trigger", constants::laben::channels, 1., constants::laben::channels + 1., 512, -256., 256.);
  
  trigref = new TH1F ("trigref", "lg of trigref", 20, 0, 20);
  lasref = new TH1F ("lasref", "lg of lasref", 20, 0, 20);

  lasref_correl_vs_lg = new TH2F("lasref_correl_vs_lg","decoded_hits_time - laser_ref",constants::laben::channels, 1., constants::laben::channels + 1., 2000, -500., 1500.);  
  
  trigref_correl_vs_lg[pulser] = new TH2F("trigref_correl_in_pulser_vs_lg","decoded_hits_time - trigger_time",constants::laben::channels, 1., constants::laben::channels + 1., 1000, -500., 500.);  
  trigref_correl_vs_lg[laser] = new TH2F("trigref_correl_in_laser_vs_lg","decoded_hits_time - trigger_time",constants::laben::channels, 1., constants::laben::channels + 1., 1000, -500., 500.);  
  trigref_correl_vs_lg[neutrino] = new TH2F("trigref_correl_in_neutrino_vs_lg","decoded_hits_time - trigger_time",constants::laben::channels, 1., constants::laben::channels + 1., 1000, -500., 500.);  
  end_gate_correl_vs_lg[pulser] = new TH2F("end_gate_correl_in_pulser_vs_lg","decoded_hits_time - end_gate",constants::laben::channels, 1., constants::laben::channels + 1., 1000, -500., 500.);  
  end_gate_correl_vs_lg[laser] = new TH2F("end_gate_correl_in_laser_vs_lg","decoded_hits_time - end_gate",constants::laben::channels, 1., constants::laben::channels + 1., 1000, -500., 500.);  
  end_gate_correl_vs_lg[neutrino] = new TH2F("end_gate_correl_in_neutrino_vs_lg","decoded_hits_time - end_gate",constants::laben::channels, 1., constants::laben::channels + 1., 1000, -500., 500.);  
  
  Nhits_vs_lg[pulser]   = new TH1F("Nhits_pulser_vs_lg","decoded_hits, pulser_trg",constants::laben::channels, 1., constants::laben::channels + 1.);  
  Nhits_vs_lg[laser]    = new TH1F("Nhits_laser_vs_lg","decoded_hits, laser_trg", constants::laben::channels, 1., constants::laben::channels + 1.);  
  Nhits_vs_lg[neutrino] = new TH1F("Nhits_neutrino_vs_lg","decoded_hits, neutrino_trg", constants::laben::channels, 1., constants::laben::channels + 1.);

  raw_Nhits_vs_lg[pulser]   = new TH1F("raw_Nhits_pulser_vs_lg","raw_hits, pulser_trg",constants::laben::channels, 1., constants::laben::channels + 1.);  
  raw_Nhits_vs_lg[laser]    = new TH1F("raw_Nhits_laser_vs_lg","raw_hits, laser_trg", constants::laben::channels, 1., constants::laben::channels + 1.);  
  raw_Nhits_vs_lg[neutrino] = new TH1F("raw_Nhits_neutrino_vs_lg","raw_hits, neutrino_trg", constants::laben::channels, 1., constants::laben::channels + 1.);

  frac_Nhits_vs_lg[pulser]   = new TH1F("frac_Nhits_pulser_vs_lg","Ndec/Nraw (-1 for Nraw = 0, -2 not checked), pulser_trg",constants::laben::channels, 1., constants::laben::channels + 1.);  
  frac_Nhits_vs_lg[laser]    = new TH1F("frac_Nhits_laser_vs_lg","Ndec/Nraw (-1 for Nraw = 0, -2 not checked), laser_trg", constants::laben::channels, 1., constants::laben::channels + 1.);  
  frac_Nhits_vs_lg[neutrino] = new TH1F("frac_Nhits_neutrino_vs_lg","Ndec/Nraw (-1 for Nraw = 0, -2 not checked), neutrino_trg", constants::laben::channels, 1., constants::laben::channels + 1.);


  retrigger_dt[pulser]  = new TH2F ("retrigger_dt_in_pulser","time 2nd_hit - 1st_hit in one lg same event", constants::laben::channels, 1., constants::laben::channels + 1., 2000, 0., 2000.);
  retrigger_dt[laser]  = new TH2F ("retrigger_dt_in_laser","time 2nd_hit - 1st_hit in one lg same event",constants::laben::channels, 1., constants::laben::channels + 1., 2000, 0., 2000.);
  retrigger_dt[neutrino]  = new TH2F ("retrigger_dt_in_neutrino","time 2nd_hit - 1st_hit in one lg same event",constants::laben::channels, 1., constants::laben::channels + 1., 2000, 0., 2000.);

  barn_interface::get ()->store (barn_interface::file, fifo_full_vs_lg, this);
  barn_interface::get ()->store (barn_interface::file, fifo_empty_vs_lg, this);

  barn_interface::get ()->store (barn_interface::file, base_vs_lg, this);
  barn_interface::get ()->store (barn_interface::file, peak_vs_lg, this);
  barn_interface::get ()->store (barn_interface::file, charge_vs_lg, this);

  barn_interface::get ()->store (barn_interface::file, charge_tt1_vs_lg, this);

  barn_interface::get ()->store (barn_interface::file, trigref, this);
  barn_interface::get ()->store (barn_interface::file, lasref, this);
  barn_interface::get ()->store (barn_interface::file, pulser_charge_vs_lg, this);
  barn_interface::get ()->store (barn_interface::file, laser_charge_vs_lg, this);


  barn_interface::get ()->store (barn_interface::file, lasref_correl_vs_lg, this);



  for(int32_t trg = 0; trg < 3; trg++) {
    barn_interface::get ()->store (barn_interface::file, trigref_correl_vs_lg[trg], this);
    barn_interface::get ()->store (barn_interface::file, end_gate_correl_vs_lg[trg], this);
    barn_interface::get ()->store (barn_interface::file, Nhits_vs_lg[trg], this);
    barn_interface::get ()->store (barn_interface::file, raw_Nhits_vs_lg[trg], this);
    barn_interface::get ()->store (barn_interface::file, frac_Nhits_vs_lg[trg], this);
    barn_interface::get ()->store (barn_interface::file, retrigger_dt[trg], this);

    Nevents_type[trg] = 0;

      // Internal vectors and arrays
    hits[trg] = new int32_t[constants::laben::channels];
    std::fill_n (hits[trg], constants::laben::channels, 0);
 
    raw_hits[trg] = new int32_t[constants::laben::channels];
    std::fill_n (raw_hits[trg], constants::laben::channels, 0);
 
  }

  get_message(bx_message::debug) << "begin" << dispatch;
}

//DOIT
bx_echidna_event* bx_calib_laben_electronics::doit (bx_echidna_event *ev) {

  double Nraw_hits = ev->get_laben ().get_raw_nhits ();
  double Ndecoded_hits = ev->get_laben ().get_decoded_nhits ();

  int32_t event_trg_type = -10;
  if (ev->get_trigger ().is_pulser ()) event_trg_type = 0;
  if (ev->get_trigger ().is_laser394 ()) event_trg_type = 1;
  if (ev->get_trigger ().is_neutrino ()) event_trg_type = 2;

    //skip for trigger types
  if (event_trg_type == -10) return ev;

    //just do not consider "muons"
  if(event_trg_type == 2 && (Ndecoded_hits > 1000 || ev->get_trigger().get_btb_inputs() != 0) ) return ev; 

  Nevents_type[event_trg_type]++;
 
    //raw hits - FIFO control for all 3 trigger types (neu, pulser, laser) together
  for (int32_t i = 0; i < ev->get_laben ().get_raw_nhits (); i++) {
    const bx_laben_raw_hit &hit = ev->get_laben ().get_raw_hit (i);
    if (hit.check_flag (bx_laben_raw_hit::fifo_full))  fifo_full_vs_lg->Fill (hit.get_logical_channel ());
    if (hit.check_flag (bx_laben_raw_hit::fifo_empty)) fifo_empty_vs_lg->Fill (hit.get_logical_channel ());
  }

    // position of pulses, time 0 is gate start
  db_run& run_info = bx_dbi::get ()->get_run ();
  double gate_width   = run_info.get_laben_gate_width ();	     
  double trigger_offset  = -1 * run_info.get_laben_gate_start ();
 
    //vectors for the calculation of multiplicity_dt
  std::vector<int32_t> lg_one_event(constants::laben::channels, 0);
  std::vector<double> time_by_lg_one_event(constants::laben::channels, 0);

    //loop in all raw hits 
  for (int32_t i = 0; i < Nraw_hits; i++) {
    int32_t lg = ev->get_laben ().get_raw_hit (i).get_logical_channel ();
    raw_hits[event_trg_type][lg - 1]++;
    raw_Nhits_vs_lg[event_trg_type]->Fill(lg);
  }    

    //loop in all decoded hits 
  for (int32_t i = 0; i < Ndecoded_hits; i++) {

      //lg properties
    int32_t lg = ev->get_laben ().get_decoded_hit (i).get_raw_hit ().get_logical_channel ();
    double charge = ev->get_laben ().get_decoded_hit (i).get_charge_bin ();
    hits[event_trg_type][lg - 1]++;
    lg_one_event[lg - 1]++;
    Nhits_vs_lg[event_trg_type]->Fill(lg);
    
    //time
    float hit_time = ev->get_laben ().get_decoded_hit (i).get_raw_time ();
    float trigger_time = ev->get_laben ().get_trigger_rawt ();
    trigref_correl_vs_lg[event_trg_type]->Fill (lg, hit_time - trigger_time); 
    end_gate_correl_vs_lg[event_trg_type]->Fill (lg, hit_time - trigger_time + trigger_offset - gate_width);
    
    //for laser triggers
    if(event_trg_type == 1){
      float laser_time = ev->get_laben ().get_laser_rawt ();
      lasref_correl_vs_lg->Fill (lg, hit_time - laser_time);
      laser_charge_vs_lg->Fill (lg, charge);
    }
    
    //retrigger
    if (lg_one_event[lg - 1] > 1) {
      double dt = ev->get_laben ().get_decoded_hit (i).get_raw_time () - time_by_lg_one_event[lg - 1];
      retrigger_dt[event_trg_type]->Fill (lg, (float)dt);
    }
    time_by_lg_one_event[lg - 1] = ev->get_laben ().get_decoded_hit (i).get_raw_time ();
    
    //pulser triggers 
    if(event_trg_type == 0){
      pulser_charge_vs_lg->Fill (lg, charge);
    }
    
    //charge values filled only for the first hit, pulser      
    if(event_trg_type == 0 && lg_one_event[lg - 1] == 1){
      double base = ev->get_laben ().get_decoded_hit (i).get_raw_hit ().get_base ();
      double peak = ev->get_laben ().get_decoded_hit (i).get_raw_hit ().get_peak ();
      base_vs_lg->Fill (lg, base);
      peak_vs_lg->Fill (lg, peak);
      //      charge_vs_lg->Fill (lg, peak - base); 
      charge_vs_lg->Fill (lg, charge); 
    }

  }
  
  //charge values filled only for the first hit, neutrino
  if(event_trg_type == 2){ //neutrino event
    if (ev->get_trigger().get_btb_inputs() == 0) { //1 cluster, btb input 0
      if(Ndecoded_hits < 100) {  // less than 100 hits
	  for(int32_t i = 0; i < Ndecoded_hits; i++){
	    int32_t lg = ev->get_laben ().get_decoded_hit (i).get_raw_hit ().get_logical_channel ();
	    double charge = ev->get_laben ().get_decoded_hit (i).get_charge_bin ();
	    charge_tt1_vs_lg->Fill (lg, charge); 
	  }
      }
    }
  }
  
  return ev;  
}

//END
void bx_calib_laben_electronics::end () {

    //vectors for later setting of visitors 
  std::vector<std::vector<ADC_status> > charge_base_status_vec(constants::laben::channels); 
  std::vector<std::vector<ADC_status> > charge_peak_status_vec(constants::laben::channels); 
  std::vector<std::vector<ADC_status> > charge_status_vec(constants::laben::channels); 
  std::vector<std::vector<timing_status> > timing_status_pl_vec(constants::laben::channels); 
  std::vector<std::vector<multiplicity> > multiplicity_pl_vec(constants::laben::channels); 
  std::vector<std::vector<timing_status> > timing_status_n_vec(constants::laben::channels); 
  std::vector<std::vector<multiplicity> > multiplicity_n_vec(constants::laben::channels);
  std::vector<std::vector<timing_status> > prev_timing_status_n_vec(constants::laben::channels); 
  std::vector<std::vector<multiplicity> > prev_multiplicity_n_vec(constants::laben::channels);

 //inizialize vectors
  for(int32_t ch = 0; ch < constants::laben::channels; ch++) {
    charge_base_status_vec.push_back (std::vector<ADC_status>());
    charge_peak_status_vec.push_back (std::vector<ADC_status>());
    charge_status_vec.push_back (std::vector<ADC_status>());

    timing_status_pl_vec.push_back (std::vector<timing_status>());
    multiplicity_pl_vec.push_back (std::vector<multiplicity>());

    timing_status_n_vec.push_back (std::vector<timing_status>());
    multiplicity_n_vec.push_back (std::vector<multiplicity>());

    prev_timing_status_n_vec.push_back (std::vector<timing_status>());
    prev_multiplicity_n_vec.push_back (std::vector<multiplicity>());
  }


    //find Nlg_x from the previous run 
  int32_t prev_Nlg_fifo_empty =0;
  int32_t prev_Nlg_fifo_full =0;
  int32_t prev_Nlg_loosing_raw_hits[3] = {};
        
  int32_t prev_Nlg_laser_ref_correl          = 0;
  int32_t prev_Nlg_trigger_ref_correl[3]     = {};
  int32_t prev_Nlg_end_of_gate_correl[3]     = {};
  int32_t prev_Nlg_bad_timing_shape_in_laser = 0;

  int32_t prev_Nlg_dead[4]         = {};
  int32_t prev_Nlg_low_eff[4]      = {};
  int32_t prev_Nlg_hot[4]          = {};
  int32_t prev_Nlg_retriggering[3] = {};
  
  int32_t prev_Nlg_base_many_zeros        = 0;
  int32_t prev_Nlg_base_many_FF           = 0;
  int32_t prev_Nlg_base_too_spread        = 0;
  int32_t prev_Nlg_base_shifted_from_mean = 0;
  int32_t prev_Nlg_base_bad_rms           = 0;
  int32_t prev_Nlg_base_very_small_rms    = 0;
   
  int32_t prev_Nlg_peak_many_zeros        = 0;
  int32_t prev_Nlg_peak_many_FF           = 0;
  int32_t prev_Nlg_peak_too_spread        = 0;
  int32_t prev_Nlg_peak_shifted_from_mean = 0;
  int32_t prev_Nlg_peak_bad_rms           = 0;
  int32_t prev_Nlg_peak_very_small_rms    = 0;
 
  int32_t prev_Nlg_charge_many_zeros           = 0;
  int32_t prev_Nlg_charge_many_FF              = 0;
  int32_t prev_Nlg_charge_many_negative_values = 0;
  int32_t prev_Nlg_charge_too_spread           = 0;
  int32_t prev_Nlg_charge_shifted_from_mean    = 0;
  int32_t prev_Nlg_charge_bad_rms              = 0;
  int32_t prev_Nlg_charge_low_gain               = 0;
  int32_t prev_Nlg_charge_high_gain               = 0;

  int32_t module_says_DB_write = 1;
             
    //vectors of lg, input 0 = does not have this characteristics, input 1 = has this characteristics, 
  std::vector<int32_t> prev_lg_fifo_empty(constants::laben::channels, 0) ;
  std::vector<int32_t> prev_lg_fifo_full (constants::laben::channels, 0);
  std::vector<int32_t> prev_lg_loosing_raw_hits_p (constants::laben::channels, 0);
  std::vector<int32_t> prev_lg_loosing_raw_hits_l (constants::laben::channels, 0);
  std::vector<int32_t> prev_lg_loosing_raw_hits_n (constants::laben::channels, 0);
  
  std::vector<int32_t> prev_lg_laser_ref_correl(constants::laben::channels, 0);
  std::vector<int32_t> prev_lg_trigger_ref_correl_p(constants::laben::channels, 0);
  std::vector<int32_t> prev_lg_trigger_ref_correl_l(constants::laben::channels, 0);
  std::vector<int32_t> prev_lg_trigger_ref_correl_n(constants::laben::channels, 0);
  std::vector<int32_t> prev_lg_end_of_gate_correl_p(constants::laben::channels, 0);
  std::vector<int32_t> prev_lg_end_of_gate_correl_l(constants::laben::channels, 0);
  std::vector<int32_t> prev_lg_end_of_gate_correl_n(constants::laben::channels, 0);
  std::vector<int32_t> prev_lg_bad_timing_shape_in_laser(constants::laben::channels, 0);

  std::vector<int32_t> prev_lg_dead_p(constants::laben::channels, 0);
  std::vector<int32_t> prev_lg_dead_l(constants::laben::channels, 0);
  std::vector<int32_t> prev_lg_dead_n(constants::laben::channels, 0);
  std::vector<int32_t> prev_lg_dead_raw(constants::laben::channels, 0);
  std::vector<int32_t> prev_lg_low_eff_p(constants::laben::channels, 0);
  std::vector<int32_t> prev_lg_low_eff_l(constants::laben::channels, 0);
  std::vector<int32_t> prev_lg_low_eff_n(constants::laben::channels, 0);
  std::vector<int32_t> prev_lg_low_eff_raw(constants::laben::channels, 0);
  std::vector<int32_t> prev_lg_hot_p(constants::laben::channels, 0);
  std::vector<int32_t> prev_lg_hot_l(constants::laben::channels, 0);
  std::vector<int32_t> prev_lg_hot_n(constants::laben::channels, 0);
  std::vector<int32_t> prev_lg_hot_raw(constants::laben::channels, 0);
  std::vector<int32_t> prev_lg_retriggering_p(constants::laben::channels, 0);
  std::vector<int32_t> prev_lg_retriggering_l(constants::laben::channels, 0);
  std::vector<int32_t> prev_lg_retriggering_n(constants::laben::channels, 0);
  
  std::vector<int32_t> prev_lg_base_many_zeros(constants::laben::channels, 0);
  std::vector<int32_t> prev_lg_base_many_FF(constants::laben::channels, 0);
  std::vector<int32_t> prev_lg_base_too_spread(constants::laben::channels, 0);
  std::vector<int32_t> prev_lg_base_shifted_from_mean(constants::laben::channels, 0);
  std::vector<int32_t> prev_lg_base_bad_rms(constants::laben::channels, 0);
  std::vector<int32_t> prev_lg_base_very_small_rms(constants::laben::channels, 0);
  
  std::vector<int32_t> prev_lg_peak_many_zeros(constants::laben::channels, 0);
  std::vector<int32_t> prev_lg_peak_many_FF(constants::laben::channels, 0);
  std::vector<int32_t> prev_lg_peak_too_spread(constants::laben::channels, 0);
  std::vector<int32_t> prev_lg_peak_shifted_from_mean(constants::laben::channels, 0);
  std::vector<int32_t> prev_lg_peak_bad_rms(constants::laben::channels, 0);
  std::vector<int32_t> prev_lg_peak_very_small_rms(constants::laben::channels, 0);
  
  std::vector<int32_t> prev_lg_charge_many_zeros(constants::laben::channels, 0);
  std::vector<int32_t> prev_lg_charge_many_FF(constants::laben::channels, 0);
  std::vector<int32_t> prev_lg_charge_many_negative_values(constants::laben::channels, 0);
  std::vector<int32_t> prev_lg_charge_too_spread(constants::laben::channels, 0);
  std::vector<int32_t> prev_lg_charge_shifted_from_mean(constants::laben::channels, 0);
  std::vector<int32_t> prev_lg_charge_bad_rms(constants::laben::channels, 0);
  std::vector<int32_t> prev_lg_charge_low_gain(constants::laben::channels, 0);
  std::vector<int32_t> prev_lg_charge_high_gain(constants::laben::channels, 0);


  db_run& run_info_prev = bx_dbi::get ()->get_run ();
  
  for(int32_t ch = 0; ch < constants::laben::channels ; ch++) {
     
    // base
    const std::vector<std::string>& laben_charge_base_status_v = run_info_prev.get_laben_charge_base_status  (ch + 1);
    for(unsigned i = 0; i < laben_charge_base_status_v.size (); i ++ ){
      switch (ADC_translation_map[laben_charge_base_status_v[i]]){
        case many_zero:         prev_Nlg_base_many_zeros ++; prev_lg_base_many_zeros[ch] = 1; break;
	case many_FF:           prev_Nlg_base_many_FF ++; prev_lg_base_many_FF[ch] = 1;  break;
        case too_spread:        prev_Nlg_base_too_spread ++; prev_lg_base_too_spread[ch] = 1; break;
        case shifted_from_mean: prev_Nlg_base_shifted_from_mean ++; prev_lg_base_shifted_from_mean[ch] = 1; break;
        case bad_rms:           prev_Nlg_base_bad_rms ++;  prev_lg_base_bad_rms[ch] = 1; break;
        case very_small_rms:    prev_Nlg_base_very_small_rms ++; prev_lg_base_very_small_rms[ch] = 1; break;
        case many_negative_values: ;
        case low_gain: ;
        case high_gain: ;
      }
    } 

    //peak
    const std::vector<std::string>& laben_charge_peak_status_v = run_info_prev.get_laben_charge_peak_status  (ch + 1);
    for(unsigned i = 0; i < laben_charge_peak_status_v.size (); i ++ ){
      switch (ADC_translation_map[laben_charge_peak_status_v[i]]) {
      case many_zero:         prev_Nlg_peak_many_zeros ++; prev_lg_peak_many_zeros[ch] = 1; break;
      case many_FF:           prev_Nlg_peak_many_FF ++; prev_lg_peak_many_FF[ch] = 1; break;
      case too_spread:        prev_Nlg_peak_too_spread ++; prev_lg_peak_too_spread[ch] = 1; break;
      case shifted_from_mean: prev_Nlg_peak_shifted_from_mean ++; prev_lg_peak_shifted_from_mean[ch] = 1; break;
      case bad_rms:           prev_Nlg_peak_bad_rms ++; prev_lg_peak_bad_rms[ch] = 1; break;
      case very_small_rms:    prev_Nlg_peak_very_small_rms ++;  prev_lg_peak_very_small_rms[ch] = 1; break;
      case many_negative_values: ;
      case low_gain: ;
      case high_gain: ;
      }
    }
    
    //charge
    const std::vector<std::string>& laben_charge_status_v = run_info_prev.get_laben_charge_status  (ch + 1);
    for(unsigned i = 0; i < laben_charge_status_v.size (); i ++ ){
      switch (ADC_translation_map[laben_charge_status_v[i]]) {
      case many_zero:         prev_Nlg_charge_many_zeros ++;  prev_lg_charge_many_zeros[ch] = 1; break;
      case many_FF:           prev_Nlg_charge_many_FF ++; prev_lg_charge_many_FF[ch] = 1; break;
      case too_spread:        prev_Nlg_charge_too_spread ++; prev_lg_charge_too_spread[ch] = 1; break;
      case shifted_from_mean: prev_Nlg_charge_shifted_from_mean ++; prev_lg_charge_shifted_from_mean[ch] = 1; break;
      case bad_rms:           prev_Nlg_charge_bad_rms ++; prev_lg_charge_bad_rms[ch] = 1; break;
      case many_negative_values: prev_Nlg_charge_many_negative_values ++;  prev_lg_charge_many_negative_values[ch] = 1; break;
      case low_gain: prev_Nlg_charge_low_gain ++;  prev_lg_charge_low_gain[ch] = 1; break;
      case high_gain: prev_Nlg_charge_high_gain ++;  prev_lg_charge_high_gain[ch] = 1; break;
      case very_small_rms:    ;
      }
    }
    
    //timing status
    const std::vector<std::string>& laben_timing_status_v = run_info_prev.get_laben_timing_status  (ch + 1);
    for(unsigned i = 0; i < laben_timing_status_v.size (); i ++ ){
      switch (timing_translation_map[laben_timing_status_v[i]]) {
      case laser_ref_correl:                 prev_Nlg_laser_ref_correl ++; break;
      case trigger_ref_correl_in_pulser:     prev_Nlg_trigger_ref_correl[0] ++; prev_lg_trigger_ref_correl_p[ch] = 1; break;
      case trigger_ref_correl_in_laser:      prev_Nlg_trigger_ref_correl[1] ++; prev_lg_trigger_ref_correl_l[ch] = 1; break;
      case trigger_ref_correl_in_neutrino:   prev_Nlg_trigger_ref_correl[2] ++; prev_lg_trigger_ref_correl_n[ch] = 1; prev_timing_status_n_vec[ch].push_back(trigger_ref_correl_in_neutrino); break;
      case end_of_gate_correl_in_pulser:     prev_Nlg_end_of_gate_correl[0] ++; prev_lg_end_of_gate_correl_p[ch] = 1; break; 
      case end_of_gate_correl_in_laser:      prev_Nlg_end_of_gate_correl[1] ++; prev_lg_end_of_gate_correl_l[ch] = 1; break; 
      case end_of_gate_correl_in_neutrino:   prev_Nlg_end_of_gate_correl[2] ++; prev_lg_end_of_gate_correl_n[ch] = 1; prev_timing_status_n_vec[ch].push_back(end_of_gate_correl_in_neutrino); break; 
      case bad_timing_shape_in_laser:        prev_Nlg_bad_timing_shape_in_laser ++; prev_lg_bad_timing_shape_in_laser[ch] = 1; break;
      }
    }
    
    //multiplicity
    const std::vector<std::string>& laben_multiplicity_status_v = run_info_prev.get_laben_multiplicity  (ch + 1);
    for(unsigned i = 0; i < laben_multiplicity_status_v.size (); i ++ ){
      switch (multiplicity_translation_map[laben_multiplicity_status_v[i]]) {
      case dead_in_pulser:            prev_Nlg_dead[0] ++;  prev_lg_dead_p[ch] = 1; break;
      case dead_in_laser:             prev_Nlg_dead[1] ++;  prev_lg_dead_l[ch] = 1; break;
      case dead_in_neutrino:          prev_Nlg_dead[2] ++;  prev_lg_dead_n[ch] = 1; prev_multiplicity_n_vec[ch].push_back(dead_in_neutrino); break;
      case dead_in_raw:               prev_Nlg_dead[3] ++;  prev_lg_dead_raw[ch] = 1; prev_multiplicity_n_vec[ch].push_back(dead_in_raw); break;
      case low_eff_in_pulser:         prev_Nlg_low_eff[0] ++; prev_lg_low_eff_p[ch] = 1; break;
      case low_eff_in_laser:          prev_Nlg_low_eff[1] ++; prev_lg_low_eff_l[ch] = 1; break;
      case low_eff_in_neutrino:       prev_Nlg_low_eff[2] ++; prev_lg_low_eff_n[ch] = 1;  prev_multiplicity_n_vec[ch].push_back(low_eff_in_neutrino); break;
      case low_eff_in_raw:            prev_Nlg_low_eff[3] ++; prev_lg_low_eff_raw[ch] = 1;  prev_multiplicity_n_vec[ch].push_back(low_eff_in_raw); break;
      case hot_in_pulser:             prev_Nlg_hot[0] ++;  prev_lg_hot_p[ch] = 1; break;
      case hot_in_laser:              prev_Nlg_hot[1] ++;  prev_lg_hot_l[ch] = 1; break;
      case hot_in_neutrino:           prev_Nlg_hot[2] ++;  prev_lg_hot_n[ch] = 1;  prev_multiplicity_n_vec[ch].push_back(hot_in_neutrino); break;
      case hot_in_raw:                prev_Nlg_hot[3] ++;  prev_lg_hot_raw[ch] = 1;  prev_multiplicity_n_vec[ch].push_back(hot_in_raw); break;
      case retriggering_in_pulser:    prev_Nlg_retriggering[0] ++; prev_lg_retriggering_p[ch] = 1; break;
      case retriggering_in_laser:     prev_Nlg_retriggering[1] ++; prev_lg_retriggering_l[ch] = 1; break; 
      case retriggering_in_neutrino:  prev_Nlg_retriggering[2] ++; prev_lg_retriggering_n[ch] = 1;  prev_multiplicity_n_vec[ch].push_back(retriggering_in_neutrino); break;
      case fifo_empty:                prev_Nlg_fifo_empty ++; prev_lg_fifo_empty[ch] = 1; break;
      case fifo_full:                 prev_Nlg_fifo_full ++; prev_lg_fifo_full[ch] = 1; break;
      case loosing_raw_hits_in_pulser:   prev_Nlg_loosing_raw_hits[0] ++; prev_lg_loosing_raw_hits_p[ch] = 1; break;
      case loosing_raw_hits_in_laser:    prev_Nlg_loosing_raw_hits[1] ++; prev_lg_loosing_raw_hits_l[ch] = 1; break;
      case loosing_raw_hits_in_neutrino: prev_Nlg_loosing_raw_hits[2] ++; prev_lg_loosing_raw_hits_n[ch] = 1;  prev_multiplicity_n_vec[ch].push_back(loosing_raw_hits_in_neutrino); break;
      }
    }
    
  }
  get_message(bx_message::log) << "prev_Nlg_base_many_zeros        " << prev_Nlg_base_many_zeros          << dispatch;
  get_message(bx_message::log) << "prev_Nlg_base_many_FF           " << prev_Nlg_base_many_FF             << dispatch;
  get_message(bx_message::log) << "prev_Nlg_base_too_spread        " << prev_Nlg_base_too_spread          << dispatch;
  get_message(bx_message::log) << "prev_Nlg_base_shifted_from_mean " << prev_Nlg_base_shifted_from_mean   << dispatch;
  get_message(bx_message::log) << "prev_Nlg_base_bad_rms           " << prev_Nlg_base_bad_rms             << dispatch;
  get_message(bx_message::log) << "prev_Nlg_base_very_small_rms    " << prev_Nlg_base_very_small_rms      << dispatch;
  				   					
  get_message(bx_message::log) << "prev_Nlg_peak_many_zeros        " << prev_Nlg_peak_many_zeros          << dispatch;
  get_message(bx_message::log) << "prev_Nlg_peak_many_FF           " << prev_Nlg_peak_many_FF             << dispatch;
  get_message(bx_message::log) << "prev_Nlg_peak_too_spread        " << prev_Nlg_peak_too_spread          << dispatch;
  get_message(bx_message::log) << "prev_Nlg_peak_shifted_from_mean " << prev_Nlg_peak_shifted_from_mean   << dispatch;
  get_message(bx_message::log) << "prev_Nlg_peak_bad_rms           " << prev_Nlg_peak_bad_rms             << dispatch;
  get_message(bx_message::log) << "prev_Nlg_peak_very_small_rms    " << prev_Nlg_peak_very_small_rms      << dispatch;
  				   
  get_message(bx_message::log) << "prev_Nlg_charge_many_zeros           " << prev_Nlg_charge_many_zeros             << dispatch;
  get_message(bx_message::log) << "prev_Nlg_charge_many_FF              " << prev_Nlg_charge_many_FF                << dispatch;
  get_message(bx_message::log) << "prev_Nlg_charge_many_negative_values " << prev_Nlg_charge_many_negative_values   << dispatch;
  get_message(bx_message::log) << "prev_Nlg_charge_too_spread           " << prev_Nlg_charge_too_spread             << dispatch;
  get_message(bx_message::log) << "prev_Nlg_charge_shifted_from_mean    " << prev_Nlg_charge_shifted_from_mean      << dispatch;
  get_message(bx_message::log) << "prev_Nlg_charge_bad_rms              " << prev_Nlg_charge_bad_rms                << dispatch;
  get_message(bx_message::log) << "prev_Nlg_charge_low_gain             " << prev_Nlg_charge_low_gain               << dispatch;
  get_message(bx_message::log) << "prev_Nlg_charge_high_gain             " << prev_Nlg_charge_high_gain               << dispatch;
  get_message(bx_message::log) << "prev_Nlg_laser_ref_correl           " <<  prev_Nlg_laser_ref_correl                 << dispatch;
  									     
  for(int32_t trg = 0; trg < 3; trg ++)  get_message(bx_message::log) << "prev_Nlg_dead_in_" <<  trg_names[trg]   << " "       << prev_Nlg_dead[trg]  << dispatch;
  for(int32_t trg = 0; trg < 3; trg ++)  get_message(bx_message::log) << "prev_Nlg_low_eff_in_" <<  trg_names[trg]   << " "    << prev_Nlg_low_eff[trg]  << dispatch;
  for(int32_t trg = 0; trg < 3; trg ++)  get_message(bx_message::log) << "prev_Nlg_hot_in_" <<  trg_names[trg]    << " "       << prev_Nlg_hot[trg]  << dispatch;

  get_message(bx_message::log) << "prev_Nlg_dead_in_raw_neutrino " << prev_Nlg_dead[3]  << dispatch;
  get_message(bx_message::log) << "prev_Nlg_low_eff_in_raw_neutrino " << prev_Nlg_low_eff[3]  << dispatch;
  get_message(bx_message::log) << "prev_Nlg_hot_in_raw_neutrino " << prev_Nlg_hot[3]  << dispatch;


  for(int32_t trg = 0; trg < 3; trg ++)  get_message(bx_message::log) << "prev_Nlg_retriggering_in_" <<  trg_names[trg] << " " << prev_Nlg_retriggering[trg]  << dispatch;
  for(int32_t trg = 0; trg < 3; trg ++)  get_message(bx_message::log) << "prev_Nlg_trigger_ref_correl_in_" << trg_names[trg]  << " " << prev_Nlg_trigger_ref_correl[trg]  << dispatch;
  for(int32_t trg = 0; trg < 3; trg ++)  get_message(bx_message::log) << "prev_Nlg_end_of_gate_correl_in_" << trg_names[trg]  << " " << prev_Nlg_end_of_gate_correl[trg]  << dispatch;
								      
  get_message(bx_message::log) << "prev_Nlg_fifo_empty          " <<  prev_Nlg_fifo_empty                 << dispatch; 						      
  get_message(bx_message::log) << "prev_Nlg_fifo_full           " <<  prev_Nlg_fifo_full                  << dispatch; 						      
  for(int32_t trg = 0; trg < 3; trg ++) get_message(bx_message::log) << "prev_Nlg_loosing_raw_hits_in_" << trg_names[trg]  << " " << prev_Nlg_loosing_raw_hits[trg] << dispatch;
  
  //get_message(bx_message::log) << "prev_Nlg_bad_timing_shape_in_laser      " << prev_Nlg_bad_timing_shape_in_laser        << dispatch;
								      
    // current run
  int32_t Nlg_fifo_empty = 0;
  int32_t Nlg_fifo_full = 0;
  int32_t Nlg_loosing_raw_hits[3] = {};
  
  int32_t Nlg_laser_ref_correl          = 0;
  int32_t Nlg_trigger_ref_correl[3]     = {};
  int32_t Nlg_end_of_gate_correl[3]     = {};
  int32_t Nlg_bad_timing_shape_in_laser = 0;
  
  int32_t Nlg_dead[4]         = {};
  int32_t Nlg_low_eff[4]      = {};
  int32_t Nlg_hot[4]          = {};
  int32_t Nlg_retriggering[3] = {};
  
  int32_t Nlg_base_many_zeros        = 0;
  int32_t Nlg_base_many_FF           = 0;
  int32_t Nlg_base_too_spread        = 0;
  int32_t Nlg_base_shifted_from_mean = 0;
  int32_t Nlg_base_bad_rms           = 0;
  int32_t Nlg_base_very_small_rms    = 0;
  
  int32_t Nlg_peak_many_zeros        = 0;
  int32_t Nlg_peak_many_FF           = 0;
  int32_t Nlg_peak_too_spread        = 0;
  int32_t Nlg_peak_shifted_from_mean = 0;
  int32_t Nlg_peak_bad_rms           = 0;
  int32_t Nlg_peak_very_small_rms    = 0;
  
  int32_t Nlg_charge_many_zeros           = 0;
  int32_t Nlg_charge_many_FF              = 0;
  int32_t Nlg_charge_many_negative_values = 0;
  int32_t Nlg_charge_too_spread           = 0;
  int32_t Nlg_charge_shifted_from_mean    = 0;
  int32_t Nlg_charge_bad_rms              = 0;
  int32_t Nlg_charge_low_gain             = 0;
  int32_t Nlg_charge_high_gain             = 0;


    //last two bins of the following vectors are used to store
    //pen-last: number of lg which did loose the characteristics
    //    last: number of lg which did "gain" the characteristics
  std::vector<int32_t> lg_fifo_empty(constants::laben::channels + 2, 0) ;
  std::vector<int32_t> lg_fifo_full (constants::laben::channels + 2, 0);
  std::vector<int32_t> lg_loosing_raw_hits_p (constants::laben::channels + 2, 0);
  std::vector<int32_t> lg_loosing_raw_hits_l (constants::laben::channels + 2, 0);
  std::vector<int32_t> lg_loosing_raw_hits_n (constants::laben::channels + 2, 0);
  
  std::vector<int32_t> lg_laser_ref_correl(constants::laben::channels + 2, 0);
  std::vector<int32_t> lg_trigger_ref_correl_p(constants::laben::channels + 2, 0);
  std::vector<int32_t> lg_trigger_ref_correl_l(constants::laben::channels + 2, 0);
  std::vector<int32_t> lg_trigger_ref_correl_n(constants::laben::channels + 2, 0);
  std::vector<int32_t> lg_end_of_gate_correl_p(constants::laben::channels + 2, 0);
  std::vector<int32_t> lg_end_of_gate_correl_l(constants::laben::channels + 2, 0);
  std::vector<int32_t> lg_end_of_gate_correl_n(constants::laben::channels + 2, 0);
  std::vector<int32_t> lg_bad_timing_shape_in_laser(constants::laben::channels + 2, 0);

  std::vector<int32_t> lg_dead_p(constants::laben::channels + 2, 0);
  std::vector<int32_t> lg_dead_l(constants::laben::channels + 2, 0);
  std::vector<int32_t> lg_dead_n(constants::laben::channels + 2, 0);
  std::vector<int32_t> lg_dead_raw(constants::laben::channels + 2, 0);
  std::vector<int32_t> lg_low_eff_p(constants::laben::channels + 2, 0);
  std::vector<int32_t> lg_low_eff_l(constants::laben::channels + 2, 0);
  std::vector<int32_t> lg_low_eff_n(constants::laben::channels + 2, 0);
  std::vector<int32_t> lg_low_eff_raw(constants::laben::channels + 2, 0);
  std::vector<int32_t> lg_hot_p(constants::laben::channels + 2, 0);
  std::vector<int32_t> lg_hot_l(constants::laben::channels + 2, 0);
  std::vector<int32_t> lg_hot_n(constants::laben::channels + 2, 0);
  std::vector<int32_t> lg_hot_raw(constants::laben::channels + 2, 0);
  std::vector<int32_t> lg_retriggering_p(constants::laben::channels + 2, 0);
  std::vector<int32_t> lg_retriggering_l(constants::laben::channels + 2, 0);
  std::vector<int32_t> lg_retriggering_n(constants::laben::channels + 2, 0);
  
  std::vector<int32_t> lg_base_many_zeros(constants::laben::channels + 2, 0);
  std::vector<int32_t> lg_base_many_FF(constants::laben::channels + 2, 0);
  std::vector<int32_t> lg_base_too_spread(constants::laben::channels + 2, 0);
  std::vector<int32_t> lg_base_shifted_from_mean(constants::laben::channels + 2, 0);
  std::vector<int32_t> lg_base_bad_rms(constants::laben::channels + 2, 0);
  std::vector<int32_t> lg_base_very_small_rms(constants::laben::channels + 2, 0);
  
  std::vector<int32_t> lg_peak_many_zeros(constants::laben::channels + 2, 0);
  std::vector<int32_t> lg_peak_many_FF(constants::laben::channels + 2, 0);
  std::vector<int32_t> lg_peak_too_spread(constants::laben::channels + 2, 0);
  std::vector<int32_t> lg_peak_shifted_from_mean(constants::laben::channels + 2, 0);
  std::vector<int32_t> lg_peak_bad_rms(constants::laben::channels + 2, 0);
  std::vector<int32_t> lg_peak_very_small_rms(constants::laben::channels + 2, 0);
  
  std::vector<int32_t> lg_charge_many_zeros(constants::laben::channels + 2, 0);
  std::vector<int32_t> lg_charge_many_FF(constants::laben::channels + 2, 0);
  std::vector<int32_t> lg_charge_many_negative_values(constants::laben::channels + 2, 0);
  std::vector<int32_t> lg_charge_too_spread(constants::laben::channels + 2, 0);
  std::vector<int32_t> lg_charge_shifted_from_mean(constants::laben::channels + 2, 0);
  std::vector<int32_t> lg_charge_bad_rms(constants::laben::channels + 2, 0);
  std::vector<int32_t> lg_charge_low_gain(constants::laben::channels + 2, 0);
  std::vector<int32_t> lg_charge_high_gain(constants::laben::channels + 2, 0);


  std::vector<double> electronics_eff_list(constants::laben::channels, -10);

  double mean_raw_nhits = 0;
  double rms_raw_nhits = 0;

  double mean_nhits[3] = {};
  double rms_nhits[3] = {};

  double mean_base_ref = 0;
  double sigma_base_ref = 0;
  double mean_peak_ref = 0;
  double sigma_peak_ref = 0;
  double mean_charge_ref = 0;
  double sigma_charge_ref = 0;

  //calculate Ndecoded_hits/Nrawhits
  for(int32_t ch = 1; ch <= constants::laben::channels ; ch++) {
    int32_t lg_to_check = 0;
    if(bx_dbi::get ()->get_channel (ch).is_ordinary()) lg_to_check = 1;
    if(get_parameter ("check_reference_lg").get_int () && bx_dbi::get ()->get_channel (ch).is_laser()) lg_to_check = 1;
    if(get_parameter ("check_reference_lg").get_int () && bx_dbi::get ()->get_channel (ch).is_trigger()) lg_to_check = 1;
    if(lg_to_check){	
      for(int32_t trg = 0;trg < 3; trg ++ ) {
	double frac_Nhits;
	if(raw_hits[trg][ch - 1]) frac_Nhits = (double) hits[trg][ch - 1] / (double) raw_hits[trg][ch - 1];
	else frac_Nhits = -1;
	frac_Nhits_vs_lg[trg]->SetBinContent(ch, frac_Nhits);
	if(frac_Nhits >= 0 && frac_Nhits < get_parameter ("good_hits_fraction").get_float ()){
	  get_message(bx_message::log) << "Lg " << ch << " only " << frac_Nhits * 100 << " % of raw hits (total of " << raw_hits[trg][ch - 1] << ") makes it to decoded hits level in " << Nevents_type[trg] << " " << trg_names[trg]  << " triggers " <<  dispatch;      
	  Nlg_loosing_raw_hits[trg] ++;
	  if(trg == 0)  {
	    lg_loosing_raw_hits_p[ch - 1] = 1;
	    multiplicity_pl_vec[ch - 1].push_back(loosing_raw_hits_in_pulser);
	  }
	  if(trg == 1) {
	    lg_loosing_raw_hits_l[ch - 1] = 1;
	    multiplicity_pl_vec[ch - 1].push_back(loosing_raw_hits_in_laser);
	  }
	  if(trg == 2)  {
	    lg_loosing_raw_hits_n[ch - 1] = 1;
	    multiplicity_n_vec[ch - 1].push_back(loosing_raw_hits_in_neutrino);
	  }
	}
      }
    }
    else {
      for(int32_t trg = 0;trg < 3; trg ++ ) {
	frac_Nhits_vs_lg[trg]->SetBinContent(ch, -2);
	if(raw_hits[trg][ch - 1] && !bx_dbi::get ()->get_channel (ch).is_laser() && !bx_dbi::get ()->get_channel (ch).is_trigger()) 
	  get_message(bx_message::warn) << "Lg " << ch << " has raw hits in trigger type " << trg_names[trg] << " and is not an ordinary, neither reference channel" << dispatch;
      }
    }
  }

    //*************************  
    //claculate mean and rms for fifo empty and full
  //(This algorithm is due to Knuth,[1] who cites Welford.[2] from Wikipedia
  double mean_full = 0;     
  double S_full = 0;  
  double mean_empty = 0;     
  double S_empty = 0;  
  int32_t nlg_to_check = 0;
  int32_t nlg_has_empty = 0;
  int32_t nlg_has_full = 0;

  for(int32_t ch = 1; ch <= constants::laben::channels ; ch++) {
    int32_t lg_to_check = 0;
    if(bx_dbi::get ()->get_channel (ch).is_ordinary()) lg_to_check = 1;
    if(get_parameter ("check_reference_lg").get_int () && bx_dbi::get ()->get_channel (ch).is_laser()) lg_to_check = 1;
    if(get_parameter ("check_reference_lg").get_int () && bx_dbi::get ()->get_channel (ch).is_trigger()) lg_to_check = 1;
    if(lg_to_check){
      nlg_to_check ++;
      double n_full_this_lg  = fifo_full_vs_lg->GetBinContent (ch);
      double n_empty_this_lg = fifo_empty_vs_lg->GetBinContent (ch);
      if (n_full_this_lg) nlg_has_full ++;
      if (n_empty_this_lg) nlg_has_empty ++;
        //full
      double delta_full = n_full_this_lg - mean_full;
      mean_full = mean_full + delta_full / nlg_to_check;
      S_full = S_full + delta_full * (n_full_this_lg - mean_full);
        //empty
      double delta_empty = n_empty_this_lg - mean_empty;
      mean_empty = mean_empty + delta_empty / nlg_to_check;
      S_empty = S_empty + delta_empty * (n_empty_this_lg - mean_empty);
    }
  }

    //fifo_empty mean (or truncated mean if necessary) 
  double  rms_empty = std::sqrt(S_empty / nlg_to_check);
  bool    flag_empty_trunc = 0;
  get_message(bx_message::log) << "Mean fifo empty flag " << mean_empty<< " +- " << rms_empty << " in " << Nevents_type[0] + Nevents_type[1] + Nevents_type[2] <<  " triggers " << dispatch;
    //if rms and sqrt(mean) too different for empty, calculate truncated mean and rms if many lg have this problem!
  if(rms_empty > 5 * sqrt(mean_empty) && nlg_has_empty/nlg_to_check > 0.75){
    get_message(bx_message::info) << "Truncated mean fifo empty is going to be calculated " << dispatch;
    mean_empty = 0;     
    S_empty = 0;
    int32_t nlg_to_check_trunc = 0;
    for(int32_t ch = 1; ch <= constants::laben::channels ; ch++) {
      int32_t lg_to_check = 0;
      if(bx_dbi::get ()->get_channel (ch).is_ordinary()) lg_to_check = 1;
      if(get_parameter ("check_reference_lg").get_int () && bx_dbi::get ()->get_channel (ch).is_laser()) lg_to_check = 1;
      if(get_parameter ("check_reference_lg").get_int () && bx_dbi::get ()->get_channel (ch).is_trigger()) lg_to_check = 1;
      if(lg_to_check){
	double n_empty_this_lg = fifo_empty_vs_lg->GetBinContent (ch);
	//nempty is 1 rms around mean value, use the value to calculate truncated mean and rms
	if(n_empty_this_lg  > (mean_empty - rms_empty) && n_empty_this_lg < (mean_empty + rms_empty)){ 
	  nlg_to_check_trunc ++;
	  double delta_empty = n_empty_this_lg - mean_empty;
	  mean_empty = mean_empty + delta_empty / nlg_to_check_trunc;
	  S_empty = S_empty + delta_empty * ( n_empty_this_lg - mean_empty);
	}
      }
    }
    if(nlg_to_check_trunc){
      rms_empty =  std::sqrt(S_empty / nlg_to_check_trunc);
      get_message(bx_message::log) << "Truncated mean fifo empty flag " << mean_empty<< " +- " << rms_empty << " in " << Nevents_type[0] + Nevents_type[1] + Nevents_type[2] <<  " triggers " << dispatch;
      flag_empty_trunc = 1;
    }
  }
  
    //fifo_full mean (or truncated mean if necessary) 
  double  rms_full = std::sqrt(S_full / nlg_to_check);
  bool    flag_full_trunc = 0;
  get_message(bx_message::log) << "Mean fifo full flag " << mean_full << " +- " << rms_full << " in " << Nevents_type[0] + Nevents_type[1] + Nevents_type[2] <<  " triggers " << dispatch;
    //if rms and sqrt(mean) too different for fifo full, calculate truncated mean and rms
  if(rms_full > 5 * sqrt(mean_full)  && nlg_has_full/nlg_to_check > 0.75){
    get_message(bx_message::info) << "Truncated mean fifo full is going to be calculated " << dispatch;
    mean_full = 0;
    S_full = 0;  
    int32_t nlg_to_check_trunc = 0;
    for(int32_t ch = 1; ch <= constants::laben::channels ; ch++) {
      int32_t lg_to_check = 0;
      if(bx_dbi::get ()->get_channel (ch).is_ordinary()) lg_to_check = 1;
      if(get_parameter ("check_reference_lg").get_int () && bx_dbi::get ()->get_channel (ch).is_laser()) lg_to_check = 1;
      if(get_parameter ("check_reference_lg").get_int () && bx_dbi::get ()->get_channel (ch).is_trigger()) lg_to_check = 1;
      if(lg_to_check){
	double n_full_this_lg = fifo_full_vs_lg->GetBinContent (ch);
	  //nfull is 1 rms around mean value, use the value to calculate truncated mean and rms
	if(n_full_this_lg  > (mean_full - rms_full) && n_full_this_lg < (mean_full + rms_full)){ 
	  nlg_to_check_trunc ++;
	  double delta_full = n_full_this_lg - mean_full;
	  mean_full = mean_full + delta_full / nlg_to_check_trunc;
	  S_full = S_full + delta_full * ( n_full_this_lg - mean_full);
	}
      }
    }
    if(nlg_to_check_trunc){
      rms_full = std::sqrt(S_full / nlg_to_check_trunc);
      get_message(bx_message::log) << "Truncated mean fifo fifo flag " << mean_full << " +- " << rms_full << " in " << Nevents_type[0] + Nevents_type[1] + Nevents_type[2] <<  " triggers " << dispatch;
      flag_full_trunc = 1;
    }
  }
  
    //*************************  
    //claculate mean and rms for Nhits_pulser, laser, neutrino, raw_neutrino excluding dead in pulser lg 
   //(This algorithm is due to Knuth,[1] who cites Welford.[2] from Wikipedia
  for(int32_t trg = 0; trg < 3; trg++) {
    double nlg_to_check = 0; 
    double mean = 0;     
    double S = 0;  
    if(Nevents_type[trg]){
      for(int32_t ch = 1; ch <= constants::laben::channels ; ch++) {
	int32_t lg_to_check = 0;
	if(bx_dbi::get ()->get_channel (ch).is_ordinary()) lg_to_check = 1; //ref channels NOT used to calculate mean
	if(lg_to_check){
	  double eff_pulser = hits[0][ch - 1]/ (float)Nevents_type[0] ;
          //not dead_in_pulser
	  if(eff_pulser > f_dead ) {
	    nlg_to_check ++;
	    double delta  = hits[trg][ch - 1] - mean;
	    mean = mean + delta / nlg_to_check;
	    S = S + delta * (hits[trg][ch - 1] - mean);
  	  }
	}
      }
      mean_nhits[trg] = mean;
      if(nlg_to_check) {
	  //calculate rms
	rms_nhits[trg] = std::sqrt(S / nlg_to_check);
	get_message(bx_message::log) << "Mean nhits per lg, trg type  " << trg_names[trg] << ":  " << mean_nhits[trg] << " +- " << rms_nhits[trg] << " in " << Nevents_type[trg] <<  " triggers " << dispatch;
	//check if the pulser mean is compatible with number of pulser triggers
	if(trg == 0 && (mean_nhits[0] > 1.2 * Nevents_type[0] || mean_nhits[0] < 0.8 * Nevents_type[0]) ) 
	  get_message(bx_message::info) << "Mean nhits per lg "<< mean_nhits[0] << " +- " << rms_nhits[0] << " in pulser differs too much from the number of pulser triggers: " << Nevents_type[0] << dispatch;
	//check if the pulser mean is compatible with number of pulser triggers
	if(trg == 1 && mean_nhits[1] < 0.05 * Nevents_type[1]) 
	  get_message(bx_message::warn) << "Mean nhits per lg "<< mean_nhits[1] << " +- " << rms_nhits[1] << " in laser is too low for the number of laser triggers: " << Nevents_type[1] << dispatch;
	if(trg == 1 && mean_nhits[1] > 0.5 * Nevents_type[1])  
	  get_message(bx_message::warn) << "Mean nhits per lg "<< mean_nhits[1] << " +- " << rms_nhits[1] << " in laser is too high for the number of laser triggers: " << Nevents_type[1] << dispatch;
       
	  //if rms and sqrt(mean) too different, calculate truncated mean and rms
	if(rms_nhits[trg] > 5 * sqrt(mean_nhits[trg])){
	  get_message(bx_message::info) << "Nhits vs lg in " <<  trg_names[trg] << " fluctuates too much, stand. dev " << rms_nhits[trg] / sqrt(mean_nhits[trg])<< " times bigger than sqrt(mean), truncated mean and rms are going to be calculated" << dispatch;
	  mean = 0;     
	  S = 0; 
	  nlg_to_check = 0;
	  for(int32_t ch = 1; ch <= constants::laben::channels ; ch++) {
	    int32_t lg_to_check = 0;
	    if(bx_dbi::get ()->get_channel (ch).is_ordinary()) lg_to_check = 1;  //ref channels NOT used to calculate mean
	    if(lg_to_check){
	      double eff_pulser = hits[0][ch - 1]/ (float)Nevents_type[0] ;
	        //not dead_in_pulser and nhits is 1 rms around mean value, use the value to calculate truncated mean and rms
	      if(eff_pulser > f_dead && hits[trg][ch - 1] > (mean_nhits[trg] - rms_nhits[trg]) && hits[trg][ch - 1] < (mean_nhits[trg] + rms_nhits[trg]) ) {
	      nlg_to_check ++;
	      double delta  = hits[trg][ch - 1] - mean;
	      mean = mean + delta / nlg_to_check;
	      S = S + delta * (hits[trg][ch - 1] - mean);
	      }
	    }
	  }
	  mean_nhits[trg] = mean;
	  rms_nhits[trg] =  std::sqrt(S / nlg_to_check);
	  get_message(bx_message::log) << "Truncated mean nhits per lg, trg type  " << trg_names[trg] << ":  " << mean_nhits[trg] << " +- " << rms_nhits[trg] << " in " << Nevents_type[trg] <<  " triggers " << dispatch;     	
	}
      }
      else { 
	get_message(bx_message::error) << "All lg appear dead in pulser" << dispatch;  
	module_says_DB_write = 0;
      } 
    }
  }  

//*************************  
    //claculate mean and rms for Nhits_raw_neutrino, raw_neutrino excluding dead in pulser lg 
   //(This algorithm is due to Knuth,[1] who cites Welford.[2] from Wikipedia
  if(Nevents_type[2]){//if neutrino events present
    double nlg_to_check = 0; 
    double mean = 0;     
    double S = 0;  
    
    for(int32_t ch = 1; ch <= constants::laben::channels ; ch++) {
      int32_t lg_to_check = 0;
      if(bx_dbi::get ()->get_channel (ch).is_ordinary()) lg_to_check = 1; //ref channels NOT used to calculate mean
      if(lg_to_check){
	double eff_pulser = hits[0][ch - 1]/ (float)Nevents_type[0] ;
	//not dead_in_pulser
	if(eff_pulser > f_dead ) {
	  nlg_to_check ++;
	  double delta  = raw_hits[2][ch - 1] - mean;
	  mean = mean + delta / nlg_to_check;
	  S = S + delta * (raw_hits[2][ch - 1] - mean);
	}
      }
    }
    mean_raw_nhits = mean;
    if(nlg_to_check) {
        //calculate rms
      rms_raw_nhits = std::sqrt(S / nlg_to_check);
      get_message(bx_message::log) << "Mean RAW nhits per lg, neutrino trg type:  " << mean_raw_nhits << " +- " << rms_raw_nhits << " in " << Nevents_type[2] <<  " neutrino triggers " << dispatch;
        //if rms and sqrt(mean) too different, calculate truncated mean and rms
      if(rms_raw_nhits > 5 * sqrt(mean_raw_nhits) ){
	get_message(bx_message::info) << "Raw Nhits vs lg in neutrino fluctuates too much, stand. dev " << rms_raw_nhits/ sqrt(mean_raw_nhits) << " times bigger than sqrt(mean), truncated mean and rms are going to be calculated" << dispatch;
	mean = 0;     
	S = 0; 
	nlg_to_check = 0;
	for(int32_t ch = 1; ch <= constants::laben::channels ; ch++) {
	  int32_t lg_to_check = 0;
	  if(bx_dbi::get ()->get_channel (ch).is_ordinary()) lg_to_check = 1;  //ref channels NOT used to calculate mean
	  if(lg_to_check){
	    double eff_pulser = hits[0][ch - 1]/ (float)Nevents_type[0] ;
	      //not dead_in_pulser and nhits is 1 rms around mean value, use the value to calculate truncated mean and rms
	    if(eff_pulser > f_dead && raw_hits[2][ch - 1] > (mean_raw_nhits - rms_raw_nhits) && raw_hits[2][ch - 1] < (mean_raw_nhits + rms_raw_nhits) ) {
	      nlg_to_check ++;
	      double delta  = raw_hits[2][ch - 1] - mean;
	      mean = mean + delta / nlg_to_check;
	      S = S + delta * (raw_hits[2][ch - 1] - mean);
	    }
	  }
	}
	mean_raw_nhits = mean;
	rms_raw_nhits =  std::sqrt(S / nlg_to_check);
	get_message(bx_message::log) << "Truncated mean RAW nhits per lg, neutrino trg:  " << mean_raw_nhits << " +- " << rms_raw_nhits << " in " << Nevents_type[2] <<  " neutrino triggers " << dispatch;     	
      }
    }
    else { 
      get_message(bx_message::error) << "All lg appear dead in pulser" << dispatch;  
      module_says_DB_write = 0;
    } 
  }
  

    //*************************  
    //mean values for charge parameters in only pulser triggers
  if(Nevents_type[0]){
  
    // Mean and rms for the base global distribution
    TH1D *proj1 = base_vs_lg->ProjectionY ("proj1");
    proj1->SetAxisRange (1,254);
    float max_base = proj1->GetBinCenter (proj1->GetMaximumBin ());
    proj1->SetAxisRange (max_base - 40, max_base + 40);
    mean_base_ref =  proj1->GetMean ();
    sigma_base_ref = proj1->GetRMS ();
    proj1->Delete ();
    get_message(bx_message::log) << "Mean base : " << mean_base_ref << " +- " << sigma_base_ref << dispatch;
    
    // Mean and rms for the peak global distribution
    TH1D *proj2 = peak_vs_lg->ProjectionY ("proj2");
    proj2->SetAxisRange (1,254);
    float max_peak = proj2->GetBinCenter (proj2->GetMaximumBin ());
    proj2->SetAxisRange (max_peak - 40, max_peak + 40);
    mean_peak_ref =  proj2->GetMean ();
    sigma_peak_ref = proj2->GetRMS ();
    proj2->Delete ();
    get_message(bx_message::log) << "Mean peak : " << mean_peak_ref << " +- " << sigma_peak_ref << dispatch;
    
        // Mean and rms for the charge global distribution
    TH1D *proj3 = charge_vs_lg->ProjectionY ("proj3");
    float max_charge = proj3->GetBinCenter (proj3->GetMaximumBin ());
    proj3->SetAxisRange (max_charge - 40, max_charge + 40);
    mean_charge_ref =  proj3->GetMean ();
    sigma_charge_ref = proj3->GetRMS ();
    proj3->Delete ();
    get_message(bx_message::log) << "Mean charge = peak - base: " << mean_charge_ref << " +- " << sigma_charge_ref << dispatch;
  }      
  
    //*************************  
    // find single electronic channels  with too many fifo empty flags
  for(int32_t ch = 1; ch <= constants::laben::channels; ch++) {
    int32_t lg_to_check = 0;
    if(bx_dbi::get ()->get_channel (ch).is_ordinary()) lg_to_check = 1;
    if(get_parameter ("check_reference_lg").get_int () && bx_dbi::get ()->get_channel (ch).is_laser()) lg_to_check = 1;
    if(get_parameter ("check_reference_lg").get_int () && bx_dbi::get ()->get_channel (ch).is_trigger()) lg_to_check = 1;
    if(lg_to_check){
      double n_empty_this_lg = fifo_empty_vs_lg->GetBinContent (ch);
      int32_t nrms_fifo_empty;
      if (flag_empty_trunc == 1) nrms_fifo_empty = 10;
      if (flag_empty_trunc == 0) nrms_fifo_empty = 3;
      if(n_empty_this_lg > mean_empty + nrms_fifo_empty  * rms_empty){
	Nlg_fifo_empty ++;
	lg_fifo_empty[ch - 1] = 1;
	get_message(bx_message::log) << "Lg " << ch << ": fifo empty flag (raw hits) " <<  n_empty_this_lg << " times in " << Nevents_type[0] + Nevents_type[1] + Nevents_type[2] <<  " triggers " << dispatch;
	multiplicity_pl_vec[ch - 1].push_back(fifo_empty);
      }
    }
  }

    // find single electronic channels  with too many fifo full flags
  for(int32_t ch = 1; ch <= constants::laben::channels; ch++) {
    int32_t lg_to_check = 0;
    if(bx_dbi::get ()->get_channel (ch).is_ordinary()) lg_to_check = 1;
    if(get_parameter ("check_reference_lg").get_int () && bx_dbi::get ()->get_channel (ch).is_laser()) lg_to_check = 1;
    if(get_parameter ("check_reference_lg").get_int () && bx_dbi::get ()->get_channel (ch).is_trigger()) lg_to_check = 1;
    if(lg_to_check){
      double n_full_this_lg = fifo_full_vs_lg->GetBinContent (ch);
      int32_t nrms_fifo_full;
      if (flag_full_trunc == 1) nrms_fifo_full = 10;
      if (flag_full_trunc == 0) nrms_fifo_full = 3;
      if(n_full_this_lg > mean_full + nrms_fifo_full * rms_full){
	Nlg_fifo_full ++;
	lg_fifo_full[ch - 1] = 1;
	get_message(bx_message::log) << "Lg " << ch << ": fifo full flag (raw hits) " <<  n_full_this_lg << " times in " << Nevents_type[0] + Nevents_type[1] + Nevents_type[2] <<  " triggers " << dispatch;
	multiplicity_pl_vec[ch - 1].push_back(fifo_full);
      }
    }
  }
  
    //fill trigref and lasref
  int32_t found_lasref = 0;
  int32_t found_trigref = 0;
  for(int32_t ch = 1; ch <= constants::laben::channels; ch++) {
    if (bx_dbi::get ()->get_channel (ch).is_laser() ) {
       found_lasref ++;
       lasref->SetBinContent(found_lasref,ch);
    }
    if (bx_dbi::get ()->get_channel (ch).is_trigger() ) {
       found_trigref ++;
       trigref->SetBinContent(found_trigref,ch);
    }
  }
  

  // Status of the single electronic channels  1 (dead, hot, low_eff)
  for(int32_t trg = 0; trg < 3; trg++) {
    for(int32_t ch = 1; ch <= constants::laben::channels; ch++) {
      int32_t lg_to_check = 0;
      int32_t is_trigger_ref = 0;
      int32_t is_laser_ref = 0;

      if(bx_dbi::get ()->get_channel (ch).is_ordinary()) lg_to_check = 1;
      if(get_parameter ("check_reference_lg").get_int () && bx_dbi::get ()->get_channel (ch).is_laser()) {
	is_laser_ref = 1;
	lg_to_check = 1;
      }
      if(get_parameter ("check_reference_lg").get_int () && bx_dbi::get ()->get_channel (ch).is_trigger()) {
	is_trigger_ref = 1;
	lg_to_check = 1;
      }
      
       if(lg_to_check){	
	if(Nevents_type[trg] && rms_nhits[trg]){
	  int32_t dead = 0;
	  int32_t hits_in_lg_for_dead = hits[trg][ch - 1];
	    //lasref in neutrino should have 0 hits, so to cheat the check, the value is set to the mean
	  if (is_laser_ref && trg == 2) hits_in_lg_for_dead =  (int32_t) mean_nhits[2];
	  //dead  ?
	  if(hits_in_lg_for_dead <=  f_dead * mean_nhits[trg]){
	    Nlg_dead[trg] ++;
	    dead = 1;
	    if(trg == 0) lg_dead_p[ch - 1] = 1;
	    if(trg == 1) lg_dead_l[ch - 1] = 1;
	    if(trg == 2) lg_dead_n[ch - 1] = 1;
	    get_message(bx_message::log) << "Lg " << ch << ": DEAD_in_" << trg_names[trg] << "  (" << hits[trg][ch - 1] << " in " << Nevents_type[trg] << "  triggers )" << dispatch;
	    if(trg == 0) multiplicity_pl_vec[ch - 1].push_back(dead_in_pulser);
	    if(trg == 1) multiplicity_pl_vec[ch - 1].push_back(dead_in_laser);
	    if(trg == 2) multiplicity_n_vec[ch - 1].push_back(dead_in_neutrino);
	  }
	    //is hot?
	  int32_t hits_in_lg_for_hot = hits[trg][ch - 1];
	    //for ref channels remove refence hits
	  if (is_trigger_ref) hits_in_lg_for_hot = hits_in_lg_for_hot - Nevents_type[trg];
	  if (is_laser_ref && trg == 1) hits_in_lg_for_hot = hits_in_lg_for_hot - Nevents_type[1];
	  if(hits_in_lg_for_hot > (mean_nhits[trg] + f_hot[trg] * rms_nhits[trg])){
	    Nlg_hot[trg] ++;
	    if(trg == 0) lg_hot_p[ch - 1] = 1;
	    if(trg == 1) lg_hot_l[ch - 1] = 1;
	    if(trg == 2) lg_hot_n[ch - 1] = 1;
	    get_message(bx_message::log) << "Lg " << ch << " Hot_in_ " << trg_names[trg] << " " << hits[trg][ch-1] << " hits in " << Nevents_type[trg] << " triggers" << dispatch;
	    if(trg == 0) multiplicity_pl_vec[ch - 1].push_back(hot_in_pulser);
	    if(trg == 1) multiplicity_pl_vec[ch - 1].push_back(hot_in_laser);
	    if(trg == 2) multiplicity_n_vec[ch - 1].push_back(hot_in_neutrino);
	  }
	    // low_eff ?
	  int32_t ref_for_low = (int32_t) mean_nhits[trg];
	    //for ref channels refernce nTriggers
	  if (is_trigger_ref) ref_for_low = Nevents_type[trg];
	  if (is_laser_ref && trg == 1) ref_for_low = Nevents_type[1];
	  if (is_laser_ref && trg != 1) ref_for_low = 0; //cannot be low eff for neutrino and pulser triggers
	  if(dead == 0 && hits[trg][ch - 1] < f_low_eff[trg] * ref_for_low){ 
	    Nlg_low_eff[trg] ++;
	    if(trg == 0) lg_low_eff_p[ch - 1] = 1;
	    if(trg == 1) lg_low_eff_l[ch - 1] = 1;
	    if(trg == 2) lg_low_eff_n[ch - 1] = 1;
	    get_message(bx_message::log) << "Lg " << ch << ": Low_eff_in_ " <<  trg_names[trg] << " " << hits[trg][ch-1] << " hits in " << Nevents_type[trg] << " triggers" << dispatch;
	    if(trg == 0) multiplicity_pl_vec[ch - 1].push_back(low_eff_in_pulser);
	    if(trg == 1) multiplicity_pl_vec[ch - 1].push_back(low_eff_in_laser);
	    if(trg == 2) multiplicity_n_vec[ch - 1].push_back(low_eff_in_neutrino);
	  }
	}		
      }
    }
  }


  // Multiplicity in Raw hits, neutrino triggers: dead, hot, low_eff)
  for(int32_t ch = 1; ch <= constants::laben::channels; ch++) {
    int32_t lg_to_check = 0;
    int32_t is_trigger_ref = 0;
    int32_t is_laser_ref = 0;

    if(bx_dbi::get ()->get_channel (ch).is_ordinary()) lg_to_check = 1;
    if(get_parameter ("check_reference_lg").get_int () && bx_dbi::get ()->get_channel (ch).is_laser()) {
      is_laser_ref = 1;
      lg_to_check = 1;
    }
    if(get_parameter ("check_reference_lg").get_int () && bx_dbi::get ()->get_channel (ch).is_trigger()) {
      is_trigger_ref = 1;
      lg_to_check = 1;
    }
      
    if(lg_to_check){	
      if(Nevents_type[2] && rms_raw_nhits){
	int32_t dead = 0;
	int32_t hits_in_lg_for_dead = raw_hits[2][ch - 1];
	  //lasref in neutrino should have 0 hits, so to cheat the check, the value is set to the mean
	if (is_laser_ref) hits_in_lg_for_dead =  (int32_t) mean_raw_nhits;
	  //dead  ?
	if(hits_in_lg_for_dead <=  f_dead * mean_raw_nhits){
	  Nlg_dead[3] ++;
	  dead = 1;
	  lg_dead_raw[ch - 1] = 1;
	  get_message(bx_message::log) << "Lg " << ch << ": DEAD_in_raw_neutrino (" << raw_hits[2][ch - 1] << " raw hits in " << Nevents_type[2] << "  neutrino triggers )" << dispatch;
	  multiplicity_n_vec[ch - 1].push_back(dead_in_raw);
	}
	  //is hot in raw?
	int32_t hits_in_lg_for_hot = raw_hits[2][ch - 1];
	  //for ref channels remove refence hits
	if (is_trigger_ref) hits_in_lg_for_hot = hits_in_lg_for_hot - Nevents_type[2];
	if(hits_in_lg_for_hot > (mean_raw_nhits + f_hot[2] * rms_raw_nhits)){
	  Nlg_hot[3] ++;
	  lg_hot_raw[ch - 1] = 1;
	  get_message(bx_message::log) << "Lg " << ch << " Hot_in_raw_neutrino (" << raw_hits[2][ch-1] << " hits in " << Nevents_type[2] << " neutrino triggers" << dispatch;
	  multiplicity_n_vec[ch - 1].push_back(hot_in_raw);
	}
	  // low_eff in raw?
	int32_t ref_for_low = (int32_t) mean_raw_nhits;
	  //for ref channels refernce nTriggers
	if (is_trigger_ref) ref_for_low = Nevents_type[2];
	if (is_laser_ref) ref_for_low = 0; //cannot be low eff for neutrino and pulser triggers
	if(dead == 0 && raw_hits[2][ch - 1] < f_low_eff[2] * ref_for_low){ 
	  Nlg_low_eff[3] ++;
	  lg_low_eff_raw[ch - 1] = 1;
	  get_message(bx_message::log) << "Lg " << ch << ": Low_eff_in_raw_neutrino (" <<  raw_hits[2][ch-1] << " hits in " << Nevents_type[2] << " neutrino triggers" << dispatch;
	  multiplicity_n_vec[ch - 1].push_back(low_eff_in_raw);
	}
      }		
    }
  }

  //*************************************low and high gain for ordinary channels and neutrino triggers

  if(  (charge_tt1_vs_lg->Integral () / 2000.)  > 50) {//we have some statistics

    for(int32_t ch = 1; ch <= constants::laben::channels; ch++) {
      int32_t lg_to_check = 0;
      
      if(bx_dbi::get ()->get_channel (ch).is_ordinary()) lg_to_check = 1;
      if(bx_dbi::get ()->get_channel (ch).is_laser() || bx_dbi::get ()->get_channel (ch).is_trigger()) lg_to_check = 0; // do not check charge for reference triggers
      
      if(lg_to_check){	

	TH1D *charge_tt1_one_lg  = charge_tt1_vs_lg->ProjectionY ("charge_tt1_one_lg",ch, ch);

	double Ntotal = charge_tt1_one_lg->Integral ();
	double Nlow = charge_tt1_one_lg->Integral (255,255+10); //0 is 255
	double Nhigh = charge_tt1_one_lg->Integral (255+40,255+255); //0 is 255


	if (Ntotal > 100) {
	  if( (Nlow/Ntotal) > 0.3) {
	
	    get_message(bx_message::log) << "Lg " << ch << "> 30% of hits in 0-10 ADC "  << dispatch;
	    Nlg_charge_low_gain ++;
	    lg_charge_low_gain[ch - 1] = 1;
	    charge_status_vec[ch - 1].push_back(low_gain);
	  }
	  
	  if( (Nhigh/Ntotal) > 0.5) {
	
	    get_message(bx_message::log) << "Lg " << ch << "> 50% of hits in 40-255 ADC"  << dispatch;
	    Nlg_charge_high_gain ++;
	    lg_charge_high_gain[ch - 1] = 1;
	    charge_status_vec[ch - 1].push_back(high_gain);
	  }
	  
	  charge_tt1_one_lg->Delete();
	}
      }
    }
  }
  
    //*************************  
    //Status of the single electronic channels 2 (retriggering, correlations for not dead lg, ADC quality)
  for(int32_t ch = 1; ch <= constants::laben::channels; ch++) {
    int32_t lg_to_check = 0;                 
    if(bx_dbi::get ()->get_channel (ch).is_ordinary()) lg_to_check = 1;
    if(get_parameter ("check_reference_lg").get_int () && bx_dbi::get ()->get_channel (ch).is_laser()) lg_to_check = 1;
    if(get_parameter ("check_reference_lg").get_int () && bx_dbi::get ()->get_channel (ch).is_trigger()) lg_to_check = 1;
  
    if(lg_to_check){	
    
      for(int32_t trg = 0; trg < 3; trg++){  
          //check if triggers of a given type present
	if (Nevents_type[trg]){
	    //if not dead, 
	  if(hits[trg][ch - 1] > f_dead * mean_nhits[trg]){

	      //retriggering?
	    TH1D* retrigger_dt_one_lg   = retrigger_dt[trg]->ProjectionY ("retrigger_dt_one_lg", ch, ch);
	    retrigger_dt_one_lg->SetAxisRange (0, 200); 
	    if((retrigger_dt_one_lg->Integral () / hits[trg][ch - 1]) > f_retriggering) { 
	      Nlg_retriggering[trg] ++;
	      if(trg == 0) lg_retriggering_p[ch - 1] = 1;
	      if(trg == 1) lg_retriggering_l[ch - 1] = 1;
	      if(trg == 2) lg_retriggering_n[ch - 1] = 1;
	      get_message(bx_message::log) << "Lg " << ch << ": Retriggering_in_" << trg_names[trg] << " (" << retrigger_dt_one_lg->Integral () << " second hits from the total of " << hits[trg][ch - 1] << " in " << Nevents_type[trg] << " triggers )" << dispatch;
	      if(trg == 0) multiplicity_pl_vec[ch - 1].push_back(retriggering_in_pulser);
	      if(trg == 1) multiplicity_pl_vec[ch - 1].push_back(retriggering_in_laser);
	      if(trg == 2) multiplicity_n_vec[ch - 1].push_back(retriggering_in_neutrino); 
	    }
	    
	      //correlation with trigger reference channels check just for ordinary channels
	    if (bx_dbi::get ()->get_channel (ch).is_ordinary()) {
	      double mean_nhits_around_trg = trigref_correl_vs_lg[trg]->Integral(1,constants::laben::channels,300,700) / (constants::laben::channels -  Nlg_dead[trg]) ;    
	      TH1D* trigref_correl_one_lg = trigref_correl_vs_lg[trg]->ProjectionY("trigref_correl_one_lg",ch,ch);
	      double nhits_around_trg_one_lg = trigref_correl_one_lg->Integral (300,700);
	      if (nhits_around_trg_one_lg > (mean_nhits_around_trg + 6 * std::max(1, (int32_t) std::sqrt(mean_nhits_around_trg)))) {
	        //almost all hits around trigger reference (+- 100 ns)
		if(trigref_correl_one_lg->Integral(400,600)/trigref_correl_one_lg->Integral(1,1000) > 0.75){
		  Nlg_trigger_ref_correl[trg] ++;
		  if(trg == 0) lg_trigger_ref_correl_p[ch - 1] = 1;
		  if(trg == 1) lg_trigger_ref_correl_l[ch - 1] = 1;
		  if(trg == 2) lg_trigger_ref_correl_n[ch - 1] = 1;
		  get_message(bx_message::log) << "Lg " << ch << " correlated with trigger ref. in " << trg_names[trg] << dispatch;
		  if(trg == 0) timing_status_pl_vec[ch - 1].push_back(trigger_ref_correl_in_pulser);
		  if(trg == 1) timing_status_pl_vec[ch - 1].push_back(trigger_ref_correl_in_laser);
		  if(trg == 2) timing_status_n_vec[ch - 1].push_back(trigger_ref_correl_in_neutrino);
		}
		else{ 
		  trigref_correl_one_lg->SetAxisRange(-200,200);
		  double rms = trigref_correl_one_lg->GetRMS ();
		  double max = trigref_correl_one_lg->GetMaximum ();
		  if( (max > 10 * std::max(1,(int32_t) mean_nhits_around_trg / 400)) && rms > 2 && rms < 100 ){
		    Nlg_trigger_ref_correl[trg] ++;
		    if(trg == 0) lg_trigger_ref_correl_p[ch - 1] = 1;
		    if(trg == 1) lg_trigger_ref_correl_l[ch - 1] = 1;
		    if(trg == 2) lg_trigger_ref_correl_n[ch - 1] = 1;
		    get_message(bx_message::log) << "Lg " << ch << " correlated with trigger ref. in " << trg_names[trg] << dispatch; 
		    if(trg == 0) timing_status_pl_vec[ch - 1].push_back(trigger_ref_correl_in_pulser);
		    if(trg == 1) timing_status_pl_vec[ch - 1].push_back(trigger_ref_correl_in_laser);
		    if(trg == 2) timing_status_n_vec[ch - 1].push_back(trigger_ref_correl_in_neutrino);
		  }
		}
	      }
	    }
	    
	    //correlation with end of gate  
	    double mean_nhits_around_end_gate = end_gate_correl_vs_lg[trg]->Integral(1,constants::laben::channels,300,700) / (constants::laben::channels -  Nlg_dead[trg]);    	  
	    TH1D* end_gate_correl_one_lg = end_gate_correl_vs_lg[trg]->ProjectionY("end_gate_correl_one_lg",ch,ch);
	    double nhits_around_end_gate_one_lg = end_gate_correl_one_lg->Integral (300,700);
	    if (nhits_around_end_gate_one_lg > (mean_nhits_around_end_gate + 6 * std::max(1, (int32_t) std::sqrt(mean_nhits_around_end_gate)))){
	        //almost all hits around end of gate
	      if(end_gate_correl_one_lg->Integral(400,600)/end_gate_correl_one_lg->Integral(1,1000) > 0.75){
		Nlg_end_of_gate_correl[trg] ++;
		if(trg == 0) lg_end_of_gate_correl_p[ch - 1] = 1;
		if(trg == 1) lg_end_of_gate_correl_l[ch - 1] = 1;
		if(trg == 2) lg_end_of_gate_correl_n[ch - 1] = 1;
		get_message(bx_message::log) << "Lg " << ch << " correlated with end_gate in " << trg_names[trg] << dispatch;
		if(trg == 0) timing_status_pl_vec[ch - 1].push_back(end_of_gate_correl_in_pulser);
		if(trg == 1) timing_status_pl_vec[ch - 1].push_back(end_of_gate_correl_in_laser);
		if(trg == 2) timing_status_n_vec[ch - 1].push_back(end_of_gate_correl_in_neutrino);
	      }
	      else { 
		end_gate_correl_one_lg->SetAxisRange(-200,200);
		double rms = end_gate_correl_one_lg->GetRMS ();
		double max = end_gate_correl_one_lg->GetMaximum ();
		if( (max > 10 * std::max(1,(int32_t) mean_nhits_around_end_gate/400)) && rms > 2 && rms < 100 ){
		  int32_t not_flat = 0;
		    //not flat before 0 inside gate ?
		  for(int32_t j = 0; j < 4; j++){
		    double diff = std::fabs(end_gate_correl_one_lg->Integral ((j+1)*100,(j+2)*100) - end_gate_correl_one_lg->Integral (1 + j*100, (j+1)*100 ));
		    double mean = (end_gate_correl_one_lg->Integral ((j+1)*100,(j+2)*100) + end_gate_correl_one_lg->Integral (1 + j*100, (j+1)*100 )) / 2.;
		    if (diff/mean > 1 && mean > 100) not_flat = 1;
		  }
		  if(not_flat == 0){//check not flat after 0 out of gate ?
		    for(int32_t j = 5; j < 9; j++){
		      double diff = std::fabs(end_gate_correl_one_lg->Integral ((j+1)*100,(j+2)*100) - end_gate_correl_one_lg->Integral (1 + j*100, (j+1)*100 ));
		      double mean = (end_gate_correl_one_lg->Integral ((j+1)*100,(j+2)*100) + end_gate_correl_one_lg->Integral (1 + j*100, (j+1)*100 )) / 2.;
		      if (diff/mean > 1 && mean > 100) not_flat = 1;
		    }
		  }
		  if(not_flat == 1){  
		    Nlg_end_of_gate_correl[trg] ++;
		    if(trg == 0) lg_end_of_gate_correl_p[ch - 1] = 1;
		    if(trg == 1) lg_end_of_gate_correl_l[ch - 1] = 1;
		    if(trg == 2) lg_end_of_gate_correl_n[ch - 1] = 1;
		    get_message(bx_message::log) << "Lg " << ch << " correlated with end_gate in " << trg_names[trg] << dispatch;
		    if(trg == 0) timing_status_pl_vec[ch - 1].push_back(end_of_gate_correl_in_pulser);
		    if(trg == 1) timing_status_pl_vec[ch - 1].push_back(end_of_gate_correl_in_laser);
		    if(trg == 2) timing_status_n_vec[ch - 1].push_back(end_of_gate_correl_in_neutrino);
		  }
		}
	      }
	    }
	  }
	}
      }
    
        //if pulser evenst present    
      if(Nevents_type[0]){
	  //if not dead in pulser, ADC channels properties
	if(hits[0][ch - 1] > f_dead * mean_nhits[0]){
	  TH1D *base_one_lg = base_vs_lg->ProjectionY ("base_one_lg",ch, ch);
	  TH1D *peak_one_lg = peak_vs_lg->ProjectionY ("peak_one_lg",ch, ch);
	  TH1D *charge_one_lg  = charge_vs_lg->ProjectionY ("charge_one_lg",ch, ch);
	    
	    //zeros (Intgeral has bin number input) 
	  float zero_base = base_one_lg->Integral (1,1);
	  float zero_peak = peak_one_lg->Integral (1,1);
	  float zero_charge = charge_one_lg->Integral (257,257);
	    
	    //last ADC channel
	  float FF_base = base_one_lg->Integral (256,256);
	  float FF_peak = peak_one_lg->Integral (256,256);
	  float FF_charge = charge_one_lg->Integral (512,512);
	    
	    //total integral including under- and overflows
	  float tot_base = base_one_lg->Integral (0,257);
	  float tot_peak = peak_one_lg->Integral (0,257);
	  float tot_charge = charge_one_lg->Integral (0,513);
	    
	    //negative charge
	  float neg_charge_hits = charge_one_lg->Integral (0,256);
	  
	  if(tot_base != tot_peak)  get_message(bx_message::error) << "Lg " << ch << " different #entries in base and peak histo: " << tot_base << " and " << tot_peak << std::endl;
	    
	    // too many zeros?  
	  if((zero_base / tot_base) > f_zero){
	    get_message(bx_message::log) << "Lg " << ch  << ": base has high 0-value fraction " << (zero_base / tot_base) << dispatch;
	    Nlg_base_many_zeros ++;
	    lg_base_many_zeros[ch - 1] = 1;
	    charge_base_status_vec[ch - 1].push_back(many_zero);
	  }
	  if((zero_peak / tot_peak) > f_zero){
	    get_message(bx_message::log) << "Lg " << ch << " : peak  has high 0-value fraction " << (zero_peak / tot_peak) << dispatch;
	    Nlg_peak_many_zeros ++;
	    lg_peak_many_zeros[ch - 1] = 1;
	    charge_peak_status_vec[ch - 1].push_back(many_zero);
	  }
	  if((zero_charge / tot_charge) > f_zero){
	    get_message(bx_message::log) << "Lg " << ch  << " : charge has high 0-value fraction " << (zero_charge / tot_charge) << dispatch;
	    Nlg_charge_many_zeros ++;
	    lg_charge_many_zeros[ch - 1] = 1;
	    charge_status_vec[ch - 1].push_back(many_zero);
	  }
	  
	    //too many FF ADC values 
	  if( (FF_base / tot_base) > f_0xFF ){
	    get_message(bx_message::log) << "Lg " << ch  << " : base has high 0xFF fraction " << (FF_base / tot_base) << dispatch;
	    Nlg_base_many_FF ++;
	    lg_base_many_FF[ch - 1] = 1;
	    charge_base_status_vec[ch - 1].push_back(many_FF);
	  }
	  if( (FF_peak / tot_peak) > f_0xFF ){
	    get_message(bx_message::log) << "Lg " << ch  << " : peak has high 0xFF fraction " << (FF_peak / tot_peak) << dispatch;
	    Nlg_peak_many_FF ++;
	    lg_peak_many_FF[ch - 1] = 1;
	    charge_peak_status_vec[ch - 1].push_back(many_FF);
	  }
	  if( (FF_charge / tot_charge) > f_0xFF ){
	    get_message(bx_message::log) << "Lg " << ch << " : charge has high 0xFF fraction " << (FF_charge / tot_charge) << dispatch;
	    Nlg_charge_many_FF ++;
	    lg_charge_many_FF[ch - 1] = 1;
	    charge_status_vec[ch - 1].push_back(many_FF);
	  }
	    
	    //too many negative--charge hits?
	  if( (neg_charge_hits / tot_charge) >  f_negative_charge){
	    get_message(bx_message::log) << "Lg " << ch  << " : high negative-charge hits fraction " << ( neg_charge_hits / tot_charge) << dispatch;
	    Nlg_charge_many_negative_values ++;
	    lg_charge_many_negative_values[ch - 1] = 1;
	    charge_status_vec[ch - 1].push_back(many_negative_values);
	  }
	  
	    //strange base distribution if non-0 and non-FF values present		  
	  if (base_one_lg->Integral (2,256) > 0){	     
	     //integral and max bin in region without 0 and FF
	    double integral_base_big = base_one_lg->Integral (2,256);
	    base_one_lg->SetAxisRange (1,255);
	    int32_t max_bin = base_one_lg->GetMaximumBin ();
	    int32_t int_from  = (int32_t) std::max(1.,max_bin - 20.);
	    int32_t int_to = (int32_t) std::min(255.,max_bin + 20.);
	    
	    double integral_base_around_max = base_one_lg->Integral (int_from, int_to);	     
	      
	      //base very spread
	    if(integral_base_around_max / integral_base_big < f_too_spread){
	      get_message(bx_message::log) << "Lg " << ch  << " : base values very spread" << dispatch;
	      Nlg_base_too_spread ++;
	      lg_base_too_spread[ch - 1] = 1;
	      charge_base_status_vec[ch - 1].push_back(too_spread);
	    }
	    else{
	      double max_base = base_one_lg->GetBinCenter (max_bin);
	      base_one_lg->SetAxisRange (std::max(1.,max_base - 20.), std::min(255.,max_base + 20.));	     
	      float mean_base = base_one_lg->GetMean ();
	      float sigma_base = base_one_lg->GetRMS ();		  
	        //base shifted
	      if(fabs(mean_base - mean_base_ref) > f_mean_offset){
		get_message(bx_message::log) << "Lg " << ch  << " : shifted base by " << (mean_base - mean_base_ref) << " ADC channels" << dispatch;
		Nlg_base_shifted_from_mean ++;
		lg_base_shifted_from_mean[ch - 1] = 1;
		charge_base_status_vec[ch - 1].push_back(shifted_from_mean);
	      }
	       //base bad rms
	      if(sigma_base > f_rms_ADC ){
		get_message(bx_message::log) << "Lg " << ch  << " : broad base, rms " << sigma_base << " channels" << dispatch;
		Nlg_base_bad_rms ++;
		lg_base_bad_rms[ch - 1] = 1;
		charge_base_status_vec[ch - 1].push_back(bad_rms);
	      }
	      //base very small rms
	      if(sigma_base < 1.) {
		get_message(bx_message::log) << "Lg " << ch << " : base has very small rms " << sigma_base << " channels" << dispatch;
		Nlg_base_very_small_rms ++;
		lg_base_very_small_rms[ch - 1] = 1;
		charge_base_status_vec[ch - 1].push_back(very_small_rms);
	      }
	    }
	  }
	  
	    //strange peak distribution non-0 and non-FF values present		  		  
	  if(peak_one_lg->Integral (2,256) > 0){
	    double integral_peak_big = base_one_lg->Integral (2,256);
	    peak_one_lg->SetAxisRange (1,255);
	    int32_t max_bin = peak_one_lg->GetMaximumBin ();
	    int32_t int_from  = (int32_t) std::max(1.,max_bin - 20.);
	    int32_t int_to = (int32_t) std::min(255.,max_bin + 20.);
	    double integral_peak_around_max = peak_one_lg->Integral (int_from, int_to);
	    
	      //peak too spread 
	    if(integral_peak_around_max / integral_peak_big < f_too_spread){
	      get_message(bx_message::log) << "Lg " << ch  << " : peak values very spread" << dispatch;
	      Nlg_peak_too_spread ++;
	      lg_peak_too_spread[ch - 1] = 1;
	      charge_peak_status_vec[ch - 1].push_back(too_spread);
	    }	      
	    else{
	      double max_peak = peak_one_lg->GetBinCenter (max_bin);
	      peak_one_lg->SetAxisRange (std::max(1.,max_peak - 20.), std::min(255.,max_peak + 20.));	     
	      float mean_peak = peak_one_lg->GetMean ();
	      float sigma_peak = peak_one_lg->GetRMS ();
	        //peak shifted
	      if(fabs(mean_peak - mean_peak_ref) > f_mean_offset) {
		get_message(bx_message::log) << "Lg " << ch << " : shifted peak by " << (mean_peak - mean_peak_ref) << " ADC channels" << dispatch;	
		Nlg_peak_shifted_from_mean ++;
		lg_peak_shifted_from_mean[ch - 1] = 1;
		charge_peak_status_vec[ch - 1].push_back(shifted_from_mean);
		}
	        //peak bad rms
	      if(sigma_peak > f_rms_ADC ){
		get_message(bx_message::log) << "Lg " << ch << " : broad peak, rms " << sigma_peak << " channels" << dispatch;
		Nlg_peak_bad_rms ++;
		lg_peak_bad_rms[ch - 1] = 1;
		charge_peak_status_vec[ch - 1].push_back(bad_rms);
	      }	  
	        //peak very small rms
	      if(sigma_peak < 1.){
		get_message(bx_message::log) << "Channel " << ch  << " : peak has very small rms " << sigma_peak << " channels" << dispatch;
		Nlg_peak_very_small_rms ++;
		lg_peak_very_small_rms[ch - 1] = 1;
		charge_peak_status_vec[ch - 1].push_back(very_small_rms);
	      }
	    }
	  }
	  
	    //strange charge distribution		  
	  if ( base_one_lg->Integral (2,254) > 0 && peak_one_lg->Integral (2,254) > 0){
	    double integral_charge_big = charge_one_lg->Integral ();
	    int32_t max_bin = charge_one_lg->GetMaximumBin ();
	    double integral_charge_around_max = charge_one_lg->Integral (max_bin - 20, max_bin + 20);	     
	    
	      //charge too spread
	    if(integral_charge_around_max / integral_charge_big < f_too_spread) {
	      get_message(bx_message::log) << "Lg " << ch << " : charge values are very spread" << dispatch;
	      Nlg_charge_too_spread ++;
	      lg_charge_too_spread[ch - 1] = 1;
	      charge_status_vec[ch - 1].push_back(too_spread);
	    }	      
	    else{
	      double max_charge = charge_one_lg->GetBinCenter (max_bin);
	      charge_one_lg->SetAxisRange (std::max(1.,max_charge - 20.), std::min(255.,max_charge + 20.));	     
	      float mean_charge = charge_one_lg->GetMean ();
	      float sigma_charge = charge_one_lg->GetRMS ();
	        //charge shifted from mean
	      if(fabs(mean_charge - mean_charge_ref) > f_mean_offset){
		get_message(bx_message::log) << "Lg " << ch << " : shifted mean charge by " << (mean_charge - mean_charge_ref) << " ADC channels" << dispatch;
		Nlg_charge_shifted_from_mean ++;
		lg_charge_shifted_from_mean[ch - 1] = 1;
		charge_status_vec[ch - 1].push_back(shifted_from_mean);
	      } 
	        //charge bad rms
	      if(sigma_charge > f_rms_charge){
		get_message(bx_message::log) << "Lg " << ch << " : broad resolution, rms " << sigma_charge << " channels" << dispatch;
		Nlg_charge_bad_rms ++;
		lg_charge_bad_rms[ch - 1] = 1;
		charge_status_vec[ch - 1].push_back(bad_rms);
	      }	  
	    }
	  }
	  base_one_lg->Delete ();	  
	  peak_one_lg->Delete ();
	  charge_one_lg->Delete ();
	}
      }  
      
        //for laser triggers and only ordinary channels
      if(Nevents_type[1] && bx_dbi::get ()->get_channel (ch).is_ordinary() ){
	  //correlation with laser reference channels in laser
	double mean_nhits_around_laser = lasref_correl_vs_lg->Integral(1,constants::laben::channels,300,700) / (constants::laben::channels -  Nlg_dead[1]) ;   
	TH1D* lasref_correl_one_lg = lasref_correl_vs_lg->ProjectionY("lasref_correl_one_lg",ch,ch);
	double nhits_around_laser_one_lg = lasref_correl_one_lg->Integral (300,700);
	if (nhits_around_laser_one_lg > (mean_nhits_around_laser + 6 * std::max(1,(int32_t) std::sqrt(mean_nhits_around_laser)))) {
	    //almost all hits around laser ref. +- 100 ns
	  if(lasref_correl_one_lg->Integral(400,600)/lasref_correl_one_lg->Integral(1,1000) > 0.75){
	    Nlg_laser_ref_correl ++;		
	    get_message(bx_message::log) << "Lg " << ch << " correlated with laser ref. in laser, ratio " << dispatch;
	    timing_status_pl_vec[ch - 1].push_back(laser_ref_correl);
	  }
	  else {
	    lasref_correl_one_lg->SetAxisRange(-200,200);
	    double rms = lasref_correl_one_lg->GetRMS ();
	    double max = lasref_correl_one_lg->GetMaximum ();
	    if( (max > 10 * std::max(1,(int32_t) mean_nhits_around_laser/400)) && rms > 2 && rms < 100 ){
	      Nlg_laser_ref_correl ++;
	      get_message(bx_message::log) << "Lg " << ch << " correlated with laser ref. in laser, ratio " << dispatch;
	      timing_status_pl_vec[ch - 1].push_back(laser_ref_correl);
	    }
	  }
	}	
      }
    }
  }

  get_message(bx_message::log) << "Nlg_fifo_empty          " <<  Nlg_fifo_empty                 << dispatch; 						      
  get_message(bx_message::log) << "Nlg_fifo_full           " <<  Nlg_fifo_full                  << dispatch; 						      
  for(int32_t trg = 0; trg < 3; trg ++) get_message(bx_message::log) << "Nlg_loosing_raw_hits_in_" << trg_names[trg]  << " " << Nlg_loosing_raw_hits[trg] << dispatch; 

  get_message(bx_message::debug) << "end" << dispatch;
  get_message(bx_message::log) << "Nlg_base_many_zeros        " << Nlg_base_many_zeros          << dispatch;
  get_message(bx_message::log) << "Nlg_base_many_FF           " << Nlg_base_many_FF             << dispatch;
  get_message(bx_message::log) << "Nlg_base_too_spread        " << Nlg_base_too_spread          << dispatch;
  get_message(bx_message::log) << "Nlg_base_shifted_from_mean " << Nlg_base_shifted_from_mean   << dispatch;
  get_message(bx_message::log) << "Nlg_base_bad_rms           " << Nlg_base_bad_rms             << dispatch;
  get_message(bx_message::log) << "Nlg_base_very_small_rms    " << Nlg_base_very_small_rms      << dispatch;
  
  get_message(bx_message::log) << "Nlg_peak_many_zeros        " << Nlg_peak_many_zeros          << dispatch;
  get_message(bx_message::log) << "Nlg_peak_many_FF           " << Nlg_peak_many_FF             << dispatch;
  get_message(bx_message::log) << "Nlg_peak_too_spread        " << Nlg_peak_too_spread          << dispatch;
  get_message(bx_message::log) << "Nlg_peak_shifted_from_mean " << Nlg_peak_shifted_from_mean   << dispatch;
  get_message(bx_message::log) << "Nlg_peak_bad_rms           " << Nlg_peak_bad_rms             << dispatch;
  get_message(bx_message::log) << "Nlg_peak_very_small_rms    " << Nlg_peak_very_small_rms      << dispatch;
  
  get_message(bx_message::log) << "Nlg_charge_many_zeros           " << Nlg_charge_many_zeros             << dispatch;
  get_message(bx_message::log) << "Nlg_charge_many_FF              " << Nlg_charge_many_FF                << dispatch;
  get_message(bx_message::log) << "Nlg_charge_many_negative_values " << Nlg_charge_many_negative_values   << dispatch;
  get_message(bx_message::log) << "Nlg_charge_too_spread           " << Nlg_charge_too_spread             << dispatch;
  get_message(bx_message::log) << "Nlg_charge_shifted_from_mean    " << Nlg_charge_shifted_from_mean      << dispatch;
  get_message(bx_message::log) << "Nlg_charge_bad_rms              " << Nlg_charge_bad_rms                << dispatch;
  get_message(bx_message::log) << "Nlg_charge_low_gain             " << Nlg_charge_low_gain               << dispatch;
  get_message(bx_message::log) << "Nlg_charge_high_gain            " << Nlg_charge_high_gain               << dispatch;
  
  get_message(bx_message::log) << "Nlg_laser_ref_correl               " << Nlg_laser_ref_correl                 << dispatch;
  
  for(int32_t trg = 0; trg < 3; trg ++)  get_message(bx_message::log) << "Nlg_dead_in_" <<  trg_names[trg] << " " << Nlg_dead[trg]  << dispatch;
  for(int32_t trg = 0; trg < 3; trg ++)  get_message(bx_message::log) << "Nlg_low_eff_in_" <<  trg_names[trg] << " " << Nlg_low_eff[trg]  << dispatch;
  for(int32_t trg = 0; trg < 3; trg ++)  get_message(bx_message::log) << "Nlg_hot_in_" <<  trg_names[trg] << " " << Nlg_hot[trg]  << dispatch;

  get_message(bx_message::log) << "Nlg_dead_in_raw_neutrino " << Nlg_dead[3]  << dispatch;
  get_message(bx_message::log) << "Nlg_low_eff_in_raw_neutrino " <<  Nlg_low_eff[3]  << dispatch;
  get_message(bx_message::log) << "Nlg_hot_in_raw_neutrino " << Nlg_hot[3]  << dispatch;

  for(int32_t trg = 0; trg < 3; trg ++)  get_message(bx_message::log) << "Nlg_retriggering_in_" <<  trg_names[trg]  << " " << Nlg_retriggering[trg]  << dispatch;
  for(int32_t trg = 0; trg < 3; trg ++)  get_message(bx_message::log) << "Nlg_trigger_ref_correl_in_" << trg_names[trg]  << " " << Nlg_trigger_ref_correl[trg]  << dispatch;
  for(int32_t trg = 0; trg < 3; trg ++)  get_message(bx_message::log) << "Nlg_end_of_gate_correl_in_" << trg_names[trg]  << " " << Nlg_end_of_gate_correl[trg]  << dispatch;

  get_message(bx_message::log) << "Nlg_bad_timing_shape_in_laser      " << Nlg_bad_timing_shape_in_laser        << dispatch;
      
  


    // to check how many channels did change status
  
  for(int32_t ch = 1; ch <= constants::laben::channels; ch ++){
        if ((lg_fifo_empty[ch - 1] - prev_lg_fifo_empty[ch - 1]) == -1) lg_fifo_empty[constants::laben::channels] += 1;
    if ((lg_fifo_empty[ch - 1] - prev_lg_fifo_empty[ch - 1]) ==  1) lg_fifo_empty[constants::laben::channels + 1] += 1;
  
    if ((lg_fifo_full[ch - 1] - prev_lg_fifo_full[ch - 1]) == -1) lg_fifo_full[constants::laben::channels] += 1;
    if ((lg_fifo_full[ch - 1] - prev_lg_fifo_full[ch - 1]) ==  1) lg_fifo_full[constants::laben::channels + 1] += 1;
           
    if (( lg_loosing_raw_hits_p[ch - 1] - prev_lg_loosing_raw_hits_p[ch - 1]) == -1)  lg_loosing_raw_hits_p[constants::laben::channels] += 1;
    if (( lg_loosing_raw_hits_p[ch - 1] - prev_lg_loosing_raw_hits_p[ch - 1]) ==  1)  lg_loosing_raw_hits_p[constants::laben::channels + 1] += 1;
       
    if (( lg_loosing_raw_hits_l[ch - 1] - prev_lg_loosing_raw_hits_l[ch - 1]) == -1)  lg_loosing_raw_hits_l[constants::laben::channels] += 1;
    if (( lg_loosing_raw_hits_l[ch - 1] - prev_lg_loosing_raw_hits_l[ch - 1]) ==  1)  lg_loosing_raw_hits_l[constants::laben::channels + 1] += 1;
        
    if (( lg_loosing_raw_hits_n[ch - 1] - prev_lg_loosing_raw_hits_n[ch - 1]) == -1)  lg_loosing_raw_hits_n[constants::laben::channels] += 1;
    if (( lg_loosing_raw_hits_n[ch - 1] - prev_lg_loosing_raw_hits_n[ch - 1]) ==  1)  lg_loosing_raw_hits_n[constants::laben::channels + 1] += 1;
         
    if (( lg_laser_ref_correl[ch - 1] - prev_lg_laser_ref_correl[ch - 1]) == -1)  lg_laser_ref_correl[constants::laben::channels] += 1;
    if (( lg_laser_ref_correl[ch - 1] - prev_lg_laser_ref_correl[ch - 1]) ==  1)  lg_laser_ref_correl[constants::laben::channels + 1] += 1;
       
    if (( lg_trigger_ref_correl_p[ch - 1] - prev_lg_trigger_ref_correl_p[ch - 1]) == -1)  lg_trigger_ref_correl_p[constants::laben::channels] += 1;
    if (( lg_trigger_ref_correl_p[ch - 1] - prev_lg_trigger_ref_correl_p[ch - 1]) ==  1)  lg_trigger_ref_correl_p[constants::laben::channels + 1] += 1;
    
    if (( lg_trigger_ref_correl_l[ch - 1] - prev_lg_trigger_ref_correl_l[ch - 1]) == -1)  lg_trigger_ref_correl_l[constants::laben::channels] += 1;
    if (( lg_trigger_ref_correl_l[ch - 1] - prev_lg_trigger_ref_correl_l[ch - 1]) ==  1)  lg_trigger_ref_correl_l[constants::laben::channels + 1] += 1;
      
    if (( lg_trigger_ref_correl_n[ch - 1] - prev_lg_trigger_ref_correl_n[ch - 1]) == -1)  lg_trigger_ref_correl_n[constants::laben::channels] += 1;
    if (( lg_trigger_ref_correl_n[ch - 1] - prev_lg_trigger_ref_correl_n[ch - 1]) ==  1)  lg_trigger_ref_correl_n[constants::laben::channels + 1] += 1;
       
    if (( lg_end_of_gate_correl_p[ch - 1] - prev_lg_end_of_gate_correl_p[ch - 1]) == -1)  lg_end_of_gate_correl_p[constants::laben::channels] += 1;
    if (( lg_end_of_gate_correl_p[ch - 1] - prev_lg_end_of_gate_correl_p[ch - 1]) ==  1)  lg_end_of_gate_correl_p[constants::laben::channels + 1] += 1;
         
    if (( lg_end_of_gate_correl_l[ch - 1] - prev_lg_end_of_gate_correl_l[ch - 1]) == -1)  lg_end_of_gate_correl_l[constants::laben::channels] += 1;
    if (( lg_end_of_gate_correl_l[ch - 1] - prev_lg_end_of_gate_correl_l[ch - 1]) ==  1)  lg_end_of_gate_correl_l[constants::laben::channels + 1] += 1;
      
    if (( lg_end_of_gate_correl_n[ch - 1] - prev_lg_end_of_gate_correl_n[ch - 1]) == -1)  lg_end_of_gate_correl_n[constants::laben::channels] += 1;
    if (( lg_end_of_gate_correl_n[ch - 1] - prev_lg_end_of_gate_correl_n[ch - 1]) ==  1)  lg_end_of_gate_correl_n[constants::laben::channels + 1] += 1;
       
    if (( lg_bad_timing_shape_in_laser[ch - 1] - prev_lg_bad_timing_shape_in_laser[ch - 1]) == -1)  lg_bad_timing_shape_in_laser[constants::laben::channels] += 1;
    if (( lg_bad_timing_shape_in_laser[ch - 1] - prev_lg_bad_timing_shape_in_laser[ch - 1]) ==  1)  lg_bad_timing_shape_in_laser[constants::laben::channels + 1] += 1;
    
    if (( lg_dead_p[ch - 1] - prev_lg_dead_p[ch - 1]) == -1)  lg_dead_p[constants::laben::channels] += 1;
    if (( lg_dead_p[ch - 1] - prev_lg_dead_p[ch - 1]) ==  1)  lg_dead_p[constants::laben::channels + 1] += 1;
   
    if (( lg_dead_l[ch - 1] - prev_lg_dead_l[ch - 1]) == -1)  lg_dead_l[constants::laben::channels] += 1;
    if (( lg_dead_l[ch - 1] - prev_lg_dead_l[ch - 1]) ==  1)  lg_dead_l[constants::laben::channels + 1] += 1;
       
    if (( lg_dead_n[ch - 1] - prev_lg_dead_n[ch - 1]) == -1)  lg_dead_n[constants::laben::channels] += 1;
    if (( lg_dead_n[ch - 1] - prev_lg_dead_n[ch - 1]) ==  1)  lg_dead_n[constants::laben::channels + 1] += 1;

    if (( lg_dead_raw[ch - 1] - prev_lg_dead_raw[ch - 1]) == -1)  lg_dead_raw[constants::laben::channels] += 1;
    if (( lg_dead_raw[ch - 1] - prev_lg_dead_raw[ch - 1]) ==  1)  lg_dead_raw[constants::laben::channels + 1] += 1;
     
    if (( lg_low_eff_p[ch - 1] - prev_lg_low_eff_p[ch - 1]) == -1)  lg_low_eff_p[constants::laben::channels] += 1;
    if (( lg_low_eff_p[ch - 1] - prev_lg_low_eff_p[ch - 1]) ==  1)  lg_low_eff_p[constants::laben::channels + 1] += 1;
   
    if (( lg_low_eff_l[ch - 1] - prev_lg_low_eff_l[ch - 1]) == -1)  lg_low_eff_l[constants::laben::channels] += 1;
    if (( lg_low_eff_l[ch - 1] - prev_lg_low_eff_l[ch - 1]) ==  1)  lg_low_eff_l[constants::laben::channels + 1] += 1;
        
    if (( lg_low_eff_n[ch - 1] - prev_lg_low_eff_n[ch - 1]) == -1)  lg_low_eff_n[constants::laben::channels] += 1;
    if (( lg_low_eff_n[ch - 1] - prev_lg_low_eff_n[ch - 1]) ==  1)  lg_low_eff_n[constants::laben::channels + 1] += 1;
    
    if (( lg_low_eff_raw[ch - 1] - prev_lg_low_eff_raw[ch - 1]) == -1)  lg_low_eff_raw[constants::laben::channels] += 1;
    if (( lg_low_eff_raw[ch - 1] - prev_lg_low_eff_raw[ch - 1]) ==  1)  lg_low_eff_raw[constants::laben::channels + 1] += 1;
    
    if (( lg_hot_p[ch - 1] - prev_lg_hot_p[ch - 1]) == -1)  lg_hot_p[constants::laben::channels] += 1;
    if (( lg_hot_p[ch - 1] - prev_lg_hot_p[ch - 1]) ==  1)  lg_hot_p[constants::laben::channels + 1] += 1;
        
    if (( lg_hot_l[ch - 1] - prev_lg_hot_l[ch - 1]) == -1)  lg_hot_l[constants::laben::channels] += 1;
    if (( lg_hot_l[ch - 1] - prev_lg_hot_l[ch - 1]) ==  1)  lg_hot_l[constants::laben::channels + 1] += 1;
         
    if (( lg_hot_n[ch - 1] - prev_lg_hot_n[ch - 1]) == -1)  lg_hot_n[constants::laben::channels] += 1;
    if (( lg_hot_n[ch - 1] - prev_lg_hot_n[ch - 1]) ==  1)  lg_hot_n[constants::laben::channels + 1] += 1;
       
    if (( lg_hot_raw[ch - 1] - prev_lg_hot_raw[ch - 1]) == -1)  lg_hot_raw[constants::laben::channels] += 1;
    if (( lg_hot_raw[ch - 1] - prev_lg_hot_raw[ch - 1]) ==  1)  lg_hot_raw[constants::laben::channels + 1] += 1;
       
    if (( lg_retriggering_p[ch - 1] - prev_lg_retriggering_p[ch - 1]) == -1)  lg_retriggering_p[constants::laben::channels] += 1;
    if (( lg_retriggering_p[ch - 1] - prev_lg_retriggering_p[ch - 1]) ==  1)  lg_retriggering_p[constants::laben::channels + 1] += 1;
    
    if (( lg_retriggering_l[ch - 1] - prev_lg_retriggering_l[ch - 1]) == -1)  lg_retriggering_l[constants::laben::channels] += 1;
    if (( lg_retriggering_l[ch - 1] - prev_lg_retriggering_l[ch - 1]) ==  1)  lg_retriggering_l[constants::laben::channels + 1] += 1;
         
    if (( lg_retriggering_n[ch - 1] - prev_lg_retriggering_n[ch - 1]) == -1)  lg_retriggering_n[constants::laben::channels] += 1;
    if (( lg_retriggering_n[ch - 1] - prev_lg_retriggering_n[ch - 1]) ==  1)  lg_retriggering_n[constants::laben::channels + 1] += 1;
    
    if (( lg_base_many_zeros[ch - 1] - prev_lg_base_many_zeros[ch - 1]) == -1)  lg_base_many_zeros[constants::laben::channels] += 1;
    if (( lg_base_many_zeros[ch - 1] - prev_lg_base_many_zeros[ch - 1]) ==  1)  lg_base_many_zeros[constants::laben::channels + 1] += 1;
   
    if (( lg_base_many_FF[ch - 1] - prev_lg_base_many_FF[ch - 1]) == -1)  lg_base_many_FF[constants::laben::channels] += 1;
    if (( lg_base_many_FF[ch - 1] - prev_lg_base_many_FF[ch - 1]) ==  1)  lg_base_many_FF[constants::laben::channels + 1] += 1;
   
    if (( lg_base_too_spread[ch - 1] - prev_lg_base_too_spread[ch - 1]) == -1)  lg_base_too_spread[constants::laben::channels] += 1;
    if (( lg_base_too_spread[ch - 1] - prev_lg_base_too_spread[ch - 1]) ==  1)  lg_base_too_spread[constants::laben::channels + 1] += 1;
    
    if (( lg_base_shifted_from_mean[ch - 1] - prev_lg_base_shifted_from_mean[ch - 1]) == -1)  lg_base_shifted_from_mean[constants::laben::channels] += 1;
    if (( lg_base_shifted_from_mean[ch - 1] - prev_lg_base_shifted_from_mean[ch - 1]) ==  1)  lg_base_shifted_from_mean[constants::laben::channels + 1] += 1;
    
    if (( lg_base_bad_rms[ch - 1] - prev_lg_base_bad_rms[ch - 1]) == -1)  lg_base_bad_rms[constants::laben::channels] += 1;
    if (( lg_base_bad_rms[ch - 1] - prev_lg_base_bad_rms[ch - 1]) ==  1)  lg_base_bad_rms[constants::laben::channels + 1] += 1;
   
    if (( lg_base_very_small_rms[ch - 1] - prev_lg_base_very_small_rms[ch - 1]) == -1)  lg_base_very_small_rms[constants::laben::channels] += 1;
    if (( lg_base_very_small_rms[ch - 1] - prev_lg_base_very_small_rms[ch - 1]) ==  1)  lg_base_very_small_rms[constants::laben::channels + 1] += 1;

    if (( lg_peak_many_zeros[ch - 1] - prev_lg_peak_many_zeros[ch - 1]) == -1)  lg_peak_many_zeros[constants::laben::channels] += 1;
    if (( lg_peak_many_zeros[ch - 1] - prev_lg_peak_many_zeros[ch - 1]) ==  1)  lg_peak_many_zeros[constants::laben::channels + 1] += 1;
    
    if (( lg_peak_many_FF[ch - 1] - prev_lg_peak_many_FF[ch - 1]) == -1)  lg_peak_many_FF[constants::laben::channels] += 1;
    if (( lg_peak_many_FF[ch - 1] - prev_lg_peak_many_FF[ch - 1]) ==  1)  lg_peak_many_FF[constants::laben::channels + 1] += 1;
  
    if (( lg_peak_too_spread[ch - 1] - prev_lg_peak_too_spread[ch - 1]) == -1)  lg_peak_too_spread[constants::laben::channels] += 1;
    if (( lg_peak_too_spread[ch - 1] - prev_lg_peak_too_spread[ch - 1]) ==  1)  lg_peak_too_spread[constants::laben::channels + 1] += 1;
  
    if (( lg_peak_shifted_from_mean[ch - 1] - prev_lg_peak_shifted_from_mean[ch - 1]) == -1)  lg_peak_shifted_from_mean[constants::laben::channels] += 1;
    if (( lg_peak_shifted_from_mean[ch - 1] - prev_lg_peak_shifted_from_mean[ch - 1]) ==  1)  lg_peak_shifted_from_mean[constants::laben::channels + 1] += 1;
    
    if (( lg_peak_bad_rms[ch - 1] - prev_lg_peak_bad_rms[ch - 1]) == -1)  lg_peak_bad_rms[constants::laben::channels] += 1;
    if (( lg_peak_bad_rms[ch - 1] - prev_lg_peak_bad_rms[ch - 1]) ==  1)  lg_peak_bad_rms[constants::laben::channels + 1] += 1;
   
    if (( lg_peak_very_small_rms[ch - 1] - prev_lg_peak_very_small_rms[ch - 1]) == -1)  lg_peak_very_small_rms[constants::laben::channels] += 1;
    if (( lg_peak_very_small_rms[ch - 1] - prev_lg_peak_very_small_rms[ch - 1]) ==  1)  lg_peak_very_small_rms[constants::laben::channels + 1] += 1;
       
    if (( lg_charge_many_zeros[ch - 1] - prev_lg_charge_many_zeros[ch - 1]) == -1)  lg_charge_many_zeros[constants::laben::channels] += 1;
    if (( lg_charge_many_zeros[ch - 1] - prev_lg_charge_many_zeros[ch - 1]) ==  1)  lg_charge_many_zeros[constants::laben::channels + 1] += 1;
   
    if (( lg_charge_many_FF[ch - 1] - prev_lg_charge_many_FF[ch - 1]) == -1)  lg_charge_many_FF[constants::laben::channels] += 1;
    if (( lg_charge_many_FF[ch - 1] - prev_lg_charge_many_FF[ch - 1]) ==  1)  lg_charge_many_FF[constants::laben::channels + 1] += 1;
   
    if (( lg_charge_many_negative_values[ch - 1] - prev_lg_charge_many_negative_values[ch - 1]) == -1)  lg_charge_many_negative_values[constants::laben::channels] += 1;
    if (( lg_charge_many_negative_values[ch - 1] - prev_lg_charge_many_negative_values[ch - 1]) ==  1)  lg_charge_many_negative_values[constants::laben::channels + 1] += 1;

     if (( lg_charge_too_spread[ch - 1] - prev_lg_charge_too_spread[ch - 1]) == -1)  lg_charge_too_spread[constants::laben::channels] += 1;
     if (( lg_charge_too_spread[ch - 1] - prev_lg_charge_too_spread[ch - 1]) ==  1)  lg_charge_too_spread[constants::laben::channels + 1] += 1;
    
     if (( lg_charge_shifted_from_mean[ch - 1] - prev_lg_charge_shifted_from_mean[ch - 1]) == -1)  lg_charge_shifted_from_mean[constants::laben::channels] += 1;
     if (( lg_charge_shifted_from_mean[ch - 1] - prev_lg_charge_shifted_from_mean[ch - 1]) ==  1)  lg_charge_shifted_from_mean[constants::laben::channels + 1] += 1;

    if (( lg_charge_bad_rms[ch - 1] - prev_lg_charge_bad_rms[ch - 1]) == -1)  lg_charge_bad_rms[constants::laben::channels] += 1;
    if (( lg_charge_bad_rms[ch - 1] - prev_lg_charge_bad_rms[ch - 1]) ==  1)  lg_charge_bad_rms[constants::laben::channels + 1] += 1;

    if (( lg_charge_low_gain[ch - 1] - prev_lg_charge_low_gain[ch - 1]) == -1)  lg_charge_low_gain[constants::laben::channels] += 1;
    if (( lg_charge_low_gain[ch - 1] - prev_lg_charge_low_gain[ch - 1]) ==  1)  lg_charge_low_gain[constants::laben::channels + 1] += 1;

    if (( lg_charge_high_gain[ch - 1] - prev_lg_charge_high_gain[ch - 1]) == -1)  lg_charge_high_gain[constants::laben::channels] += 1;
    if (( lg_charge_high_gain[ch - 1] - prev_lg_charge_high_gain[ch - 1]) ==  1)  lg_charge_high_gain[constants::laben::channels + 1] += 1;
   }

  get_message(bx_message::log) << "lg_fifo_empty: " <<  lg_fifo_empty[constants::laben::channels] << " lg lost and " << lg_fifo_empty[constants::laben::channels + 1] << " lg gained it" << dispatch;
  get_message(bx_message::log) << "lg_fifo_full: " <<  lg_fifo_full[constants::laben::channels] << " lg lost and " << lg_fifo_full[constants::laben::channels + 1] << " lg gained it" << dispatch;

  get_message(bx_message::log) << " lg_loosing_raw_hits_p: " <<   lg_loosing_raw_hits_p[constants::laben::channels] << " lg lost and " <<  lg_loosing_raw_hits_p[constants::laben::channels + 1] << " lg gained it" << dispatch;

  get_message(bx_message::log) << " lg_loosing_raw_hits_l: " <<   lg_loosing_raw_hits_l[constants::laben::channels] << " lg lost and " <<  lg_loosing_raw_hits_l[constants::laben::channels + 1] << " lg gained it" << dispatch;
  
  get_message(bx_message::log) << " lg_loosing_raw_hits_n: " <<   lg_loosing_raw_hits_n[constants::laben::channels] << " lg lost and " <<  lg_loosing_raw_hits_n[constants::laben::channels + 1] << " lg gained it" << dispatch;

  get_message(bx_message::log) << " lg_laser_ref_correl: " <<   lg_laser_ref_correl[constants::laben::channels] << " lg lost and " <<  lg_laser_ref_correl[constants::laben::channels + 1] << " lg gained it" << dispatch;

  get_message(bx_message::log) << " lg_trigger_ref_correl_p: " <<   lg_trigger_ref_correl_p[constants::laben::channels] << " lg lost and " <<  lg_trigger_ref_correl_p[constants::laben::channels + 1] << " lg gained it" << dispatch;

  get_message(bx_message::log) << " lg_trigger_ref_correl_l: " <<   lg_trigger_ref_correl_l[constants::laben::channels] << " lg lost and " <<  lg_trigger_ref_correl_l[constants::laben::channels + 1] << " lg gained it" << dispatch;

  get_message(bx_message::log) << " lg_trigger_ref_correl_n: " <<   lg_trigger_ref_correl_n[constants::laben::channels] << " lg lost and " <<  lg_trigger_ref_correl_n[constants::laben::channels + 1] << " lg gained it" << dispatch;

  get_message(bx_message::log) << " lg_end_of_gate_correl_p: " <<   lg_end_of_gate_correl_p[constants::laben::channels] << " lg lost and " <<  lg_end_of_gate_correl_p[constants::laben::channels + 1] << " lg gained it" << dispatch;
  
  get_message(bx_message::log) << " lg_end_of_gate_correl_l: " <<   lg_end_of_gate_correl_l[constants::laben::channels] << " lg lost and " <<  lg_end_of_gate_correl_l[constants::laben::channels + 1] << " lg gained it" << dispatch;
 
  get_message(bx_message::log) << " lg_end_of_gate_correl_n: " <<   lg_end_of_gate_correl_n[constants::laben::channels] << " lg lost and " <<  lg_end_of_gate_correl_n[constants::laben::channels + 1] << " lg gained it" << dispatch;

  get_message(bx_message::log) << " lg_bad_timing_shape_in_laser: " <<   lg_bad_timing_shape_in_laser[constants::laben::channels] << " lg lost and " <<  lg_bad_timing_shape_in_laser[constants::laben::channels + 1] << " lg gained it" << dispatch;

  get_message(bx_message::log) << " lg_dead_p: " <<   lg_dead_p[constants::laben::channels] << " lg lost and " <<  lg_dead_p[constants::laben::channels + 1] << " lg gained it" << dispatch;

  get_message(bx_message::log) << " lg_dead_l: " <<   lg_dead_l[constants::laben::channels] << " lg lost and " <<  lg_dead_l[constants::laben::channels + 1] << " lg gained it" << dispatch;

  get_message(bx_message::log) << " lg_dead_n: " <<   lg_dead_n[constants::laben::channels] << " lg lost and " <<  lg_dead_n[constants::laben::channels + 1] << " lg gained it" << dispatch;

  get_message(bx_message::log) << " lg_dead_raw_n: " <<   lg_dead_raw[constants::laben::channels] << " lg lost and " <<  lg_dead_raw[constants::laben::channels + 1] << " lg gained it" << dispatch;

  get_message(bx_message::log) << " lg_low_eff_p: " <<   lg_low_eff_p[constants::laben::channels] << " lg lost and " <<  lg_low_eff_p[constants::laben::channels + 1] << " lg gained it" << dispatch;

  get_message(bx_message::log) << " lg_low_eff_l: " <<   lg_low_eff_l[constants::laben::channels] << " lg lost and " <<  lg_low_eff_l[constants::laben::channels + 1] << " lg gained it" << dispatch;

  get_message(bx_message::log) << " lg_low_eff_n: " <<   lg_low_eff_n[constants::laben::channels] << " lg lost and " <<  lg_low_eff_n[constants::laben::channels + 1] << " lg gained it" << dispatch;

  get_message(bx_message::log) << " lg_low_eff_raw: " <<   lg_low_eff_raw[constants::laben::channels] << " lg lost and " <<  lg_low_eff_raw[constants::laben::channels + 1] << " lg gained it" << dispatch;

  get_message(bx_message::log) << " lg_hot_p: " <<   lg_hot_p[constants::laben::channels] << " lg lost and " <<  lg_hot_p[constants::laben::channels + 1] << " lg gained it" << dispatch;

  get_message(bx_message::log) << " lg_hot_l: " <<   lg_hot_l[constants::laben::channels] << " lg lost and " <<  lg_hot_l[constants::laben::channels + 1] << " lg gained it" << dispatch;

  get_message(bx_message::log) << " lg_hot_n: " <<   lg_hot_n[constants::laben::channels] << " lg lost and " <<  lg_hot_n[constants::laben::channels + 1] << " lg gained it" << dispatch;

  get_message(bx_message::log) << " lg_hot_raw_n: " <<   lg_hot_raw[constants::laben::channels] << " lg lost and " <<  lg_hot_raw[constants::laben::channels + 1] << " lg gained it" << dispatch;

  get_message(bx_message::log) << " lg_retriggering_p: " <<   lg_retriggering_p[constants::laben::channels] << " lg lost and " <<  lg_retriggering_p[constants::laben::channels + 1] << " lg gained it" << dispatch;

  get_message(bx_message::log) << " lg_retriggering_l: " <<   lg_retriggering_l[constants::laben::channels] << " lg lost and " <<  lg_retriggering_l[constants::laben::channels + 1] << " lg gained it" << dispatch;

  get_message(bx_message::log) << " lg_retriggering_n: " <<   lg_retriggering_n[constants::laben::channels] << " lg lost and " <<  lg_retriggering_n[constants::laben::channels + 1] << " lg gained it" << dispatch;

  get_message(bx_message::log) << " lg_base_many_zeros: " <<   lg_base_many_zeros[constants::laben::channels] << " lg lost and " <<  lg_base_many_zeros[constants::laben::channels + 1] << " lg gained it" << dispatch;

  get_message(bx_message::log) << " lg_base_many_FF: " <<   lg_base_many_FF[constants::laben::channels] << " lg lost and " <<  lg_base_many_FF[constants::laben::channels + 1] << " lg gained it" << dispatch;

  get_message(bx_message::log) << " lg_base_too_spread: " <<   lg_base_too_spread[constants::laben::channels ] << " lg lost and " <<  lg_base_too_spread[constants::laben::channels + 1] << " lg gained it" << dispatch;

  get_message(bx_message::log) << " lg_base_shifted_from_mean: " <<   lg_base_shifted_from_mean[constants::laben::channels] << " lg lost and " <<  lg_base_shifted_from_mean[constants::laben::channels + 1] << " lg gained it" << dispatch;

  get_message(bx_message::log) << " lg_base_bad_rms: " <<   lg_base_bad_rms[constants::laben::channels] << " lg lost and " <<  lg_base_bad_rms[constants::laben::channels + 1] << " lg gained it" << dispatch;

  get_message(bx_message::log) << " lg_base_very_small_rms: " <<   lg_base_very_small_rms[constants::laben::channels] << " lg lost and " <<  lg_base_very_small_rms[constants::laben::channels + 1] << " lg gained it" << dispatch;

  get_message(bx_message::log) << " lg_peak_many_zeros: " <<   lg_peak_many_zeros[constants::laben::channels] << " lg lost and " <<  lg_peak_many_zeros[constants::laben::channels + 1] << " lg gained it" << dispatch;

  get_message(bx_message::log) << " lg_peak_many_FF: " <<   lg_peak_many_FF[constants::laben::channels] << " lg lost and " <<  lg_peak_many_FF[constants::laben::channels + 1] << " lg gained it" << dispatch;

  get_message(bx_message::log) << " lg_peak_too_spread: " <<   lg_peak_too_spread[constants::laben::channels] << " lg lost and " <<  lg_peak_too_spread[constants::laben::channels + 1] << " lg gained it" << dispatch;

  get_message(bx_message::log) << " lg_peak_shifted_from_mean: " <<   lg_peak_shifted_from_mean[constants::laben::channels] << " lg lost and " <<  lg_peak_shifted_from_mean[constants::laben::channels + 1] << " lg gained it" << dispatch;

  get_message(bx_message::log) << " lg_peak_bad_rms: " <<   lg_peak_bad_rms[constants::laben::channels] << " lg lost and " <<  lg_peak_bad_rms[constants::laben::channels + 1] << " lg gained it" << dispatch;

  get_message(bx_message::log) << " lg_peak_very_small_rms: " <<   lg_peak_very_small_rms[constants::laben::channels] << " lg lost and " <<  lg_peak_very_small_rms[constants::laben::channels + 1] << " lg gained it" << dispatch;  

  get_message(bx_message::log) << " lg_charge_many_zeros: " <<   lg_charge_many_zeros[constants::laben::channels] << " lg lost and " <<  lg_charge_many_zeros[constants::laben::channels + 1] << " lg gained it" << dispatch;
 
  get_message(bx_message::log) << " lg_charge_many_FF: " <<   lg_charge_many_FF[constants::laben::channels] << " lg lost and " <<  lg_charge_many_FF[constants::laben::channels + 1] << " lg gained it" << dispatch;

  get_message(bx_message::log) << " lg_charge_many_negative_values: " << lg_charge_many_negative_values[constants::laben::channels] << " lg lost and " << lg_charge_many_negative_values[constants::laben::channels + 1] << " lg gained it" << dispatch;
  
  get_message(bx_message::log) << " lg_charge_too_spread: " <<   lg_charge_too_spread[constants::laben::channels] << " lg lost and " <<  lg_charge_too_spread[constants::laben::channels + 1] << " lg gained it" << dispatch;
  
  get_message(bx_message::log) << " lg_charge_shifted_from_mean: " <<   lg_charge_shifted_from_mean[constants::laben::channels] << " lg lost and " <<  lg_charge_shifted_from_mean[constants::laben::channels + 1] << " lg gained it" << dispatch;
  
  get_message(bx_message::log) << " lg_charge_bad_rms: " <<   lg_charge_bad_rms[constants::laben::channels] << " lg lost and " <<  lg_charge_bad_rms[constants::laben::channels + 1] << " lg gained it" << dispatch;
   
  get_message(bx_message::log) << " lg_charge_low_gain: " <<   lg_charge_low_gain[constants::laben::channels] << " lg lost and " <<  lg_charge_low_gain[constants::laben::channels + 1] << " lg gained it" << dispatch;
   
  get_message(bx_message::log) << " lg_charge_high_gain: " <<   lg_charge_high_gain[constants::laben::channels] << " lg lost and " <<  lg_charge_high_gain[constants::laben::channels + 1] << " lg gained it" << dispatch;
   
      // Delete vectors
  delete [] hits[0];
  delete [] hits[1];
  delete [] hits[2];
    
    //check if the situation did not change with respect to the previous run
  int32_t i_Nlg_status_change   = get_parameter ("Nlg_status_change").get_int ();

     //for base
  if( abs(Nlg_base_many_zeros - prev_Nlg_base_many_zeros) >  i_Nlg_status_change)
    get_message(bx_message::log) << "Nlg_base_many_zeros changed by " << (Nlg_base_many_zeros - prev_Nlg_base_many_zeros) << " from the previous run " <<  dispatch;
  if( abs(Nlg_base_many_FF - prev_Nlg_base_many_FF) >  i_Nlg_status_change)
    get_message(bx_message::log) << "Nlg_base_many_FF changed by " << (Nlg_base_many_FF - prev_Nlg_base_many_FF) << " from the previous run " << dispatch;
  if( abs(Nlg_base_too_spread - prev_Nlg_base_too_spread) >  i_Nlg_status_change)
    get_message(bx_message::log) << "Nlg_base_too_spread changed by " << (Nlg_base_too_spread - prev_Nlg_base_too_spread) << " from the previous run " << dispatch;
  if( abs(Nlg_base_shifted_from_mean - prev_Nlg_base_shifted_from_mean) >  i_Nlg_status_change)
    get_message(bx_message::log) << "Nlg_base_shifted_from_mean changed by " << (Nlg_base_shifted_from_mean - prev_Nlg_base_shifted_from_mean) << " from the previous run " << dispatch;
  if( abs(Nlg_base_bad_rms - prev_Nlg_base_bad_rms) >  i_Nlg_status_change)
    get_message(bx_message::log) << "Nlg_base_bad_rms changed by " << (Nlg_base_bad_rms - prev_Nlg_base_bad_rms) << " from the previous run "  << dispatch;
  if( abs(Nlg_base_very_small_rms - prev_Nlg_base_very_small_rms) >  i_Nlg_status_change)
    get_message(bx_message::log) << "Nlg_base_very_small_rms changed by " << (Nlg_base_very_small_rms - prev_Nlg_base_very_small_rms) << " from the previous run "  << dispatch;

    //for peak
  if( abs(Nlg_peak_many_zeros - prev_Nlg_peak_many_zeros) >  i_Nlg_status_change)
    get_message(bx_message::log) << "Nlg_peak_many_zeros changed by " << (Nlg_peak_many_zeros - prev_Nlg_peak_many_zeros) << " from the previous run " <<  dispatch;
  if( abs(Nlg_peak_many_FF - prev_Nlg_peak_many_FF) >  i_Nlg_status_change)
    get_message(bx_message::log) << "Nlg_peak_many_FF changed by " << (Nlg_peak_many_FF - prev_Nlg_peak_many_FF) << " from the previous run " << dispatch;
  if( abs(Nlg_peak_too_spread - prev_Nlg_peak_too_spread) >  i_Nlg_status_change)
    get_message(bx_message::log) << "Nlg_peak_too_spread changed by " << (Nlg_peak_too_spread - prev_Nlg_peak_too_spread) << " from the previous run " << dispatch;
  if( abs(Nlg_peak_shifted_from_mean - prev_Nlg_peak_shifted_from_mean) >  i_Nlg_status_change)
    get_message(bx_message::log) << "Nlg_peak_shifted_from_mean changed by " << (Nlg_peak_shifted_from_mean - prev_Nlg_peak_shifted_from_mean) << " from the previous run " << dispatch;
  if( abs(Nlg_peak_bad_rms - prev_Nlg_peak_bad_rms) >  i_Nlg_status_change)
    get_message(bx_message::log) << "Nlg_peak_bad_rms changed by " << (Nlg_peak_bad_rms - prev_Nlg_peak_bad_rms) << " from the previous run "  << dispatch;
  if( abs(Nlg_peak_very_small_rms - prev_Nlg_peak_very_small_rms) >  i_Nlg_status_change)
    get_message(bx_message::log) << "Nlg_peak_very_small_rms changed by " << (Nlg_peak_very_small_rms - prev_Nlg_peak_very_small_rms) << " from the previous run "  << dispatch;

   //for charge
  if( abs(Nlg_charge_many_zeros - prev_Nlg_charge_many_zeros) >  i_Nlg_status_change)
    get_message(bx_message::log) << "Nlg_charge_many_zeros changed by " << (Nlg_charge_many_zeros - prev_Nlg_charge_many_zeros) << " from the previous run " <<  dispatch;
  if( abs(Nlg_charge_many_FF - prev_Nlg_charge_many_FF) >  i_Nlg_status_change)
    get_message(bx_message::log) << "Nlg_charge_many_FF changed by " << (Nlg_charge_many_FF - prev_Nlg_charge_many_FF) << " from the previous run " << dispatch;
  if( abs(Nlg_charge_many_negative_values - prev_Nlg_charge_many_negative_values) >  i_Nlg_status_change)
    get_message(bx_message::log) << "Nlg_charge_many_negative_values changed by " << (Nlg_charge_many_negative_values - prev_Nlg_charge_many_negative_values) << " from the previous run " << dispatch;
  if( abs(Nlg_charge_too_spread - prev_Nlg_charge_too_spread) >  i_Nlg_status_change)
    get_message(bx_message::log) << "Nlg_charge_too_spread changed by " << (Nlg_charge_too_spread - prev_Nlg_charge_too_spread) << " from the previous run " << dispatch;
  if( abs(Nlg_charge_shifted_from_mean - prev_Nlg_charge_shifted_from_mean) >  i_Nlg_status_change)
    get_message(bx_message::log) << "Nlg_charge_shifted_from_mean changed by " << (Nlg_charge_shifted_from_mean - prev_Nlg_charge_shifted_from_mean) << " from the previous run " << dispatch;
  if( abs(Nlg_charge_bad_rms - prev_Nlg_charge_bad_rms) >  i_Nlg_status_change)
    get_message(bx_message::log) << "Nlg_charge_bad_rms changed by " << (Nlg_charge_bad_rms - prev_Nlg_charge_bad_rms) << " from the previous run "  << dispatch;
  if( abs(Nlg_charge_low_gain - prev_Nlg_charge_low_gain) >  i_Nlg_status_change)
    get_message(bx_message::log) << "Nlg_charge_low_gain changed by " << (Nlg_charge_low_gain - prev_Nlg_charge_low_gain) << " from the previous run "  << dispatch;
  if( abs(Nlg_charge_high_gain - prev_Nlg_charge_high_gain) >  i_Nlg_status_change)
    get_message(bx_message::log) << "Nlg_charge_high_gain changed by " << (Nlg_charge_high_gain - prev_Nlg_charge_high_gain) << " from the previous run "  << dispatch;


  if( abs(Nlg_laser_ref_correl - prev_Nlg_laser_ref_correl) >  i_Nlg_status_change)
    get_message(bx_message::log) << "Nlg_laser_ref_correl changed by " << (Nlg_laser_ref_correl - prev_Nlg_laser_ref_correl) << " from the previous run "  << dispatch;
    
  for(int32_t trg = 0; trg < 3; trg ++){
    if( abs(Nlg_dead[trg] - prev_Nlg_dead[trg]) >  i_Nlg_status_change)
      get_message(bx_message::log) << "Nlg_dead_in_" << trg_names[trg] << " changed by " << (Nlg_dead[trg] - prev_Nlg_dead[trg]) << " from the previous run "  << dispatch;
  }
  for(int32_t trg = 0; trg < 3; trg ++){
    if( abs(Nlg_low_eff[trg] - prev_Nlg_low_eff[trg]) >  i_Nlg_status_change)
      get_message(bx_message::log) << "Nlg_low_eff_in_" << trg_names[trg] << " changed by " << (Nlg_low_eff[trg] - prev_Nlg_low_eff[trg]) << " from the previous run "  << dispatch;
  } 
  for(int32_t trg = 0; trg < 3; trg ++){
    if( abs(Nlg_hot[trg] - prev_Nlg_hot[trg]) >  i_Nlg_status_change)
      get_message(bx_message::log) << "Nlg_hot_in_" << trg_names[trg] << " changed by " << (Nlg_hot[trg] - prev_Nlg_hot[trg]) << " from the previous run "  << dispatch;
  }
  for(int32_t trg = 0; trg < 3; trg ++){
    if( abs(Nlg_retriggering[trg] - prev_Nlg_retriggering[trg]) >  i_Nlg_status_change)
      get_message(bx_message::log) << "Nlg_retriggering_in_" << trg_names[trg] << " changed by " << (Nlg_retriggering[trg] - prev_Nlg_retriggering[trg]) << " from the previous run "  << dispatch;
  }
  for(int32_t trg = 0; trg < 3; trg ++){
    if( abs(Nlg_trigger_ref_correl[trg] - prev_Nlg_trigger_ref_correl[trg]) >  i_Nlg_status_change)
      get_message(bx_message::log) << "Nlg_trigger_ref_correl_in_" << trg_names[trg] << " changed by " << (Nlg_trigger_ref_correl[trg] - prev_Nlg_trigger_ref_correl[trg]) << " from the previous run "  << dispatch;
  } 
  for(int32_t trg = 0; trg < 3; trg ++){
    if( abs(Nlg_end_of_gate_correl[trg] - prev_Nlg_end_of_gate_correl[trg]) >  i_Nlg_status_change)
      get_message(bx_message::log) << "Nlg_end_of_gate_correl_in_" << trg_names[trg] << " changed by " << (Nlg_end_of_gate_correl[trg] - prev_Nlg_end_of_gate_correl[trg]) << " from the previous run "  << dispatch;
  } 
  
  if(abs(prev_Nlg_fifo_empty - Nlg_fifo_empty) > i_Nlg_status_change)
    get_message(bx_message::log) << "Nlg_fifo_empty changed by " << prev_Nlg_fifo_empty - Nlg_fifo_empty << " channels "<< dispatch; 						      
  if(abs(prev_Nlg_fifo_full - Nlg_fifo_full) > i_Nlg_status_change)
    get_message(bx_message::log) << "Nlg_fifo_full changed by " << prev_Nlg_fifo_full - Nlg_fifo_full << " channels " << dispatch; 						      
 
  for(int32_t trg = 0; trg < 3; trg ++){
    if(abs(prev_Nlg_loosing_raw_hits[trg] - Nlg_loosing_raw_hits[trg]) > i_Nlg_status_change)
      get_message(bx_message::log) << "Nlgl_loosing_raw_hits_in_" << trg_names[trg] << " changed by " << Nlg_loosing_raw_hits[trg] - prev_Nlg_loosing_raw_hits[trg] << " channels " << dispatch; 
  }					      
	
  if(abs(Nlg_bad_timing_shape_in_laser - prev_Nlg_bad_timing_shape_in_laser) >  i_Nlg_status_change)
    get_message(bx_message::log) << "Nlg_bad_timing_shape_in_laser changed by " << (Nlg_bad_timing_shape_in_laser - prev_Nlg_bad_timing_shape_in_laser) << " from the previous run "  << dispatch;
  

    //set visitors and write to the DB (mean neutrino > 400) 
  if(get_parameter ("db_write").get_bool () &&  Nlg_dead[2] < get_parameter ("Nlg_dead_nu_thresh").get_float () &&  mean_nhits[2] > get_parameter ("mean_nhits_nu_thresh").get_float ()  && module_says_DB_write){ 
    db_run& run_info = bx_dbi::get ()->get_run ();
    std::vector<std::string> vstr;
    for(int32_t ch = 0; ch < constants::laben::channels; ch++) {
      
      vstr.clear ();
      for (unsigned i = 0; i < charge_base_status_vec[ch].size (); i++) vstr.push_back (ADC_translation_map.rfind (charge_base_status_vec[ch][i])->first);
      run_info.set_laben_charge_base_status  (ch + 1, vstr, this);
     
      vstr.clear ();
      for (unsigned i = 0; i < charge_peak_status_vec[ch].size (); i++) vstr.push_back (ADC_translation_map.rfind ( charge_peak_status_vec[ch][i])->first);
      run_info.set_laben_charge_peak_status (ch + 1, vstr, this);
     
      vstr.clear ();     
      for (unsigned i = 0; i < charge_status_vec[ch].size (); i++) vstr.push_back (ADC_translation_map.rfind (charge_status_vec[ch][i])->first);
      run_info.set_laben_charge_status  (ch + 1, vstr, this); 
      
      vstr.clear ();
      for (unsigned i = 0; i < timing_status_pl_vec[ch].size (); i++) vstr.push_back (timing_translation_map.rfind (timing_status_pl_vec[ch][i])->first);
        //enough neutrino triggers using current run
      if(mean_nhits[2] > 100) 
	for (unsigned i = 0; i < timing_status_n_vec[ch].size (); i++) vstr.push_back (timing_translation_map.rfind (timing_status_n_vec[ch][i])->first);
        //not enough neutrino triggers using previous run
      if(mean_nhits[2] < 100) 
	for (unsigned i = 0; i < prev_timing_status_n_vec[ch].size (); i++) vstr.push_back (timing_translation_map.rfind (prev_timing_status_n_vec[ch][i])->first);      
      run_info.set_laben_timing_status  (ch + 1, vstr, this);
     
      vstr.clear ();
      for (unsigned i = 0; i < multiplicity_pl_vec[ch].size (); i++) vstr.push_back (multiplicity_translation_map.rfind (multiplicity_pl_vec[ch][i])->first);
      //enough neutrino triggers using current run
      if(mean_nhits[2] > 100) 
	for (unsigned i = 0; i < multiplicity_n_vec[ch].size (); i++) vstr.push_back (multiplicity_translation_map.rfind (multiplicity_n_vec[ch][i])->first);
      //not enough neutrino triggers using previous run
      if(mean_nhits[2] < 100) 
	for (unsigned i = 0; i < prev_multiplicity_n_vec[ch].size (); i++) vstr.push_back (multiplicity_translation_map.rfind (prev_multiplicity_n_vec[ch][i])->first);
      run_info.set_laben_multiplicity  (ch + 1, vstr, this);
    }

     //write to the DB
    run_info.write_laben_electronic_channel (true, this);
    get_message(bx_message::log) << "Writing to DB" << dispatch; 
  }
}




  /*   
  double laser_peaks_ratio = 0;
  double error_laser_peaks_ratio = 0;
    //time shift of reflected peak in ns
    //FIXME READ n from database
  double refr_indx = 1;
  double time_shift =  13.4 / 3E+8 * 1E+9 * refr_indx;
  int32_t bin_of_maximum;
  double rms_main_laser_peak;


    //mean values for laser-hit time distribution
    //ratio of Nhits in direct and reflected peak for all Lg  
  TH1D* lasref_correl_all_lg = lasref_correl_vs_lg->ProjectionY("lasref_correl_one_lg",1,2240);
  bin_of_maximum = lasref_correl_all_lg->GetMaximumBin ();
  if(bin_of_maximum > 0){
    double time_of_maximum = bin_of_maximum - 500.; 
    lasref_correl_all_lg->SetAxisRange(time_of_maximum - 50, time_of_maximum + 50);
    rms_main_laser_peak = lasref_correl_all_lg->GetRMS ();
    std::cout << "main laser peak " <<  time_of_maximum << " ns and rms of " << rms_main_laser_peak << "time shift is  " << time_shift << std::endl;
    if(rms_main_laser_peak > 20) 
      get_message(bx_message::warn) << "Too braod laser--peak resolution of rms " << rms_main_laser_peak << " ns " << dispatch;  
    lasref_correl_all_lg->SetAxisRange(-500, 1500);
      //check if we are within histo bounderies and calculate ratio of hits in the direct and reflected peak
    if((bin_of_maximum - 2 * (int32_t) rms_main_laser_peak) > 0 && (bin_of_maximum + (int32_t) time_shift + 2 * (int32_t) rms_main_laser_peak) < 1500) {  
      double Nhits_main_peak_all_lg = lasref_correl_all_lg->Integral(bin_of_maximum - 2 *(int32_t) rms_main_laser_peak, bin_of_maximum + 2 * (int32_t) rms_main_laser_peak );
      double Nhits_reflected_peak_all_lg = lasref_correl_all_lg->Integral(bin_of_maximum + (int32_t) time_shift - 3 * (int32_t) rms_main_laser_peak, bin_of_maximum + (int32_t) time_shift + 3 *(int32_t) rms_main_laser_peak);
      if(Nhits_main_peak_all_lg && Nhits_reflected_peak_all_lg){
	laser_peaks_ratio = Nhits_reflected_peak_all_lg/ Nhits_main_peak_all_lg;
	error_laser_peaks_ratio = std::sqrt( pow((laser_peaks_ratio / sqrt(Nhits_reflected_peak_all_lg)), 2) + pow ((laser_peaks_ratio/sqrt(Nhits_main_peak_all_lg)),2)); 
	get_message(bx_message::log) << "Laser_peaks_ratio " << laser_peaks_ratio << " +- " << error_laser_peaks_ratio << ", maximum bin is " << bin_of_maximum << " corresponding to " << bin_of_maximum - 500 << " ns after laser trigger" << dispatch;
      }
      else get_message(bx_message::warn) << "Zero hits in the main peaks region  for laser triggers" << dispatch;
    }    
    else get_message(bx_message::warn) << "Laser peak out oh histogram bounders" << dispatch;
  }
  else get_message(bx_message:: warn) << "Laser hits time distribution has maximum before laser_time!" << dispatch;
  */


	    /*	      //strange hit-time distribution in laser triggers ?
	    lasref_correl_one_lg->SetAxisRange(-500,1500);
	    double Nhits_main_peak_one_lg = lasref_correl_one_lg->Integral(bin_of_maximum - 2 *(int32_t) rms_main_laser_peak , bin_of_maximum + 2 *(int32_t) rms_main_laser_peak);
	    double Nhits_reflected_peak_one_lg = lasref_correl_one_lg->Integral(bin_of_maximum + (int32_t) time_shift - 3 *(int32_t) rms_main_laser_peak, bin_of_maximum + (int32_t) time_shift + 3 *(int32_t) rms_main_laser_peak);
	    if( std::fabs(lasref_correl_one_lg->GetMaximumBin () - bin_of_maximum) < 90){ 
	      std::cout << "Lg " <<  ch << " Nhits main peak " << Nhits_main_peak_one_lg << " reflected " << Nhits_reflected_peak_one_lg << std::endl;
	      if(Nhits_main_peak_one_lg && Nhits_reflected_peak_one_lg){
		double laser_peaks_ratio_one_lg = Nhits_reflected_peak_one_lg/ Nhits_main_peak_one_lg;
		double error_laser_peaks_ratio_one_lg = std::sqrt( pow((laser_peaks_ratio_one_lg / sqrt(Nhits_reflected_peak_one_lg)), 2) + pow ((laser_peaks_ratio_one_lg/sqrt(Nhits_main_peak_one_lg)),2)); 
				
		if( std::fabs(laser_peaks_ratio_one_lg - laser_peaks_ratio) > 5 * (error_laser_peaks_ratio_one_lg + error_laser_peaks_ratio)){
		Nlg_bad_timing_shape_in_laser ++;
		get_message(bx_message::log) << "Lg " << ch << " bad time distribution  of laser hits" << dispatch;
		timing_status_vec[ch - 1].push_back(bad_timing_shape_in_laser);
		}
	      }
	      else {
		get_message(bx_message::log) << "Lg " << ch << " Zero hits in the direct or reflected peak region for laser " << dispatch;
		Nlg_bad_timing_shape_in_laser ++;
		timing_status_vec[ch - 1].push_back(bad_timing_shape_in_laser);
	      }
	    }
	    else {
	      Nlg_bad_timing_shape_in_laser ++;
	      timing_status_vec[ch - 1].push_back( bad_timing_shape_in_laser);
	      get_message(bx_message::log) << " Lg " << ch << " Bin_of_max has shifted position in laser hits to bin " << lasref_correl_one_lg->GetMaximumBin () << ", lg  has " << lasref_correl_one_lg->Integral (0,2001) <<  " total hits " << dispatch;
	    }
	    */
   

     
    
