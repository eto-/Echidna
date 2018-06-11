/* BOREXINO Reconstruction program
 *
 * Author: Livia Ludhova<livia.ludhova@mi.infn.it>
 * Maintainer: Livia Ludhova<livia.ludhova@mi.infn.it>
 * 
 *
 */
#include <algorithm>
#include "bx_detector_monitor.hh"
#include "bx_echidna_event.hh"
#include "db_run.hh"
#include "db_channel.hh"
#include "db_calib.hh"
#include "bx_dbi.hh"
#include "constants.hh"
#include "barn_interface.hh"
#include "db_profile.hh"
#include "constants.hh"
#include "laben_event_shape.hh"

#include <assert.h>

#include <cmath> 

// ctor
bx_detector_monitor::bx_detector_monitor() : bx_base_module("bx_detector_monitor", bx_base_module::main_loop) {
	require_event_stage (bx_detector::laben, bx_base_event::decoded);
 
  module_type_map[crate]   = "crate";
  module_type_map[hvb]     = "hvb";
  module_type_map[feb]     = "feb";
  module_type_map[lbnb]    = "lbnb";
  module_type_map[channel] = "channel";
 
  trigger_type_map[pulser]   = "pulser";
  trigger_type_map[random]   = "random";
  trigger_type_map[neutrino] = "neutrino";
  trigger_type_map[laser]    = "laser";

} 

// begin
void bx_detector_monitor::begin () {

  calibration_mode  = get_parameter ("calibration_mode").get_int ();
  disabling_lg  = get_parameter ("disabling_lg").get_int ();
  if (get_parameter_other ("bx_reco_framework", "skip_events").get_int () > 0) disabling_lg = false;
  if (get_parameter_other ("bx_filter_trigger", "skip_trigger_types").get_vector ().size () > 0) disabling_lg = false;

  check_after_precalib_done = 0;
 
  n_trgrate_cycle = 0;
  n_events_in_this_cycle = 0;
  n_nhits_in_cycle = 0;

  nhits_pulser_sum = 0;
  
  n[crate] = constants::laben::ncrates; //14
  n[hvb] = (int32_t) n[crate] *  constants::laben::frontend::board_per_rack/2;
  n[feb] = n[crate] *  constants::laben::frontend::board_per_rack; // 14 * 14 = 196
  n[lbnb] = n[crate] * constants::laben::board_per_rack; // 14 * 20 = 280
  n[channel] =   constants::laben::channels; //2240
  start_lg_on = 0;
  start_lg_off = 0;

  if(calibration_mode){
     
      // charge histos
    neutrino_charge_vs_lg = new TH2F("neutrino_charge_vs_lg","neutrino_charge_vs_lg", n[channel], 1., 1. + n[channel], 512, -256., 256.);
    barn_interface::get ()->store (barn_interface::file, neutrino_charge_vs_lg, this);

    charge_pulser_vs_lg = new TH2F("charge_pulser_vs_lg","charge_pulser_vs_lg", n[channel], 1., 1. + n[channel], 512, -256., 256.);
    barn_interface::get ()->store (barn_interface::file, charge_pulser_vs_lg, this);

      // how many time lg changed status to on (off) during the run for pulser and "HV" = neutrino + random + laser triggers
      // 0 = pulser, 1 = random+neutrino+laser 0= off 1= on
    h_lg_Nx_changed[0][0] = new TH1F("h_lg_Nx_changed_to_off_p", "pulser trigger", n[channel], 1., 1. + n[channel]);
    h_lg_Nx_changed[0][0]->SetXTitle("Logical channel number ");
    h_lg_Nx_changed[0][0]->SetYTitle("How many times changed to off");
  
    h_lg_Nx_changed[0][1] = new TH1F("h_lg_Nx_changed_to_on_p", "pulser trigger", n[channel], 1., 1. + n[channel]);
    h_lg_Nx_changed[0][1]->SetXTitle("Logical channel number ");
    h_lg_Nx_changed[0][1]->SetYTitle("How many times changed to on");
    
    h_lg_Nx_changed[1][0] = new TH1F("h_lg_Nx_changed_to_off_HV", "random+neutrino+laser trigger", n[channel], 1., 1. + n[channel]);
    h_lg_Nx_changed[1][0]->SetXTitle("Logical channel number ");
    h_lg_Nx_changed[1][0]->SetYTitle("How many times changed to off");
    
    h_lg_Nx_changed[1][1] = new TH1F("h_lg_Nx_changed_to_on_HV", "random+neutrino+laser trigger", n[channel], 1., 1. + n[channel]);
    h_lg_Nx_changed[1][1]->SetXTitle("Logical channel number ");
    h_lg_Nx_changed[1][1]->SetYTitle("How many times changed to on");
    
    for(int32_t i = 0; i < 2; i++) {
      for(int32_t j = 0; j < 2; j++) barn_interface::get ()->store (barn_interface::file, h_lg_Nx_changed[i][j] , this);
    }
  }

  //for online mode only  
  if(!calibration_mode){
      //neutrino trigger rate (histo ready up to 10 hours)
    neutrino_trigger_rate = new TH1F("neutrino_trigger_rate","neutrino_trigger_rate",600,0,36000);
    neutrino_trigger_rate->SetXTitle("time [seconds, 60 s/bin] ");
    neutrino_trigger_rate->SetYTitle("neutrino triggers [s-1]");
    barn_interface::get ()->store (barn_interface::file, neutrino_trigger_rate, this);
    
      //neutrino nhits rate (histo ready up to 10 hours)
    neutrino_nhits_rate = new TH1F("neutrino_nhits_rate","neutrino_nhits_rate",600,0,36000);
    neutrino_nhits_rate->SetXTitle("time [seconds, 60 s/bin] ");
    neutrino_nhits_rate->SetYTitle("nhits, neutrino trg [s-1]");
    barn_interface::get ()->store (barn_interface::file, neutrino_nhits_rate, this);

      //random nhits rate (histo ready up to 10 hours)
    random_nhits_rate = new TH1F("random_nhits_rate","random_nhits_rate",600/4,0,36000);
    random_nhits_rate->SetXTitle("time [seconds, 240 s/bin] ");
    random_nhits_rate->SetYTitle("nhits, random trg [s-1]");
    barn_interface::get ()->store (barn_interface::file, random_nhits_rate, this);

      //neutrino nhits per event (histo ready up to 10 hours)
    neutrino_nhits_per_event = new TH1F("neutrino_nhits_per_event","neutrino_nhits_per_event",600,0,36000);
    neutrino_nhits_per_event->SetXTitle("time [seconds, 60 s/bin] ");
    neutrino_nhits_per_event->SetYTitle("nhits per event, neutrino trg");
    barn_interface::get ()->store (barn_interface::file, neutrino_nhits_per_event, this);
     
     //random nhits per event (histo ready up to 10 hours)
    random_nhits_per_event = new TH1F("random_nhits_per_event","random_nhits_per_event",600/4,0,36000);
    random_nhits_per_event->SetXTitle("time [seconds, 240 s/bin] ");
    random_nhits_per_event->SetYTitle("nhits per event, random trg ");
    barn_interface::get ()->store (barn_interface::file, random_nhits_per_event, this);

    //pulser nhits rate (histo ready up to 10 hours)
    pulser_bins = 800000;
    max_events  = 800000;
    pulser_nhits = new TH1F("pulser_nhits","pulser_nhits",pulser_bins,1,max_events);
    pulser_nhits->SetXTitle("event number ");
    pulser_nhits->SetYTitle("nhits, pulser trg ");
    pulser_nhits->SetMarkerStyle(20);
    //pulser_nhits->SetMarkerSize(1);
  

    //occupancy
    //crates
    dec_occupancy[crate][0] = new TH1F ("dec_occupancy_cr_p", "crate in pulser" , n[crate], 1., 1. + n[crate]);
    dec_occupancy[crate][0]->SetXTitle("crate number");
    
    dec_occupancy[crate][1] = new TH1F ("dec_occupancy_cr_r", "crate in random" , n[crate], 1., 1. + n[crate]);
    dec_occupancy[crate][1]->SetXTitle("crate number");
    
    dec_occupancy[crate][2] = new TH1F ("dec_occupancy_cr_n", "crate in neutrino" , n[crate], 1., 1. + n[crate]);
    dec_occupancy[crate][2]->SetXTitle("crate number");
    
    dec_occupancy[crate][3] = new TH1F ("dec_occupancy_cr_l", "crate in laser" , n[crate], 1., 1. + n[crate]);
    dec_occupancy[crate][3]->SetXTitle("crate number");  
  
      //hvb
    dec_occupancy[hvb][0] = new TH1F ("dec_occupancy_hvb_p", "hvb in pulser" , n[hvb], 1., 1. + n[hvb]);
    dec_occupancy[hvb][0]->SetXTitle("hvb number");
    
    dec_occupancy[hvb][1] = new TH1F ("dec_occupancy_hvb_r", "hvb in random" , n[hvb], 1., 1. + n[hvb]);
    dec_occupancy[hvb][1]->SetXTitle("hvb number");
    
    dec_occupancy[hvb][2] = new TH1F ("dec_occupancy_hvb_n", "hvb in neutrino" , n[hvb], 1., 1. + n[hvb]);
    dec_occupancy[hvb][2]->SetXTitle("hvb number");
    
    dec_occupancy[hvb][3] = new TH1F ("dec_occupancy_hvb_l", "hvb in laser" , n[hvb], 1., 1. + n[hvb]);
    dec_occupancy[hvb][3]->SetXTitle("hvb number");
    
     //feb
    dec_occupancy[feb][0] = new TH1F ("dec_occupancy_feb_p", "feb in pulser" , n[feb], 1., 1. + n[feb]);
    dec_occupancy[feb][0]->SetXTitle("feb number");
    
    dec_occupancy[feb][1] = new TH1F ("dec_occupancy_feb_r", "feb in random" , n[feb], 1., 1. + n[feb]);
    dec_occupancy[feb][1]->SetXTitle("feb number");
    
    dec_occupancy[feb][2] = new TH1F ("dec_occupancy_feb_n", "feb in neutrino" , n[feb], 1., 1. + n[feb]);
    dec_occupancy[feb][2]->SetXTitle("feb number");
    
    dec_occupancy[feb][3] = new TH1F ("dec_occupancy_feb_l", "feb in laser" , n[feb], 1., 1. + n[feb]);
    dec_occupancy[feb][3]->SetXTitle("feb number");

     //lbnb 
    dec_occupancy[lbnb][0] = new TH1F ("dec_occupancy_lbnb_p", "lbnb in pulser" , n[lbnb], 1., 1. + n[lbnb]);
    dec_occupancy[lbnb][0]->SetXTitle("lbnb number");
    
    dec_occupancy[lbnb][1] = new TH1F ("dec_occupancy_lbnb_r", "lbnb in random" , n[lbnb], 1., 1. + n[lbnb]);
    dec_occupancy[lbnb][1]->SetXTitle("lbnb number");
    
    dec_occupancy[lbnb][2] = new TH1F ("dec_occupancy_lbnb_n", "lbnb in neutrino" , n[lbnb], 1., 1. + n[lbnb]);
    dec_occupancy[lbnb][2]->SetXTitle("lbnb number");
    
    dec_occupancy[lbnb][3] = new TH1F ("dec_occupancy_lbnb_l", "lbnb in laser" , n[lbnb], 1., 1. + n[lbnb]);
    dec_occupancy[lbnb][3]->SetXTitle("lbnb number");

      //channel
    dec_occupancy[channel][0] = new TH1F ("dec_occupancy_ch_p", "channel in pulser" , n[channel], 1., 1. + n[channel]);
    dec_occupancy[channel][0]->SetXTitle("channel number");
    
    dec_occupancy[channel][1] = new TH1F ("dec_occupancy_ch_r", "channel in random" , n[channel], 1., 1. + n[channel]);
    dec_occupancy[channel][1]->SetXTitle("channel number");
    
    dec_occupancy[channel][2] = new TH1F ("dec_occupancy_ch_n", "channel in neutrino" , n[channel], 1., 1. + n[channel]);
    dec_occupancy[channel][2]->SetXTitle("channel number");
    
    dec_occupancy[channel][3] = new TH1F ("dec_occupancy_ch_l", "channel in laser" , n[channel], 1., 1. + n[channel]);
    dec_occupancy[channel][3]->SetXTitle("channel number");
    
    //barn_interface + initialisation
    for(int32_t module_type = 0; module_type < 5; module_type ++) {
      for(int32_t trg = 0; trg < 4; trg ++){ 
	barn_interface::get ()->store (barn_interface::file, dec_occupancy[module_type][trg], this);
	for(int32_t i = 0; i < n[module_type]; i++){
	  times_found_noisy[module_type][trg].push_back (0); 
	  times_found_empty[module_type][trg].push_back(0); 
	}
      }
    }
  }//end of online mode
        
  total_triggers = 0;
  for(int32_t i = 0; i < 4; i++)  n_triggers[i] = 0; 


  if(calibration_mode){
    //  *************************************************
    // vectors for finding OFF channels 
    //  **************************************************

    //nhits and dead lg in 1000 precalibration events
    nhits_in_precalib.resize(n[channel],0);
    
      //list of all disabled lg
    cummulative_off_lg.resize(n[channel],0);
    cummulative_charge_off_lg.resize(n[channel],0);

      
    //how many time the lg changed its status to off (on) in pulser and HV triggers
    lg_Nx_to_off_p.resize(n[channel],0);
    lg_Nx_to_on_p.resize(n[channel],0);
    lg_Nx_to_off_HV.resize(n[channel],0);
    lg_Nx_to_on_HV.resize(n[channel],0);
    
    //the event number of the LATEST change
    ev_lg_to_off_p.resize(n[channel],0);
    ev_lg_to_on_p.resize(n[channel],0);
    ev_lg_to_off_HV.resize(n[channel],0);
    ev_lg_to_on_HV.resize(n[channel],0);
    
    latest_off_lg_p.resize(n[channel],0);
    latest_off_lg_HV.resize(n[channel],0);
    cummulative_occupancy_p.resize(n[channel],0);
    cummulative_lbnb_occupancy_p.resize(n[lbnb],0);
    cummulative_occupancy_HV.resize(n[channel],0);
    cummulative_hvb_occupancy_HV.resize(n[hvb],0);
  }

}

//DO IT

bx_echidna_event* bx_detector_monitor::doit (bx_echidna_event *ev) {

    //disable channels disabling for simulations
  if (ev->is_mctruth_enabled ()) return ev; 
  
  total_triggers ++;

  int32_t evnum = ev->get_event_number ();
  int32_t first_time_hot  = get_parameter ("first_time_hot").get_int ();
  int32_t event_trg_type = -10;

  if (ev->get_trigger ().is_pulser   ()) {
    event_trg_type = 0;
    n_triggers[0] ++;
  } 
  if (ev->get_trigger ().is_random   ()) {
    event_trg_type = 1;
    n_triggers[1] ++; 
  }
  if (ev->get_trigger ().is_neutrino ()) {
    event_trg_type = 2;
    n_triggers[2] ++; 
  }
  if (ev->get_trigger ().is_laser394 ()) {
    event_trg_type = 3;
    n_triggers[3] ++; 
  } 

  if (event_trg_type == -10) return ev; 
      
    //  *************************************************
    // crates, hvb, feb, lbnb and channels occupancy
    //  **************************************************
    
    //fill for decoded hits
     
  std::vector<int32_t> dec_lg_list;
  int32_t dec_nhits_one_event = 0;  //WE CONSIDER ONLY DECODED HITS FROM ORDINARY CHANNELS 
    
  //loop on hits 
  int32_t nhits = ev->get_laben ().get_decoded_nhits ();
  for(int32_t i = 0; i < nhits; i ++) {
    int32_t lg = ev->get_laben ().get_decoded_hit (i).get_raw_hit ().get_logical_channel ();  
      //charge for neutrino events and "not muons"
    if(event_trg_type == 2 && nhits < 1000) {
      double charge = ev->get_laben ().get_decoded_hit (i).get_charge_bin ();
      neutrino_charge_vs_lg->Fill (lg, charge);
    }
      //charge for pulser events
    if(event_trg_type == 0){
      double charge = ev->get_laben ().get_decoded_hit (i).get_charge_bin ();
      charge_pulser_vs_lg->Fill (lg, charge);
    }
    //count ndecoded hits ALWAYS without ref channels, even if in decoded they are enabled
    if (bx_dbi::get ()->get_channel (lg).is_ordinary ())  dec_nhits_one_event ++; 
    
    dec_lg_list.push_back (lg);
  }
  
  //getting occupancies using laben_event_shape  
  laben_event_shape dec_ev_shape (dec_lg_list);

    //laben_event_shape does not consider hits from ref channels, occupancy is divided by the number of nhits in ordinary channels allowed by decoder
  const std::vector<float>& dec_crate_occupancy = dec_ev_shape.get_crate_occupancy ();
  if(!calibration_mode) for (int32_t i = 0; i < n[crate]; i++)  dec_occupancy[crate][event_trg_type]->Fill (i + 1, dec_crate_occupancy[i]);
 
  const std::vector<float>& dec_hvb_occupancy = dec_ev_shape.get_hvb_occupancy ();
  const std::vector<float>& dec_hvb_nhits = dec_ev_shape.get_nhits_in_hvb ();
  if(!calibration_mode) for (int32_t i = 0; i < n[hvb]; i++)  dec_occupancy[hvb][event_trg_type]->Fill (i + 1, dec_hvb_occupancy[i]);
  
  const std::vector<float>& dec_feb_occupancy = dec_ev_shape.get_feb_occupancy ();
  if(!calibration_mode) for (int32_t i = 0; i < n[feb] ; i++)  dec_occupancy[feb][event_trg_type]->Fill (i + 1, dec_feb_occupancy[i]);
  
  const std::vector<float>& dec_lbnb_nhits = dec_ev_shape.get_nhits_in_lbnb ();
  const std::vector<float>& dec_lbnb_occupancy = dec_ev_shape.get_lbnb_occupancy ();
  if(!calibration_mode) for (int32_t i = 0; i < n[lbnb]; i++)  dec_occupancy[lbnb][event_trg_type]->Fill (i + 1, dec_lbnb_occupancy[i]);

  const std::vector<float>& dec_channel_occupancy = dec_ev_shape.get_channel_occupancy ();
  if(!calibration_mode) for (int32_t i = 0; i < n[channel]; i++)  dec_occupancy[channel][event_trg_type]->Fill (i + 1, dec_channel_occupancy[i]);

  const std::vector<float>& dec_nhits_in_channel = dec_ev_shape.get_nhits_in_channel ();

  if(calibration_mode){
  //  *************************************************
  // filling vectors for finding OFF channels 
  //  **************************************************
      
      //channels not having precalibration hits
    if(event_trg_type == 0 && evnum <= 1000){
      for (int32_t i = 0; i < n[channel]; i++) {
	if (bx_dbi::get ()->get_channel (i+1).is_ordinary ()){
	  if (dec_nhits_in_channel[i] != 0) nhits_in_precalib[i] ++;
	}
      }
    }

    if(evnum > 1000 &&  check_after_precalib_done == 0) {

      check_after_precalib_done = 1;

        //A check for dead channels
      for (int32_t i = 0; i < n[channel]; i++) {
	if (bx_dbi::get ()->get_channel (i+1).is_ordinary ()){
	  if (nhits_in_precalib[i] < 10){
	    get_message(bx_message::log) << "Channel " << i+1 << " is off for the whole run, has only " << nhits_in_precalib[i] << " hits in precalib events" << dispatch;
	    cummulative_off_lg[i] = 1;
	    cummulative_charge_off_lg[i] = 1; //no charge for dead channels
	    if(disabling_lg) detector_interface::get ()->add_disabled_channel (i + 1, evnum, db_run::timing, this);
	    start_lg_off ++;
	  }
	  else if (bx_dbi::get ()->get_channel (i+1).is_ordinary ()) start_lg_on ++;
	}
      }
        //check if the off channels belong to a single lbnb
      int32_t N_off_already = 0;
      int32_t N_ordinary_in_lbnb = 0;
      for (int32_t i = 0; i < n[channel]; i++) {
	  //when starting a new board, set to 0 the counter of off events in one lbnb
	if(!(i%8)) {
	  N_off_already = 0;
	  N_ordinary_in_lbnb = 0;
	}
	if(cummulative_off_lg[i]) N_off_already ++;
	if (bx_dbi::get ()->get_channel (i+1).is_ordinary ()) N_ordinary_in_lbnb ++;
	if((i%8) == 7 && i > 0 && N_off_already == N_ordinary_in_lbnb) get_message(bx_message::warn) << "Laben board " << i/constants::laben::channels_per_board + 1<< " (channels " << i + 1 - 7 << "-" <<  i + 1 << ") is off during the whole run" << dispatch; 
      }
      get_message(bx_message::log) << start_lg_off << " ordinary channels OFF after precalib " << dispatch;
      get_message(bx_message::log) << start_lg_on  << " ordinary channels ON after precalib " << dispatch;
      
        //B check for in-pulser-charge-dead channels (among all alive channels)
      for (int32_t i = 0; i < n[channel]; i++) {
	if (bx_dbi::get ()->get_channel (i+1).is_ordinary () && cummulative_charge_off_lg[i] == 0){
	  int32_t ch = i + 1;
	  TH1D* pulser_charge_one_lg  = charge_pulser_vs_lg->ProjectionY ("pulser_charge_one_lg", ch, ch);
	  //0-bin 257
	  double n_ok_charge_hits = pulser_charge_one_lg->Integral (257+5, 257+100);
	  if (n_ok_charge_hits < 10.) {
	    cummulative_charge_off_lg[i] = 1;
	    if(disabling_lg) detector_interface::get ()->add_disabled_channel (i + 1, evnum, db_run::charge, this);
	    get_message(bx_message::log) << "Channel " << i+1 << " after precalib is on, but OFF IN CHARGE (pulser)" << dispatch;
	  }
	  pulser_charge_one_lg->Delete ();
	}
      }
        //check how many off and on ordinary channels in charge
      start_lg_off_charge = 0;
      start_lg_on_charge = 0;
      for (int32_t i = 0; i < n[channel]; i++) {
	if (bx_dbi::get ()->get_channel (i+1).is_ordinary ()){
	  if (cummulative_charge_off_lg[i]) {
	    start_lg_off_charge ++;
	  }
	  else  start_lg_on_charge ++;
	}
      }
      get_message(bx_message::log) << start_lg_off_charge << " ordinary channels OFF IN CHARGE after precalib " << dispatch;
      get_message(bx_message::log) << start_lg_on_charge  << " ordinary channels ON IN CHARGE after precalib " << dispatch;  
   
    }//end of loop after precalib
  
      //check for channels with bad charge in neutrino only
    if( n_triggers[2] == 10000) {
      for (int32_t i = 0; i < n[channel]; i++) {
	if (bx_dbi::get ()->get_channel (i+1).is_ordinary () && cummulative_charge_off_lg[i] == 0){
	  int32_t ch = i + 1;
	  TH1D* neutrino_charge_one_lg  = neutrino_charge_vs_lg->ProjectionY ("neutrino_charge_one_lg", ch, ch);
	  //0-bin 257
	  double n_ok_charge_hits = neutrino_charge_one_lg->Integral (257+5, 257+100);
	  if (n_ok_charge_hits < 10.) {
	    cummulative_charge_off_lg[i] = 1;
	    get_message(bx_message::log) << "Channel " << i+1 << " after 10.000 neutrino triggers, OFF IN CHARGE (neutrino)" << dispatch;
	    if(disabling_lg) detector_interface::get ()->add_disabled_channel (i + 1, evnum, db_run::charge, this);
	  }
	  neutrino_charge_one_lg->Delete ();
	}
      }
        //check how many off and on ordinary channels in charge
      start_lg_off_charge = 0;
      start_lg_on_charge = 0;
      for (int32_t i = 0; i < n[channel]; i++) {
	if (bx_dbi::get ()->get_channel (i+1).is_ordinary ()){
	  if (cummulative_charge_off_lg[i])  start_lg_off_charge ++;
	  else  start_lg_on_charge ++;
	}
      }
      get_message(bx_message::log) << start_lg_off_charge << " ordinary channels OFF IN CHARGE after 10.000 neutrino triggers" << dispatch;
      get_message(bx_message::log) << start_lg_on_charge  << " ordinary channels ON IN CHARGE after 10.000 neutrino trigger" << dispatch;  
   
    }//end of if after 10.000 neutrino triggers




    //search for disabled/off electronics channels
    if(event_trg_type == 0 && evnum > 1000){

        //search for laben boards off every 3 pulser events
          //accummulate channel occupnacy for Nevents
      for(int32_t i = 0; i < n[lbnb]; i++) {
	  cummulative_lbnb_occupancy_p[i] += dec_lbnb_nhits[i];
      }
      int32_t Nevents_lbnb_p = 3;	
      if(!(n_triggers[0] % Nevents_lbnb_p)) {
	for (int32_t i = 0; i < n[lbnb]; i++) {
	  //lbnb is now off 
	  if(cummulative_lbnb_occupancy_p[i] == 0){
	    int32_t lg_in_lbnb = i * constants::laben::channels_per_board + 1;
	    //check if the lbnb is already off
	    int32_t N_is_off_already = 0;
	    int32_t N_is_ordinary_in_lbnb = 0;
	    for (int32_t j = 0; j < 8; j++) {
	      if(cummulative_off_lg[lg_in_lbnb - 1 + j] && bx_dbi::get ()->get_channel (lg_in_lbnb + j).is_ordinary () ) N_is_off_already ++;
	      if (bx_dbi::get ()->get_channel (lg_in_lbnb +j).is_ordinary ()) N_is_ordinary_in_lbnb ++;
	    }
	    //if the whole board is not yet off, disbale it
	    if(N_is_off_already != N_is_ordinary_in_lbnb ){
	      get_message(bx_message::warn) << "Laben board " << i + 1 << " (channels " << lg_in_lbnb << "-" << lg_in_lbnb + 7 << ") is off from event " << evnum << dispatch; 
	      for (int32_t j = 0; j < 8; j++) {
		cummulative_off_lg[lg_in_lbnb - 1 + j] = 1;
		cummulative_charge_off_lg[lg_in_lbnb - 1 + j] = 1;
		if(disabling_lg) detector_interface::get ()->add_disabled_channel (lg_in_lbnb + j, evnum, db_run::timing, this);
	      }
	    }
	  }
	}
	cummulative_lbnb_occupancy_p.clear ();
	cummulative_lbnb_occupancy_p.resize(n[lbnb],0);
      }
	 
        //search for single channels off in pulser every Nevents_p pulser events
        //sometimes a pulser hit does miss at a single channel         
	 int32_t Nevents_p = 10;	
        //accummulate channel occupnacy for Nevents
      for(int32_t i = 0; i < n[channel]; i++) {
	if (bx_dbi::get ()->get_channel (i+1).is_ordinary ()){
	  cummulative_occupancy_p[i] += dec_nhits_in_channel[i];
	}
      }
    
        //if I accummulated already Nevents from pulser triggers
      if(!(n_triggers[0] % Nevents_p)) {	
	for (int32_t i = 0; i < n[channel]; i++) {
	  if (bx_dbi::get ()->get_channel (i+1).is_ordinary ()){
	    
	    //if lg is now off 
	    if (cummulative_occupancy_p[i] == 0) {
	        //if channel was not yet found as off
	      if(cummulative_off_lg[i] == 0) {
		get_message(bx_message::warn) << "Channel " << i + 1 << " is off (pulser) from event " << evnum << dispatch; 
	        cummulative_off_lg[i] = 1;
	        cummulative_charge_off_lg[i] = 1;
		if(disabling_lg) detector_interface::get ()->add_disabled_channel (i + 1, evnum, db_run::timing, this);	    
	      }
	      //it was off also before, so no change
	      //lg was on (could have some previous off periods though) 
	      if(latest_off_lg_p[i] == 0) {
		latest_off_lg_p[i] = 1;
		lg_Nx_to_off_p[i] ++;
		h_lg_Nx_changed[0][0]->Fill(i+1);
		ev_lg_to_off_p[i] = evnum;
	      }
	    }
	    else if (cummulative_occupancy_p[i] > 0) {//channel is ON
	      //lg was on also before -> no change
	      //but when it was off before:
	      if (latest_off_lg_p[i] == 1){
		latest_off_lg_p[i] = 0;
		lg_Nx_to_on_p[i] ++;
		h_lg_Nx_changed[0][1]->Fill(i+1);
		ev_lg_to_on_p[i] = evnum;
	      }
	    }
	  }
	}
	cummulative_occupancy_p.clear ();
	//2nd input of resize gives the value of newly added inputs ONLY
	cummulative_occupancy_p.resize (n[channel],0);
      }
          
    }   

    //for HV triggers (random, laser, neutrino) we check every Nevents_HV triggers 
    if (event_trg_type > 0 && evnum > 1000){
   
      //accummulate hvb occupancy for Nevents_HVB
      for(int32_t i = 0; i < n[hvb]; i++) {
        cummulative_hvb_occupancy_HV[i] += dec_hvb_nhits[i];
      }
      
      int32_t Nevents_HVB = 100;  
        //if I accummulated already Nevents_HVB from HV triggers
      if(!((n_triggers[1] + n_triggers[2]+ n_triggers[3]) % Nevents_HVB)) {
          //search for HV boards off each  Nevents_HVB  events
	for (int32_t i = 0; i < n[hvb]; i++) {
	  //hvb is now off
	  if(cummulative_hvb_occupancy_HV[i] == 0){
	    int32_t crate = i/7;
	    int32_t board_in_crate = i%7;
	    int32_t lg_in_hvb = crate*160 + board_in_crate*24  + 1;
	    //check if the hvb is already off
	    int32_t N_off_already = 0;
	    int32_t N_ordinary_in_hvb = 0;
	     for (int32_t j = 0; j < (board_in_crate == 6 ? 16 : 24); j++) {
	      if(cummulative_off_lg[lg_in_hvb - 1 + j] && bx_dbi::get ()->get_channel (lg_in_hvb + j).is_ordinary ()) N_off_already ++;
	      if (bx_dbi::get ()->get_channel (lg_in_hvb + j).is_ordinary ()) N_ordinary_in_hvb ++;
	    }
	    //if the whole board is not yet off, disbale it
	    if(N_off_already != N_ordinary_in_hvb){
	      get_message(bx_message::warn) << "HVB board " << i + 1 << " (channels " << lg_in_hvb << "-" << lg_in_hvb + 23 << ") is off from event " << evnum << dispatch; 
	      for (int32_t j = 0; j < (board_in_crate == 6 ? 16 : 24); j++) {
		//FIX ?? 
		cummulative_off_lg[lg_in_hvb - 1 + j] = 1;
		cummulative_charge_off_lg[lg_in_hvb - 1 + j] = 1;
		if(disabling_lg) detector_interface::get ()->add_disabled_channel (lg_in_hvb + j, evnum, db_run::timing, this);
	      }
	    }
	  }
	}
	cummulative_hvb_occupancy_HV.clear ();
	cummulative_hvb_occupancy_HV.resize(n[hvb],0);		
      }

        //accummulate channel occupancy for Nevents
      for(int32_t i = 0; i < n[channel]; i++) {
	if (bx_dbi::get ()->get_channel (i+1).is_ordinary ()) cummulative_occupancy_HV[i] += dec_nhits_in_channel[i];
      }
      
      int32_t Nevents_HV = 12000; //with the rate of 30 Hz about 2000 events/minute
      double dead_fraction = 0.05;

         //if I accummulated already Nevents from HV triggers
      if(!((n_triggers[1] + n_triggers[2]+ n_triggers[3]) % Nevents_HV)) {
	
        //calculate mean for cummulative_occupancy HV
	double mean_cummulative_occupancy_HV = 0;
	int32_t N_lg_on = 0;
	for (int32_t i = 0; i < n[channel]; i++) {
	  if (bx_dbi::get ()->get_channel (i+1).is_ordinary ()){	   
	    if(cummulative_occupancy_HV[i]) {
	      N_lg_on ++;
	      mean_cummulative_occupancy_HV += cummulative_occupancy_HV[i];
	    }
	  }
	}
	if(N_lg_on) mean_cummulative_occupancy_HV = mean_cummulative_occupancy_HV/ N_lg_on;
	  
	for (int32_t i = 0; i < n[channel]; i++) {
	  if (bx_dbi::get ()->get_channel (i+1).is_ordinary ()){
	    //if lg is now off (less than dead_fraction of the mean
	    if (cummulative_occupancy_HV[i] <= dead_fraction * mean_cummulative_occupancy_HV) {
	      if(cummulative_off_lg[i] == 0) {
		get_message(bx_message::warn) << "Channel " << i + 1 << " is off (HV) from event " << evnum << dispatch; 
	        cummulative_off_lg[i] = 1;
	        cummulative_charge_off_lg[i] = 1;
		if(disabling_lg) detector_interface::get ()->add_disabled_channel (i + 1, evnum, db_run::timing, this);
	      }
 
	      //lg was off also before -> no change
	      //lg was on (could have some previous off periods though) 
	      if(latest_off_lg_HV[i] == 0) {
		latest_off_lg_HV[i] = 1;
		lg_Nx_to_off_HV[i] ++;
		h_lg_Nx_changed[1][0]->Fill(i+1);
		ev_lg_to_off_HV[i] = evnum - Nevents_HV;
	      }
	    }
	    else if (cummulative_occupancy_HV[i] > dead_fraction * mean_cummulative_occupancy_HV) {//channel is ON, has hits in this event
	      //lg was on also before -> no change
	      //but when it was off before:
	      if (latest_off_lg_HV[i] == 1){
		latest_off_lg_HV[i] = 0;
		lg_Nx_to_on_HV[i] ++;
		h_lg_Nx_changed[1][1]->Fill(i+1);
		ev_lg_to_on_HV[i] = evnum - Nevents_HV;
	      }
	    }
	  }
	}
	cummulative_occupancy_HV.clear ();
	//2nd input of resize gives the value of newly added inputs ONLY
	cummulative_occupancy_HV.resize (n[channel],0);
      }
    }
  } 
      
      
  if(!calibration_mode){

    //  **************************************************************
    //  crates, fe and laben boards occupancy check each n_events
    //  **************************************************************
     
    int32_t  n_events = get_parameter ("n_events").get_int ();;
  
    //check if it is time to analyse occupancy for the trigger type of the actual event
    if ((evnum < 1050 && n_triggers[0] == 1000 && event_trg_type == 0) || (evnum > 1000 && n_triggers[event_trg_type] && !(n_triggers[event_trg_type] % n_events)))  {
      get_message(bx_message::log) << "Total triggers " << total_triggers <<  ", considered triggers " << n_triggers[event_trg_type] << " in trg type " << trigger_type_map[event_trg_type] << dispatch;
      
      for(int32_t module_type = 0; module_type < 5; module_type++) n_on[module_type] = 0; 
      
      //loop crates = 0, hvb = 1, feb = 2, lbnb = 3, channel = 4)
      for(int32_t module_type = 0; module_type < 5; module_type ++){
        
          //mean and rms of the occupancy distributon
	int32_t     n_cycles = get_parameter ("n_cycles").get_int ();
	double  n_rms = get_parameter ("n_rms").get_float ();
	double mean = 0;
	double rms = 0;
	double S = 0;
	double this_modules_on = 0; 	
          // loop on individual crates etc.
	for(int32_t i_module = 1; i_module <= n[module_type]; i_module++){
	  double bin_content =  dec_occupancy[module_type][event_trg_type]->GetBinContent (i_module);
	  if(bin_content) n_on[module_type] ++;
	}
	
	if (n_on[module_type] < (int32_t) 0.5 *  n[module_type]) 
	  get_message(bx_message::error) << "Only " << n_on[module_type] << " modules on for module_type " << module_type_map[module_type] << dispatch;
	else { //if more than 1/2 of modules on
	  for(int32_t i_module = 1; i_module <= n[module_type]; i_module++){
	    double bin_content =  dec_occupancy[module_type][event_trg_type]->GetBinContent (i_module);
	    if(bin_content) {
	      this_modules_on ++; 	
	      double delta  = bin_content - mean;
	      mean = mean + delta / this_modules_on;
	      S = S + delta * (bin_content - mean);
	    }
	  }
	  
	  
	  rms = std::sqrt(S /  (double) this_modules_on);
	  get_message(bx_message::log) << "Mean occupation module type " << module_type_map[module_type]  << " in trigger type  " << trigger_type_map[event_trg_type] << " is " << mean << "  +- " << rms  << dispatch;
	  
	  // truncated mean and rms of the occupancy distributon
	  double trun_mean = 0;
	  double trun_rms = 0;
	  S = 0;
	  this_modules_on = 0; 	
	  // loop on individual crates etc.
	  for(int32_t i_module = 1; i_module <= n[module_type]; i_module++){
	    double bin_content =  dec_occupancy[module_type][event_trg_type]->GetBinContent (i_module);
	    if(bin_content && bin_content < (mean + 4 * rms)) {
	      this_modules_on ++; 	
	      double delta  = bin_content - trun_mean;
	      trun_mean = trun_mean + delta / this_modules_on;
	      S = S + delta * (bin_content - trun_mean);
	    }
	  }
	  
	  trun_rms = std::sqrt(S /  (double) this_modules_on);
	  get_message(bx_message::log) << "Truncated mean occupation module type " << module_type_map[module_type]  << " in trigger type  " << trigger_type_map[event_trg_type] << " is " << trun_mean << "  +- " << trun_rms  << dispatch;
	  
	  mean = trun_mean;
	  rms = trun_rms;
	  
	  // check for noisy individual feb and lbnb
	  if(module_type == 2 || module_type == 3){
	    for(int32_t i_module = 1; i_module <= n[module_type]; i_module ++){
	      double bin_content = dec_occupancy[module_type][event_trg_type]->GetBinContent (i_module);
	      if(bin_content > (mean + n_rms * rms)) {//is noisy
		(times_found_noisy[module_type][event_trg_type])[i_module - 1] ++;
		//write warning if found for the first time
		if (first_time_hot && (times_found_noisy[module_type][event_trg_type])[i_module - 1] == 1)
		  get_message(bx_message::warn) << "dec hits:  noisy  " <<  module_type_map[module_type] << ", number "<< i_module << "  in trg_type  " << trigger_type_map[event_trg_type] << ", " << (bin_content - mean)/ rms  << " rms above mean, found first time, last event analysed " << evnum  << dispatch;
		//write warning if found again after n_cycles times
		if ( (times_found_noisy[module_type][event_trg_type])[i_module - 1]  &&  !((times_found_noisy[module_type][event_trg_type])[i_module - 1] % n_cycles))
		  get_message(bx_message::warn) << "WAKE UP: dec hits:  noisy " <<  module_type_map[module_type] << ", number "<< i_module << "  in trg_type  " << trigger_type_map[event_trg_type] << ", " << (bin_content - mean)/ rms  << " rms above mean, found already " << (times_found_noisy[module_type][event_trg_type])[i_module - 1] << " times, last event analysed " << evnum <<  dispatch;
	      }	
	    }
	  }
	  
	  // check for off crates and hvb 
	  if(module_type == 0 || module_type == 1){
	    for(int32_t i_module = 1; i_module <= n[module_type]; i_module ++){
	      double bin_content = dec_occupancy[module_type][event_trg_type]->GetBinContent (i_module);
	      if(bin_content < 0.001 * mean ) {//is empty
		(times_found_empty[module_type][event_trg_type])[i_module - 1] ++;
		//write warning if found for the first time
		if ((times_found_empty[module_type][event_trg_type])[i_module - 1] == 1)
		  get_message(bx_message::warn) << "NO HITS in " <<  module_type_map[module_type] << ", number "<< i_module << " in trg_type " << trigger_type_map[event_trg_type] << " found first time, last event analysed " << evnum  << dispatch;
		//write warning if found again after n_cycles times
		if ( (times_found_empty[module_type][event_trg_type])[i_module - 1]  &&  !((times_found_empty[module_type][event_trg_type])[i_module - 1] % n_cycles))
		  get_message(bx_message::error) << "WAKE UP: NO  HITS in " <<  module_type_map[module_type] << ", number "<< i_module << "  in trg_type  " << trigger_type_map[event_trg_type]  << " found already " << (times_found_noisy[module_type][event_trg_type])[i_module - 1] << " times, last event analysed " << evnum <<  dispatch;
	      }	
	    }
	  }
	  
	}
    
	barn_interface::get ()->network_send (dec_occupancy[module_type][event_trg_type], this);
	dec_occupancy[module_type][event_trg_type]->Reset ();  
      }
    }
  

    //neutrino trigger rate and nhits rate for neutrino triggers
    
    double mean_trigger_rate = -10;
    double mean_neutrino_nhits = -10;
 
    if (event_trg_type == 2 && n_triggers[2] == 1) { //first neutrino event
      uint32_t bo;
      ev->get_trigger().get_gps_time ( n_very_first_event_tsec, bo);
      n_events_in_this_cycle = 1;
      n_nhits_in_cycle = dec_nhits_one_event;
      n_time_last_previous_cycle = 0;
    } else if  (event_trg_type == 2 && n_triggers[2] > 1){
      
      n_events_in_this_cycle ++;
      
      uint32_t tsec, tnsec;
      ev->get_trigger().get_gps_time (tsec, tnsec);
      tsec -=  n_very_first_event_tsec;
      double time_this_event =  tsec + (double) tnsec * 1E-9 ; //event time in sec with respect to to the first enutrino event
      //create dynamic tlist of times of the last 300 neutrino events (gps time with respect to the first neutrino event) 
      if (tlist.size () < 300) {
        //time list
	tlist.push_back (time_this_event);
        //nhits list
	neutrino_nhits_list.push_back (dec_nhits_one_event);
	n_nhits_in_cycle += dec_nhits_one_event;
      } else {
        //time list
	tlist.push_back (time_this_event);
	tlist.pop_front ();
	//nhits list
	n_nhits_in_cycle = n_nhits_in_cycle + dec_nhits_one_event - *(neutrino_nhits_list.begin ());
	neutrino_nhits_list.push_back (dec_nhits_one_event);
	neutrino_nhits_list.pop_front ();
      }
      
      // when at least 100 new neutrino events and minimum 60 s
      if ( (time_this_event - n_time_last_previous_cycle) > 61 &&  n_events_in_this_cycle > 100) {
	
        //time list
	double first_time = *(tlist.begin ());
	double last_time = *(tlist.rbegin ());
	mean_trigger_rate = tlist.size ()/(last_time - first_time);
	
	//neutrino trigger rate
	neutrino_trigger_rate->Fill(time_this_event - (time_this_event - n_time_last_previous_cycle)/2, mean_trigger_rate);
	neutrino_trigger_rate->SetAxisRange (0,1.2*time_this_event);
	barn_interface::get ()->network_send (neutrino_trigger_rate, this);
	
	//neutrino nhits per event
	double mean_neutrino_nhits_per_event;
	mean_neutrino_nhits_per_event =  n_nhits_in_cycle /  (double) tlist.size ();
	neutrino_nhits_per_event->Fill(time_this_event - (time_this_event - n_time_last_previous_cycle)/2 , mean_neutrino_nhits_per_event);
	neutrino_nhits_per_event->SetAxisRange (0,1.2*time_this_event);
	barn_interface::get ()->network_send (neutrino_nhits_per_event, this);

	//neutrino nhits rate
	mean_neutrino_nhits =  n_nhits_in_cycle / (last_time - first_time) ;
	neutrino_nhits_rate->Fill(time_this_event - (time_this_event - n_time_last_previous_cycle)/2, mean_neutrino_nhits);
	neutrino_nhits_rate->SetAxisRange (0,1.2*time_this_event);
	barn_interface::get ()->network_send (neutrino_nhits_rate, this);
	n_events_in_this_cycle = 0;      
	
	n_time_last_previous_cycle = last_time;
	n_trgrate_cycle ++;
      }
    }
    
    //random nhits rate
    if (event_trg_type == 1 && n_triggers[1] == 1) { //first random event
      uint32_t bo;
      ev->get_trigger().get_gps_time ( r_very_first_event_tsec, bo);
      r_events_in_this_cycle = 1;
      r_nhits_in_cycle = dec_nhits_one_event;
      r_time_last_previous_cycle = 0;
    } else if  (event_trg_type == 1 && n_triggers[1] > 1){
      
      r_events_in_this_cycle ++;
      uint32_t tsec, tnsec;
      ev->get_trigger().get_gps_time (tsec, tnsec);
      tsec -=  r_very_first_event_tsec;
      double time_this_event =  tsec + (double) tnsec * 1E-9 ; //event time in sec with respect to to the first random event
      //create dynamic tlist of times of the last 300 random events (gps time with respect to the first random event) 
      if (random_tlist.size () < 300) {
	//time list
	random_tlist.push_back (time_this_event);
        //nhits list
	random_nhits_list.push_back (dec_nhits_one_event);
	r_nhits_in_cycle += dec_nhits_one_event;
      } else {
	//time list
	random_tlist.push_back (time_this_event);
	random_tlist.pop_front ();
        //nhits list
	r_nhits_in_cycle = r_nhits_in_cycle + dec_nhits_one_event - *(random_nhits_list.begin ());
	random_nhits_list.push_back (dec_nhits_one_event);
	random_nhits_list.pop_front ();
      }
      
      // when at least 100 new random events and minimum 240 s
      if ( (time_this_event - r_time_last_previous_cycle) > 240 &&  r_events_in_this_cycle > 100) {
	
	//time list
	double r_first_time = *(random_tlist.begin ());
	double r_last_time = *(random_tlist.rbegin ());
	
	//nhits per event
	double mean_random_nhits_per_event;
	mean_random_nhits_per_event =  r_nhits_in_cycle /  (double) random_tlist.size ();
	random_nhits_per_event->Fill(time_this_event - (time_this_event - r_time_last_previous_cycle)/2 , mean_random_nhits_per_event);
	random_nhits_per_event->SetAxisRange (0,1.2*time_this_event);
	barn_interface::get ()->network_send (random_nhits_per_event, this);
	
	//random nhits rate
	double mean_random_nhits =  r_nhits_in_cycle / (r_last_time - r_first_time);
	random_nhits_rate->Fill(time_this_event  - (time_this_event - r_time_last_previous_cycle)/2, mean_random_nhits);
	random_nhits_rate->SetAxisRange (0,1.2*time_this_event);
	barn_interface::get ()->network_send (random_nhits_rate, this);
      
	r_events_in_this_cycle = 0;
	r_time_last_previous_cycle =  r_last_time;
	r_trgrate_cycle ++;
      }
    }
  
    
    //pulser nhits rate
    
    if (event_trg_type == 0 ) { 
      //int32_t bin = evnum/(max_events/pulser_bins);
      nhits_pulser_sum += dec_nhits_one_event;
      int32_t pulser_events_in_cycle = 1;
      if(!(n_triggers[0] %  pulser_events_in_cycle)) {
	pulser_nhits->SetAxisRange (0,evnum+100);
	pulser_nhits->SetMinimum (0.95 * nhits_pulser_sum/ pulser_events_in_cycle);
	pulser_nhits->Fill(evnum, nhits_pulser_sum/ pulser_events_in_cycle);
	for(int32_t i = 1; i < pulser_bins; i ++) pulser_nhits->SetBinError(i,1);
	(barn_interface::get ()->network_send (pulser_nhits, this));
	nhits_pulser_sum = 0;
      }
    }
  }
  
  return ev;
  }

//END
void bx_detector_monitor::end () {
  get_message(bx_message::debug) << "end" << dispatch;    
  
  
  //check how many ordinary lg are off and on at the end
  end_lg_on = 0;
  end_lg_off = 0;
  end_lg_on_charge = 0;
  end_lg_off_charge = 0;
  for (int32_t i = 0; i < n[channel]; i++) {
    if (bx_dbi::get ()->get_channel (i+1).is_ordinary ()){
      if (cummulative_off_lg[i] == 1) end_lg_off ++;
      else if (cummulative_off_lg[i] == 0) end_lg_on ++;
      if (cummulative_charge_off_lg[i] == 1) end_lg_off_charge ++;
      else if (cummulative_charge_off_lg[i] == 0) end_lg_on_charge ++;
    }
  }
  
  get_message(bx_message::info) << "At the end, ordinary channels: " << end_lg_on << " on, " << end_lg_off << " off" << dispatch;
  get_message(bx_message::info) << "During the run " << start_lg_on - end_lg_on << " ordinary channels lost " << dispatch;
  get_message(bx_message::info) << "At the end, ordinary channels IN CHARGE: " << end_lg_on_charge << " on, " << end_lg_off_charge << " off" << dispatch;
  get_message(bx_message::info) << "During the run " << start_lg_on_charge - end_lg_on_charge << " ordinary channels lost IN CHARGE" << dispatch;
  

    //if allowed from echidna.cfg, write to DB (if not already present, checked in db_run)
    db_run& run_info = bx_dbi::get()->get_run ();
    if (get_parameter ("db_write").get_bool ()) run_info.write_disabled_channels (true, this);


  /*  if(calibration_mode){   
    
      //to print off lg 
    for (int32_t i = 0; i < n[channel]; i++){
      int32_t indx = 0;
      
      //lg with missing precalibration hits
      

      //lg turned off once in both pulser and HV triggers (probably disabled from DAQ) 
      if(latest_off_lg_p[i] && latest_off_lg_HV[i] && lg_Nx_to_off_p[i] == 1 &&  lg_Nx_to_off_HV[i] == 1){
	get_message(bx_message::warn) << "Lg " << i+1 << " is off from the event " <<  ev_lg_to_off_p[i] << " (pulser) and event " <<  ev_lg_to_off_HV[i] << " in HV triggers" << dispatch;
	indx = 1;
      }
      
      //lg turned off once in only pulser (should not happen!)
      if(latest_off_lg_p[i] && lg_Nx_to_off_p[i]== 1 &&  lg_Nx_to_on_p[i]== 0 && latest_off_lg_HV[i]==0 &&  lg_Nx_to_off_HV[i] == 0){
	get_message(bx_message::warn) << "Lg " << i+1 << " is off only in pulser from event " <<  ev_lg_to_off_p[i] << dispatch;
	indx = 1;
      }
      
      //lg turned off once in HV trigger only (HV tripped?)
    if(latest_off_lg_HV[i] && lg_Nx_to_off_HV[i] == 1 && lg_Nx_to_on_HV[i] == 0 && latest_off_lg_p[i] == 0 && lg_Nx_to_off_p[i] == 0){
      get_message(bx_message::warn) << "Lg " << i+1 << " is off only in HV triggers from event " <<  ev_lg_to_off_HV[i] << dispatch;
      indx = 1;
    }
    
      //is now on in p and HV, but had an off period in pulser only (temporary lost pulser events)- should not really  happen
    if(latest_off_lg_p[i] == 0 && latest_off_lg_HV[i] == 0 && lg_Nx_to_off_p[i] > 0 && lg_Nx_to_off_HV[i] == 0){
      get_message(bx_message::warn) << "Lg " << i+1 << " is on in pulser and HV, but was off in pulser " <<  lg_Nx_to_off_p[i] << " times (last time from  event " << ev_lg_to_off_p[i] << " ) and got on " << lg_Nx_to_on_p[i] << " times (last time at event "  << ev_lg_to_on_p[i] << ")" << dispatch;
      indx = 1;
    }

    //is now on in p and HV, but had an off period in HV only (temporary HV failure)
    if(latest_off_lg_p[i] == 0 && latest_off_lg_HV[i] == 0 && lg_Nx_to_off_p[i] == 0 && lg_Nx_to_off_HV[i] > 0){
      get_message(bx_message::warn) << "Lg " << i+1 << " is on in pulser and HV, but was off in HV triggers " <<  lg_Nx_to_off_HV[i] << " times (last time at event " << ev_lg_to_off_HV[i] << " ) and got on " << lg_Nx_to_on_HV[i] << " times (last time at event "  << ev_lg_to_on_HV[i] << ")" << dispatch;
      indx = 1;
    }
    
    //is now on in p and HV, but had an off period in both pulser and HV 
    if(latest_off_lg_p[i] == 0 && latest_off_lg_HV[i] == 0 && lg_Nx_to_off_p[i] > 0 && lg_Nx_to_off_HV[i] > 0){
       get_message(bx_message::warn) << "Lg " << i+1 << " is on in pulser and HV, but was off in HV triggers" <<  lg_Nx_to_off_HV[i] << " times (last time at event " << ev_lg_to_off_HV[i] << " ) and got on " << lg_Nx_to_on_HV[i] << " times (last time at event "  << ev_lg_to_on_HV[i] << ") and in pulser " <<  lg_Nx_to_off_p[i] << " times (last time from  event " << ev_lg_to_off_p[i] << " ) and got on " << lg_Nx_to_on_p[i] << " times (last time at event "  << ev_lg_to_on_p[i] << ")" << dispatch;
      indx = 1;
    }
    
    //some other combinations
    if((lg_Nx_to_off_p[i] > 0 || lg_Nx_to_on_p[i] > 0 ||  lg_Nx_to_off_HV[i] > 0 || lg_Nx_to_on_HV[i] > 0) && indx == 0)
      get_message(bx_message::warn) << "Lg " << i+1 << " strange: Now: off in pulser ?: " << latest_off_lg_p[i] << " , Is off is HV ?: " << latest_off_lg_HV[i] << " changed to off in pulser " << lg_Nx_to_off_p[i] << " times and to on in pulser " <<  lg_Nx_to_on_p[i] << " and changed to off in HV triggers " << lg_Nx_to_off_HV[i] << " times and to on in HV triggers " <<  lg_Nx_to_on_HV[i] << dispatch;
    
    }    
  }
  */
}
  
    


