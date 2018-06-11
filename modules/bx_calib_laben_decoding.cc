/* BOREXINO Reconstruction program
 *
 * Author: Livia Ludhova <livia.ludhova@mi.infn.it>
 * 	   Alessandro Fontana <alessandro.fontana@mi.infn.it>
 * 	   Davide Franco <davide.franco@mi.infn.it>:(
 * Maintainer: Livia Ludhova <livia.ludhova@mi.infn.it>
 *
 * $Id: bx_calib_laben_decoding.cc,v 1.35 2012/05/09 17:49:34 ludhova Exp $
 */

#include "bx_calib_laben_decoding.hh"
#include "bx_echidna_event.hh"
#include "barn_interface.hh"
#include "db_run.hh"
#include "db_profile.hh"
#include "bx_dbi.hh"
#include "TH1F.h"
#include "TH2F.h"

// ctor
bx_calib_laben_decoding::bx_calib_laben_decoding (): bx_base_module("bx_calib_laben_decoding", bx_base_module::main_loop) {
  require_event_stage (bx_detector::laben, bx_base_event::decoded);
  require_trigger_type (bx_trigger_event::pulser);
}


// BEGIN
void bx_calib_laben_decoding::begin () {

  //parameters from echidna.cfg
  check_reference_lg  = get_parameter ("check_reference_lg").get_int ();
  limit_nlg_precalib_not_ok  =  get_parameter("limit_nlg_precalib_not_ok").get_int ();
  mean_limit_single = get_parameter ("mean_limit_single").get_float ();
  rms_limit_single =  get_parameter ("rms_limit_single").get_float (); 
  mean_limit_total = get_parameter ("mean_limit_total").get_float ();
  rms_limit_total =  get_parameter ("rms_limit_total").get_float (); 
  monitor_mode =  get_parameter ("monitor_mode").get_bool (); 
  

  offset  = new TH1F ("offset", "offset of individual channels ", 300, -150, 150); 
  barn_interface::get ()->store (barn_interface::file, offset, this);

  rms  = new TH1F ("rms", "rms of individual channels", 500, 0, 50); 
  barn_interface::get ()->store (barn_interface::file, rms, this);

  laben_decoding  = new TH1F ("laben_decoding", "times", 2000, -1000, 1000); //range influences the offset value!
  barn_interface::get ()->store (barn_interface::junk, laben_decoding, this);
    
  charge_vs_channel = new TH2F ("charge_vs_channel","charge decoded hits", constants::laben::channels, 1, constants::laben::channels + 1, 512, -256, 256);
  barn_interface::get ()->store (barn_interface::file, charge_vs_channel, this);
  charge_vs_channel->SetXTitle("logical channel number");
  charge_vs_channel->SetYTitle("ADC channel ");


  time_vs_channel = new TH2F ("time_vs_channel","pulser triggers, decoded hits", constants::laben::channels, 1, constants::laben::channels + 1, 600, -300, 300);
  time_all_vs_channel = new TH2F ("time_all_vs_channel","pulser triggers, all decoded hits", constants::laben::channels, 1, constants::laben::channels + 1, 2650, -300, 5000);
  barn_interface::get ()->store (barn_interface::file, time_vs_channel, this);
  time_vs_channel->SetXTitle("logical channel number");
  time_vs_channel->SetYTitle("raw_time - mean - mode ");

  barn_interface::get ()->store (barn_interface::file, time_all_vs_channel, this);
  time_all_vs_channel->SetXTitle("logical channel number");
  time_all_vs_channel->SetYTitle("raw_time - mean - mode ");

  trigger_sum = 0;
  b_found_mctruth = false;
}


//DOIT
bx_echidna_event* bx_calib_laben_decoding::doit (bx_echidna_event *ev) {
  const bx_laben_event& er = ev->get_laben ();
  if (ev->is_mctruth_enabled ()) b_found_mctruth = true;

  trigger_sum++;

  if(trigger_sum == 1000 && monitor_mode)
    end_1000ev (0);

  double mean = 0;
  int32_t hits_in_gate = 0;
  int32_t nhits = er.get_decoded_nhits ();
  int32_t i;

   for(i = 0; i < nhits; i++) {
    if(!er.get_decoded_hit (i).is_out_of_gate ()){
      mean += er.get_decoded_hit (i).get_raw_time ();
      hits_in_gate ++;
    }
   }
   mean = mean / hits_in_gate;

   for(i = 0; i < nhits; i++){
     if(!er.get_decoded_hit (i).is_out_of_gate ()){
       double t = er.get_decoded_hit (i).get_raw_time () - mean;
       laben_decoding->Fill (t);
     }
   }

   double offset = laben_decoding->GetMaximumBin () - 1000.;
   laben_decoding->Reset (); 
  
   for(i = 0; i < nhits; i++){
     if(!er.get_decoded_hit (i).is_out_of_gate ()){
      int32_t fill = 0;
      double normalized_time = er.get_decoded_hit (i).get_raw_time () - mean - offset ;
      double charge = er.get_decoded_hit (i).get_uncorrected_charge_bin ();
      double charge1 = er.get_decoded_hit (i).get_charge_bin ();
      int32_t lg = er.get_decoded_hit (i).get_raw_hit ().get_logical_channel ();
      charge_vs_channel->Fill(lg,charge1);
         //for ordinary channels filling all hits
         //for laser and trigger reference filling only hits with charge correpsonding to pulser hits
      if (bx_dbi::get ()->get_channel (lg).is_ordinary ()) fill = 1;
      if (bx_dbi::get ()->get_channel (lg).is_laser() && charge > 30 && charge < 75) fill = 1;
      if (bx_dbi::get ()->get_channel (lg).is_trigger() && charge > 30 && charge < 75) fill = 1;
      //cngs reference channels
      if (!bx_dbi::get ()->get_channel (lg).is_ordinary () && !bx_dbi::get ()->get_channel (lg).is_empty () && !bx_dbi::get ()->get_channel (lg).is_trigger () && !bx_dbi::get ()->get_channel (lg).is_laser () && charge > 30 && charge < 75) fill = 1;

      
      if(fill) time_vs_channel->Fill (lg, normalized_time);
      time_all_vs_channel->Fill (lg, normalized_time);
     }
  }
  
   if(!(int32_t(trigger_sum) % 200)) barn_interface::get ()->network_send (time_vs_channel, this);


  return ev;
   

}

//END_1000Ev
void bx_calib_laben_decoding::end_1000ev (int32_t indx) {
  
  if(trigger_sum == 0) return;
 
  
  int32_t n_off_channels = 0;
  int32_t n_bad_channels = 0;
  
  int32_t n_trigref = 0;
  int32_t n_off_trigref = 0;
  int32_t n_bad_trigref = 0;
  
  int32_t n_lasref = 0;
  int32_t n_off_lasref = 0;
  int32_t n_bad_lasref = 0;

  int32_t n_cngsref = 0;
  int32_t n_off_cngsref = 0;
  int32_t n_bad_cngsref = 0;
  
  //if precalib not yet present, read the list of bad and off lg from previous run 
  // if ( !bx_dbi::get ()->get_run ().is_precalib_quality_present ()){
  //  db_run& run_info = bx_dbi::get ()->get_run ();
  //  prev_bad_channels_list = run_info.get_laben_precalib_bad_channels ();
  //  prev_off_channels_list = run_info.get_laben_precalib_off_channels (); 
  // }
  
  //this run 
  bad_channels_list.clear ();
  off_channels_list.clear ();
    
    //copy the time_vs_channel to time_vs_channel_good, later it will be cleaned from bad lg
  time_vs_channel_good = (TH2F*) time_vs_channel->Clone (); 
  time_vs_channel_good->SetName ("time_vs_channel_good");
  if (indx == 1) barn_interface::get ()->store (barn_interface::file, time_vs_channel_good, this);
   
    // loop for all logical channels
  for(int32_t ch = 1;ch < (constants::laben::channels + 1);ch++){
    double total_sum = 0; //sum in 602 time bins,including underflow and overflow   
    double peak_sum = 0; //number of events in the peak around 0 (-55,55)ns 
    peak_mean = 0;
    peak_rms = 0; 
    double outer_sum = 0;  //events in time channels dt=(-inf,-54) + (200,inf), including overflow, underflow 
    int32_t indx_bad = 0;  //set to 1 if channel is found as problematic
    int32_t indx_off = 0;  //set to 1 if channel is found as problematic
    //int32_t was_off = 0;
    //int32_t was_bad = 0;
    int32_t lg_to_check = 0;
    int32_t lg_ordinary = 0;
    int32_t lg_laser = 0;    
    int32_t lg_trigger = 0;    
    int32_t lg_cngs = 0;    
    

    if(bx_dbi::get ()->get_channel (ch).is_ordinary()) {
      lg_to_check = 1;
      lg_ordinary = 1;
    }

    if(check_reference_lg && bx_dbi::get ()->get_channel (ch).is_laser()){
      lg_to_check = 1;
      lg_laser = 1;
      n_lasref ++; 
    }
    
    if(check_reference_lg && bx_dbi::get ()->get_channel (ch).is_trigger()){
      lg_to_check = 1;
      lg_trigger = 1;
      n_trigref ++;
    }

    //cngs reference channels are considered as not empty
    if(check_reference_lg && !(bx_dbi::get ()->get_channel (ch).is_empty()) && !(bx_dbi::get ()->get_channel (ch).is_ordinary ()) && !(bx_dbi::get ()->get_channel (ch).is_laser()) && !(bx_dbi::get ()->get_channel (ch).is_trigger ())  )  {
      lg_to_check = 1;
      lg_cngs = 1;
      n_cngsref ++;
    }

    if(lg_to_check) {


        //integral of entries for one logical channels, including overflow and underflow
      total_sum = time_vs_channel->Integral (ch,ch,0,601); 
 
        //integral of the central part for 1 log. channel, time channels time =(-55,55)ns around 0
      peak_sum = time_vs_channel->Integral (ch,ch,245,355); 
    
        //integral in time channels (-inf,-54)ns + (200,inf)ns, including overflow, underflow (excluding retriggering region)
      outer_sum = time_vs_channel->Integral (ch,ch,0,244) + time_vs_channel->Integral (ch, ch, 500, 601);
       
      if((total_sum/trigger_sum) < 0.01) {//channel has  no statistics 

	if(lg_ordinary) get_message(bx_message::log) << "ordinary logical channel " << ch << " has less hits than 1% of triggers" << dispatch;
	if(lg_laser) {
	  get_message(bx_message::log) << "laser reference channel " << ch << " has less hits than 1% of triggers" << dispatch;
	  n_off_lasref ++;
	}
	if(lg_trigger){
	  get_message(bx_message::log) << "trigger reference channel " << ch << " has less hits than 1% of triggers" << dispatch;
	  n_off_trigref ++;
	}	
	if(lg_cngs){
	  get_message(bx_message::log) << "cngs reference channel " << ch << " has less hits than 1% of triggers" << dispatch;
	  n_off_cngsref ++;
	}	

	off_channels_list.push_back (ch);
	indx_off = 1;
      }
      
      else{//channel has statistics 
    
	  //to find the mean and rms of the events in the time interval (-55,55) ns
        central_one_channel = time_vs_channel->ProjectionY ("central_one_channel",ch,ch); 
        central_one_channel->SetAxisRange (-55,55); 
        peak_mean = central_one_channel->GetMean ();
	peak_rms = central_one_channel->GetRMS ();
	offset->Fill(peak_mean);
	rms->Fill(peak_rms);
 
       	  // to find the channels with not goog mean or offset (not centered, or too braod, or too thin == with few hits)
	if( fabs(peak_mean) > mean_limit_single || peak_rms > rms_limit_single || peak_rms < 0.2){ 
  
	  if(lg_ordinary) get_message (bx_message::log) << "problem with ordinary logical channel " << ch << " , peak_mean is " << peak_mean << " and peak_rms is " << peak_rms << dispatch;	  
	  if(lg_laser){
	    get_message (bx_message::log) << "problem with laser reference channel " << ch << " , peak_mean is " << peak_mean << " and peak_rms is " << peak_rms << dispatch;	
	    n_bad_lasref ++;
	  }
	  if(lg_trigger){
	    get_message (bx_message::log) << "problem with trigger reference channel " << ch << " , peak_mean is " << peak_mean << " and peak_rms is " << peak_rms << dispatch;	  
	    n_bad_trigref ++;
	  }
	  if(lg_cngs){
	    get_message (bx_message::log) << "problem with cngs reference channel " << ch << " , peak_mean is " << peak_mean << " and peak_rms is " << peak_rms << dispatch;	  
	    n_bad_cngsref ++;
	  }
	  bad_channels_list.push_back (ch);
	  indx_bad = 1; 
        }
	 
	else { //too many events outside the central region?
	  if( (peak_sum/trigger_sum < 0.60 && outer_sum/trigger_sum > 0.20) ||  peak_sum/(peak_sum + outer_sum) < 0.5 ) {
	    if(lg_ordinary) get_message (bx_message::log) << "ordinary channel " << ch << " too many events outside central region" << dispatch;
	    if(lg_laser){
	      get_message (bx_message::log) << "laser reference channel " << ch << " too many events outside central region" << dispatch;
	      n_bad_lasref ++;
	    }
	    if(lg_trigger){
	      get_message (bx_message::log) << "trigger reference channel " << ch << " too many events outside central region" << dispatch;
	      n_bad_trigref ++;
	    }
	    if(lg_cngs){
	      get_message (bx_message::log) << "cngs reference channel " << ch << " too many events outside central region" << dispatch;
	      n_bad_cngsref ++;
	    }
	    bad_channels_list.push_back (ch);
	    indx_bad = 1;          
	  }
	}
      }
    } 

    //we clean time_vs_channel_good histo, for bad and not checked lg the content is set to 0
    if(indx_bad == 1 || indx_off == 1) for(int32_t i = 0;i < 602;i++) time_vs_channel_good->SetBinContent (ch,i,0);

    
    
    if(indx_off == 1){
      n_off_channels ++;
      /*
	if ( !bx_dbi::get ()->get_run ().is_precalib_quality_present ()){
	//look if found off channel was off also in prev run; if not put it to vector of new off lg
	unsigned i = 0; 
	while (i < prev_off_channels_list.size () &&  prev_off_channels_list[i] <= ch) {
	  if (prev_off_channels_list[i] == ch) was_off = 1;
	  i ++;
	}      
	if (was_off == 0) new_off_channels_list.push_back (ch);
      }
      */
    }
    if(indx_bad == 1){
      n_bad_channels ++;
      /*
      if ( !bx_dbi::get ()->get_run ().is_precalib_quality_present ()){
	  //look if found bad channel was bad also in prev run; if not put it to vector of new bad lg
	unsigned i = 0; 
	while (i < prev_bad_channels_list.size () &&  prev_bad_channels_list[i] <= ch) {
	  if (prev_bad_channels_list[i] == ch) was_bad = 1;
	  i ++;
	}      
	if (was_bad == 0) new_bad_channels_list.push_back (ch);
      }
      */
    } 
  }
  
  time_vs_channel_good_central = time_vs_channel_good->ProjectionY ("all_channels_good",1,constants::laben::channels); 
  time_vs_channel_good_central->SetAxisRange (-55,55); 
  rms_clean = time_vs_channel_good_central->GetRMS ();
  mean_clean = time_vs_channel_good_central->GetMean ();
   
  time_vs_channel_central = time_vs_channel->ProjectionY ("all_channels",1,constants::laben::channels); 
  time_vs_channel_central->SetAxisRange (-55,55); 
  rms_old = time_vs_channel_central->GetRMS ();
  mean_old = time_vs_channel_central->GetMean ();
    
    //warning if offset or rms of the detector (considering precalib-ok channels ) is bad
  if(fabs(mean_clean) > mean_limit_total) get_message(bx_message::warn) << "mean_offset (for channels with ok_precalib) is too bad as " << rms_clean << " ns" << dispatch;
  if(rms_clean > rms_limit_total) get_message(bx_message::warn) << "rms_clean (for channels with ok_precalib) is too bad as " << rms_clean << " ns" << dispatch;  
  
  //log about mean and rms
  get_message(bx_message::log) << "rms with all channels: " << rms_old << " ns, rms with only good channels:  " << rms_clean << " ns" << dispatch;
  get_message(bx_message::log) << "mean with all channels: " << mean_old << " ns, mean with only good channels:" << mean_clean << " ns" << dispatch;


    //warning if too many channels have no or bad precalibration
  if ( (n_bad_channels +  n_off_channels) > limit_nlg_precalib_not_ok)
    get_message(bx_message::warn) << "Too many channels with failed or missing precalib: " <<  n_bad_channels + n_off_channels << " channels " << dispatch;
    
  if(check_reference_lg){
      //warnings if too many trigref channels have no or bad precalib
    if ((n_trigref - n_bad_trigref - n_off_trigref) < 1)
      get_message(bx_message::error) << "No trigger reference channel with a good precalibration! (n_off = " << n_off_trigref << ", n_bad = " << n_bad_trigref << ")" << dispatch;
    else if ((n_trigref - n_bad_trigref - n_off_trigref) < (int32_t) n_trigref /2)
      get_message(bx_message::warn) << "More than half of trigref channels have no or bad precalibration (n_off = " << n_off_trigref << ", n_bad = " << n_bad_trigref << ")" << dispatch;
    
 
      //warnings if too many lasref channels have no or bad precalib
    if ((n_lasref - n_bad_lasref - n_off_lasref) < 1)
      get_message(bx_message::error) << "No laser reference channel with a good precalibration! (n_off = " << n_off_lasref << ", n_bad = " << n_bad_lasref << ")" << dispatch;
    else if ((n_lasref - n_bad_lasref - n_off_lasref) < (int32_t) n_lasref /2)
      get_message(bx_message::warn) << "More than half of lasref channels have no or bad precalibration (n_off = " << n_off_lasref << ", n_bad = " << n_bad_lasref << ")" << dispatch;
    
  }
  
  //warning for new bad and off channels
  if ( !bx_dbi::get ()->get_run ().is_precalib_quality_present ()){
    
    /*
    get_message(bx_message::warn) << "Comparing to previous run, new off channels are: " << dispatch;
    for(unsigned i = 0; i < new_off_channels_list.size (); i ++)  get_message(bx_message::warn) << new_off_channels_list[i] << dispatch;
    
    get_message(bx_message::warn) << "Comparing to previous run, new bad channels are: " << dispatch;
    for(unsigned i = 0; i < new_bad_channels_list.size (); i ++)  get_message(bx_message::warn) << new_bad_channels_list[i] << dispatch;
    

    get_message(bx_message::log) << "previous run: number of bad electronics channels: " << prev_bad_channels_list.size () << dispatch;
    get_message(bx_message::log) << "previous run: number of off electronics channels: " << prev_off_channels_list.size () << dispatch;
    */
 }
    
    //logs about number of bad and off channels
   get_message(bx_message::log) << "number of bad electronic channels " << n_bad_channels << ", from checked ordinary and reference (if required) channels " << dispatch;
  get_message(bx_message::log) << "number of off electronic channels " << n_off_channels <<  ", from checked ordinary and reference (if required) channels " << dispatch;

  //logs about number of bad and off reference channels
  if(check_reference_lg){
    get_message(bx_message::log) << "number of off: " << n_off_trigref << " and bad: " << n_bad_trigref << " trigger reference channels" << dispatch;
    get_message(bx_message::log) << "number of off: " << n_off_lasref << " and bad: " << n_bad_lasref << " laser reference channels" << dispatch;
  }

     //after 1000 events
  if (indx == 0) {
    time_vs_channel_good->Delete ();
    if(time_vs_channel->Integral () > 2240 * 100){
      if (::fabs (mean_clean) <= mean_limit_total && rms_clean < rms_limit_total )     
	get_message (bx_message::info) << "precalibrations ok: mean " << std::setprecision (2) << mean_clean << "+-" << std::setprecision (2) << rms_clean << "ns. # bad lg: " <<  n_bad_channels + n_off_channels << dispatch;
      else get_message (bx_message::warn) << "ERROR: precalibrations failed: " << std::setprecision (2) << mean_clean << "+-" << std::setprecision (2) << rms_clean << " ns. # bad lg: " <<  n_bad_channels + n_off_channels << dispatch;
    }
  }

}


//END
void bx_calib_laben_decoding::end () {
  
  end_1000ev (1);
  
  if( time_vs_channel->Integral () > 2240 * 100){
    db_run& run_info = bx_dbi::get ()->get_run ();
    if (::fabs (mean_clean) <= mean_limit_total && rms_clean < rms_limit_total && get_parameter ("db_write").get_bool ()) {
      if (b_found_mctruth) get_message (bx_message::warn) << "precalib not written for simulation" << dispatch;
      else {
	run_info.write_laben_precalib (true, this);
        //if(!run_info.is_precalib_quality_present () && get_parameter ("db_write").get_bool ()){
	run_info.set_laben_precalib_mean_time ( mean_clean, this);
	run_info.set_laben_precalib_sigma_time (rms_clean, this);
	run_info.set_laben_precalib_bad_channels (bad_channels_list, this);
	run_info.set_laben_precalib_off_channels (off_channels_list, this);

	run_info.write_laben_precalib_quality (true, this);
      }
    }
  }
}
/*
 * $Log: bx_calib_laben_decoding.cc,v $
 * Revision 1.35  2012/05/09 17:49:34  ludhova
 * chnage in range of time histo for all hits
 *
 * Revision 1.34  2012-05-09 15:13:10  ludhova
 * time histo for all hits
 *
 * Revision 1.33  2012-04-26 11:17:23  ludhova
 * charge control for cngs referebnce channels
 *
 * Revision 1.32  2012-04-01 14:38:35  ludhova
 * treatment of cngs reference channels
 *
 * Revision 1.31  2008-11-19 11:22:56  ludhova
 * change cut on charge of precalib pulses of reference hits
 *
 * Revision 1.30  2007-09-13 15:01:25  ddangelo
 * modified msg for precalib status to comply with omon needs
 *
 * Revision 1.29  2007-09-13 12:58:56  ludhova
 * cosmetics for online
 *
 * Revision 1.28  2007-09-13 12:32:52  ludhova
 * new status messages after 1000 events
 *
 * Revision 1.27  2007-05-25 11:16:10  razeto
 * bx_calib_laben_decoding will validate precalibrations
 *
 * Revision 1.26  2007/02/09 17:28:17  ludhova
 * comparison with previous run removed...
 *
 * Revision 1.25  2007-01-24 15:07:19  ludhova
 * warning for new off and bad precalibration channels (with respect to previous run, online mode
 *
 * Revision 1.24  2007-01-24 12:17:47  ludhova
 * new histos for mean and rms; small changes in bad-channel condition
 *
 * Revision 1.23  2007-01-24 09:49:07  ludhova
 * additional control of reference channels; cut on charge (30-60) for rerence channels only added
 *
 * Revision 1.22  2006-11-23 16:13:15  ludhova
 * change in outer sum for trigger reference channels
 *
 * Revision 1.21  2006-10-17 15:36:10  ludhova
 * monitor mode enabled
 *
 * Revision 1.20  2006-09-18 17:27:50  ludhova
 * bug correction
 *
 * Revision 1.19  2006-09-18 12:06:50  ludhova
 * better quality controll
 *
 * Revision 1.18  2006-09-14 13:05:34  ludhova
 * new quality controll added
 *
 * Revision 1.17  2006/08/28 16:09:23  ludhova
 * db_write condition added, default 0 (with Ale)
 *
 * Revision 1.16  2006/08/21 11:23:17  razeto
 * Updated to new barn_interface + do not return 0
 *
 * Revision 1.15  2006/07/18 12:11:11  ludhova
 * check if data written in DB added
 *
 * Revision 1.14  2006/07/18 10:24:50  razeto
 * Fixed a bug
 *
 * Revision 1.13  2006/07/13 13:32:22  ludhova
 * a check if data present added
 *
 * Revision 1.12  2006/06/21 12:34:21  ludhova
 * added visitors
 *
 * Revision 1.11  2006/06/20 14:19:55  ludhova
 * output histogram to be sent in Gviewer
 *
 * Revision 1.10  2006/05/17 14:39:31  ludhova
 * warnings changed to log's
 *
 * Revision 1.9  2006/05/09 08:00:08  ludhova
 * variable name changed on Ale suggestion
 *
 * Revision 1.8  2006/05/05 14:50:25  ludhova
 * precalibration checks for reference (trigger, laser) channels added
 *
 * Revision 1.7  2006/03/13 13:09:59  ludhova
 * correction of warnings
 *
 * Revision 1.6  2006/03/13 13:07:32  ludhova
 * correction of warnings
 *
 * Revision 1.5  2006/03/02 14:09:28  ludhova
 * output in th form ready for database
 *
 * Revision 1.4  2006/02/27 11:19:28  ludhova
 * small changes concerning hit out of gate
 *
 * Revision 1.3  2006/01/16 14:30:08  ludhova
 * module checking the precalibration alignment edited
 *
 * Revision 1.2  2006/01/11 13:28:07  fontana
 * *** empty log message ***
 *
 * Revision 1.1  2006/01/11 11:54:22  dfranco
 * added a new module: bx_calib_laben_decoding.hh
 *
 */
