/* BOREXINO Reconstruction program
 *
 * Author: Maria Elena Monzani <monzani@mi.infn.it>
 * Maintainer: Livia Ludhova <ludhova@gmail.com>
 *
 * $Id: bx_calib_laben_time_align.cc,v 1.12 2008/07/01 13:48:38 ludhova Exp $ 
 *
 * Implementation of bx_calib_laben_time_align
 *
 */
#include "bx_calib_laben_time_align.hh"
#include "messenger.hh"
#include "bx_dbi.hh"
#include "db_profile.hh"
#include "db_channel.hh"
#include "db_run.hh"
#include "barn_interface.hh"
#include <fstream>

// ctor
bx_calib_laben_time_align::bx_calib_laben_time_align (): bx_base_module("bx_calib_laben_time_align", bx_base_module::main_loop) {
  require_event_stage (bx_detector::laben, bx_base_event::decoded);
  require_trigger_type (bx_trigger_event::laser394);
}

void bx_calib_laben_time_align::begin () {
  get_message(bx_message::debug) << "begin" << dispatch;
  // creation of the 2-dim histogram 
  int32_t nch = constants::laben::channels;
  bx_time_calib = new TH2S ("bx_time_calib", "Time Calibration Histogram",nch, 1, nch + 1, 1000, 500., 1000.);
  barn_interface::get ()->store (barn_interface::file, bx_time_calib, this);
  // creation of the result vectors
  f_time_offset = new float[nch]; 
  f_time_sigma  = new float[nch]; 
  std::fill_n (f_time_offset,nch,0.);
  std::fill_n (f_time_sigma,nch,0.);
  const db_run& run_info = bx_dbi::get()->get_run ();
  if (get_parameter ("write_calib").get_bool () && run_info.is_laben_laser_present ()) 
    get_message (bx_message::warn) << "Laser calibrations already present for present run !!!" << dispatch; 
}

//DO IT
// filling of the 2-dim histogram
bx_echidna_event* bx_calib_laben_time_align::doit (bx_echidna_event *ev) {
  const bx_laben_event& el = ev->get_laben();
  // Update has data field
  if (el.get_raw_nhits ()) b_has_data = true;
  int32_t nhits = el.get_decoded_nhits();
  for(int32_t i = 0; i < nhits; i++){
    const bx_laben_decoded_hit& hit = el.get_decoded_hit(i);
    const bx_laben_raw_hit& rawhit = hit.get_raw_hit();
    uint16_t ch = rawhit.get_logical_channel();
    double rawt = hit.get_raw_time();
    rawt -= el.get_laser_rawt();
    bx_time_calib->Fill (ch, rawt);
  }
  return ev;
}

void bx_calib_laben_time_align::end () {
  if (!b_has_data) {
    get_message(bx_message::error) << "No laser data are present in Run: time calibration not performed." << dispatch;
    // deallocate the vectors
    delete [] f_time_offset;
    delete [] f_time_sigma;
    return;
  }  
  db_run& run_info = bx_dbi::get()->get_run ();

  // mean value of the laser peak
  TH1D *mean_proj = bx_time_calib->ProjectionY("mean_proj");
  int32_t max_bin = mean_proj->GetMaximumBin();
  float offset = mean_proj->GetBinCenter(max_bin);
  int32_t entries_run = (int32_t) mean_proj->Integral();
  float mean_entries = (float) entries_run / (constants::laben::channels);

  //number of events in the main and reflected peak
  int32_t peak_entries_run = (int32_t) mean_proj->Integral(max_bin - 20, max_bin + 20);

  mean_proj->Delete();
  get_message(bx_message::log) << "Mean value of the Laser Peak: " << offset << " ns" << dispatch;
  get_message(bx_message::log) << "Mean entries on each Channel: " << mean_entries << dispatch;
  get_message(bx_message::log) << "Peak contains about " <<  (double) peak_entries_run/ (double) entries_run * 100. << "% of all hits" << dispatch;
 
  // fit of the single channel histos
  TF1* g0 = new TF1("g0","gaus",500.,1000.);
  Double_t tpar[3];
  int32_t not_connected = 0;
  int32_t no_light = 0;
  int32_t not_conv = 0;
  int32_t large_drift = 0;

  // main loop
  for (int32_t i=0; i<constants::laben::channels; i++) {
   
    // discard non Pmt channels and disconnected PMTs
    const db_channel_laben &ch_info = dynamic_cast <const db_channel_laben&> (bx_dbi::get ()->get_channel (i + 1));
    if (!ch_info.is_ordinary() || ch_info.is_pmt_disconnected()){
      get_message(bx_message::log) << "Channel " << i+1 << " empty or disconnected" << dispatch;
      not_connected++;
      continue;
    }

    TH1D *proj = bx_time_calib->ProjectionY("proj", i+1, i+1);
    int32_t entries_channel = (int32_t) proj->Integral();
    int32_t max_bin_channel = proj->GetMaximumBin();
    int32_t max_used = max_bin_channel;
    if (max_bin_channel < 240 || max_bin_channel > 340) max_used = max_bin;
    int32_t peak_entries_channel = (int32_t) proj->Integral(max_used - 20, max_used + 20);

    // discard channels with no light 
    if (entries_channel == 0 ){
      proj->Delete();
      get_message(bx_message::log) << "Pmt " << i+1 << " has no entries " <<  dispatch;
      no_light++;  
        //check if the previous values were ok and if yes set the current vector top these values 
        //may be in this run just the precalibration failed!
      if(fabs(run_info.get_laben_time_offset (i+1)) > 0 && run_info.get_laben_time_sigma(i+1) > 0){
	f_time_offset[i] = run_info.get_laben_time_offset (i+1);
	f_time_sigma[i] = run_info.get_laben_time_sigma (i+1);
	get_message(bx_message::log) << "Pmt " << i+1 << " used values from the previous run, offset  " << f_time_offset[i] << ", sigma " << 	f_time_sigma[i] <<  dispatch;
      }
      continue;
    }

    //discard PMTs with some entries but very few in the main peak
    if (peak_entries_channel < 100 || (double) peak_entries_channel / (double) entries_channel < 0.10) {
      get_message(bx_message::log) << "Pmt " << i+1 << ": low statistics, total entries = " << entries_channel << ", entries in the peak " << peak_entries_channel << " in the region " << proj->GetBinCenter(max_used - 20) << "-" << proj->GetBinCenter(max_used + 20) <<  " ns " << dispatch;
      proj->Delete();
      no_light++; 
      //check if the previous values were ok and if yes set the current vector top these values 
      if(fabs(run_info.get_laben_time_offset (i+1)) > 0 && run_info.get_laben_time_sigma(i+1) > 0){
	f_time_offset[i] = run_info.get_laben_time_offset (i+1);
	f_time_sigma[i] = run_info.get_laben_time_sigma (i+1);
	get_message(bx_message::log) << "Pmt " << i+1 << " used values from the previous run, offset  " << f_time_offset[i] << ", sigma " << 	f_time_sigma[i] <<  dispatch;
      }
      continue;
    }
    
    tpar[0] = peak_entries_channel;
    tpar[1] = proj->GetBinCenter(max_used);
    tpar[2] = 3.;
    g0->SetParameters(&tpar[0]);
    g0->SetRange(tpar[1] - 12., tpar[1] + 12.);
    int32_t FitResult = proj->Fit("g0","QRL0");
    g0->GetParameters(tpar);
    proj->Delete();
    
    // discard channels with bad fit
    if (FitResult != 0 || fabs (tpar[1] - offset) > 20 || tpar[2] <= 0 || tpar[2] > 5 ){
      get_message(bx_message::log) << "Pmt " << i+1 << ": bad fit; entries = " 
	<< entries_channel << "; peak = " << (int32_t) (tpar[1] - offset) << "; sigma = " << (int32_t) tpar[2] << dispatch;
      not_conv++;  
      //check if the previous values were ok and if yes set the current vector top these values 
      if(fabs(run_info.get_laben_time_offset (i+1)) > 0 && run_info.get_laben_time_sigma(i+1) > 0){
	f_time_offset[i] = run_info.get_laben_time_offset (i+1);
	f_time_sigma[i] = run_info.get_laben_time_sigma (i+1);
	get_message(bx_message::log) << "Pmt " << i+1 << " used values from the previous run, offset  " << f_time_offset[i] << ", sigma " << 	f_time_sigma[i] <<  dispatch;
      }
      continue;
    }

    // warn for channels with large peak drift compared to previous run
    if ( fabs (run_info.get_laben_time_offset (i+1) - (tpar[1] - offset)) > 10 && run_info.get_laben_time_sigma (i+1) != 0 ) {
      get_message(bx_message::log)<<"Pmt "<<i+1<<": large peak drift; new = " << tpar[1] - offset <<
	", old = " << run_info.get_laben_time_offset(i+1) << "; sigma = " << tpar[2] << dispatch;
      large_drift++;
      // if requested upon configuration, discard these data
      if (get_parameter ("large_drift").get_bool ()) continue;
    }
    
    //fill the vectors which will set the visitors
    f_time_offset[i] = (float) tpar[1] - offset;
    f_time_sigma[i]  = (float) tpar[2];
    get_message(bx_message::log) << "Pmt " << i+1 << ": well calibrated; entries = " << entries_channel
      << " ; peak = " << (int32_t) (tpar[1] - offset) << " ; sigma = " << (int32_t) tpar[2] << dispatch;
  }
  g0->Delete();

  int32_t good_ch = constants::laben::channels - not_connected - no_light - not_conv;
  
  get_message(bx_message::info) << not_connected << " channels not connected to any Pmt" << dispatch;
  get_message(bx_message::info) << no_light << " channels with no direct light" << dispatch;
  get_message(bx_message::info) << not_conv << " channels with bad calibration fit" << dispatch;
  get_message(bx_message::info) << large_drift << " channels have large peak drift" << dispatch;
  get_message(bx_message::info) << good_ch << " channels correctly calibrated" << dispatch;
 
  //fill the visitors to db_run
  if ((float)good_ch/(constants::laben::channels-not_connected)<=0.7 && get_parameter ("write_calib").get_bool ()) {
    get_message(bx_message::warn) << "Too may channels have bad time calibration: DB tables will not be updated!" << dispatch;
  }
  else if (get_parameter ("write_calib").get_bool ()){//set visitors
    for (int32_t i=0; i<constants::laben::channels; i++) {
      run_info.set_laben_time_offset (i + 1, f_time_offset[i], this);
      run_info.set_laben_time_sigma (i + 1, f_time_sigma[i], this);
      get_message(bx_message::log) << "Values in DB: lG: "<< i + 1 << " peak " <<  f_time_offset[i] << " rms " <<  f_time_sigma[i]  << dispatch;
    }
    run_info.write_laben_laser_time (true, this);//wtite to DB
  }

  get_message(bx_message::debug) << "end" << dispatch;
  // deallocate the vectors
  delete [] f_time_offset;
  delete [] f_time_sigma;
}

/*
 * $Log: bx_calib_laben_time_align.cc,v $
 * Revision 1.12  2008/07/01 13:48:38  ludhova
 * check of write_calib before writing to DB
 *
 * Revision 1.11  2007-12-17 15:55:38  ludhova
 * cut on absolute number of hits in direct peak added
 *
 * Revision 1.10  2007-11-23 08:43:01  ludhova
 * small bug in the check of the previous values
 *
 * Revision 1.9  2007-11-22 14:52:34  ludhova
 * changing fit region
 *
 * Revision 1.8  2007-11-22 14:03:57  ludhova
 * setting previous values if the fit fails
 *
 * Revision 1.7  2007-11-21 17:26:08  ludhova
 * new logics in setting visitors
 *
 * Revision 1.6  2007-11-21 13:33:42  ludhova
 * disconnected pmts skipped
 *
 * Revision 1.5  2006-11-17 17:40:08  razeto
 * Never return 0
 *
 * Revision 1.4  2006/08/21 11:20:08  razeto
 * Updated to new barn_interface
 *
 * Revision 1.3  2006/01/02 21:27:18  razeto
 * Changed test target to file target for db_barn (missing maintainer)
 *
 * Revision 1.2  2005/05/05 17:13:40  monzani
 * S/Getters names updated according to the interface modification. Debugging.
 *
 * Revision 1.1  2005/05/05 10:02:22  monzani
 * Laser time calibrations added.
 *
 */
