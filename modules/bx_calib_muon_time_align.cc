/* BOREXINO Reconstruction program
 *
 * Author: Davide D'Angelo <davide.dangelo@mi.infn.it>, Michael Wurm <mwurm@ph.tum.de>, Quirin Meindl <qmeindl@ph.tum.de>
 * Based on code from Maria Elena Monzani
 * Maintainer: Michael Wurm <mwurm@ph.tum.de>
 *
 * $Id: bx_calib_muon_time_align.cc,v 1.7 2008/10/16 19:48:16 ddangelo Exp $
 *
 * Implementation of bx_calib_muon_time_align *
 * 
 */
#include "bx_calib_muon_time_align.hh"
#include "messenger.hh"
#include "bx_dbi.hh"
#include "db_profile.hh"
#include "db_run.hh"
#include "barn_interface.hh"
#include "TF1.h"
#include <cmath>


inline bool bigger (float a, float b, float min){
return a > b + min;
}

inline bool smaller (float a, float b, float min){
return a < b - min;
}

inline bool bigger_error (float a, float b, float c){
float err = :: sqrt(a+b);
if (err < 2) err = 2;
err = err * c;
return bigger (a, b , err);
}

inline bool smaller_error (float a, float b, float c){
float err = :: sqrt(a+b);
if (err < 2) err = 2;
err = err * c;
return smaller (a, b , err);
}




// ctor
bx_calib_muon_time_align::bx_calib_muon_time_align (): bx_base_module("bx_calib_muon_time_align", bx_base_module::main_loop) {
  require_event_stage (bx_detector::muon, bx_base_event::decoded);
  require_trigger_type (bx_trigger_event::laser394);
}

void bx_calib_muon_time_align::begin () {
  get_message(bx_message::debug) << "begin" << dispatch;

  // creation of the 2-dim histogram 
  int32_t nch = constants::muon::channels;
  muon_time_calib   = new TH2S ("muon_time_calib"  , "Time Calibration Histogram"  ,nch,0,nch, 1200,500,1100);
  barn_interface::get ()->store (barn_interface::file, muon_time_calib  , this);

  b_has_data = false;
  i4_nevents = 0;

  if ( get_parameter ("write_calib_time").get_bool () && bx_dbi::get()->get_run ().is_muon_laser_present () ) 
    get_message (bx_message::warn) << "Laser calibrations already present for this run." << dispatch; 
}

bx_echidna_event* bx_calib_muon_time_align::doit (bx_echidna_event *ev) {
  const bx_muon_event& er = ev->get_muon();

  int32_t nhits = er.get_decoded_nhits();

  // Update has data field
  if (!nhits) return ev;
  i4_nevents++;
  b_has_data |= true;

  for(int32_t i = 0; i < nhits; i++){
    const bx_muon_decoded_hit& hit = er.get_decoded_hit(i);
    uint16_t mch = hit.get_raw_hit().get_muon_channel();
    float time   = hit.get_time  ();
    muon_time_calib  ->Fill(mch,time);
  }

  return ev;
}

void bx_calib_muon_time_align::end () {
  if (!b_has_data) {
    get_message(bx_message::error) << "No laser data are present in Run: time calibration not performed." << dispatch;
    return;
  }  

  // retrieve parameters
  float minimum_efficiency = get_parameter ("minimum_efficiency").get_float () ;
  float time_allowance     = get_parameter ("time_allowance"    ).get_float () ;
  float max_time_offset    = get_parameter ("max_time_offset"   ).get_float () ;
  float min_time_rms       = get_parameter ("min_time_rms"      ).get_float () ;
  float max_time_rms       = get_parameter ("max_time_rms"      ).get_float () ;
  float time_drift         = get_parameter ("time_drift"        ).get_float () ;

  const db_profile& profile_info = bx_dbi::get()->get_profile ();
            db_run&     run_info = bx_dbi::get()->get_run     ();

  // mean value of the laser peak
  TH1D *mean_proj_time   = muon_time_calib  ->ProjectionY("mean_proj_time"  );
  int32_t avg_entries = int32_t(mean_proj_time  ->Integral()/208.);
  float avg_offset  = mean_proj_time  ->GetBinCenter(mean_proj_time  ->GetMaximumBin());

  mean_proj_time  ->Delete();
  get_message(bx_message::log) << "Average number of entries: " << avg_entries             << dispatch;
  get_message(bx_message::log) << "Average time offset: "       << avg_offset  << " ns"    << dispatch;

  // creation of the result vectors
  float *time_offset   = new float[constants::muon::channels]; 
  float *time_sigma    = new float[constants::muon::channels]; 
  std::fill_n (time_offset ,constants::muon::channels,       0.);
  std::fill_n (time_sigma  ,constants::muon::channels,       0.);

  // some ctrs
  int32_t not_connected = 0, no_light = 0;
  int32_t not_conv_time = 0, no_mean_time = 0;
  int32_t large_drift_time = 0;
  int32_t good_in_time = 0;

  // main loop
  for (int32_t i=0; i<constants::muon::channels; i++) {
    // discard non Pmt channels
    if (profile_info.logical_channel_description (i+constants::muon::channel_offset+1) != db_profile::ordinary){
      get_message(bx_message::log) << "Mch " << i << " is service channel." << dispatch;
      not_connected++;
      continue;
    }
    TH1D *proj_time = muon_time_calib->ProjectionY("proj_time", i+1, i+1);
    int32_t entries = int32_t(proj_time->Integral());
    // discard channels with no direct light
    if (float(entries)/float(i4_nevents) < minimum_efficiency) {
      get_message(bx_message::log) << "Mch " << i << ": low statistics; entries = " << entries << dispatch;
      proj_time->Delete();
      no_light++;  
      continue;
    }

    // *** TIME CALIBRATION ***
    TF1* func_time = new TF1("func_time","gaus");
    Double_t tpar[3];
    tpar[0] = proj_time->GetBinContent(proj_time->GetMaximumBin());
    tpar[1] = proj_time->GetBinCenter (proj_time->GetMaximumBin());
    tpar[2] = max_time_rms/2.;
    func_time->SetParameters(tpar);
    func_time->SetRange(tpar[1]-time_allowance,tpar[1]+time_allowance);
    proj_time->Fit("func_time","QRL0"); // Q=quit; R=use_function_range; L=loglikelihood ; 0=do_not_draw;
    func_time->GetParameters(tpar);
    delete func_time;
    // discard channels with bad fit
    if ( fabs (tpar[1]-avg_offset) > max_time_offset || tpar[2] <= min_time_rms || tpar[2] > max_time_rms ){

      get_message(bx_message::log) << "Mch " << i << ": time : entries = " << entries
                                    << "; fitted peak = " << tpar[1] - avg_offset << "; fitted sigma = " << tpar[2] << dispatch;

      if (fabs (tpar[1]-avg_offset) > max_time_offset)
        get_message(bx_message::log) << "Mch " << i << ": time : bad fit due to bad peak position." << dispatch;

      if (tpar[2] <= min_time_rms || tpar[2] > max_time_rms)
        get_message(bx_message::log) << "Mch " << i << ": time : bad fit due to bad sigma." << dispatch;

      not_conv_time++;  
      tpar[1] = proj_time->GetMean();
      tpar[2] = proj_time->GetRMS();

      if ( fabs (tpar[1]-avg_offset) > max_time_offset || tpar[2] <= min_time_rms || tpar[2] > max_time_rms ){
 
        if (fabs (tpar[1]-avg_offset) > max_time_offset)
          get_message(bx_message::log) << "Mch " << i << ": time : bad mean due to bad peak position." << dispatch;

        if (tpar[2] <= min_time_rms || tpar[2] > max_time_rms)
          get_message(bx_message::log) << "Mch " << i << ": time : bad rms due to bad sigma." << dispatch;

	no_mean_time++;
	
	tpar[1] = avg_offset;
	tpar[2] = 0.;
	continue;
      }
    }
    proj_time->Delete();
    // warn for channels with large peak drift compared to previous run
    if ( fabs (run_info.get_muon_time_offset (i+constants::muon::channel_offset+1) - (tpar[1] - avg_offset)) > time_drift && run_info.get_muon_time_sigma (i+constants::muon::channel_offset+1) != 0 ) {
      get_message(bx_message::log)<<"Mch "<<i<<": time : large peak drift; new peak = " << tpar[1] - avg_offset << ", new sigma = " << tpar[2] << "; old peak = " << run_info.get_muon_time_offset(i+constants::muon::channel_offset+1) << dispatch;
      large_drift_time++;
      // if requested upon configuration, discard these data
      if (get_parameter ("large_drift").get_bool ()) continue;
    }
    time_offset[i] = (float) tpar[1] - avg_offset;
    time_sigma [i] = (float) tpar[2];
    good_in_time++;

    get_message(bx_message::log) << "Mch " << i << ": time : calibrated; entries = " << int32_t(entries) 
				 << " ; peak = " << /*(int32_t) */(tpar[1] - avg_offset) << " ; sigma = " << /*(int32_t)*/ tpar[2] << dispatch;
  }

  get_message(bx_message::info) << not_connected      << " channels not connected to a Pmt"  << dispatch;
  get_message(bx_message::info) << no_light           << " channels with no direct light"    << dispatch;
  get_message(bx_message::info) << not_conv_time      << " channels with bad time fit"       << dispatch;
  get_message(bx_message::info) << no_mean_time       << " channels with bad mean time"      << dispatch;
  get_message(bx_message::info) << large_drift_time   << " channels have large time drift"   << dispatch;
  get_message(bx_message::info) << good_in_time       << " channels with time calibration"   << dispatch;

  // fill the visitors to db_run
  // forward mean time peak to calib_muon_charge_peak via mch 0
  run_info.set_muon_time_offset ( 3001, avg_offset, this );
  for (int32_t i=1; i<constants::muon::channels; i++) {
    int32_t lg = i + 1 + constants::muon::channel_offset;
    if (profile_info.logical_channel_description (lg) != db_profile::ordinary) continue;
    run_info.set_muon_time_offset  (lg , time_offset [i], this);
    run_info.set_muon_time_sigma   (lg , time_sigma  [i], this);
  }

  // deallocate the vectors
  delete [] time_offset;
  delete [] time_sigma;

  // decide whether to update DB tables
  if ( get_parameter ("write_calib_time").get_bool () ) {
    if ( float(good_in_time) / (constants::muon::channels-not_connected) > get_parameter ("update_threshold").get_float() )
      run_info.write_muon_laser_time   (true, this);
    else
      get_message(bx_message::warn) << "Too may channels have bad time calibration: DB time table will not be updated." << dispatch;
  }

  get_message(bx_message::debug) << "end" << dispatch;
}

/*
 * $Log: bx_calib_muon_time_align.cc,v $
 * Revision 1.7  2008/10/16 19:48:16  ddangelo
 * channel efficiency check improved,
 * time cut
 * some parameters tuning
 *
 * Revision 1.6  2008-10-14 15:25:58  wurm
 * debugging
 *
 * Revision 1.5  2008-10-14 14:33:49  wurm
 * split module off calib_muon_channel
 *
 * Revision 1.6  2008-10-07 14:42:54  wurm
 * omit charge calibration if time calibration files
 *
 * Revision 1.5  2008-10-07 14:16:27  wurm
 * modified time fit
 *
 * Revision 1.4  2008-10-01 09:44:23  wurm
 * adjusted histogram range of charge
 *
 * Revision 1.3  2008-09-29 13:32:47  wurm
 * rms for bad charge fits, further minor adjustments
 *
 * Revision 1.2  2006-09-11 20:56:40  ddangelo
 * added drift check, params readout, charge fit.
 * First params tuning, first run tests.
 * More to be done
 *
 * Revision 1.1  2006/09/11 14:12:38  ddangelo
 * bx_calib_muon_charge_peak and bx_calib_muon_time_align
 * replaced by
 * bx_calib_muon_channel
 *
 *
 */
