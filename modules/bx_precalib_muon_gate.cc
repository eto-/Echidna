/* BOREXINO Reconstruction program
 *
 * Author: Davide D'Angelo <davide.dangelo@mi.infn.it>
 * Maintainer: Davide D'Angelo <davide.dangelo@mi.infn.it>
 *
 * $Id: bx_precalib_muon_gate.cc,v 1.4 2016/07/29 00:49:43 ddangelo Exp $
 *
 * Implementation of bx_precalib_muon_gate.hh
 *
 */
#include <stdlib.h>
#include "bx_precalib_muon_gate.hh"
#include "messenger.hh"
#include "constants.hh"
#include "bx_dbi.hh"
#include "db_run.hh"
#include "db_profile.hh"
#include "bx_detector.hh"
#include "barn_interface.hh"

// ctor
bx_precalib_muon_gate::bx_precalib_muon_gate (): bx_base_module("bx_precalib_muon_gate", bx_base_module::precalib_cycle2) {
  require_event_stage (bx_detector::muon, bx_base_event::raw);
  require_trigger_type (bx_trigger_event::pulser);
}

void bx_precalib_muon_gate::begin () {
  // retrieve parameters
  int run_number =bx_dbi::get ()->get_run ().get_number();
  bool is_new_tdc = run_number > 17500 ? true : false;
  bool is_new_trg = run_number > 26579 ? true : false;
  min_evnum      = (unsigned short)(get_parameter("min_evnum"     ).get_int());
  max_evnum      = (unsigned short)(get_parameter("max_evnum"     ).get_int());
  min_gate_width = (unsigned short)(get_parameter(is_new_tdc ? "min_new_gate_width" : "min_gate_width").get_int());
  max_gate_width = (unsigned short)(get_parameter(is_new_tdc ? "max_new_gate_width" : "max_gate_width").get_int());
  precalib_pulse_tolerance = (unsigned short)(get_parameter(is_new_trg ? "pulse_tolerance_new_trg" : "pulse_tolerance_old_trg").get_int());

  if (get_parameter("use_fixed_time").get_bool())
    pulse_time = (unsigned short)(get_parameter("fixed_pulse_time").get_int());
  else { 
    pulse_time =  bx_dbi::get ()->get_run ().get_muon_precalib_pulse_time();
    if (!pulse_time) {
      get_message(bx_message::warn) << "No precalib pulse time from DB. Using fixed time." << dispatch;
      pulse_time = (unsigned short)(get_parameter("fixed_pulse_time").get_int());
    }
  }
  get_message(bx_message::info) << "Pulse time " << pulse_time << " pulse tolerance " << precalib_pulse_tolerance << dispatch;

  gates = new TH2S("muon_gate_entries", "Muon Gate Entries", constants::muon::channels, 0, constants::muon::channels, 
		       10*(max_gate_width-min_gate_width), min_gate_width, max_gate_width);
  barn_interface::get ()->store (barn_interface::file, gates, this);

  means = new TH1F("muon_gate_means", "Muon Gate Means", max_gate_width-min_gate_width, min_gate_width, max_gate_width);
  sigmas = new TH1F("muon_gate_sigmas", "Muon Gate Sigmas", 80, 0, 8);
  barn_interface::get ()->store (barn_interface::file, means, this);
  barn_interface::get ()->store (barn_interface::file, sigmas, this);

  b_has_data = false;

  get_message(bx_message::debug) << "begin" << dispatch;
}

bx_echidna_event* bx_precalib_muon_gate::doit (bx_echidna_event *ev) {
  if (ev->get_event_number() < min_evnum || ev->get_event_number() > max_evnum) return ev;
  const bx_muon_event& er = ev->get_muon();
  if (!er.get_raw_nhits()) return ev;

  b_has_data = true;

  const db_profile& profile_info = bx_dbi::get ()->get_profile();

  for (int i = 0; i <er.get_raw_nhits(); i++) { 
    const bx_muon_raw_hit& h = er.get_raw_hit(i);

    // skip special channels
    if (profile_info.logical_channel_description(h.get_logical_channel()) != db_profile::ordinary)
      continue;

    // skip hits unrelated to calibration pulse (time filter)
    if ( ::abs( h.get_lead_time() - pulse_time ) > precalib_pulse_tolerance )
      continue;

    // check if width is ok
    if ( h.get_time_diff() > min_gate_width && h.get_time_diff() <= max_gate_width )
      gates->Fill(h.get_muon_channel(), h.get_time_diff());
    else
      get_message(bx_message::warn) << "Unusual gate width (skipping); ev " << ev->get_event_number() 
				    << " mch " << h.get_muon_channel() << " width " << h.get_time_diff() << dispatch; 
  } // end ov loop over hits
  return ev;
}

void bx_precalib_muon_gate::end () {
  if(!b_has_data) return;

//  db_run& run_info = bx_dbi::get ()->get_run ();
  const db_profile& profile_info = bx_dbi::get ()->get_profile();

  for(int i=0; i<constants::muon::channels; i++) { 
    int lg = constants::muon::channel_offset + i + 1;
    if (profile_info.logical_channel_description(lg) == db_profile::ordinary) {
      TH1D* proj = gates->ProjectionY("proj_gate", i+1, i+1);
      double entries = proj->Integral();

      if (entries) {
	double mode = proj->GetBinCenter(proj->GetMaximumBin());
	double mean = proj->GetMean();
	double rms  = proj->GetRMS();
	means->Fill(mean);
	sigmas->Fill(rms);

//	run_info.set_muon_precalib_gate(i+constants::muon::channel_offset+1, mean, this);	
//	run_info.set_muon_precalib_pedsigma(i+constants::muon::channel_offset+1, rms, this);	
	get_message(bx_message::info) << "mch " << i << " entries " << entries << " mode " << mode << " mean " << mean << " rms " << rms << dispatch;
      }
      else {
//	run_info.set_muon_precalib_gate(i+constants::muon::channel_offset+1, 0., this);	
//	run_info.set_muon_precalib_pedsigma(i+constants::muon::channel_offset+1, 0., this);	
	get_message(bx_message::warn) << "mch " << i << " gate not available" << dispatch;
      }	
    }
    else 
      ;//      get_message(bx_message::info) << "mch " << i << " is not ordinary channel" << dispatch;

  } /* end of loop on channels */

//  if (get_parameter ("write_precalib").get_bool ()) bx_dbi::get ()->get_run ().write_muon_precalib (true, this);

  get_message(bx_message::debug) << "end" << dispatch;
}  

/*
 * $Log: bx_precalib_muon_gate.cc,v $
 * Revision 1.4  2016/07/29 00:49:43  ddangelo
 * dynamic handling of muon precalib pulse
 *
 * Revision 1.3  2012/01/16 12:21:00  ddangelo
 * modified to handle different tdc resolutions
 *
 * Revision 1.2  2008-10-26 12:04:36  ddangelo
 * debug
 *
 * Revision 1.1  2008-10-25 09:59:54  ddangelo
 * added. cloned from pedestal module
 *
 *
 */
