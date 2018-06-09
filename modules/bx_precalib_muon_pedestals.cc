/* BOREXINO Reconstruction program
 *
 * Author: Davide D'Angelo <davide.dangelo@lngs.infn.it>
 * Maintainer: Davide D'Angelo <davide.dangelo@lngs.infn.it>
 *
 * $Id: bx_precalib_muon_pedestals.cc,v 1.23 2016/07/29 00:49:43 ddangelo Exp $
 *
 * Implementation of bx_precalib_muon_pedestals.hh
 *
 */
#include <stdlib.h>
#include "bx_precalib_muon_pedestals.hh"
#include "messenger.hh"
#include "constants.hh"
#include "bx_dbi.hh"
#include "db_run.hh"
#include "db_profile.hh"
#include "bx_detector.hh"
#include "barn_interface.hh"

// ctor
bx_precalib_muon_pedestals::bx_precalib_muon_pedestals (): bx_base_module("bx_precalib_muon_pedestals", bx_base_module::precalib_cycle2) {
  require_event_stage (bx_detector::muon, bx_base_event::raw);
  require_trigger_type (bx_trigger_event::pulser);
}

void bx_precalib_muon_pedestals::begin () {
  // retrieve parameters
  int run_number =bx_dbi::get ()->get_run ().get_number();
  bool is_new_tdc = run_number > 17500 ? true : false;
  bool is_new_trg = run_number > 26579 ? true : false;
  min_evnum      = (unsigned short)(get_parameter("min_evnum"     ).get_int());
  max_evnum      = (unsigned short)(get_parameter("max_evnum"     ).get_int());
  min_pedestal_width = (unsigned short)(get_parameter(is_new_tdc ? "min_new_pedestal_width" : "min_pedestal_width").get_int());
  max_pedestal_width = (unsigned short)(get_parameter(is_new_tdc ? "max_new_pedestal_width" : "max_pedestal_width").get_int());
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

  pedestals = new TH2S("muon_pedestals_entries", "Muon Pedestal Entries", constants::muon::channels, 0, constants::muon::channels, 
		       10*(max_pedestal_width-min_pedestal_width), min_pedestal_width, max_pedestal_width);
  barn_interface::get ()->store (barn_interface::file, pedestals, this);

  means = new TH1F("muon_pedestal_means", "Muon Pedestal Means", max_pedestal_width-min_pedestal_width, min_pedestal_width, max_pedestal_width);
  sigmas = new TH1F("muon_pedestal_sigmas", "Muon Pedestal Sigmas", 80, 0, 8);
  barn_interface::get ()->store (barn_interface::file, means, this);
  barn_interface::get ()->store (barn_interface::file, sigmas, this);

  b_has_data = false;

  get_message(bx_message::debug) << "begin" << dispatch;
}

bx_echidna_event* bx_precalib_muon_pedestals::doit (bx_echidna_event *ev) {
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
    if ( h.get_time_diff() > min_pedestal_width && h.get_time_diff() <= max_pedestal_width )
      pedestals->Fill(h.get_muon_channel(), h.get_time_diff());
    else
      get_message(bx_message::warn) << "Unusual pedestal width (skipping); ev " << ev->get_event_number() 
				    << " mch " << h.get_muon_channel() << " width " << h.get_time_diff() << dispatch; 
  } // end ov loop over hits
  return ev;
}

void bx_precalib_muon_pedestals::end () {
  if(!b_has_data) return;

  db_run& run_info = bx_dbi::get ()->get_run ();
  const db_profile& profile_info = bx_dbi::get ()->get_profile();

  for(int i=0; i<constants::muon::channels; i++) { 
    int lg = constants::muon::channel_offset + i + 1;
    if (profile_info.logical_channel_description(lg) == db_profile::ordinary) {
      TH1D* proj = pedestals->ProjectionY("proj_ped", i+1, i+1);
      double entries = proj->Integral();

      if (entries) {
	double mode = proj->GetBinCenter(proj->GetMaximumBin());
	double mean = proj->GetMean();
	double rms  = proj->GetRMS();
	means->Fill(mean);
	sigmas->Fill(rms);

	run_info.set_muon_precalib_pedestal(i+constants::muon::channel_offset+1, mean, this);	
	run_info.set_muon_precalib_pedsigma(i+constants::muon::channel_offset+1, rms, this);	
	get_message(bx_message::info) << "mch " << i << " entries " << entries << " mode " << mode << " mean " << mean << " rms " << rms << dispatch;
      }
      else {
	run_info.set_muon_precalib_pedestal(i+constants::muon::channel_offset+1, 0., this);	
	run_info.set_muon_precalib_pedsigma(i+constants::muon::channel_offset+1, 0., this);	
	get_message(bx_message::warn) << "mch " << i << " pedestal not available" << dispatch;
      }	
    }
    else 
      ;//      get_message(bx_message::info) << "mch " << i << " is not ordinary channel" << dispatch;

  } /* end of loop on channels */

  if (get_parameter ("write_precalib").get_bool ()) bx_dbi::get ()->get_run ().write_muon_precalib (true, this);

  get_message(bx_message::debug) << "end" << dispatch;
}  

/*
 * $Log: bx_precalib_muon_pedestals.cc,v $
 * Revision 1.23  2016/07/29 00:49:43  ddangelo
 * dynamic handling of muon precalib pulse
 *
 * Revision 1.22  2012/01/16 12:21:00  ddangelo
 * modified to handle different tdc resolutions
 *
 * Revision 1.21  2008-10-26 12:04:36  ddangelo
 * debug
 *
 * Revision 1.20  2008-10-25 09:59:32  ddangelo
 * added check on event number
 *
 * Revision 1.19  2008-06-20 16:21:36  razeto
 * Added an include to compile with gcc 4.3
 *
 * Revision 1.18  2007-02-23 14:47:04  ddangelo
 * fixed data types
 *
 * Revision 1.17  2006/09/05 12:06:03  ddangelo
 * cleaned up
 * histograms filled
 *
 * Revision 1.16  2006/08/21 15:37:47  ddangelo
 * added a few test histograms (new barn interface).
 * introduced the possibility to use a fixed pulse time (parameter switch and constant setting).
 * Retrieving of find pulse results moved in begin().
 *
 * Revision 1.15  2006/01/25 13:12:00  misiaszek
 * Added a missing include (auth from maintainer)
 *
 * Revision 1.14  2005/03/02 16:02:06  ddangelo
 * included trg type check to follow framework changes
 *
 * Revision 1.13  2004/11/27 16:42:49  ddangelo
 * sigma writing to db_run added.
 * command to write to db added.
 * Adapted to use logical channel instead of muon channel.
 *
 * Revision 1.12  2004/11/26 15:25:11  razeto
 * Thanks to Davide who so kindly underlined a syntactical error, I fixed the "Maintainer" word
 *
 * Revision 1.11  2004/11/26 14:06:23  razeto
 * Added Mantainer field
 *
 * Revision 1.10  2004/09/22 13:29:16  ddangelo
 * updated bx_reco_event to bx_echidna_event
 *
 * Revision 1.9  2004/09/22 11:29:58  ddangelo
 * added require_event_stage in module's constructor
 *
 * Revision 1.8  2004/06/01 11:36:27  ddangelo
 * updated to new constants syntax
 *
 * Revision 1.7  2004/05/20 15:45:56  ddangelo
 *  check on has_data introduced
 *  check on chann descr implemented with db data (db_profile for the moment)
 *  instead of static event members.
 *  computetion done with algotithms instead of c like math.
 *  Data stored in a std::vector* (histo use).
 *
 * Revision 1.5  2004/04/28 11:04:49  ddangelo
 * Modules adapted to call the right const/nonconst event portion getter
 *
 * Revision 1.4  2004/04/28 10:48:21  ddangelo
 * Pedestal calculation reworked.
 *
 * Revision 1.3  2004/04/27 09:46:50  ddangelo
 * modifications to match new event structure.
 * Other debugging
 *
 * Revision 1.1  2004/04/20 11:36:25  ddangelo
 * added 2 modules for outer muon precalibration
 *
 *
 */
