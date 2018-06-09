/* BOREXINO Reconstruction program
 *
 * Author: Davide D'Angelo <davide.dangelo@lngs.infn.it>
 * Maintainer: Davide D'Angelo <davide.dangelo@lngs.infn.it>
 *
 * $Id: bx_precalib_muon_findpulse.cc,v 1.21 2016/07/29 00:49:43 ddangelo Exp $
 *
 * Implementation of bx_precalib_muon_findpulse.hh
 *
 */

#include "bx_precalib_muon_findpulse.hh"
#include "messenger.hh"
#include "constants.hh"
#include "bx_dbi.hh"
#include "db_run.hh"
#include "db_profile.hh"
#include "bx_detector.hh"
#include "barn_interface.hh"

bx_precalib_muon_findpulse::bx_precalib_muon_findpulse(): bx_base_module("bx_precalib_muon_findpulse", bx_base_module::precalib_cycle1) {
  require_event_stage (bx_detector::muon, bx_base_event::raw);
  require_trigger_type (bx_trigger_event::pulser);
}

void bx_precalib_muon_findpulse::begin () {
  const db_profile& profile_info = bx_dbi::get ()->get_profile ();

  pulses = new TH1F("muon_precalib_pulse_time", "Muon precalib pulse time", profile_info.muon_tdc_range()/16 , 0, profile_info.muon_tdc_range());
  barn_interface::get ()->store (barn_interface::file, pulses, this);

  b_has_data = false;
  get_message(bx_message::debug) << "begin" << dispatch;
}

bx_echidna_event* bx_precalib_muon_findpulse::doit (bx_echidna_event *ev) {

  const bx_muon_event& er = ev->get_muon();

  b_has_data |= bool(er.get_raw_nhits());

  const db_profile& profile_info = bx_dbi::get ()->get_profile ();

  // start looping over hits
  for (int i = 0; i <er.get_raw_nhits(); i++) { 
    const bx_muon_raw_hit& h = er.get_raw_hit(i);

    // skip special channels
    if (profile_info.logical_channel_description(h.get_logical_channel()) != db_profile::ordinary)
      continue;

    pulses->Fill(h.get_lead_time());

  } // end of loop over hits

  return ev; 
}

void bx_precalib_muon_findpulse::end () {
  if (!b_has_data) return; 

  db_run& run_info = bx_dbi::get ()->get_run ();

  double entries = pulses->GetEntries();
  if (entries < 50000 || entries > 600000) // We expect 1 or 2 entries per channel x 1000 events 
	get_message(bx_message::error) << "entries " << entries << " is out of range (50k-600k)." << dispatch;

 
  TSpectrum *s = new TSpectrum();
  int n_peaks = s->Search(pulses);
  float *positions = s->GetPositionX();
  if (n_peaks < 1 || n_peaks > 2) 
	get_message(bx_message::error) << "npeaks " << entries << " is out of range (1-2)." << dispatch;
  int precalib_peak_index = 0;
  if (n_peaks == 2 && positions[1] < positions[0]) 
	precalib_peak_index = 1;

  get_message(bx_message::info) << "entries " << entries << " n_peaks " << n_peaks << " pos0 " << positions[0] << " pos1 " << positions[1] << " selecting peak " << precalib_peak_index << dispatch;

  unsigned long value_for_db = positions[precalib_peak_index];
  run_info.set_muon_precalib_pulse_time( value_for_db, this);	

  get_message(bx_message::debug) << "end" << dispatch;
}

/*
 * $Log: bx_precalib_muon_findpulse.cc,v $
 * Revision 1.21  2016/07/29 00:49:43  ddangelo
 * dynamic handling of muon precalib pulse
 *
 * Revision 1.20  2007/12/07 17:33:59  ddangelo
 * compliant with new dynamically loaded tdc clock range
 *
 * Revision 1.19  2007-02-23 14:47:04  ddangelo
 * fixed data types
 *
 * Revision 1.18  2006/09/05 12:06:03  ddangelo
 * cleaned up
 * histograms filled
 *
 * Revision 1.17  2006/08/21 15:34:01  ddangelo
 * test histogram added (new barn interface).
 * introduced the possibility to use mode instead of mean (parameter switch).
 *
 * Revision 1.16  2006/01/25 13:12:00  misiaszek
 * Added a missing include (auth from maintainer)
 *
 * Revision 1.15  2005/03/04 14:21:09  ddangelo
 * fixed a typo
 *
 * Revision 1.14  2005/03/02 16:02:06  ddangelo
 * included trg type check to follow framework changes
 *
 * Revision 1.13  2004/11/26 15:25:11  razeto
 * Thanks to Davide who so kindly underlined a syntactical error, I fixed the "Maintainer" word
 *
 * Revision 1.12  2004/11/26 14:06:23  razeto
 * Added Mantainer field
 *
 * Revision 1.11  2004/09/22 13:29:16  ddangelo
 * updated bx_reco_event to bx_echidna_event
 *
 * Revision 1.10  2004/09/22 11:29:58  ddangelo
 * added require_event_stage in module's constructor
 *
 * Revision 1.9  2004/06/01 11:36:27  ddangelo
 * updated to new constants syntax
 *
 * Revision 1.8  2004/05/25 12:59:33  ddangelo
 * std::vector<>::at() replaced by std::vector<>::operator[]()
 * to compile with gcc 2.95 as well
 *
 * Revision 1.7  2004/05/20 15:41:53  ddangelo
 * check on has_data introduced
 * check on chann descr implemented with db data (db_profile for the moment)
 * instead of static event members.
 * computetion done with algotithms instead of c like math.
 * Data stored in a std::vector* (histo use).
 *
 * Revision 1.5  2004/04/28 11:04:49  ddangelo
 * Modules adapted to call the right const/nonconst event portion getter
 *
 * Revision 1.4  2004/04/27 15:05:53  ddangelo
 * fixed the bug that was causing the crash with no muon data.
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
