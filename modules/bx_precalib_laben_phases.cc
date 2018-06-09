/* BOREXINO Reconstruction program
 *
 * Author: Alessandro Razeto <Alessandro.Razeto@ge.infn.it>
 * Maintainer: Alessandro Razeto <Alessandro.Razeto@ge.infn.it>
 *
 * $Id: bx_precalib_laben_phases.cc,v 1.21 2006/08/21 11:16:46 razeto Exp $
 *
 * Implementation of bx_precalibrator
 *
 */
#include "bx_precalib_laben_phases.hh"
#include "messenger.hh"
#include "laben_time_hit.hh"
#include "constants.hh"
#include "bx_dbi.hh"
#include "db_run.hh"
#include "barn_interface.hh"

// ctor
bx_precalib_laben_phases::bx_precalib_laben_phases (): bx_base_module("bx_precalib_laben_phases", bx_base_module::precalib_cycle2) {
  require_event_stage (bx_detector::laben, bx_base_event::raw);
  require_trigger_type (bx_trigger_event::pulser);
}

// module interface
void bx_precalib_laben_phases::begin () {
  right_phase_v = new int[constants::laben::channels];
  inverse_phase_v = new int[constants::laben::channels];
  error_phase_v = new int[constants::laben::channels];

  std::fill_n (right_phase_v, constants::laben::channels, 0);
  std::fill_n (inverse_phase_v, constants::laben::channels, 0);
  std::fill_n (error_phase_v, constants::laben::channels, 0);

  phase_map = new TH1S("phase_map", "bx_precalib_laben_phases TDC phase map", 2400, 1, 2401);
  phase_channel_distribution = new TH2S ("phase_channel_distribution", "bx_precalib_laben_phases TDC phase distribution", constants::laben::channels, 1, constants::laben::channels + 1, 3, 0, 3);
  barn_interface::get ()->store (barn_interface::file, phase_map, this);
  barn_interface::get ()->store (barn_interface::file, phase_channel_distribution, this);
  get_message(bx_message::debug) << "begin" << dispatch;
}

bx_echidna_event* bx_precalib_laben_phases::doit (bx_echidna_event *ev) {
  const bx_laben_event& er = ev->get_laben();

    // Update has data field
  if (er.get_raw_nhits ()) b_has_data = true;

    // Use only first hit
  int fired_channels[constants::laben::channels];
  std::fill_n (fired_channels, constants::laben::channels, 0);
  
  laben_time_hit t_hit;
  for (int i = 0; i < er.get_raw_nhits (); i++) {
    const bx_laben_raw_hit &hit = er.get_raw_hit (i);
    if (fired_channels[hit.get_logical_channel () - 1]++) continue;
    t_hit.init (hit);
    laben_time_hit::ramp_slope guessed = t_hit.guess_slope ();
    if (!t_hit.is_good ()) {
      error_phase_v[hit.get_logical_channel () - 1] ++;
      phase_channel_distribution->Fill (hit.get_logical_channel (), 2);
    }
    if (guessed != t_hit.hw_slope ()) {
      inverse_phase_v[hit.get_logical_channel () - 1] ++;
      phase_channel_distribution->Fill (hit.get_logical_channel (), 0);
    } else {
      right_phase_v[hit.get_logical_channel () - 1] ++;
      phase_channel_distribution->Fill (hit.get_logical_channel (), 1);
    }
  }
  return ev;
}

void bx_precalib_laben_phases::end () {
  if (b_has_data) {
  
    db_run& run_info = bx_dbi::get ()->get_run ();
  
    for (int i = 0; i < constants::laben::channels; i++) {
      bool rising_on_even = true;
      if (error_phase_v[i] * 2 >= (inverse_phase_v[i] + right_phase_v[i]) && error_phase_v[i]) 
        get_message(bx_message::info) << "channel " << i + 1 << " too many errors " << error_phase_v[i] << " against " << inverse_phase_v[i] << " + " << right_phase_v[i] << dispatch;

      if (inverse_phase_v[i] > right_phase_v[i]) rising_on_even = false;

      int v = rising_on_even;
      if (error_phase_v[i] >= (inverse_phase_v[i] + right_phase_v[i])) {
        if (error_phase_v[i]) v = 2;
        else v = -1;
      }
      phase_map->SetBinContent(i + 1, v);

      if (run_info.check_laben_precalib_rising_on_even (i + 1) && run_info.get_laben_precalib_rising_on_even (i + 1) != rising_on_even)
	get_message(bx_message::log) << "replacing rising on even from " << run_info.get_laben_precalib_rising_on_even (i + 1) << " to " << rising_on_even << " for channel " << i + 1 << dispatch;
      run_info.set_laben_precalib_rising_on_even (i + 1, rising_on_even, this);
      get_message(bx_message::debug) << "rising_on_even for channel " << i + 1 << " is " << rising_on_even << dispatch;
    }
  }

    // deallocate the vectors
  delete [] right_phase_v;
  delete [] inverse_phase_v;
  delete [] error_phase_v;

  get_message(bx_message::debug) << "end" << dispatch;
}
/*
 * $Log: bx_precalib_laben_phases.cc,v $
 * Revision 1.21  2006/08/21 11:16:46  razeto
 * Updated to new barn_interface
 *
 * Revision 1.20  2006/01/02 21:24:39  razeto
 * Changed test target to file target for db_barn
 *
 * Revision 1.19  2005/12/03 15:22:35  razeto
 * Added a message when recalculated precalib values differ from the database ones
 *
 * Revision 1.18  2005/08/03 13:24:12  razeto
 * Fixed the lenght distribution histograms
 *
 * Revision 1.17  2005/03/22 12:10:19  razeto
 * Fixed a message
 *
 * Revision 1.16  2005/03/02 17:33:35  razeto
 * Updated to new laben_time_hit
 *
 * Revision 1.15  2005/03/02 15:45:43  razeto
 * Fixed a test (buggy for older runs) from framework to precalib modules
 *
 * Revision 1.14  2004/11/26 15:25:11  razeto
 * Thanks to Davide who so kindly underlined a syntactical error, I fixed the "Maintainer" word
 *
 * Revision 1.13  2004/11/26 14:01:21  razeto
 * Added Mantainer field
 *
 * Revision 1.12  2004/09/28 13:43:27  razeto
 * Removed the ifdef for root barn histos
 *
 * Revision 1.11  2004/09/22 13:25:39  razeto
 * Changed bx_reco_event to bx_echidna_event
 *
 * Revision 1.10  2004/09/22 10:34:04  razeto
 * Updated to follow sub_detector enum in bx_detector
 *
 * Revision 1.9  2004/07/12 10:47:10  razeto
 * Updated to use only the first hit on every channel
 *
 * Revision 1.8  2004/06/01 11:50:36  razeto
 * Updated to follow new constants design
 *
 * Revision 1.7  2004/05/31 15:00:55  razeto
 * Updated to use the base module require statements.
 * Upgraded to avoid doing stuff at the end if there is no data (has_data is false)
 *
 * Revision 1.6  2004/05/21 11:28:23  razeto
 * Updated conforming to the new conventions of event::get_logical_channel
 *
 * Revision 1.5  2004/05/21 08:39:13  razeto
 * Updated to use barn_interface
 *
 * Revision 1.4  2004/05/18 15:01:54  razeto
 * A lot of development done
 *
 * Revision 1.3  2004/04/28 11:04:49  ddangelo
 * Modules adapted to call the right const/nonconst event portion getter
 *
 * Revision 1.2  2004/04/27 14:25:18  ddangelo
 * Modifications to match new event structure.
 * Useless dynamic_cast in event reading replaced by getter.
 * safer const reference is used instead of non-const pointer.
 *
 * Revision 1.1  2004/04/12 15:59:49  razeto
 * Added a firsto core of laben precalibration
 *
 */
