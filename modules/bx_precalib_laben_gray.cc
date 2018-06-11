/* BOREXINO Reconstruction program
 *
 * Author: Alessandro Razeto <Alessandro.Razeto@ge.infn.it>
 * Maintainer: Alessandro Razeto <Alessandro.Razeto@ge.infn.it>
 *
 * $Id: bx_precalib_laben_gray.cc,v 1.29 2007/02/06 14:37:43 razeto Exp $
 *
 * Implementation of bx_precalibrator
 *
 */
#include "bx_precalib_laben_gray.hh"
#include "messenger.hh"
#include "laben_time_hit.hh"
#include "constants.hh"
#include "bx_dbi.hh"
#include "db_run.hh"
#include "db_profile.hh"
#include "barn_interface.hh"

#include <algorithm>
#include <numeric>


// ctor
bx_precalib_laben_gray::bx_precalib_laben_gray (): bx_base_module("bx_precalib_laben_gray", bx_base_module::precalib_cycle1) {
  require_event_stage (bx_detector::laben, bx_base_event::raw);
  require_trigger_type (bx_trigger_event::pulser);
}

// module interface
void bx_precalib_laben_gray::begin () {

  i_channel_skip = get_parameter ("channel_reference_skip").get_int ();
  f_mean_bound = get_parameter ("mean_bound").get_float ();
  f_rms_bound = get_parameter ("rms_bound").get_float ();

  i_nreferences = int32_t(constants::laben::channels / i_channel_skip);
  gray_diffs = new TH2F* [i_nreferences];
  for (int32_t i = 0; i < i_nreferences; i++) {
    std::ostringstream name, title;
    name << "gray_diff_" << (i + 1) * i_channel_skip;
    title << "bx_precalib_laben_gray TDC gray counter diffs for lg " << (i + 1) * i_channel_skip;
    gray_diffs[i] = new TH2F (name.str ().c_str (), title.str ().c_str (), constants::laben::channels, 1, constants::laben::channels + 1, 256, -128, 128);
    barn_interface::get ()->store (barn_interface::file, gray_diffs[i], this);
  }
  gray_shifts = new TH1S ("gray_shifts", "bx_precalib_laben_gray TDC gray counter shift", 2400, 1, 2401);
  barn_interface::get ()->store (barn_interface::file, gray_shifts, this);

  gray_map = new int32_t[i_nreferences];
}

bx_echidna_event* bx_precalib_laben_gray::doit (bx_echidna_event *ev) {
  const bx_laben_event& er = ev->get_laben();
  bool count[100] = { 0, };
  

    // Update has data field
  if (er.get_raw_nhits ()) b_has_data = true;

    // First set the refernce channels
  for (int32_t i = 0; i < er.get_raw_nhits (); i++) {
    const bx_laben_raw_hit &hit = er.get_raw_hit (i);
    for (int32_t j = 0; j < constants::laben::channels / i_channel_skip; j++) 
      if (hit.get_logical_channel () == i_channel_skip * (j + 1) && hit.get_time_1 () != 0xff && hit.get_time_2 () != 0xff)
	if (!count[j]) {
          gray_map[j] = hit.get_gray_counter ();
	}
  }
  
    // Then fill the shift histogram using only first hit
  int32_t fired_channels[constants::laben::channels];
  std::fill_n (fired_channels, constants::laben::channels, 0);
  for (int32_t i = 0; i < er.get_raw_nhits (); i++) {
    const bx_laben_raw_hit &hit = er.get_raw_hit (i);
    if (hit.get_time_1 () == 0xff && hit.get_time_2 () == 0xff) continue;
    int32_t lg = hit.get_logical_channel ();
    if (fired_channels[lg - 1]++) continue;
    int32_t gray = hit.get_gray_counter ();
    for (int32_t j = 0; j < i_nreferences; j++) gray_diffs[j]->Fill (lg, gray - gray_map[j]);
  }
  return ev;
}

void bx_precalib_laben_gray::end () {
  if(b_has_data) {

    db_run& run_info = bx_dbi::get ()->get_run ();
    bool found = false;

    for (int32_t j = 0; j < constants::laben::channels / i_channel_skip; j++) {
      float mean, rms;
      mean = gray_diffs[j]->GetMean (2);
      rms = gray_diffs[j]->GetRMS (2);
      if (::fabs (mean) < f_mean_bound && rms < f_rms_bound && rms > 0.1) { 
        get_message (bx_message::log) << "found channel reference " << i_channel_skip * (j + 1) << " with shift " << mean << "+-" << rms << dispatch;
        found = true;

        for (int32_t i = 0; i < constants::laben::channels; i++) {
	  TH1D * p = gray_diffs[j]->ProjectionY ("base_one_lg", i + 1, i + 1);
	  int32_t gray_shift = p->GetMaximumBin () - 129; // 129 is the center
	  if (p->Integral () < 100) gray_shift = 0; // ROOT SUCKS!!! GetEntries does not work with older root version
	  if (gray_shift < -20 || gray_shift > 20) get_message (bx_message::info) << "ch " << i + 1 << " shift " << gray_shift << " out of bounds" << dispatch;
	
	  run_info.set_laben_precalib_gray_shift (i + 1, gray_shift, this);
  

	  get_message (bx_message::debug) << "gray shift for channel " << i + 1 << " is " << gray_shift << dispatch;
	  gray_shifts->SetBinContent (i + 1, gray_shift); // Bins start at 1
        }
        break;
      } else 
        get_message (bx_message::log) << "skipped channel reference " << i_channel_skip * (j + 1) << " with shift " << mean << "+-" << rms << dispatch;
    }
    if (!found) get_message (bx_message::critic) << "no gray shift reference channel found" << dispatch;    

  }

    // delete the matrix
  delete [] gray_map;
}
/*
 * $Log: bx_precalib_laben_gray.cc,v $
 * Revision 1.29  2007/02/06 14:37:43  razeto
 * ROOT SUCKSvi modules/bx_precalib_laben_gray.cc ! OLD root does not calculate correctly entries for projections. Fixed
 *
 * Revision 1.28  2007/02/02 16:51:39  razeto
 * Check that the histogram contains data
 *
 * Revision 1.27  2007/01/29 19:03:04  razeto
 * Use root histogram to calulate mean and rms. This will produce shift in barn
 *
 * Revision 1.26  2007/01/27 15:35:32  razeto
 * Handle reference channels in the same way of ordinary ones
 *
 * Revision 1.25  2006-08-21 11:16:46  razeto
 * Updated to new barn_interface
 *
 * Revision 1.24  2006/01/02 21:24:39  razeto
 * Changed test target to file target for db_barn
 *
 * Revision 1.23  2005/12/03 15:22:35  razeto
 * Added a message when recalculated precalib values differ from the database ones
 *
 * Revision 1.22  2005/08/03 13:24:12  razeto
 * Fixed the lenght distribution histograms
 *
 * Revision 1.21  2005/07/13 12:34:41  razeto
 * Fixed a typo
 *
 * Revision 1.20  2005/03/22 12:10:19  razeto
 * Fixed a message
 *
 * Revision 1.19  2005/03/02 15:45:43  razeto
 * Fixed a test (buggy for older runs) from framework to precalib modules
 *
 * Revision 1.18  2004/11/26 15:25:11  razeto
 * Thanks to Davide who so kindly underlined a syntactical error, I fixed the "Maintainer" word
 *
 * Revision 1.17  2004/11/26 14:01:21  razeto
 * Added Mantainer field
 *
 * Revision 1.16  2004/09/30 13:03:42  razeto
 * Fixed a printout syntax
 *
 * Revision 1.15  2004/09/30 12:54:24  razeto
 * Changed a printout level from info to log
 *
 * Revision 1.14  2004/09/28 13:43:27  razeto
 * Removed the ifdef for root barn histos
 *
 * Revision 1.13  2004/09/22 13:25:39  razeto
 * Changed bx_reco_event to bx_echidna_event
 *
 * Revision 1.12  2004/09/22 10:34:04  razeto
 * Updated to follow sub_detector enum in bx_detector
 *
 * Revision 1.11  2004/08/10 11:05:53  razeto
 * Fixed a bug in finding valid reference channel
 *
 * Revision 1.10  2004/07/12 10:47:10  razeto
 * Updated to use only the first hit on every channel
 *
 * Revision 1.9  2004/06/01 11:50:36  razeto
 * Updated to follow new constants design
 *
 * Revision 1.8  2004/05/31 15:00:55  razeto
 * Updated to use the base module require statements.
 * Upgraded to avoid doing stuff at the end if there is no data (has_data is false)
 *
 * Revision 1.7  2004/05/21 11:28:23  razeto
 * Updated conforming to the new conventions of event::get_logical_channel
 *
 * Revision 1.6  2004/05/21 08:39:13  razeto
 * Updated to use barn_interface
 *
 * Revision 1.5  2004/05/18 15:01:54  razeto
 * A lot of development done
 *
 * Revision 1.4  2004/04/28 11:04:49  ddangelo
 * Modules adapted to call the right const/nonconst event portion getter
 *
 * Revision 1.3  2004/04/27 14:25:18  ddangelo
 * Modifications to match new event structure.
 * Useless dynamic_cast in event reading replaced by getter.
 * safer const reference is used instead of non-const pointer.
 *
 * Revision 1.2  2004/04/26 13:53:00  razeto
 * Updated to follow new name/syntax of db_run and db_profile
 *
 * Revision 1.1  2004/04/12 15:59:49  razeto
 * Added a firsto core of laben precalibration
 *
 */
