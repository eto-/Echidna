/* BOREXINO Reconstruction program
 *
 * Author: Alessandro Razeto <Alessandro.Razeto@ge.infn.it>
 * Maintainer: Alessandro Razeto <Alessandro.Razeto@ge.infn.it>
 *
 * $Id: bx_precalib_laben_adc.cc,v 1.22 2006/08/21 11:16:46 razeto Exp $
 *
 * Implementation of bx_precalibrator
 *
 */
#include "bx_precalib_laben_adc.hh"
#include "messenger.hh"
#include "laben_time_hit.hh"
#include "constants.hh"
#include "bx_dbi.hh"
#include "db_run.hh"
#include "barn_interface.hh"

#include <algorithm>

// ctor
bx_precalib_laben_adc::bx_precalib_laben_adc (): bx_base_module("bx_precalib_laben_adc", bx_base_module::precalib_cycle1) {
  require_event_stage (bx_detector::laben, bx_base_event::raw);
  require_trigger_type (bx_trigger_event::pulser);
}

// module interface
void bx_precalib_laben_adc::begin () {
    // allocate the matrix
  adc_sample_map = new int32_t*[constants::laben::channels]; 
  for (int32_t i = 0; i < constants::laben::channels; i++) {
    adc_sample_map[i] = new int32_t[256];
    std::fill_n (adc_sample_map[i], 256, 0);
  }

  u1_maxima_bin_add = get_parameter ("maxima_bin_add").get_int ();
  
  adc_samples = new TH1F ("adc_samples", "bx_precalib_laben_adc TDC RAMP samples", 256, 0, 256);
  adc_limits = new TH1F ("adc_limits", "bx_precalib_laben_adc TDC RAMP limits", 256, 0, 256);
  barn_interface::get ()->store (barn_interface::file, adc_samples, this);
  barn_interface::get ()->store (barn_interface::file, adc_limits, this);

  get_message(bx_message::debug) << "begin" << dispatch;
}

bx_echidna_event* bx_precalib_laben_adc::doit (bx_echidna_event *ev) {
  const bx_laben_event& er = ev->get_laben();

    // Update has data field
  if (er.get_raw_nhits ()) b_has_data = true;

    // Use only first hit
  int32_t fired_channels[constants::laben::channels];
  std::fill_n (fired_channels, constants::laben::channels, 0);
      
  for (int32_t i = 0; i < er.get_raw_nhits (); i++) {
    const bx_laben_raw_hit &hit = er.get_raw_hit (i);
    if (fired_channels[hit.get_logical_channel () - 1]++) continue;
    adc_sample_map[hit.get_logical_channel () - 1][hit.get_time_1()] ++;
    adc_sample_map[hit.get_logical_channel () - 1][hit.get_time_2()] ++;
    adc_samples->Fill (hit.get_time_1());
    adc_samples->Fill (hit.get_time_2());
  }
  return ev;
}

void bx_precalib_laben_adc::end () {
  if (b_has_data) {

    db_run& run_info = bx_dbi::get ()->get_run ();

    for (int32_t i = 0; i < constants::laben::channels; i++) {
        // calculate local maxima in the first 20 and in the last 20 bins
      int32_t low_bin = std::max_element (adc_sample_map[i], adc_sample_map[i] + 20) - adc_sample_map[i];
      int32_t high_bin = std::max_element (adc_sample_map[i] + 236, adc_sample_map[i] + 256) - adc_sample_map[i];
    
        // add/remove a costant offset
      if (low_bin >= u1_maxima_bin_add) low_bin -= u1_maxima_bin_add;
      if (high_bin + u1_maxima_bin_add <= 255) high_bin += u1_maxima_bin_add;
    
        // if adc_sample_map[i][high_bin] == 0 no hit present
      if (! adc_sample_map[i][high_bin]) high_bin = 255; // low_bin is already at 0 considering max_element engine

        // set the values
      run_info.set_laben_precalib_low_bin (i + 1, low_bin, this);
      run_info.set_laben_precalib_high_bin (i + 1, high_bin, this);

      adc_limits->Fill (low_bin);
      adc_limits->Fill (high_bin);
    
      get_message(bx_message::debug) << "adc limits for channel " << i + 1 << " are " <<  low_bin << ":" << high_bin << dispatch;
    } 
  } 
  
    // deallocate the matrix
  for (int32_t i = 0; i < constants::laben::channels; i++) delete [] adc_sample_map[i];
  delete [] adc_sample_map;

  get_message(bx_message::debug) << "end" << dispatch;
}
/*
 * $Log: bx_precalib_laben_adc.cc,v $
 * Revision 1.22  2006/08/21 11:16:46  razeto
 * Updated to new barn_interface
 *
 * Revision 1.21  2006/01/02 21:28:47  razeto
 * Changed test target to file target for db_barn
 *
 * Revision 1.20  2005/05/16 20:05:39  razeto
 * Fixed a bug (thanks to marcin)
 *
 * Revision 1.19  2005/03/22 12:10:19  razeto
 * Fixed a message
 *
 * Revision 1.18  2005/03/15 09:25:16  razeto
 * Fixed a bug
 *
 * Revision 1.17  2005/03/02 15:45:43  razeto
 * Fixed a test (buggy for older runs) from framework to precalib modules
 *
 * Revision 1.16  2004/11/26 15:25:11  razeto
 * Thanks to Davide who so kindly underlined a syntactical error, I fixed the "Maintainer" word
 *
 * Revision 1.15  2004/11/26 14:01:20  razeto
 * Added Mantainer field
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
 * Revision 1.11  2004/07/12 10:47:10  razeto
 * Updated to use only the first hit on every channel
 *
 * Revision 1.10  2004/06/01 11:50:36  razeto
 * Updated to follow new constants design
 *
 * Revision 1.9  2004/05/31 15:00:55  razeto
 * Updated to use the base module require statements.
 * Upgraded to avoid doing stuff at the end if there is no data (has_data is false)
 *
 * Revision 1.8  2004/05/21 11:28:23  razeto
 * Updated conforming to the new conventions of event::get_logical_channel
 *
 * Revision 1.7  2004/05/21 08:39:13  razeto
 * Updated to use barn_interface
 *
 * Revision 1.6  2004/05/18 15:01:54  razeto
 * A lot of development done
 *
 * Revision 1.5  2004/04/28 11:04:49  ddangelo
 * Modules adapted to call the right const/nonconst event portion getter
 *
 * Revision 1.4  2004/04/27 14:25:18  ddangelo
 * Modifications to match new event structure.
 * Useless dynamic_cast in event reading replaced by getter.
 * safer const reference is used instead of non-const pointer.
 *
 * Revision 1.3  2004/04/26 13:53:00  razeto
 * Updated to follow new name/syntax of db_run and db_profile
 *
 * Revision 1.2  2004/04/13 12:09:17  razeto
 * Fixed a bug: the values were not set
 *
 * Revision 1.1  2004/04/12 15:59:49  razeto
 * Added a firsto core of laben precalibration
 *
 */
