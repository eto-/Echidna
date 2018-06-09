/* BOREXINO Reconstruction program
 *
 * Author: Alessandro Razeto <Alessandro.Razeto@ge.infn.it>
 * Maintainer: Alessandro Razeto <Alessandro.Razeto@ge.infn.it>
 *
 * $Id: bx_precalib_laben_d80.cc,v 1.21 2006/08/21 11:16:46 razeto Exp $
 *
 * Implementation of bx_precalibrator
 * D80 is istogrammed in a bin large 600 bin with the following method:
 * each bin has 0.1 resolution and is centered on the 80 value; so the
 * range is 50-109.
 *
 */
#include "bx_precalib_laben_d80.hh"
#include "messenger.hh"
#include "laben_time_hit.hh"
#include "constants.hh"
#include "bx_dbi.hh"
#include "db_run.hh"
#include "barn_interface.hh"

#include <algorithm>

// ctor
bx_precalib_laben_d80::bx_precalib_laben_d80 (): bx_base_module("bx_precalib_laben_d80", bx_base_module::precalib_cycle3) {
  require_event_stage (bx_detector::laben, bx_base_event::raw);
  require_trigger_type (bx_trigger_event::pulser);
}

// module interface
void bx_precalib_laben_d80::begin () {
    // allocate the matrix
  d80_delay_map = new int*[constants::laben::channels]; 
  for (int i = 0; i < constants::laben::channels; i++) {
    d80_delay_map[i] = new int[600];
    std::fill_n (d80_delay_map[i], 600, 0);
  }

  f_d80_low_limit = get_parameter ("d80_low_limit").get_float ();
  f_d80_high_limit = get_parameter ("d80_high_limit").get_float ();
  
  d80 = new TH1F ("d80", "bx_precalib_laben_d80 TDC D80 values", 600, 50, 110);
  d80_channel_distrubution = new TH2F ("d80_channel_distrubution", "bx_precalib_laben_d80 TDC D80 channel distribution", constants::laben::channels, 1, constants::laben::channels + 1, 600, 50, 110);
  barn_interface::get ()->store (barn_interface::file, d80, this);
  barn_interface::get ()->store (barn_interface::file, d80_channel_distrubution, this);
      
  get_message(bx_message::debug) << "begin" << dispatch;
}

bx_echidna_event* bx_precalib_laben_d80::doit (bx_echidna_event *ev) {
  const bx_laben_event& er = ev->get_laben ();

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
    float d80 = t_hit.get_d80 ();
    if (!t_hit.is_good (0.25)) continue;
    if (d80 < 50.) d80 = 50.;
    else if (d80 > 109.) d80 = 109.;
    d80_delay_map[hit.get_logical_channel () - 1][int ((d80 - 50.) * 10.)] ++;
    d80_channel_distrubution->Fill (hit.get_logical_channel (), d80);
  }
  return ev;
}

void bx_precalib_laben_d80::end () {
  if (b_has_data) {

    db_run& run_info = bx_dbi::get ()->get_run ();
  
    for (int i = 0; i < constants::laben::channels; i++) {
        // calculate the histogram mean
      double sum = 0.;
      long int counts = 0;
      for (int j = 0; j < 600; j++) {
         counts += d80_delay_map[i][j];
         sum += (float(j) / 10 + 50) * d80_delay_map[i][j];
      }
      double mean = sum / counts;

      if (counts < 20) mean = 80.1235;
      else if (mean < f_d80_low_limit || mean > f_d80_high_limit) {
        get_message (bx_message::info) << "d80 not found for channel " << i + 1 << " between limits (" << f_d80_low_limit << ":" << f_d80_high_limit << ")" << dispatch;
        mean = 56.;
      }
    
        // set the value
      run_info.set_laben_precalib_delta80 (i + 1, mean, this);

      d80->Fill (mean);
    
      get_message (bx_message::debug) << "d80 for channel " << i + 1 << " is " << mean << dispatch;
    }
  }

    // deallocate the matrix
  for (int i = 0; i < constants::laben::channels; i++) delete [] d80_delay_map[i];
  delete [] d80_delay_map;

  get_message(bx_message::debug) << "end" << dispatch;
}
/*
 * $Log: bx_precalib_laben_d80.cc,v $
 * Revision 1.21  2006/08/21 11:16:46  razeto
 * Updated to new barn_interface
 *
 * Revision 1.20  2006/01/02 21:24:39  razeto
 * Changed test target to file target for db_barn
 *
 * Revision 1.19  2005/08/03 13:24:12  razeto
 * Fixed the lenght distribution histograms
 *
 * Revision 1.18  2005/03/22 12:12:59  razeto
 * Fixed a message and changed the default value for no data from 55 to 80.1235 (this helps bx_elec)
 *
 * Revision 1.17  2005/03/02 17:33:35  razeto
 * Updated to new laben_time_hit
 *
 * Revision 1.16  2005/03/02 15:45:43  razeto
 * Fixed a test (buggy for older runs) from framework to precalib modules
 *
 * Revision 1.15  2004/11/26 15:25:11  razeto
 * Thanks to Davide who so kindly underlined a syntactical error, I fixed the "Maintainer" word
 *
 * Revision 1.14  2004/11/26 14:01:21  razeto
 * Added Mantainer field
 *
 * Revision 1.13  2004/10/19 18:45:11  razeto
 * Improved handling of unplugged channel
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
 * Revision 1.1  2004/04/26 13:50:51  razeto
 * Added
 *
 */
