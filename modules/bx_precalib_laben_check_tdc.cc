/* BOREXINO Reconstruction program
 *
 * Author: Alessandro Razeto <Alessandro.Razeto@ge.infn.it>
 * Maintainer: Alessandro Razeto <Alessandro.Razeto@ge.infn.it>
 *
 * $Id: bx_precalib_laben_check_tdc.cc,v 1.29 2007/03/28 14:16:16 razeto Exp $
 *
 * Implementation of bx_precalibrator
 *
 */
#include "bx_precalib_laben_check_tdc.hh"
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
bx_precalib_laben_check_tdc::bx_precalib_laben_check_tdc (): bx_base_module("bx_precalib_laben_check_tdc", bx_base_module::precalib_cycle4) {
  require_event_stage (bx_detector::laben, bx_base_event::raw);
  require_trigger_type (bx_trigger_event::pulser);
}

// module interface
void bx_precalib_laben_check_tdc::begin () {

  time_channel_distribution = new TH2F ("time_channel_distribution", "bx_precalib_laben_check_tdc TDC time distribution", constants::laben::channels, 1, constants::laben::channels + 1, 600, -300, 300);
  error_channel_distribution = new TH2F ("error_channel_distribution", "bx_precalib_laben_check_tdc TDC error distribution", constants::laben::channels, 1, constants::laben::channels + 1, 100, -1, 1);
  time_bits = new TH2F ("time_bits", "bx_precalib_laben_check_tdc TDC time_bit/slope distribution", 4, 0, 4, 4, 0, 4);
  barn_interface::get ()->store (barn_interface::file, time_channel_distribution, this);
  barn_interface::get ()->store (barn_interface::file, error_channel_distribution, this);
  barn_interface::get ()->store (barn_interface::file, time_bits, this);
      
  i_ref_channel = get_parameter ("ref_channel").get_int ();
  b_found_mctruth = false;
  get_message(bx_message::debug) << "begin" << dispatch;
}

bx_echidna_event* bx_precalib_laben_check_tdc::doit (bx_echidna_event *ev) {
  if (ev->is_mctruth_enabled ()) b_found_mctruth = true;

  const bx_laben_event& er = ev->get_laben();
  
  if (er.get_raw_nhits ()) b_has_data = true;
  
  double reference_time = 0;  // The assegnation is just to avoid warns from compiler
  laben_time_hit t_hit;
    // First set the refernce channels
  for (int i = 0; i < er.get_raw_nhits (); i++) {
    const bx_laben_raw_hit &hit = er.get_raw_hit (i);
    if (hit.get_logical_channel () != i_ref_channel) continue;
    t_hit.init (hit);
    reference_time = t_hit.get_time ();
    if (t_hit.is_good (0.45)) break;
  }

  //const db_profile &profile_data = bx_dbi::get ()->get_profile ();
  const db_run &run_data = bx_dbi::get ()->get_run ();
    // Then fill the shift histogram
  for (int i = 0; i < er.get_raw_nhits (); i++) {
    float dt;
    const bx_laben_raw_hit &hit = er.get_raw_hit (i);
    int lg = hit.get_logical_channel ();
    if (hit.get_time_1 () == 0xff && hit.get_time_2 () == 0xff) continue; // check faster than fist converting hit to t_hit
 //   if (profile_data.logical_channel_description (lg) != db_profile::ordinary) continue;
    t_hit.init (hit);
    dt = t_hit.get_time () - reference_time;
    if ((dt > -110 && dt < -90) || (dt > 90 && dt < 110)) {
      int rising_on_even = run_data.get_laben_precalib_rising_on_even (lg);
      int slope = t_hit.guess_slope ();
      if (slope != t_hit.hw_slope ()) {
        int bit = hit.get_flags_ch () & 0x3;
        time_bits->Fill (bit, slope + 2 * rising_on_even);
      }
    }
    if (t_hit.is_valid ()) {
      if (dt < -300) dt = -290;
      else if (dt > 300) dt = 290;
      time_channel_distribution->Fill (lg, dt);
    }
    error_channel_distribution->Fill (lg, t_hit.get_error ());
  }
  return ev;
}

void bx_precalib_laben_check_tdc::end () {
  if (b_has_data) {
    TH1D *p = time_channel_distribution->ProjectionY ("bx_precalib_laben_check_tdc_unneded");
    if (!p) get_message (bx_message::critic) << "time_channel_distribution->ProjectionY failed" << dispatch;
    p->SetAxisRange(-120,120);
    float mean = fabsf (p->GetMean ());
    float rms = p->GetRMS ();
    if (mean <= get_parameter ("time_distrib_max_mean").get_float () && rms <= get_parameter ("time_distrib_max_rms").get_float () && rms > 0.3) {
      if (get_parameter ("write_precalib").get_bool () && !b_found_mctruth) bx_dbi::get ()->get_run ().write_laben_precalib (true, this);
      get_message (bx_message::info) << "laben precalibrations converged (" << mean << "+-" << rms << ")" << dispatch;
    } else {
      get_message (bx_message::error) << "laben precalibrations failed (" << mean << "+-" << rms << ")" << dispatch;
    }
    p->Delete ();
  }
  get_message(bx_message::debug) << "end" << dispatch;
}
/*
 * $Log: bx_precalib_laben_check_tdc.cc,v $
 * Revision 1.29  2007/03/28 14:16:16  razeto
 * Quiter
 *
 * Revision 1.28  2007-01-27 15:35:32  razeto
 * Handle reference channels in the same way of ordinary ones
 *
 * Revision 1.27  2007-01-24 16:59:13  razeto
 * Use info to say that precalibs are fine
 *
 * Revision 1.26  2006/10/13 15:31:41  razeto
 * Updated to new getter name
 *
 * Revision 1.25  2006-08-21 11:16:46  razeto
 * Updated to new barn_interface
 *
 * Revision 1.24  2006/01/02 21:24:39  razeto
 * Changed test target to file target for db_barn
 *
 * Revision 1.23  2005/08/03 13:24:12  razeto
 * Fixed the lenght distribution histograms
 *
 * Revision 1.22  2005/06/28 12:33:26  razeto
 * Fixed a bug with no entries in histogram
 *
 * Revision 1.21  2005/06/27 17:05:57  razeto
 * Fixed a bug and updated the time_bits handling
 *
 * Revision 1.20  2005/03/22 12:11:33  razeto
 * Do not use channel other than ordinary to fill the time distribution histo
 *
 * Revision 1.19  2005/03/17 16:07:41  razeto
 * Added check for mctruth presence before writing precalib data
 *
 * Revision 1.18  2005/03/15 09:26:00  razeto
 * Added a log message
 *
 * Revision 1.17  2005/03/03 15:17:03  razeto
 * Updated to the new error range of laben_time_hit
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
 * Revision 1.13  2004/11/26 15:06:13  razeto
 * Added test of convergence of the laben precalib
 *
 * Revision 1.12  2004/11/26 14:01:20  razeto
 * Added Mantainer field
 *
 * Revision 1.11  2004/10/19 18:45:51  razeto
 * Added enable of visitor writing
 *
 * Revision 1.10  2004/09/28 13:43:27  razeto
 * Removed the ifdef for root barn histos
 *
 * Revision 1.9  2004/09/22 13:25:39  razeto
 * Changed bx_reco_event to bx_echidna_event
 *
 * Revision 1.8  2004/09/22 10:34:04  razeto
 * Updated to follow sub_detector enum in bx_detector
 *
 * Revision 1.7  2004/08/10 11:12:42  razeto
 * Added a parameter
 *
 * Revision 1.6  2004/05/31 15:00:55  razeto
 * Updated to use the base module require statements.
 * Upgraded to avoid doing stuff at the end if there is no data (has_data is false)
 *
 * Revision 1.5  2004/05/25 17:16:49  razeto
 * Removed a gcc warning
 *
 * Revision 1.4  2004/05/21 11:28:23  razeto
 * Updated conforming to the new conventions of event::get_logical_channel
 *
 * Revision 1.3  2004/05/21 08:39:13  razeto
 * Updated to use barn_interface
 *
 * Revision 1.2  2004/05/20 11:02:05  razeto
 * Changed the lenght of reference time to support the time full range
 *
 * Revision 1.1  2004/05/18 15:01:54  razeto
 * A lot of development done
 *
 */
