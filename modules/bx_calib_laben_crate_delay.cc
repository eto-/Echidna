/* BOREXINO Reconstruction program
 *
 * Author: Alessandro Razeto <Alessandro.Razeto@ge.infn.it>
 * Maintainer: Alessandro Razeto <Alessandro.Razeto@ge.infn.it>
 *
 * $Id: bx_calib_laben_crate_delay.cc,v 1.4 2007/10/30 17:33:31 ddangelo Exp $
 *
 * Implementation of bx_laben_allign_crates
 *
 */
#include "bx_calib_laben_crate_delay.hh"
#include "messenger.hh"
#include "db_channel.hh"
#include "constants.hh"
#include "TH1F.h"

// ctor
bx_calib_laben_crate_delay::bx_calib_laben_crate_delay (): bx_base_module("bx_calib_laben_crate_delay", bx_base_module::main_loop) {
  require_event_stage (bx_detector::laben, bx_base_event::decoded);
}

// module interface
void bx_calib_laben_crate_delay::begin () {
  get_message(bx_message::debug) << "begin" << dispatch;
  
  p_crate_time_sums = new double[constants::laben::ncrates];
  std::fill_n (p_crate_time_sums, constants::laben::ncrates, 0.);
  i4_event_count = 0;
}


bx_echidna_event* bx_calib_laben_crate_delay::doit (bx_echidna_event *ev) {
  const bx_laben_event& er = ev->get_laben ();

  double ev_crate_time_sums[constants::laben::ncrates];
  int32_t ev_crate_hit_count[constants::laben::ncrates];
  int32_t fired_channels[constants::laben::channels];
  std::fill_n (ev_crate_time_sums, constants::laben::ncrates, 0.);
  std::fill_n (ev_crate_hit_count, constants::laben::ncrates, 0);
  std::fill_n (fired_channels, constants::laben::channels, 0);

  for (int32_t i = 0; i < er.get_decoded_nhits (); i++) {
    const bx_laben_decoded_hit &hit = er.get_decoded_hit (i);
    int32_t lg = hit.get_raw_hit ().get_logical_channel ();
    if (fired_channels[lg - 1]++ || !hit.get_db_channel ()->is_ordinary () || !hit.is_timing_good (0.4) || hit.is_out_of_gate ()) continue;
    int32_t cr = constants::crate (lg);
    double dt = hit.get_raw_time () - er.get_trigger_rawt ();
    if (::fabs (dt) > 20) continue;
    ev_crate_time_sums[cr - 1] += (hit.get_raw_time () - er.get_trigger_rawt ()); 
    ev_crate_hit_count[cr - 1] ++;
  }

  if (!ev_crate_hit_count[0]) return ev;
  ev_crate_time_sums[0] /= double(ev_crate_hit_count[0]);
  
  for (int32_t i = 1; i < constants::laben::ncrates; i++) { 
    if (ev_crate_hit_count[i]) {
      ev_crate_time_sums[i] /= double(ev_crate_hit_count[i]);
      ev_crate_time_sums[i] -= ev_crate_time_sums[0];
      get_message(bx_message::debug) << ev_crate_time_sums[i] << dispatch;
      p_crate_time_sums[i] += ev_crate_time_sums[i];
    }
  }

  i4_event_count ++;
  
  return ev;
}

void bx_calib_laben_crate_delay::end () {
  for (int32_t i = 0; i < constants::laben::ncrates; i++) get_message(bx_message::warn) << p_crate_time_sums[i] / double(i4_event_count) << " " << i << dispatch;

  delete [] p_crate_time_sums;

  get_message(bx_message::debug) << "end" << dispatch;
}

/*
 * $Log: bx_calib_laben_crate_delay.cc,v $
 * Revision 1.4  2007/10/30 17:33:31  ddangelo
 * just method renamed
 *
 * Revision 1.3  2004-11-26 15:25:11  razeto
 * Thanks to Davide who so kindly underlined a syntactical error, I fixed the "Maintainer" word
 *
 * Revision 1.2  2004/11/26 14:01:20  razeto
 * Added Mantainer field
 *
 * Revision 1.1  2004/10/05 13:50:59  razeto
 * Changed name to conform to the new bx_calib_* standard
 *
 * Revision 1.3  2004/09/22 13:27:34  razeto
 * Changed bx_reco_event to bx_echidna_event
 *
 * Revision 1.2  2004/09/22 10:36:31  razeto
 * Updated to follow sub_detector enum in bx_detector
 *
 * Revision 1.1  2004/08/10 11:14:03  razeto
 * Added a module for calculating the time differences from crate to crate
 *
 */
