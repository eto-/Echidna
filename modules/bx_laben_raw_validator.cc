/* BOREXINO Reconstruction program
 *
 * Author: Alessandro Razeto <Alessandro.Razeto@ge.infn.it>
 * Maintainer: Alessandro Razeto <Alessandro.Razeto@ge.infn.it>
 *
 * $Id: bx_laben_raw_validator.cc,v 1.15 2010/07/01 18:22:43 razeto Exp $
 *
 */

#include "bx_laben_raw_validator.hh"
#include "bx_echidna_event.hh"
#include "barn_interface.hh"
#include "TH1F.h"
#include "constants.hh"
#include "laben_time_hit.hh"
#include <algorithm>

// ctor
bx_laben_raw_validator::bx_laben_raw_validator (): bx_base_module("bx_laben_raw_validator", bx_base_module::main_loop) {
  require_event_stage (bx_detector::laben, bx_base_event::raw);
}

// module interface
void bx_laben_raw_validator::begin () {
  flag_lg[bx_laben_raw_hit::good] = new TH1F("good_flag_lg", "Map of good hits", constants::laben::channels, 1, constants::laben::channels + 1);
  flag_lg[bx_laben_raw_hit::fifo_full] = new TH1F("fifo_full_flag_lg", "Map of FIFO FULL hits", constants::laben::channels, 1, constants::laben::channels + 1);
  flag_lg[bx_laben_raw_hit::fifo_empty] = new TH1F("fifo_empty_flag_lg", "Map of FIFO EMPTY hits", constants::laben::channels, 1, constants::laben::channels + 1);
  flag_lg[bx_laben_raw_hit::counter] = new TH1F("counter_flag_lg", "Map of FIFO EMPTY hits", constants::laben::channels, 1, constants::laben::channels + 1);
  flag_lg[bx_laben_raw_hit::trg_jump] = new TH1F("trg_jump_flag_lg", "Map of TRG JUMP hits", constants::laben::channels, 1, constants::laben::channels + 1);
  flag_lg[bx_laben_raw_hit::trg_jump_large] = new TH1F("trg_jump_large_flag_lg", "Map of TRG JUMP LARGE hits", constants::laben::channels, 1, constants::laben::channels + 1);
  flag_lg[bx_laben_raw_hit::trg_in_busy] = new TH1F("trg_in_busy_flag_lg", "Map of TRG IN BUSY hits", constants::laben::channels, 1, constants::laben::channels + 1);
  flag_lg[bx_laben_raw_hit::invalid] = new TH1F("invalid_lg", "Map of invalid hits", constants::laben::channels, 1, constants::laben::channels + 1);
  for (int32_t i = 0; i < 8; i++)
    barn_interface::get ()->store (barn_interface::file, flag_lg[i], this);

  good_flag_lg_c = (TH1F*)flag_lg[bx_laben_raw_hit::good]->Clone ("good_flag_lg_c");
  board_occupancy = new char[constants::laben::nboards];
  count = 0;
}


bx_echidna_event* bx_laben_raw_validator::doit (bx_echidna_event *ev) {
  const bx_laben_event& er = ev->get_laben ();
  bx_laben_raw_event& e = dynamic_cast<bx_laben_raw_event&>(ev->get_laben ());
  int32_t n_flags[8] = {};

  ::memset (board_occupancy, 0, constants::laben::nboards);
  //std::fill_n (board_occupancy, constants::laben::nboards, 0); memset is much faster

  int32_t last_lg = -1;
  int32_t nhits_lg = 0;
  bool invalid_lg = false;
  for (int32_t i = 0; i < e.get_raw_nhits (); i++) {
    const bx_laben_raw_hit &hit = er.get_raw_hit (i);
    if (nhits_lg == hit.get_logical_channel ()) nhits_lg ++;
    else {
      if (nhits_lg > 0 && invalid_lg) get_message (bx_message::error) << "event " << ev->get_event_number () << "laben channel " << last_lg << " have " << nhits_lg << " hits with an invalid hit" << dispatch;
      last_lg = hit.get_logical_channel ();
      invalid_lg = false; 
      nhits_lg = 0;
    }
    for (int32_t flag = 0; flag < 8; flag++) {
      if (hit.check_flag (bx_laben_raw_hit::flags(flag))) {
	n_flags[flag] ++;
	flag_lg[flag]->Fill (hit.get_logical_channel ());
      }
    }
    if (!hit.check_flag (bx_laben_raw_hit::invalid)) { // ignore invalid hit
      if (!hit.check_flag (bx_laben_raw_hit::counter)) board_occupancy[(hit.get_logical_channel () - 1)/ 8] = 1;  
      else e.i4_nhits_fw += int32_t(hit.get_base ()) + (int32_t(hit.get_peak ()) >> 8);
    } else invalid_lg = true;
  }

  int32_t zero_count = 0;
  int32_t max_board = constants::laben::nboards;
  if (ev->get_run_number () >= 12000) max_board = constants::laben::board_per_rack * 13;
  for (int32_t i = 0; i < max_board; i++) if (!board_occupancy[i]) zero_count++;
  e.i4_empty_boards = zero_count;


  if (n_flags[bx_laben_raw_hit::good]) {
    if (n_flags[bx_laben_raw_hit::trg_jump]) get_message (bx_message::log) << "event " << ev->get_event_number () << "laben jump with " << n_flags[bx_laben_raw_hit::good] << " good hits and " << n_flags[bx_laben_raw_hit::trg_jump] << " jump hits" << dispatch;
    if (n_flags[bx_laben_raw_hit::trg_jump_large]) get_message (bx_message::warn) << "event " << ev->get_event_number () << "laben jump large with " << n_flags[bx_laben_raw_hit::good] << " good hits and " << n_flags[bx_laben_raw_hit::trg_jump_large] << " jump large hits" << dispatch;
    //if (n_flags[bx_laben_raw_hit::trg_in_busy]) get_message (bx_message::debug) << "event " << ev->get_event_number () << "laben trigger in busy with " << n_flags[bx_laben_raw_hit::good] << " good hits and " << n_flags[bx_laben_raw_hit::trg_in_busy] << " trigger in busy hits" << dispatch;
  } else {
    if (n_flags[bx_laben_raw_hit::trg_jump]) get_message (bx_message::log) << "event " << ev->get_event_number () << "laben jump with " << n_flags[bx_laben_raw_hit::trg_jump] << " jump hits" << dispatch;
    if (n_flags[bx_laben_raw_hit::trg_jump_large]) get_message (bx_message::warn) << "event " << ev->get_event_number () << "laben jump large with " << n_flags[bx_laben_raw_hit::trg_jump_large] << " jump large hits" << dispatch;
    //if (n_flags[bx_laben_raw_hit::trg_in_busy])  get_message (bx_message::debug) << "event " << ev->get_event_number () << "laben trigger in busy with "  << n_flags[bx_laben_raw_hit::trg_in_busy] << " trigger in busy hits" << dispatch;
  }
  std::copy (n_flags, n_flags + 8, e.nhits_flag); 

  if (!(++count % 600)) for (int32_t i = 0; i < 8; i++) barn_interface::get ()->network_send (flag_lg[i], this);
  if (!(count % 1200)) {
    good_flag_lg_c->Add (flag_lg[bx_laben_raw_hit::good]);
    flag_lg[bx_laben_raw_hit::good]->Reset ();
  }

  return ev;
}

void bx_laben_raw_validator::end () {
  flag_lg[bx_laben_raw_hit::good]->Reset ();
  flag_lg[bx_laben_raw_hit::good]->Add (good_flag_lg_c);
}

/*
 * $Log: bx_laben_raw_validator.cc,v $
 * Revision 1.15  2010/07/01 18:22:43  razeto
 * Removed a debug printout
 *
 * Revision 1.14  2010-07-01 18:21:35  razeto
 * Added counter hit flag for new laben fw and nhits_fw from laben boards
 *
 * Revision 1.13  2010-05-21 12:35:42  razeto
 * Simpler checks
 *
 * Revision 1.12  2010-05-21 12:34:01  ddangelo
 * bx_laben_tracker renamed as bx_laben_energy_tracker
 * added bx_laben_tof_tracker and bx_cmt_tracker
 *
 * Revision 1.11  2010-05-21 11:43:34  razeto
 * Ignore crate 14 for empty boards
 *
 * Revision 1.10  2009-10-23 12:41:22  razeto
 * Muted
 *
 * Revision 1.9  2009-10-06 13:29:22  razeto
 * Check added for invalid hits
 *
 * Revision 1.8  2009-01-30 11:14:34  razeto
 * Added non cumulative plot
 *
 * Revision 1.7  2008-10-07 14:03:56  razeto
 * Using new flag
 *
 * Revision 1.6  2008-10-04 15:56:26  razeto
 * Send histos to viewer
 *
 * Revision 1.5  2007-11-05 11:03:09  razeto
 * Added static is_laben_invalid to laben_time_hit
 *
 * Revision 1.4  2007-11-03 23:42:54  razeto
 * Ignore invalid hits for empty board calculation
 *
 * Revision 1.3  2007-10-30 18:12:59  ddangelo
 * just variable renamed
 *
 * Revision 1.2  2007-10-30 17:35:14  razeto
 * Reduced verbosity, added board occupancy
 *
 * Revision 1.1  2006-10-23 14:41:16  razeto
 * New module bx_laben_raw_validator to validate event using laben raw hit flags
 *
 */
