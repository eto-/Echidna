/* BOREXINO Reconstruction program
 *
 * Author: Davide D'Angelo <davide.dangelo@lngs.infn.it>
 * Maintainer: Davide D'Angelo <davide.dangelo@lngs.infn.it>
 *
 * $Id: bx_muon_decoder.cc,v 1.43 2012/04/27 16:02:56 ddangelo Exp $
 *
 * Implementation of bx_muon_decoder.hh
 * 
*/

#include "bx_echidna_event.hh"
#include "bx_muon_decoder.hh"
#include "bx_dbi.hh"
#include "db_run.hh"
#include "db_profile.hh"
#include "db_channel.hh"
#include "bx_trigger_event.hh"
#include "bx_detector.hh"
#include <algorithm>
 
// trivial predicates for the find_if algorithms (see below)
bool is_integrity_channel (const bx_muon_raw_hit& h) {
  const db_profile& profile_info = bx_dbi::get ()->get_profile ();
  return profile_info.logical_channel_description(h.get_logical_channel()) == db_profile::integrity;
}

bool is_laser_reference_channel (const bx_muon_raw_hit& h) {
  const db_profile& profile_info = bx_dbi::get ()->get_profile ();
  return profile_info.logical_channel_description(h.get_logical_channel()) == db_profile::laser;
}

/*bool is_trigger_reference_channel (const bx_muon_raw_hit& h) {
  const db_profile& profile_info = bx_dbi::get ()->get_profile ();
  return profile_info.logical_channel_description(h.get_logical_channel()) == db_profile::trigger;
}*/

// ctor
bx_muon_decoder::bx_muon_decoder() : bx_base_module("bx_muon_decoder", bx_base_module::main_loop) {
  require_event_stage (bx_detector::muon, bx_base_event::raw);
}


// public interface methods 
     
void bx_muon_decoder::begin () {
  b_discard_disabled_lg = get_parameter ("discard_disabled_lg").get_bool ();
  b_apply_calibration   = get_parameter ("apply_calibration"  ).get_bool ();
  p_disabled_mch = new bool[constants::muon::channels];
  std::fill_n (p_disabled_mch, constants::muon::channels, false);
  u4_n_disabled_channels = 0;
  if (!bx_dbi::get ()->get_run().is_muon_alignment_present())
    get_message(bx_message::warn) << "Muon alignment not present in db" << dispatch;
}

bx_echidna_event* bx_muon_decoder::doit (bx_echidna_event *ev) {
  const std::vector<int32_t>& v = detector_interface::get ()->get_disabled_channels ();
  if (v.size () != u4_n_disabled_channels) {
    for (unsigned i = 0; i < v.size (); i++) {
      if (v[i] > constants::muon::channel_offset && v[i] <= constants::muon::channel_offset + constants::laben::channels)
        p_disabled_mch[v[i]-constants::muon::channel_offset-1] = true;
    }
    u4_n_disabled_channels = v.size ();
  }

  // cast to a muon event  
  bx_muon_event& er = ev->get_muon();
  bx_muon_decoded_event &ew = dynamic_cast<bx_muon_decoded_event&>(ev->get_muon());

  const db_run&     run_info     = bx_dbi::get ()->get_run     ();
  const db_profile& profile_info = bx_dbi::get ()->get_profile ();

  if ( !run_info.is_muon_alignment_present() || (ev->get_event_number() <= (unsigned)run_info.get_muon_aligned_nevents()) ) 
      ew.b_is_aligned = true; 
  else return ev;

  bool is_new_tdc = (er.get_error_flag() == 42) ? true : false; 

  float time_buffer_ns = is_new_tdc ? 20000 : profile_info.muon_tdc_range() * constants::muon::tdc::ns_per_clock;

  // find integrity times for the 2 tdc boards (used for NEW tdcs only)
  double integrity_time      [constants::muon::tdc::max_boards]; // will stay 0 if not a laser event
  int32_t    n_integrity_channels[constants::muon::tdc::max_boards];
  for (int32_t i=0; i < constants::muon::tdc::max_boards; i++) {
    integrity_time[i] = 0.;
    n_integrity_channels[i] = 0;
  }
  if (is_new_tdc) {
    for (int32_t iHit=0 ; iHit<er.get_raw_nhits(); iHit++) {    
      const bx_muon_raw_hit& h = er.get_raw_hit(iHit);
      if (!is_integrity_channel(h)) continue;

      integrity_time[h.get_muon_channel() < constants::muon::tdc::channels_per_board] += h.get_lead_time();
      n_integrity_channels[h.get_muon_channel() < constants::muon::tdc::channels_per_board]++;
    }
    for (int32_t i=0; i < constants::muon::tdc::max_boards; i++) {
      if (n_integrity_channels[i] != 7)
        get_message(bx_message::warn) << "Strange number of integrity channels (" << n_integrity_channels[i] << ") on board " << i << dispatch;
      if (n_integrity_channels[i]) integrity_time[i] /= n_integrity_channels[i];
     }
  }

  // find laser reference time
  std::vector<bx_muon_raw_hit>::const_iterator laser_ref;
  double laser_ref_time = 0.; // will stay 0 if not a laser event
  if (ev->get_trigger().is_laser394()) {
    laser_ref = find_if(er.get_raw_hits().begin(), er.get_raw_hits().end(), is_laser_reference_channel);
    if(laser_ref == er.get_raw_hits().end()) {
      get_message(bx_message::warn) << "Laser event without reference pulse" << dispatch;
      return ev;
    }
    laser_ref_time = laser_ref->get_lead_time();
    if (is_new_tdc) laser_ref_time -= integrity_time[laser_ref->get_muon_channel()<128];
  }

  // find trigger reference time
/*  std::vector<bx_muon_raw_hit>::const_iterator trigger_ref;
  double trigger_ref_time = 0.; // will stay 0 if old TDC or if it is a laser event
  if (is_new_tdc && !ev->get_trigger().is_laser394()) {
    trigger_ref = find_if(er.get_raw_hits().begin(), er.get_raw_hits().end(), is_trigger_reference_channel);
    if(trigger_ref == er.get_raw_hits().end()) {
      get_message(bx_message::error) << "No trigger reference found" << dispatch;
      return ev;
    }
    trigger_ref_time = trigger_ref->get_lead_time();
  }*/

  double time, charge;
  for (int32_t i=0 ; i<er.get_raw_nhits(); i++) {    
    const bx_muon_raw_hit& h = er.get_raw_hit(i);

    if (profile_info.logical_channel_description(h.get_logical_channel()) != db_profile::ordinary) continue;

    if (b_discard_disabled_lg && !ev->get_trigger().is_pulser() && p_disabled_mch[h.get_muon_channel()]) continue;

    // decode time 
    if (!is_new_tdc) { // OLD tdc 
      time = (double)h.get_lead_time();  // pasitive, backward referred to common stop
      time -= laser_ref_time; // 0 for non laser events
      time += constants::muon::tdc::ns_per_clock; // convert to ns
      time = -time; //turn to negative
    }
    else { // NEW tdc
      time = (double)h.get_lead_time(); // positive, referred to beginning of gate
      time -= integrity_time[h.get_muon_channel() < constants::muon::tdc::channels_per_board]; // negative, referred to sync signal
      time -= laser_ref_time; // 0 for non laser events, negative referred to integrity of board 2 otherwise.
      time *= constants::muon::tdc::new_ns_per_clock; // convert to ns
    }



    // skip shit, non-laser events only      
    if (!ev->get_trigger().is_laser394() && (time > 0. || time < - time_buffer_ns)) continue;

    charge = h.get_time_diff() - run_info.get_muon_precalib_pedestal(h.get_logical_channel());
    if (charge < 0) {
      charge = 0.;
      //get_message(bx_message::debug) << "ev " << ev->get_event_number() << " mch " << h.get_muon_channel() << "; negative charge." << dispatch;
    }

    // apply corrections to time and charge for non-laser events
    if (b_apply_calibration && !ev->get_trigger().is_laser394()) {
      time -= run_info.get_muon_time_offset(h.get_logical_channel());
      charge /= run_info.get_muon_charge_peak(h.get_logical_channel());
    }

    bx_muon_decoded_hit decoded_hit(time, charge, &h);

    // Associate the db_channel
    decoded_hit.p_db_channel = &dynamic_cast<const db_channel_muon&>(bx_dbi::get ()->get_channel (h.get_logical_channel()));

    ew.decoded_hits.push_back(decoded_hit);
    ew.f4_decoded_charge+=charge;
    ew.nhits_per_channel_dec[h.get_muon_channel()]++;
  }  //end of loop over raw hits

  for (int32_t i=0 ; i< constants::muon::channels; i++) { if (ew.nhits_per_channel_dec[i]>0) ew.i4_decoded_npmts++; }

  er.mark_stage (bx_base_event::decoded);

  return ev;  
}

void bx_muon_decoder::end () {
}


/*
 * $Log: bx_muon_decoder.cc,v $
 * Revision 1.43  2012/04/27 16:02:56  ddangelo
 * patched to use integrity signals as trg refs on independent boards. minimally tested.
 *
 * Revision 1.42  2012-01-16 13:05:01  ddangelo
 * cleaned up and commented
 *
 * Revision 1.41  2012-01-16 09:48:31  ddangelo
 * decoding upgraded to handle new tdc format
 *
 * Revision 1.40  2009-10-23 09:27:30  wurm
 * commented debug messages
 *
 * Revision 1.39  2009-07-30 16:32:11  ddangelo
 * re-established negative charge correction
 * channel disabling not in pulser events
 * (commit already done in c11)
 *
 * Revision 1.38  2008-12-11 17:42:11  ddangelo
 * reversed the logic to set is_aligned flag
 *
 * Revision 1.37  2008-10-09 13:07:35  ddangelo
 * debug
 *
 * Revision 1.36  2008-10-09 09:26:20  ddangelo
 * apply calibration made after parameter
 *
 * Revision 1.35  2008-10-01 16:44:05  ddangelo
 * correction for negative charge hits removed (tmp to study soem problems)
 *
 * Revision 1.34  2008-10-01 09:45:17  wurm
 * removed pe_per_tick constant
 *
 * Revision 1.33  2008-09-29 13:31:58  wurm
 * charge and time offset corrections to decoded hits
 *
 * Revision 1.32  2008-09-04 15:31:44  ddangelo
 * checking muon alignment data presence before returning on misalaigned events
 *
 * Revision 1.31  2008-08-29 16:55:45  ddangelo
 * debugging
 *
 * Revision 1.30  2008-08-26 13:46:05  ddangelo
 * debug
 *
 * Revision 1.29  2008-08-26 13:39:53  ddangelo
 * writing is_aligned flag
 *
 * Revision 1.28  2008-08-20 16:20:23  ddangelo
 * disabling of channels introduced
 *
 * Revision 1.27  2007-12-20 18:40:16  ddangelo
 * fixed a bad bug
 *
 * Revision 1.26  2007-12-07 17:33:59  ddangelo
 * compliant with new dynamically loaded tdc clock range
 *
 * Revision 1.25  2007-11-29 15:01:01  ddangelo
 * debugging
 *
 * Revision 1.24  2007-05-25 15:09:55  ddangelo
 * minor things
 *
 * Revision 1.23  2007-05-07 13:40:49  ddangelo
 * applying patch to flag TObjects with cycle numbers
 *
 * Revision 1.22  2007-03-28 17:43:15  ddangelo
 * code restyled.
 * dec charge/npmts added.
 *
 * Revision 1.21  2007-02-23 19:22:18  ddangelo
 * added a casting to double
 *
 * Revision 1.20  2007/02/23 14:17:10  ddangelo
 * hit time expressed as negative with ref to trigger
 * random events processed too
 *
 * Revision 1.19  2007/02/22 15:50:44  ddangelo
 * tmp removed an error message.
 * fixed a few data types
 *
 * Revision 1.18  2006/09/08 12:40:00  ddangelo
 * db channel association done
 *
 * Revision 1.17  2006/01/25 13:12:00  misiaszek
 * Added a missing include (auth from maintainer)
 *
 * Revision 1.16  2004/11/27 16:41:48  ddangelo
 * adapted to look precalib pedestals with logical channel instead of muon channel
 *
 * Revision 1.15  2004/11/26 15:25:11  razeto
 * Thanks to Davide who so kindly underlined a syntactical error, I fixed the "Maintainer" word
 *
 * Revision 1.14  2004/11/26 14:06:23  razeto
 * Added Mantainer field
 *
 * Revision 1.13  2004/09/22 13:29:16  ddangelo
 * updated bx_reco_event to bx_echidna_event
 *
 * Revision 1.12  2004/09/22 12:13:52  ddangelo
 * updated a getter name
 *
 * Revision 1.11  2004/09/22 11:29:58  ddangelo
 * added require_event_stage in module's constructor
 *
 * Revision 1.10  2004/06/25 14:41:30  ddangelo
 * added require and mark event stage
 *
 * Revision 1.9  2004/06/07 12:48:57  ddangelo
 * special channels removed from decoded hits vector
 *
 * Revision 1.8  2004/06/02 17:15:00  ddangelo
 * fixed some stuff
 *
 * Revision 1.7  2004/06/01 11:36:27  ddangelo
 * updated to new constants syntax
 *
 * Revision 1.6  2004/05/20 11:02:30  ddangelo
 * trigger type selection introduced.
 * Job splitted in 2 private methods (physics/led)
 * 2nd one finds led ref pulse with a pre-loop.
 * chann descr asked to db_profile (tmp)
 *
 * Revision 1.5  2004/05/18 18:21:03  ddangelo
 * added an include
 *
 * Revision 1.4  2004/04/28 11:04:49  ddangelo
 * Modules adapted to call the right const/nonconst event portion getter
 *
 * Revision 1.3  2004/04/27 09:46:50  ddangelo
 * modifications to match new event structure.
 * Other debugging
 *
 * Revision 1.2  2004/04/06 12:42:17  razeto
 * Added module role to base_module, and the code to handle everywhere. Some other minor fixes
 *
 * Revision 1.1  2004/04/02 09:56:03  ddangelo
 * created bx_muon_decoder module.
 * precalibration and calibration/db info are still commented out
 *
 *
 */

