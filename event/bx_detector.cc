/* BOREXINO Reconstruction program
 *
 * Author: Alessandro Razeto <Alessandro.Razeto@ge.infn.it>
 * Maintainer: Alessandro Razeto <Alessandro.Razeto@ge.infn.it>
 *
 * $Id: bx_detector.cc,v 1.42 2013/01/25 16:30:13 razeto Exp $
 *
 * bx_detector implementation
 *
 */
#include "bx_detector.hh"
ClassImp(bx_detector);

#ifndef _ECHIDNA_ROOTLIB_
#include <algorithm>
#include "bx_base_module.hh"
#include "constants.hh"
#include "bx_dbi.hh"
#include "db_run.hh"
#include "db_channel.hh"
#include "bx_echidna_event.hh"

detector_interface * detector_interface::me = 0;

detector_interface::detector_interface (): bx_named ("detector_interface") {
  detector = new bx_detector;
  detector->total_skipped_events = 0;
  detector->last_evnum_read = detector->next_evnum_to_read = 0;
  detector->use_db = false;

  detector->b_laben_enabled   = get_parameter ("laben")  .get_bool ();
  detector->b_muon_enabled    = get_parameter ("muon")   .get_bool ();
  detector->b_mctruth_enabled = get_parameter ("mctruth").get_bool ();
  detector->scaling_factor    = get_parameter ("scaling_factor").get_float ();

  const vdt::vdt_vector& dlg = get_parameter("disable_lg").get_vector ();
  for (unsigned i = 0; i < dlg.size (); i++) detector->user_disabled_channels_v.push_back(dlg[i].get_int ());
  const vdt::vdt_vector& blcp = get_parameter("disabled_laben_channel_properties").get_vector ();
  for (unsigned i = 0; i < blcp.size (); i++) detector->disabled_laben_channel_properties_v.push_back(blcp[i].get_string ());
  const vdt::vdt_vector& bmcp = get_parameter("disabled_muon_channel_properties").get_vector ();
  for (unsigned i = 0; i < bmcp.size (); i++) detector->disabled_muon_channel_properties_v.push_back(bmcp[i].get_string ());

}

void detector_interface::set_event_detector_status (const bx_echidna_event &e) {
  detector->compute_disabled_channels (e.get_run_number ());

  detector->b_laben_enabled   &= e.is_laben_enabled   ();
  detector->b_muon_enabled    &= e.is_muon_enabled    ();
  detector->b_mctruth_enabled &= e.is_mctruth_enabled ();

  bx_message &msg =  get_message (bx_message::info);
  const char *status[] = { "off", "on" };
  msg << newline << "laben detector is " << status[int(detector->b_laben_enabled)] << newline;
  msg << "muon detector is " << status[int(detector->b_muon_enabled)] << newline;
}

void bx_detector::skip_event (bx_base_module *module, int trgtype) {
  skipped_events[trgtype][module->get_name ()]++;
  total_skipped_events++;
}

void bx_detector::disable_lg (int lg, const std::string& reason) {
  bx_message msg (bx_message::log, "bx_detector: ");
  msg << "channel " << lg << " because of " << reason << dispatch;
  disabled_channels_v.push_back (lg);
}

void bx_detector::disable_charge (int lg, const std::string& reason) {
  bx_message msg (bx_message::log, "bx_detector: ");
  msg << "charge channel " << lg << " because of " << reason << dispatch;
  disabled_charge_v.push_back (lg);
}

void bx_detector::compute_disabled_channels (int run_number) {
  disabled_channels_v.clear ();

    // First disable user specified channels
  if (check_laben_property ("user")) 
    for (uint32_t i = 0; i < user_disabled_channels_v.size (); i++)
      disable_lg (user_disabled_channels_v[i], "user");

  if (check_laben_property ("db")) use_db = true;

    // Then check each channel quality
  for (int i = 1; i <= constants::laben::channels; i++) {
    const db_channel_laben& ch_info = dynamic_cast<const db_channel_laben&>(bx_dbi::get ()->get_channel (i));
    if (check_laben_property ("empty") && ch_info.is_empty ())
      disable_lg (i, "empty");
    else if (check_laben_property ("disconnected") && ch_info.is_pmt_disconnected () && ch_info.is_ordinary ())
      disable_lg (i, "disconnected");
    else if (check_laben_property ("bad_precalib") && (ch_info.is_precalib_bad () || ch_info.is_precalib_off ()))
      disable_lg (i, "bad_precalib");
    else if (check_laben_property ("bad_pmt_timing") && (::fabs(ch_info.time_offset ()) > 20. || ch_info.time_sigma () > 10))
      disable_lg (i, "bad_pmt_timing");
    else if (check_laben_property ("dead_pmt") && ch_info.pmt_status () == "dead")
      disable_lg (i, "dead_pmt");
    else if (check_laben_property ("bad_elec_timing") && ch_info.timing_status ().size ())
      disable_lg (i, "bad_elec_timing");
    else if (check_laben_property ("hot_pmt") && ch_info.dark_noise () > 3 )
      disable_lg (i, "hot_pmt");
    else if (check_laben_property ("hot_in_neutrino", ch_info.multiplicity ()) && ch_info.dark_noise () < 2 && !check_vector_property ("retriggering_in_neutrino", ch_info.multiplicity ()))
      disable_lg (i, "hot_in_neutrino");
    else if (check_laben_property ("hot_in_laser", ch_info.multiplicity ()))
      disable_lg (i, "hot_in_laser");
    else if (check_laben_property ("dead_in_neutrino", ch_info.multiplicity ()))
      disable_lg (i, "dead_in_neutrino");
    else if (check_laben_property ("dead_in_raw", ch_info.multiplicity ()))
      disable_lg (i, "dead_in_raw");
    else if (check_laben_property ("low_eff_in_neutrino", ch_info.multiplicity ()))
      disable_lg (i, "low_eff_in_neutrino");
    else if (check_laben_property ("low_gain", ch_info.charge_status ()))
      disable_lg (i, "low_gain");
    else if (check_laben_property ("high_gain", ch_info.charge_status ()))
      disable_charge (i, "high_gain");
  }
  
  for (int i = 1+constants::muon::channel_offset; i <= constants::muon::channel_offset+constants::muon::channels; i++) {
    const db_channel_muon& ch_info = dynamic_cast<const db_channel_muon&>(bx_dbi::get ()->get_channel (i));
    if (check_muon_property ("empty") && ch_info.is_empty ())
      disable_lg (i, "empty");
    else if (check_muon_property ("disconnected") && ch_info.is_disconnected () && ch_info.is_ordinary ())
      disable_lg (i, "disconnected");
//    else if (check_muon_property ("bad_precalib") && (ch_info.is_precalib_bad () || ch_info.is_precalib_off ()))
//      disable_lg (i, "bad_precalib");
//    else if (check_muon_property ("bad_pmt_timing") && (::fabs(ch_info.time_offset ()) > 20. || ch_info.time_sigma () > 10))
//      disable_lg (i, "bad_pmt_timing");
//    else if (check_muon_property ("dead_pmt") && ch_info.pmt_status () == "dead")
//      disable_lg (i, "dead_pmt");
//    else if (check_muon_property ("bad_elec_timing") && ch_info.timing_status ().size ())
//      disable_lg (i, "bad_elec_timing");
    else if (check_muon_property ("hot_in_muon", ch_info.multiplicity ()) && !check_vector_property ("retriggering_in_muon", ch_info.multiplicity ())) //check if the last condition applies 
      disable_lg (i, "hot_in_muon");
    else if (check_muon_property ("hot_in_neutrino", ch_info.multiplicity ()) && !check_vector_property ("retriggering_in_neutrino", ch_info.multiplicity ())) //check if the last condition applies 
      disable_lg (i, "hot_in_neutrino");
    else if (check_muon_property ("hot_in_laser", ch_info.multiplicity ()))
      disable_lg (i, "hot_in_laser");
    else if (check_muon_property ("dead_in_neutrino", ch_info.multiplicity ()))
      disable_lg (i, "dead_in_neutrino");
    else if (check_muon_property ("low_eff_in_neutrino", ch_info.multiplicity ()))
      disable_lg (i, "low_eff_in_neutrino");
    else if (check_muon_property ("dead_in_muon", ch_info.multiplicity ()))
      disable_lg (i, "dead_in_muon");
    else if (check_muon_property ("low_eff_in_muon", ch_info.multiplicity ()))
      disable_lg (i, "low_eff_in_muon");
  }
    // Strips copies
  std::sort (disabled_channels_v.begin (), disabled_channels_v.end ());
  std::vector<int>::iterator new_end = std::unique (disabled_channels_v.begin (), disabled_channels_v.end ());
  disabled_channels_v.erase (new_end, disabled_channels_v.end ());
  bx_message msg (bx_message::info, "bx_detector: ");
  msg << "found " << disabled_channels_v.size () << " bad channels" << dispatch;
}

bool bx_detector::check_laben_property (const std::string& p) {
  return std::find (disabled_laben_channel_properties_v.begin (), disabled_laben_channel_properties_v.end (), p) != disabled_laben_channel_properties_v.end ();
}
bool bx_detector::check_muon_property (const std::string& p) {
  return std::find (disabled_muon_channel_properties_v.begin (), disabled_muon_channel_properties_v.end (), p) != disabled_muon_channel_properties_v.end ();
}
bool bx_detector::check_vector_property (const std::string& p, const std::vector<std::string> &v) {
  return std::find (v.begin (), v.end (), p) != v.end ();
}

bool bx_detector::check_laben_property (const std::string& p, const std::vector<std::string> &v) {
  if (std::find (disabled_laben_channel_properties_v.begin (), disabled_laben_channel_properties_v.end (), p) == disabled_laben_channel_properties_v.end ()) return false;
  return std::find (v.begin (), v.end (), p) != v.end ();
}
bool bx_detector::check_muon_property (const std::string& p, const std::vector<std::string> &v) {
  if (std::find (disabled_muon_channel_properties_v.begin (), disabled_muon_channel_properties_v.end (), p) == disabled_muon_channel_properties_v.end ()) return false;
  return std::find (v.begin (), v.end (), p) != v.end ();
}


void bx_detector::add_disabled_channel (int lg, int evnum, int type, const bx_named* obj) {
  if (obj->get_name () != "bx_detector_monitor") {
    bx_message msg (bx_message::error, "bx_detector: ");
    msg << obj->get_name () << " is not authorized to call add_disabled_channel" << dispatch;
  } else if (type != db_run::timing && type != db_run::charge) {
    bx_message msg (bx_message::error, "bx_detector: ");
    msg << obj->get_name () << " type is not timing or charge" << dispatch;
  } else if (!bx_dbi::get ()->get_run ().is_disabled_channels_present () || !use_db) {
    bx_dbi::get ()->get_run ().add_disabled_channel (lg, evnum, db_run::disabled_type(type), obj);
    if (type == db_run::timing && std::find (disabled_channels_v.begin (), disabled_channels_v.end (), lg) == disabled_channels_v.end ()) {
      disabled_channels_v.push_back (lg);
      std::sort (disabled_channels_v.begin (), disabled_channels_v.end ());
      std::vector<int>::iterator new_end = std::unique (disabled_channels_v.begin (), disabled_channels_v.end ());
      disabled_channels_v.erase (new_end, disabled_channels_v.end ());
    } else if (type == db_run::charge && std::find (disabled_charge_v.begin (), disabled_charge_v.end (), lg) == disabled_charge_v.end ()) {
      disabled_charge_v.push_back (lg);
      std::sort (disabled_charge_v.begin (), disabled_charge_v.end ());
      std::vector<int>::iterator new_end = std::unique (disabled_charge_v.begin (), disabled_charge_v.end ());
      disabled_charge_v.erase (new_end, disabled_charge_v.end ());
    }
  }
}

void bx_detector::read_disabled_channels (int evnum) {
  if (!use_db) return;
  evnum = int(evnum / scaling_factor);
  if (evnum < next_evnum_to_read) return;

  const std::map<int, std::map<int, db_run::disabled_type> >& dc = bx_dbi::get ()->get_run ().get_disabled_channels ();
  for (std::map<int, std::map<int, db_run::disabled_type> >::const_iterator it = dc.begin (); it != dc.end (); it++) {
    int ev = it->first;
    if (ev > evnum) {
      next_evnum_to_read = ev;
      break;
    }
    if (ev > last_evnum_read || (last_evnum_read == 0 && ev == 0)) {
      last_evnum_read = ev;
      for (std::map<int, db_run::disabled_type>::const_iterator jt = it->second.begin (); jt != it->second.end (); jt++) {
	if (jt->second == db_run::timing) disabled_channels_v.push_back (jt->first);
	else disabled_charge_v.push_back (jt->first);
      }
    }
  }
  std::sort (disabled_channels_v.begin (), disabled_channels_v.end ());
  std::vector<int>::iterator new_end = std::unique (disabled_channels_v.begin (), disabled_channels_v.end ());
  disabled_channels_v.erase (new_end, disabled_channels_v.end ());
  std::sort (disabled_charge_v.begin (), disabled_charge_v.end ());
  new_end = std::unique (disabled_charge_v.begin (), disabled_charge_v.end ());
  disabled_charge_v.erase (new_end, disabled_charge_v.end ());
}

#endif
/*
 * $Log: bx_detector.cc,v $
 * Revision 1.42  2013/01/25 16:30:13  razeto
 * do not use db unless requested
 *
 * Revision 1.41  2012-10-09 15:15:02  ludhova
 * hot_pmt option for disabling hot channels with DR > 3 kHz
 *
 * Revision 1.40  2011-03-01 19:21:20  razeto
 * Faster method for get_channel
 *
 * Revision 1.39  2011-03-01 18:12:17  razeto
 * low gain pmt are disabled in timing
 *
 * Revision 1.38  2011-03-01 12:56:21  razeto
 * Low/high gain log separated
 *
 * Revision 1.37  2011-02-28 12:33:42  razeto
 * Added bad_gain disabling
 *
 * Revision 1.36  2011-02-18 17:10:05  ddangelo
 * major code cleanup: removed fadc throughout the program
 *
 * Revision 1.35  2008-11-20 09:22:57  ludhova
 * change on limits on dark noise for hot_in_neutrino disabling
 *
 * Revision 1.34  2008-09-25 13:36:19  ludhova
 * disable dead_in_raw
 *
 * Revision 1.33  2008-09-25 13:12:35  razeto
 * Scaling factor added for disabling channels (intended for mc)
 *
 * Revision 1.32  2008-09-25 11:52:17  razeto
 * Added db property for bx_detector::disabled_laben_channel_properties
 *
 * Revision 1.31  2008-09-24 17:14:16  razeto
 * Added reading runtime disabled channels from db
 *
 * Revision 1.30  2008-08-20 14:03:19  ddangelo
 * finilizing with muon db channel
 *
 * Revision 1.29  2008-08-20 12:45:10  ddangelo
 * fixing old and new bugs...
 *
 * Revision 1.28  2008-08-19 18:12:03  ddangelo
 * debugging
 *
 * Revision 1.27  2008-08-19 17:53:00  ddangelo
 * added muon channel disabling
 *
 * Revision 1.26  2008-04-10 09:47:13  ludhova
 * check first disconnected, then bad precalib
 *
 * Revision 1.25  2007-12-15 12:15:44  ludhova
 * disabled pmts hot_in_neutrino, dark_rate < 3kHz, not retriggering_in_neutrino
 *
 * Revision 1.24  2007-12-14 19:13:33  ludhova
 * disabled pmts hot_in_neutrino and dark < 3 kHz
 *
 * Revision 1.23  2007-12-10 15:59:27  ludhova
 * enlarge limit for bad_pmt_timing
 *
 * Revision 1.22  2007-12-07 14:26:22  ddangelo
 * added vector for bad charge channels, getter, initialization and add() method (auth by maintainer)
 *
 * Revision 1.21  2007-11-29 17:17:59  razeto
 * Do not disable reference channels
 *
 * Revision 1.20  2007-11-27 15:08:30  razeto
 * Log the disable reason
 *
 * Revision 1.19  2007-11-20 12:26:53  razeto
 * discard empty and disconnected channels
 *
 * Revision 1.18  2007-10-30 17:36:53  razeto
 * Removed a warn
 *
 * Revision 1.17  2007-05-25 13:43:44  razeto
 * Be quiter
 *
 * Revision 1.16  2007-05-07 10:49:35  razeto
 * Bug fixing
 *
 * Revision 1.15  2007-05-07 10:46:03  razeto
 * Using right unique algorithm
 *
 * Revision 1.14  2007-05-04 11:08:46  razeto
 * add_disabled_channel implemented
 *
 * Revision 1.13  2007-05-02 15:58:12  razeto
 * Upgraded bx_detector to support variable disabled channel list. Few names changed
 *
 * Revision 1.12  2007-03-10 15:18:22  ddangelo
 * Now reader init the bx_detector (even with the event detector status).
 *
 * Revision 1.11  2007-01-25 17:21:37  razeto
 * Hot in neutrino are even hot channels, do not disable hot channels with hot in neutrino request
 *
 * Revision 1.10  2006/11/27 10:59:01  razeto
 * Upgraded bx_detector to have more detailed check on bad channels multiplicity:
 * now bad_multiplicity do not exist any more but there are other fields.
 * Upgraded laben_decoder: the fix_bad_channel flag do not exists anymore. To
 * avoid bad channels skipping, set detector_interface.bad_laben_channel_properties
 * to {}.
 * Configurations updated.
 *
 * Revision 1.9  2006/09/10 14:41:22  razeto
 * Removed inline since creates some problems with gcc 2.95
 *
 * Revision 1.8  2006/09/09 18:59:30  razeto
 * Upgraded to handle configuration for bad channel
 *
 * Revision 1.7  2006/09/09 14:06:56  razeto
 * Upgraded to have a draft calculation of bad channels
 *
 * Revision 1.6  2006/08/21 11:11:59  razeto
 * Improved bx_detector and created detector_interface
 *
 * Revision 1.5  2004/11/26 15:25:10  razeto
 * Thanks to Davide who so kindly underlined a syntactical error, I fixed the "Maintainer" word
 *
 * Revision 1.4  2004/11/26 14:01:19  razeto
 * Added Mantainer field
 *
 * Revision 1.3  2004/10/19 16:23:40  razeto
 * Updated to use vdt::get_bool
 *
 * Revision 1.2  2004/09/22 11:18:52  razeto
 * Added mctruth to bx_detector::sub_detector
 *
 * Revision 1.1  2004/09/22 10:17:16  razeto
 * Added
 *
 */
