/* BOREXINO Reconstruction program
 *
 * Author: Alessandro Razeto <Alessandro.Razeto@ge.infn.it>
 * Maintainer: Alessandro Razeto <Alessandro.Razeto@ge.infn.it>
 *
 * $Id: db_channel.cc,v 1.22 2011/02/18 17:10:05 ddangelo Exp $
 *
 * The database interface for channel informations. READONLY
 *
 */
#include "db_channel.hh"
#include "constants.hh"
#include <algorithm>
#include "db_run.hh"

ClassImp(db_channel)
ClassImp(db_channel_laben)
ClassImp(db_channel_muon)

#ifndef _ECHIDNA_ROOTLIB_
db_channel::db_channel (int lg, const db_run& run, const db_profile& profile): TObject(), i4_lg(lg) {
  channel_description = profile.logical_channel_description (lg);
}

db_channel_laben::db_channel_laben (int lg, const db_run& run, const db_profile& profile): db_channel (lg, run, profile) {
  if (is_ordinary ()) {
    b_pmt_has_cone = profile.pmt_has_cone (lg);
    i4_pmt_hole_id = profile.pmt_hole_id (lg);
    const db_profile::pmt_coordinates &c = profile.get_pmt_coordinates (lg);
    f_x = c.x;
    f_y = c.y;
    f_z = c.z;
    f_theta = c.theta;
    f_phi = c.phi;
    i4_pmt_fiber_bundle = profile.pmt_fiber_bundle (lg);
    b_disconnected = run.is_pmt_disconnected (lg);
  }

  const std::vector<int>& bad_ch = run.get_laben_precalib_bad_channels ();
  if (std::find (bad_ch.begin (), bad_ch.end (), lg) != bad_ch.end ()) b_precalib_bad = true;
  else b_precalib_bad = false;

  const std::vector<int>& off_ch = run.get_laben_precalib_off_channels ();
  if (std::find (off_ch.begin (), off_ch.end (), lg) != off_ch.end ()) b_precalib_off = true;
  else b_precalib_off = false;

  f_time_offset = run.get_laben_time_offset (lg);
  f_time_sigma = run.get_laben_time_sigma (lg);
  f_charge_peak = run.get_laben_charge_peak (lg);
  f_charge_sigma = run.get_laben_charge_sigma (lg);
  f_dark_noise = run.get_laben_dark_noise (lg);
  f_dark_sigma = run.get_laben_dark_sigma (lg);
  s_pmt_status = run.get_laben_pmt_status (lg);

  v_charge_base_status = run.get_laben_charge_base_status (lg);
  v_charge_peak_status = run.get_laben_charge_peak_status (lg);
  v_charge_status = run.get_laben_charge_status (lg);
  v_timing_status = run.get_laben_timing_status (lg);
  v_multiplicity = run.get_laben_multiplicity (lg);
}

db_channel_muon::db_channel_muon (int lg, const db_run& run, const db_profile& profile): db_channel (lg, run, profile) {
  if (is_ordinary ()) {
    i4_hole_id = profile.pmt_hole_id (lg);
    const db_profile::pmt_coordinates &c = profile.get_pmt_coordinates (lg);
    f4_x = c.x;
    f4_y = c.y;
    f4_z = c.z;
    b_disconnected = run.is_pmt_disconnected (lg);
  }

  f4_time_offset  = run.get_muon_time_offset  (lg);
  f4_time_sigma   = run.get_muon_time_sigma   (lg);
  f4_charge_peak  = run.get_muon_charge_peak  (lg);
  f4_charge_sigma = run.get_muon_charge_sigma (lg);
  f4_dark_noise   = run.get_muon_dark_noise   (lg);
  f4_dark_sigma   = run.get_muon_dark_sigma   (lg);
  s_pmt_status    = run.get_muon_pmt_status   (lg);
  v_multiplicity  = run.get_muon_multiplicity (lg);
}

db_channel* db_channel_builder::build (int lg, const db_run& run, const db_profile& profile) { 
  if (constants::laben::is_laben (lg)) return new db_channel_laben (lg, run, profile);
  if (constants::muon::is_muon (lg)) return new db_channel_muon (lg, run, profile);
  return 0;
}
#endif

/*  
 *  $Log: db_channel.cc,v $
 *  Revision 1.22  2011/02/18 17:10:05  ddangelo
 *  major code cleanup: removed fadc throughout the program
 *
 *  Revision 1.21  2008-08-22 17:16:25  ddangelo
 *  fixed the bug
 *
 *  Revision 1.20  2008-08-21 17:12:57  ddangelo
 *  tmp commit.
 *  commenting muon db channel filling to avid crash after a bug in led calibrations
 *
 *  Revision 1.19  2008-08-20 14:03:42  ddangelo
 *  muon db channel completed
 *
 *  Revision 1.18  2008-08-19 17:52:18  ddangelo
 *  debugging
 *
 *  Revision 1.17  2008-08-19 17:51:31  ddangelo
 *  added muon multiplicity and pmt disconnected
 *
 *  Revision 1.16  2007-11-09 18:49:38  razeto
 *  Disconnected added to laben_channel
 *
 *  Revision 1.15  2007-03-27 11:00:09  razeto
 *  Fixed a *conceptual* bug
 *
 *  Revision 1.14  2007-03-24 15:59:12  ddangelo
 *  muon db channel improved
 *
 *  Revision 1.13  2007-03-23 18:56:56  ddangelo
 *  loading of fadc channel tmp commented out:
 *  it causes run time exception.
 *
 *  Revision 1.12  2007-03-22 19:36:58  ddangelo
 *  implemented muon db channel
 *
 *  Revision 1.11  2006/09/09 18:51:32  razeto
 *  Fixed the default value for few fields
 *
 *  Revision 1.10  2006/09/08 10:11:47  razeto
 *  Upgraded db_channels and added muon db_channel
 *
 *  Revision 1.9  2006/05/08 17:31:33  razeto
 *  Added db_channel patch for fadc (sent to the mailing list)
 *
 *  Revision 1.8  2006/01/25 13:24:34  misiaszek
 *  Added ifndef _ECHIDNA_ROOTLIB_
 *
 *  Revision 1.7  2005/11/04 00:58:41  misiaszek
 *
 *  ----------------------------------------------------------------------
 *   Modified Files:
 *   	db_channel.cc
 *
 *  TObject ctor called.
 *  ----------------------------------------------------------------------
 *
 *  Revision 1.6  2004/11/26 15:25:10  razeto
 *  Thanks to Davide who so kindly underlined a syntactical error, I fixed the "Maintainer" word
 *
 *  Revision 1.5  2004/11/26 14:01:20  razeto
 *  Added Mantainer field
 *
 *  Revision 1.4  2004/09/27 13:31:06  razeto
 *  Fibre bundle, theta and phi read from database
 *
 *  Revision 1.3  2004/06/01 11:51:48  razeto
 *  Temporary removed a FADC copy since FADC data are not read in db_profile
 *
 *  Revision 1.2  2004/05/19 11:06:27  razeto
 *  Updated to a working version of db_channel
 *
 *  Revision 1.1  2004/05/19 10:25:50  razeto
 *  Added db_channel as discussed at last cycle meeting
 *
 */
