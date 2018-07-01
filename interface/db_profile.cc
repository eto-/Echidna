/* BOREXINO Reconstruction program
 *
 * Author: Alessandro Razeto <Alessandro.Razeto@ge.infn.it>
 * Maintainer: Alessandro Razeto <Alessandro.Razeto@ge.infn.it>
 *
 * $Id: db_profile.cc,v 1.26 2012/03/22 19:08:36 razeto Exp $
 *
 * The database interface for profile informations
 *
 */
#include "db_profile.hh"
#include "bx_dbi.hh"
#include "messenger.hh"
#include "messenger.hh"

#include <sstream>

ClassImp(db_profile)

#ifndef _ECHIDNA_ROOTLIB_
	
db_profile::db_profile (int profile_id) : db_acl(), TObject() {

  channel_description_type_map["Ordinary"] = db_profile::ordinary;
  channel_description_type_map["Trigger"] = db_profile::trigger;
  channel_description_type_map["Laser"] = db_profile::laser;
  
  channel_description_type_map["Pulser"] = db_profile::pulser;

  channel_description_type_map["Cngs-InnerDetector"] = db_profile::cngs_id;
  channel_description_type_map["Cngs-OuterDetector"] = db_profile::cngs_od;
  channel_description_type_map["Cngs-Trigger"] = db_profile::cngs_trg;
  
  // muon channel descriptions
  channel_description_type_map["integrity"] = db_profile::integrity;
  channel_description_type_map["ordinary"] = db_profile::ordinary;
  channel_description_type_map["trigger"] = db_profile::trigger;
  channel_description_type_map["laser"] = db_profile::laser;
  channel_description_type_map["omt"] = db_profile::omt;
  channel_description_type_map["clock"] = db_profile::clock;
  channel_description_type_map["empty"] = db_profile::empty;

  
  m_load_laben_channel_mapping (profile_id);  
  m_load_muon_channel_mapping (profile_id);  
  m_load_muon_general_setting (profile_id);  

  m_load_profile_table (profile_id);

  m_load_detector_geometry (profile_id);

  m_load_physics_constants (profile_id);
  if (profile_id <= 10) {
    laben_crate_delay_v[1]  = 4.8;
    laben_crate_delay_v[2]  = 4.7;
    laben_crate_delay_v[3]  = 3.7;
    laben_crate_delay_v[4]  = 3.7;
    laben_crate_delay_v[5]  = 1.9;
    laben_crate_delay_v[6]  = 2.7;
    laben_crate_delay_v[7]  = 0.;
    laben_crate_delay_v[8]  = 3.1;
    laben_crate_delay_v[9]  = 4.5;
    laben_crate_delay_v[10] = 4.5;
    laben_crate_delay_v[11] = -7.7;
    laben_crate_delay_v[12] = 2.4;
    laben_crate_delay_v[13] = 4.7;
    laben_crate_delay_v[14] = 4.;
  } else {
    laben_crate_delay_v[1]  = 0.;
    laben_crate_delay_v[2]  = 0.;
    laben_crate_delay_v[3]  = 0.;
    laben_crate_delay_v[4]  = 0.;
    laben_crate_delay_v[5]  = -2.48;
    laben_crate_delay_v[6]  = -1.85;
    laben_crate_delay_v[7]  = -4.53;
    laben_crate_delay_v[8]  = -1.41;
    laben_crate_delay_v[9]  = -0.20;
    laben_crate_delay_v[10] = -0.18 ;
    laben_crate_delay_v[11] = -12.65;
    laben_crate_delay_v[12] = -2.22;
    laben_crate_delay_v[13] = -0.09;
    laben_crate_delay_v[14] = -0.83;
  }
  bx_dbi::get ()->close_db_connections ();
  bx_dbi::get ()->get_message (bx_message::info) << "db_profile: db_profile (" << profile_id << ") initialized" << dispatch;
}

void db_profile::m_load_laben_channel_mapping (int profile_id) {
  bx_dbi *dbi = bx_dbi::get ();
  bx_message &msg = bx_dbi::get ()->get_message (bx_message::critic);
  msg << "db_profile: ";

  std::ostringstream t;
  t << "\"ProfileID\"=" << profile_id << " AND \"ChannelID\">0";

  const bx_dbi::table &table = dbi->query (bx_dbi::daq_config, "\"LabenChannelMapping\"", t.str (), "\"ChannelID\",\"ChannelDescription\"", "\"ChannelID\"", 0);

  if (!table.size ()) msg << "no LabenChannelMapping entries for profile " << profile_id << dispatch;

  const bx_dbi::column& ch_id = table["ChannelID"];
  const bx_dbi::column& ch_descr = table["ChannelDescription"];
  
  channel_description_type_map["Ordinary"] = db_profile::empty; // this is an hack since the db has ordinary even for
 								// channels without pmt. The channel will be set to ordinary
								// later while loading bx_geometry.
  for (size_t i = 0; i < ch_id.size (); i++) {
    if (!map_check(ch_descr[i].get_string (), channel_description_type_map)) msg << "unknown channel description " << ch_descr[i] << dispatch;
    else logical_channel_description_v[ch_id[i].get_int ()] = channel_description_type_map[ch_descr[i].get_string ()];
  }  
  channel_description_type_map["Ordinary"] = db_profile::ordinary; // restore the map


  const bx_dbi::table &table2 = dbi->query (bx_dbi::bx_geometry, "\"HolesMapping\"", t.str (), "*", "\"ChannelID\"", 0);

  if (!table2.size ()) msg << "no HolesMapping entries for profile " << profile_id << dispatch;
  const bx_dbi::column& ch_id2 = table2["ChannelID"];
  const bx_dbi::column& hole = table2["HoleLabel"];
  const bx_dbi::column& conc = table2["Conc"];
  const bx_dbi::column& x = table2["X"];
  const bx_dbi::column& y = table2["Y"];
  const bx_dbi::column& z = table2["Z"];
  const bx_dbi::column& theta = table2["Theta"];
  const bx_dbi::column& phi = table2["Phi"];
  const bx_dbi::column& fiber_bundle = table2["FiberBundle"];

  i_number_of_internal_pmt = 0;
  for (size_t i = 0; i < ch_id2.size (); i++) {
    int lg = ch_id2[i].get_int ();
    pmt_has_cone_v[lg] = conc[i].get_int ();
    pmt_hole_id_v[lg] = hole[i].get_int ();
    pmt_coordinates c;
    c.x = x[i].get_float ();
    c.y = y[i].get_float ();
    c.z = z[i].get_float ();
    c.theta = theta[i].get_float ();
    c.phi = phi[i].get_float ();
    pmt_coordinates_v[lg] = c;
    pmt_fiber_bundle_v[lg] = fiber_bundle[i].get_int ();
    if (logical_channel_description_v[lg] == db_profile::empty) {
      logical_channel_description_v[lg] = db_profile::ordinary; // this ends the hack described above
      i_number_of_internal_pmt ++;
    }
  }

  dbi->get_message (bx_message::debug) << "db_profile: load_laben_channel_mapping OK" << dispatch;
}


void db_profile::m_load_muon_general_setting (int profile_id) {
  bx_dbi *dbi = bx_dbi::get ();

  bx_message &msg = bx_dbi::get ()->get_message (bx_message::critic);
  msg << "db_profile: ";

  std::ostringstream t;
  t << "\"ProfileID\"=" << profile_id;

  const bx_dbi::table &table = dbi->query (bx_dbi::daq_config, "\"MuonGeneralSetting\"", t.str (), "\"TDCClockRange\"", "", 1);

  if (!table.size ()) msg << "no MuonGeneralSetting entry for profile " << profile_id << dispatch;

  const bx_dbi::column& range = table["TDCClockRange"];
  
  i4_muon_tdc_range = range[0].get_int ();

  dbi->get_message (bx_message::debug) << "db_profile: load_muon_general_setting OK" << dispatch;
}

void db_profile::m_load_muon_channel_mapping (int profile_id) {
  bx_dbi *dbi = bx_dbi::get ();

  bx_message &msg = bx_dbi::get ()->get_message (bx_message::critic);
  msg << "db_profile: ";

  std::ostringstream t;
  t << "\"ProfileID\"=" << profile_id << " AND \"ChannelID\">0";

  const bx_dbi::table &table = dbi->query (bx_dbi::daq_config, "\"MuonChannelMapping\"", t.str (), "\"ChannelID\",\"ChannelDescription\"", "\"ChannelID\"", 0);

  if (!table.size ()) msg << "no MuonChannelMapping entries for profile " << profile_id << dispatch;

  const bx_dbi::column& ch_id = table["ChannelID"];
  const bx_dbi::column& ch_descr = table["ChannelDescription"];
  
  for (size_t i = 0; i < ch_id.size (); i++) {
    if (!map_check(ch_descr[i].get_string (), channel_description_type_map)) msg << "unknown channel description " << ch_descr[i] << dispatch;
    else logical_channel_description_v[ch_id[i].get_int ()] = channel_description_type_map[ch_descr[i].get_string ()];
  }  

  const bx_dbi::table &table2 = dbi->query (bx_dbi::bx_geometry, "\"MuonHolesMapping\"", t.str (), "\"ChannelID\",\"HoleLabel\",\"X\",\"Y\",\"Z\"", "\"ChannelID\"", 0);

  if (!table2.size ()) msg << "no MuonHolesMapping entries for profile " << profile_id << dispatch;
  const bx_dbi::column& ch_id2 = table2["ChannelID"];
  const bx_dbi::column& hole = table2["HoleLabel"];
  const bx_dbi::column& x = table2["X"];
  const bx_dbi::column& y = table2["Y"];
  const bx_dbi::column& z = table2["Z"];
  for (size_t i = 0; i < ch_id2.size (); i++) {
    int lg = ch_id2[i].get_int ();
    pmt_hole_id_v[lg] = hole[i].get_int ();
    pmt_coordinates c;
    c.x = x[i].get_float ();
    c.y = y[i].get_float ();
    c.z = z[i].get_float ();
    c.phi = c.theta = 0;
    pmt_coordinates_v[lg] = c;
  }
  i_number_of_external_pmt = ch_id2.size ();

  dbi->get_message (bx_message::debug) << "db_profile: load_muon_channel_mapping OK" << dispatch;
}

void db_profile::m_load_detector_geometry (int profile_id) {
  bx_dbi *dbi = bx_dbi::get ();
  bx_message &msg = dbi->get_message (bx_message::critic);
  msg << "db_profile: ";

  std::ostringstream t;
  t << "\"ProfileID\"=" << profile_id;
  const bx_dbi::table &detector_params = dbi->query (bx_dbi::bx_geometry, "\"DetectorGeometry\"", t.str (), "*", "", -1);

  f_stainless_steel_sphere_radius = detector_params["SSSRadius"][0].get_float ();
  f_inner_vessel_radius = detector_params["InnerVesselRadius"][0].get_float ();
  f_light_collector_radius = f_stainless_steel_sphere_radius - detector_params["LGEntranceOffset"][0].get_float ();
  f_light_guide_exit_aperture = detector_params["LGExitRadius"][0].get_float ();
  f_light_guide_entry_aperture = detector_params["LGEntryRadius"][0].get_float ();
  f_light_guide_lenght = detector_params["LGLength"][0].get_float ();
  f_pmt_chathode_radius = detector_params["CathodeRadius"][0].get_float ();
  f_pmt_bulb_radius = detector_params["BulbRadius"][0].get_float ();

  dbi->get_message (bx_message::debug) << "db_profile: m_load_detector_geometry OK" << dispatch;
}


void db_profile::m_load_profile_table (int profile_id) {
  bx_dbi *dbi = bx_dbi::get ();
  bx_message &msg = dbi->get_message (bx_message::critic);
  msg << "db_profile: ";

  std::ostringstream t;
  t << "\"ProfileID\"=" << profile_id;
  const bx_dbi::table &profile_params = dbi->query (bx_dbi::daq_config, "\"Profiles\"", t.str (), "*", "", -1);

  i4_neutrino_trigger_tag = profile_params["NeutrinoTrigger"][0].get_int ();
  i4_muon_trigger_tag = profile_params["MuonMTBTrigger"][0].get_int ();
  i4_neutron_trigger_tag = profile_params["MuonTCTTrigger"][0].get_int ();
  i4_laser266_trigger_tag = profile_params["Laser266Trigger"][0].get_int ();
  i4_laser355_trigger_tag = profile_params["Laser355Trigger"][0].get_int ();
  i4_laser394_trigger_tag = profile_params["Laser394Trigger"][0].get_int ();
  i4_pulser_trigger_tag = profile_params["CalibPulserTrigger"][0].get_int ();
  i4_random_trigger_tag = profile_params["RandomTrigger"][0].get_int ();
  f_trigger_start_gate = profile_params["StartGateTrigger"][0].get_float ();
  f_trigger_end_gate = profile_params["EndGateTrigger"][0].get_float ();
  f_laser394_trigger_start_gate = profile_params["StartGateLaser"][0].get_float ();
  f_laser394_trigger_end_gate = profile_params["EndGateLaser"][0].get_float ();

  dbi->get_message (bx_message::debug) << "db_profile: load_profile_table OK" << dispatch;
}

void db_profile::m_load_physics_constants (int profile_id) {
  bx_dbi *dbi = bx_dbi::get ();

  std::ostringstream t;
  t << "\"ProfileID\"=" << profile_id;
  const bx_dbi::table &physics_params = dbi->query (bx_dbi::bx_physics, "\"FixedParameters\"", t.str (), "*", "", -1);

  f_pmt_time_jitter = physics_params["PMTTimeJitter"][0].get_float ();
  f_pmt_reflectivity = physics_params["PMTReflCoeff"][0].get_float ();
  f_light_guide_reflectivity = physics_params["LGExterCoeff"][0].get_float ();
  f_sss_specular_reflectivity = physics_params["SSSSpecRefl"][0].get_float ();
  f_sss_diffusion_reflectivity = physics_params["SSSDiffRefl"][0].get_float ();
  f_sss_uni_reflectivity = physics_params["SSSUnifRefl"][0].get_float ();

  dbi->get_message (bx_message::debug) << "db_profile: m_load_physics_constants OK" << dispatch;
}

#endif

/*  
 *  $Log: db_profile.cc,v $
 *  Revision 1.26  2012/03/22 19:08:36  razeto
 *  Added cngs reference channels
 *
 *  Revision 1.25  2011-02-18 17:10:05  ddangelo
 *  major code cleanup: removed fadc throughout the program
 *
 *  Revision 1.24  2009-10-26 11:19:48  ddangelo
 *  trgtype 'muon_tct' (or 'muon_totc') renamed as 'neutron' after real use
 *  trgtype 'muon_mtb' renamed as 'muon' as duplicity with the above one disappeared and 'mtb' is associated to btb_inputs flag
 *
 *  Revision 1.23  2007-12-10 09:40:44  ludhova
 *  buggy name of the DB column fixed
 *
 *  Revision 1.22  2007-12-07 17:33:02  ddangelo
 *  added loading of muon_general_stting table
 *
 *  Revision 1.21  2007-06-03 16:13:30  razeto
 *  Do not keep db connection opened (else bxdb has lots of problems)
 *
 *  Revision 1.20  2007-05-07 10:13:08  razeto
 *  Count correctly number of pmts
 *
 *  Revision 1.19  2006/07/20 14:53:34  ludhova
 *  Added FadcChannelMapping
 *
 *  Revision 1.18  2006/01/25 13:23:54  misiaszek
 *  Moved from cmap to simple map (to work with root)
 *
 *  Revision 1.17  2005/12/11 16:36:07  razeto
 *  Fixed a warning
 *
 *  Revision 1.16  2004/11/26 15:25:10  razeto
 *  Thanks to Davide who so kindly underlined a syntactical error, I fixed the "Maintainer" word
 *
 *  Revision 1.15  2004/11/26 14:01:20  razeto
 *  Added Mantainer field
 *
 *  Revision 1.14  2004/11/24 13:06:03  razeto
 *  Upgraded to the new cmap ctor with name
 *
 *  Revision 1.13  2004/09/27 13:31:06  razeto
 *  Fibre bundle, theta and phi read from database
 *
 *  Revision 1.12  2004/08/10 11:24:33  razeto
 *  Fixed some constants for crate delays
 *
 *  Revision 1.11  2004/08/06 10:30:28  razeto
 *  cycle_1 branch merged in the main trunk.
 *
 *  Revision 1.10  2004/08/05 09:05:20  razeto
 *  Used postgres EXTRACT command to convert date to time_t
 *
 *  Revision 1.9.2.2  2004/07/29 08:54:29  razeto
 *  A dirty patch to support new crate delays
 *
 *  Revision 1.9.2.1  2004/07/26 18:04:54  razeto
 *  Fixed a bug in reading the channel mapping
 *
 *  Revision 1.9  2004/06/10 10:40:01  razeto
 *  Fixed some bugs
 *
 *  Revision 1.8  2004/06/10 10:07:42  razeto
 *  Added FADC handling
 *
 *  Revision 1.7  2004/05/26 10:32:13  razeto
 *  Fixed a bug (using the wrong table name)
 *
 *  Revision 1.6  2004/05/26 09:21:59  razeto
 *  Upgraded to last version of vdt
 *
 *  Revision 1.5  2004/05/25 17:15:47  razeto
 *  Upgraded to use name in vdt::get_*
 *
 *  Revision 1.4  2004/05/20 10:49:08  razeto
 *  Applied Davide's patch for loading muon tables
 *
 *  Revision 1.3  2004/04/26 13:48:29  razeto
 *  Added db_acl to check set calls for calling module to have the right privileges.
 *  Fixed the names of set/get methods.
 *  Modifications decided at the software Paris meeting.
 *
 *  Revision 1.2  2004/04/12 16:08:19  razeto
 *  Added crate delay fixed map
 *
 *  Revision 1.1  2004/04/09 08:05:00  razeto
 *  Added
 *
 */
