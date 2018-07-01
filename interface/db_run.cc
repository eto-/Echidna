/* BOREXINO Reconstruction program
 *
 * Author: Alessandro Razeto <Alessandro.Razeto@ge.infn.it>
 * Maintainers: Alessandro Razeto <Alessandro.Razeto@ge.infn.it>,
 *
 * $Id: db_run.cc,v 1.68 2015/07/24 12:37:27 ilia.drachnev Exp $
 *
 * The database interface for run informations
 *
 */
#include "db_run.hh"
#include "bx_dbi.hh"
#include "bx_detector.hh"
#include "messenger.hh"
#include "constants.hh"

#include <string>
#include <sstream>

ClassImp(db_run)

#ifndef _ECHIDNA_ROOTLIB_

  //namespace {
  void copy_vdt_vector (std::vector<std::string>& destination, const vdt::vdt_vector& source) {
    for (unsigned i = 0; i < source.size (); i++) destination.push_back (source [i].get_string ());
  }
//};

db_run::db_run (int run_number): db_acl(), TObject(), b_precalib_present(0), b_precalib_quality_present(0), 
				 b_write_laben_precalib(0), b_write_laben_precalib_quality(0), b_write_muon_precalib(0), 
				 b_laben_laser_present(0), b_write_laben_laser_time(0), b_write_laben_laser_charge(0), 
				 b_laben_calib_tt1_present(0), b_write_laben_tt1_charge(0),
				 b_muon_laser_present(0), b_write_muon_laser_time(0), b_write_muon_laser_charge(0),
				 b_trigger_parameters_present(0), b_write_trigger_parameters(0), b_laben_dark_rates_present(0), 
				 b_write_laben_dark_rates(0),  b_muon_dark_rates_present(0), b_write_muon_dark_rates(0), 
				 b_laben_electronic_channel_present(0), b_write_laben_electronic_channel(0),
				 b_muon_electronic_channel_present(0), b_write_muon_electronic_channel(0), 
				 b_muon_alignment_present(0), b_write_muon_alignment(0), b_disabled_channels_present(0), b_write_disabled_channels(0), 
                                 b_laben_dark_rates_parametrisation_present(0), b_write_laben_dark_rates_parametrisation(0) {


  set_acl ("set_laben_precalib_low_bin",             "bx_precalib_laben_adc");
  set_acl ("set_laben_precalib_high_bin",            "bx_precalib_laben_adc");
  set_acl ("set_laben_precalib_gray_shift",          "bx_precalib_laben_gray");
  set_acl ("set_laben_precalib_rising_on_even",      "bx_precalib_laben_phases");
  set_acl ("set_laben_precalib_delta80",             "bx_precalib_laben_d80");
  set_acl ("set_muon_precalib_pulse_time",           "bx_precalib_muon_findpulse");
  set_acl ("set_muon_precalib_pedestal", 	     "bx_precalib_muon_pedestals");
  set_acl ("set_muon_precalib_pedsigma", 	     "bx_precalib_muon_pedestals");
  set_acl ("set_laben_precalib_mean_time",           "bx_calib_laben_decoding");
  set_acl ("set_laben_precalib_sigma_time",          "bx_calib_laben_decoding");
  set_acl ("set_laben_precalib_bad_channels",        "bx_calib_laben_decoding");
  set_acl ("set_laben_precalib_off_channels",        "bx_calib_laben_decoding");
  set_acl ("set_laben_time_offset",                  "bx_calib_laben_time_align");
  set_acl ("set_laben_time_sigma",                   "bx_calib_laben_time_align");
  set_acl ("set_laben_charge_peak",                  "bx_calib_laben_charge_peak");
  set_acl ("set_laben_charge_sigma",                 "bx_calib_laben_charge_peak");
  set_acl ("set_laben_charge_tt1_peak",              "bx_calib_laben_charge_tt1");
  set_acl ("set_laben_charge_tt1_sigma",             "bx_calib_laben_charge_tt1");
  set_acl ("set_laben_charge_tt1_mean",              "bx_calib_laben_charge_tt1");
  set_acl ("set_laben_charge_tt1_rms",               "bx_calib_laben_charge_tt1");
  set_acl ("set_laben_charge_tt1_p0",                "bx_calib_laben_charge_tt1");
  set_acl ("set_muon_time_offset",                   "bx_calib_muon_time_align");
  set_acl ("set_muon_time_sigma",                    "bx_calib_muon_time_align");
  set_acl ("set_muon_charge_peak",                   "bx_calib_muon_charge_peak");
  set_acl ("set_muon_charge_sigma",                  "bx_calib_muon_charge_peak");
  set_acl ("set_laben_gate_width",                   "bx_calib_gate");
  set_acl ("set_laben_gate_start",                   "bx_calib_gate");
  set_acl ("set_laben_laser_offset",                 "bx_calib_gate");
  set_acl ("set_laben_pulser_offset",                "bx_calib_gate");
  set_acl ("set_laben_cluster_offset",               "bx_calib_gate");
  set_acl ("set_laben_mean_dark_noise" ,             "bx_calib_laben_dark_rates");
  set_acl ("set_laben_mean_dark_sigma",              "bx_calib_laben_dark_rates");
  set_acl ("set_laben_dead_cone",                    "bx_calib_laben_dark_rates");
  set_acl ("set_laben_dead_no_cone",                 "bx_calib_laben_dark_rates");
  set_acl ("set_laben_hot_cone",                     "bx_calib_laben_dark_rates");
  set_acl ("set_laben_hot_no_cone",                  "bx_calib_laben_dark_rates");
  set_acl ("set_laben_dark_noise", 	             "bx_calib_laben_dark_rates");
  set_acl ("set_laben_dark_sigma",                   "bx_calib_laben_dark_rates");
  set_acl ("set_laben_pmt_status",                   "bx_calib_laben_dark_rates");
  set_acl ("set_muon_mean_dark_noise" ,              "bx_calib_muon_dark_rates");
  set_acl ("set_muon_mean_dark_sigma",               "bx_calib_muon_dark_rates");
  set_acl ("set_muon_dead",                          "bx_calib_muon_dark_rates");
  set_acl ("set_muon_hot",                           "bx_calib_muon_dark_rates");
  set_acl ("set_muon_dark_noise", 	             "bx_calib_muon_dark_rates");
  set_acl ("set_muon_dark_sigma",                    "bx_calib_muon_dark_rates");
  set_acl ("set_muon_pmt_status",                    "bx_calib_muon_dark_rates");
  set_acl ("set_laben_charge_base_status",           "bx_calib_laben_electronics");
  set_acl ("set_laben_charge_peak_status",           "bx_calib_laben_electronics");
  set_acl ("set_laben_charge_status",                "bx_calib_laben_electronics");
  set_acl ("set_laben_timing_status",                "bx_calib_laben_electronics");
  set_acl ("set_laben_multiplicity",                 "bx_calib_laben_electronics");
  set_acl ("set_laben_status",                       "bx_calib_laben_electronics");
  set_acl ("set_muon_multiplicity",                  "bx_calib_muon_electronics");
  set_acl ("set_muon_alignment",                     "bx_calib_muon_alignment");
  set_acl ("write_laben_precalib", 	             "bx_precalib_laben_check_tdc");
  set_acl ("set_mean_dark_rate_per_used_pmt",        "bx_laben_dark_rates");
  set_acl ("set_error_mean_dark_rate_per_used_pmt",  "bx_laben_dark_rates");
  set_acl ("set_mu_win1",                            "bx_laben_dark_rates");
  set_acl ("set_mu_win2",                            "bx_laben_dark_rates");
  set_acl ("set_pois_const_win1",                    "bx_laben_dark_rates");
  set_acl ("set_error_pois_const_win1",              "bx_laben_dark_rates");
  set_acl ("set_mu_fit_win1",                        "bx_laben_dark_rates");
  set_acl ("set_error_mu_fit_win1",                  "bx_laben_dark_rates");
  set_acl ("set_exp_const_win1",                     "bx_laben_dark_rates");
  set_acl ("set_error_exp_const_win1",               "bx_laben_dark_rates");
  set_acl ("set_tau_win1",                           "bx_laben_dark_rates");
  set_acl ("set_error_tau_win1",                     "bx_laben_dark_rates");
  set_acl ("set_pois_const_win2",                    "bx_laben_dark_rates");
  set_acl ("set_error_pois_const_win2",              "bx_laben_dark_rates");
  set_acl ("set_mu_fit_win2",                        "bx_laben_dark_rates");
  set_acl ("set_error_mu_fit_win2",                  "bx_laben_dark_rates");
  set_acl ("set_exp_const_win2",                     "bx_laben_dark_rates");
  set_acl ("set_error_exp_const_win2",               "bx_laben_dark_rates");
  set_acl ("set_tau_win2",                           "bx_laben_dark_rates");
  set_acl ("set_error_tau_win2",                     "bx_laben_dark_rates");
  set_acl ("write_laben_precalib", 	             "bx_calib_laben_decoding");
  set_acl ("write_muon_precalib",                    "bx_precalib_muon_pedestals");
  set_acl ("write_laben_precalib_quality", 	     "bx_calib_laben_decoding");
  set_acl ("write_laben_laser_time", 	             "bx_calib_laben_time_align");
  set_acl ("write_laben_laser_charge",               "bx_calib_laben_charge_peak");
  set_acl ("write_laben_tt1_charge",                 "bx_calib_laben_charge_tt1");
  set_acl ("write_muon_laser_time", 	             "bx_calib_muon_time_align");
  set_acl ("write_muon_laser_charge",                "bx_calib_muon_charge_peak");
  set_acl ("write_laben_dark_rates",                 "bx_calib_laben_dark_rates");
  set_acl ("write_muon_dark_rates",                  "bx_calib_muon_dark_rates");
  set_acl ("write_trigger_parameters",               "bx_calib_gate");
  set_acl ("write_laben_electronic_channel",         "bx_calib_laben_electronics");
  set_acl ("write_muon_electronic_channel",          "bx_calib_muon_electronics");
  set_acl ("write_muon_alignment",                   "bx_calib_muon_alignment");
  set_acl ("write_disabled_channels",                "bx_detector_monitor");
  set_acl ("write_laben_dark_rates_parametrisation", "bx_laben_dark_rates");
  set_acl ("reset_precalibrations",		     "bx_reco_framework");
  set_acl ("add_disabled_channel",		     "bx_detector_monitor");

  run_type_map["normal"     ] = db_run::normal;
  run_type_map["calibration"] = db_run::calibration;
  run_type_map["laser394"   ] = db_run::laser394;
  run_type_map["laser355"   ] = db_run::laser355;
  run_type_map["laser266"   ] = db_run::laser266;
  run_type_map["random"     ] = db_run::random;
  run_type_map["pulser"     ] = db_run::pulser;
  run_type_map["custom"     ] = db_run::custom;
  run_type_map["no_file"    ] = db_run::nofile;
  run_type_map["test"       ] = db_run::test;
  run_type_map["source"     ] = db_run::random;
  run_type_map["water"      ] = db_run::water;

  bx_dbi *dbi = bx_dbi::get ();
  bx_message &msg = bx_dbi::get ()->get_message (bx_message::critic);
  msg << "db_run: ";
  
  m_read_runinfo                  (run_number, dbi);
  m_read_laben_precalib           (run_number, dbi);
  m_read_laben_precalib_quality   (run_number, dbi);
  m_read_muon_precalib            (run_number, dbi);
  m_read_laben_laser_calib        (run_number, dbi);
  m_read_laben_tt1_calib          (run_number, dbi);
  m_read_muon_laser_calib         (run_number, dbi);
  m_read_laben_dark_rates         (run_number, dbi); 
  m_read_muon_dark_rates          (run_number, dbi); 
  m_read_trigger_parameters       (run_number, dbi);
  m_read_laben_electronic_channel (run_number, dbi); 
  m_read_muon_electronic_channel  (run_number, dbi); 
  m_read_muon_alignment           (run_number, dbi); 
  m_read_laben_dark_rates_parametrisation(run_number, dbi);
  m_read_disconnected_pmts        (run_number, dbi);
  m_read_disabled_channels	  (run_number, dbi);
  m_read_laben_qe                 (run_number, dbi);
  bx_dbi::get ()->close_db_connections ();
  bx_dbi::get ()->get_message (bx_message::info) << "db_run: db_run (" << run_number << ") initialized" << dispatch;
}

void db_run::reset_precalibrations (const bx_named* obj) {
  if (!check_acl ("reset_precalibrations", obj)) return;

  laben_precalib_low_bin_v.clear ();
  laben_precalib_high_bin_v.clear ();
  laben_precalib_rising_on_even_v.clear ();
  laben_precalib_gray_shift_v.clear ();
  laben_precalib_delta80_v.clear ();
  muon_precalib_pedestal_v.clear ();
  muon_precalib_pedsigma_v.clear ();
  u4_muon_precalib_pulse_time = 0;
}


/******************************/
/********* Read methods *******/

void db_run::m_read_runinfo (int run_number, bx_dbi *dbi) {
  std::ostringstream where_str;
  where_str << "\"RunNumber\"=" << run_number;
  const bx_dbi::table &table = dbi->query (bx_dbi::daq_config, "\"Run\"", where_str.str (), "*", "", -1);

  i_run_number = table["RunNumber"][0].get_int ();
  if (i_run_number != run_number)  {
    bx_message &msg = bx_dbi::get ()->get_message (bx_message::critic);
    msg << "db_run: ";
    msg << "database corruption asked run " << run_number << " received run " << i_run_number << dispatch;
  }
  i_profile_id	= table["ProfileID"][0].get_int ();

  if (!map_check(table["Type"][0].get_string(), run_type_map)) {
    bx_dbi::get ()->get_message (bx_message::error) << "db_run: unknown run type " << table["Type"][0] << dispatch;
  } else i_run_type = run_type_map[table["Type"][0].get_string()];

  i_max_time	= table["MaxTime"][0].get_int ();
  i_max_events	= table["MaxEvent"][0].get_int ();
  i_stop_reason	= table["StopReason"][0].get_int ();
  i_number_of_files = table["NumberOfFiles"][0].get_int ();
  i_duration	= table["Duration"][0].get_int ();
  i_events	= table["Events"][0].get_int ();
  std::string t = std::string("SELECT EXTRACT(EPOCH FROM TIMESTAMP \'") + table["StartTime"][0].get_string () + "\')";
  u4_start_time	= dbi->query (bx_dbi::daq_config, t, -1)["date_part"][0].get_int ();
  t = std::string("SELECT EXTRACT(EPOCH FROM TIMESTAMP \'") + table["StopTime"][0].get_string () + "\')";
  u4_stop_time = dbi->query (bx_dbi::daq_config, t, -1)["date_part"][0].get_int ();
  
  where_str.str ("");
  where_str << "\"RunNumber\" = (SELECT MAX(\"RunNumber\") FROM \"RunToCalibProfile\" where \"RunNumber\"<=" << run_number << ")";
  i_calib_profile = dbi->query (bx_dbi::bx_calib, "\"RunToCalibProfile\"", where_str.str (), "\"CalibProfile\"", "", -1)["CalibProfile"][0].get_int ();
}


void db_run::m_read_laben_qe (int run_number, bx_dbi *dbi)
{
  std::ostringstream where_str;
  where_str.str ("");
  where_str << "\"RunNumber\"<=" << run_number; 
  int currentMarkerRun = dbi->query (bx_dbi::bx_calib, "\"QuantumEfficiency\"", where_str.str (), "MAX(\"RunNumber\")", "", -1)["max"][0].get_int ();
  where_str.str ("");
  where_str << "\"RunNumber\"=" <<  currentMarkerRun;//to be firsr run dst period
 
  const bx_dbi::table &table5_qe_lg = dbi->query (bx_dbi::bx_calib, "\"QuantumEfficiency\"", where_str.str ());
  if (table5_qe_lg.size() == 0) {
      bx_message &msg = bx_dbi::get ()->get_message (bx_message::critic);
      msg << "db_run: internal unconsistency in db. RUN " << currentMarkerRun << ": fields filled in  InnerDetectorDarkRate but not in InnerPmtsDarkRate" << dispatch ;
  };
  const bx_dbi::column& effqe = table5_qe_lg["Qe"];
  const bx_dbi::column& effqe_nc = table5_qe_lg["QeWoc"];
  const bx_dbi::column& effqe_noise = table5_qe_lg["Noise"];
  const bx_dbi::column& ch_id5     = table5_qe_lg["ChannelID"]; 
 
  for (size_t i = 0; i < ch_id5.size (); i++) {
    int lg = ch_id5[i].get_int ();
    laben_qe_v[lg] = effqe[i].get_float ();
    laben_qe_nocorr_v[lg] = effqe_nc[i].get_float ();
    laben_qe_dark_noise_v[lg] = effqe_noise[i].get_float ();
  }

}

void db_run::m_read_laben_precalib  (int run_number, bx_dbi *dbi) {
  std::ostringstream where_str;
  where_str << "\"RunNumber\"=" << run_number;
  const bx_dbi::table &table3 = dbi->query (bx_dbi::bx_precalib, "\"LabenPrecalibData\"", where_str.str ());
  if (table3.size () > 0) {
    b_precalib_present = true;
    const bx_dbi::column& ch_id = table3["ChannelID"];
    const bx_dbi::column& low = table3["LowADC"];
    const bx_dbi::column& high = table3["HighADC"];
    const bx_dbi::column& gshift = table3["GrayShift"];
    const bx_dbi::column& roe = table3["RisingOnEven"];
    const bx_dbi::column& d80 = table3["Delta80"];
    for (size_t i = 0; i < ch_id.size (); i++) {
      int lg = ch_id[i].get_int ();
      laben_precalib_low_bin_v[lg] = low[i].get_int ();
      laben_precalib_high_bin_v[lg] = high[i].get_int ();
      laben_precalib_rising_on_even_v[lg] = roe[i].get_bool ();
      laben_precalib_gray_shift_v[lg] = gshift[i].get_int ();
      laben_precalib_delta80_v[lg] = d80[i].get_float ();
    }
  }
}

void db_run::m_read_laben_precalib_quality  (int run_number, bx_dbi *dbi) {
  std::ostringstream where_str;
  where_str << "\"RunNumber\"=" << run_number;
  //  where_str << "\"RunNumber\"<=" << run_number;
  //int last_precalib_quality_run = dbi->query (bx_dbi::bx_precalib, "\"LabenPrecalibDecodingQuality\"", where_str.str (), "MAX(\"RunNumber\")", "", -1)["max"][0].get_int ();
  //where_str.str ("");
  //where_str << "\"RunNumber\"=" << last_precalib_quality_run;
  const bx_dbi::table &table3quality = dbi->query (bx_dbi::bx_precalib, "\"LabenPrecalibDecodingQuality\"", where_str.str ());
  if (table3quality.size () > 0) {
    b_precalib_quality_present = true;
    f_laben_precalib_mean_time  = table3quality ["MeanTime"] [0].get_float();
    f_laben_precalib_sigma_time = table3quality ["SigmaTime"][0].get_float();
    const vdt& bad_channel_cell = table3quality ["BadChannelsList"][0];
    const vdt& off_channel_cell = table3quality ["OffChannelsList"][0];
    for (unsigned iq = 0; iq < off_channel_cell.get_vector().size(); iq++ ) laben_precalib_off_channels_v.push_back(off_channel_cell.get_vector()[iq].get_int());
    for (unsigned iq = 0; iq < bad_channel_cell.get_vector().size(); iq++ ) laben_precalib_bad_channels_v.push_back(bad_channel_cell.get_vector()[iq].get_int());    
  } else if (b_precalib_present) {
      bx_message &msg = bx_dbi::get ()->get_message (bx_message::error);
      msg << "Precalib present while LabenPrecalibDecodingQuality is absent" << dispatch;
  }

}

void db_run::m_read_muon_precalib (int run_number, bx_dbi *dbi) {
  std::ostringstream where_str;
  where_str << "\"RunNumber\"=" << run_number;
  const bx_dbi::table &table_4 = dbi->query (bx_dbi::bx_precalib, "\"MuonPrecalibData\"", where_str.str ());
  if (table_4.size () > 0) {
    b_precalib_present = true;
    const bx_dbi::column& ch_id = table_4["ChannelID"];
    const bx_dbi::column& ped = table_4["Pedestal"];
    const bx_dbi::column& sig = table_4["PedSigma"];
    for (size_t i = 0; i < ch_id.size (); i++) {
      int lg = ch_id[i].get_int ();
      muon_precalib_pedestal_v[lg] = ped[i].get_float ();
      muon_precalib_pedsigma_v[lg] = sig[i].get_float ();
    }
  }
}
    
void db_run::m_read_laben_laser_calib (int run_number, bx_dbi *dbi) {
  std::ostringstream where_str;
  where_str << "\"RunNumber\"<=" << run_number;
  int last_laser_pmt_calib_run = dbi->query (bx_dbi::bx_calib, "\"LaserPmtCalibration\"", where_str.str (), "MAX(\"RunNumber\")", "", -1)["max"][0].get_int ();
  if (last_laser_pmt_calib_run == run_number) b_laben_laser_present = true;
  where_str.str ("");
  where_str << "\"RunNumber\"=" << last_laser_pmt_calib_run;
  const bx_dbi::table &table2 = dbi->query (bx_dbi::bx_calib, "\"LaserPmtCalibration\"", where_str.str ());
  const bx_dbi::column& ch_id2 = table2["ChannelID"];
  const bx_dbi::column& time_offset = table2["TimeOffset"];
  const bx_dbi::column& time_sigma = table2["TimeSigma"];
  const bx_dbi::column& charge_peak = table2["ChargePeak"];
  const bx_dbi::column& charge_sigma = table2["ChargeSigma"];
  for (size_t i = 0; i < ch_id2.size (); i++) {
    int lg = ch_id2[i].get_int ();
    laben_time_offset_v[lg] = time_offset[i].get_float ();
    laben_time_sigma_v[lg] = time_sigma[i].get_float ();
    laben_charge_peak_v[lg] = charge_peak[i].get_float ();
    laben_charge_sigma_v[lg] = charge_sigma[i].get_float ();
  }
}

void db_run::m_read_laben_tt1_calib (int run_number, bx_dbi *dbi) {
  std::ostringstream where_str;
  where_str << "\"RunNumber\"<=" << run_number;
  int last_neutrino_pmt_calib_run = dbi->query (bx_dbi::bx_calib, "\"NeutrinoPmtCalibration\"", where_str.str (), "MAX(\"RunNumber\")", "", -1)["max"][0].get_int ();
  if (last_neutrino_pmt_calib_run == run_number) b_laben_calib_tt1_present = true;
  where_str.str ("");
  where_str << "\"RunNumber\"=" << last_neutrino_pmt_calib_run;
  const bx_dbi::table &table2 = dbi->query (bx_dbi::bx_calib, "\"NeutrinoPmtCalibration\"", where_str.str ());
  const bx_dbi::column& ch_id2 = table2["ChannelID"];
  const bx_dbi::column& charge_peak = table2["ChargePeak"];
  const bx_dbi::column& charge_sigma = table2["ChargeSigma"];
  const bx_dbi::column& charge_mean = table2["ChargeMean"];
  const bx_dbi::column& charge_rms = table2["ChargeRms"];
  const bx_dbi::column& charge_p0 = table2["P0"];
  for (size_t i = 0; i < ch_id2.size (); i++) {
    int lg = ch_id2[i].get_int ();
    laben_charge_tt1_peak_v[lg] = charge_peak[i].get_float ();
    laben_charge_tt1_sigma_v[lg] = charge_sigma[i].get_float ();
    laben_charge_tt1_mean_v[lg] = charge_mean[i].get_float ();
    laben_charge_tt1_rms_v[lg] = charge_rms[i].get_float ();
    laben_charge_tt1_p0_v[lg] = charge_p0[i].get_float ();
  }
}


void db_run::m_read_muon_laser_calib (int run_number, bx_dbi *dbi) {
  std::ostringstream where_str;
  where_str << "\"RunNumber\"<=" << run_number;
  int last_laser_pmt_calib_run = dbi->query (bx_dbi::bx_calib, "\"MuonPmtCalibration\"", where_str.str (), "MAX(\"RunNumber\")", "", -1)["max"][0].get_int ();
  if (last_laser_pmt_calib_run == run_number) b_muon_laser_present = true;
  where_str.str ("");
  where_str << "\"RunNumber\"=" << last_laser_pmt_calib_run;
  const bx_dbi::table & table = dbi->query (bx_dbi::bx_calib, "\"MuonPmtCalibration\"", where_str.str ());
  const bx_dbi::column& ch_id        = table["ChannelID"];
  const bx_dbi::column& time_offset  = table["TimeOffset"];
  const bx_dbi::column& time_sigma   = table["TimeSigma"];
  const bx_dbi::column& charge_peak  = table["ChargePeak"];
  const bx_dbi::column& charge_sigma = table["ChargeSigma"];
  for (size_t i = 0; i < ch_id.size (); i++) {
    int lg = ch_id[i].get_int ();
    muon_time_offset_v [lg] = time_offset [i].get_float ();
    muon_time_sigma_v  [lg] = time_sigma  [i].get_float ();
    muon_charge_peak_v [lg] = charge_peak [i].get_float ();
    muon_charge_sigma_v[lg] = charge_sigma[i].get_float ();
  }
}

void db_run::m_read_laben_dark_rates (int run_number, bx_dbi *dbi) {
  std::ostringstream where_str;
  where_str << "\"RunNumber\"<=" << run_number;
  int last_detector_parameters_calib_run = dbi->query (bx_dbi::bx_calib, "\"InnerDetectorDarkRate\"", where_str.str (), "MAX(\"RunNumber\")", "", -1)["max"][0].get_int ();
  if (last_detector_parameters_calib_run == run_number) b_laben_dark_rates_present = true;
  where_str.str ("");
  where_str << "\"RunNumber\"=" << last_detector_parameters_calib_run;
  const bx_dbi::table &table4detector = dbi->query (bx_dbi::bx_calib, "\"InnerDetectorDarkRate\"", where_str.str (), "*","" ,-1);
  f_laben_mean_dark_noise   = table4detector["MeanDarkNoise"] [0].get_float();
  f_laben_mean_dark_sigma   = table4detector["MeanDarkSigma"][0].get_float();
  i_laben_dead_cone         = table4detector["DeadCone"][0].get_int();
  i_laben_dead_no_cone      = table4detector["DeadNoCone"][0].get_int();
  i_laben_hot_cone          = table4detector["HotCone"][0].get_int();
  i_laben_hot_no_cone       = table4detector["HotNoCone"][0].get_int();

  
  const bx_dbi::table &table5dark_lg = dbi->query (bx_dbi::bx_calib, "\"InnerPmtsDarkRate\"", where_str.str ());
  if (table5dark_lg.size() == 0) {
      bx_message &msg = bx_dbi::get ()->get_message (bx_message::critic);
      msg << "db_run: internal unconsistency in db. RUN " << last_detector_parameters_calib_run << ": fields filled in  InnerDetectorDarkRate but not in InnerPmtsDarkRate" << dispatch ;
  };
  const bx_dbi::column& dark_noise = table5dark_lg["DarkNoise"];
  const bx_dbi::column& dark_sigma = table5dark_lg["DarkSigma"];
  const bx_dbi::column& pmt_status = table5dark_lg["PmtStatus"];
  const bx_dbi::column& ch_id5     = table5dark_lg["ChannelID"];
  
  for (size_t i = 0; i < ch_id5.size (); i++) {
    int lg = ch_id5[i].get_int ();
    laben_dark_noise_v[lg] = dark_noise[i].get_float ();
    laben_dark_sigma_v[lg] = dark_sigma[i].get_float ();
    laben_pmt_status_v[lg] = pmt_status[i].get_string ();
  }
}  

void db_run::m_read_muon_dark_rates (int run_number, bx_dbi *dbi) {
  std::ostringstream where_str;
  where_str << "\"RunNumber\"<=" << run_number;
  int last_detector_parameters_calib_run = dbi->query (bx_dbi::bx_calib, "\"OuterDetectorDarkRate\"", where_str.str (), "MAX(\"RunNumber\")", "", -1)["max"][0].get_int ();
  if (last_detector_parameters_calib_run == run_number) b_muon_dark_rates_present = true;
  where_str.str ("");
  where_str << "\"RunNumber\"=" << last_detector_parameters_calib_run;
  const bx_dbi::table &table_detector = dbi->query (bx_dbi::bx_calib, "\"OuterDetectorDarkRate\"", where_str.str (), "*","" ,-1);
  f_muon_mean_dark_noise   = table_detector["MeanDarkNoise"] [0].get_float();
  f_muon_mean_dark_sigma   = table_detector["MeanDarkSigma"][0].get_float();
  i_muon_dead              = table_detector["Dead"][0].get_int();
  i_muon_hot               = table_detector["Hot"][0].get_int();

  
  const bx_dbi::table &table_dark_lg = dbi->query (bx_dbi::bx_calib, "\"OuterPmtsDarkRate\"", where_str.str ());
  if (table_dark_lg.size() == 0) {
      bx_message &msg = bx_dbi::get ()->get_message (bx_message::critic);
      msg << "db_run: internal unconsistency in db. RUN " << last_detector_parameters_calib_run << ": fields filled in  OuterDetectorDarkRate but not in OuterPmtsDarkRate" << dispatch ;
  };
  const bx_dbi::column& dark_noise = table_dark_lg["DarkNoise"];
  const bx_dbi::column& dark_sigma = table_dark_lg["DarkSigma"];
  const bx_dbi::column& pmt_status = table_dark_lg["PmtStatus"];
  const bx_dbi::column& ch_id      = table_dark_lg["ChannelID"];
  
  for (size_t i = 0; i < ch_id.size (); i++) {
    int lg = ch_id[i].get_int ();
    muon_dark_noise_v[lg] = dark_noise[i].get_float ();
    muon_dark_sigma_v[lg] = dark_sigma[i].get_float ();
    muon_pmt_status_v[lg] = pmt_status[i].get_string ();
  }
}  

void db_run::m_read_trigger_parameters (int run_number, bx_dbi *dbi) {
  std::ostringstream where_str;
  where_str << "\"RunNumber\"<=" << run_number;
  int last_trigger_parameters_calib_run = dbi->query (bx_dbi::bx_calib, "\"TriggerParameters\"", where_str.str (), "MAX(\"RunNumber\")", "", -1)["max"][0].get_int ();
  if (last_trigger_parameters_calib_run == run_number) b_trigger_parameters_present = true;
  where_str.str ("");
  where_str << "\"RunNumber\"=" << last_trigger_parameters_calib_run;
  const bx_dbi::table &table5detector = dbi->query (bx_dbi::bx_calib, "\"TriggerParameters\"", where_str.str (), "*","" ,-1);
  f_laben_gate_width      = table5detector["GateWidth"] [0].get_float();
  f_laben_gate_start      = table5detector["GateStart"] [0].get_float();    
  f_laben_laser_offset    = table5detector["LaserOffset"] [0].get_float();    
  f_laben_pulser_offset   = table5detector["PulserOffset"] [0].get_float();    
  f_laben_cluster_offset  = table5detector["ClusterOffset"] [0].get_float();    
}

void db_run::m_read_laben_electronic_channel (int run_number, bx_dbi *dbi) {
  std::ostringstream where_str;
  where_str << "\"RunNumber\"<=" << run_number;
  int last_laben_electronic_channel_run = dbi->query (bx_dbi::bx_calib, "\"LabenChannelsProperties\"", where_str.str (), "MAX(\"RunNumber\")", "", -1)["max"][0].get_int ();
  if (last_laben_electronic_channel_run  == run_number) b_laben_electronic_channel_present = true;
  where_str.str ("");
  where_str << "\"RunNumber\"=" << last_laben_electronic_channel_run;
  const bx_dbi::table &table = dbi->query (bx_dbi::bx_calib, "\"LabenChannelsProperties\"", where_str.str ());
  const bx_dbi::column& ch_id              = table["ChannelID"       ];
  const bx_dbi::column& charge_base_status = table["ChargeBaseStatus"];
  const bx_dbi::column& charge_peak_status = table["ChargePeakStatus"];
  const bx_dbi::column& charge_status      = table["ChargeStatus"    ];
  const bx_dbi::column& timing_status      = table["TimingStatus"    ];
  const bx_dbi::column& multiplicity       = table["Multiplicity"    ];
  for (size_t i = 0; i < ch_id.size (); i++) {
    int lg = ch_id[i].get_int ();
    copy_vdt_vector (laben_charge_base_status_v[lg], charge_base_status[i].get_vector ());
    copy_vdt_vector (laben_charge_peak_status_v[lg], charge_peak_status[i].get_vector ());
    copy_vdt_vector (laben_charge_status_v[lg], charge_status[i].get_vector ());
    copy_vdt_vector (laben_timing_status_v[lg], timing_status[i].get_vector ());
    copy_vdt_vector (laben_multiplicity_v[lg], multiplicity[i].get_vector ());
  }
}

void db_run::m_read_muon_electronic_channel (int run_number, bx_dbi *dbi) {
  std::ostringstream where_str;
  where_str << "\"RunNumber\"<=" << run_number;
  int last_muon_electronic_channel_run = dbi->query (bx_dbi::bx_calib, "\"MuonChannelsProperties\"", where_str.str (), "MAX(\"RunNumber\")", "", -1)["max"][0].get_int ();
  if (last_muon_electronic_channel_run == run_number) b_muon_electronic_channel_present = true;
  where_str.str ("");
  where_str << "\"RunNumber\"=" << last_muon_electronic_channel_run;
  const bx_dbi::table &table = dbi->query (bx_dbi::bx_calib, "\"MuonChannelsProperties\"", where_str.str ());
  const bx_dbi::column& ch_id       = table["ChannelID"   ];
  const bx_dbi::column& multiplicity= table["Multiplicity"];
  for (size_t i = 0; i < ch_id.size (); i++) {
    int mch = ch_id[i].get_int () - 3001;
    copy_vdt_vector (muon_multiplicity_v[mch], multiplicity[i].get_vector ());
  }
}

void db_run::m_read_muon_alignment (int run_number, bx_dbi *dbi) {
  std::ostringstream where_str;
  where_str << "\"RunNumber\"=" << run_number;
  const bx_dbi::table &table = dbi->query (bx_dbi::bx_calib, "\"MuonAlignment\"", where_str.str ());
  if (table.size() >0 ) {
    const bx_dbi::column& n  = table["LastEvent"       ];
    const bx_dbi::column& t  = table["LastTime"        ];
    const bx_dbi::column& un = table["LastUpEvent"     ];
    const bx_dbi::column& ut = table["LastUpTime"      ];
    const bx_dbi::column& an = table["LastAlignedEvent"];
    const bx_dbi::column& at = table["LastAlignedTime" ];
    b_muon_alignment_present = true;
    i_muon_nevents         =  n[0].get_int();
    std::string st = std::string("SELECT EXTRACT(EPOCH FROM TIMESTAMP \'") + t[0].get_string () + "\')";
    t_muon_time	    = dbi->query (bx_dbi::daq_config, st, -1)["date_part"][0].get_int ();
    i_muon_up_nevents      = un[0].get_int();
    std::string sut = std::string("SELECT EXTRACT(EPOCH FROM TIMESTAMP \'") + ut[0].get_string () + "\')";
    t_muon_up_time	   = dbi->query (bx_dbi::daq_config, sut, -1)["date_part"][0].get_int ();
    i_muon_aligned_nevents = an[0].get_int();
    std::string sat = std::string("SELECT EXTRACT(EPOCH FROM TIMESTAMP \'") + at[0].get_string () + "\')";
    t_muon_aligned_time	    = dbi->query (bx_dbi::daq_config, sat, -1)["date_part"][0].get_int ();
  } 
}

void db_run::m_read_disconnected_pmts (int run_number, bx_dbi *dbi) {
  std::ostringstream query; 
  query << "select \"ChannelID\" from \"DisconnectedPmts\" p where \"Date\" = ( select max(\"Date\") from \"DisconnectedPmts\" where \"ChannelID\" = p.\"ChannelID\") and \"Disconnected\" is true and \"RunNumber\" <= " << run_number;
  const bx_dbi::table &table = dbi->query (bx_dbi::bx_calib, query.str (), 0);

  const bx_dbi::column& ch_id = table["ChannelID"];
  for (size_t i = 0; i < ch_id.size (); i++) disconnected_pmts_v[ch_id[i].get_int ()] = true;
}

void db_run::m_read_disabled_channels (int run_number, bx_dbi *dbi) {
  std::ostringstream where_str;
  int cycle = bx_dbi::get ()->get_parameter ("disabled_channel_cycle").get_int ();
  if (cycle <= 0) cycle = CYCLE_NUMBER;
  where_str << "\"RunNumber\"=" << run_number << " AND \"Cycle\">=" << cycle;
  const bx_dbi::table &ctest = dbi->query (bx_dbi::bx_calib, "\"DisabledChannels\"", where_str.str (), "MIN(\"Cycle\")");
  if (!ctest.size () || ctest["min"][0].get_type () != vdt::int_vdt) {
    bx_message msg(bx_message::warn);
    msg << "no disabled channels for cycle " << cycle << dispatch;
    return;
  }
  cycle = ctest["min"][0].get_int ();
  where_str.str ("");
  where_str << "\"RunNumber\"=" << run_number << " AND \"Cycle\"=" << cycle;
  const bx_dbi::table &table = dbi->query (bx_dbi::bx_calib, "\"DisabledChannels\"", where_str.str ());
  if(table.size () > 0) {
    b_disabled_channels_present = true;
    const bx_dbi::column& ch_id = table["ChannelID"];
    const bx_dbi::column& evnum = table["EvNum"];
    const bx_dbi::column& t_on = table["Timing"];
    const bx_dbi::column& c_on = table["Charge"];

    for (size_t i = 0; i < ch_id.size (); i++) {
      int lg = ch_id[i].get_int ();
      int e = evnum[i].get_int ();
      bool t = t_on[i].get_bool ();
      bool c = c_on[i].get_bool ();
      if (t) disabled_channels_v[e][lg] = timing;
      else if (c && disabled_channels_v[e][lg] != timing) disabled_channels_v[e][lg] = charge;
    }
  }
}

void db_run::m_read_laben_dark_rates_parametrisation(int run_number, bx_dbi* dbi) {
  std::ostringstream where_str;
  where_str << "\"RunNumber\"=" << run_number << " AND \"EchidnaCycle\" = " << CYCLE_NUMBER;
  const bx_dbi::table &table = dbi->query (bx_dbi::bx_calib, "\"InnerDetectorDarkRateParametrisation\"", where_str.str (), "*", "", 1);
  if(table.size () > 0) {
    f_mean_dark_rate_per_used_pmt        = table["MeanDarkRatePerUsedPmt"     ] [0].get_float ();
    f_error_mean_dark_rate_per_used_pmt  = table["ErrorMeanDarkRatePerUsedPmt"] [0].get_float ();
    f_mu_win1                            = table["MuWin1"                     ] [0].get_float ();
    f_mu_win2                            = table["MuWin2"                     ] [0].get_float ();
    f_pois_const_win1                    = table["PoisConstWin1"              ] [0].get_float ();
    f_error_pois_const_win1              = table["ErrorPoisConstWin1"         ] [0].get_float ();
    f_mu_fit_win1                        = table["MuFitWin1"                  ] [0].get_float ();
    f_error_mu_fit_win1                  = table["ErrorMuFitWin1"             ] [0].get_float ();
    f_exp_const_win1                     = table["ExpConstWin1"               ] [0].get_float ();
    f_error_exp_const_win1               = table["ErrorExpConstWin1"          ] [0].get_float ();
    f_tau_win1                           = table["TauWin1"                    ] [0].get_float ();
    f_error_tau_win1                     = table["ErrorTauWin1"               ] [0].get_float ();
    f_pois_const_win2                    = table["PoisConstWin2"              ] [0].get_float ();
    f_error_pois_const_win2              = table["ErrorPoisConstWin2"         ] [0].get_float ();
    f_mu_fit_win2                        = table["MuFitWin2"                  ] [0].get_float ();
    f_error_mu_fit_win2                  = table["ErrorMuFitWin2"             ] [0].get_float ();
    f_exp_const_win2                     = table["ExpConstWin2"               ] [0].get_float ();
    f_error_exp_const_win2               = table["ErrorExpConstWin2"          ] [0].get_float ();
    f_tau_win2                           = table["TauWin2"                    ] [0].get_float ();
    f_error_tau_win2                     = table["ErrorTauWin2"               ] [0].get_float ();
    b_laben_dark_rates_parametrisation_present = true;
  }
}

/**************************************/
/********** write methods *************/

void db_run::m_write_laben_precalib () {
  if (b_write_laben_precalib && !b_precalib_present) {
    bx_dbi::table laben_precalib_write_table("laben_precalib_write_table");
    laben_precalib_write_table["RunNumber"].reserve (constants::laben::channels);
    laben_precalib_write_table["ChannelID"].reserve (constants::laben::channels);
    laben_precalib_write_table["LowADC"].reserve (constants::laben::channels);
    laben_precalib_write_table["HighADC"].reserve (constants::laben::channels);
    laben_precalib_write_table["GrayShift"].reserve (constants::laben::channels);
    laben_precalib_write_table["RisingOnEven"].reserve (constants::laben::channels);
    laben_precalib_write_table["Delta80"].reserve (constants::laben::channels);
    bx_dbi::column& run_number = laben_precalib_write_table["RunNumber"];
    bx_dbi::column& ch_id = laben_precalib_write_table["ChannelID"];
    bx_dbi::column& low_bin = laben_precalib_write_table["LowADC"];
    bx_dbi::column& high_bin = laben_precalib_write_table["HighADC"];
    bx_dbi::column& gray_shift = laben_precalib_write_table["GrayShift"];
    bx_dbi::column& rising_on_even = laben_precalib_write_table["RisingOnEven"];
    bx_dbi::column& d80 = laben_precalib_write_table["Delta80"];
    for (int i = 0; i < constants::laben::channels; i++) {
      int lg = i + 1;
      if (!check_laben_precalib_delta80 (lg)) continue; // Assume precalib are present for this channel if d80 is calculated (d80 is the last precalib done)
      run_number.push_back (i_run_number);
      ch_id.push_back (lg);
      low_bin.push_back (laben_precalib_low_bin_v[lg]);
      high_bin.push_back (laben_precalib_high_bin_v[lg]);
      gray_shift.push_back (laben_precalib_gray_shift_v[lg]);
      rising_on_even.push_back (laben_precalib_rising_on_even_v[lg]);
      d80.push_back (laben_precalib_delta80_v[lg]);
    }
    bx_dbi::column_names index_columns;
    index_columns.push_back ("RunNumber");
    bx_dbi::get ()->insert (bx_dbi::bx_precalib, "\"LabenPrecalibData\"", laben_precalib_write_table, index_columns);
  }
}

void db_run::m_write_laben_precalib_quality () {
  if (b_write_laben_precalib_quality && !b_precalib_quality_present) {
    bx_dbi::table laben_precalib_quality_write_table("laben_precalib_quality_write_table");
    laben_precalib_quality_write_table["RunNumber"].push_back (i_run_number);
    laben_precalib_quality_write_table["MeanTime"].push_back (f_laben_precalib_mean_time);
    laben_precalib_quality_write_table["SigmaTime"].push_back (f_laben_precalib_sigma_time);
    laben_precalib_quality_write_table["BadChannelsList"].push_back (laben_precalib_bad_channels_v);
    laben_precalib_quality_write_table["BadChannelsNumber"].push_back (int (laben_precalib_bad_channels_v.size ()));
    laben_precalib_quality_write_table["OffChannelsList"].push_back (laben_precalib_off_channels_v);
    laben_precalib_quality_write_table["OffChannelsNumber"].push_back (int (laben_precalib_off_channels_v.size ()));
    bx_dbi::column_names index_columns;
    index_columns.push_back ("RunNumber");
    bx_dbi::get ()->insert (bx_dbi::bx_precalib, "\"LabenPrecalibDecodingQuality\"", laben_precalib_quality_write_table, index_columns);
  }
}

void db_run::m_write_muon_precalib () {
  if (b_write_muon_precalib && !b_precalib_present) {
    bx_dbi::table muon_precalib_write_table("muon_precalib_write_table");
    muon_precalib_write_table["RunNumber"].reserve (constants::muon::channels);
    muon_precalib_write_table["ChannelID"].reserve (constants::muon::channels);
    muon_precalib_write_table["Pedestal"].reserve (constants::muon::channels);
    muon_precalib_write_table["PedSigma"].reserve (constants::muon::channels);
    bx_dbi::column& run_number = muon_precalib_write_table["RunNumber"];
    bx_dbi::column& ch_id = muon_precalib_write_table["ChannelID"];
    bx_dbi::column& ped = muon_precalib_write_table["Pedestal"];
    bx_dbi::column& sig = muon_precalib_write_table["PedSigma"];
    for (int i = 0; i < constants::muon::channels; i++) {
      int lg = i + 1 + constants::muon::channel_offset;
      if (!check_muon_precalib_pedestal (lg)) continue;
      run_number.push_back ((int32_t)i_run_number);
      ch_id.push_back ((int32_t)lg);
      ped.push_back ((double)muon_precalib_pedestal_v[lg]);
      sig.push_back ((double)muon_precalib_pedsigma_v[lg]);
    }
    bx_dbi::column_names index_columns;
    index_columns.push_back ("RunNumber");
    bx_dbi::get ()->insert (bx_dbi::bx_precalib, "\"MuonPrecalibData\"", muon_precalib_write_table, index_columns);
  }
}

void db_run::m_write_laben_laser_calib () {
  if (b_write_laben_laser_time && b_write_laben_laser_charge && !b_laben_laser_present) {
    bx_dbi::table laben_laser_write_table("laben_laser_write_table");
    laben_laser_write_table["RunNumber"].reserve (constants::laben::channels);
    laben_laser_write_table["ChannelID"].reserve (constants::laben::channels);
    laben_laser_write_table["TimeOffset"].reserve (constants::laben::channels);
    laben_laser_write_table["TimeSigma"].reserve (constants::laben::channels);
    laben_laser_write_table["ChargePeak"].reserve (constants::laben::channels);
    laben_laser_write_table["ChargeSigma"].reserve (constants::laben::channels);
    bx_dbi::column& run_number = laben_laser_write_table["RunNumber"];
    bx_dbi::column& ch_id = laben_laser_write_table["ChannelID"];
    bx_dbi::column& tioff = laben_laser_write_table["TimeOffset"];
    bx_dbi::column& tisig = laben_laser_write_table["TimeSigma"];
    bx_dbi::column& chpea = laben_laser_write_table["ChargePeak"];
    bx_dbi::column& chsig = laben_laser_write_table["ChargeSigma"];
    for (int i = 0; i < constants::laben::channels; i++) {
      int lg = i + 1;
      run_number.push_back ((int32_t)i_run_number);
      ch_id.push_back ((int32_t)lg);
      tioff.push_back ((double)laben_time_offset_v[lg]);
      tisig.push_back ((double)laben_time_sigma_v[lg]);
      chpea.push_back ((double)laben_charge_peak_v[lg]);
      chsig.push_back ((double)laben_charge_sigma_v[lg]);
    }
    bx_dbi::column_names index_columns;
    index_columns.push_back ("RunNumber");
    bx_dbi::get ()->insert (bx_dbi::bx_calib, "\"LaserPmtCalibration\"", laben_laser_write_table, index_columns);
  }
}

void db_run::m_write_laben_tt1_calib () {
  if (b_write_laben_tt1_charge && !b_laben_calib_tt1_present) {
    bx_dbi::table laben_tt1_write_table("laben_tt1_write_table");
    laben_tt1_write_table["RunNumber"].reserve (constants::laben::channels);
    laben_tt1_write_table["ChannelID"].reserve (constants::laben::channels);
    laben_tt1_write_table["ChargePeak"].reserve (constants::laben::channels);
    laben_tt1_write_table["ChargeSigma"].reserve (constants::laben::channels);
    laben_tt1_write_table["ChargeMean"].reserve (constants::laben::channels);
    laben_tt1_write_table["ChargeRms"].reserve (constants::laben::channels);
    laben_tt1_write_table["P0"].reserve (constants::laben::channels);
    bx_dbi::column& run_number = laben_tt1_write_table["RunNumber"];
    bx_dbi::column& ch_id = laben_tt1_write_table["ChannelID"];
    bx_dbi::column& chpea = laben_tt1_write_table["ChargePeak"];
    bx_dbi::column& chsig = laben_tt1_write_table["ChargeSigma"];
    bx_dbi::column& chmean = laben_tt1_write_table["ChargeMean"];
    bx_dbi::column& chrms = laben_tt1_write_table["ChargeRms"];
    bx_dbi::column& chp0 = laben_tt1_write_table["P0"];
    for (int i = 0; i < constants::laben::channels; i++) {
      int lg = i + 1;
      run_number.push_back ((int32_t)i_run_number);
      ch_id.push_back ((int32_t)lg);
      chpea.push_back ((double)laben_charge_tt1_peak_v[lg]);
      chsig.push_back ((double)laben_charge_tt1_sigma_v[lg]);
      chmean.push_back ((double)laben_charge_tt1_mean_v[lg]);
      chrms.push_back ((double)laben_charge_tt1_rms_v[lg]);
      chp0.push_back ((double)laben_charge_tt1_p0_v[lg]);
    }
    bx_dbi::column_names index_columns;
    index_columns.push_back ("RunNumber");
    bx_dbi::get ()->insert (bx_dbi::bx_calib, "\"NeutrinoPmtCalibration\"", laben_tt1_write_table, index_columns);
  }
}


void db_run::m_write_muon_laser_calib () {
  if (b_write_muon_laser_time && b_write_muon_laser_charge && !b_muon_laser_present) {
    bx_dbi::table muon_laser_write_table("muon_laser_write_table");
    muon_laser_write_table["RunNumber"  ].reserve (constants::muon::channels);
    muon_laser_write_table["ChannelID"  ].reserve (constants::muon::channels);
    muon_laser_write_table["TimeOffset" ].reserve (constants::muon::channels);
    muon_laser_write_table["TimeSigma"  ].reserve (constants::muon::channels);
    muon_laser_write_table["ChargePeak" ].reserve (constants::muon::channels);
    muon_laser_write_table["ChargeSigma"].reserve (constants::muon::channels);
    bx_dbi::column& run_number = muon_laser_write_table["RunNumber"];
    bx_dbi::column& ch_id = muon_laser_write_table["ChannelID"];
    bx_dbi::column& tioff = muon_laser_write_table["TimeOffset"];
    bx_dbi::column& tisig = muon_laser_write_table["TimeSigma"];
    bx_dbi::column& chpea = muon_laser_write_table["ChargePeak"];
    bx_dbi::column& chsig = muon_laser_write_table["ChargeSigma"];
    for (int i = 0; i < constants::muon::channels; i++) {
      int lg = i + 1 + constants::muon::channel_offset;
      run_number.push_back ((int32_t)i_run_number);
      ch_id.push_back ((int32_t)lg);
      tioff.push_back ((double)muon_time_offset_v[lg]);
      tisig.push_back ((double)muon_time_sigma_v[lg]);
      chpea.push_back ((double)muon_charge_peak_v[lg]);
      chsig.push_back ((double)muon_charge_sigma_v[lg]);
    }
    bx_dbi::column_names index_columns;
    index_columns.push_back ("RunNumber");
    bx_dbi::get ()->insert (bx_dbi::bx_calib, "\"MuonPmtCalibration\"", muon_laser_write_table, index_columns);
  }
}

void db_run::m_write_trigger_parameters () {
  if (b_write_trigger_parameters && !b_trigger_parameters_present) {
    bx_dbi::table trigger_parameters_write_table("trigger_parameters_write_table");
    trigger_parameters_write_table["RunNumber"].push_back (i_run_number);
    trigger_parameters_write_table["GateWidth"].push_back (f_laben_gate_width);
    trigger_parameters_write_table["GateStart"].push_back (f_laben_gate_start);
    trigger_parameters_write_table["LaserOffset"].push_back (f_laben_laser_offset);
    trigger_parameters_write_table["PulserOffset"].push_back (f_laben_pulser_offset);
    trigger_parameters_write_table["ClusterOffset"].push_back (f_laben_cluster_offset);
    bx_dbi::column_names index_columns;
    index_columns.push_back ("RunNumber");
    bx_dbi::get ()->insert (bx_dbi::bx_calib, "\"TriggerParameters\"", trigger_parameters_write_table, index_columns);
  }
}

void db_run::m_write_muon_alignment () {
  if (b_write_muon_alignment && !b_muon_alignment_present) {
    bx_dbi::table muon_alignments_write_table("muon_alignments_write_table");
    muon_alignments_write_table["RunNumber"       ].push_back (i_run_number          );
    muon_alignments_write_table["LastEvent"       ].push_back (i_muon_nevents        );
    muon_alignments_write_table["LastTime"        ].push_back (std::string(ctime(&t_muon_time))    );
    muon_alignments_write_table["LastUpEvent"     ].push_back (i_muon_up_nevents     );
    muon_alignments_write_table["LastUpTime"      ].push_back (std::string(ctime(&t_muon_up_time)) );
    muon_alignments_write_table["LastAlignedEvent"].push_back (i_muon_aligned_nevents);
    muon_alignments_write_table["LastAlignedTime" ].push_back (std::string(ctime(&t_muon_aligned_time)));
    bx_dbi::column_names index_columns;
    index_columns.push_back ("RunNumber");
    bx_dbi::get ()->insert (bx_dbi::bx_calib, "\"MuonAlignment\"", muon_alignments_write_table, index_columns);
  }
}

void db_run::m_write_laben_dark_rates () {
  if (b_write_laben_dark_rates && !b_laben_dark_rates_present) {
    bx_dbi::table laben_detector_dark_rates_write_table("laben_detector_dark_rates_write_table");
    laben_detector_dark_rates_write_table["RunNumber"].push_back (i_run_number);
    laben_detector_dark_rates_write_table["MeanDarkNoise"].push_back (f_laben_mean_dark_noise);
    laben_detector_dark_rates_write_table["MeanDarkSigma"].push_back (f_laben_mean_dark_sigma);
    laben_detector_dark_rates_write_table["DeadCone"].push_back (i_laben_dead_cone);	 
    laben_detector_dark_rates_write_table["DeadNoCone"].push_back (i_laben_dead_no_cone);   
    laben_detector_dark_rates_write_table["HotCone"].push_back (i_laben_hot_cone);	 
    laben_detector_dark_rates_write_table["HotNoCone"].push_back (i_laben_hot_no_cone);	 
    bx_dbi::column_names index_columns;
    index_columns.push_back ("RunNumber");
    bx_dbi::get ()->insert (bx_dbi::bx_calib, "\"InnerDetectorDarkRate\"", laben_detector_dark_rates_write_table, index_columns);

    bx_dbi::table laben_pmt_dark_rates_write_table("laben_pmt_dark_rates_write_table");
    laben_pmt_dark_rates_write_table["RunNumber"].reserve (constants::laben::channels);
    laben_pmt_dark_rates_write_table["ChannelID"].reserve (constants::laben::channels);
    laben_pmt_dark_rates_write_table["DarkNoise"].reserve (constants::laben::channels);
    laben_pmt_dark_rates_write_table["DarkSigma"].reserve (constants::laben::channels);
    laben_pmt_dark_rates_write_table["PmtStatus"].reserve (constants::laben::channels);
    bx_dbi::column& run_number  = laben_pmt_dark_rates_write_table["RunNumber"];
    bx_dbi::column& ch_id       = laben_pmt_dark_rates_write_table["ChannelID"];
    bx_dbi::column& dark_noise  = laben_pmt_dark_rates_write_table["DarkNoise"];
    bx_dbi::column& dark_sigma  = laben_pmt_dark_rates_write_table["DarkSigma"];
    bx_dbi::column& pmt_status  = laben_pmt_dark_rates_write_table["PmtStatus"];
    for (int i = 0; i < constants::laben::channels; i++) {
      int lg = i + 1;
      run_number.push_back (i_run_number);
      ch_id.push_back (lg);
      dark_noise.push_back (laben_dark_noise_v[lg]);
      dark_sigma.push_back (laben_dark_sigma_v[lg]);
      pmt_status.push_back (laben_pmt_status_v[lg]);	  
    }
    bx_dbi::get ()->insert (bx_dbi::bx_calib, "\"InnerPmtsDarkRate\"", laben_pmt_dark_rates_write_table, index_columns);
  }
}

void db_run::m_write_muon_dark_rates () {
  if (b_write_muon_dark_rates && !b_muon_dark_rates_present) {
    bx_dbi::table muon_detector_dark_rates_write_table("muon_detector_dark_rates_write_table");
    muon_detector_dark_rates_write_table["RunNumber"].push_back (i_run_number);
    muon_detector_dark_rates_write_table["MeanDarkNoise"].push_back (f_muon_mean_dark_noise);
    muon_detector_dark_rates_write_table["MeanDarkSigma"].push_back (f_muon_mean_dark_sigma);
    muon_detector_dark_rates_write_table["Dead"].push_back (i_muon_dead);	 
    muon_detector_dark_rates_write_table["Hot"].push_back (i_muon_hot);	 
    bx_dbi::column_names index_columns;
    index_columns.push_back ("RunNumber");
    bx_dbi::get ()->insert (bx_dbi::bx_calib, "\"OuterDetectorDarkRate\"", muon_detector_dark_rates_write_table, index_columns);

    bx_dbi::table muon_pmt_dark_rates_write_table("muon_pmt_dark_rates_write_table");
    muon_pmt_dark_rates_write_table["RunNumber"].reserve (constants::muon::channels);
    muon_pmt_dark_rates_write_table["ChannelID"].reserve (constants::muon::channels);
    muon_pmt_dark_rates_write_table["DarkNoise"].reserve (constants::muon::channels);
    muon_pmt_dark_rates_write_table["DarkSigma"].reserve (constants::muon::channels);
    muon_pmt_dark_rates_write_table["PmtStatus"].reserve (constants::muon::channels);
    bx_dbi::column& run_number  = muon_pmt_dark_rates_write_table["RunNumber"];
    bx_dbi::column& ch_id       = muon_pmt_dark_rates_write_table["ChannelID"];
    bx_dbi::column& dark_noise  = muon_pmt_dark_rates_write_table["DarkNoise"];
    bx_dbi::column& dark_sigma  = muon_pmt_dark_rates_write_table["DarkSigma"];
    bx_dbi::column& pmt_status  = muon_pmt_dark_rates_write_table["PmtStatus"];
    for (int i = 0; i < constants::muon::channels; i++) {
      int lg = i + constants::muon::channel_offset + 1;
      run_number.push_back (i_run_number);
      ch_id.push_back (lg);
      dark_noise.push_back (muon_dark_noise_v[lg]);
      dark_sigma.push_back (muon_dark_sigma_v[lg]);
      pmt_status.push_back (muon_pmt_status_v[lg]);	  
    }
    bx_dbi::get ()->insert (bx_dbi::bx_calib, "\"OuterPmtsDarkRate\"", muon_pmt_dark_rates_write_table, index_columns);
  }
}

void db_run::m_write_laben_electronic_channel() {
  if (b_write_laben_electronic_channel && !b_laben_electronic_channel_present) {
    bx_dbi::table electronic_channel_write_table("electronic_channel_write_table");
    electronic_channel_write_table["RunNumber"].reserve (constants::laben::channels);
    electronic_channel_write_table["ChannelID"].reserve (constants::laben::channels); 
    electronic_channel_write_table["ChargeBaseStatus"].reserve (constants::laben::channels);
    electronic_channel_write_table["ChargePeakStatus"].reserve (constants::laben::channels);
    electronic_channel_write_table["ChargeStatus"].reserve (constants::laben::channels);
    electronic_channel_write_table["TimingStatus"].reserve (constants::laben::channels);
    electronic_channel_write_table["Multiplicity"].reserve (constants::laben::channels);
    bx_dbi::column_names index_columns;
    index_columns.push_back ("RunNumber");
    bx_dbi::column& run_number  = electronic_channel_write_table["RunNumber"];
    bx_dbi::column& ch_id       = electronic_channel_write_table["ChannelID"];
    bx_dbi::column& charge_base_status = electronic_channel_write_table["ChargeBaseStatus"];
    bx_dbi::column& charge_peak_status = electronic_channel_write_table["ChargePeakStatus"];
    bx_dbi::column& charge_status = electronic_channel_write_table["ChargeStatus"];
    bx_dbi::column& timing_status = electronic_channel_write_table["TimingStatus"]; 
    bx_dbi::column& multiplicity = electronic_channel_write_table["Multiplicity"];
    for (int i = 0; i < constants::laben::channels; i++) {
      int lg = i + 1;
      run_number.push_back (i_run_number);
      ch_id.push_back (lg);
      charge_base_status.push_back (laben_charge_base_status_v[lg]);
      charge_peak_status.push_back (laben_charge_peak_status_v[lg]);
      charge_status.push_back (laben_charge_status_v[lg]);
      timing_status.push_back  (laben_timing_status_v[lg]);
      multiplicity.push_back (laben_multiplicity_v[lg]);
    }
    bx_dbi::get ()->insert (bx_dbi::bx_calib, "\"LabenChannelsProperties\"", electronic_channel_write_table, index_columns);
  } 
}

void db_run::m_write_muon_electronic_channel() {
  if (b_write_muon_electronic_channel && !b_muon_electronic_channel_present) {
    bx_dbi::table muon_electronic_channel_write_table("muon_electronic_channel_write_table");
    muon_electronic_channel_write_table["RunNumber"].reserve (constants::muon::channels);
    muon_electronic_channel_write_table["ChannelID"].reserve (constants::muon::channels); 
    muon_electronic_channel_write_table["Multiplicity"].reserve (constants::muon::channels);
    bx_dbi::column_names index_columns;
    index_columns.push_back ("RunNumber");
    bx_dbi::column& run_number  = muon_electronic_channel_write_table["RunNumber"];
    bx_dbi::column& ch_id       = muon_electronic_channel_write_table["ChannelID"];
    bx_dbi::column& multiplicity = muon_electronic_channel_write_table["Multiplicity"];
    for (int i = 0; i < constants::muon::channels; i++) {
      int lg = i + 3001;
      run_number.push_back (i_run_number);
      ch_id.push_back (lg);
      multiplicity.push_back (muon_multiplicity_v[i]);
    }
    bx_dbi::get ()->insert (bx_dbi::bx_calib, "\"MuonChannelsProperties\"", muon_electronic_channel_write_table, index_columns);
  } 
}

void db_run::m_write_disabled_channels () {
  if(b_write_disabled_channels && !b_disabled_channels_present) {
    bx_dbi::table disabled_channels("disabled_channels");
    
    disabled_channels["RunNumber"].reserve (constants::laben::channels);
    disabled_channels["ChannelID"].reserve (constants::laben::channels);
    disabled_channels["EvNum"].reserve (constants::laben::channels);
    disabled_channels["Timing"].reserve (constants::laben::channels);
    disabled_channels["Charge"].reserve (constants::laben::channels);
    disabled_channels["Cycle"].reserve (constants::laben::channels);
    bx_dbi::column& run_number  = disabled_channels["RunNumber"];
    bx_dbi::column& ch_id       = disabled_channels["ChannelID"];
    bx_dbi::column& evnum       = disabled_channels["EvNum"];
    bx_dbi::column& t           = disabled_channels["Timing"];
    bx_dbi::column& c           = disabled_channels["Charge"];
    bx_dbi::column& cy          = disabled_channels["Cycle"];
    
    bx_dbi::column_names index_columns;
    index_columns.push_back ("RunNumber");
    index_columns.push_back ("Cycle");
    
    for (std::map<int, std::map<int, db_run::disabled_type> >::const_iterator it = disabled_channels_v.begin (); it != disabled_channels_v.end (); it++) {
      int ev = it->first;
      for (std::map<int, db_run::disabled_type>::const_iterator jt = it->second.begin (); jt != it->second.end (); jt++) {
	int lg = jt->first;
	run_number.push_back (i_run_number);
	ch_id.push_back (lg);
	evnum.push_back (ev);
	t.push_back (jt->second == timing);
	c.push_back (jt->second == charge);
	cy.push_back (CYCLE_NUMBER);
      }
    }

    bx_dbi::get ()->insert (bx_dbi::bx_calib, "\"DisabledChannels\"", disabled_channels, index_columns);
  }
}

void db_run::m_write_laben_dark_rates_parametrisation() {
  if(b_write_laben_dark_rates_parametrisation && !b_laben_dark_rates_parametrisation_present) {
    bx_dbi::table inner_detector_dark_rates_parametrisation_write_table("inner_detector_dark_rates_parametrisation_write_table");
    inner_detector_dark_rates_parametrisation_write_table["RunNumber"                  ].push_back (i_run_number);
    inner_detector_dark_rates_parametrisation_write_table["EchidnaCycle"               ].push_back (CYCLE_NUMBER);
    inner_detector_dark_rates_parametrisation_write_table["MeanDarkRatePerUsedPmt"     ].push_back (f_mean_dark_rate_per_used_pmt      );
    inner_detector_dark_rates_parametrisation_write_table["ErrorMeanDarkRatePerUsedPmt"].push_back (f_error_mean_dark_rate_per_used_pmt);
    inner_detector_dark_rates_parametrisation_write_table["MuWin1"                     ].push_back (f_mu_win1                          );
    inner_detector_dark_rates_parametrisation_write_table["MuWin2"                     ].push_back (f_mu_win2                          );
    inner_detector_dark_rates_parametrisation_write_table["PoisConstWin1"              ].push_back (f_pois_const_win1                  );
    inner_detector_dark_rates_parametrisation_write_table["ErrorPoisConstWin1"         ].push_back (f_error_pois_const_win1            );
    inner_detector_dark_rates_parametrisation_write_table["MuFitWin1"                  ].push_back (f_mu_fit_win1                      );
    inner_detector_dark_rates_parametrisation_write_table["ErrorMuFitWin1"             ].push_back (f_error_mu_fit_win1                );
    inner_detector_dark_rates_parametrisation_write_table["ExpConstWin1"               ].push_back (f_exp_const_win1                   );
    inner_detector_dark_rates_parametrisation_write_table["ErrorExpConstWin1"          ].push_back (f_error_exp_const_win1             );
    inner_detector_dark_rates_parametrisation_write_table["TauWin1"                    ].push_back (f_tau_win1                         );
    inner_detector_dark_rates_parametrisation_write_table["ErrorTauWin1"               ].push_back (f_error_tau_win1                   );
    inner_detector_dark_rates_parametrisation_write_table["PoisConstWin2"              ].push_back (f_pois_const_win2                  );
    inner_detector_dark_rates_parametrisation_write_table["ErrorPoisConstWin2"         ].push_back (f_error_pois_const_win2            );
    inner_detector_dark_rates_parametrisation_write_table["MuFitWin2"                  ].push_back (f_mu_fit_win2                      );
    inner_detector_dark_rates_parametrisation_write_table["ErrorMuFitWin2"             ].push_back (f_error_mu_fit_win2                );
    inner_detector_dark_rates_parametrisation_write_table["ExpConstWin2"               ].push_back (f_exp_const_win2                   );
    inner_detector_dark_rates_parametrisation_write_table["ErrorExpConstWin2"          ].push_back (f_error_exp_const_win2             );
    inner_detector_dark_rates_parametrisation_write_table["TauWin2"                    ].push_back (f_tau_win2                         );
    inner_detector_dark_rates_parametrisation_write_table["ErrorTauWin2"               ].push_back (f_error_tau_win2                   );

    bx_dbi::column_names index_columns;
    index_columns.push_back ("RunNumber");
    index_columns.push_back ("EchidnaCycle");

    bx_dbi::get ()->insert (bx_dbi::bx_calib, "\"InnerDetectorDarkRateParametrisation\"", inner_detector_dark_rates_parametrisation_write_table, index_columns);
 }
}

void db_run::flush () {
  if (!((detector_interface::get()->is_laben_enabled() && !b_write_laben_precalib) || 
	(detector_interface::get()->is_muon_enabled() && !b_write_muon_precalib))) {
    m_write_laben_precalib         ();
    m_write_laben_precalib_quality ();
    m_write_muon_precalib          ();
  }
  m_write_laben_laser_calib        ();
  m_write_laben_tt1_calib          ();
  m_write_muon_laser_calib         ();
  m_write_laben_dark_rates         ();
  m_write_muon_dark_rates          ();
  m_write_trigger_parameters       ();
  m_write_laben_electronic_channel (); 
  m_write_muon_electronic_channel  (); 
  m_write_muon_alignment           (); 
  m_write_disabled_channels	   ();
  m_write_laben_dark_rates_parametrisation();
} 

#endif

/*  
 *  $Log: db_run.cc,v $
 *  Revision 1.68  2015/07/24 12:37:27  ilia.drachnev
 *  added effective QE vectors
 *
 *  Revision 1.67  2013/01/17 17:03:23  misiaszek
 *  Method m_read_laben_dark_rates_parametrisation added
 *
 *  Revision 1.66  2013-01-16 22:47:14  misiaszek
 *  Visitors & getters for InnerDetectorDarkRateParametrisation added
 *
 *  Revision 1.65  2011-04-04 09:33:21  razeto
 *  Added disabled_channel_cycle option to choose the disabled channels list from a specified cycle
 *
 *  Revision 1.64  2011-03-21 11:22:08  razeto
 *  Choose the lowest cycle
 *
 *  Revision 1.63  2011-03-09 13:14:30  razeto
 *  Bugfix: add cycle to index column
 *
 *  Revision 1.62  2011-03-07 17:35:00  razeto
 *  Added cycle index in DisabledChannels
 *
 *  Revision 1.61  2009-01-27 22:57:06  razeto
 *  Fixed name for bx_calib_laben_dark_rates
 *
 *  Revision 1.60  2008-11-27 11:39:37  ludhova
 *  correct acl for tt1_rms
 *
 *  Revision 1.59  2008-11-26 17:57:13  ludhova
 *  debug, bad column name
 *
 *  Revision 1.58  2008-11-26 17:37:46  ludhova
 *  change name of internal variable of m_read_laben_tt1_calib
 *
 *  Revision 1.57  2008-11-26 13:53:42  ludhova
 *  mean, rms and P0 columns in NeutrinoPmtCalibration table
 *
 *  Revision 1.56  2008-11-20 12:18:59  ludhova
 *  debug in m_read_laben_tt1_calib
 *
 *  Revision 1.55  2008-11-17 14:46:29  ludhova
 *  visitors for NeutrinoPmtCalibration DB table
 *
 *  Revision 1.54  2008-10-14 13:23:50  ddangelo
 *  module bx_calib_muon_channel replaced by bx_calib_muon_time_align and bx_calib_muon_charge_peak
 *
 *  Revision 1.53  2008-09-26 13:34:20  ddangelo
 *  fixed timestamp handeling in write muon alignmente method
 *
 *  Revision 1.52  2008-09-25 11:45:40  ludhova
 *  visitors for DisabledChannels
 *
 *  Revision 1.51  2008-09-24 17:13:26  razeto
 *  Added disabled_channels runtime map
 *
 *  Revision 1.50  2008-09-04 15:32:26  ddangelo
 *  added a getter
 *  fixed a typo
 *
 *  Revision 1.49  2008-08-29 16:55:22  ddangelo
 *  debugging
 *
 *  Revision 1.48  2008-08-26 18:05:07  ddangelo
 *  muon alignement: times implemented as timestamps
 *
 *  Revision 1.47  2008-08-25 16:37:58  ddangelo
 *  muon alignment flagging implemented. varaibles, sgetters, read/write methods, flags.
 *  to be tested
 *
 *  Revision 1.46  2008-08-22 17:15:56  ddangelo
 *  restored db_channel_muon filling
 *
 *  Revision 1.45  2008-08-11 12:41:47  ddangelo
 *  added visitors for muon electronics calibration (multiplicity only for the moment)
 *  data, s/getters, read/write methods, initializiation, control flags and flushing
 *  b_electronic_channel_present modfied to b_laben_electronic_channel_present to avoid conflicts
 *
 *  Revision 1.44  2008-06-24 15:36:18  ddangelo
 *  fixed a module's name in acl settings
 *
 *  Revision 1.43  2008-06-20 16:25:06  razeto
 *  put parentheses to avoid confusion
 *
 *  Revision 1.42  2007-11-09 18:48:48  razeto
 *  Typo fixed
 *
 *  Revision 1.41  2007-11-09 18:38:57  razeto
 *  Added disabled pmts (thanks to Marcin for SQL support)
 *
 *  Revision 1.40  2007-06-03 16:13:30  razeto
 *  Do not keep db connection opened (else bxdb has lots of problems)
 *
 *  Revision 1.39  2007-05-30 13:48:55  ludhova
 *  correct acl setting for write_laben_precalib
 *
 *  Revision 1.38  2007-05-25 11:16:09  razeto
 *  bx_calib_laben_decoding will validate precalibrations
 *
 *  Revision 1.37  2007-02-28 16:17:22  ddangelo
 *  adding a check to align ID and OD precalib writing
 *
 *  Revision 1.36  2007-02-23 19:21:33  ddangelo
 *  visitor upgraded to includ muon dark rates data. both r/w implemented.
 *  some other variables renamed to avoid name clashes: interface unbroken.
 *
 *  Revision 1.35  2007/02/19 15:12:06  razeto
 *  Added new run types
 *
 *  Revision 1.34  2007/02/09 17:14:27  ludhova
 *  removed reading of previous run for precalib quality
 *
 *  Revision 1.33  2007-02-09 15:54:22  ludhova
 *  in case of missing precalibrationm m_read_laben_precalib_quality reads previous run
 *
 *  Revision 1.32  2007-01-30 11:03:07  razeto
 *  When doing precalib alway start from scratch
 *
 *  Revision 1.31  2006/12/06 14:15:16  ludhova
 *  Correct names for columns in table InnerDetectorDarkRate in write method
 *
 *  Revision 1.30  2006-11-28 13:11:29  razeto
 *  Fixed table name
 *
 *  Revision 1.29  2006/11/22 13:38:26  ludhova
 *  bugs in settings visitors corrected
 *
 *  Revision 1.28  2006-11-17 11:42:36  razeto
 *  Added water run type
 *
 *  Revision 1.27  2006/09/11 20:55:24  ddangelo
 *  corrected acl for muon calib modules
 *
 *  Revision 1.26  2006/09/07 14:11:41  ddangelo
 *  added muon laser calibration parameters, sgetters and acl.
 *
 *  Revision 1.25  2006/08/30 09:43:09  ludhova
 *  new visitors
 *
 *  Revision 1.24  2006/07/19 10:29:23  dfranco
 *  Added visitors for dark rates and detector efficiency
 *
 *  Revision 1.23  2006/07/18 14:48:20  razeto
 *  Added support for PrecalibDecodingQuality writing
 *
 *  Revision 1.22  2006/07/18 10:42:26  razeto
 *  Merged all reading code from Ludhova
 *
 *  Revision 1.21  2006/07/18 09:56:42  dfranco
 *  Reorganization of the methods
 *
 *  Revision 1.20  2006/06/21 12:36:43  ludhova
 *  added precalibration_quality (missing write)
 *
 *  Revision 1.19  2006/01/25 13:24:01  misiaszek
 *  Moved from cmap to simple map (to work with root)
 *
 *  Revision 1.18  2005/05/05 17:12:34  monzani
 *  Changed the names of the visitors for laser calibrations (Laben specified).
 *
 *  Revision 1.17  2005/05/05 10:06:01  monzani
 *  DB update from laser calibration modules added.
 *
 *  Revision 1.16  2004/11/30 14:30:30  razeto
 *  Added calib profile reading from database
 *
 *  Revision 1.15  2004/11/26 17:42:56  razeto
 *  Added muon precalib handling
 *
 *  Revision 1.14  2004/11/26 15:25:10  razeto
 *  Thanks to Davide who so kindly underlined a syntactical error, I fixed the "Maintainer" word
 *
 *  Revision 1.13  2004/11/26 14:01:20  razeto
 *  Added Mantainer field
 *
 *  Revision 1.12  2004/11/24 13:06:03  razeto
 *  Upgraded to the new cmap ctor with name
 *
 *  Revision 1.11  2004/11/16 16:01:38  razeto
 *  Upgraded database writing with transaction. Overrite support almost done
 *
 *  Revision 1.10  2004/10/19 18:41:58  razeto
 *  Integrated visitors writing in the framework
 *
 *  Revision 1.9  2004/10/19 16:30:41  razeto
 *  Added laben predecoding reading.
 *  Added predecoding writing to the database (still experimental).
 *
 *  Revision 1.8  2004/08/14 21:21:08  razeto
 *  Added some fields used from calibrations (and some improvements)
 *
 *  Revision 1.7  2004/08/05 09:05:20  razeto
 *  Used postgres EXTRACT command to convert date to time_t
 *
 *  Revision 1.6  2004/06/10 10:06:48  razeto
 *  Added some run type, and changed an error condition
 *
 *  Revision 1.5  2004/05/26 09:21:59  razeto
 *  Upgraded to last version of vdt
 *
 *  Revision 1.4  2004/05/25 17:15:47  razeto
 *  Upgraded to use name in vdt::get_*
 *
 *  Revision 1.3  2004/05/18 14:26:59  razeto
 *  Updated
 *
 *  Revision 1.2  2004/04/26 13:48:29  razeto
 *  Added db_acl to check set calls for calling module to have the right privileges.
 *  Fixed the names of set/get methods.
 *  Modifications decided at the software Paris meeting.
 *
 *  Revision 1.1  2004/04/09 08:05:00  razeto
 *  Added
 *
 */
