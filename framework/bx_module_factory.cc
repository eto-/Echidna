/* BOREXINO Reconstruction program
 *
 * Author: Alessandro Razeto <Alessandro.Razeto@ge.infn.it>
 * 	   Marco Pallavicini <pallas@ge.infn.it>
 * Maintainer: Alessandro Razeto <Alessandro.Razeto@ge.infn.it>
 *
 * $Id: bx_module_factory.cc,v 1.96 2015/08/06 12:24:38 koun Exp $
 *
 * Implementation of bx_module_factory
 *
 */
#include <fstream>
#include <iostream>
#include <iterator>
#include <algorithm>

#include "bx_module_factory.hh"
#include "bx_base_module.hh"
#include "bx_test_module.hh"
#include "bx_trigger_decoder.hh"
#include "bx_reader.hh"
#include "bx_writer.hh"
#include "bx_muon_decoder.hh"
#include "bx_muon_findcluster.hh"
#include "bx_muon_tracker.hh"
#include "bx_laben_energy_tracker.hh"
#include "bx_laben_tof_tracker.hh"
#include "bx_global_tracker.hh"
#include "bx_reader.hh"
#include "bx_writer.hh"
#include "bx_precalib_laben_adc.hh"
#include "bx_precalib_laben_gray.hh"
#include "bx_precalib_laben_phases.hh"
#include "bx_precalib_laben_d80.hh"
#include "bx_precalib_laben_check_tdc.hh"
#include "bx_precalib_muon_gate.hh"
#include "bx_precalib_muon_findpulse.hh"
#include "bx_precalib_muon_pedestals.hh"
#include "bx_laben_findcluster.hh"
#include "bx_laben_dark_rates.hh"
#include "bx_laben_decoder.hh"
#include "bx_detector_monitor.hh"
#include "bx_calib_laben_electronics.hh"
#include "bx_calib_muon_electronics.hh"
#include "bx_calib_laben_crate_delay.hh"
#include "bx_calib_gate.hh"
#include "bx_baricentrator.hh"
#include "bx_calib_laben_dark_rates.hh"
#include "bx_calib_fiber_bundle.hh"
#include "bx_calib_laser_transparency.hh"
#include "bx_calib_laben_time_align.hh"
#include "bx_calib_laben_charge_peak.hh"
#include "bx_calib_laben_charge_tt1.hh"
#include "bx_calib_muon_time_align.hh"
#include "bx_calib_muon_charge_peak.hh"
#include "bx_calib_muon_dark_rates.hh"
#include "bx_calib_muon_pulser.hh"
#include "bx_calib_muon_alignment.hh"
#include "bx_calib_monitor.hh"
#include "bx_calib_laben_decoding.hh"
#include "bx_position_reco_mi.hh"
#include "bx_position_reco_lngs.hh"
#include "bx_position_reco_noavg.hh"
#include "bx_splitting_filter.hh"
#include "bx_calib_laben_retrigger.hh"
#include "bx_position_reco_msk.hh"
#include "bx_position_reco_dbn.hh"
#include "bx_position_reco_mach4.hh"
#include "bx_position_reco.hh"
#include "bx_energy_reco_mc.hh"
#include "bx_energy_reco_lik.hh"
#include "bx_energy_reco_msk.hh"
#include "bx_cmt_tracker.hh"
#include "bx_pid_shape.hh"
#include "bx_pid_alphabeta.hh"
#include "bx_pid_positron.hh"
#include "bx_energy_reco_dbn.hh"
#include "bx_omon.hh"
#include "bx_filter_trigger.hh"
#include "bx_laben_raw_validator.hh"
#include "bx_v1731sys.hh"
#include "bx_pid_ab_mach4.hh"
#include "messenger.hh"
#include "bx_snews.hh"
#include "bx_new_charge_weight.hh"

bx_module_factory *bx_module_factory::me = 0;

// ctor
bx_module_factory::bx_module_factory () {
  modules.clear ();
}

bx_module_factory::~bx_module_factory () {
  std::for_each (modules.begin (), modules.end (), bx_base_module::delete_operation ());
  modules.clear ();
}

// singleton
bx_module_factory* bx_module_factory::get () {
  if (!me) 
    me = new bx_module_factory;

  return me;
}

bx_base_module *bx_module_factory::get_module (const std::string &module_name) {

  bx_base_module::bx_base_module_vector::iterator item = std::find_if (modules.begin (), modules.end (), bx_base_module::compare_name_operation (module_name));
  if (item != modules.end ()) return *item;
  
  bx_message msg(bx_message::critic, "bx_module_factory: ");
  msg << "asked for module " << module_name << " which is not present. Present modules:";
  std::copy (modules.begin (), modules.end (), std::ostream_iterator<bx_base_module *>(msg, "\n"));
  msg << dispatch;

  return *item;
}


bx_base_module *bx_module_factory::get_module (bx_base_module::module_role role) {
 bx_base_module::bx_base_module_vector::iterator item = std::find_if (modules.begin (), modules.end (), bx_base_module::compare_role_operation (role));
  if (item != modules.end ()) return *item;
  
  bx_message msg(bx_message::critic, "bx_module_factory: ");
  msg << "asked for module with role " << role << " which is not present. Present modules:";
  std::copy (modules.begin (), modules.end (), std::ostream_iterator<bx_base_module *>(msg, "\n"));
  msg << dispatch;

  return *item;
}

void bx_module_factory::delete_module (const std::string &module_name) {
  
  if (module_name == "ALL") {
    std::for_each (modules.begin (), modules.end (), bx_base_module::delete_operation ());
    modules.clear ();
  } else {
    bx_base_module::bx_base_module_vector::iterator item = std::find_if (modules.begin (), modules.end (), bx_base_module::compare_name_operation (module_name));
    if (item != modules.end ()) {
      delete *item;
      modules.erase (item);
    } else {
      bx_message msg(bx_message::critic, "bx_module_factory: ");
      msg << "asked for module " << module_name << " which is not present. Present modules:";
      std::copy (modules.begin (), modules.end (), std::ostream_iterator<bx_base_module *>(msg, "\n"));
      msg << dispatch;
    }
  }
}

void bx_module_factory::create_modules (const std::vector<std::string>& module_name_list) {
  bx_message msg(bx_message::critic, "bx_module_factory: ");

  if (modules.size ()) msg << "already inited" << dispatch;

  for (unsigned long i = 0; i < module_name_list.size (); i++) {
    std::string m_name = module_name_list[i];
    if (m_name == "bx_reader") {
      modules.push_back (new bx_reader);
    } else if (m_name == "bx_writer") {
      modules.push_back (new bx_writer);
    } else if (m_name == "bx_precalib_laben_adc") {
      modules.push_back (new bx_precalib_laben_adc);
    } else if (m_name == "bx_precalib_laben_gray") {
      modules.push_back (new bx_precalib_laben_gray);
    } else if (m_name == "bx_precalib_laben_phases") {
      modules.push_back (new bx_precalib_laben_phases);
    } else if (m_name == "bx_precalib_laben_d80") {
      modules.push_back (new bx_precalib_laben_d80);
    } else if (m_name == "bx_precalib_laben_check_tdc") {
      modules.push_back (new bx_precalib_laben_check_tdc);
    } else if (m_name == "bx_precalib_muon_findpulse") {
      modules.push_back (new bx_precalib_muon_findpulse);
    } else if (m_name == "bx_precalib_muon_pedestals") {
      modules.push_back (new bx_precalib_muon_pedestals);
    } else if (m_name == "bx_precalib_muon_gate") {
      modules.push_back (new bx_precalib_muon_gate);
    } else if (m_name == "bx_test_module") {
      modules.push_back (new bx_test_module);
    } else if (m_name == "bx_trigger_decoder") {
      modules.push_back (new bx_trigger_decoder);
    } else if (m_name == "bx_muon_decoder") {
      modules.push_back (new bx_muon_decoder);
    } else if (m_name == "bx_laben_decoder") {
      modules.push_back (new bx_laben_decoder);
    } else if (m_name == "bx_laben_findcluster") {
      modules.push_back (new bx_laben_findcluster);
    } else if (m_name == "bx_laben_dark_rates") {
      modules.push_back (new bx_laben_dark_rates);
    } else if (m_name == "bx_muon_findcluster") {
      modules.push_back (new bx_muon_findcluster);
    } else if (m_name == "bx_muon_tracker") {
      modules.push_back (new bx_muon_tracker);
    } else if (m_name == "bx_global_tracker") {
      modules.push_back (new bx_global_tracker);
    } else if (m_name == "bx_detector_monitor") {
      modules.push_back (new bx_detector_monitor);
    } else if (m_name == "bx_calib_laben_electronics") {
      modules.push_back (new bx_calib_laben_electronics);
    } else if (m_name == "bx_calib_muon_electronics") {
      modules.push_back (new bx_calib_muon_electronics);
    } else if (m_name == "bx_calib_laben_crate_delay") {
      modules.push_back (new bx_calib_laben_crate_delay);
    } else if (m_name == "bx_calib_gate") {
      modules.push_back (new bx_calib_gate);
    } else if (m_name == "bx_calib_laben_time_align") {
      modules.push_back (new bx_calib_laben_time_align);
    } else if (m_name == "bx_calib_laben_charge_peak") {
      modules.push_back (new bx_calib_laben_charge_peak);
    } else if (m_name == "bx_calib_laben_charge_tt1") {
      modules.push_back (new bx_calib_laben_charge_tt1);
    } else if (m_name == "bx_calib_muon_time_align") {
      modules.push_back (new bx_calib_muon_time_align);
    } else if (m_name == "bx_calib_muon_charge_peak") {
      modules.push_back (new bx_calib_muon_charge_peak);
    } else if (m_name == "bx_calib_muon_pulser") {
      modules.push_back (new bx_calib_muon_pulser);
    } else if (m_name == "bx_calib_muon_alignment") {
      modules.push_back (new bx_calib_muon_alignment);
    } else if (m_name == "bx_calib_muon_dark_rates") {
      modules.push_back (new bx_calib_muon_dark_rates);
    } else if (m_name == "bx_calib_laben_decoding") {
      modules.push_back (new bx_calib_laben_decoding);
    } else if (m_name == "bx_baricentrator") {
      modules.push_back (new bx_baricentrator);
    } else if (m_name == "bx_calib_laben_dark_rates") {
      modules.push_back (new bx_calib_laben_dark_rates);
    } else if (m_name == "bx_calib_fiber_bundle") {
      modules.push_back (new bx_calib_fiber_bundle); 
    } else if (m_name == "bx_calib_laser_transparency") {
      modules.push_back (new bx_calib_laser_transparency);        
    } else if (m_name == "bx_position_reco_mi") {
      modules.push_back (new bx_position_reco_mi);
    } else if (m_name == "bx_position_reco_lngs") {
      modules.push_back (new bx_position_reco_lngs);
    } else if (m_name == "bx_position_reco_noavg") {
      modules.push_back (new bx_position_reco_noavg);
    } else if (m_name == "bx_position_reco_mach4") {
      modules.push_back (new bx_position_reco_mach4);
    } else if (m_name == "bx_splitting_filter") {
      modules.push_back (new bx_splitting_filter);
    } else if (m_name == "bx_calib_laben_retrigger") {
      modules.push_back (new bx_calib_laben_retrigger);
    } else if (m_name == "bx_position_reco_msk") {
      modules.push_back (new bx_position_reco_msk);
    } else if (m_name == "bx_position_reco_dbn") {
      modules.push_back (new bx_position_reco_dbn);
    } else if (m_name == "bx_position_reco") {
      modules.push_back (new bx_position_reco);
    } else if (m_name == "bx_energy_reco_mc") {
      modules.push_back (new bx_energy_reco_mc);      
    } else if (m_name == "bx_energy_reco_lik") {
      modules.push_back (new bx_energy_reco_lik);      
    } else if (m_name == "bx_energy_reco_msk") {
      modules.push_back (new bx_energy_reco_msk);      
    } else if (m_name == "bx_energy_reco_dbn") {
      modules.push_back (new bx_energy_reco_dbn);      
    } else if (m_name == "bx_laben_energy_tracker") {
      modules.push_back (new bx_laben_energy_tracker);      
    } else if (m_name == "bx_laben_tof_tracker") {
      modules.push_back (new bx_laben_tof_tracker);      
    } else if (m_name == "bx_pid_alphabeta") {
      modules.push_back (new bx_pid_alphabeta);
    } else if (m_name == "bx_pid_shape") {
      modules.push_back (new bx_pid_shape);
    } else if (m_name == "bx_pid_positron") {
      modules.push_back (new bx_pid_positron);
    } else if (m_name == "bx_omon") {
      modules.push_back (new bx_omon);
    } else if (m_name == "bx_filter_trigger") {
      modules.push_back (new bx_filter_trigger);
    } else if (m_name == "bx_snews") {
      modules.push_back (new bx_snews);
    } else if (m_name == "bx_laben_raw_validator") {
      modules.push_back (new bx_laben_raw_validator);
    } else if (m_name == "bx_calib_monitor") {
      modules.push_back (new bx_calib_monitor);
    } else if (m_name == "bx_v1731sys") {
      modules.push_back (new bx_v1731sys);
    } else if (m_name == "bx_pid_ab_mach4") {
      modules.push_back (new bx_pid_ab_mach4);
    } else if (m_name == "bx_cmt_tracker") {
      modules.push_back (new bx_cmt_tracker);
    } else if (m_name == "bx_new_charge_weight") {
      modules.push_back (new bx_new_charge_weight);
    } else {
      msg << "unknown module name: " << m_name << dispatch;
    }
  }

  // Sort the array
  std::sort (modules.begin (), modules.end (), bx_base_module::sort_operation ());

  bx_message info(bx_message::info, "bx_module_factory: module list:");
  for (unsigned long i = 0; i < modules.size (); i++) {
    if (modules[i]->is_enabled ()) {
      info << std::endl << "\t" <<  modules[i]->get_name () << " role " << modules[i]->get_role ();
      info << " priority " << modules[i]->get_priority () << " enabled";
    }
  }
  info << dispatch;
  
  // A consinstency test for the modules array: reader have to be present in
  // only 1 enabled module.
  m_check_single_presence (bx_base_module::reader);
}


void bx_module_factory::m_check_single_presence (bx_base_module::module_role role) {
  bx_message msg(bx_message::critic, "bx_module_factory: ");

  // compare_role_operation only takes in account enabled modules.
  bx_base_module::bx_base_module_vector::iterator item1 = std::find_if (modules.begin (), modules.end (), bx_base_module::compare_role_operation (role));
  if (item1 == modules.end ()) msg << "no module for role " << role << ", aborting" << dispatch;

  bx_base_module::bx_base_module_vector::iterator item2 = std::find_if (item1 + 1, modules.end (), bx_base_module::compare_role_operation (role));
  if (item2 != modules.end ()) msg << "several modules for role " << role << " (" << (*item1)->get_name () << ", " << (*item2)->get_name () << "), aborting" << dispatch;
}
  
/*
 * $Log: bx_module_factory.cc,v $
 * Revision 1.96  2015/08/06 12:24:38  koun
 * added bx_new_charge_weight
 *
 * Revision 1.95  2015/07/29 09:48:58  ilia.drachnev
 * *** empty log message ***
 *
 * Revision 1.94  2015/07/29 09:38:38  ilia.drachnev
 * added noavg pos reco
 *
 * Revision 1.93  2015/07/27 09:39:03  misiaszek
 * module bx_laben_cluster_noavg removed
 *
 * Revision 1.92  2015/07/14 10:12:45  misiaszek
 * bx_laben_cluster_noavg module added to factory
 *
 * Revision 1.91  2013/01/10 18:13:33  ludhova
 * new modul for DR parametrisation added
 *
 * Revision 1.90  2011-02-18 17:10:05  ddangelo
 * major code cleanup: removed fadc throughout the program
 *
 * Revision 1.89  2011-02-18 15:30:16  ddangelo
 * added bx_pid_positron module
 *
 * Revision 1.88  2010-05-21 12:33:34  ddangelo
 * bx_laben_tracker renamed as bx_laben_energy_tracker
 * added bx_laben_tof_tracker and bx_cmt_tracker
 *
 * Revision 1.87  2009-10-23 14:00:03  koshio
 * Add the lngs postion reconstruction
 *
 * Revision 1.86  2009-07-16 15:52:12  razeto
 * Added mach4 a/b
 *
 * Revision 1.85  2009-07-16 12:40:21  razeto
 * Added Mach4 position reco
 *
 * Revision 1.84  2008-11-13 09:06:07  ludhova
 * bx_calib_laben_charge_tt1 module
 *
 * Revision 1.83  2008-10-25 10:00:14  ddangelo
 * added bx_precalib_muon_gate module
 *
 * Revision 1.82  2008-10-23 09:15:53  dicienzo
 * Added bx_snews
 *
 * Revision 1.81  2008-10-14 15:26:51  wurm
 * debugging
 *
 * Revision 1.80  2008-10-14 13:23:43  ddangelo
 * module bx_calib_muon_channel replaced by bx_calib_muon_time_align and bx_calib_muon_charge_peak
 *
 * Revision 1.79  2008-10-01 15:56:38  guardi
 * bx_calib_monitor added
 *
 * Revision 1.78  2008-08-20 16:21:27  ddangelo
 * module bx_calib_dark_rates renamed as bx_calib_laben_dark_rates
 *
 * Revision 1.77  2008-08-11 12:49:38  ddangelo
 * added bx_calib_muon_electronics module
 *
 * Revision 1.76  2008-05-13 12:42:10  ddangelo
 * added bx_calib_muon_alignment
 *
 * Revision 1.75  2008-04-29 14:05:34  ddangelo
 * debugging
 *
 * Revision 1.74  2008-04-29 13:45:09  ddangelo
 * added bx_global_tracker
 *
 * Revision 1.73  2008-04-02 13:40:31  pallas
 * Adding new module hook in factory
 *
 * Revision 1.72  2008-02-21 17:17:29  ddangelo
 * added laben tracker
 *
 * Revision 1.71  2007-10-25 17:22:39  razeto
 * Added bx_energy_reco_msk
 *
 * Revision 1.70  2007-07-10 16:31:17  razeto
 * Added energy_dbn
 *
 * Revision 1.69  2007-05-21 15:05:53  ddangelo
 * added module bx_calib_muon_pulser
 *
 * Revision 1.68  2007-05-03 15:47:59  ddangelo
 * added bx_muon_tracker module
 *
 * Revision 1.67  2007-02-22 19:57:41  ddangelo
 * added bx_calib_muon_dark_rates
 *
 * Revision 1.66  2007/02/21 15:48:37  ddangelo
 * added bx_muon_findcluster
 *
 * Revision 1.65  2006/10/23 14:41:17  razeto
 * New module bx_laben_raw_validator to validate event using laben raw hit flags
 *
 * Revision 1.64  2006/10/19 11:50:49  ludhova
 * Added module
 *
 * Revision 1.63  2006-09-11 14:12:38  ddangelo
 * bx_calib_muon_charge_peak and bx_calib_muon_time_align
 * replaced by
 * bx_calib_muon_channel
 *
 * Revision 1.62  2006/08/23 15:56:18  ddangelo
 * added 2 muon calibration modules
 *
 * Revision 1.61  2006/08/21 11:14:44  razeto
 * Added new bx_filter_trigger module
 *
 * Revision 1.60  2006/07/17 13:49:25  razeto
 * Added omon for Andrew
 *
 * Revision 1.59  2006/05/10 12:14:39  razeto
 * Removed some useless include
 *
 * Revision 1.58  2006/03/08 18:01:16  razeto
 * Added bx_calib_laser_transparency
 *
 * Revision 1.57  2006/03/08 17:59:38  pepem
 * Adding bx_pid_alphabeta module
 *
 * Revision 1.56  2006/03/02 16:37:42  razeto
 * Added new bx_calib_laben_retrigger
 *
 * Revision 1.55  2006/01/11 11:54:20  dfranco
 * added a new module: bx_calib_laben_decoding.hh
 *
 * Revision 1.54  2005/12/07 21:36:21  misiaszek
 *
 * bx_energy_reco_lik module added to factory
 *
 * Revision 1.53  2005/12/01 12:47:57  misiaszek
 *
 * bx_energy_reco_mc module added (auth from Alessandro Razeto)
 *
 * Revision 1.52  2005/09/20 13:49:35  razeto
 * Added bx_pid_shape
 *
 * Revision 1.51  2005/07/13 12:31:51  razeto
 * Merged the the 2 laben clustering modules: now find_cluster_time is
 * a subroutine of the bx_laben_findcluster module
 *
 * Revision 1.50  2005/06/27 15:14:04  razeto
 * Added bx_position_reco_dbn to framework
 *
 * Revision 1.49  2005/06/22 14:06:43  razeto
 * Added new module
 *
 * Revision 1.48  2005/06/20 14:17:27  razeto
 * Added position reco
 *
 * Revision 1.47  2005/05/05 09:26:01  razeto
 * Added laben laser calibration modules
 *
 * Revision 1.46  2005/03/04 11:50:22  razeto
 * Removed a useless include
 *
 * Revision 1.45  2004/12/19 00:02:38  razeto
 * Modified to list only the enabled modules, since the full configuration is already logged from CM
 *
 * Revision 1.44  2004/12/13 15:39:01  razeto
 * Added Moscow spatial reconstruction
 *
 * Revision 1.43  2004/12/09 17:21:35  razeto
 * Added splitting
 *
 * Revision 1.42  2004/11/26 15:25:10  razeto
 * Thanks to Davide who so kindly underlined a syntactical error, I fixed the "Maintainer" word
 *
 * Revision 1.41  2004/11/26 14:01:19  razeto
 * Added Mantainer field
 *
 * Revision 1.40  2004/11/24 13:29:44  razeto
 * Changed name of position reco to bx_position_reco_mi
 *
 * Revision 1.39  2004/11/24 13:00:07  razeto
 * Added bx_fadc_findcluster module
 *
 * Revision 1.38  2004/11/17 11:53:55  razeto
 * Added position reco
 *
 * Revision 1.37  2004/10/05 13:51:43  razeto
 * Updated some module names
 *
 * Revision 1.36  2004/10/01 10:28:19  razeto
 * Added bx_calib_fiber_bundle
 *
 * Revision 1.35  2004/09/30 14:35:36  razeto
 * Added delete_module to module_factory to destroy one or all modules
 *
 * Revision 1.34  2004/09/28 12:28:23  razeto
 * added bx_laben_detector_status
 *
 * Revision 1.33  2004/09/23 10:09:24  razeto
 * Changed internal reader name (to follow writer)
 *
 * Revision 1.32  2004/09/17 13:35:42  razeto
 * Added bx_baricentrator
 *
 * Revision 1.31  2004/08/10 11:14:03  razeto
 * Added a module for calculating the time differences from crate to crate
 *
 * Revision 1.30  2004/07/22 10:16:29  razeto
 * Added bx_laben_electronics_status
 *
 * Revision 1.29  2004/07/22 10:13:01  razeto
 * Jumbo patch for the configuration manager applied, read the documentation for usage
 *
 * Revision 1.28  2004/07/16 07:34:38  razeto
 * Add fadc decoder
 *
 * Revision 1.27  2004/07/12 10:05:30  razeto
 * Added some clustering modules
 *
 * Revision 1.26  2004/06/10 12:40:19  razeto
 * Added some needed includes
 *
 * Revision 1.25  2004/06/07 10:59:20  razeto
 * Upgraded to handle multiple writer at the same time
 *
 * Revision 1.24  2004/06/01 15:36:14  razeto
 * Upgrade to full support vdt vectors
 *
 * Revision 1.23  2004/05/31 14:59:03  razeto
 * Updated to new method is_enabled in base module
 *
 * Revision 1.22  2004/05/21 08:41:59  razeto
 * Added laben decoder
 *
 * Revision 1.21  2004/05/18 14:54:27  razeto
 * Added few modules
 *
 * Revision 1.20  2004/04/27 16:58:00  ddangelo
 * removed old includes
 *
 * Revision 1.19  2004/04/27 09:50:43  ddangelo
 * test module precalibrator removed.
 * module bx_trigger_decode_module renamed to bx_trigger_decoder.
 *
 * Revision 1.18  2004/04/24 17:36:45  razeto
 * Added a check if the parameter already exist in setting (added an init method too)
 *
 * Revision 1.17  2004/04/20 16:08:02  ddangelo
 * added handling for a few muon models
 *
 * Revision 1.16  2004/04/18 10:30:14  razeto
 * Fixed to compile wit g++ 2.95 present on the cluster
 *
 * Revision 1.15  2004/04/12 16:01:15  razeto
 * Added the precalib modules and a message
 *
 * Revision 1.14  2004/04/09 07:50:19  razeto
 * Added some more messages
 *
 * Revision 1.13  2004/04/06 12:43:55  razeto
 * Added module role to base_module, and the code to handle everywhere. Some other minor fixes
 *
 * Revision 1.12  2004/04/03 09:23:24  razeto
 * Fixed a bug. Corrected the test module name
 *
 * Revision 1.11  2004/04/01 12:08:53  razeto
 * Moved the base_module_operation operators in the bx_base_module class
 * scope; moved the bx_base_module_vector typedef in the bx_base_module
 * class.
 * Added a event delete to the frame.
 *
 * Revision 1.10  2004/03/26 16:36:11  razeto
 * Introduced bx_options; a lot of code modified to read options and parameters.
 * Bx_event_reader interface changed: now there is a standard constructor and
 * the file opening is done at begin using the parameters for the file name.
 *
 * Revision 1.9  2004/03/26 16:31:14  razeto
 * Fixed a bug: reader, writer and precalibrator were pointers, but the could
 * be even not initialized since the initialization were left to a procedure
 * which would allocate them if they where present on the modules.cfs file.
 * This could not be true.
 * Now these special modules are in a dedicated vector which however has
 * the same paradigm of the standard module vector.
 * This solves the problem; maybe a role for the module could be introduced
 * in future, allowing to have just one vector.
 *
 * Revision 1.8  2004/03/24 16:22:08  razeto
 * Moved reader,writer,precalib creation from framework to factory
 *
 * Revision 1.7  2004/03/24 14:23:35  razeto
 * Moved module list ownership from framewrok to module_factory, this allows:
 *   - a better syntax for bx_options which can store parameters even if the
 * framework is not created
 *   - inter module comunications
 * Since the module factory was already a singleton this does not affect
 * the modules array lifecycle.
 *
 * Revision 1.6  2004/03/23 14:26:47  pallas
 * Debugging trigger module
 * GPS clock OK
 *
 * Revision 1.5  2004/03/22 14:29:25  razeto
 * Some cosmetic changes
 *
 * Revision 1.4  2004/03/21 20:19:17  razeto
 * Added std:: to string
 *
 * Revision 1.3  2004/03/20 18:46:45  pallas
 * debugging
 *
 * Revision 1.2  2004/03/20 17:50:54  pallas
 * Debugging
 *
 * Revision 1.1  2004/03/20 17:39:45  pallas
 * First release
 *
 * Revision 1.1  2004/03/20 14:16:47  razeto
 * Added sources
 *
 *
 */
