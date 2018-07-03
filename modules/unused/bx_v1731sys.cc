/* BOREXINO Reconstruction program
 *
 * Author: Stefano Davini (stefano.davini@ge.infn.it)
 * Maintainer: Stefano Davini (stefano.davini@ge.infn.it)
 *
 * $Id: bx_v1731sys.cc,v 1.19 2009/11/20 10:55:24 davini Exp $
 *
 * Module for analysis of V1731 data (28 channels system for neutron detection and other)
 * 
*/

#include "bx_v1731sys.hh"
#include "bx_echidna_event.hh"
#include "bx_neutron_event.hh"
#include "bx_dbi.hh"
#include "db_profile.hh"
#include <stdlib.h>
#include <iomanip>
#include "zlib.h"

bx_v1731sys::bx_v1731sys() : bx_base_module("bx_v1731sys", bx_base_module::main_loop) {
 }


void bx_v1731sys::begin () {
  transfer_opt.s_filepath = "/bxstorage/neutron/";  /* default value */
  transfer_opt.s_wget_cmd = "wget -nv --timeout=60 -O- http://bxmaster-data.lngs.infn.it";
  transfer_opt.transfer_mode = v1731_event_correlator::http;

  if (char *p = ::getenv("NEUTRON_FILE_PATH")) { /* environment variable for file path OUTSIDE bxmaster */
    transfer_opt.s_filepath = p;
    transfer_opt.transfer_mode = v1731_event_correlator::file;
    get_message(bx_message::info)<<"File path for neutron files set to " << transfer_opt.s_filepath << dispatch;
  }

  b_map_created = false; /* status var; changed in void m_load_file_map() */
  m_load_file_map();     /* load map of aviable file*/

  v1731event_index.set_transfer_param(transfer_opt);
  
  currNfile = 0;          /* current N (raw_data) gzipped file; opened/loaded/read/decoded in v1731_event_getter*/
  p_connect = 0;          /* current pipe */
  s_curr_file_name = "";  /* current file_name and postion are updated in v1731_event_getter*/
  i4_curr_file_position = 0;

  last_neutron_gpsstime = 0;
  last_neutron_file_position  = 0;
  last_neutron_laben_nclusters= 0;

}


bx_echidna_event* bx_v1731sys::doit (bx_echidna_event* ev) {
  
  bx_neutron_event& ew = dynamic_cast<bx_neutron_event&> (ev->get_neutron());
  m_clear_ev(ew);

  if (ev->get_run_number() < 7838) return ev;  /* runs before NeutronDAQ operative */
  if (ev->get_run_number() < 8080) return ev;   /* some runs in june 08 have daq desincronized */
  if ((ev->get_run_number()>=8080) && (ev->get_run_number()<=8300)) return ev; /* runs in july 08: daq desincrinized*/
  
  //const bool b_skip_after_run8504 = (bool) get_parameter("skip_after_run8504").get_int();
  if ((ev->get_run_number()> 8504) && (ev->get_run_number()<9313)) return ev; 

  if (!b_map_created) return ev;

  /* trigger128 event selection */
  if (!(ev->get_trigger().is_neutron())) return ev;
  
  uint32_t u4_sec, u4_nsec;
  ev->get_trigger().get_gps_time(u4_sec,u4_nsec);    /* get event gpstime */

  /* add new event's file to the index (if file exists, and not already indexed) */
  v1731event_index.new_event(*map_time_file, u4_sec);
  ew.b_is_enabled = v1731event_index.daq_active(u4_sec);

  /* search most appropriate event in these indexed */
  v1731_event_correlator::event_finder find_event(v1731event_index, u4_sec);
  //find_event.m_show();  /* debug mode */

  if (!find_event.get_found()) return ev;

  /* candidate v1731 events have been found */

  const bool b_trgid_acquired = (ev->get_run_number() >= 8346 ? true : false);
  
  if (b_trgid_acquired){
    std::vector <v1731_event::v1731event> v1731evs;
    v1731evs.reserve(find_event.get_n_founds());
    unsigned u4_min_dtrg_elem = 0;
    int32_t      i4_min_dtrg_value= 500;
    const unsigned u4_trgid = (0xffff & ev->get_event_number());
    for (unsigned u4=0; u4<find_event.get_n_founds(); ++u4){
      v1731_event::v1731event v1731ev_tmp;
      v1731_event::v1731event_decoder v1731ev_decoder(find_event, v1731ev_tmp, u4, currNfile, p_connect, s_curr_file_name, i4_curr_file_position);
      v1731evs.push_back(v1731ev_tmp);
      const int32_t i4_trg_diff = static_cast <int32_t> (u4_trgid - v1731ev_tmp.get_trigger_id());
      if ((i4_trg_diff<=1) && (abs(i4_trg_diff) <= abs(i4_min_dtrg_value))){
	i4_min_dtrg_value = i4_trg_diff;
	u4_min_dtrg_elem  = u4;
      }
    }
    if ((i4_min_dtrg_value>1) || (i4_min_dtrg_value<-7)) return ev;
    /* association done! */
    // get_message(bx_message::debug)<<" dtrg = "<<i4_min_dtrg_value<<dispatch;
    ew.b_is_associated = true;
    v1731evs.at(u4_min_dtrg_elem).clusterize(100, 130, 800, 800);
    v1731_neutron_analyzer::neutron_analyzer neutron_analysis(v1731evs[u4_min_dtrg_elem], true);
    if (neutron_analysis.get_has_neutron()) m_write_ev(ew, neutron_analysis);
    return ev;
  }

  /* procedure for runs without trgid acquired*/

  ew.b_is_associated = true;
  for (unsigned u4f=0; u4f<find_event.get_n_founds(); ++u4f){

    if ((last_neutron_gpsstime == find_event.get_v1731_timer(u4f)) && 
	(last_neutron_file_position == find_event.get_position(u4f)) &&
	(last_neutron_laben_nclusters>0)) continue;    /* avoids some wrong associations*/

    v1731_event::v1731event v1731ev;
    v1731_event::v1731event_decoder v1731ev_decoder(find_event, v1731ev, u4f, currNfile, p_connect, s_curr_file_name, i4_curr_file_position); //reads event, fills v1731ev
    
    if (v1731ev_decoder.error_occurred()){
      get_message(bx_message::error)<<"Error reading/loading "<<find_event.get_file_name(u4f)<<dispatch;
      break;
    }
    else {  
      v1731_neutron_analyzer::neutron_analyzer neutron_analysis(v1731ev, true);
      if (neutron_analysis.get_has_neutron()){
	last_neutron_gpsstime        = find_event.get_v1731_timer(u4f);
	last_neutron_file_position   = find_event.get_position(u4f);
	last_neutron_laben_nclusters = ev->get_laben().get_nclusters(); 
	m_write_ev(ew, neutron_analysis);
	break;
      }
    }
  }
  
  return ev;  
}

  
void bx_v1731sys::end () {
  if (b_map_created) delete map_time_file;
  if (currNfile) gzclose(currNfile);
  if (p_connect) pclose(p_connect);
}


/*********** private methods ****************/


void bx_v1731sys::m_load_file_map() {
    
  /* read file list of directory and store map with start and end times */
 
  std::vector <std::string> s_file_names;  // N_.dat.gz names
  std::vector <std::string> s_aux_names;   // T_.dat names
  s_file_names.reserve(1024);
  s_aux_names.reserve(1024);
  bool b_error_occurred = false;
  v1731_event_correlator::make_filename_vectors(s_file_names, s_aux_names, transfer_opt, b_error_occurred);

  if (!b_error_occurred){ /* create map */
    using v1731_event_correlator::mapper_time_file;
    map_time_file = new mapper_time_file(s_file_names, s_aux_names);
    b_map_created = true;
    get_message(bx_message::info)<<"Time-File map for Neutron system created;"<<dispatch;
  }
  else if (b_error_occurred){
    get_message(bx_message::warn)<<"Time-File map for Neutron system not created;"<<dispatch;
  }
}

void bx_v1731sys::m_clear_ev(bx_neutron_event& ew){
  ew.b_is_enabled = false;
  ew.b_is_associated = false;
  ew.i4_n_neutrons = 0;
  ew.pulses.clear();
}

void bx_v1731sys::m_write_ev(bx_neutron_event& ew, const v1731_neutron_analyzer::neutron_analyzer& neutron_analysis){
  const unsigned u4_n_neutrons = neutron_analysis.get_n_neutrons(); // in fact it's number of pulses 
  ew.i4_n_neutrons = u4_n_neutrons;  // now canditates neutron are all pulses; to be changed in future;
  bx_neutron_pulse tmp_neutron_pulse;
  for (unsigned u4n=0; u4n<u4_n_neutrons; ++u4n){
    tmp_neutron_pulse.f4_peak_time = neutron_analysis.get_time(u4n);
    tmp_neutron_pulse.f4_amplitude = neutron_analysis.get_peak(u4n);
    tmp_neutron_pulse.f4_charge    = neutron_analysis.get_charge(u4n);
    tmp_neutron_pulse.f4_rise_time = neutron_analysis.get_rise(u4n);
    tmp_neutron_pulse.f4_fall_time = neutron_analysis.get_fall(u4n);
    tmp_neutron_pulse.f4_x = 100.;
    tmp_neutron_pulse.f4_y = 100.;
    tmp_neutron_pulse.f4_z = 100.;
    ew.pulses.push_back(tmp_neutron_pulse);
  }
  //if (u4_n_neutrons>0) get_message(bx_message::info)<<"Aligned V1731 neutron(s);"<<dispatch;
}


/*
 * $Log: bx_v1731sys.cc,v $
 * Revision 1.19  2009/11/20 10:55:24  davini
 * clusterize v1731ev
 *
 * Revision 1.18  2009-10-26 11:19:37  ddangelo
 * trgtype 'muon_tct' (or 'muon_totc') renamed as 'neutron' after real use
 * trgtype 'muon_mtb' renamed as 'muon' as duplicity with the above one disappeared and 'mtb' is associated to btb_inputs flag
 *
 * Revision 1.17  2008-12-17 15:31:22  davini
 * daq repaired and online since run 9313
 *
 * Revision 1.16  2008-12-15 14:50:34  davini
 * get skip_after_run8504 to disable data processing for some runs; junk data are due to hardware problems in sept08
 *
 * Revision 1.14  2008-12-15 13:37:16  davini
 * event index
 *
 * Revision 1.13  2008-11-27 14:18:42  davini
 * addes vars and lines for sequential readout; gzFile and pipe are private vars;
 *
 * Revision 1.12  2008-10-20 13:19:58  davini
 * repository path changed to bxmaster-data;
 *
 * Revision 1.11  2008-09-16 10:12:33  davini
 * removed debug messages
 *
 * Revision 1.10  2008-08-26 15:50:12  davini
 * added daq enabled and event associated flag; removed i4_n_pulses;
 *
 * Revision 1.9  2008-08-10 19:49:21  davini
 * trgid alignment improved
 *
 * Revision 1.8  2008-08-06 14:59:31  davini
 * added TRGID alignement; added fall times;
 *
 * Revision 1.7  2008-07-30 10:49:15  davini
 * debug
 *
 * Revision 1.6  2008-07-30 08:54:59  davini
 * code ANSI C++ standardized
 *
 *
 */

