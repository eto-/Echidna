/* BOREXINO Reconstruction program
 * 
 * Author: Marcin Misiaszek
 * Maintainer: Marcin Misiaszek
 *
 * $Id: bx_barn.cc,v 1.18 2011/03/01 19:21:19 razeto Exp $
 *
 * Implementation of db_barn
 *
 */


#include <algorithm>
#include "bx_barn.hh"
#include "bx_dbi.hh"
#include "constants.hh"
#include "db_run.hh"
#include "db_profile.hh"
#include "db_calib.hh"
#include "db_channel.hh"
#include "bx_detector.hh"

ClassImp(bx_barn);


#ifndef _ECHIDNA_ROOTLIB_
bx_barn::bx_barn (const char* name, const char* title) :  TNamed(name, title) {
  bx_dbi* dbi = bx_dbi::get();
  int i_run_number = dbi->get_run().get_number();

  v_run_list.push_back(i_run_number);

  m_run_profiles[i_run_number] = &dbi->get_profile ();
  m_run_runs[i_run_number] = &dbi->get_run ();
  m_run_calibs  [i_run_number] = &dbi->get_calib();
  map_channels *channels = new map_channels;
  m_run_channels[i_run_number] = channels;
  for(int lg = 1 ; lg <= constants::laben::channels; lg++) {
    const db_channel* channel = &dbi->get_channel(lg);
    (*channels)[lg] = (db_channel*)channel;
  };
  for(int lg = constants::muon::channel_offset + 1; lg <= constants::muon::channel_offset + constants::muon::channels; lg++) {
    const db_channel* channel = &dbi->get_channel(lg);
    (*channels)[lg] = (db_channel*)channel;
  };

  m_run_detector[i_run_number] = detector_interface::get ();
};
#endif

bx_barn::~bx_barn () {
  for (std::map<int, map_channels*>::iterator it = m_run_channels.begin (); it != m_run_channels.end (); it++) {
    for (map_channels::iterator jt = it->second->begin (); jt != it->second->end (); jt++) delete jt->second;
    delete it->second;
  }
  for (std::map<int, db_profile*>::iterator it = m_run_profiles.begin (); it != m_run_profiles.end (); it++) delete it->second;
  for (std::map<int, db_run*>::iterator it = m_run_runs.begin (); it != m_run_runs.end (); it++) delete it->second;
  for (std::map<int, db_calib*>::iterator it = m_run_calibs.begin (); it != m_run_calibs.end (); it++) delete it->second;
  for (std::map<int, bx_detector*>::iterator it = m_run_detector.begin (); it != m_run_detector.end (); it++) delete it->second;
}

void bx_barn::merge (const bx_barn* that) {
  v_run_list.insert (v_run_list.end (), that->v_run_list.begin (), that->v_run_list.end ());
  std::sort (v_run_list.begin (), v_run_list.end ());
  std::vector<int>::iterator new_end = std::unique (v_run_list.begin (), v_run_list.end ());
  v_run_list.erase (new_end, v_run_list.end ());

  for (std::map<int, map_channels*>::const_iterator it = that->m_run_channels.begin (); it != that->m_run_channels.end (); it++) 
    if (!m_run_channels[it->first]) {
      map_channels *mc = new map_channels;
      for (map_channels::const_iterator jt = it->second->begin (); jt != it->second->end (); jt++)
        if(typeid(*(jt->second)) == typeid(db_channel_laben) ) (*mc)[jt->first] = new db_channel_laben (static_cast<db_channel_laben&>(*(jt->second)));
        else (*mc)[jt->first] = new db_channel (*(jt->second));
      
      m_run_channels[it->first] = mc;
    }
  for (std::map<int, db_profile*>::const_iterator it = that->m_run_profiles.begin (); it != that->m_run_profiles.end (); it++) if (!m_run_profiles[it->first]) m_run_profiles[it->first] = new db_profile(*(it->second));
  for (std::map<int, db_run*>::const_iterator it = that->m_run_runs.begin (); it != that->m_run_runs.end (); it++) if (!m_run_runs[it->first]) m_run_runs[it->first] = new db_run(*(it->second));
  for (std::map<int, db_calib*>::const_iterator it = that->m_run_calibs.begin (); it != that->m_run_calibs.end (); it++) if (!m_run_calibs[it->first]) m_run_calibs[it->first] = new db_calib(*(it->second));
  for (std::map<int, bx_detector*>::const_iterator it = that->m_run_detector.begin (); it != that->m_run_detector.end (); it++) if (!m_run_detector[it->first]) m_run_detector[it->first] = new bx_detector(*(it->second));
}

int bx_barn::live_time () const {
  double sum = 0;
  for (std::map<int, db_run*>::const_iterator it = m_run_runs.begin (); it != m_run_runs.end (); it++) sum += difftime(it->second->get_stop_time (), it->second->get_start_time ());
  return int(sum);
}

int bx_barn::astro_time () const {
  if (!m_run_runs.size ()) return 0;
  double sum = ::difftime (m_run_runs.rbegin ()->second->get_stop_time (), m_run_runs.begin ()->second->get_start_time ());
  return int(sum);
}
