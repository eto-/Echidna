/* BOREXINO Reconstruction program
 * 
 * Author: Marcin Misiaszek
 * Maintainer: Marcin Misiaszek
 *
 * $Id: bx_barn.hh,v 1.25 2015/01/09 15:03:07 misiaszek Exp $
 *
 * A repository for db objects, to be used from root file. 
 * Store db_calib, db_profile, db_run and db_channel.
 * 
 * 
 *
 */

#ifndef _BX_BARN_H
#define _BX_BARN_H
#ifndef CYCLE_NUMBER
#define CYCLE_NUMBER 18
#endif

#include <TNamed.h>
#include <vector>
#include <string>
#include <sstream>
#include <map>
#include <stdexcept>
class db_run;
class db_profile;
class db_calib;
class db_channel;
class bx_detector;

namespace std {
  template <> class pair<int, db_channel*> {
    pair(const pair &) {};
  };
}

class bx_barn: public TNamed {
   public:
    bx_barn (): TNamed () {}
    ~bx_barn ();

      // List of run number present in bx_barn    
    const std::vector<int>& get_run_list () const { return v_run_list; }

      // DB visitors
    const db_run *get_db_run (int i_run) const { return map_get (i_run, m_run_runs, "runs"); }
    const db_profile *get_db_profile (int i_run) const { return map_get (i_run, m_run_profiles, "profiles"); } 
    const db_calib *get_db_calib (int i_run) const { return map_get (i_run, m_run_calibs, "calibs"); }
    typedef std::map<int, db_channel*> map_channels;
    const map_channels& get_db_channels_map (int i_run) const { return *(map_get (i_run, m_run_channels, "channels")); };
    const db_channel *get_db_channel (int i_run, int i_lg) const { return map_get (i_lg, get_db_channels_map (i_run), "channels"); }

      // BxDetector
    const bx_detector* get_bx_detector (int i_run) const { return map_get (i_run, m_run_detector, "detectors"); }

    void merge (const bx_barn*);
    int live_time () const;
    int astro_time () const;
   private:   
    std::vector<int> v_run_list;
    std::map<int, map_channels*> m_run_channels;
    std::map<int, db_profile*> m_run_profiles;
    std::map<int, db_run*> m_run_runs;
    std::map<int, db_calib*> m_run_calibs;
    std::map<int, bx_detector*> m_run_detector;

    bx_barn (const char* name, const char* title);

    template <typename key_t, typename value_t> const value_t& map_get (key_t k, const std::map<key_t, value_t>& map, const std::string& map_name) const { 
      typename std::map<key_t, value_t>::const_iterator item = map.find (k);
      if (item == map.end ()) {
	std::ostringstream msg;
	msg << "map: element " << k << " not found";
	if (map_name.size ()) msg << " in map \"" << map_name << "\"";
#ifdef CMAP_ABORT
	abort ();
#else
	throw std::runtime_error(msg.str ());
#endif	
      }
      return item->second;
    }

  friend class barn_interface;
  ClassDef(bx_barn,CYCLE_NUMBER)
};

#endif

