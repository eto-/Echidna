/* BOREXINO Reconstruction program
 *
 * Author: Stefano Davini <stefano.davini@ge.infn.it>
 * Maintainer: Stefano Davini <stefano.davini@ge.infn.it>
 *
 * $Id:
 *
 * Definitons of classes and functions used in bx_v1731sys 
 * to correlate v1731 data to laben events  
 *
 *
 * 
 */

#ifndef _EVENT_CORRELATOR_HH
#define _EVENT_CORRELATOR_HH

#include <iostream>
#include <cstdio>
#include <string>
#include <cctype>
#include <ctime>
#include <cmath>
#include <vector>
#include <stdlib.h>
#include <algorithm>
#include "zlib.h"
#include "bx_named.hh"
#include "messenger.hh"

namespace v1731_event_correlator{
  class mapper_time_file;
  class v1731event_indexer;
  class event_finder;

  typedef std::vector<std::string> str_vec;
  enum e_transfer {http, file};

  struct transfer{
    e_transfer transfer_mode;
    std::string s_filepath;
    std::string s_wget_cmd;
  };

  void make_filename_vector(str_vec& s_filename_vector, const transfer&, bool& b_error_occurred);
  void make_filename_vectors(str_vec& s_filename_vector, str_vec& s_auxfilename_vector, const transfer&, bool& b_error_occurred);

};



/************** mapper_file_time*********** */

class v1731_event_correlator::mapper_time_file:public bx_named{
public:
  mapper_time_file (const str_vec& s_in_file_names, const str_vec& s_in_auxfile_names);
  virtual ~mapper_time_file() {}
  
  unsigned get_number_of_files() const {return struct_map_ft.size();}
  unsigned get_key_if_found(unsigned long u4_gpstime_sec) const;

  const std::string& get_Nfile_name(unsigned u4_key) const {return struct_map_ft.at(u4_key).s_Nfile_name;}
  const std::string& get_Tfile_name(unsigned u4_key) const {return struct_map_ft.at(u4_key).s_Tfile_name;}
  unsigned long get_time_start(unsigned u4_key) const {return struct_map_ft.at(u4_key).u4_time_start;}
  unsigned long get_time_stop (unsigned u4_key) const {return struct_map_ft.at(u4_key).u4_time_stop;}
  bool get_has_auxfile(unsigned u4_key) const {return struct_map_ft.at(u4_key).b_has_auxfile;}
  void show(unsigned u4_key);

private:
  struct struct_file_time{
    std::string s_Nfile_name;
    std::string s_Tfile_name;
    unsigned long u4_time_start;
    unsigned long u4_time_stop;
    bool b_has_auxfile;
  };
  std::vector<struct_file_time> struct_map_ft;

  unsigned long u4_date_to_sec(const std::string& s_in_date);
  tm tm_ref_2000;  // time ref (Jan 01, 2000, 00:00:00)
  time_t tt_ref_2000;
};


/**************** v1731event_indexer ********************/
class v1731_event_correlator::v1731event_indexer: public bx_named{
  struct event_info{
    std::string  s_filename;
    long  i4_position;   // row of T_*.dat file (and event position in N_*.dat.gz file)
  };
  struct file_borders{
    unsigned long u4_start;
    unsigned long u4_stop;
  };

public:
  v1731event_indexer(): bx_named("v1731event_indexer") {index_borders.reserve(3);}
  virtual ~v1731event_indexer() {}

  void set_transfer_param(const transfer& tr_par) {transfer_param = tr_par;}
  const transfer& get_transfer_param() const {return transfer_param;} 

  void new_event(const mapper_time_file& map_tf, unsigned long u4_gpstime);
  bool daq_active(unsigned long u4_gpstime) const;
  void get_event_info(unsigned long u4_gpstime, std::vector<std::string>& s_filenames, std::vector<long>& i4_poss) const;
  void show(unsigned long u4_gpstime);

private:
  transfer transfer_param;
  std::multimap<unsigned long, event_info> event_index;
  std::vector<file_borders> index_borders;
  void add_to_index(const mapper_time_file&, unsigned u4_key);
};


/********** event_finder *******************/

class v1731_event_correlator::event_finder: public bx_named{
 public:
  event_finder(const v1731event_indexer&, const unsigned long u4_sec);
  virtual ~event_finder() {}
  bool     get_found()    const {return (events_found.size()>0);}
  unsigned get_n_founds() const {return events_found.size();}

  std::string get_filepath() const {return transfer_param.s_filepath;}
  std::string get_wget_cmd() const {return transfer_param.s_wget_cmd;}
  e_transfer get_load_mode() const {return transfer_param.transfer_mode;}

  std::string get_file_name(unsigned u4) const {return ((u4<events_found.size()) ? events_found[u4].s_file_name: "");}  
  unsigned long get_v1731_timer(unsigned u4) const {return ((u4<events_found.size()) ? events_found[u4].u4_v1731_timer : 0);}
  long get_position(unsigned u4) const {return ((u4<events_found.size()) ? events_found[u4].i4_position : 0);}     
  short get_difference(unsigned u4) const {return ((u4<events_found.size()) ? events_found[u4].i2_difference : 10);}  

  void m_show();

 private: 
  struct struct_event{
    std::string s_file_name;
    unsigned long u4_v1731_timer;
    long  i4_position;   // row of T_*.dat file
    short i2_difference; // time difference between searched and found event
  };
  std::vector <struct_event> events_found;
  std::vector <struct_event> event_found_old;
  transfer transfer_param;

  void fill_found_events(unsigned long u4_gpstime, int i4_difftime, const v1731event_indexer&);
};

    
#endif

/*
 * $Log: v1731_event_correlator.hh,v $
 * Revision 1.5  2008/12/15 13:37:16  davini
 * event index
 *
 * Revision 1.4  2008-08-18 15:57:35  davini
 * removed warnings (thanks to A.R.)
 *
 * Revision 1.3  2008-07-30 08:55:12  davini
 * code ANSI C++ standardized
 *
 *
*/
