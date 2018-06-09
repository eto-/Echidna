/* BOREXINO Reconstruction program
 * 
 * Author: Alessandro Razeto <Alessandro.Razeto@ge.infn.it>
 * Maintainer: Alessandro Razeto <Alessandro.Razeto@ge.infn.it>
 *
 * $Id: bx_dbi.hh,v 1.27 2014/11/05 13:41:28 marcocci Exp $
 *
 * The database interface
 *
 */
#ifndef _BX_DBI_H
#define _BX_DBI_H
#include "vdt.hh"
#include "cmap.hh"
#include "messenger.hh"
#include "bx_named.hh"
#include "db_run.hh"
#include "db_profile.hh"
#include "db_calib.hh"
#include "db_channel.hh"

#include <string>
#include <map>
#include <vector>

class bx_reco_framework;

class bx_dbi: public bx_named {
  public:
    enum database_id {
      daq_config = 0,
      bx_geometry,
      bx_calib,
      bx_precalib,
      bx_physics,
      bx_slow,
      channelhistory,
    };
    typedef std::vector<vdt> column;
    typedef std::vector<std::string> column_names;
    typedef std::cmap<std::string, column> table;

    static bx_dbi* get () { if (!me) me = new bx_dbi; return me; } 

    // Do a query: the first method allow not to type a valid SQL query, while the second is present
    // if the user has to make some complicated query, not foreseen by the first.
    // max_lines is the maximum number of line to read and has the following syntax:
    //   0 means unlimited (allow 0 touples)
    //   > 0 do not read more the indicated lines (soft limit)
    //   < 0 if the query returns more than the indicated lines (abs value) there is an exception (hard limit)
    // Every field has to be quoted with the SQL syntax, ie table_fields has to be "\"ProfileID\""
    // For best use a const reference has to be assigned to the return value of query, as follow:
    // const bx_dbi::table &table = dbi->query (...);
    const table query (database_id db_id, const std::string& from_table, const std::string& where,
		  const std::string &table_fields = "*", const std::string& order = "", int max_lines = 0);
    const table query (database_id db_id, const std::string& sql_query, int max_lines);
    
    // Insert a table in the database
    void insert (database_id db_id, const std::string& to_table, const table& t, const column_names& index_columns, bool overwrite = 0);


    // Real write data
    void flush_visitors (bx_reco_framework *);
    void close_db_connections (); // close all db connections
    
    int get_current_run_number () const { return i4_current_run_number; }
    void set_current_run_number (int run_number) { i4_current_run_number = run_number; }

    // Get db_run/db_profile table for the requested run/profile.
    db_run& get_run () { if (!run_p) run_p = new db_run(i4_current_run_number); return *run_p; }
    db_profile& get_profile () { if (!profile_p) profile_p = new db_profile(get_run ().get_profile_id ()); return *profile_p; }
    db_calib& get_calib () { if (!calib_p) calib_p = new db_calib(get_run ().get_calib_profile ()); return *calib_p; }                   
    const db_channel& get_channel (int lg) { if (!channel_p[lg]) channel_p[lg] = db_channel_builder::build (lg, get_run (), get_profile ()); return *channel_p[lg]; }                   

    bool is_production () const { return b_production; }
  private:
    bx_dbi ();
    bx_dbi (const bx_dbi& dbi): bx_named ("illegal") { } // I do not want a copy constructor
    static bx_dbi *me;

    bool b_production;
    long int i4_current_run_number;

    void m_open_connection (database_id id);
    void m_close_connection (database_id id);
    bool m_check_priv (database_id db_id, std::string to_table, const std::string& type);
    bool m_transaction (database_id db_id, bool open);
    bool m_check_data_presence (database_id db_id, const std::string& to_table, const table& t, const column_names& index_columns);

      // The database are accessed by index of database_id type; since, unfortunatelly,
      // the database name is mutable a double reference process has to be done using
      // db_name_hash_index which indexes formal DB names using standards names.
      // Then the params map contain a lookup from these standard names and the 
      // db name. It is important that the paramenter names in the config file
      // (the name before the colon) are excactly the formal names.
      // BTW, no one should change the format of the configuration file.
    typedef std::map<database_id, std::string> db_name_hash_index;
    typedef std::map<database_id, void *> db_connection;
    db_name_hash_index hash_name_indexes;
    db_connection connections;
    
      // since db_run has private constructors and is friend to bx_dbi a map
      // of db_run objects can not be done (since the map allocator fails in
      // initializing the stored object). An array of pointers is instead possible.
      // Idem for profile map.
    db_run* run_p;
    db_profile* profile_p;
    db_calib* calib_p;
    db_channel** channel_p;

    bx_message message;
};

#endif
/*
 * $Log#
 */
