/* BOREXINO Reconstruction program
 *
 * Author: Alessandro Razeto <Alessandro.Razeto@ge.infn.it>
 * Maintainer: Alessandro Razeto <Alessandro.Razeto@ge.infn.it>
 *
 * $Id: bx_dbi.cc,v 1.36 2014/11/05 13:41:28 marcocci Exp $
 *
 * Implementation of bx_dbi
 *
 */
#include "bx_dbi.hh"
#include "db_run.hh"
#include "db_profile.hh"
#include "db_calib.hh"
#include "db_channel.hh"

#include <iostream>
#include <fstream>
#include <sstream>
#include <libpq-fe.h>
#include <pg_config.h>
#include <sys/types.h>
#include <pwd.h>
#include <string.h>

bx_dbi *bx_dbi::me = 0;

bx_dbi::bx_dbi (): bx_named("bx_dbi"), run_p(0), profile_p(0), calib_p(0) {
  channel_p=new db_channel*[4000]();
    // Set the enumerator/name corrispondence
  hash_name_indexes[daq_config] = "DAQCONFIGDBNAME";
  hash_name_indexes[bx_geometry] = "BXGEOMETRYDBNAME";
  hash_name_indexes[bx_calib] = "BXCALIBDBNAME";
  hash_name_indexes[bx_precalib] = "BXPRECALIBDBNAME";
  hash_name_indexes[bx_physics] = "BXPHYSICSDBNAME";
  hash_name_indexes[bx_slow] = "BXSLOWDBNAME";
  hash_name_indexes[channelhistory] = "CHHISTORYDBNAME";
  
  endpwent ();
  b_production = false;
  while (struct passwd *p = getpwent ()) {
    if (p->pw_uid == getuid () && std::string (p->pw_name) == "production") {
      b_production = true;
      break;
    }
  }
  endpwent ();

//  connections[daq_config] = connections[bx_geometry] = connections[bx_calib] = connections[bx_precalib] = connections[bx_physics] = connections[bx_slow] = connections[channelhistory] = 0;
}

void bx_dbi::close_db_connections () {
  m_close_connection (daq_config);
  m_close_connection (bx_geometry);
  m_close_connection (bx_calib);
  m_close_connection (bx_precalib);
  m_close_connection (bx_physics);
  m_close_connection (bx_slow);
  m_close_connection (channelhistory);
}


void bx_dbi::m_close_connection (database_id id) {
  if (!connections[id]) return;
  PQfinish ((PGconn*)connections[id]);
  connections[id] = 0;
  get_message (bx_message::debug) << "closed db connection for " << get_parameter (hash_name_indexes[id]) << dispatch;
}

void bx_dbi::m_open_connection (database_id id) {

  const char *db = get_parameter (hash_name_indexes[id]).get_string ().c_str ();
  const char *user = get_parameter ("PGUSER").get_string ().c_str (); 
  const char *passwd = "";
  if (check_parameter ("PGPASSWD")) passwd = get_parameter("PGPASSWD").get_string ().c_str ();

  const char *ip = get_parameter ("SERVER1").get_string ().c_str ();
  int port_i = get_parameter ("PORT1").get_int ();
  const char *ip2 = get_parameter ("SERVER2").get_string ().c_str ();
  int port_i2 = get_parameter ("PORT2").get_int ();

  std::ostringstream s_port;
  s_port << port_i;
  
  get_message (bx_message::debug) << "PQsetdbLogin (" << ip << "," << port_i << ", 0, 0, " << db << ", " << user << "," << passwd << ")" << dispatch;
  PGconn *connection = PQsetdbLogin (ip, s_port.str ().c_str (), 0, 0, db, user, passwd);


  if (PQstatus (connection) == CONNECTION_BAD) {
    get_message (bx_message::warn) << "connection with " << ip << " failed: " << PQerrorMessage (connection) << dispatch;

    std::ostringstream s_port2;
    s_port2 << port_i2;
    get_message (bx_message::debug) << "PQsetdbLogin (" << ip2 << "," << port_i2 << ", 0, 0, " << db << ", " << user << "," << passwd << ")" << dispatch;
    connection = PQsetdbLogin (ip2, s_port2.str ().c_str (), 0, 0, db, user, passwd);
  }
  
  if (PQstatus (connection) == CONNECTION_BAD) {
    get_message (bx_message::warn) << "connection with " << ip2 << " failed: " << PQerrorMessage (connection) << dispatch;
    get_message (bx_message::critic) << "cannot connect to both databases " <<  get_parameter ("SERVER1") << " and " << get_parameter ("SERVER2") << dispatch;
  }

  connections[id] = (void *)connection;
  get_message (bx_message::debug) << "opened db connection for " << get_parameter (hash_name_indexes[id]) << dispatch;
}


const bx_dbi::table bx_dbi::query (database_id db_id, const std::string& from_table, const std::string& where,
		      const std::string &table_fields, const std::string& order, int max_lines) {
  std::ostringstream sql_query;
  
  sql_query << "SELECT " << table_fields << " FROM " << from_table;
  if (where.size ()) sql_query << " WHERE " << where;
  if (order.size ()) sql_query << " ORDER BY " << order;
  
  return query (db_id, sql_query.str (), max_lines);
}

const bx_dbi::table bx_dbi::query (database_id db_id, const std::string& sql_query, int max_lines) {
  const std::string &db_name = get_parameter (hash_name_indexes[db_id]).get_string ();

  table query_result("query_result");
  query_result.clear ();

  if (!connections[db_id]) m_open_connection (db_id);
  
  PGresult *result = PQexec ((PGconn *)connections[db_id], sql_query.c_str ());

  if (PQresultStatus (result) != PGRES_TUPLES_OK) 
    get_message (bx_message::critic) << "query \"" << sql_query << "\" on database \"" << db_name << "\" failed with error \"" << PQresultErrorMessage (result) << "\"" << dispatch;

  int tuples = PQntuples (result);
  if (!tuples) {
    if (max_lines != 0) get_message (bx_message::warn) << "no tuples found for query \"" << sql_query << "\" on database \"" << db_name << "\"" << dispatch;
    return query_result;
  }

  if (max_lines < 0) {
    max_lines *= -1;
    if (tuples < max_lines) 
      get_message (bx_message::critic) << "query \"" <<  sql_query << "\" on database \"" << db_name + "\" returned " << tuples << " tuples while only " << max_lines << " were expected" << dispatch;
  }


  int lines = max_lines < tuples ? max_lines : tuples;
  if (!max_lines) lines = tuples;

  int fields = PQnfields (result);
  for (int i = 0; i < fields; i++) {
    const char *field = PQfname (result, i);
    query_result[field].reserve (lines);
    for (int j = 0; j < lines; j++) {
      const char *value = PQgetvalue (result, j, i);
      query_result[field].push_back (vdt(value, field));
    }
  }

  PQclear (result);
  return query_result;
}

bool bx_dbi::m_transaction (database_id db_id, bool open) {
  std::string cmd("BEGIN;");
  if (!open) cmd = "COMMIT;";
  
  PGresult *result = PQexec ((PGconn *)connections[db_id], cmd.c_str ());
  if (PQresultStatus (result) != PGRES_COMMAND_OK) {
    const std::string &db_name = get_parameter (hash_name_indexes[db_id]).get_string ();
    if (open) get_message (bx_message::critic) << "starting a tranaction block on on database \"" << db_name << "\" failed with error \"" << PQresultErrorMessage (result) << "\"" << dispatch;
    else get_message (bx_message::critic) << "committing a tranaction block on on database \"" << db_name << "\" failed with error \"" << PQresultErrorMessage (result) << "\"" << dispatch;
    return false;
  }

  PQclear (result);
  return true;
}

bool bx_dbi::m_check_priv (database_id db_id, std::string to_table, const std::string& type) {
    // Strip \" on table name since 7.2.2 does not support it
  if (std::string ("7.2.2") == std::string (PG_VERSION)) {
    to_table.erase (to_table.begin ());
    to_table.erase (to_table.end () - 1);
  }
  
  std::string check_privilege_query = std::string("SELECT has_table_privilege(\'") + to_table + "\', '" + type + "')";
  if (query (db_id, check_privilege_query, -1)["has_table_privilege"][0].get_string () != std::string ("t")) {
    const std::string &db_name = get_parameter (hash_name_indexes[db_id]).get_string ();
    get_message (bx_message::error) << "you do not have permission to write on table " << to_table << " of database " << db_name << dispatch;
    return false;
  }

  return true;
}

bool bx_dbi::m_check_data_presence (database_id db_id, const std::string& to_table, const table& t, const column_names& index_columns) {
  if (!index_columns.size ()) return false;
  
  std::ostringstream where;
  where << "\"" << index_columns[0] << "\"=" << t[index_columns[0]][0];
  for (unsigned int i = 1; i < index_columns.size (); i++) where << " and \"" << index_columns[i] << "\"=" << t[index_columns[i]][0];

  const table &j = query (db_id, to_table, where.str ());
  if (j.size ()) return true;
  
  return false;
}
  
void bx_dbi::insert (database_id db_id, const std::string& to_table, const table& t, const column_names& index_columns, bool overwrite) {
  if ((!is_production () && !get_parameter ("force_write").get_bool ()) || !get_parameter ("write_enabled").get_bool ()) {
    get_message (bx_message::error) << "db in readonly mode (set force_write or become production)" << dispatch;
    return;
  }

  if (!t.begin ()->second.size ()) {
    get_message (bx_message::error) << "no data to insert in " << to_table << dispatch;
    return;
  }

  const std::string &db_name = get_parameter (hash_name_indexes[db_id]).get_string ();

  if (connections[db_id]) m_open_connection (db_id);

    // Check write privileges
  //if (!m_check_priv (db_id, to_table, "INSERT")) return;
  
  ///bool do_overrite = false;
  if (m_check_data_presence (db_id, to_table, t, index_columns)) {
//    if (overwrite) do_overrite = true;
//    else {
      get_message (bx_message::error) << "data already present in table \"" << to_table << "\"" << dispatch;
//      return;
//    }
  }

    // Open the transaction
  if (!m_transaction (db_id, true)) return;
  
    // Prepare the query common string
  std::ostringstream str;
  table::const_iterator it = t.begin ();
  str << "INSERT INTO " << to_table << " (\"" << it->first;
  for (it++; it != t.end (); it++) str << "\", \"" << it->first;
  str << "\") VALUES (";
  
    // Send queries  
  for (unsigned int i = 0; i < t.begin ()->second.size (); i++) {
    std::ostringstream full_query;
    it = t.begin ();
  
    if (it->second[i].get_type () == vdt::int_vdt || it->second[i].get_type () == vdt::float_vdt) full_query << str.str () << it->second[i];
    else full_query << str.str () << "'" << it->second[i] << "'";
    for (it++; it != t.end (); it++) {
      if (it->second[i].get_type () == vdt::int_vdt || it->second[i].get_type () == vdt::float_vdt) full_query << ", " << it->second[i];
      else full_query << ", '" << it->second[i] << "'";
    }
    full_query << ")";
  
    PGresult *result = PQexec ((PGconn *)connections[db_id], full_query.str ().c_str ());
    if (PQresultStatus (result) != PGRES_COMMAND_OK) 
      get_message (bx_message::critic) << "query \"" << full_query.str () << "\" on database \"" << db_name << "\" failed with error \"" << PQresultErrorMessage (result) << "\"" << dispatch;
    PQclear (result);
  }
  
    // Commit the transaction
  if (!m_transaction (db_id, false)) return;
}

void bx_dbi::flush_visitors (bx_reco_framework *f) {
  run_p->flush ();
//  for (db_profile_map::iterator it = profile_map.begin (); it != profile_map.end (); it++) *it->flush ();
//  for (db_calib_map::iterator it = calib_map.begin (); it != calib_map.end (); it++) *it->flush ();
}
/*
 * $Log: bx_dbi.cc,v $
 * Revision 1.36  2014/11/05 13:41:28  marcocci
 * fixed a bad memory allocation for g4bx2 compatibility
 *
 * Revision 1.35  2012/04/02 11:40:09  razeto
 * Added write_enabled parameter + cngs configuration
 *
 * Revision 1.34  2011-03-02 17:21:33  razeto
 * Included added for memset
 *
 * Revision 1.33  2011-03-01 19:21:19  razeto
 * Faster method for get_channel
 *
 * Revision 1.32  2008-10-01 16:17:52  razeto
 * Error if writing nothing into the database
 *
 * Revision 1.31  2008-09-26 14:02:41  razeto
 * DB write authorized only to production user
 *
 * Revision 1.30  2007-06-03 16:13:30  razeto
 * Do not keep db connection opened (else bxdb has lots of problems)
 *
 * Revision 1.29  2007-03-29 09:01:51  razeto
 * Clear query result (fix small leak)
 *
 * Revision 1.28  2006/08/30 09:43:44  ludhova
 * removed useless and buggy check
 *
 * Revision 1.27  2006/07/18 14:47:41  razeto
 * Fixed quoting for composite fields
 *
 * Revision 1.26  2006/07/18 14:35:43  razeto
 * Added data type check before insert
 *
 * Revision 1.25  2006/05/08 17:31:33  razeto
 * Added db_channel patch for fadc (sent to the mailing list)
 *
 * Revision 1.24  2005/08/04 15:00:44  razeto
 * Added a parameter to dbi
 *
 * Revision 1.23  2005/06/29 11:48:34  razeto
 * Removed libpq-fe.h include to allow rootcint on bx_dbi; this was done using
 * void * intead of PGconn * for connection pointers (with some casting).
 *
 * Revision 1.22  2005/05/05 08:53:42  razeto
 * Now bx_dbi is a bx_named and gets params from parameter_broker (ie echidna|user.cfg)
 *
 * Revision 1.21  2005/04/04 13:23:28  razeto
 * Added a workaround for a 7.2.2 postgresql anomaly.
 *
 * Revision 1.20  2004/12/07 16:01:31  razeto
 * Fixed a memory leak bug
 *
 * Revision 1.19  2004/11/26 15:25:10  razeto
 * Thanks to Davide who so kindly underlined a syntactical error, I fixed the "Maintainer" word
 *
 * Revision 1.18  2004/11/26 14:01:20  razeto
 * Added Mantainer field
 *
 * Revision 1.17  2004/11/24 13:03:58  razeto
 * Upgraded to the new cmap ctor with name
 *
 * Revision 1.16  2004/11/16 16:01:38  razeto
 * Upgraded database writing with transaction. Overrite support almost done
 *
 * Revision 1.15  2004/10/19 18:41:58  razeto
 * Integrated visitors writing in the framework
 *
 * Revision 1.14  2004/10/19 16:21:12  razeto
 * Added bx_precalib database stuff
 * Upgraded some params with named vdt
 * Upgraded query method to not complain on empty tables
 * Added insert method to write the database. First attempt, a lot to be done.
 *
 * Revision 1.13  2004/09/16 11:45:21  razeto
 * Added db_calib
 *
 * Revision 1.12  2004/09/14 18:24:34  razeto
 * Fixed a bug (DDangelo)
 *
 * Revision 1.11  2004/05/26 11:02:44  razeto
 * Changed the query return type from reference to const value, which even if
 * requires more copies, allows to have different result tables at the same time.
 *
 * Revision 1.10  2004/05/26 10:31:06  razeto
 * removed unusefull stuff
 *
 * Revision 1.9  2004/05/26 09:21:29  razeto
 * Upgraded to set column name to vdt
 *
 * Revision 1.8  2004/05/19 11:07:05  razeto
 * Updated to manage an internal hash of db_channel, in a way similar to db_profile
 *
 * Revision 1.7  2004/05/18 14:25:10  razeto
 * Fixed a bug. Some updates
 *
 * Revision 1.6  2004/04/26 14:18:05  razeto
 * Bug still not found, but fixed using port_i
 *
 * Revision 1.5  2004/04/26 13:49:17  razeto
 * Removed db_run/db_profile header include in bx_dbi.hh
 *
 * Revision 1.4  2004/04/09 07:56:48  razeto
 * Changed the enum names for some dbs to a more familiar name: moved the enum in the class local namespace.
 * Fixed the sql query build syntax in query method.
 * Exported get_message.
 * Added db_run and db_profile getters.
 *
 * Revision 1.3  2004/04/05 13:14:45  razeto
 * Used messages for expcetions
 *
 * Revision 1.2  2004/04/03 09:25:46  razeto
 * Added messenger
 *
 * Revision 1.1  2004/04/01 13:02:04  razeto
 * Added bx_dbi
 *
 */
