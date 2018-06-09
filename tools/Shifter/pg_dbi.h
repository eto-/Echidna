#ifndef _PG_DBI_H
#define _PG_DBI_H

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "/usr/include/postgresql/libpq-fe.h"
//#include "/usr/include/postgresql/pg_config.h"

// DATABASE description ***************

const char user[]             = "borex_guest";
const char passwd[]           = "xyz";
const char ip1[]              = "bxdb.lngs.infn.it";
const char ip2[]              = "localhost";
const char port1[]            = "5432";
const char port2[]            = "5432";

const int  NDB                = 8;
const char db_names[NDB][127] = {"daq_config", "bx_geometry", "bx_calib", "bx_precalib", "bx_physics", "bx_runvalidation", "bx_slow", "channelhistory"};

enum database_id {
  DAQ_CONFIG = 0,
  BX_GEOMETRY,
  BX_CALIB,
  BX_PRECALIB,
  BX_PHYSICS,
  BX_RUNVALIDATION,
  BX_SLOW,
  CHANNELHISTORY
};

// ************************************

class pg_dbi {

  public:
    pg_dbi () {
      if (!me) me = this;
      else printf("pg_dbi:: ATTENTION: more than one instance of the PG database manager created.\n");
      for (int i = 0; i < NDB; i++) connections[i] = NULL;
      q_result = NULL; q_tuples = q_fields = 0;
    }
    ~pg_dbi() { close_connections(); }

    static pg_dbi*   get ()                                        { if (!me) me = new pg_dbi; return me; }

    const PGconn*    open_connection(int id);
    int              open_connections();                           //{ for (int i = 0; i < NDB; i++) open_connection(i); }
    void             close_connection(int id);
    void             close_connections()                           { for (int i = 0; i < NDB; i++) close_connection(i); }
    int              is_opened(int id)                             { if (connections[id]) return 1; else return 0; }

    const PGresult*  query_result() const                          { return q_result; }
    static void      test_query(int db_id = -1, const char* sql_query = NULL);

    // Do a query: the first method allow not to type a valid SQL query, while the second is present
    // if the user has to make some complicated query, not foreseen by the first.
    // max_lines is the maximum number of line to read and has the following syntax:
    //   <= 0 means unlimited (allow 0 touples)
    //   >  0 do not read more the indicated lines
    // Every field has to be quoted with the SQL syntax, ie table_fields has to be "\"ProfileID\""
    // Number of successfully read lines is returned (or -1 if the query is failed)

    int              query (int db_id, const char* table_fields,  const char* from_table,
                            const char* where = NULL, const char* order = NULL, const char *max_lines = NULL);
    int              query (int db_id, const char* sql_query, int wr = 0);

    const char*      query_value(const int line, const int column);

    // !! query_table(...): you >> MUST << guarantee that the allocation table (lines x columns)
    // is large enough to store the query table:
    //
    // pg_dbi::get()->query_dimensions(lines, columns);
    // char ***table = new char**[lines];
    // for (i = 0; i < lines; i++){
    //   table[i] = new char*[columns];
    //   for (j = 0; j < columns; j++) table[i][j] = new char[2048];
    // }
    //

    int              query_table (char ***table, const int lines, const int columns, int max_lines);

    void             query_dimensions(int &lines, int &columns)     { lines = q_tuples; columns = q_fields; }
    const char*      query_column_name(const int column);
    int              query_value_length (const int line, const int column);

    // not tested

    int              check_privilege(int db_id, const char* to_table, const char* type);
    int              transaction(int db_id, int open);
    //int              check_data_presence (int db_id, const char* to_table, const char* t, const char* index_columns);

    // Insert a table in the database
    //void insert (int db_id, const std::string& to_table, const table& t, const column_names& index_columns, int overwrite = 0);

  private:
    static pg_dbi *me;

    void     *connections[NDB];
    PGresult *q_result;
    int  q_tuples, q_fields;
};

#endif
