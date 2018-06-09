// PG database interface
// (modified standalone version of bx_dbi class written by A. Razzeto)
// author and maintainer: K. Fomenko
#include "pg_dbi.h"

pg_dbi *pg_dbi::me = 0;

void pg_dbi::test_query(int id, const char *the_query){
  if (!me) return;

  int lines, columns;
  if (id >= 0 && id < NDB && the_query){
    if (me->query(id, the_query) < 0) return;
  }
  me->query_dimensions(lines, columns);
  for (int ii = 0; ii < lines; ii++) {
    for (int jj = 0; jj < columns; jj++) {
      printf("%s ", me->query_value(ii, jj));
    }
    printf("\n");
  }
}

const PGconn* pg_dbi::open_connection (int id) {

  if (connections[id]) return (PGconn *) connections[id];

  char log_str[2048];
  PGconn *connection;

  sprintf(log_str, "host = %s port = %s dbname = %s user = %s password = %s", ip1, port1, db_names[id], user, passwd);
  connection = PQconnectdb (log_str);

  if (PQstatus (connection) == CONNECTION_BAD) {
    fprintf(stderr, "pg_dbi::open_connection: Connection to %s failed: %s\n", ip1, PQerrorMessage (connection));

    sprintf(log_str, "host = %s port = %s dbname = %s user = %s password = %s", ip2, port2, db_names[id], user, passwd);
    connection = PQconnectdb (log_str);
  }

  if (PQstatus (connection) == CONNECTION_BAD) {
    fprintf(stderr, "pg_dbi::open_connection: Connection to %s failed: %s\n", ip2, PQerrorMessage (connection));
    return NULL;
  }

  connections[id] = (void *) connection;

  return connection;
}

int pg_dbi::open_connections() {
  int n_connections = 0;

  for (int i = 0; i < NDB; i++)
    if (open_connection(i)) ++n_connections;

  return n_connections;
}


void pg_dbi::close_connection(int id){
  if (connections[id]){
    PQfinish((PGconn *) connections[id]);
    connections[id] = NULL;
  }
}


int pg_dbi::query (int db_id, const char* table_q_fields, const char* from_table,
			       const char* where, const char* order, const char *max_lines)
{
  char sql_query[2048], limit[128];

  if (max_lines > 0) strcpy(limit, max_lines);
  else strcpy(limit, "ALL");

  if (order)
    sprintf(sql_query, "SELECT \"%s\" FROM \"%s\" WHERE %s ORDER BY \"%s\" LIMIT %s", table_q_fields, from_table, where, order, limit);
  else if (where)
    sprintf(sql_query, "SELECT \"%s\" FROM \"%s\" WHERE %s LIMIT %s", table_q_fields, from_table, where, limit);
  else
    sprintf(sql_query, "SELECT \"%s\" FROM \"%s\" LIMIT %s", table_q_fields, from_table, limit);

  //printf("DEBUG: query was: %s\n", sql_query);

  return query(db_id, sql_query);
}


int pg_dbi::query (int db_id, const char* sql_query, int wr) {

  if (!connections[db_id]){
    fprintf(stderr, "pg_dbi::query: database \"%s\" is not opened, ignoring request\n", db_names[db_id]);
    return -1;
  }

  q_result = PQexec((PGconn*) connections[db_id], sql_query);

  if (PQresultStatus(q_result) != PGRES_TUPLES_OK && !wr){
    fprintf(stderr, "pg_dbi::query: \"%s\" on database \"%s\" failed: %s", sql_query, db_names[db_id], PQresultErrorMessage(q_result));
    q_tuples = q_fields = 0;
    q_result = NULL;
    return -1;
  }

  q_tuples = PQntuples(q_result);
  q_fields = PQnfields(q_result);

  return q_tuples;
}


const char* pg_dbi::query_value (const int line, const int column) {
  if (!q_result){
    fprintf(stderr, "pg_dbi::query_value: resulting query table does not exist\n");
    return NULL;
  }

  if (column >= q_fields || line >= q_tuples){
    fprintf(stderr, "pg_dbi::query_value: requested value is out of resulting query table dimensions\n");
    return NULL;
  }

  return PQgetvalue(q_result, line, column);
}


int pg_dbi::query_value_length (const int line, const int column) {
    if (!q_result){
	fprintf(stderr, "pg_dbi::query_value: resulting query table does not exist\n");
	return 0;
    }

    if (column >= q_fields || line >= q_tuples){
	fprintf(stderr, "pg_dbi::query_value: requested value is out of resulting query table dimensions\n");
	return 0;
    }

    return PQgetlength(q_result, line, column);
}


int pg_dbi::query_table (char ***table, const int lines, const int columns, int max_lines) {
  if (!q_result){
    fprintf(stderr, "pg_dbi::query_table: resulting query table does not exist\n");
    return -1;
  }

  if (lines > q_tuples || columns > q_fields){
    fprintf(stderr, "pg_dbi::query_table: dimensions mismatch: query_tbl: %d x %d <-> allocation_tbl: %d x %d\n", q_tuples, q_fields, lines, columns);
    return -1;
  }

  if (max_lines > 0) max_lines = max_lines > q_tuples ? q_tuples : max_lines;
  else max_lines = q_tuples;

  int i, j;
  //table = new char**[lines];
  for (i = 0; i < max_lines; i++){
    //table[i] = new char*[columns];
    for (j = 0; j < columns; j++){
      //printf("DEBUG: line %d column %d value %s\n", i, j, PQgetvalue(q_result, i, j));
      //table[i][j] = new char[256];
      strcpy(table[i][j], PQgetvalue(q_result, i, j));
    }
  }

  return max_lines;
}


const char* pg_dbi::query_column_name (const int column) {
  if (!q_result){
    printf("pg_dbi::query_column_name: query table does not exist\n");
    return NULL;
  }

  if (column >= q_fields){
    printf("pg_dbi::query_column_name: requested value is out of query table dimensions\n");
    return NULL;
  }

  return PQfname(q_result, column);
}


int pg_dbi::transaction (int db_id, int open) {

  char cmd[128];
  if (open) sprintf(cmd, "BEGIN;");
  else      sprintf(cmd, "COMMIT;");

  PGresult *result = PQexec((PGconn*) connections[db_id], cmd);
  if (PQresultStatus(result) != PGRES_COMMAND_OK) {
    if (open) printf("starting transaction block on database \"%s\" failed with %s\n",  db_names[db_id], PQresultErrorMessage(result));
    else      printf("commiting transaction block on database \"%s\" failed with %s\n", db_names[db_id], PQresultErrorMessage(result));
    return 0;
  }

  return 1;
}

int pg_dbi::check_privilege (int db_id, const char* to_table, const char* type) {
  char check_privilege_query[2048], field_name[2048], value[2048];

  sprintf(check_privilege_query, "SELECT has_table_privilege(\'%s\', \'%s\')", to_table, type);
  printf("DEBUG: check_privilege_query string: %s\n", check_privilege_query);

  PGresult *i_result = PQexec((PGconn*) connections[db_id], check_privilege_query);
  int i_tuples   = PQntuples(i_result);
  int i_fields   = PQnfields(i_result);

  printf("DEBUG: check_privilege_query result: tuples %d, fields %d\n", i_tuples, i_fields);

  for (int i = 0; i < i_fields; i++){
    strcpy(field_name, PQfname(i_result, i));
    printf("DEBUG: field %s\n", field_name);
    if (!strcmp(field_name, "has_table_privilege")){
      strcpy(value, PQgetvalue(i_result, 0, i));
      printf("DEBUG: ... value %s\n", value);
      break;
    }
  }

  if (strcmp(value, "t")){
    printf("write access denied: table \"%s\" database \"%s\"\n", to_table, db_names[db_id]);
    return false;
  }

  return true;
}

//   for (int i = 0; i < q_fields; i++) {
//     const char *field = PQfname(q_result, i);
//     query_result[field].reserve(q_tuples);
//     for (int j = 0; j < q_tuples; j++) {
//       const char *value = PQgetvalue(q_result, j, i);
//       query_result[field].push_back(vdt(value, field));
//     }
//   }

// bool bx_dbi::check_data_presence (int db_id, const char* to_table, const table& t, const column_names& index_columns) {
//   if (!index_columns.size ()) return false;
//
//   std::ostringstream where;
//   where << "\"" << index_columns[0] << "\"=" << t[index_columns[0]][0];
//   for (unsigned int i = 1; i < index_columns.size (); i++) where << " and \"" << index_columns[i] << "\"=" << t[index_columns[i]][0];
//
//   query (db_id, to_table, where.str());
//   if (j.size ()) return true;
//
//   return false;
// }
//
// void bx_dbi::insert (int db_id, const std::string& to_table, const table& t, const column_names& index_columns, bool overwrite) {
//   const std::string &db_name = get_parameter (hash_name_indexes[db_id]).get_string ();
//
//   // Check write privileges
//   if (!m_check_priv (db_id, to_table, "INSERT")) return;
//
//   bool do_overrite = false;
//   if (m_check_data_presence (db_id, to_table, t, index_columns)) {
//     if (overwrite) do_overrite = true;
//     else {
//       get_message (bx_message::error) << "data already present in table \"" << to_table << "\"" << dispatch;
//       return;
//     }
//   }
//
//   // Open the transaction
//   if (!m_transaction (db_id, true)) return;
//
//   // Prepare the query common string
//   std::ostringstream str;
//   table::const_iterator it = t.begin ();
//   str << "INSERT INTO " << to_table << " (\"" << it->first;
//   for (it++; it != t.end (); it++) str << "\", \"" << it->first;
//   str << "\") VALUES (";
//
//   // Send queries
//   for (unsigned int i = 0; i < t.begin ()->second.size (); i++) {
//     std::ostringstream full_query;
//     it = t.begin ();
//
//     full_query << str.str () << it->second[i];
//     for (it++; it != t.end (); it++) full_query << ", " << it->second[i];
//     full_query << ")";
//
//     PGresult *result = PQexec ((PGconn *)connections[db_id], full_query.str ().c_str ());
//     if (PQresultStatus (result) != PGRES_COMMAND_OK)
//       get_message (bx_message::critic) << "query \"" << full_query.str () << "\" on database \"" << db_name << "\" failed with error \"" << PQresultErrorMessage (result) << "\"" << dispatch;
//   }
//
//   // Commit the transaction
//   if (!m_transaction (db_id, false)) return;
// }
