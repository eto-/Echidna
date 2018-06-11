/* BOREXINO Reconstruction program
 *
 * Author: Stefano Davini <stefano.davini@ge.infn.it>
 * Maintainer: Stefano Davini <stefano.davini@ge.infn.it>
 *
 * $Id:
 *
 * Definitons of classes declared in v1731_event_correlator; 
 *
 *
 * 
 */

#include "v1731_event_correlator.hh"
using namespace v1731_event_correlator;

const unsigned u4_NO_KEY = 0xFFFF; /* map-key meaning "no-found" */


void v1731_event_correlator::make_filename_vector(str_vec& s_filename_vector, const transfer& file_transfer, bool& b_error_occurred){
   
  gzFile gzfile = 0;
  FILE* p_open_connection = 0;

  switch (file_transfer.transfer_mode){
  case (http):{
    p_open_connection=::popen(file_transfer.s_wget_cmd.c_str(), "r");
    if (!p_open_connection) {
      std::cerr<<"bx_v1731sys: error in void v1731_event_correlator::make_filename_vector(): !p_open_connection; \n";
      b_error_occurred=true;
    }
    gzfile=::gzdopen(::fileno(p_open_connection), "r");
    break;
  }
  case (file):{
    gzfile=::gzopen(file_transfer.s_filepath.c_str(), "r");
    break;
  }
  }

  if (!gzfile){
    std::cerr<<"bx_v1731sys: error in void v1731_event_correlator::make_filename_vector(): !gzfile \n";
    b_error_occurred=true;
  }
  else {
    char c_tmp_filename [50];               // char buffer of gzgets
    std::string s_current_filename;
    while (!gzeof(gzfile)) {
      gzgets(gzfile, c_tmp_filename, 50);   // reads gzfile string
      s_current_filename = c_tmp_filename;  // char* -> std::string
      const int32_t i2_name_len= s_current_filename.length();
      if (i2_name_len>5){                   // avoids eof or other
	const std::string::iterator last_char = -- s_current_filename.end();
	const std::string::value_type c = *last_char;
	if (isspace(c)) s_current_filename.erase(last_char); // erase endl or spaces char at end
	s_filename_vector.push_back(s_current_filename);
      }
    }
  }

  gzclose(gzfile);
  if (file_transfer.transfer_mode==http) pclose(p_open_connection);
}

////////////////////////////////////////////////////////////////////

void v1731_event_correlator::make_filename_vectors(str_vec& s_filename_vector, str_vec& s_auxfilename_vector, const transfer& transfer_opt, bool& b_error_occurred){

  const std::string s_wget  = transfer_opt.s_wget_cmd + transfer_opt.s_filepath;
  const std::string s_fileN = "FileListN.list";
  const std::string s_fileT = "FileListT.list";

  const std::string s_wget_fileT = s_wget + s_fileT;
  const std::string s_wget_fileN = s_wget + s_fileN;
  const std::string s_fileN_path = transfer_opt.s_filepath + s_fileN;
  const std::string s_fileT_path = transfer_opt.s_filepath + s_fileT;

  transfer fileN_transfer;
  fileN_transfer.transfer_mode = transfer_opt.transfer_mode;
  fileN_transfer.s_wget_cmd = s_wget_fileN;
  fileN_transfer.s_filepath = s_fileN_path;

  transfer fileT_transfer;
  fileT_transfer.transfer_mode = transfer_opt.transfer_mode;
  fileT_transfer.s_wget_cmd = s_wget_fileT;
  fileT_transfer.s_filepath = s_fileT_path;

  v1731_event_correlator::make_filename_vector(s_filename_vector, fileN_transfer, b_error_occurred);
  v1731_event_correlator::make_filename_vector(s_auxfilename_vector, fileT_transfer, b_error_occurred);
}


/***********************  mapper_file_time ****************/

v1731_event_correlator::mapper_time_file::mapper_time_file(const str_vec& s_in_file_names,const str_vec& s_in_auxfile_names): bx_named("v1731_event_correlator: mapper_file_time"){

  struct_map_ft.reserve(1024);
  struct_map_ft.clear();

  /* creates a time reference (Jan 01, 2000 00:00:00) */
  tm_ref_2000.tm_year = 100;// 100 = 2000-1900
  tm_ref_2000.tm_mon  = 0;  // 0 = Jan
  tm_ref_2000.tm_mday = 1; 
  tm_ref_2000.tm_hour = 0;
  tm_ref_2000.tm_min  = 0;
  tm_ref_2000.tm_sec  = 0;
  tm_ref_2000.tm_isdst=-1;  // daylight saving enabled
  tt_ref_2000=mktime(&tm_ref_2000); // convert tm into time_t (in order to use difftime(time_t,time_t) in mapper_file_time::u4_date_to_sec(string&))

  /* read input string vector in order to create the map */
  const std::string s_Nfile_name = "N_yyMMddhhmm_yyMMddhhmm.dat.gz";
  const std::string s_Tfile_name = "T_yyMMddhhmm_yyMMddhhmm.dat";
  const int32_t i2_Ngood_length = s_Nfile_name.length();

  for (str_vec::const_iterator ci= s_in_file_names.begin(); ci!=s_in_file_names.end(); ++ci){
    str_vec::value_type s_current_name = *ci;
    const std::string::iterator last_char = -- s_current_name.end();
    const std::string::value_type c = *last_char;
    if (isspace(c)) s_current_name.erase(last_char);
    
    struct_file_time current;
    const int32_t i2_Nname_len = s_current_name.length();
    if (i2_Nname_len == i2_Ngood_length){  // length control validates the name;
      current.s_Nfile_name = s_current_name; 
      
      const std::string s_sub_date_begin=s_current_name.substr(2,10); //start chars
      const std::string s_sub_date_end=s_current_name.substr(13,10); //stop chars
      const uint32_t u4_current_time_start=u4_date_to_sec(s_sub_date_begin);
      const uint32_t u4_current_time_stop=u4_date_to_sec(s_sub_date_end);
      current.u4_time_start = u4_current_time_start;
      current.u4_time_stop  = u4_current_time_stop;
      
      const std::string s_current_auxfile_name='T'+s_current_name.substr(1,i2_Nname_len-4);   // cuts .gz
      const str_vec::const_iterator ci_aux = std::find(s_in_auxfile_names.begin(), s_in_auxfile_names.end(), s_current_auxfile_name);
      const bool b_current_has_aux=(ci_aux != s_in_auxfile_names.end() ? true : false);
      
      current.b_has_auxfile  = b_current_has_aux;
      current.s_Tfile_name = (b_current_has_aux ? s_current_auxfile_name : "");
      struct_map_ft.push_back(current);
    }
    else { // file name has wrong lenght 
    }
  }
}


uint32_t v1731_event_correlator::mapper_time_file::u4_date_to_sec(const std::string& s_in_date){
  tm tm_current;      //date and time of current begin/end file
  time_t tt_current;  //1 jan 2000 00:00:00

  const int32_t i2_current_year=atoi(s_in_date.substr(0,2).c_str());
  const int32_t i2_current_mon =atoi(s_in_date.substr(2,2).c_str());
  const int32_t i2_current_day =atoi(s_in_date.substr(4,2).c_str());
  const int32_t i2_current_hour=atoi(s_in_date.substr(6,2).c_str());
  const int32_t i2_current_min =atoi(s_in_date.substr(8,2).c_str());
  
  tm_current.tm_year= 100 + i2_current_year;
  tm_current.tm_mon = i2_current_mon-1;
  tm_current.tm_mday= i2_current_day;
  tm_current.tm_hour= i2_current_hour;
  tm_current.tm_min = i2_current_min;
  tm_current.tm_sec = 0;
  tm_current.tm_isdst=-1;  //save daylight enabled (it: ora legale/solare attiva)
  tt_current=mktime(&tm_current);

  const uint32_t u4_sec_since_2000 = static_cast<uint32_t> (difftime(tt_current,tt_ref_2000));
  return u4_sec_since_2000;
}


unsigned v1731_event_correlator::mapper_time_file::get_key_if_found(uint32_t u4_gpstime) const{
  for (unsigned u4=0; u4<struct_map_ft.size(); u4++){
    if ((u4_gpstime>=struct_map_ft[u4].u4_time_start) && (u4_gpstime<=struct_map_ft[u4].u4_time_stop)){
      return u4;
    }
  }
  return (u4_NO_KEY);  /* if no file is found, then return NO_KEY key (0xFFFF) */ 
}


void v1731_event_correlator::mapper_time_file::show(unsigned u4_key){
  get_message(bx_message::info)<<"V1731 time-file map content for key "<<u4_key<<dispatch;
  if (u4_key==u4_NO_KEY) {
    get_message(bx_message::info)<<"no file found;"<<dispatch;
    return;
  }
  if (u4_key>=struct_map_ft.size()) {
    get_message(bx_message::info)<<"key out of range;"<<dispatch;
    return;
  } 
  get_message(bx_message::info)<<"file: "<<struct_map_ft[u4_key].s_Nfile_name<<"\t start: "<<struct_map_ft[u4_key].u4_time_start<<"\t stop "<<struct_map_ft[u4_key].u4_time_stop<<dispatch;
  if (struct_map_ft[u4_key].b_has_auxfile) get_message(bx_message::info)<<"auxfile: "<<struct_map_ft[u4_key].s_Tfile_name<<dispatch;
}


/******************** v1731event_indexer ****************/

void v1731_event_correlator::v1731event_indexer::new_event(const mapper_time_file& map_tf, uint32_t u4_gpstime){
  
  const bool b_event_already_indexed = ((daq_active(u4_gpstime)) && (!event_index.empty()));
  const unsigned u4_margin = 60; // margin: look for new file to index

  if (!b_event_already_indexed){
    const unsigned u4_map_key = map_tf.get_key_if_found(u4_gpstime);
    add_to_index(map_tf, u4_map_key);
    
    /* indexing "next" file if we are close to the edge of next file */
    const uint32_t u4_next_gpstimes    = u4_gpstime + u4_margin;
    const bool b_nextevents_already_indexed = ((daq_active(u4_next_gpstimes)));
    if (!b_nextevents_already_indexed){
      const unsigned u4_map_nextkey = map_tf.get_key_if_found(u4_next_gpstimes);
      add_to_index(map_tf, u4_map_nextkey); 
    }
  }

  //AGGIUNGI: SE SIAMO TROPPO VICINO AI BORDI E "L'EVENTO SUCCESSIVO" NON E' TRA I BORDI, INDICIZZALO
  
}

void v1731_event_correlator::v1731event_indexer::add_to_index(const mapper_time_file& map_tf, unsigned u4_key){
  
  if ( (u4_key == u4_NO_KEY) || (!map_tf.get_has_auxfile(u4_key)) )
    return;
  if (u4_key>= map_tf.get_number_of_files()){
    get_message(bx_message::error)<<" map-key out of range;"<<dispatch;
    return;
  }

  const std::string s_tfile_name = map_tf.get_Tfile_name(u4_key);
  
  gzFile Tfile = 0;
  FILE* p_popen_connection = 0;

  switch (transfer_param.transfer_mode){
  case (http):{
    const std::string s_wget = transfer_param.s_wget_cmd + transfer_param.s_filepath + s_tfile_name; //" wget -nv --timeout=0 -O- http://bxmaster-data.lngs.infn.it/bxstorage/neutron/"+s_tfile_name;
    p_popen_connection=::popen(s_wget.c_str(),"r");
    if (!p_popen_connection){
      get_message(bx_message::error)<<"in void event_indexer::add_to_index: cannot do "<<s_wget<<dispatch;
      return;
    }
    
    Tfile=::gzdopen(::fileno(p_popen_connection), "r"); 
    break;
  }
  case (file):{
    Tfile=::gzopen((transfer_param.s_filepath + s_tfile_name).c_str(), "r");
    break;
  }
  }
  
  if (!Tfile) { 
    get_message(bx_message::error)<<"in void event_indexer::add_to_index: cannot open "<<s_tfile_name<<dispatch;
    pclose(p_popen_connection);
    return;
  }
  else {
    file_borders curr_borders;
    long i4_curr_position = 1; // row in T.dat (-> event num. in N_*.dat); iterated; START FROM 1
 
    while (!gzeof(Tfile)){
      event_info curr_event;
      char c_curr_daqtime[15]; // char buffer (-> timeref in T_*.dat)
      gzgets(Tfile, c_curr_daqtime, 15);
      const uint32_t u4_curr_daqtime = atol(c_curr_daqtime);
      curr_event.s_filename  = map_tf.get_Nfile_name(u4_key);
      curr_event.i4_position = i4_curr_position;
      if (u4_curr_daqtime){    /* prevents eof lines or other bad stuffs*/
	event_index.insert(std::pair<uint32_t, event_info> (u4_curr_daqtime, curr_event));
	if (i4_curr_position == 1) curr_borders.u4_start = u4_curr_daqtime;
	curr_borders.u4_stop = u4_curr_daqtime;
      }
      ++i4_curr_position;
    }
    
    curr_borders.u4_start = std::min(curr_borders.u4_start, map_tf.get_time_start(u4_key));
    curr_borders.u4_stop  = std::max(curr_borders.u4_stop, map_tf.get_time_stop(u4_key));
    index_borders.push_back(curr_borders);
    
    //get_message(bx_message::info)<<"Current Borders: start "<<curr_borders.u4_start<<" stop "<<curr_borders.u4_stop<<dispatch;
  }
  get_message(bx_message::info)<<"File "<<s_tfile_name<<" indexed; "<<event_index.size()<<" event indexed till now;"<<dispatch;

  gzclose(Tfile);
  pclose(p_popen_connection);

}


bool v1731_event_correlator::v1731event_indexer::daq_active(uint32_t u4_gpstime) const{
  for (std::vector<file_borders>::const_iterator cit = index_borders.begin();  cit != index_borders.end(); ++cit){
    if ((u4_gpstime>=(*cit).u4_start) && (u4_gpstime<=(*cit).u4_stop)) return true;
  }
  return false;
}

void v1731_event_correlator::v1731event_indexer::get_event_info(uint32_t u4_gpstime, std::vector<std::string>& s_filenames, std::vector<long>& i4_poss) const{

  s_filenames.clear();
  i4_poss.clear();

  std::multimap<uint32_t, event_info>::const_iterator cit;
  std::pair<std::multimap<uint32_t, event_info>::const_iterator, std::multimap<uint32_t, event_info>::const_iterator> ret;

  ret = event_index.equal_range(u4_gpstime);
  for (cit=ret.first; cit!=ret.second; ++cit){
    s_filenames.push_back(((*cit).second).s_filename);
    i4_poss.push_back(((*cit).second).i4_position);
  }
}

void v1731_event_correlator::v1731event_indexer::show(uint32_t u4_gpstime){
  std::multimap<uint32_t, event_info>::const_iterator cit;
  std::pair<std::multimap<uint32_t, event_info>::const_iterator, std::multimap<uint32_t, event_info>::const_iterator> ret;

  ret = event_index.equal_range(u4_gpstime);
  get_message(bx_message::info)<<"List of indexed event with gpstime_sec "<<u4_gpstime<<" : "<<dispatch;
  for (cit=ret.first; cit!=ret.second; ++cit){
    get_message(bx_message::info)<<"GpsTimeSec (key): "<<(*cit).first <<"\t| File: "<<((*cit).second).s_filename<<" | Position (event number): "<< ((*cit).second).i4_position<<dispatch;
  }
}


/************ event_finder ****************************/

v1731_event_correlator::event_finder::event_finder(const v1731event_indexer& v1731event_index, const uint32_t u4_sec):bx_named("v1731_event_correlator::event_finder"){
  
  transfer_param =  v1731event_index.get_transfer_param();
  
  std::vector <std::string> s_curr_filenames;
  std::vector <long> i4_curr_poss;
  struct_event curr_struct_event;
  events_found.reserve(4);
  events_found.clear();
  
  fill_found_events(u4_sec, 0, v1731event_index);
  fill_found_events(u4_sec, 1, v1731event_index);

}

void v1731_event_correlator::event_finder::fill_found_events(uint32_t u4_gpstime, int32_t i4_difftime, const v1731event_indexer& v1731event_index){
  
  std::vector <std::string> s_curr_filenames;
  std::vector <long> i4_curr_poss;
  struct_event curr_struct_event;
  const uint32_t u4_curr_gpstime = u4_gpstime + i4_difftime;
  
  v1731event_index.get_event_info(u4_curr_gpstime, s_curr_filenames, i4_curr_poss);
  for (unsigned u4=0; u4<i4_curr_poss.size(); ++u4){
    curr_struct_event.s_file_name   = s_curr_filenames.at(u4);
    curr_struct_event.i4_position   = i4_curr_poss.at(u4);
    curr_struct_event.u4_v1731_timer= u4_curr_gpstime; 
    curr_struct_event.i2_difference = i4_difftime;
    events_found.push_back(curr_struct_event);
  }

}


void v1731_event_correlator::event_finder::m_show() { /* debugging method */
  std::string s_found="Found ";
  std::string s_events=" event(s): ";
  std::string s_not_found="Event Not Found";
  std::string s_timer=" at (time ref.) ";
  std::string s_position=" | Position (event number) ";
  std::string s_difference=" | Difference (seconds) : ";
  
  bx_message &msg= get_message (bx_message::info);
  if (const unsigned u4_nfounds = events_found.size()) {
    msg<<s_found<<u4_nfounds<<s_events<<'\n';
    for (std::vector<struct_event>::const_iterator cit=events_found.begin(); cit!=events_found.end(); cit++ ){
      msg<<(*cit).s_file_name;
      msg<<s_timer<<(*cit).u4_v1731_timer;
      msg<<s_position<<(*cit).i4_position;
      msg<<s_difference<<(*cit).i2_difference<<'\n';
    }
  }
  else msg<<s_not_found;
  msg<<dispatch;

}

/*
 * $Log: v1731_event_correlator.cc,v $
 * Revision 1.7  2008/12/15 13:50:42  davini
 * quiet
 *
 * Revision 1.6  2008-12-15 13:37:16  davini
 * event index
 *
 * Revision 1.5  2008-07-30 10:49:02  davini
 * debug
 *
 * Revision 1.4  2008-07-30 08:55:16  davini
 * code ANSI C++ standardized
 *
 */
