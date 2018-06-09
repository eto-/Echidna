/* BOREXINO Reconstruction program
 *
 * Author: Stefano Davini <stefano.davini@ge.infn.it>
 * Maintainer: Stefano Davini <stefano.davini@ge.infn.it>
 *
 * $Id:
 *
 * Definitons of v1731_event namespace (used in bx_v1731sys)  
 * 
 *
 *
 * 
 */

#include <algorithm>
#include <TH1F.h>
#include "v1731_event.hh"

using v1731_event_correlator::event_finder;
using namespace v1731_event;


/***************** v1731event_decoder ********************/

v1731_event::v1731event_decoder::v1731event_decoder(const event_finder& ef_find_event, v1731event& v1731ev, unsigned u4, gzFile& currNfile, FILE* p_popen_connection, std::string& s_curr_file_name, long& i4_curr_file_position): bx_named("v1731event_decoder"){
  
  const std::string s_file_name = ef_find_event.get_filepath() + ef_find_event.get_file_name(u4);
  const long i4_ev_position = ef_find_event.get_position(u4);
  set_time_reference(v1731ev, ef_find_event.get_v1731_timer(u4));

  b_error_occurred = false;
  int i2_deco_err  = 0;

  const bool b_seq_rd = ((s_file_name==s_curr_file_name) && (i4_ev_position>=i4_curr_file_position) && (currNfile));

  /*get_message(bx_message::info)<<"File name: "<<s_file_name<<"\t CF name: "<<s_curr_file_name<<"\t Position: "<<i4_ev_position<<"\t C Position: "<<i4_curr_file_position<<dispatch;
  if (b_seq_rd) get_message(bx_message::info)<<"Sequential"<<dispatch;
  if (!b_seq_rd) get_message(bx_message::info)<<"Not Sequential"<<dispatch;*/

  switch (ef_find_event.get_load_mode()){
  
  case (v1731_event_correlator::http):{ /* "Gran Sasso" */
    if (!b_seq_rd){          /* Nfile has to be closed and reopened */
      if (currNfile) {
	gzclose(currNfile); 
	get_message(bx_message::info)<<"Nfile "<<s_curr_file_name<<" closed;"<<dispatch;
      }
      if (p_popen_connection){
	pclose(p_popen_connection);
      }
      const std::string s_wget = ef_find_event.get_wget_cmd() + s_file_name;//"wget -nv --timeout=60 -O- http://bxmaster-data.lngs.infn.it/bxstorage/neutron/"+s_file_name;
      p_popen_connection=::popen(s_wget.c_str(),"r");
      if (!p_popen_connection) {
	get_message(bx_message::error)<<"Error in v1731event_decoder::v1731event_decoder (FILE* p_popen_connection) "<<s_file_name<<dispatch;
	b_error_occurred = true;
	i2_deco_err = -2;
      }
      currNfile=::gzdopen(::fileno(p_popen_connection),"r");
      s_curr_file_name = s_file_name;
      i4_curr_file_position = 0;
      if (currNfile) get_message(bx_message::info)<<"Nfile "<<s_curr_file_name<<" opened;"<<dispatch;
    }
    break;  
  }

  case (v1731_event_correlator::file):{ /* "Genova" */
    if (!b_seq_rd){
      if (currNfile) gzclose(currNfile);
      currNfile=::gzopen(s_file_name.c_str(), "r");
      s_curr_file_name = s_file_name;
      i4_curr_file_position = 0;
    }
    break;
  }
  }
  
  if (!currNfile){
    get_message(bx_message::error)<<"in v1731event_decoder::v1731event_decoder gzFile "<<s_file_name<<" is zombie;"<<dispatch;
    b_error_occurred=true;
    i2_deco_err = -3;
  }
 
  if (!i2_deco_err){
    const long i4_start_position = ( b_seq_rd ? (i4_curr_file_position-1): 0); // positions from v1731_event_finder start from 1
    i2_deco_err = load_event(currNfile, i4_start_position, i4_ev_position);    // load_event return (-n) if error occurred
    if ((b_error_occurred = i2_deco_err)) {
      get_message(bx_message::error)<<"in v1731event_decoder::load_event "<<s_file_name<<" event "<<i4_ev_position<<dispatch;
    }
    else{
      i4_curr_file_position = i4_ev_position;
      i2_deco_err = read_event(currNfile, v1731ev); // read_event return (-n) if error occurred
      if (!(b_error_occurred = i2_deco_err)){
	++i4_curr_file_position;
      }
      else {  /* error in reading */
	get_message(bx_message::error)<<"in v1731event_decoder::read_event "<<s_file_name<<" event "<<i4_ev_position<<dispatch;
      }
    }
  }

  if (b_error_occurred){   /* if error, close Nfile and pipe connection, reset status vars*/
    switch (i2_deco_err){
    case (-3): {
      get_message(bx_message::error)<<"gzFile or pipe was not opened correcly;"<<dispatch;
      break;
    }
    case (-33):{
      get_message(bx_message::error)<<"in v1731_event_decoder::load_event, fist word is not 3333;"<<dispatch;
      break;
    }
    case (-66):{
      get_message(bx_message::error)<<"in v1731_event_decoder::load_event, last word is not 6666;"<<dispatch;
      break;
    }
    case (-330):{
      get_message(bx_message::error)<<"in v1731_event_decoder::read_event, first word is not 3333;"<<dispatch; 
      break;
    }
    case (-400):{
      get_message(bx_message::error)<<"in v1731_event_decoder::read_event, ZLE control word is not valid;"<<dispatch; 
      break;
    }
    case (-660):{
      get_message(bx_message::error)<<"in v1731_event_decoder::read_event, last  word is not 6666;"<<dispatch; 
      break;
    }
    default:{
      get_message(bx_message::error)<<"unknown error;"<<dispatch;
    }
    }
    
    gzclose(currNfile);
    currNfile = 0;
    s_curr_file_name = "";
    i4_curr_file_position = 0;
    if (ef_find_event.get_load_mode() == v1731_event_correlator::http) pclose(p_popen_connection);
    get_message(bx_message::warn)<<"current gzFile (and pipe) closed;"<<dispatch;
    
  }
}


v1731_event::v1731event_decoder::v1731event_decoder(gzFile & gzfileN, v1731event & v1731ev):bx_named("v1731event_decoder"){
  b_error_occurred=false;

  if (!gzfileN) {
    b_error_occurred=true;
    get_message(bx_message::error)<<"Error occurred in v1731event_decoder::v1731event_decoder(gzFile&, v1731event&)"<<dispatch;
  }
  else{
    const int i2_read_err = read_event(gzfileN, v1731ev);
    if ((b_error_occurred = i2_read_err)) {
      get_message(bx_message::error)<<"in v1731event_decoder::read_event"<<dispatch;
    }
  }
}


int v1731_event::v1731event_decoder::load_event(gzFile & Nfile, long i4_start_position, long i4_event_position){ 
  for (long i4_i=i4_start_position; i4_i<i4_event_position-1; ++i4_i){// start_position starts from 0; event_position starts from 1 (see class constructor); 
    char c_buff [20];  // char buffer
    int i2_chk;        // first and last word
    gzgets(Nfile, c_buff, 20);
    i2_chk=atoi(c_buff);
    if (i2_chk!=3333){ 
      return (-33); // error occurred    
    }
    gzgets(Nfile, c_buff, 20);
    const long i4_length=atol(c_buff);
    for (long i4_j=0; i4_j<i4_length; ++i4_j){
      gzgets(Nfile, c_buff, 20);   // streams the event
    }
    gzgets(Nfile, c_buff, 20);
    i2_chk=atoi(c_buff);
    if (i2_chk!=6666) {
      return (-66);  // error occurred
    }
  }
  return 0;
}


int v1731_event::v1731event_decoder::check_zle(short i1_zle){
  switch (i1_zle){
  case (1):{
    b_zle_enabled = true;
    return 0;
    break;
  }
  case (0):{
    b_zle_enabled = false;
    return 0;
    break;
  }
  default:
    return (-400);
  }
}

short v1731_event::v1731event_decoder::check_mask(short i1_num, int i2_mask, bool* b_active){
  short i1_num_active=0;
  for (short i1_i=0; i1_i<i1_num; ++i1_i){
    b_active[i1_i]=i2_mask%2;
    if (b_active[i1_i]) ++i1_num_active;
    i2_mask=i2_mask/2;
  }
  return i1_num_active;
}

v1731_event::v1731event_decoder::zle_action v1731_event::v1731event_decoder::check_zle_control_word(unsigned long u4_control_word, unsigned long& u4_num_of_data_to){
  const unsigned long u4_2pow31=0x80000000;
  if (u4_control_word>u4_2pow31) {
    u4_num_of_data_to=(u4_control_word-u4_2pow31)*4; //control word contains the number of L32 to be written; 1L32=4byte -> 4datum (v1731 is 8bit digitizer)
    return (good);
  }
  else {
    u4_num_of_data_to=(u4_control_word*4);
    return (skip);
  }
} 

int v1731_event::v1731event_decoder::read_event(gzFile& Nfile, v1731event & v1731ev){
  #define NBUFF 20
  char c_buff [NBUFF]; // char buffer
  
  gzgets(Nfile, c_buff, NBUFF);
  int i2_tmp=atoi(c_buff);// first word:3333
  if (i2_tmp!=3333) return (-330); //error has occurred in reading 
  gzgets(Nfile, c_buff, NBUFF);
  //const long i4_event_len=atol(c_buff); //length (number of lines) of event
  gzgets(Nfile, c_buff, NBUFF);
  const unsigned u4_trgid=atol(c_buff); //TRGID; 
  gzgets(Nfile, c_buff, NBUFF);
  const unsigned long u4_LV_time=atol(c_buff); //event time reference (seconds) provided by Neutron DAQ (PC & LabView clock time)
  gzgets(Nfile, c_buff, NBUFF);
  //const int i2_ev_num=atoi(c_buff); //event number; not in use now
  //get_message(bx_message::info)<<"Event Number read: "<<i2_ev_num<<dispatch;
  gzgets(Nfile, c_buff, NBUFF);
  const short i1_zle=atoi(c_buff); // 0 if not zle, 1 if zle
  gzgets(Nfile, c_buff, NBUFF);
  const int i2_board_mask=atoi(c_buff); // board mask (read it in binary repr.)
  gzgets(Nfile, c_buff, NBUFF);
  //const long i4_data_sum=atol(c_buff); // number of data (all board, all channel); useless?

  const int i2_zle_err = check_zle(i1_zle);  //v1731::check_zle changes private bool b_zle_enabled; returns (-400) if error occurred
  if (i2_zle_err) return (i2_zle_err);
  if (b_zle_enabled) set_zle_enabled(v1731ev,true);
  else set_zle_enabled(v1731ev,false);

  set_trigger_id(v1731ev, u4_trgid);
  set_time_reference(v1731ev, u4_LV_time);

  i1_num_active_boards=check_mask(4,i2_board_mask,b_active_board);

  for (short i1_b=0; i1_b<4; ++i1_b){ //goes through v1731 boards 

    if (b_active_board[i1_b]) {
      gzgets(Nfile, c_buff, NBUFF);
      //const long i4_ev_counter=atol(c_buff); // event counter; useless?
      gzgets(Nfile, c_buff, NBUFF);
      const int i2_channel_mask=atoi(c_buff);  // channel mask
      gzgets(Nfile, c_buff, NBUFF);
      //const double f8_trigger_time_tag=atof(c_buff); //trigger time tag
      //const unsigned long u4_trigger_time_tag= static_cast<unsigned long> (f8_trigger_time_tag);
      gzgets(Nfile, c_buff, NBUFF);
      const long i4_ev_size=atol(c_buff); //event size
      
      i1_num_active_channel[i1_b]=check_mask(8,i2_channel_mask,b_active_channel[i1_b]);

      if (b_zle_enabled){
	long u4_current_rel_pos=0;  //relative data position (in #samples)
	for (short i1_ch=0; i1_ch<8; ++i1_ch){  //goes through channels of single v1731 board
	  
	  u4_current_rel_pos=0; // event start: relative data position = 0 
	 
	  const bool b_channel_active=b_active_channel[i1_b][i1_ch];
	  set_active(v1731ev,i1_b,i1_ch,b_channel_active);

	  if (b_channel_active){  //active channel: read it
	    
	    gzgets(Nfile, c_buff, NBUFF);
	    const long i4_ch_size=atol(c_buff); //channel size: number of data to stream for THIS channel
	    long i4_s=0;  //i4_s increase every stream/gzgets (after channel size)
	    while (i4_s<i4_ch_size){

	      gzgets(Nfile, c_buff, NBUFF);
	      const double f8_control_word=atof(c_buff);   //WARNING: unsigned long
	      const unsigned long u4_control_word= static_cast<unsigned long> (f8_control_word);
	      ++i4_s;

	      unsigned long u4_num_of_data; // number of data to read (or number of data skipped)
	      const zle_action action_to_do=check_zle_control_word(u4_control_word,u4_num_of_data);
	      
	      switch (action_to_do){
	      case (skip): {
		u4_current_rel_pos+=(u4_num_of_data);  
	      }
		break;
	      case (good): {
		const unsigned long u4_final_rel_pos= u4_current_rel_pos + u4_num_of_data;
		for (unsigned long u4_g=u4_current_rel_pos; u4_g<u4_final_rel_pos; ++u4_g ){
		  gzgets(Nfile, c_buff, NBUFF);
		  float f_data=atof(c_buff); //data (i.e. v1731 samples)
		  ++i4_s;
		  fill(v1731ev, i1_b, i1_ch, 2.*u4_g, f_data);
		  if (!((i1_b==0) && (i1_ch==0)))digi_sum_fill(v1731ev, 2.*u4_g, f_data);
		}
		set_good_zones(v1731ev,i1_b,i1_ch, 2*u4_current_rel_pos, 2*u4_num_of_data);
		increase_num_of_good_zones(v1731ev,i1_b,i1_ch);
		u4_current_rel_pos+=(u4_num_of_data);
	      }
		break;
	      }	    
	    }
	  }
	  else { //channel not active (in zle enabled)
	  }

	  set_length(v1731ev, i1_b, i1_ch, 2*u4_current_rel_pos);
	  set_bins(v1731ev, i1_b, i1_ch, u4_current_rel_pos);
	}
      }
      else {   // zle not enabled
	if (i1_num_active_channel[i1_b]==0) return (-1); //avoids division by zero
	const long i4_ch_size=i4_ev_size/i1_num_active_channel[i1_b]; // channel size: number of data to stream for ALL channels
	for (short i1_ch=0; i1_ch<8; ++i1_ch){

	  const bool b_channel_active=b_active_channel[i1_b][i1_ch];
	  set_active(v1731ev, i1_b, i1_ch, b_channel_active);
	  
	  if (b_channel_active){
	    set_length(v1731ev, i1_b, i1_ch, 2*i4_ch_size); // bug fixed: (i4_ch_size)
	    for (long i4_n=0; i4_n<i4_ch_size; ++i4_n){
	      gzgets(Nfile, c_buff, NBUFF);
	      float f_data=atof(c_buff); //data (v1731 samples)
	      fill(v1731ev, i1_b, i1_ch, 2.*i4_n, f_data);
	      if (!((i1_b==0) && (i1_ch==0))) digi_sum_fill(v1731ev, 2.*i4_n, f_data);
	    }
	  }
	  else if (!b_channel_active){
	    set_length(v1731ev, i1_b, i1_ch, 0);
	  }       
	  set_bins(v1731ev, i1_b, i1_ch,i4_ch_size);
	}
      }
    }
    //goto next board
  }

  gzgets(Nfile, c_buff, NBUFF);
  i2_tmp=atoi(c_buff);         // last word: 6666
  if (i2_tmp==6666) return 0;  // correct readout
  else return (-660);            // incorrect readout
}





void v1731_event::v1731event_decoder::set_time_reference(v1731event& v1731ev, unsigned long u4_time_ref){
  v1731ev.u4_time_reference = u4_time_ref;
}

void v1731_event::v1731event_decoder::set_trigger_id(v1731event& v1731ev, unsigned u4_trgid){
  v1731ev.u4_trigger_id = u4_trgid;
}

void v1731_event::v1731event_decoder::set_active(v1731event& v1731ev, short i1_b, short i1_ch, bool b_in_active){
  v1731ev.b_active[i1_b][i1_ch] = b_in_active;
}

void v1731_event::v1731event_decoder::set_zle_enabled(v1731event& v1731ev, bool b_in_zle_enabled){
  v1731ev.b_zle_enabled = b_in_zle_enabled;
}

void v1731_event::v1731event_decoder::set_hf_samples(v1731event& v1731ev, short i1_b, short i1_ch, TH1F hf_in_samples){
  v1731ev.hf_samples[i1_b][i1_ch] = hf_in_samples;
}

void v1731_event::v1731event_decoder::set_length(v1731event& v1731ev, short i1_b, short i1_ch, long i4_in_length){
  v1731ev.i4_event_length[i1_b][i1_ch] = i4_in_length;
}

void v1731_event::v1731event_decoder::increase_num_of_good_zones(v1731event& v1731ev, short i1_b, short i1_ch){
  ++v1731ev.i2_num_of_good_zones[i1_b][i1_ch];
}

void v1731_event::v1731event_decoder::set_good_zones(v1731event& v1731ev, short i1_b, short i1_ch, unsigned long u4_in_begin, unsigned long u4_in_length){
  v1731ev.i4_begin_of_good[i1_b][i1_ch].push_back(u4_in_begin);
  v1731ev.i4_length_of_good[i1_b][i1_ch].push_back(u4_in_length);
  v1731ev.i4_end_of_good[i1_b][i1_ch].push_back(u4_in_begin+ u4_in_length);
}

void v1731_event::v1731event_decoder::set_bins(v1731event& v1731ev, short i1_b, short i1_ch, long i4_ch_size){
  v1731ev.hf_samples[i1_b][i1_ch].SetBins(i4_ch_size, -1, 2.* i4_ch_size -1);
}

void v1731_event::v1731event_decoder::fill(v1731event& v1731ev, short i1_b, short i1_ch, float f_position, float f_weigth){
  v1731ev.hf_samples[i1_b][i1_ch].Fill(f_position, f_weigth);
}

void v1731_event::v1731event_decoder::digi_sum_fill(v1731event& v1731ev, float f_position, float f_weigth){
  v1731ev.hf_digi_sum.Fill(f_position, f_weigth);
}



/******************** v1731event ******************/

v1731_event::v1731event::v1731event(): bx_named("v1731event"){
  const long i4_max_range = static_cast<long>(1.1e6);// 1.1 ms
  const long i4_bin_number= i4_max_range/2;          // 1bin=2 ns
  
  b_clusterized = false;
  
  hf_digi_sum.SetName("v1731event_digi_sum");
  hf_digi_sum.SetTitle("v1731event Digital Sum; Time (ns);;");
  hf_digi_sum.SetBins(i4_bin_number, -1, i4_max_range -1);

  for (short i1_b=0; i1_b<4;++i1_b){
    for (short i1_ch=0; i1_ch<8; ++i1_ch){
      b_active[i1_b][i1_ch]=false;
      i2_num_of_good_zones[i1_b][i1_ch]=0;

      std::ostringstream oss_histo_name;
      std::ostringstream oss_histo_title;
      oss_histo_name<<"board:"<<i1_b<<"ch:"<<i1_ch;
      oss_histo_title<<"v1731 Samples board "<<i1_b<<" channel "<<i1_ch<<"; Time (ns); ADC bin;";
      const std::string s_histo_name=oss_histo_name.str();
      const std::string s_histo_title=oss_histo_title.str();
      hf_samples[i1_b][i1_ch].SetName(s_histo_name.c_str());
      hf_samples[i1_b][i1_ch].SetTitle(s_histo_title.c_str());
      hf_samples[i1_b][i1_ch].SetBins(i4_bin_number, -1, i4_max_range -1);
    }
  }
}


void v1731_event::v1731event::clusterize(int i4_base, int i4_thr, int i4_bck, int i4_fwd){
	// CAEN v1731 board include a zero suppression algorithm
	// however, this algorithm works only if the numer of transitions good/zero
	// is <14 : all samples belonging over the 14th transitions are acquired
	// 14 transitions -> num of good zones <8, 8th zone is the long last one
	// see S.Davini Master Thesis (2008) for detalis
	// this method is trying to "simulate" v1731 Zero Length Encoding algorithm

	if (b_clusterized)
		return;

	const int nshortMAX = 7;
	// max number of "short" good zone; note that .at(nshortMAX) refers to the "long last zone"
	// remind for future: arrays start from 0: 1 zone => .size =1, .at(0) [not .at(1)]

	if (i2_num_of_good_zones[0][0]<= nshortMAX){
		b_clusterized = true;
		return;
	}

	// check if this zone is long?
	if (i4_length_of_good[0][0].at(nshortMAX)<5000){
		b_clusterized = true;
		return;
	}
		
	// here is the last "long" good zone
	
	const long i4_lbegin = i4_begin_of_good[0][0].at(nshortMAX);
	const long i4_lend   = i4_end_of_good[0][0].at(nshortMAX);
	
	// now i delete long zone info: i will fill them later
	i4_begin_of_good [0][0].pop_back();
	i4_end_of_good   [0][0].pop_back();
	i4_length_of_good[0][0].pop_back();
	i2_num_of_good_zones[0][0]--;

	long i4_tmp_begin = i4_lbegin, i4_tmp_end =i4_lend;

	for (long i4t=i4_lbegin; i4t<=i4_lend-2; i4t+=2){
		const int fcurr_sample = int(get_sample_at_time(0, 0, float(i4t)));
		if (fcurr_sample> i4_thr){
			i4_tmp_begin = std::min(i4_tmp_begin, i4t - i4_bck);
			i4_tmp_end   = i4t + i4_fwd;
		}

		if (i4t > i4_tmp_end){
			// cluster found
			i4_tmp_begin = std::max(i4_tmp_begin, i4_lbegin);
			i4_tmp_end   = std::min(i4_tmp_end, i4_lend); 
			i4_begin_of_good [0][0].push_back(i4_tmp_begin);
			i4_end_of_good   [0][0].push_back(i4_tmp_end);
			i4_length_of_good[0][0].push_back(i4_tmp_end - i4_tmp_begin);
			i2_num_of_good_zones[0][0]++;
			i4_tmp_begin = i4_lend;
			i4_tmp_end   = i4_lend;
		}
		
	}
	b_clusterized = true;
	
}

/*
 * $Log: v1731_event.cc,v $
 * Revision 1.5  2009/11/20 10:47:23  davini
 * added method clusterize to emulate CAEN v1731 ZLE sero suppression algorithm for last long sample
 *
 * Revision 1.4  2008-11-27 16:38:06  davini
 * new error code in v1731event_decoder;
 *
 * Revision 1.3  2008-11-27 14:19:35  davini
 * sequential readout in v1731event_decoder
 *
 *
 */
