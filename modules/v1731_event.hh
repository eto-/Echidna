/* BOREXINO Reconstruction program
 *
 * Author: Stefano Davini <stefano.davini@ge.infn.it>
 * Maintainer: Stefano Davini <stefano.davini@ge.infn.it>
 *
 * $Id:
 *
 * Declarations of v1731 event namespace (used in bx_v1731sys)  
 * v1731event class manages v1731 samples;
 * v1731event_getter loads and reads events from N_*.dat.gz files;
 *
 * 
 */

#ifndef _V1731_EVENT_HH
#define _V1731_EVENT_HH

#include <cstdio>
#include <string>
#include <sstream>
#include <cmath>
#include <vector>
#include <TH1F.h>
#include <TAxis.h>
#include "zlib.h"
#include "v1731_event_correlator.hh"

namespace v1731_event{
  class v1731event_decoder;
  class v1731event;
  
  using v1731_event_correlator::e_transfer;
  using v1731_event_correlator::transfer;
  using v1731_event_correlator::event_finder;
};

/************* v1731event_decoder *******************/

class v1731_event::v1731event_decoder: public bx_named{
public:
  v1731event_decoder(const event_finder&, v1731event&, unsigned u4, gzFile&, FILE*, std::string&, long&); // "random" readout (u4 is an offset in case of >1 event found in event_finder)
  v1731event_decoder(gzFile& gzfileN, v1731event&); // sequential readout
  virtual ~v1731event_decoder() {}
  bool error_occurred() const {return b_error_occurred;}

  //friend functions of v1731event
public:
  void set_time_reference(v1731event&, unsigned long u4_time_ref);
  void set_trigger_id(v1731event&, unsigned u4_trgid);
  void set_active(v1731event&, short i1_b, short i1_ch, bool b_in_active);
  void set_zle_enabled(v1731event&, bool b_in_zle_enabled);
  void set_hf_samples(v1731event&, short i1_b, short i1_ch, TH1F hf_in_samples);
  void set_length(v1731event&, short i1b, short i1_ch, long i4_in_lenght);
  void increase_num_of_good_zones(v1731event&, short i1_b, short i1_ch);
  void set_good_zones(v1731event&, short i1_b,short i1_ch, unsigned long u4_in_begin, unsigned long u4_in_length);
  void set_bins(v1731event&, short i1_b, short i1_ch, long i4_ch_size);
  void fill(v1731event&, short i1_b, short i1_ch, float f_position, float f_weigth);
  void digi_sum_fill(v1731event&, float f_position, float f_weigth);

private:
  bool b_error_occurred;
  int load_event(gzFile&, long i4_start_position, long i4_event_position);
  int read_event(gzFile&, v1731event& v1731);
  int check_zle(short i1_zle);
  short check_mask(short i1_num, int i2_mask, bool* b_active);
  enum zle_action {skip, good};
  zle_action check_zle_control_word(unsigned long u4_control_word, unsigned long& u4_num_of_data_to);

  bool b_zle_enabled;  //true if zero lenght suppression is enabled (~always)
  bool b_active_board[4];
  bool b_active_channel[4][8];
  short i1_num_active_boards;
  short i1_num_active_channel[4];
  unsigned long u4_trigger_time_tag[4];

};


/****************** v1731event ********************/


class v1731_event::v1731event: public bx_named{
  public:
  v1731event();
  virtual ~v1731event() {}
  bool get_active(short i2b, short i2ch) const {return ( ((i2b<4)&&(i2ch<8)) ? b_active[i2b][i2ch] : false);}
  bool get_zle_enabled() const {return b_zle_enabled;}
  long get_event_length(short i2b, short i2ch) const {return ( ((i2b<4)&&(i2ch<8)) ? i4_event_length[i2b][i2ch] : 0);}
  unsigned long get_time_reference() const {return u4_time_reference;}
  unsigned long get_trigger_id() const {return u4_trigger_id;}
  float get_sample_at_time(short i1_b, short i1_ch, float f_time) const 
  {return hf_samples[i1_b][i1_ch].GetBinContent(((long)(f_time/2.))+1);};
  float get_sample_at_bin(short i1_b, short i1_ch, long i4_bin) const 
  {return hf_samples[i1_b][i1_ch].GetBinContent(i4_bin);};
  void get_sample_array(short i1_b, short i1_ch, float* f_array, long i4_begin, long i4_size) const
  {for (long i4s=0; i4s<i4_size; ++i4s){
    f_array[i4s]=hf_samples[i1_b][i1_ch].GetBinContent(i4s + i4_begin +1);
  }
  };
  int get_number_of_good_zones(short i2b, short i2ch) const 
  {return ( ((i2b<4) && (i2ch<8)) ? i2_num_of_good_zones[i2b][i2ch] : 0 );};
  long get_begin_of_good_zone(short i2b, short i2ch, int i4_zone) const
  {return i4_begin_of_good[i2b][i2ch].at(i4_zone);};
  long get_length_of_good_zone(short i2b, short i2ch, int i4_zone) const
  {return i4_length_of_good[i2b][i2ch].at(i4_zone);};
  void clusterize(int base, int thr, int bck, int fwd);
  bool is_clusterized() const {return b_clusterized;}


private:
  bool b_active[4][8];
  bool b_zle_enabled;
  bool b_clusterized;
  TH1F hf_samples[4][8];
  TH1F hf_digi_sum;
  long i4_event_length[4][8];
  int i2_num_of_good_zones[4][8];
  std::vector <long> i4_begin_of_good[4][8];
  std::vector <long> i4_end_of_good[4][8];
  std::vector <long> i4_length_of_good[4][8];
  unsigned long u4_time_reference;
  unsigned long u4_trigger_id;

public:
  friend void v1731_event::v1731event_decoder::set_time_reference(v1731event&, unsigned long u4_time_ref);
  friend void v1731_event::v1731event_decoder::set_trigger_id(v1731event&, unsigned u4_trgid);
  friend void v1731_event::v1731event_decoder::set_active(v1731event&, short i1_b, short i1_ch, bool b_in_active);  
  friend void v1731_event::v1731event_decoder::set_zle_enabled(v1731event&, bool b_in_zle_enabled);
  friend void v1731_event::v1731event_decoder::set_hf_samples(v1731event&, short i1_b, short i1_ch, TH1F hf_in_samples);
  friend void v1731_event::v1731event_decoder::set_length(v1731event&, short i1b, short i1_ch, long i4_in_lenght);
  friend void v1731_event::v1731event_decoder::increase_num_of_good_zones(v1731event&, short i1_b, short i1_ch);
  friend void v1731_event::v1731event_decoder::set_good_zones(v1731event&, short i1_b, short i1_ch, unsigned long u4_in_begin, unsigned long u4_in_length);
  friend void v1731_event::v1731event_decoder::set_bins(v1731event&, short i1_b, short i1_ch, long i4_ch_size);
  friend void v1731_event::v1731event_decoder::fill(v1731event&, short i1_b, short i1_ch, float f_position, float f_weigth);
  friend void v1731_event::v1731event_decoder::digi_sum_fill(v1731event&, float f_position, float f_weigth);

};

///////////////////////////////////////////

#endif

/*
 * $Log: v1731_event.hh,v $
 * Revision 1.5  2009/11/20 10:47:23  davini
 * added method clusterize to emulate CAEN v1731 ZLE sero suppression algorithm for last long sample
 *
 * Revision 1.4  2008-11-27 14:19:35  davini
 * sequential readout in v1731event_decoder
 *
 * Revision 1.3  2008-07-30 08:55:21  davini
 * code ANSI C++ standardized
 *
 */
