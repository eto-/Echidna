/* BOREXINO Reconstruction program
 *
 * Author: Livia Ludhova <livia.ludhova@mi.infn.it>
 * 	   Alessandro Fontana <alessandro.fontana@mi.infn.it>
 * 	   Davide Franco <davide.franco@mi.infn.it>
 * Maintainer: Livia Ludhova <livia.ludhova@mi.infn.it>
 *
 * $Id: bx_calib_laben_decoding.hh,v 1.15 2012/05/09 15:13:12 ludhova Exp $
 * 
 * module to check the precalibration for each logial channel
 * to find those channels which are misaligned, too broad or having too much events outside the central interval 
 * logics:
 * 1)  channel has less hits than 1% of pulser triggers  
 * 2) if not 1) control of peak position and rms within -55 and +55 ns
 * 3) if peak and rms are ok ( check in 2) , control if there are not too many hits outside the central -55 +55 ns region
*/
 
#ifndef _BX_CALIB_LABEN_DECODING_HH
#define _BX_CALIB_LABEN_DECODING_HH

#include "bx_base_module.hh"
#include "bx_rec_general.hh"
#include "bx_laben_event.hh"
#include  "bx_dbi.hh"
#include "db_channel.hh"
#include <vector>

class bx_laben_event;
class TH1F;
class TH2F;
class TH1D;

class bx_calib_laben_decoding: public bx_base_module {
  public:
  // this section if fixed; do not edit
    bx_calib_laben_decoding ();
    virtual ~bx_calib_laben_decoding () {}

    virtual void begin ();
    virtual bx_echidna_event* doit (bx_echidna_event *ev);
    virtual void end ();
  private:
  // this section is free; add here members, methods, histo pointers
    int check_reference_lg;
    int  limit_nlg_precalib_not_ok;
    double trigger_sum;
    double peak_mean,peak_rms;
    double mean_limit_single,rms_limit_single;
    double mean_limit_total,rms_limit_total;
    bool monitor_mode;
    bool b_found_mctruth;
    
    TH1F* offset;
    TH1F* rms; 

    TH1F* laben_decoding; 
    
    TH2F* charge_vs_channel;

    TH2F* time_vs_channel;
    TH2F* time_all_vs_channel;
    TH1D* time_vs_channel_central;
    TH1D* central_one_channel;
    
    TH2F* time_vs_channel_good;
    TH1D* time_vs_channel_good_central;

  //std::vector<int> new_bad_channels_list;
  // std::vector<int> new_off_channels_list;   
  // std::vector<int> prev_bad_channels_list;
  // std::vector<int> prev_off_channels_list;   
    std::vector<int> bad_channels_list;
    std::vector<int> off_channels_list;   
    double rms_old,rms_clean;  
    double mean_old,mean_clean;

    void end_1000ev(int indx);
  
};

#endif
