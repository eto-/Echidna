/* BOREXINO Reconstruction program
 *
 * Author: Davide D'Angelo <davide.dangelo@mi.infn.it>
 * Maintainer: Davide D'Angelo <davide.dangelo@mi.infn.it>
 *             inspired by bx_calib_dark_rates
 *
 * $Id: bx_calib_muon_dark_rates.hh,v 1.2 2007/02/23 14:16:08 ddangelo Exp $
 *
 * study of the dark rates in the outer detector.
 * 
 * HOT and DEAD PMTs are identified based on the user defined thresholds 
 * The mean dark rate / PMT is calculated with and without bad lg's.  
 * Input: 
 *    the decoded hits
 * Output: 
 *    2 TH1F to barn 
 *    list and number of HOT and DEAD PMTs 
 *    DB writing
 * 
 * 2x2D histograms for monitoring the random hits, "illumination" of channels as function of theta, phi (hits_map)
 * one integral version which is sent to gviewer every 200 random triggers, is cumulative
 * one version which is reset every 1000 random triggers
 * 
*/

#ifndef _BX_CALIB_MUON_DARK_RATES_HH
#define _BX_CALIB_MUON_DARK_RATES_HH

#include "bx_rec_general.hh"
#include "bx_base_module.hh"

#include "TH1F.h"

#include <map>

class bx_echidna_event;

class bx_calib_muon_dark_rates : public bx_base_module {
  public:
    bx_calib_muon_dark_rates ();
    virtual ~bx_calib_muon_dark_rates () {}
      
    virtual void begin ();
    virtual bx_echidna_event* doit (bx_echidna_event *ev);
    virtual void end ();

  private:  
    // internal ev/hit ctrs
    int i4_trg_ctr; 
    std::vector<int> v_nhits;

    //to be read from echidna.cfg
    float f4_dark_rate_thresh_high, f4_dark_rate_thresh_low;        
    float f4_gate_start, f4_gate_end, f4_gate_width;
    
    // histograms for barn
    TH1F* dark_rate_vs_channel;
    TH1F* ok_dark_rates;
        
    // histograms for detector monitor. "illumination" of channels as function of theta, phi
  //    TH2F* random_hits_map_cumulative;   
  //    TH2F* random_hits_map;   
   
    enum pmt_status_t {
      no_pmt = 0,
      good = 1,
      hot = 2,
      dead = 3, 
      new_no_pmt = 10,
      new_good = 11,
      new_hot = 12,
      new_dead = 13,
    };
  
   std::cmap <std::string, pmt_status_t> translation_map;

};

#endif
/*
 * $Log: bx_calib_muon_dark_rates.hh,v $
 * Revision 1.2  2007/02/23 14:16:08  ddangelo
 * finished, tested, polished.
 * still missing: detector monitor histos
 *                db interface
 *
 * Revision 1.1  2007/02/22 19:57:13  ddangelo
 * added. untested.
 *
 */
