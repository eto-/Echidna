/* BOREXINO Reconstruction program
 *
 * Author: Daniela Manuzio <dmanuzio@ge.infn.it> and Livia Ludhova<livia.ludhova@mi.infn.it> 
 * Maintainer: Livia Ludhova<livia.ludhova@mi.infn.it>
 *
 * $Id: bx_calib_laben_dark_rates.hh,v 1.2 2008/08/20 14:14:02 ddangelo Exp $
 *
 * study of the dark rates in the inner detector.
 * 
 * The dark rate for each lg (and for the whole detector, e.g. mean value per PMT) is
 * calculated by counting the number of hits considering individual trigger types and their sum 
 * (which trigger types contribute to sum is given by index_trg_type set in echidna.cfg)
 * random (whole gate)
 * pulser (before the pulser peak)
 * laser (before the laser peak)
 * 1-cluster neutrino events (after the first cluster hits)
 * trigger types given by index_trg_type
 * for the whole detector the 3-sigma compatibility of individual trigger types with "all trigger types" is tested
 *
 *
 * HOT and DEAD PMTs (with/without cone) are identified based on the trigger sum-value 
 * The mean dark rate / PMT is calculated with and without bad lg's.  
 * Input: the decoded hits (after clustering, 1-cluster events have to be found) 
 * Output: 
 * 1x2D histo: dark rates vs lg (x-axis), y-values, trigger type:  1: random, 2: pulser, 3: laser 4: neutrino-1cluster events, 5:sum of all trigger types
 * list and number of (with/without cone) HOT and DEAD PMTs 
 * number of different trigger types and corresponding dark time
 * 
 * 2x2D histograms for monitoring the random hits, "illumination" of channels as function of theta, phi (hits_map)
 * one integral version which is sent to gviewer every 200 random triggers, is cumulative
 * one version which is reset every 1000 random triggers
 * 
*/

#ifndef _BX_CALIB_LABEN_DARK_RATES_HH
#define _BX_CALIB_LABEN_DARK_RATES_HH
#include "bx_rec_general.hh"
#include "bx_base_module.hh"
#include "TROOT.h"
#include "TH1F.h"
#include "TH2F.h"
#include <map>

class bx_echidna_event;

class bx_calib_laben_dark_rates : public bx_base_module {
  public:
    bx_calib_laben_dark_rates ();
    virtual ~bx_calib_laben_dark_rates () {}
      
    // Operations for the framework 
    virtual void begin ();
    virtual bx_echidna_event* doit (bx_echidna_event *ev);
    virtual void end ();

  private:
    
    int count_random_triggers; 
    int count_pulser_triggers; 
    int count_laser_triggers; 
    int count_neutrino_triggers; 
    

      // used to define the range od dark-times histos
    double time_after_gate; 
    double time_before_gate; 
    int time_bins;  //number of bins in dark-times histos (50 ns binning)
          
      //to be read from echidna.cfg
    int index_trg_type;
    float min_dark_time; 
    float dark_rate_thresh_high, dark_rate_thresh_low;  
      
      //to be read from DB
    float gate_width;
    float trigger_offset;
    float time_pulser;   //edge in ns, 0 - start of the gate
    float time_laser;   //edge in ns, 0 - start of the gate
    float time_neutrino;   //edge in 1cluster events - decoded hits time distribution, in ns, 0 - start of the gate
    float width_neutrino;   //width of 1-cluster events 
    
    //on y axis:
       //1: random, 2: pulser, 3: laser 4: neutrino-1cluster events, 5:sum of all trigger types
    TH2F* ID_dark_rate_vs_channel;         //for non ordinary channels set to -1
    TH1F* ID_ok_dark_rates;                //dark rates of good PMTs
    TH1F* ID_all_dark_rates;                //dark rates of good PMTs
  
      
    //histograms for monitoring the random hits, "illumination" of channels as function of theta, phi
    TH2F* random_hits_map_cumulative;   
    TH2F* random_hits_map;   


    std::vector<int> nhits_random;
    std::vector<int> nhits_pulser;
    std::vector<int> nhits_laser;
    std::vector<int> nhits_neutrino;
   
    float fit_rate(TH1D* h, int npmts);
    float fit_error_rate(TH1D* h, float rate);
  
    enum pmt_status {
      no_pmt = 0,
      good = 1,
      hot = 2,
      dead = 3, 
      new_no_pmt = 10,
      new_good = 11,
      new_hot = 12,
      new_dead = 13,
    };
  
   std::cmap <std::string, pmt_status> translation_map;

};

#endif
/*
 * $Log: bx_calib_laben_dark_rates.hh,v $
 * Revision 1.2  2008/08/20 14:14:02  ddangelo
 * class and methods renamed
 *
 * Revision 1.1  2008-08-20 14:06:36  ddangelo
 * module bx_laben_dark_rates renamed as bx_calib_laben_dark_rates
 *
 * Revision 1.10  2007-05-31 15:08:07  ludhova
 * histo new names + 1 new histo
 *
 * Revision 1.9  2006-12-12 11:07:40  ludhova
 * some structural changes
 *
 * Revision 1.8  2006-08-30 09:43:09  ludhova
 * new visitors
 *
 * Revision 1.7  2006/07/18 17:38:13  ludhova
 * added writing to DB
 *
 * Revision 1.6  2006/03/22 12:43:54  ludhova
 * dark noise is calcluated based on 4 trigger types
 *
 * Revision 1.5  2006/03/02 16:45:51  ludhova
 * new outputs prepared for the db tables
 *
 * Revision 1.4  2006/02/28 15:25:42  ludhova
 * New header file
 *  ----------------------------------------------------------------------
 *
 * Revision 1.3  2004/11/26 15:25:11  razeto
 * Thanks to Davide who so kindly underlined a syntactical error, I fixed the "Maintainer" word
 *
 * Revision 1.2  2004/11/26 14:16:54  razeto
 * Added Mantainer field
 *
 * Revision 1.1  2004/10/05 13:33:28  dmanuzio
 * Removed a bug and changed the name of 2 modules (for consistency with
 * bx_calib_laben_fiber_bundle module)
 *
 * Revision 1.2  2004/09/28 13:31:13  dmanuzio
 * Removed the ifdef for root barn histos
 *
 * Revision 1.1  2004/09/28 09:51:32  dmanuzio
 * Added the first version of the module to compute the dead pmts and the
 * pmts' dark rates
 *
 *
 */
