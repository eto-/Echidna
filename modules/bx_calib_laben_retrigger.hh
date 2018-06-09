/* BOREXINO Reconstruction program
 *
 * Author: Livia Ludhova <Livia.Ludhova@mi.infn.it>
 * 	   
 * Maintainer: Livia Ludhova <Livia.Ludhova@mi.infn.it>
 *   
 * 1-cluster events are considered
 *
 * 1) Module that finds those lg which do retrigger (more than 1 hit from the same lg in one event)
 * it happens from time to time in all lg -> mean is calculated
 * lg where this happens "factor" (echidna.cfg) times more are identidied and the 2nd bit of "index" is set (value 10)  
 *
 * 2) retriggering happnes about 140 ns after the fisrt hit.... so in the clustered_hit_time distribution there is
 * a peak in time region 150-200 ns
 * the mean total_ratio of N-bad-events (in this retrigger region) / N-good-events (first 30ns) is calculated
 * this ratio and error_ratio is calculated also for each lg individually
 * lg where this ratio is too high (ratio - total_ratio)/error_ratio > n_sigma (from echidna.cfg) identified
 * and 1st bit of the  "index" is set (01) 
 * INPUT: clusterd hits
 * OUTPUT: list of lg with index > 0
 * index = 11 (retriggering and bad hit--time distribution)
 * index = 10 (retriggering)
 * index = 01 (bad hit--time distribution) (it happens that some lg have the main peak shifted (even if precalib is ok), * or the ratio of the peaks due to the light reflection is anomalous etc....)
 */

#ifndef _BX_CALIB_LABEN_RETRIGGER_HH
#define _BX_CALIB_LABEN_RETRIGGER_HH

#include "bx_base_module.hh"
#include "bx_rec_general.hh"
#include "TROOT.h"
#include "TH1F.h"
#include "TH2F.h"
#include <map>


class bx_laben_event;
class TH1F;
class TH2F;

class bx_calib_laben_retrigger: public bx_base_module {
  public:
  // this section if fixed; do not edit
    bx_calib_laben_retrigger ();
  virtual ~bx_calib_laben_retrigger () {}

    virtual void begin ();
    virtual bx_echidna_event* doit (bx_echidna_event *ev);
    virtual void end ();
  private:
    
   //to be read from echidna.cfg
  double n_sigma;
  double factor;
  int    index;

    //for each lg time distribution of clustered hits, time 0 is the beginning of the cluster 
  TH2F*  time_vs_lg; 
    // lg distrubution for event by event (is reset for each event)
  TH1F*  lg_junk;    
    // lg distribution of retriggering 
  TH1F*  lg_retrig;  
  };

#endif
 






