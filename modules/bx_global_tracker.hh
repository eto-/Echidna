/* BOREXINO Reconstruction program
 *
 * Authors: Werner Maneschg <werner.maneschg@mpi-hd.mpg.de>, Michael Wurm 
 *
 * $Id: bx_global_tracker.hh,v 1.5 2010/06/28 10:44:03 wurm Exp $
 *
 * Module to fit the muon track among the 4 points found by ID and OD trackers.
 *
 *  last modification: 07.05.08
*/

#ifndef _BX_GLOBAL_TRACKER_HH
#define _BX_GLOBAL_TRACKER_HH

#include "bx_rec_general.hh"
#include "bx_base_module.hh"
#include "TF1.h"
#include "TGraphErrors.h"
#include "bx_track.hh"


class bx_echidna_event;

class bx_global_tracker : public bx_base_module {
  public:
    bx_global_tracker ();
    virtual ~bx_global_tracker () {}
      
    // Operations for the framework 
    virtual void begin ();
    virtual bx_echidna_event* doit (bx_echidna_event *ev);
    virtual void end ();

  private:
    int32_t   i4_global_events;		  // counter for all muon events
    int32_t   i4_tracked_events;		// counter for reconstructed events

    int32_t   i4_enable_histos;

    int32_t   i4_min_nhits;
    float f4_max_dev;
    float f4_omega;
    float f4_chi2max44;
    float f4_max_par;

    float *m_fit_track(float *px, float *py, float *pz, float *pdx, float *pdy, float *pdz, int32_t npts, int32_t bpt); 
};

#endif

/* 
 * $Log: bx_global_tracker.hh,v $
 * Revision 1.5  2010/06/28 10:44:03  wurm
 * new laben tof tracking, modifications to muon tracking (introducing phi,theta,impact to the event) and modification of global tracker to prefer laben tof over laben energy (and conditions)
 *
 * Revision 1.4  2008-12-04 15:38:47  wurm
 * added necessary variables
 *
 * Revision 1.3  2008-07-11 15:31:26  maneschg
 * add new module for global fitting of muontrack
 *
 * Revision 1.2  2008-04-29 14:05:54  ddangelo
 * fixed the name
 *
 * Revision 1.1  2008-04-29 13:45:55  ddangelo
 * added empty
 *
 *
 */
