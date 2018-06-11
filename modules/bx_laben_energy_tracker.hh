/* BOREXINO Reconstruction program
 *
 * Authors: Michael Wurm <mwurm@ph.tum.de>, Davide Franco <davide.franco@mi.infn.it>, 
 * Maintainer: Michael Wurm <mwurm@ph.tum.de>
 *
 * $Id: bx_laben_energy_tracker.hh,v 1.1 2010/05/21 12:34:01 ddangelo Exp $
 *
 * Implemenentation of bx_laben_energy_tracker
 *
*/

#ifndef _BX_LABEN_ENERGY_TRACKER_HH
#define _BX_LABEN_ENERGY_TRACKER_HH

#include "bx_rec_general.hh"
#include "bx_base_module.hh"

class bx_echidna_event;

class bx_laben_energy_tracker : public bx_base_module {
  public:
    bx_laben_energy_tracker ();
    virtual ~bx_laben_energy_tracker () {}
      
    // Operations for the framework 
    virtual void begin ();
    virtual bx_echidna_event* doit (bx_echidna_event *ev);
    virtual void end ();

  private:
    int32_t   i4_laben_events;		// counter for all ID laben events
    int32_t   i4_tracked_events;		// counter for reconstructed events

    int32_t   i4_enable_histos;

    float f4_tau;                       // weight
    float f4_time_limit;                // take only the first hits within this limit from the beginning of cluster
    float f4_long_limit;                // limit for the hits regarded for track orientation
    int32_t   i4_nsteps;
};

#endif

/* 
 * $Log: bx_laben_energy_tracker.hh,v $
 * Revision 1.1  2010/05/21 12:34:01  ddangelo
 * bx_laben_tracker renamed as bx_laben_energy_tracker
 * added bx_laben_tof_tracker and bx_cmt_tracker
 *
 *
 */
