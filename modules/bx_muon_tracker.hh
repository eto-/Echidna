/* BOREXINO Reconstruction program
 *
 * Authors: Michael Wurm <mwurm@ph.tum.de>, Davide D'Angelo <davide.dangelo@mi.infn.it>, 
 * Maintainer: Davide D'Angelo <davide.dangelo@mi.infn.it>
 *
 * $Id: bx_muon_tracker.hh,v 1.8 2010/06/28 10:44:03 wurm Exp $
 *
 * Implemenentation of bx_laben_tracker
 *
*/

#ifndef _BX_MUON_TRACKER_HH
#define _BX_MUON_TRACKER_HH

#include "bx_rec_general.hh"
#include "bx_base_module.hh"
#include <vector>
#include "TH1F.h"

class bx_echidna_event;

struct coo {
  float x;
  float y;
  float z;
  float theta;
  float phi;
  float rad;
  float time;
  float dx;
  float dy;
  float dz;
  float dtheta;
  float dphi;
};



class bx_muon_tracker : public bx_base_module {
  public:
    bx_muon_tracker ();
    virtual ~bx_muon_tracker () {}
      
    // Operations for the framework 
    virtual void begin ();
    virtual bx_echidna_event* doit (bx_echidna_event *ev);
    virtual void end ();

  private:
    int32_t   i4_muon_events;		// counter for all ID muon events
    int32_t   i4_rec_events;		// counter for reconstructed events

    float f4_tau;			// time constant for clustering

    int32_t   i4_enable_histos;
    float f4_hit_charge_threshold;	// minimum charge for hit to be considered as "cluster"
    float f4_entry_tau;                 // suppresion of late entry points
    float f4_entry_mean;                // best time relative to ID
    float f4_entry_sigma;               // accepted deviation
    float f4_exit_incline;              // distance to time correlation of entry/exit points
    float f4_exit_dt_min;               // minimum distance
    float f4_exit_sigma;                // accepted deviation

    // methods for coordinate trafos
    coo m_rotate_x(coo input, float phi);
    coo m_rotate_y(coo input, float phi);
    coo m_rotate_z(coo input, float phi);
    coo m_construct_xyz(float rad, float phi, float theta, float dphi, float dtheta);
    coo m_construct_rpt(float x, float y, float z);
    coo m_normalize_vector(coo input);
    float m_limit_angle(float phi, float min_angle, float max_angle);
    float m_get_radius(coo input);
    float m_get_radius(float x, float y, float z);
    
    // Histograms for track variables
    TH1F *tracks_costheta;		// displays cos(theta) distribution of muons
    TH1F *tracks_phi;			// displays phi distribution
    TH1F *tracks_distance;		// displays minimum distance to the center

};

#endif

/* 
 * $Log: bx_muon_tracker.hh,v $
 * Revision 1.8  2010/06/28 10:44:03  wurm
 * new laben tof tracking, modifications to muon tracking (introducing phi,theta,impact to the event) and modification of global tracker to prefer laben tof over laben energy (and conditions)
 *
 * Revision 1.7  2008-07-21 15:00:33  wurm
 * added variables for time dependence
 *
 * Revision 1.6  2008-02-02 15:29:52  wurm
 *
 *
 * added private variables
 *
 * Revision 1.5  2008-01-30 17:56:14  wurm
 *
 * added variable "tau" to header
 *
 * Revision 1.4  2007-11-26 18:49:59  ddangelo
 * cleaned up and commented
 *
 * Revision 1.3  2007-11-26 14:06:27  ddangelo
 * completely redone
 *
 * Revision 1.2  2007-11-14 19:03:46  ddangelo
 * added michi's version.
 * writing to event
 * much more job to come
 *
 * Revision 1.1  2007-05-03 15:48:30  ddangelo
 * just added. empty skeleton
 *
 */
