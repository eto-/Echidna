/* BOREXINO Reconstruction program
 *
 * Author: Michael Wurm <michael.wurm@ph.tum.de>
 * MaInt_tainer: Davide.Dangelo <Davide.Dangelo@lngs.infn.it>
 *
 * $Id: bx_laben_tof_tracker.hh,v 1.4 2010/08/10 16:18:24 wurm Exp $
 *
 * Sample module to explain echidna module writing
 * to new developers.
 * Please read the "Module programmer's guide" before 
 * or while looking at the examples here 
 * 
 */

#ifndef _BX_TOF_TRACKER_HH
#define _BX_TOF_TRACKER_HH

#include "bx_base_module.hh"
#include "bx_rec_general.hh"
#include "bx_muon_tracker.hh"

class bx_laben_event;
class TH1F;
class TH2F;
class TF1;
class TGraph;
class TGraphErrors;

/*struct coo {
	Float_t x;
	Float_t y;
	Float_t z;
	Float_t theta;
	Float_t phi;
	Float_t rad;
	Float_t time;
	Float_t dx;
	Float_t dy;
	Float_t dz;
	Float_t dtheta;
	Float_t dphi;
};*/

class bx_laben_tof_tracker: public bx_base_module {
  public:
  // this section if fixed; do not edit
    bx_laben_tof_tracker ();
    virtual ~bx_laben_tof_tracker () {}

    virtual void begin ();
    virtual bx_echidna_event* doit (bx_echidna_event *ev);
    virtual void end ();
  private:
  // this section is free; add here members, methods, histo poInt_ters
    	Int_t   i4_laben_events;
	Int_t   i4_tracked_events;
	Int_t   i4_low_nhits;
	Float_t ep_time_cut;
	Float_t pi;
	Float_t ep_max_theta_low;
	Float_t ep_max_theta_high;
	Float_t ep_max_phi_low;
	Float_t ep_max_phi_high;	
	Int_t   ep_diff_hits;
	Float_t ep_par_threshold;
	Float_t ep_theta_pos_threshold;
	Int_t   xp_plains;
	Float_t xp_theta_diff;
	Int_t   xp_theta_bins;
	Int_t   xp_nhits_diff;

	coo m_rotate_x(coo input, Float_t phi);
    	coo m_rotate_y(coo input, Float_t phi);
	coo m_rotate_z(coo input, Float_t phi);
	coo m_construct_xyz(Float_t rad, Float_t phi, Float_t theta, Float_t dphi, Float_t dtheta);
	coo m_construct_rpt(Float_t x, Float_t y, Float_t z);
	coo m_normalize_vector(coo input);
	Float_t m_limit_angle(Float_t phi, Float_t min_angle, Float_t max_angle);
	Float_t m_get_radius(coo input);
	Float_t m_get_radius(Float_t x, Float_t y, Float_t z);
};

#endif
/*
 * $Log: bx_laben_tof_tracker.hh,v $
 * Revision 1.4  2010/08/10 16:18:24  wurm
 * adjusted the xp phi fit for virginia
 *
 * Revision 1.3  2010-06-28 10:44:03  wurm
 * new laben tof tracking, modifications to muon tracking (Int_troducing phi,theta,impact to the event) and modification of global tracker to prefer laben tof over laben energy (and conditions)
 *
 * Revision 1.2  2010-05-21 13:17:37  ddangelo
 * adding laben_tof_tracker and cmt_tracker
 *
 * Revision 1.1  2010-05-21 12:34:01  ddangelo
 * bx_laben_tracker renamed as bx_laben_energy_tracker
 * added bx_laben_tof_tracker and bx_cmt_tracker
 *
 * Revision 1.6  2005-02-03 19:00:18  ddangelo
 * maInt_tainer changed (Razeto->D'Angelo)
 * Module really implemented with many possible examples for developers.
 * It features all aspects described in "Module's programmers guide" in the docs.
 *
 * Revision 1.5  2004/11/26 15:25:11  razeto
 * Thanks to Davide who so kindly underlined a syntactical error, I fixed the "MaInt_tainer" word
 *
 * Revision 1.4  2004/11/26 14:01:21  razeto
 * Added Mantainer field
 *
 * Revision 1.3  2004/09/22 13:28:37  razeto
 * Changed bx_reco_event to bx_echidna_event
 *
 * Revision 1.2  2004/03/21 18:55:05  razeto
 * Some cosmetic changes
 *
 * Revision 1.1  2004/03/20 17:39:45  pallas
 * First release
 *
 * Revision 1.1  2004/03/20 14:16:47  razeto
 * Added sources
 *
 */
