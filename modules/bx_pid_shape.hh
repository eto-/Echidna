/* BOREXINO Reconstruction program
 *
 * Author: Sandra Zavatarelli <zavatare@ge.infn.it> and Livia Ludhova <livia.ludhova@mi.infn.it>
 * Maintainer: Livia Ludhova <livia.ludhova@mi.infn.it>
 *
 * Event's asymmetry tester
 * 1) looks for strange single-cluster hit distributions in crates, fe-boards (feb) and laben-boards, as:
 * 1a) only small portion of crates (feb, laben) has triggered
 * 1b) some crates, boards, channels are too hot
 */

#ifndef _BX_PID_SHAPE_HH
#define _BX_PID_SHAPE_HH

#include <vector>
#include "bx_base_module.hh"
#include "TVector3.h"
#include "TF1.h"
#include "TF2.h"
#include "TH1F.h"
#include "TH2F.h"

class bx_echidna_event;

class bx_pid_shape: public bx_base_module {
  public:
    bx_pid_shape ();
    virtual ~bx_pid_shape () {}
    
    // Operations for the framework 
    virtual void begin ();
    virtual bx_echidna_event* doit (bx_echidna_event *ev);
    virtual void end ();

  private:
    int i4_min_nhits; // from echidna.cfg
    int i4_nbins;  // from echidna.cfg
    int event_debug;  // from echidna.cfg
    float f4_meantime;  // from echidna.cfg
    float f4_peaktime;  // from echidna.cfg


      // number of:
    int n_good_lg;
    int n_on_crate;
    int n_on_feb;
    int n_on_laben;

      //number of good_lg in individual crate, feb, laben, lg...
    std::vector<double> crate_weight;
    std::vector<double> feb_weight;
    std::vector<double> laben_weight;
    std::vector<double> lg_weight;

      //event numbers of identified events with strange hit distrib. in carte, feb and laben
    std::vector<double> ev_el;
       
      //histo for electronics
    TH2F *percentage_of_nhits_crate;
    TH2F *percentage_of_nhits_feb;
    TH2F *percentage_of_nhits_laben;
    TH2F *percentage_of_nhits_lg;

    TH2F *nhitted_crates_vs_nhits;
    TH2F *nhitted_feb_vs_nhits;
    TH2F *nhitted_laben_vs_nhits;
    TH2F *nhitted_lg_vs_nhits;

    TH1F* suspicious_events;
  
    TH1F* h_crate_weight;
    TH1F* h_feb_weight; 
    TH1F* h_laben_weight;
    TH1F* h_lg_weight;

    TH1F* crate_all_weighted;
    TH1F* feb_all_weighted;
    TH1F* laben_all_weighted;

    TH1F* crate_all;
    TH1F* feb_all;
    TH1F* laben_all;

    TH1F* crate_1cl;
    TH1F* feb_1cl;
    TH1F* laben_1cl;
    TH1F* lg_1cl;

    TH1F* bad_crate;
    TH1F* bad_feb;
    TH1F* bad_laben;
    TH1F* bad_lg;
    
       //phi - cos(theta) & company
    TH2F* theta_vs_phi_PMT;
    TH2F* theta_vs_phi_PMT_ev_alive;
    TH2F* theta_vs_phi_PMT_ev_tot;
    TH2F* theta_vs_phi_all;
    TH2F* theta_vs_phi;
    TH1F* virt_pmt_rel_var;

    TH1F* ns_asym;

    TH2F* plane_cos_chi2;  
    TH1F* h_plane_chi2;  

    TH1F* sphere_chi2;
    TH1F* sphere_lkl;

    TH1F* num_bad_crates;
    TH1F* num_bad_boards;
    TH1F* aYp0;
    TH1F* aYp1;
    TH1F* aYp2;
    TH1F* aYp3; 
  
    bool *p_disabled_lg;

   int n_crate;
  int n_feb;
  int n_laben;
  int nlg;
  
};
#endif
/* 
 * $Log: bx_pid_shape.hh,v $
 * Revision 1.10  2009/10/08 17:35:23  ludhova
 * new variables
 *
 * Revision 1.9.6.1  2009-09-28 11:56:02  ludhova
 * distribution in cosTheta-Phi, charge and position corrections
 *
 * Revision 1.9  2008-02-26 18:31:56  ddangelo
 * added filling of is_pointlike variable based on meantime and peaktime values
 * added paramaters with limits and their readout from config
 * Events without any peak are considered pointlike.
 *
 * Commit done while maintainer was out of her mind.
 *
 * Revision 1.8  2007-10-26 13:21:55  ludhova
 * new variable names
 *
 * Revision 1.7  2007-10-25 10:38:50  ludhova
 * some cosmetics
 *
 * Revision 1.6  2007-10-24 13:53:32  ludhova
 * new search for events with strange electronics distribution
 *
 * Revision 1.5  2006-08-30 15:51:05  ludhova
 * identificattion of bad crates and boards changed
 *
 * Revision 1.4  2005/09/21 14:00:43  razeto
 * Added binning in spheric harmonic calculation
 *
 * Revision 1.3  2005/09/21 12:49:33  razeto
 * Fixed normalization coefficient for spherical harmonic
 *
 * Revision 1.2  2005/09/20 17:19:46  razeto
 * Updated to a just working code
 */
