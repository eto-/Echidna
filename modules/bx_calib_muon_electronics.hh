/* BOREXINO Reconstruction program
 *
 * Author: Davide D'Angelo <davide.dangelo@mi.infn.it> (starting from bx_calib_muon_ectronics)
 * Maintainer: Davide D'Angelo <davide.dangelo@mi.infn.it> 
 * 
 * $Id: bx_calib_muon_electronics.hh,v 1.5 2008/08/19 17:53:48 ddangelo Exp $
 *
 * Definition of bx_calib_muon_electronics
 * Module that studies the status of the inner detector electronics channels
 * In detail: 
 * list and number of dead channels (pulser, laser, neutrino triggers separately, threshold from echidna.cfg)
 * list and number of low efficient channels  (pulser, laser, neutrino triggers separately, threshold from echidna.cfg)
 * list and number of hot channels  (pulser, laser, neutrino triggers separately, threshold from echidna.cfg)
 * list and number of retriggering channels  == second hits 0-200 ns after the 1st hit (pulser, laser, neutrino triggers separately, threshold from echidna.cfg)
 *
 * Input: the decoded hits
 * Output: some histos and a list and number of bad channels
 * It writes to DB, bx_calib, MuonChannelProperties
 * 
 */

#ifndef _BX_CALIB_MUON_ELECTRONICS_HH
#define _BX_CALIB_MUON_ELECTRONICS_HH
#include "bx_rec_general.hh"
#include "bx_base_module.hh"
#include "TROOT.h"
#include "TH1F.h"
#include "TH2F.h"
#include <map>

class bx_echidna_event;

class bx_calib_muon_electronics : public bx_base_module {
  public:
    bx_calib_muon_electronics ();
    virtual ~bx_calib_muon_electronics () {}
      
    // Operations for the framework 
    virtual void begin ();
    virtual bx_echidna_event* doit (bx_echidna_event *ev);
    virtual void end ();

  private:
    int  nevents_per_tt[4];
    int* nhits_per_tt_mch[4];
      
    float f_dead;
    float f_low_eff[4];
    float f_hot[4];
    float f_retriggering;
    
    TH1F *h_nhits_vs_mch[4];	
    TH2F *h_dt_vs_mch   [4];
   
  enum multiplicity {
    dead_in_pulser           = 10,
    dead_in_laser            = 11,
    dead_in_neutrino         = 12,
    dead_in_muon             = 13,
    low_eff_in_pulser        = 20,
    low_eff_in_laser         = 21,
    low_eff_in_neutrino      = 22,
    low_eff_in_muon          = 23,
    hot_in_pulser            = 30,
    hot_in_laser             = 31,
    hot_in_neutrino          = 32,
    hot_in_muon              = 33,
    retriggering_in_pulser   = 40,
    retriggering_in_laser    = 41,
    retriggering_in_neutrino = 42,
    retriggering_in_muon     = 43,
  };
  
  std::cmap <std::string, multiplicity> multiplicity_translation_map;

  enum trg_type{
    pulser   = 0,
    laser    = 1,
    neutrino = 2,
    muon     = 3,
  };

  std::cmap <int, std::string> trg_names;
};

#endif
/*
 * $Log: bx_calib_muon_electronics.hh,v $
 * Revision 1.5  2008/08/19 17:53:48  ddangelo
 * debugging
 * some parameter tuning
 *
 * Revision 1.4  2008-08-13 12:25:48  ddangelo
 * added "muon" tt checks
 *
 * Revision 1.3  2008-08-12 17:13:53  ddangelo
 * starting to converge
 *
 * Revision 1.2  2008-08-11 16:36:43  ddangelo
 * it compiles.
 *
 * Revision 1.1  2008-08-11 12:47:52  ddangelo
 * added, still junk
 *
 */
