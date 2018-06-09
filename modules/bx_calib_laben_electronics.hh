/* BOREXINO Reconstruction program
 *
 * Author: Daniela Manuzio <dmanuzio@ge.infn.it> and Maria Elena Monzani <monzani@mi.infn.it>
 * Maintainer: Livia Ludhova <livia.ludhova@mi.infn.it>
 *
 * $Id: bx_calib_laben_electronics.hh,v 1.10 2011/02/25 19:02:33 ludhova Exp $
 *
 * Definition of bx_calib_laben_electronics
 * Module that studies the status of the inner detector electronics channels
 * In detail: 
 * checking channels with too many FIFO full and empty raw hits
 * gives mean + error for base, peak, charge for the whole detector using pulser triggers
 * list and number of dead channels (pulser, laser, neutrino triggers separately, threshold from echidna.cfg)
 * list and number of low efficient channels  (pulser, laser, neutrino triggers separately, threshold from echidna.cfg)
 * list and number of hot channels  (pulser, laser, neutrino triggers separately, threshold from echidna.cfg)
 * list and number of retriggering channels  == second hits 0-200 ns after the 1st hit (pulser, laser, neutrino triggers separately, threshold from echidna.cfg)
 * list and number of negative-charge-hits  channels  (threshold from echidna.cfg)
 * identifies channels with too many 0-values for base, peak, charge
 * identifies channels with too many FF-values for base and peak
 * identifies channels with strange distribution of base, peak, charge:
 *    too much spread values
 *    shifted values from mean
 *    too big (rms_ADC threshold from echidna.cfg) or too thin rms (less than 1 ADC channel)
 * list and number of channels correlated with
 *     trigger reference (pulser, laser, neutrino triggers separately)
 *     end of the gate (pulser, laser, neutrino triggers separately)
 *     laser reference (laser triggers separately)
 *
 * Input: the decoded hits
 * Output: some histos and a list and number of bad channels
 * 
*/

#ifndef _BX_CALIB_LABEN_ELECTRONICS_HH
#define _BX_CALIB_LABEN_ELECTRONICS_HH
#include "bx_rec_general.hh"
#include "bx_base_module.hh"
#include "TROOT.h"
#include "TH1F.h"
#include "TH2F.h"
#include <map>

class bx_echidna_event;

class bx_calib_laben_electronics : public bx_base_module {
  public:
    bx_calib_laben_electronics ();
    virtual ~bx_calib_laben_electronics () {}
      
    // Operations for the framework 
    virtual void begin ();
    virtual bx_echidna_event* doit (bx_echidna_event *ev);
    virtual void end ();

  private:
    int  Nevents_type[3];
    int* hits[3];
    int* raw_hits[3];
      
    float f_dead;
    float f_low_eff[3];
    float f_hot[3];
    float f_retriggering;
    
    float f_zero;
    float f_0xFF;
    float f_too_spread;
    float f_mean_offset;
    float f_rms_ADC, f_rms_charge;
    float f_negative_charge;
      
    TH1F *fifo_full_vs_lg;
    TH1F *fifo_empty_vs_lg;

    TH2F *base_vs_lg;
    TH2F *peak_vs_lg;
    TH2F *charge_vs_lg;

    TH2F *charge_tt1_vs_lg;

    TH2F *pulser_charge_vs_lg;
    TH2F *laser_charge_vs_lg;
  
    TH1F *trigref;
    TH1F *lasref;
    
    TH2F *lasref_correl_vs_lg; 
    TH2F *tmp_vs_lg;
    TH2F *trigref_correl_vs_lg[3];
    TH2F *end_gate_correl_vs_lg[3];
        
    
    TH1F *raw_Nhits_vs_lg[3];	
    TH1F *Nhits_vs_lg[3];	
    TH1F *frac_Nhits_vs_lg[3];	
    
    TH2F *retrigger_dt[3];
   
  enum  ADC_status {
    many_zero            = 1,
    many_FF              = 2,
    too_spread           = 3,
    shifted_from_mean    = 4,
    bad_rms              = 5,
    very_small_rms       = 6,
    many_negative_values = 7,
    low_gain             = 8,
    high_gain            = 9,
    };
  
  enum timing_status {
    laser_ref_correl                = 11,
    trigger_ref_correl_in_pulser    = 20,
    trigger_ref_correl_in_laser     = 21,
    trigger_ref_correl_in_neutrino  = 22,
    end_of_gate_correl_in_pulser    = 30,
    end_of_gate_correl_in_laser     = 31,
    end_of_gate_correl_in_neutrino  = 32,
    bad_timing_shape_in_laser       = 41,
  };
  
  enum multiplicity {
    dead_in_pulser           = 10,
    dead_in_laser            = 11,
    dead_in_neutrino         = 12,
    dead_in_raw              = 13,
    low_eff_in_pulser        = 20,
    low_eff_in_laser         = 21,
    low_eff_in_neutrino      = 22,
    low_eff_in_raw           = 23,
    hot_in_pulser            = 30,
    hot_in_laser             = 31,
    hot_in_neutrino          = 32,
    hot_in_raw               = 33,
    retriggering_in_pulser   = 40,
    retriggering_in_laser    = 41,
    retriggering_in_neutrino = 42,
    fifo_empty               = 5012,
    fifo_full                = 6012,
    loosing_raw_hits_in_pulser   = 50,
    loosing_raw_hits_in_laser    = 51,
    loosing_raw_hits_in_neutrino = 52,  
};
  
  std::cmap <std::string, ADC_status> ADC_translation_map;
  std::cmap <std::string, timing_status> timing_translation_map;
  std::cmap <std::string, multiplicity> multiplicity_translation_map;

  enum trg_type{
    pulser = 0,
    laser = 1,
    neutrino = 2,
  };


  std::cmap <int, std::string> trg_names;


  float fit_line(TH1F* h, int level_guess);

};

#endif
/*
 * $Log: bx_calib_laben_electronics.hh,v $
 * Revision 1.10  2011/02/25 19:02:33  ludhova
 * low and high gain during el calib
 *
 * Revision 1.9  2007-11-12 14:08:28  ludhova
 * search for dead, low_eff, hot in raw hits in neutrino triggers - write it in DB
 *
 * Revision 1.8  2007-06-18 14:13:56  ludhova
 * histos for ref channels and lg of raw hits added
 *
 * Revision 1.7  2006-11-07 15:56:17  ludhova
 * control of the ratio N_raw_hits/N_decoded_hits added
 *
 * Revision 1.6  2006-10-30 11:47:35  ludhova
 * fifo full and fifo empty controll added
 *
 * Revision 1.5  2006-09-19 14:26:33  ludhova
 * substantial changes in the electronics channels quality controll
 *
 * Revision 1.4  2006/03/22 13:22:43  ludhova
 * some new diagnostics for electronic channels added
 *
 * Revision 1.3  2004/11/26 15:25:11  razeto
 * Thanks to Davide who so kindly underlined a syntactical error, I fixed the "Maintainer" word
 *
 * Revision 1.2  2004/11/26 14:16:54  razeto
 * Added Mantainer field
 *
 * Revision 1.1  2004/10/05 13:33:28  dmanuzio
 * Removed a bug and changed the name of 2 modules (for consistency with
 * bx_calib_fiber_bundle module)
 *
 * Revision 1.3  2004/09/28 13:31:13  dmanuzio
 * Removed the ifdef for root barn histos
 *
 * Revision 1.2  2004/09/22 14:20:47  dmanuzio
 * Changed bx_reco_event to bx_echidna_event
 *
 * Revision 1.1  2004/07/22 09:57:50  dmanuzio
 * New module to test the electronics (laben and front end problems)
 *
 *
 */
