/* BOREXINO Reconstruction program
 *
 * Author: Evgeny Litvinovich <litvinov@lngs.infn.it>
 * Maintainer: Evgeny Litvinovich <litvinov@lngs.infn.it>
 *
 * $Id: 
 *
 * Moscow event's energy reconstruction
*/

#ifndef _BX_ENERGY_RECO_MSK_HH
#define _BX_ENERGY_RECO_MSK_HH
#include "bx_base_module.hh"
#include "TVector3.h"
#include "TH1F.h"

class bx_echidna_event;
class TMinuit;

class bx_energy_reco_msk : public bx_base_module {
  public:
    bx_energy_reco_msk();
    virtual ~bx_energy_reco_msk() {}
  
    // Operations for the framework 
    virtual void begin();
    virtual bx_echidna_event* doit(bx_echidna_event *ev);
    virtual void end();
     double msk_fcn_e(double *x);

  private:
     double f8_attenuation_length;
     double f8_coef_to_return_npe_cone;
     double f8_coef_to_return_npe_nocone;

      float f4_ref_index;
      float f4_cathode_efficiency;
      float f4_lg_entry_radius;
      float f4_cathode_radius;
    
     std::vector<double> f8_omega;
     std::vector<double> f8_attenuation_factor;

     std::vector<TVector3> *f4_pmt_positions;
     std::vector<bool> *b_pmt_cone;
     std::vector<float> f4_collected_charge;
     std::vector<int> i4_disabled_channels;

     TMinuit *p_minuit;
     bx_echidna_event *p_fit_ev;
     int i4_fit_cluster;
     int i4_n_minuit_warns;

     TH1F *msk_nphotons;
};

#endif

/*
 * $Log: bx_energy_reco_msk.hh,v $
 * Revision 1.3  2008/10/17 10:09:02  litvinov
 * no more msk_charge, only msk_nphotons
 *
 * Revision 1.2  2007-11-07 19:32:15  litvinov
 * using "dubna" reconstructed position; using "charge" as input instead of "npe"
 *
 */
