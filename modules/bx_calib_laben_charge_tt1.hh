/* BOREXINO Reconstruction program
 * 
* Author: Livia Ludhova and Gemma Testera
 *
 *PMT charge calibration using tt1
 */
#ifndef _BX_CALIB_LABEN_CHARGE_TT1_H
#define _BX_CALIB_LABEN_CHARGE_TT1_H

#include "bx_base_module.hh"
#include "bx_echidna_event.hh"
#include "bx_rec_general.hh"
#include "TH2S.h"
#include "TF1.h"

class TH2S;
class TF1;

class bx_calib_laben_charge_tt1: public bx_base_module {
  public:
    bx_calib_laben_charge_tt1 ();
    virtual ~bx_calib_laben_charge_tt1 () {}

    virtual void begin ();
    virtual bx_echidna_event* doit (bx_echidna_event *ev);
    virtual void end ();
  private:
    TH2S* bx_adc_charge_tt1;
    TH1F  *Amp_vs_lg, *C1_vs_lg, *Sig1_vs_lg, *R12_vs_lg, *R13_vs_lg, *Chi2_vs_lg;
    TH1F  *Mean_vs_lg, *Rms_vs_lg, *P0_vs_lg;
    TH1F  *hC1, *hSig1, *hR12, *hR13, *hChi2, *hP0;
    TH1F  *hMean, *hRms;

    int32_t current_run;
    int32_t nTriggers;
    float *f_charge_peak, *f_charge_sigma;
    float *f_charge_mean, *f_charge_rms, *f_p0;
};

#endif
