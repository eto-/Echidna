/* BOREXINO Reconstruction program
 *
 * Author: Maria Elena Monzani <monzani@mi.infn.it>
 * Maintainer: Maria Elena Monzani <monzani@mi.infn.it>
 *
 * $Id: bx_calib_laben_charge_peak.hh,v 1.5 2008/11/19 15:50:40 ludhova Exp $
 *
 * Get the time allignement for crates
 *
 */
#ifndef _BX_CALIB_LABEN_CHARGE_PEAK_H
#define _BX_CALIB_LABEN_CHARGE_PEAK_H

#include "bx_base_module.hh"
#include "bx_echidna_event.hh"
#include "bx_rec_general.hh"
#include "TH2S.h"
#include "TF1.h"

class TH2S;
class TF1;
class TH1F;

class bx_calib_laben_charge_peak: public bx_base_module {
  public:
    bx_calib_laben_charge_peak ();
    virtual ~bx_calib_laben_charge_peak () {}

    virtual void begin ();
    virtual bx_echidna_event* doit (bx_echidna_event *ev);
    virtual void end ();
  private:
    TH2S* bx_adc_charge_calib, *bx_q_charge_calib;
    TH1F  *p0_vs_lg, *mu_vs_lg;
    TH1F  *fit1_vs_lg, *mean_vs_lg, *fit2_vs_lg, *rms_vs_lg;
    TH1F  *h_mean, *h_rms, *h_fit1, *h_fit2;
    int32_t nTriggers;
    float* f_charge_peak;
    float* f_charge_sigma;
};

#endif
/*
 * $Log: bx_calib_laben_charge_peak.hh,v $
 * Revision 1.5  2008/11/19 15:50:40  ludhova
 * some new histos + new binning
 *
 * Revision 1.4  2007-11-13 14:12:51  smirnov
 * Gaussian fit tuned. Default values in the case of the fit failure will be set to the adc mean value
 *
 * Revision 1.3  2007-04-13 13:53:45  razeto
 * Added real charge histogram (for displaying only)
 *
 * Revision 1.2  2005-05-05 17:13:40  monzani
 * S/Getters names updated according to the interface modification. Debugging.
 *
 * Revision 1.1  2005/05/05 10:03:48  monzani
 * Laser charge calibrations added (some debug still needed).
 *
 */
