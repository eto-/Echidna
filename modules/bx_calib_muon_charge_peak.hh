/* BOREXINO Reconstruction program
 *
 * Author: Davide D'Angelo <davide.dangelo@mi.infn.it>
 * Based on code from Maria Elena Monzani
 * Maintainer: Michael Wurm <mwurm@ph.tum.de>
 *
 * $Id: bx_calib_muon_charge_peak.hh,v 1.6 2008/12/15 14:32:01 wurm Exp $
 *
 * Outer Muon calibration.
 * Performes channel time alignment and charge peak computation 
 * Works on timing laser runs (LED data) 
 *
 */

#ifndef _BX_CALIB_MUON_CHARGE_PEAK_HH
#define _BX_CALIB_MUON_CHARGE_PEAK_HH

#include "bx_base_module.hh"
#include "bx_echidna_event.hh"
#include "bx_rec_general.hh"
#include "TH2S.h"

class TH2S;

class bx_calib_muon_charge_peak: public bx_base_module {
  public:
    bx_calib_muon_charge_peak ();
    virtual ~bx_calib_muon_charge_peak () {}

    virtual void begin ();
    virtual bx_echidna_event* doit (bx_echidna_event *ev);
    virtual void end ();
  private:
    TH2S* muon_charge_calib;
    int32_t i4_nevents;
    float f4_time_offset;
    float f4_time_allowance;
    float f4_chi2_max;
    float f4_chi2_min;
};

#endif
/*
 * $Log: bx_calib_muon_charge_peak.hh,v $
 * Revision 1.6  2008/12/15 14:32:01  wurm
 * added parameters for chi2
 *
 * Revision 1.5  2008-10-16 19:48:16  ddangelo
 * channel efficiency check improved,
 * time cut
 * some parameters tuning
 *
 * Revision 1.4  2008-10-14 15:25:44  wurm
 * debugging
 *
 * Revision 1.3  2008-10-14 14:32:50  wurm
 * split module off calib_muon_channel
 *
 * Revision 1.1  2006-09-11 14:12:38  ddangelo
 * bx_calib_muon_charge_peak and bx_calib_muon_time_align
 * replaced by
 * bx_calib_muon_channel
 *
 *
 */
