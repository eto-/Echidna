/* BOREXINO Reconstruction program
 *
 * Author: Davide D'Angelo <davide.dangelo@mi.infn.it>, Michael Wurm <mwurm@ph.tum.de>, Quirin Meindl <qmeindl@ph.tum.de>
 * Based on code from Maria Elena Monzani
 * Maintainer: Michael Wurm <mwurm@ph.tum.de>
 *
 * $Id: bx_calib_muon_time_align.hh,v 1.5 2008/10/16 19:48:16 ddangelo Exp $
 *
 * Outer Muon calibration.
 * Performes channel time alignment 
 * Works on timing laser runs (LED data) 
 *
 */

#ifndef _BX_CALIB_MUON_TIME_ALIGN_HH
#define _BX_CALIB_MUON_TIME_ALIGN_HH

#include "bx_base_module.hh"
#include "bx_echidna_event.hh"
#include "bx_rec_general.hh"
#include "TH2S.h"

class TH2S;

class bx_calib_muon_time_align: public bx_base_module {
  public:
    bx_calib_muon_time_align ();
    virtual ~bx_calib_muon_time_align () {}

    virtual void begin ();
    virtual bx_echidna_event* doit (bx_echidna_event *ev);
    virtual void end ();
  private:
    int32_t i4_nevents;
    TH2S* muon_time_calib;
};

#endif
