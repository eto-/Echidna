/* BOREXINO Reconstruction program
 *
 * Author: Barbara Caccianiga <Barbara.Caccianiga@mi.infn.it>
 * Maintainer: Barbara Caccianiga <Barbara.Caccianiga@mi.infn.it>
 *
 * $Id: bx_calib_fiber_bundle.hh,v 1.3 2004/11/26 15:25:11 razeto Exp $
 *
 * Calculates some useful parameters of fiber bundles (luminosity, dead fibers...)
 * 
 */
#ifndef _BX_CALIB_FIBER_BUNDLE_H
#define _BX_CALIB_FIBER_BUNDLE_H

#include "bx_base_module.hh"
#include "bx_rec_general.hh"
#include "bx_echidna_event.hh"
#include "TH1F.h"
#include "TH2F.h"

class bx_calib_fiber_bundle: public bx_base_module {
  public:
    bx_calib_fiber_bundle ();
    virtual ~bx_calib_fiber_bundle () {}

    virtual void begin ();
    virtual bx_echidna_event* doit (bx_echidna_event *ev);
    virtual void end ();
  private:
    int32_t i4_times;
    float f_laser_low;
    float f_laser_high;
    float f_dark_low;
    float f_dark_high;
    int32_t n_bundles;
    TH1F *p_timetot;
    TH2F *p_laser;
    TH2F *p_dark_noise;
};

#endif
/*
 * $Log: bx_calib_fiber_bundle.hh,v $
 * Revision 1.3  2004/11/26 15:25:11  razeto
 * Thanks to Davide who so kindly underlined a syntactical error, I fixed the "Maintainer" word
 *
 * Revision 1.2  2004/11/26 14:16:54  razeto
 * Added Mantainer field
 *
 * Revision 1.1  2004/10/01 10:10:16  bcaccian
 * Added
 *
 */
