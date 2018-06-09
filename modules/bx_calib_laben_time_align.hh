/* BOREXINO Reconstruction program
 *
 * Author: Maria Elena Monzani <monzani@mi.infn.it>
 * Maintainer: Maria Elena Monzani <monzani@mi.infn.it>
 *
 * $Id: bx_calib_laben_time_align.hh,v 1.2 2005/05/05 17:13:40 monzani Exp $
 *
 * Get the time allignement for crates
 *
 */
#ifndef _BX_CALIB_LABEN_TIME_ALIGN_H
#define _BX_CALIB_LABEN_TIME_ALIGN_H

#include "bx_base_module.hh"
#include "bx_echidna_event.hh"
#include "bx_rec_general.hh"
#include "TH2S.h"
#include "TF1.h"

class TH2S;
class TF1;

class bx_calib_laben_time_align: public bx_base_module {
  public:
    bx_calib_laben_time_align ();
    virtual ~bx_calib_laben_time_align () {}

    virtual void begin ();
    virtual bx_echidna_event* doit (bx_echidna_event *ev);
    virtual void end ();
  private:
    TH2S* bx_time_calib;
    float* f_time_offset;
    float* f_time_sigma;
};

#endif
/*
 * $Log: bx_calib_laben_time_align.hh,v $
 * Revision 1.2  2005/05/05 17:13:40  monzani
 * S/Getters names updated according to the interface modification. Debugging.
 *
 * Revision 1.1  2005/05/05 10:02:22  monzani
 * Laser time calibrations added.
 *
 */
