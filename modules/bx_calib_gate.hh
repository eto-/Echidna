/* BOREXINO Reconstruction program
 *
 * Author: Barbara Caccianiga <Barbara.Caccianiga@mi.infn.it>
 * Maintainer: Barbara Caccianiga <Barbara.Caccianiga@mi.infn.it>
 *
 * 
 */
#ifndef _BX_CALIB_GATE_H
#define _BX_CALIB_GATE_H

#include "bx_base_module.hh"
#include "bx_rec_general.hh"
#include "bx_echidna_event.hh"
#include "TH1F.h"
#include "TH2F.h"

class bx_calib_gate: public bx_base_module {
  public:
    bx_calib_gate ();
    virtual ~bx_calib_gate () {}

    virtual void begin ();
    virtual bx_echidna_event* doit (bx_echidna_event *ev);
    virtual void end ();
  private:
    int i4_times;
    int i4_bin;
    double up_limit;
    double low_limit;
    double media[1000];
    TH1F *p_tempo_random;
    TH1F *p_tempo_laser;
    TH1F *p_tempo_pulser;
    TH1F *p_tempo_neutrino;
    TH1F *p_tempo_random_only_gate;
    TH2F *p_random_ch;
};

#endif
/*
 */
