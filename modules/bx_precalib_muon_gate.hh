/* BOREXINO Reconstruction program
 *
 * Author: Davide D'Angelo <davide.dangelo@mi.infn.it>
 * Maintainer: Davide D'Angelo <davide.dangelo@mi.infn.it>
 *
 * $Id: bx_precalib_muon_gate.hh,v 1.1 2008/10/25 09:59:54 ddangelo Exp $
 *
 * Special module to do outer muon precalibration
 * cycle 2: compute QTC gate
 */

#ifndef _BX_PRECALIB_MUON_GATE_HH
#define _BX_PRECALIB_MUON_GATE_HH

#include "bx_base_module.hh"
#include "bx_echidna_event.hh"
#include "bx_rec_general.hh"
#include "TH1F.h"
#include "TH2F.h"

class bx_precalib_muon_gate: public bx_base_module {
  public:
    bx_precalib_muon_gate ();
    virtual ~bx_precalib_muon_gate () {}

    virtual void begin ();
    virtual bx_echidna_event* doit (bx_echidna_event *ev);
    virtual void end ();
  private:
    unsigned short min_evnum;
    unsigned short max_evnum;
    unsigned short min_gate_width; 
    unsigned short max_gate_width; 
    unsigned short precalib_pulse_tolerance; 
    unsigned short pulse_time;

    TH2S* gates;
    TH1F* means;
    TH1F* sigmas;
};

#endif

/*
 * $Log: bx_precalib_muon_gate.hh,v $
 * Revision 1.1  2008/10/25 09:59:54  ddangelo
 * added. cloned from pedestal module
 *
 *
 */
