/* BOREXINO Reconstruction program
 *
 * Author: Davide D'Angelo <davide.dangelo@lngs.infn.it>
 * Maintainer: Davide D'Angelo <davide.dangelo@lngs.infn.it>
 *
 * $Id: bx_precalib_muon_findpulse.hh,v 1.11 2016/07/29 00:49:43 ddangelo Exp $
 *
 * Special module to do outer muon precalibration
 * Cycle 1: find the precalibration signal time (raw)
 *
 */

#ifndef _BX_PRECALIB_MUON_FINDPULSE_HH
#define _BX_PRECALIB_MUON_FINDPULSE_HH

#include "bx_base_module.hh"
#include "bx_echidna_event.hh"
#include "bx_rec_general.hh"
#include "TH1F.h"
#include "TSpectrum.h"

class bx_precalib_muon_findpulse: public bx_base_module {
  public:
    bx_precalib_muon_findpulse ();
    virtual ~bx_precalib_muon_findpulse () {}

    virtual void begin ();
    virtual bx_echidna_event* doit (bx_echidna_event *ev);
    virtual void end ();
  private:
    TH1F* pulses;

};

#endif

/*
 * $Log: bx_precalib_muon_findpulse.hh,v $
 * Revision 1.11  2016/07/29 00:49:43  ddangelo
 * dynamic handling of muon precalib pulse
 *
 * Revision 1.10  2006/09/05 12:06:03  ddangelo
 * cleaned up
 * histograms filled
 *
 * Revision 1.9  2006/08/21 15:34:01  ddangelo
 * test histogram added (new barn interface).
 * introduced the possibility to use mode instead of mean (parameter switch).
 *
 * Revision 1.8  2004/11/26 15:25:11  razeto
 * Thanks to Davide who so kindly underlined a syntactical error, I fixed the "Maintainer" word
 *
 * Revision 1.7  2004/11/26 14:06:23  razeto
 * Added Mantainer field
 *
 * Revision 1.6  2004/09/22 13:29:16  ddangelo
 * updated bx_reco_event to bx_echidna_event
 *
 * Revision 1.5  2004/05/30 13:51:19  ddangelo
 * b_has_data removed since now bx_base_module provides it
 *
 * Revision 1.4  2004/05/20 15:41:53  ddangelo
 * check on has_data introduced
 * check on chann descr implemented with db data (db_profile for the moment)
 * instead of static event members.
 * computetion done with algotithms instead of c like math.
 * Data stored in a std::vector* (histo use).
 *
 * Revision 1.3  2004/04/27 15:05:53  ddangelo
 * fixed the bug that was causing the crash with no muon data.
 *
 * Revision 1.2  2004/04/27 09:46:50  ddangelo
 * modifications to match new event structure.
 * Other debugging
 *
 * Revision 1.1  2004/04/20 11:36:25  ddangelo
 * added 2 modules for outer muon precalibration
 *
 *
 */
