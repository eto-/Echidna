/* BOREXINO Reconstruction program
 *
 * Author: Davide D'Angelo <davide.dangelo@lngs.infn.it>
 * Maintainer: Davide D'Angelo <davide.dangelo@lngs.infn.it>
 *
 * $Id: bx_precalib_muon_pedestals.hh,v 1.10 2008/10/25 09:59:32 ddangelo Exp $
 *
 * Special module to do outer muon precalibration
 * cycle 2: compute QTC pedestals and their spreads
 */

#ifndef _BX_PRECALIB_MUON_PEDESTALS_HH
#define _BX_PRECALIB_MUON_PEDESTALS_HH

#include "bx_base_module.hh"
#include "bx_echidna_event.hh"
#include "bx_rec_general.hh"
#include "TH1F.h"
#include "TH2F.h"

class bx_precalib_muon_pedestals: public bx_base_module {
  public:
    bx_precalib_muon_pedestals ();
    virtual ~bx_precalib_muon_pedestals () {}

    virtual void begin ();
    virtual bx_echidna_event* doit (bx_echidna_event *ev);
    virtual void end ();
  private:
    unsigned short min_evnum;
    unsigned short max_evnum;
    unsigned short min_pedestal_width; 
    unsigned short max_pedestal_width; 
    unsigned short precalib_pulse_tolerance; 
    unsigned short pulse_time;

    TH2S* pedestals;
    TH1F* means;
    TH1F* sigmas;
};

#endif

/*
 * $Log: bx_precalib_muon_pedestals.hh,v $
 * Revision 1.10  2008/10/25 09:59:32  ddangelo
 * added check on event number
 *
 * Revision 1.9  2006-09-05 12:06:03  ddangelo
 * cleaned up
 * histograms filled
 *
 * Revision 1.8  2006/08/21 15:37:47  ddangelo
 * added a few test histograms (new barn interface).
 * introduced the possibility to use a fixed pulse time (parameter switch and constant setting).
 * Retrieving of find pulse results moved in begin().
 *
 * Revision 1.7  2004/11/26 15:25:11  razeto
 * Thanks to Davide who so kindly underlined a syntactical error, I fixed the "Maintainer" word
 *
 * Revision 1.6  2004/11/26 14:06:23  razeto
 * Added Mantainer field
 *
 * Revision 1.5  2004/09/22 13:29:16  ddangelo
 * updated bx_reco_event to bx_echidna_event
 *
 * Revision 1.4  2004/05/30 13:51:19  ddangelo
 * b_has_data removed since now bx_base_module provides it
 *
 * Revision 1.3  2004/05/20 15:45:56  ddangelo
 *  check on has_data introduced
 *  check on chann descr implemented with db data (db_profile for the moment)
 *  instead of static event members.
 *  computetion done with algotithms instead of c like math.
 *  Data stored in a std::vector* (histo use).
 *
 * Revision 1.2  2004/04/28 10:48:21  ddangelo
 * Pedestal calculation reworked.
 *
 * Revision 1.1  2004/04/20 11:36:25  ddangelo
 * added 2 modules for outer muon precalibration
 *
 *
 */
