/* BOREXINO Reconstruction program
 *
 * Author: Alessandro Razeto <Alessandro.Razeto@ge.infn.it>
 * 	   Marco Pallavicini <pallas@ge.infn.it>
 * Maintainer: Davide.Dangelo <Davide.Dangelo@lngs.infn.it>
 *
 * $Id: bx_cmt_tracker.hh,v 1.2 2010/05/21 14:53:37 bick Exp $
 *
 * Sample module to explain echidna module writing
 * to new developers.
 * Please read the "Module programmer's guide" before 
 * or while looking at the examples here 
 * 
 */

#ifndef _BX_CMT_TRACKER_HH
#define _BX_CMT_TRACKER_HH

#include "bx_base_module.hh"
#include "bx_rec_general.hh"

class bx_laben_event;
class TH1F;
class TH2F;

class bx_cmt_tracker: public bx_base_module {
  public:
  // this section if fixed; do not edit
    bx_cmt_tracker ();
    virtual ~bx_cmt_tracker () {}

    virtual void begin ();
    virtual bx_echidna_event* doit (bx_echidna_event *ev);
    virtual void end ();
  private:
  // this section is free; add here members, methods, histo pointers
    int i4_times;
    std::vector<double> *my_vector;
    bool m_check_this_and_that(const bx_laben_event& er);
    TH1F* my_histo_check;
    TH2F* my_histo_compute;
};

#endif
/*
 * $Log: bx_cmt_tracker.hh,v $
 * Revision 1.2  2010/05/21 14:53:37  bick
 * Just getting started
 *
 * Revision 1.1  2010-05-21 12:34:01  ddangelo
 * bx_laben_tracker renamed as bx_laben_energy_tracker
 * added bx_laben_tof_tracker and bx_cmt_tracker
 *
 * Revision 1.6  2005-02-03 19:00:18  ddangelo
 * maintainer changed (Razeto->D'Angelo)
 * Module really implemented with many possible examples for developers.
 * It features all aspects described in "Module's programmers guide" in the docs.
 *
 * Revision 1.5  2004/11/26 15:25:11  razeto
 * Thanks to Davide who so kindly underlined a syntactical error, I fixed the "Maintainer" word
 *
 * Revision 1.4  2004/11/26 14:01:21  razeto
 * Added Mantainer field
 *
 * Revision 1.3  2004/09/22 13:28:37  razeto
 * Changed bx_reco_event to bx_echidna_event
 *
 * Revision 1.2  2004/03/21 18:55:05  razeto
 * Some cosmetic changes
 *
 * Revision 1.1  2004/03/20 17:39:45  pallas
 * First release
 *
 * Revision 1.1  2004/03/20 14:16:47  razeto
 * Added sources
 *
 */
