/* BOREXINO Reconstruction program
 *
 * Author: Alessandro Razeto <Alessandro.Razeto@ge.infn.it>
 * Maintainer: Alessandro Razeto <Alessandro.Razeto@ge.infn.it>
 *
 * $Id: bx_precalib_laben_d80.hh,v 1.7 2004/11/26 15:25:11 razeto Exp $
 *
 * Precalib for laben phases (3st cycle): calcolate d80
 *
 */
#ifndef _BX_PRECALIB_LABEN_D80_H
#define _BX_PRECALIB_LABEN_D80_H
#include "bx_base_module.hh"
#include "bx_echidna_event.hh"
#include "bx_rec_general.hh"


#include "TROOT.h"
#include "TH1F.h"
#include "TH2F.h"


class bx_precalib_laben_d80: public bx_base_module {
  public:
    bx_precalib_laben_d80 ();
    virtual ~bx_precalib_laben_d80 () {}

    virtual void begin ();
    virtual bx_echidna_event* doit (bx_echidna_event *ev);
    virtual void end ();
  private:
    int32_t **d80_delay_map;
    float f_d80_low_limit;
    float f_d80_high_limit;
    TH1F *d80;
    TH2F *d80_channel_distrubution;
};

#endif
/*
 * $Log: bx_precalib_laben_d80.hh,v $
 * Revision 1.7  2004/11/26 15:25:11  razeto
 * Thanks to Davide who so kindly underlined a syntactical error, I fixed the "Maintainer" word
 *
 * Revision 1.6  2004/11/26 14:01:21  razeto
 * Added Mantainer field
 *
 * Revision 1.5  2004/09/28 13:43:27  razeto
 * Removed the ifdef for root barn histos
 *
 * Revision 1.4  2004/09/22 13:25:24  razeto
 * Changed bx_reco_event to bx_echidna_event
 *
 * Revision 1.3  2004/05/21 08:39:13  razeto
 * Updated to use bx_root_barn
 *
 * Revision 1.2  2004/05/18 15:01:54  razeto
 * A lot of development done
 *
 * Revision 1.1  2004/04/26 13:50:51  razeto
 * Added
 *
 */
