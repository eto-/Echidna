/* BOREXINO Reconstruction program
 *
 * Author: Alessandro Razeto <Alessandro.Razeto@ge.infn.it>
 * Maintainer: Alessandro Razeto <Alessandro.Razeto@ge.infn.it>
 *
 * $Id: bx_precalib_laben_check_tdc.hh,v 1.8 2005/03/17 16:07:41 razeto Exp $
 *
 * Precalib for laben TDC (4st cycle): check the time allignement
 *
 */
#ifndef _BX_PRECALIB_LABEN_CHECK_TDC_H
#define _BX_PRECALIB_LABEN_CHECK_TDC_H
#include "bx_base_module.hh"
#include "bx_echidna_event.hh"
#include "bx_rec_general.hh"
#include "cmap.hh"

#include "TROOT.h"
#include "TH2F.h"


class bx_precalib_laben_check_tdc: public bx_base_module {
  public:
    bx_precalib_laben_check_tdc ();
    virtual ~bx_precalib_laben_check_tdc () {}

    virtual void begin ();
    virtual bx_echidna_event* doit (bx_echidna_event *ev);
    virtual void end ();
  private:
    int32_t i_ref_channel;
    bool b_found_mctruth;
    TH2F *time_channel_distribution;
    TH2F *error_channel_distribution;
    TH2F *time_bits;
};

#endif
/*
 * $Log: bx_precalib_laben_check_tdc.hh,v $
 * Revision 1.8  2005/03/17 16:07:41  razeto
 * Added check for mctruth presence before writing precalib data
 *
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
 * Revision 1.3  2004/08/10 11:12:42  razeto
 * Added a parameter
 *
 * Revision 1.2  2004/05/21 08:39:13  razeto
 * Updated to use bx_root_barn
 *
 * Revision 1.1  2004/05/18 15:01:54  razeto
 * A lot of development done
 *
 */
