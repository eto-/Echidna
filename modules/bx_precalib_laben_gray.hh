/* BOREXINO Reconstruction program
 *
 * Author: Alessandro Razeto <Alessandro.Razeto@ge.infn.it>
 * Maintainer: Alessandro Razeto <Alessandro.Razeto@ge.infn.it>
 *
 * $Id: bx_precalib_laben_gray.hh,v 1.8 2007/01/29 19:03:04 razeto Exp $
 *
 * Precalib for laben gray counter (1st cycle): allign gray counter shift 
 *
 */
#ifndef _BX_PRECALIB_LABEN_GRAY_H
#define _BX_PRECALIB_LABEN_GRAY_H
#include "bx_base_module.hh"
#include "bx_echidna_event.hh"
#include "bx_rec_general.hh"
#include "cmap.hh"

#include "TROOT.h"
#include "TH2F.h"
#include "TH1D.h"


class bx_precalib_laben_gray: public bx_base_module {
  public:
    bx_precalib_laben_gray ();
    virtual ~bx_precalib_laben_gray () {}

    virtual void begin ();
    virtual bx_echidna_event* doit (bx_echidna_event *ev);
    virtual void end ();
  private:
    int i_channel_skip;
    float f_mean_bound;
    float f_rms_bound;
    int i_nreferences;
    int *gray_map;

    TH2F **gray_diffs;
    TH1S *gray_shifts;
	
};

#endif
/*
 * $Log: bx_precalib_laben_gray.hh,v $
 * Revision 1.8  2007/01/29 19:03:04  razeto
 * Use root histogram to calulate mean and rms. This will produce shift in barn
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
 * Revision 1.3  2004/05/21 08:39:13  razeto
 * Updated to use bx_root_barn
 *
 * Revision 1.2  2004/05/18 15:01:54  razeto
 * A lot of development done
 *
 * Revision 1.1  2004/04/12 15:59:49  razeto
 * Added a firsto core of laben precalibration
 *
 */
