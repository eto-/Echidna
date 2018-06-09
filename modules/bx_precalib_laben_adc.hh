/* BOREXINO Reconstruction program
 *
 * Author: Alessandro Razeto <Alessandro.Razeto@ge.infn.it>
 * Maintainer: Alessandro Razeto <Alessandro.Razeto@ge.infn.it>
 *
 * $Id: bx_precalib_laben_adc.hh,v 1.7 2004/11/26 15:25:11 razeto Exp $
 *
 * Precalib for laben phases (1st cycle): calcolate the ramp limits and update them to
 * db_run. 
 *
 */
#ifndef _BX_PRECALIB_LABEN_ADC_H
#define _BX_PRECALIB_LABEN_ADC_H
#include "bx_base_module.hh"
#include "bx_echidna_event.hh"
#include "bx_rec_general.hh"

#include "TROOT.h"
#include "TH1F.h"

class bx_precalib_laben_adc: public bx_base_module {
  public:
    bx_precalib_laben_adc ();
    virtual ~bx_precalib_laben_adc () {}

    virtual void begin ();
    virtual bx_echidna_event* doit (bx_echidna_event *ev);
    virtual void end ();
  private:
    int **adc_sample_map;
    unsigned char u1_maxima_bin_add;
    TH1F *adc_limits, *adc_samples;
};

#endif
/*
 * $Log: bx_precalib_laben_adc.hh,v $
 * Revision 1.7  2004/11/26 15:25:11  razeto
 * Thanks to Davide who so kindly underlined a syntactical error, I fixed the "Maintainer" word
 *
 * Revision 1.6  2004/11/26 14:01:20  razeto
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
