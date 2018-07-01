/* BOREXINO Reconstruction program
 *
 * Author: Alessandro Razeto <Alessandro.Razeto@ge.infn.it>
 * Maintainer: Alessandro Razeto <Alessandro.Razeto@ge.infn.it>
 *
 * $Id: bx_calib_laben_crate_delay.hh,v 1.3 2004/11/26 15:25:11 razeto Exp $
 *
 * Get the time allignement for crates
 *
 */
#ifndef _BX_CALIB_LABEN_CRATE_DELAY_H
#define _BX_CALIB_LABEN_CRATE_DELAY_H

#include "bx_base_module.hh"
#include "bx_echidna_event.hh"
#include "bx_rec_general.hh"
#include "constants.hh"

class bx_calib_laben_crate_delay: public bx_base_module {
  public:
    bx_calib_laben_crate_delay ();
    virtual ~bx_calib_laben_crate_delay () {}

    virtual void begin ();
    virtual bx_echidna_event* doit (bx_echidna_event *ev);
    virtual void end ();
  private:
    double *p_crate_time_sums;
    int32_t i4_event_count;
};

#endif
/*
 * $Log: bx_calib_laben_crate_delay.hh,v $
 * Revision 1.3  2004/11/26 15:25:11  razeto
 * Thanks to Davide who so kindly underlined a syntactical error, I fixed the "Maintainer" word
 *
 * Revision 1.2  2004/11/26 14:01:20  razeto
 * Added Mantainer field
 *
 * Revision 1.1  2004/10/05 13:50:59  razeto
 * Changed name to conform to the new bx_calib_* standard
 *
 * Revision 1.2  2004/09/22 13:27:34  razeto
 * Changed bx_reco_event to bx_echidna_event
 *
 * Revision 1.1  2004/08/10 11:14:03  razeto
 * Added a module for calculating the time differences from crate to crate
 *
 */
