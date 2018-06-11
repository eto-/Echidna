/* BOREXINO Reconstruction program
 *
 * Author: Alessandro Razeto <Alessandro.Razeto@ge.infn.it>
 * Maintainer: Alessandro Razeto <Alessandro.Razeto@ge.infn.it>
 *
 * $Id: bx_laben_raw_validator.hh,v 1.5 2010/07/01 18:21:35 razeto Exp $
 *
 * Validate raw laben events, checking for patterns
 * 
 */

#ifndef _BX_LABEN_RAW_VALIDATOR_HH
#define _BX_LABEN_RAW_VALIDATOR_HH

#include "bx_base_module.hh"

class TH1F;

class bx_laben_raw_validator: public bx_base_module {
  public:
    bx_laben_raw_validator ();
    virtual ~bx_laben_raw_validator () {}

    virtual void begin ();
    virtual bx_echidna_event* doit (bx_echidna_event *ev);
    virtual void end ();
  private:
    TH1F *flag_lg[8], *good_flag_lg_c; // 8 are the possible flags (see bx_laben_raw_hits::flags)
    char *board_occupancy;
    int32_t count;
};

#endif
/*
 * $Log: bx_laben_raw_validator.hh,v $
 * Revision 1.5  2010/07/01 18:21:35  razeto
 * Added counter hit flag for new laben fw and nhits_fw from laben boards
 *
 * Revision 1.4  2009-01-30 11:14:34  razeto
 * Added non cumulative plot
 *
 * Revision 1.3  2008-10-07 14:03:56  razeto
 * Using new flag
 *
 * Revision 1.2  2007-10-30 17:35:14  razeto
 * Reduced verbosity, added board occupancy
 *
 * Revision 1.1  2006-10-23 14:41:16  razeto
 * New module bx_laben_raw_validator to validate event using laben raw hit flags
 *
 */
