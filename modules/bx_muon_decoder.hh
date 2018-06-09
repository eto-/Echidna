/* BOREXINO Reconstruction program
 *
 * Author: Davide D'Angelo <davide.dangelo@lngs.infn.it>
 * Maintainer: Davide D'Angelo <davide.dangelo@lngs.infn.it>
 *
 * $Id: bx_muon_decoder.hh,v 1.8 2008/10/09 09:26:20 ddangelo Exp $
 *
 * Modules for decoding of muon infos
 * 
*/

#ifndef _BX_MUON_DECODER_HH
#define _BX_MUON_DECODER_HH

#include "bx_rec_general.hh"
#include "bx_base_module.hh"

class bx_muon_decoder : public bx_base_module {
  public:
    bx_muon_decoder ();
    virtual ~bx_muon_decoder () {}
      
    // Operations for the framework 
    virtual void begin ();
    virtual bx_echidna_event* doit (bx_echidna_event *ev);
    virtual void end ();
  private:
    bool b_discard_disabled_lg, b_apply_calibration;
    bool* p_disabled_mch;
    unsigned int u4_n_disabled_channels;
};

#endif
/*
 * $Log: bx_muon_decoder.hh,v $
 * Revision 1.8  2008/10/09 09:26:20  ddangelo
 * apply calibration made after parameter
 *
 * Revision 1.7  2008-08-20 16:20:23  ddangelo
 * disabling of channels introduced
 *
 * Revision 1.6  2007-03-28 17:43:15  ddangelo
 * code restyled.
 * dec charge/npmts added.
 *
 * Revision 1.5  2004-11-26 15:25:11  razeto
 * Thanks to Davide who so kindly underlined a syntactical error, I fixed the "Maintainer" word
 *
 * Revision 1.4  2004/11/26 14:06:23  razeto
 * Added Mantainer field
 *
 * Revision 1.3  2004/09/22 13:29:16  ddangelo
 * updated bx_reco_event to bx_echidna_event
 *
 * Revision 1.2  2004/05/20 11:02:30  ddangelo
 * trigger type selection introduced.
 * Job splitted in 2 private methods (physics/led)
 * 2nd one finds led ref pulse with a pre-loop.
 * chann descr asked to db_profile (tmp)
 *
 * Revision 1.1  2004/04/02 09:56:03  ddangelo
 * created bx_muon_decoder module.
 * precalibration and calibration/db info are still commented out
 *
 *
 */
