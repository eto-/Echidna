/* BOREXINO Reconstruction program
 *
 * Author: Marco Pallavicini <pallas@ge.infn.it>
 * Maintainer: Marco Pallavicini <pallas@ge.infn.it>
 *
 * $Id: bx_trigger_decoder.hh,v 1.11 2011/02/22 16:24:11 razeto Exp $
 *
 * Modules for decoding of trigger infos
 * 
*/

#ifndef _BX_TRIGGER_DECODER_HH
#define _BX_TRIGGER_DECODER_HH

#include "bx_rec_general.hh"
#include "bx_base_module.hh"
#include <time.h>

class bx_echidna_event;
class bx_trigger_event;
class bx_trigger_decoded_event;
class TH1F;

class bx_trigger_decoder : public bx_base_module {
  public:
    bx_trigger_decoder ();
    virtual ~bx_trigger_decoder () {}
      
    // Operations for the framework 
    virtual void begin ();
    virtual bx_echidna_event* doit (bx_echidna_event *ev);
    virtual void end ();

  private:
    bool b_was_a_muon;
    uint32_t u4_prev_trgid;
    int32_t i4_missing_events;
    int32_t i4_missing_neutron_triggers;
    int32_t i4_year, i4_mon, i4_day, i4_hour, i4_min, i4_sec;                                  
    int32_t i4_prev_year,i4_prev_mon,i4_prev_day,i4_prev_hour,i4_prev_min,i4_prev_sec;                                                                            //
    uint32_t prev_gps_times[2];
    
    TH1F *build_dt;
    TH1F *ppc0_dt;

    void m_decode_trgtype(const bx_trigger_event&, bx_trigger_decoded_event&);
    // decode GPS clock infos
    void gpsclock_decode(uint32_t, uint32_t, uint32_t, 
	                 int32_t, 
			 int32_t&, int32_t&, int32_t&, int32_t&, int32_t&, int32_t&, int32_t&, int32_t&, 
			 uint32_t&, uint32_t&, time_t&);
};

#endif
/*
 * $Log: bx_trigger_decoder.hh,v $
 * Revision 1.11  2011/02/22 16:24:11  razeto
 * GPS start from 2000 UTC + added leap seconds + timet
 *
 * Revision 1.10  2009-10-23 13:34:31  ghiano
 * debugging
 *
 * Revision 1.9  2009-10-19 10:30:47  ghiano
 * gps decoding fixed for bad data
 *
 * Revision 1.8  2009-10-15 15:49:55  ghiano
 * gps decoding fixed for bad data
 *
 * Revision 1.7  2009-07-17 14:52:42  ddangelo
 * added histo of dt (ppc0 - gps)
 *
 * Revision 1.6  2008-05-08 19:07:24  razeto
 * Added time checking
 *
 * Revision 1.5  2008-04-14 17:04:20  ddangelo
 * added checks for trgid consecutiveness and for trg128 presence after muon
 *
 * Revision 1.4  2004-11-29 13:21:23  razeto
 * Added Mantainer field
 *
 * Revision 1.3  2004/09/22 14:00:07  ddangelo
 * updated bx_reco_event into bx_echidna_event (on behalf of Marco)
 *
 * Revision 1.2  2004/05/20 12:22:34  ddangelo
 * added a method to decode trgtype.
 * Fixed include problem.
 * commit with delegation from module coordinator
 *
 * Revision 1.1  2004/04/27 10:02:48  ddangelo
 * module bx_trigger_decode_module deleted and re-added as bx_trigger_decoder
 *
 *
 */
