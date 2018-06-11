/* BOREXINO Reconstruction program
 *
 * Author: Alessandro Razeto <Alessandro.Razeto@ge.infn.it>
 * Maintainer: Alessandro Razeto <Alessandro.Razeto@ge.infn.it>
 *
 * $Id: bx_laben_decoder.hh,v 1.25 2017/05/12 14:32:38 misiaszek Exp $
 *
 * Decoder for the laben event
 *
 */
#ifndef _BX_LABEN_DECODER_H
#define _BX_LABEN_DECODER_H

#include "bx_base_module.hh"
#include "bx_echidna_event.hh"
#include "bx_rec_general.hh"
#include "TH1F.h"
#include "TH2F.h"

#include <vector>
#include <list>
class db_channel_laben;

class bx_laben_decoder: public bx_base_module {
  public:
    bx_laben_decoder ();
    virtual ~bx_laben_decoder () {}

    virtual void begin ();
    virtual bx_echidna_event* doit (bx_echidna_event *ev);
    virtual void end ();
  private:
    float m_check_gray_cross (const bx_laben_event& er);
    double m_get_moda_mean (const std::vector<double>& v);
    double m_search_moda (const double *times, uint16_t size);

    typedef std::list<bx_laben_decoded_hit> bx_laben_decoded_hit_list;
    void remove_retrigger_hits (bx_laben_decoded_hit_list&);
    void remove_cable_reflection_hits (bx_laben_decoded_hit_list&);

    bool b_discard_out_of_gate_hits, b_discard_reference_hits, b_discard_calib_data, b_warn_empy_channel, b_mc_use_charge_calib_data;
    bool b_discard_retrigger_hits, b_discard_reflection_hits, b_discard_disabled_lg;
    int32_t i4_nhits_threshold;
    float f4_trigger_start_;
    float f4_shift_min, f4_shift_max;
    uint8_t *p_disabled_lg;
    bool *p_empty_lg;
    uint32_t i4_ordinary_pmt, i4_n_disabled_channels, i4_n_disabled_charge, i4_n_disabled_pmts, i4_n_disabled_pmts_charge;
    const db_channel_laben **ch_info_v; 
    TH1F *trigger_charge, *laser_charge;
    TH2F *gray_cross_h;
};

#endif
/*
 * $Log: bx_laben_decoder.hh,v $
 * Revision 1.25  2017/05/12 14:32:38  misiaszek
 * denis patch
 *
 *
 * Revision 1.24  2012/03/22 19:19:57  razeto
 * Added nhits_threshold
 *
 * Revision 1.23  2011-03-24 15:10:19  razeto
 * Fixed hit name
 *
 * Revision 1.22  2011-01-19 10:42:44  davini
 * added the possibility to use charge calibration data for channels in montecarlo
 *
 * Revision 1.21  2008-09-25 13:11:32  razeto
 * Zeroing charge for charge disabled events
 *
 * Revision 1.20  2007-12-10 11:17:56  razeto
 * use only ordinary channels as base for n_live_pmts
 *
 * Revision 1.19  2007-11-05 23:39:53  razeto
 * hits_on_empty filled
 *
 * Revision 1.18  2007-10-30 17:34:26  razeto
 * Experimental code:
 * - better algorithm for finding tigger time
 * - handle better gray window
 * - handle large gate for neutron events
 *
 * Revision 1.17  2007-05-08 23:54:23  razeto
 * Live pmt calculation integrated
 *
 * Revision 1.16  2007-05-02 16:38:57  razeto
 * Updated to new bx_detector: now channels can be disabled runtime (to be tested)
 *
 * Revision 1.15  2006/11/27 10:59:01  razeto
 * Upgraded bx_detector to have more detailed check on bad channels multiplicity:
 * now bad_multiplicity do not exist any more but there are other fields.
 * Upgraded laben_decoder: the fix_bad_channel flag do not exists anymore. To
 * avoid bad channels skipping, set detector_interface.bad_laben_channel_properties
 * to {}.
 * Configurations updated.
 *
 * Revision 1.14  2006/11/16 10:08:15  razeto
 * Experimental version of decoder with new peaks suppression (to be debugged)
 *
 * Revision 1.13  2006/09/09 14:08:21  razeto
 * Upgraded to remove retrigger hits and added a fix_bad_channels parameter
 *
 * Revision 1.12  2006/07/13 14:27:30  razeto
 * Upgraded gray cross handling (enlarged again window, lowered probability requirement and add debug histo)
 *
 * Revision 1.11  2006/03/21 15:03:45  razeto
 * Now time resolution is better calculated in bx_calib_laben_decoding
 *
 * Revision 1.10  2005/06/27 16:15:02  razeto
 * Added a new histogram to have the detector time resolution witout recalculating the precalibrations
 *
 * Revision 1.9  2005/03/11 15:29:48  razeto
 * Added calibration time data usage in the decoder
 *
 * Revision 1.8  2005/03/03 09:30:36  razeto
 * Added a parameter to even copy reference channels hits and
 * fixed the relevant bug.
 * Changed the window for out of gate hits
 *
 * Revision 1.7  2004/11/26 15:25:11  razeto
 * Thanks to Davide who so kindly underlined a syntactical error, I fixed the "Maintainer" word
 *
 * Revision 1.6  2004/11/26 14:01:20  razeto
 * Added Mantainer field
 *
 * Revision 1.5  2004/09/22 13:26:19  razeto
 * Changed bx_reco_event to bx_echidna_event
 *
 * Revision 1.4  2004/09/09 11:53:18  razeto
 * Doit method mainly rewritten; now several loops on the hit matrix
 * are done now to allow a better (and hopefully cleaner) decoding.
 * New decoded event format introduced (to fully have time reference
 * hits).
 * Decoder is now almost ready to recognise time reference hits from
 * predecoding hits using the charge samples (for runid > 1800).
 *
 * Revision 1.3  2004/07/14 11:49:00  razeto
 * Fixed a bug (array overflow) and upgraded some routines
 *
 * Revision 1.2  2004/06/01 15:38:25  razeto
 * A lot of development to support time references, out_of_gate hits and much more
 *
 * Revision 1.1  2004/05/21 08:40:33  razeto
 * Added
 *
 */
