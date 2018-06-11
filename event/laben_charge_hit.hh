/* BOREXINO Reconstruction program
 *
 * Author: Alessandro Razeto <Alessandro.Razeto@ge.infn.it>
 * Maintainer: Alessandro Razeto <Alessandro.Razeto@ge.infn.it>
 *
 * $Id: laben_charge_hit.hh,v 1.9 2008/12/10 11:40:24 razeto Exp $
 *
 * A class to decode the laben hit charge, the interface (ctor and 
 * init) is the same than laben_time_hit.
 *
 */
#ifndef _LABEN_CHARGE_HIT
#define _LABEN_CHARGE_HIT

#include "bx_laben_event.hh"
#include "bx_named.hh"

class laben_charge_hit: public bx_named {
  public:
    laben_charge_hit (bool use_calib_data = false);

    void init (const bx_laben_decoded_hit& hit);
    void init (const bx_laben_decoded_hit& hit, const bx_laben_decoded_hit& previus_hit);

    float get_uncorrected_charge_bin () { return u1_peak - u1_base; }
    float get_charge_bin ();

    float get_uncorrected_charge () { return get_uncorrected_charge_bin () / f4_channel_bin_to_charge; }
    float get_charge () { return get_charge_bin () / f4_channel_bin_to_charge; }
    int get_npe ();

    float get_charge_mean () { return get_charge_bin () / f4_channel_bin_to_charge_mean; }
  private:
    int b_use_calib_data;
      // Constants from parameters
    float f4_baseline_fast_drift;
    float f4_zero_charge_limit;

      // Hit data
    uint8_t u1_peak, u1_base;
    float f4_prev_dt;
    uint8_t u1_prev_peak, u1_prev_base;
    float f4_channel_bin_to_charge, f4_channel_bin_to_charge_mean;

      // Some real constants
    static const float default_bin_to_charge;
    static const float integrator_decay_time;
    static const float restoration_time, d80;

    void get_calib (int lg);
};


#endif
/*
 * $Log: laben_charge_hit.hh,v $
 * Revision 1.9  2008/12/10 11:40:24  razeto
 * Added mean charge and peak from tt1
 *
 * Revision 1.8  2005-11-20 18:12:05  razeto
 * Updated to a new scheme to use calibration and added calibration code
 *
 * Revision 1.7  2005/11/18 16:43:33  razeto
 * Added calculation of integer charge
 *
 * Revision 1.6  2005/03/02 17:32:49  razeto
 * Modified laben_{charge|time}_hit to speedup echidna avoiding to create too many objects
 *
 * Revision 1.5  2004/11/26 15:25:10  razeto
 * Thanks to Davide who so kindly underlined a syntactical error, I fixed the "Maintainer" word
 *
 * Revision 1.4  2004/11/26 14:01:19  razeto
 * Added Mantainer field
 *
 * Revision 1.3  2004/09/09 11:44:28  razeto
 * Upgraded to be a bx_named; added pileup calculation
 *
 * Revision 1.2  2004/05/25 17:18:20  razeto
 * Upgraded
 *
 * Revision 1.1  2004/05/20 10:51:34  razeto
 * Added sketch version
 *
 */
