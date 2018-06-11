/* BOREXINO Reconstruction program
 *
 * Author: Alessandro Razeto <Alessandro.Razeto@ge.infn.it>
 * Maintainer: Alessandro Razeto <Alessandro.Razeto@ge.infn.it>
 *
 * $Id: laben_time_hit.hh,v 1.15 2008/10/07 14:03:20 razeto Exp $
 *
 * A class to decode the laben hit time: the obj is created empty,
 * with just some parameter, than it is initialized to raw hit by the
 * init method.
 * The second field in laben_time_hit::init (gray_cross) is intended
 * to be used when an event is in the gray counter cross reagion.
 *
 */
#ifndef _LABEN_TIME_HIT
#define _LABEN_TIME_HIT

#include "bx_laben_event.hh"
#include "bx_named.hh"

class laben_time_hit: public bx_named {
  public:
    enum ramp_slope { // Do not change the values!!!
      rising = 0,
      falling = 1,
    };

    laben_time_hit (bool use_calib_data = false);
    void init (const bx_laben_raw_hit& hit, bool gray_cross = false);
    
    ramp_slope guess_slope ();
    ramp_slope hw_slope () const { return i_hw_slope; } 
    double get_time ();	// Relative to a gray counter cycle window
    float get_d80 ();
    float get_error 	() const { return f_error; }
    bool is_valid 	() const { return b_valid; }
    bool is_good	() const { if (b_valid && f_error <= f_error_threshold) return true; return false; }
    bool is_good	(float threshold) const { if (b_valid && f_error <= threshold) return true; return false; }
    bool is_gray_crossing_window () { return (u4_gray_count < gray_crossing_window_low) || (u4_gray_count > gray_crossing_window_high); }
    
  private:
    bool b_use_calib_data;
    float f_d80;
    int i_low, i_high;
    float crate_time_correction;
    bool rising_on_even;
    int i_gr_shift;
    float f4_calib_offset;
      // some constants
    float f_ramp_slope, f_t20, f_t30, f_error_threshold;

    void m_calc_ramp_times (ramp_slope slope, float& t1, float& t2);

    uint8_t u1_ramp_1, u1_ramp_2;
    uint8_t u1_time_bits;
    uint32_t u4_gray_count;
    ramp_slope i_hw_slope;

    float f_error;
    bool b_valid;

    static const float ramp_length_ns;
    static const float window_upper_bound;
    static const float window_lower_bound;
    static const uint8_t spire_threshold;
    static const float d80;
    static const uint16_t gray_crossing_window_low, gray_crossing_window_high;
};


#endif
/*
 * $Log: laben_time_hit.hh,v $
 * Revision 1.15  2008/10/07 14:03:20  razeto
 * Added invalid flag to event (and removed checker from laben_time_hit)
 *
 * Revision 1.14  2007-11-05 11:03:08  razeto
 * Added static is_laben_invalid to laben_time_hit
 *
 * Revision 1.13  2005-11-20 18:11:40  razeto
 * Updated to a new scheme to use calibration
 *
 * Revision 1.12  2005/03/11 15:29:49  razeto
 * Added calibration time data usage in the decoder
 *
 * Revision 1.11  2005/03/03 15:11:27  razeto
 * Fixed a bug in handling gray counter circularity
 *
 * Revision 1.10  2005/03/02 17:32:49  razeto
 * Modified laben_{charge|time}_hit to speedup echidna avoiding to create too many objects
 *
 * Revision 1.9  2004/11/26 15:25:10  razeto
 * Thanks to Davide who so kindly underlined a syntactical error, I fixed the "Maintainer" word
 *
 * Revision 1.8  2004/11/26 14:01:19  razeto
 * Added Mantainer field
 *
 * Revision 1.7  2004/09/09 11:44:39  razeto
 * Upgraded to be a bx_named
 *
 * Revision 1.6  2004/05/21 08:42:23  razeto
 * Fixed a typo
 *
 * Revision 1.5  2004/05/20 10:51:00  razeto
 * Updated to support gray counter cycle crossing and to support the full time range
 *
 * Revision 1.4  2004/05/18 14:21:05  razeto
 * Fixed some bugs.
 * Added support for using laben hit time bin. Now the resolution is very
 * fine and no 50 or 100 ns mirror peaks are present.
 *
 * Revision 1.3  2004/04/27 09:41:09  ddangelo
 * Event structure redisegned.
 * bx_reco_event now handles event header and has object for different detector portions.
 * Inheritance from different reconstruction status is kept within each detector portion obj.
 * Hit lists are linked by const ptr (maybe a temporary solution).
 *
 * Revision 1.2  2004/04/24 17:37:40  razeto
 * Added a method for calculating d80
 *
 * Revision 1.1  2004/04/15 13:00:14  razeto
 * Added
 *
 */
