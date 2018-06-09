/* BOREXINO Reconstruction program
 *
 * Author: Alessandro Razeto <Alessandro.Razeto@ge.infn.it>
 * Maintainer: Alessandro Razeto <Alessandro.Razeto@ge.infn.it>
 *
 * $Id: laben_time_hit.cc,v 1.26 2009/07/16 16:09:43 razeto Exp $
 *
 * laben_time_hit implementation
 *
 */
#include "bx_dbi.hh"
#include "db_run.hh"
#include "db_profile.hh"
#include "constants.hh"
#include "laben_time_hit.hh"
#include "messenger.hh"

#include <math.h>

const float laben_time_hit::ramp_length_ns = 50.;
const float laben_time_hit::window_upper_bound = 17.;	// These values correspond to an error of 0.45
const float laben_time_hit::window_lower_bound = 14.;	// 		" "
const unsigned char laben_time_hit::spire_threshold = 10;
const float laben_time_hit::d80 = 80.;
  // Define the window bounds for a gray counter to be close to a cycle end/start.
const unsigned short laben_time_hit::gray_crossing_window_low = 800;
const unsigned short laben_time_hit::gray_crossing_window_high = (1UL << 16) - 800;

laben_time_hit::laben_time_hit (bool use_calib_data): bx_named ("laben_time_hit"), b_use_calib_data(use_calib_data) { }

void laben_time_hit::init (const bx_laben_raw_hit& hit, bool gray_cross) {
  f_error = 1.;
  if (hit.check_flag (bx_laben_raw_hit::invalid)) {
    b_valid = false;
    return;
  }
  else b_valid = true;
  
    // Get some parameters from the db
  int lg = hit.get_logical_channel ();
  const db_run& run_info = bx_dbi::get()->get_run ();
  const db_profile& profile_info = bx_dbi::get()->get_profile ();

  crate_time_correction = profile_info.laben_crate_delay (constants::crate (lg));

  if (!run_info.check_laben_precalib_low_bin (lg) || !run_info.check_laben_precalib_high_bin (lg))
    get_message (bx_message::critic) << "trying to run without ADC precalib data for channel " << lg << dispatch;
  i_low = run_info.get_laben_precalib_low_bin (lg);
  i_high = run_info.get_laben_precalib_high_bin (lg);
  if (i_low >= i_high) b_valid = false;
  
  if (!run_info.check_laben_precalib_gray_shift (lg)) 
    get_message (bx_message::critic) << "trying to run without GRAY SHIFT precalib data for channel " << lg << dispatch;
  i_gr_shift = run_info.get_laben_precalib_gray_shift (lg);

  if (b_use_calib_data) f4_calib_offset = run_info.get_laben_time_offset (lg);
  else f4_calib_offset = 0.;

    // The following parameters can not be in bd_run since laben_time_hit::guess_slope can 
    // be called during precalib. If not present use the default.
  if (run_info.check_laben_precalib_delta80 (lg)) f_d80 = run_info.get_laben_precalib_delta80 (lg);
  else f_d80 = d80;
  if (run_info.check_laben_precalib_rising_on_even (lg)) rising_on_even = run_info.get_laben_precalib_rising_on_even (lg);
  else rising_on_even = true;
 
    // Store the time data from the hit
  u1_ramp_1 = hit.get_time_1 ();
  u1_ramp_2 = hit.get_time_2 ();
  u1_time_bits = hit.get_flags_ch () & 0x3;
  u4_gray_count = hit.get_gray_counter ();
  if (gray_cross && u4_gray_count < gray_crossing_window_low) u4_gray_count += 1 << 16;

    // Calculate the hw slope  
  i_hw_slope = laben_time_hit::rising;			// suppose rising (avoid to much typing below)
  if (u4_gray_count & 0x1) {  	// odd
    if (rising_on_even) i_hw_slope = laben_time_hit::falling;
  } else {			// even
    if (!rising_on_even) i_hw_slope = laben_time_hit::falling;
  }

    // Some constants used later
  f_ramp_slope = (i_high - i_low) / ramp_length_ns;
  f_t20 = 2 * ramp_length_ns - f_d80;
  f_t30 = f_d80 - ramp_length_ns;
  f_error_threshold = window_lower_bound / (window_lower_bound + window_upper_bound);
  
  if (!b_valid) f_error = .5;
} 

double laben_time_hit::get_time () {
  
    // If no predecoding info are loaded for this channel no calculations can be done, return 0.
  if (!is_valid ()) return 0.;

    // Calculate the gray counter time
  double time_gray = ramp_length_ns * (double(u4_gray_count) - i_gr_shift);
  
    // Use guess_slope returned slope instead of hw_slope. 
  laben_time_hit::ramp_slope slope = guess_slope ();

    // Calculate the 50ns relative time of the the 2 hits
  float t1, t2;
  m_calc_ramp_times (slope, t1, t2);

    // If guessed slope and hw slope does not coincide suppose the gray counter is late.
  if (slope != i_hw_slope) {
      // Use the information of slope bits that even if not well understood allows a 100ns peak of 99%.
      // The code is only an IF with a complicated conditions: the flags can be 0 or 3 for the 100ns
      // peak condition and there are 2 combinations dependig on rising_on_even flag.
      // if rising_on_even {
      //   if (rising && 3) -50ns
      //   else if (falling && 0) -50ns
      //   else +50ns
      // } else {
      //   if (falling && 3) -50ns
      //   else if (rising && 0) -50ns
      //   else +50ns
      // }
      // The following if just reduces the number of the if intructions shown above with just 1 test 
      // (and a lot of bolean math).
    if ((!rising_on_even && 
	  ( ((slope == laben_time_hit::rising) && (u1_time_bits == 3)) ||
	    ((slope == laben_time_hit::falling) && (u1_time_bits == 0)) )
	) || (rising_on_even &&
	  ( ((slope == laben_time_hit::falling) && (u1_time_bits == 3)) ||
	    ((slope == laben_time_hit::rising) && (u1_time_bits == 0)) )
	))
      time_gray -= ramp_length_ns;
    else time_gray += ramp_length_ns;
  }
  
    // If sample 1 is to close to spires use sample 2
  if (u1_ramp_1 < (i_low + spire_threshold) || u1_ramp_1 > (i_high + spire_threshold))
    return time_gray + t2 - f_d80 + crate_time_correction - f4_calib_offset;
    // Else use sample 1 (which is more accurate than d80)
  return time_gray + t1 + crate_time_correction - f4_calib_offset;
}

void laben_time_hit::m_calc_ramp_times (laben_time_hit::ramp_slope slope, float& t1, float& t2) {
  if (slope == laben_time_hit::rising) {
    t1 = (u1_ramp_1 - i_low) / f_ramp_slope;
    if (t1 < f_t20) t2 = ramp_length_ns + (i_high - u1_ramp_2) / f_ramp_slope;
    else t2 = 2 * ramp_length_ns + (u1_ramp_2 - i_low) / f_ramp_slope;
  } else {
    t1 = (i_high - u1_ramp_1) / f_ramp_slope;
    if (t1 < f_t20) t2 = ramp_length_ns + (u1_ramp_2 - i_low) / f_ramp_slope;
    else t2 = 2 * ramp_length_ns + (i_high - u1_ramp_2) / f_ramp_slope;
  }
}

float laben_time_hit::get_d80 () {
  
    // If no predecoding info are loaded for this channel no calculations can be done, return 0.
  if (!is_valid ()) return 0.;

    // Use guess_slope returned slope instead of hw_slope. 
  laben_time_hit::ramp_slope slope = guess_slope ();

    // Calculate the 50ns relative time of the the 2 hits
  float t1, t2;
  m_calc_ramp_times (slope, t1, t2);

  return t2 - t1;
}

laben_time_hit::ramp_slope laben_time_hit::guess_slope () {

    // If no predecoding info are loaded for this channel use the graycounte parity
  if (!is_valid ()) return i_hw_slope;

    // Sistances from the second sample for the 2 projections
  float d_rising, d_falling;

    // Rising ramp  
  float t = (u1_ramp_1 - i_low) / f_ramp_slope; 	// the time of the first sample
  if (t > ramp_length_ns) t = ramp_length_ns;
  else if (t < 0) t = 0;

  if (t < f_t20) d_rising = ::fabs (i_high - (t + f_t30) * f_ramp_slope - u1_ramp_2);
  else d_rising = ::fabs (i_low + (t - f_t20) * f_ramp_slope - u1_ramp_2);

    // Falling ramp
  t = ramp_length_ns - t;				// the time of the first sample
  if (t < f_t20) d_falling = ::fabs (i_low + (t + f_t30) * f_ramp_slope - u1_ramp_2);
  else d_falling = ::fabs (i_high - (t - f_t20) * f_ramp_slope - u1_ramp_2);

    // Calculate the error and the slope
  float d_min, d_max;
  laben_time_hit::ramp_slope guess_slope;
  if (d_falling > d_rising) { 
    d_max = d_falling;
    d_min = d_rising;
    guess_slope = laben_time_hit::rising;
  } else {
    d_min = d_falling;
    d_max = d_rising;
    guess_slope = laben_time_hit::falling;
  }
  if (d_max == 0.) f_error = 0.5;	   // if max = 0 even min = 0 so the error is maximal
  else f_error = d_min / ( d_min + d_max); // else use the pondered distance
  
    // If the one and only if min is in the lower bound and max is > the the upper bound use
    // the found slope
  if (is_good () && d_max > window_upper_bound) {
    if (guess_slope == i_hw_slope) f_error -= 0.5;
    return guess_slope;
  } else {
    f_error += 0.5;
    return i_hw_slope;
  }
} 

/*
 * $Log: laben_time_hit.cc,v $
 * Revision 1.26  2009/07/16 16:09:43  razeto
 * Removed a wrong assumption in the heuristic code for guessing (does not affects results)
 *
 * Revision 1.25  2008-10-07 14:03:20  razeto
 * Added invalid flag to event (and removed checker from laben_time_hit)
 *
 * Revision 1.24  2007-11-05 11:03:08  razeto
 * Added static is_laben_invalid to laben_time_hit
 *
 * Revision 1.23  2007-10-30 17:36:18  razeto
 * Gray cross only really close to edge
 *
 * Revision 1.22  2006-10-13 15:31:18  razeto
 * Added variables from raw data
 *
 * Revision 1.21  2006-07-13 14:27:30  razeto
 * Upgraded gray cross handling (enlarged again window, lowered probability requirement and add debug histo)
 *
 * Revision 1.20  2006/07/13 12:54:27  razeto
 * Enlarged gray cross limits from 16 to 25us
 *
 * Revision 1.19  2005/11/20 18:11:40  razeto
 * Updated to a new scheme to use calibration
 *
 * Revision 1.18  2005/08/24 14:35:34  razeto
 * Fixed a typo
 *
 * Revision 1.17  2005/05/05 16:58:49  razeto
 * Adapted to the new interface
 *
 * Revision 1.16  2005/03/11 15:29:49  razeto
 * Added calibration time data usage in the decoder
 *
 * Revision 1.15  2005/03/03 15:29:27  razeto
 * Better handling of anomalous zero boundary
 *
 * Revision 1.14  2005/03/03 15:19:38  razeto
 * Moved validity test on the top, changed the error range
 *
 * Revision 1.13  2005/03/02 17:32:49  razeto
 * Modified laben_{charge|time}_hit to speedup echidna avoiding to create too many objects
 *
 * Revision 1.12  2005/03/01 15:12:41  razeto
 * Merged with cycle_2
 *
 * Revision 1.11.2.1  2005/02/25 17:00:00  razeto
 * Upgraded adding some sanity checks on guess_slope (required for bx_elec)
 *
 * Revision 1.11  2004/11/26 15:25:10  razeto
 * Thanks to Davide who so kindly underlined a syntactical error, I fixed the "Maintainer" word
 *
 * Revision 1.10  2004/11/26 14:01:19  razeto
 * Added Mantainer field
 *
 * Revision 1.9  2004/09/09 11:44:39  razeto
 * Upgraded to be a bx_named
 *
 * Revision 1.8  2004/08/06 10:30:28  razeto
 * cycle_1 branch merged in the main trunk.
 *
 * Revision 1.7.2.1  2004/08/06 10:01:22  razeto
 * Fixed the for the gray shift correction
 *
 * Revision 1.7  2004/05/28 13:48:11  razeto
 * Added some test for precalib data
 *
 * Revision 1.6  2004/05/21 11:28:42  razeto
 * Updated conforming to the new conventions of event::get_logical_channel
 *
 * Revision 1.5  2004/05/20 10:51:00  razeto
 * Updated to support gray counter cycle crossing and to support the full time range
 *
 * Revision 1.4  2004/05/18 14:21:05  razeto
 * Fixed some bugs.
 * Added support for using laben hit time bin. Now the resolution is very
 * fine and no 50 or 100 ns mirror peaks are present.
 *
 * Revision 1.3  2004/04/26 13:53:22  razeto
 * Updated to follow new name/syntax of db_run and db_profile
 *
 * Revision 1.2  2004/04/24 17:37:40  razeto
 * Added a method for calculating d80
 *
 * Revision 1.1  2004/04/15 13:00:14  razeto
 * Added
 *
 */
