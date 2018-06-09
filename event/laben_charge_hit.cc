/* BOREXINO Reconstruction program
 *
 * Author: Alessandro Razeto <Alessandro.Razeto@ge.infn.it>
 * Maintainer: Alessandro Razeto <Alessandro.Razeto@ge.infn.it>
 *
 * $Id: laben_charge_hit.cc,v 1.16 2011/01/19 15:13:58 davini Exp $
 *
 * laben_charge_hit implementation
 *
 */
#include "laben_charge_hit.hh"
#include "messenger.hh"
#include "bx_dbi.hh"
#include "db_run.hh"
#include <math.h>

const float laben_charge_hit::default_bin_to_charge = 26.; //era 22
const float laben_charge_hit::integrator_decay_time = 500.;
const float laben_charge_hit::restoration_time = 2700;  // e^-(2800/500) < 1/255
const float laben_charge_hit::d80 = 80;

laben_charge_hit::laben_charge_hit (bool use_calib_data): bx_named ("laben_charge_hit"), b_use_calib_data(use_calib_data) { 
  f4_baseline_fast_drift = get_parameter ("baseline_fast_drift").get_float (); 
  f4_zero_charge_limit = get_parameter ("zero_charge_limit").get_float (); 
}

void laben_charge_hit::get_calib (int lg) {
  const db_run& run_info = bx_dbi::get()->get_run ();
  float calib_peak = run_info.get_laben_charge_tt1_peak (lg);
  float calib_mean = run_info.get_laben_charge_tt1_mean (lg);
  if (b_use_calib_data && calib_peak > 0) {
    f4_channel_bin_to_charge = calib_peak;
    f4_channel_bin_to_charge_mean = calib_mean;
  } else f4_channel_bin_to_charge_mean = f4_channel_bin_to_charge = default_bin_to_charge;
}

void laben_charge_hit::init (const bx_laben_decoded_hit& hit) {
  u1_peak = hit.get_raw_hit ().get_peak ();
  u1_base = hit.get_raw_hit ().get_base ();
  f4_prev_dt = restoration_time;
  get_calib (hit.get_raw_hit ().get_logical_channel ());  
}	

void laben_charge_hit::init (const bx_laben_decoded_hit& hit, const bx_laben_decoded_hit& previus_hit) {
  int lg = hit.get_raw_hit ().get_logical_channel ();
  int lg_prev = previus_hit.get_raw_hit ().get_logical_channel ();

  if (lg != lg_prev) 
    get_message (bx_message::critic) << "laben_charge_hit constructed with hit of different channels " << lg << "!=" << lg_prev << dispatch;

  f4_prev_dt = hit.get_raw_time () - previus_hit.get_raw_time ();
  if (f4_prev_dt < 0) 
    get_message (bx_message::critic) << "laben_charge_hit construced with wrong time sorted hit " << hit.get_raw_time () << "<" << previus_hit.get_raw_time () << dispatch;
  
  u1_peak = hit.get_raw_hit ().get_peak ();
  u1_base = hit.get_raw_hit ().get_base ();
  u1_prev_peak = previus_hit.get_raw_hit ().get_peak ();
  u1_prev_base = previus_hit.get_raw_hit ().get_base ();
  get_calib (lg);
}

float laben_charge_hit::get_charge_bin () {
  float simple_charge = u1_peak - u1_base;
  if (simple_charge < 0) simple_charge = 0;

    // Some simple check to see if a overlapping routine is to be used
  if (f4_prev_dt >= restoration_time) return simple_charge; 	// keep >= since one hit ctor sets f4_prev_dt = restoration_time
  if (f4_prev_dt < d80) return simple_charge; 			// It should not happen since dt for good channels should be > 140ns
  if (u1_prev_peak <= u1_prev_base) return simple_charge; 	// No info can be found on prev hit
  
    // A more detailed check to see if the actual base sample is too distant from where it should
    // be (simply considering integrator discarge. If not the baseline is too noisy.
  float prev_simple_charge = u1_prev_peak - u1_prev_base;
  float base_check = prev_simple_charge * ::expf ((d80 - f4_prev_dt) / integrator_decay_time) + u1_prev_base;
  if (::fabsf(base_check - u1_base) > f4_baseline_fast_drift) return simple_charge; 

    // If here all check were sucesfull so calculate pile-up
  float prev_base_delayed = prev_simple_charge * ::expf (-f4_prev_dt / integrator_decay_time) + u1_prev_base;
  return (u1_peak > prev_base_delayed) ? u1_peak - prev_base_delayed : 0;
}

int laben_charge_hit::get_npe () {
  if (get_charge () < f4_zero_charge_limit) return 0;
  if (f4_zero_charge_limit < 0.5 && get_charge () < 0.5) return 1; // For values lower than 0.5 do not return 0
  return int (::roundf (get_charge ()));
}
/*
 * $Log: laben_charge_hit.cc,v $
 * Revision 1.16  2011/01/19 15:13:58  davini
 * undo last change
 *
 * Revision 1.15  2011-01-19 10:46:11  davini
 * integrator_decay_time = 235.
 *
 * Revision 1.14  2008-12-10 11:40:24  razeto
 * Added mean charge and peak from tt1
 *
 * Revision 1.13  2008-07-11 08:40:30  razeto
 * Now default charge is 26 (according to data and simulations)
 *
 * Revision 1.12  2007-03-29 09:00:54  razeto
 * Never have a negative charge, npe always >= 1
 *
 * Revision 1.11  2006/02/28 12:54:26  testera
 * Set the default charge bin to the same value as bx_elec (by Alessandro)
 *
 * Revision 1.10  2005/11/20 18:12:05  razeto
 * Updated to a new scheme to use calibration and added calibration code
 *
 * Revision 1.9  2005/11/18 16:43:33  razeto
 * Added calculation of integer charge
 *
 * Revision 1.8  2005/03/02 17:32:49  razeto
 * Modified laben_{charge|time}_hit to speedup echidna avoiding to create too many objects
 *
 * Revision 1.7  2004/11/26 15:25:10  razeto
 * Thanks to Davide who so kindly underlined a syntactical error, I fixed the "Maintainer" word
 *
 * Revision 1.6  2004/11/26 14:01:19  razeto
 * Added Mantainer field
 *
 * Revision 1.5  2004/09/12 12:57:15  razeto
 * Added math.h include
 *
 * Revision 1.4  2004/09/09 11:44:28  razeto
 * Upgraded to be a bx_named; added pileup calculation
 *
 * Revision 1.3  2004/06/08 10:25:19  razeto
 * Updated to new hit name (rawt -> raw_time)
 *
 * Revision 1.2  2004/05/25 17:18:20  razeto
 * Upgraded
 *
 * Revision 1.1  2004/05/20 10:51:34  razeto
 * Added sketch version
 *
 */
