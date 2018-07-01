/* BOREXINO Reconstruction program
 *
 * Author: Davide D'Angelo <davide.dangelo@lngs.infn.it>
 * based on work by Razeto&Pallas 
 * Maintainer: Davide D'Angelo <davide.dangelo@lngs.infn.it>
 *
 * $Id: bx_laben_event.cc,v 1.63 2015/07/28 09:14:54 misiaszek Exp $
 *
 * Implementation of bx_laben_event
 *
 */
#include "bx_laben_event.hh"
#include "bx_mctruth_event.hh"
#include "bx_event_disk_format.h"
#include "messenger.hh"
#include "laben_time_hit.hh"
#include "laben_charge_hit.hh"
#include "bx_detector.hh"

bx_laben_raw_event::bx_laben_raw_event (const char *disk_event): i4_nhits_fw(0) {
  if (!detector_interface::get ()->is_laben_enabled ()) return;

  laben_header_disk_format *head = (laben_header_disk_format *)disk_event;

  if (head->entries) u4_errors = head->errors;

  raw_hits.reserve (head->entries);
  int last_channel = -1;
  uint8_t hit_count = 1;
  for (uint32_t i = 0; i < head->entries; i++) {
    bx_laben_raw_hit hit (disk_event + sizeof (laben_header_disk_format) + i * sizeof (laben_hit_disk_format), hit_count);
    raw_hits.push_back (hit);
    if (hit.get_logical_channel () == last_channel) hit_count++;
    else { hit_count = 1; last_channel = hit.get_logical_channel (); }
  }		
}

uint16_t bx_laben_raw_hit::flags_bits[__max__] = { 0x8770, 0x40, 0x20, 0x10, 0x100, 0x200, 0x400, 0x8000 }; // KEEP alligned to flags enum !!!

bx_laben_raw_hit::bx_laben_raw_hit (const char *disk_event, uint8_t order_in_channel) {
  
  laben_hit_disk_format *hit = (laben_hit_disk_format *)disk_event;

  u2_channel = hit->channel;
  u1_time_1 = hit->time_s[0];
  u1_time_2 = hit->time_s[1];
  u2_gray_counter = hit->time_l;
  u1_base = hit->charge[0];
  u1_peak = hit->charge[1];
  u2_flags = hit->flags & 0x760; // VALUE FROM ENUM, only possible good values of hw bits
  if (u2_flags & flags_bits[fifo_full] && u2_flags & flags_bits[fifo_empty]) 
    u2_flags = (u2_flags & ~(flags_bits[fifo_full] | flags_bits[fifo_empty])) | flags_bits[counter];
  else if (u1_time_1 == 0xFF && u1_time_2 == 0xFF && u2_gray_counter == 0xFFFF) u2_flags |= flags_bits[invalid];
  u2_flags |= hit->flags & 0x3;
//  u2_errors = hit->errors; UNUSED
  u1_order_in_channel = order_in_channel;
}

bx_laben_decoded_event::bx_laben_decoded_event () {
}

const bx_laben_decoded_hit& bx_laben_decoded_hit::operator= (laben_time_hit& t_hit) {
  f8_raw_time = t_hit.get_time ();
  f4_d80 = t_hit.get_d80 ();
  f4_time_error = t_hit.get_error ();

  return *this;
}

const bx_laben_decoded_hit& bx_laben_decoded_hit::operator= (laben_charge_hit& c_hit) {
  f4_charge_bin = c_hit.get_charge_bin ();
  f4_uncorrected_charge_bin = c_hit.get_uncorrected_charge_bin ();
  i4_charge_npe = c_hit.get_npe ();
  f4_charge_pe = c_hit.get_charge ();
  f4_uncorrected_charge_pe = c_hit.get_uncorrected_charge ();
  f4_charge_mean_pe = c_hit.get_charge_mean ();

  return *this;
}

bx_laben_clustered_event::bx_laben_clustered_event() : i4_nclusters_found (0), i4_nclusters_old(0), i4_nclusters_neutron(0){  
}

bx_laben_cluster::bx_laben_cluster () : i4_clustered_nhits_conc (0), 
                                        i4_clustered_nhits_thresh (0), 
					i4_clustered_nhits_short (0),
					i4_clustered_nhits_400 (0),
					f4_clustered_nhits_bkg (0.),
					f4_charge (0.), 
					f4_charge_conc (0.), 
					f4_charge_mean (0.), 
					f4_charge_thresh (0.), 
					f4_charge_short (0.),
					f4_charge_400 (0.),
					f4_charge_dt1(0),
					f4_charge_dt2(0),
					f4_charge_npmts (0.),
					f4_charge_clean (0.),
					f4_charge_noavg_dt1 (0.),
					f4_charge_noavg_dt2 (0.),
					f4_charge_noavg (0.),
					f4_charge_noavg_short (0.),
					i4_npe (0), 
					i4_npe_conc (0),
					i4_npmts (0),
					i4_npmts_conc (0),
					i4_npmts_thresh (0),
					i4_npmts_short (0),
					i4_npmts_400 (0),
					i4_npmts_dt1(0), 
					i4_npmts_dt2(0), 
					f8_start_time (0.),
//					f8_rough_time (0.),
					f4_mean_time (0.),
					f4_mean_time_short (0),
					f4_rms_time (0.),
					f4_rms_time_short (0.),
					f4_duration_short (0.),
					u1_flag (0),
			                b_is_neutron (0) {
}

const bx_base_position& bx_position::operator= (const bx_mctruth_frame& frame) {
  f4_x = frame.get_position (0);
  f4_y = frame.get_position (1);
  f4_z = frame.get_position (2);
  f4_user = f4_t = i4_matrix = 0;
  b_converged = false;

  return *(bx_base_position *)this;
}

bx_laben_shaped_cluster::bx_laben_shaped_cluster (const bx_laben_cluster &cluster): bx_laben_rec_cluster(cluster), 
										    f4_ns_asymmetry(0.),
										    f4_sphere_chi2(0.),
										    f4_sphere_lkl(0.),
										    f4_sphere_rel_var(0.),
										    f4_plane_cos(0.),
										    f4_plane_chi2(0.),
										    f4_h_plane_chi2(0.),
										    i1_quality_flags(0) {
  for (int i = 0; i<4; i++) v_sh_power[i] = 0.;
}

float bx_laben_ab_cluster::get_tailtot(int tail) const {
  if ( tail < 40 || tail > 130 )
    return -1.;
  int index = (tail - 40)/10;
  return v_tailtot[index];
}

float bx_laben_ab_cluster::get_tailtot_ab_mlp(int tail_index) const {
  if ( tail_index < 0 || tail_index > 9 )
    return -1.;
  return v_tailtot_ab_mlp[tail_index];
}

float bx_laben_ab_cluster::get_tailtot_c11_mva(int tail_index) const {
  if ( tail_index < 0 || tail_index > 9 )
    return -1.;
  return v_tailtot_c11_mva[tail_index];
}



/*float bx_laben_ab_mach4_cluster::get_tailtot_mach4(int tail) const {
  if ( tail < 30 || tail > 110 )
    return -1.;
  int index = (tail - 30)/5;
  return v_tailtot_mach4[index];
}*/

bx_laben_positron_cluster::bx_laben_positron_cluster (const bx_laben_cluster &cluster) : bx_laben_ab_cluster(cluster),
											f4_gatti_ops_beta(-100.),
											f4_gatti_c11_beta(-100.),
											f4_gatti_ops_nops(-100.){
}



/*
 * $Log: bx_laben_event.cc,v $
 * Revision 1.63  2015/07/28 09:14:54  misiaszek
 * tailtot for c11/b discrimination added
 *
 * Revision 1.62  2015/07/14 09:42:21  misiaszek
 * f4_charge_noavg_dt1,  f4_charge_noavg_dt2,  f4_charge_noavg , f4_charge_noavg_short added
 *
 * Revision 1.61  2015/07/13 15:07:14  misiaszek
 * tailtot_ab_mlp for a/b discrimination with MLP algorithm added
 *
 * Revision 1.60  2014/12/11 21:27:12  wurm
 * added getters and variables to BxEvent that concern C-11 likelihood tagging
 * * for tracks: new distance-to- and projection-on-track getters, dE/dx calculation
 * * counting of neutron-like clusters in tt128, new flag of neutron clusters
 *
 * Revision 1.59  2013/01/22 09:22:17  mosteiro
 * dt1,2 vars calculated for each cluster
 *
 * Revision 1.58  2012-10-30 15:39:17  ddangelo
 * fixed typos
 *
 * Revision 1.57  2011-04-13 10:37:35  ddangelo
 * added variable nclusters_old
 *
 * Revision 1.56  2011-03-23 16:04:48  ddangelo
 * added charge_clean, zeroing, copy. internal and root event
 *
 * Revision 1.55  2011-02-19 13:29:11  davini
 * bx_pid_positron stuffs
 *
 * Revision 1.54  2011-02-18 18:22:01  ddangelo
 * event support to mach4 a/b discrimination commented out. variable writing with the module commented out.
 * added event support for gatti ops. to be improved.
 *
 * Revision 1.53  2011-02-18 16:04:27  davini
 * added npmts_400 nhits_400 charge_400 on Oleg's request; added charge_npmts based on Alessandro, Stefano and Livia idea;
 *
 * Revision 1.52  2011-01-15 20:57:40  razeto
 * Use label for max
 *
 * Revision 1.51  2010-07-23 21:28:58  razeto
 * time bits restored
 *
 * Revision 1.50  2010-07-01 18:21:36  razeto
 * Added counter hit flag for new laben fw and nhits_fw from laben boards
 *
 * Revision 1.49  2010-07-01 13:28:25  ddangelo
 * nhits_bkg changed to float to accomodate info for tt1
 *
 * Revision 1.48  2009-11-18 11:43:31  ddangelo
 * added rms_time, duration and npmts short version for back compatibility.
 * npmt variables renamed to npmts also in internal event.
 * mach4_n1700 renamed as mach4_fixed throughout the event classes.
 *
 * Revision 1.47  2009-10-26 19:17:23  ddangelo
 * added bx_position_mach4_n1700 class (internal and root event, copy and getters)
 * in (base)postion class variable likelihood renamed as user. parallel getters to use it as likelihood or refraction index
 *
 * Revision 1.46  2009-10-22 16:06:46  ddangelo
 * in bx_laben_cluster added nhits_bkg for evaluate of bkg contribution in case of piled-up events.
 * internal and root event, copy and getters
 *
 * Revision 1.45  2009-10-08 15:45:59  ddangelo
 * implemented pid event variables removal and addition.
 * Internal and root, inizialization, copy and getters.
 *
 * Revision 1.44  2009-10-06 17:05:58  razeto
 * Use unused bit for invalid (0x8000) and mask. Check_flag fixed
 *
 * Revision 1.43  2009-09-18 22:20:12  ddangelo
 * added sphere_rel_var to laben shaped cluster. Internal and root event, inizialization and copy.
 *
 * Revision 1.42  2009-09-17 15:54:29  razeto
 * Initialize the short variables (non maintainer commit)
 *
 * Revision 1.41  2009-07-22 10:41:22  ddangelo
 * in position classes, both internal and root event:
 * - n_iterations removed
 * + matrix and converged variables added with getters.
 *
 * Revision 1.40  2009-07-17 15:39:51  ddangelo
 * laben cluster:
 * + added npmts_thresh, nhits_thresh and charge_thresh, to be computed with >0.2pe hits only (internal and root event).
 * - removed cluster rough time (internal and root).
 * - removed npe and npe_conc (root only)
 *
 * Revision 1.39  2009-07-16 15:17:57  ddangelo
 * mach a/b ported to root event
 * debugging tailtot getter for m4
 * other debugging
 *
 * Revision 1.38  2009-07-16 10:50:50  ddangelo
 * new class for mach4 a/b discrimination
 *
 * Revision 1.37  2008-12-15 12:12:44  razeto
 * Added rms_time to cluster (to improve PID)
 *
 * Revision 1.36  2008-12-11 17:12:25  razeto
 * Do not remove low level flags (not to break time_decondig)
 *
 * Revision 1.35  2008-12-10 11:40:24  razeto
 * Added mean charge and peak from tt1
 *
 * Revision 1.34  2008-11-19 11:27:24  razeto
 * Copy only flags (fixes a bug)
 *
 * Revision 1.33  2008-10-07 14:03:20  razeto
 * Added invalid flag to event (and removed checker from laben_time_hit)
 *
 * Revision 1.32  2008-06-20 16:23:41  razeto
 * Initialize the right variable (small bug fix)
 *
 * Revision 1.31  2007-11-12 15:19:52  razeto
 * Added cluster flags and removed broad variable (davide auth)
 *
 * Revision 1.30  2007-11-05 23:39:02  razeto
 * Added new hits_on_empty variable to laben data
 *
 * Revision 1.29  2007-10-30 15:45:02  ddangelo
 * added # of laben cluster found by algorythm (different from saved one for high multiplicity events)
 * added end hit time of a cluster
 * internal and root event. getters, copy, initialization, etc...
 *
 * Revision 1.28  2007-10-25 15:46:03  ddangelo
 * added a bit field to shaped event, meaning to be assigned. smart getters to be added at that time. internal and external event.
 * all variables zeroed in ctor for shaped event.
 *
 * Revision 1.27  2007-10-11 11:23:53  ddangelo
 * new mctruth format (internal AND root event)
 *
 * Revision 1.26  2007-05-07 13:40:26  ddangelo
 * applying patch to flag TObjects with cycle numbers
 *
 * Revision 1.25  2007-05-04 16:33:24  pallas
 * Changed variables for alpha beta
 * Now tailtot is an array of 10 values
 *
 * Revision 1.24  2007-03-27 15:18:20  ddangelo
 * variables npe_conc, charge_conc, nhits_conc added to laben cluster
 * f4_pe renamed as f4_charge in bx_laben_event.hh
 * decoded_charge and decoded_npmts added tu muon event
 *
 * Revision 1.23  2006-10-23 15:34:34  ddangelo
 * applied ale's patch to inlcude laben integrity flags
 *
 * Revision 1.22  2006/10/13 15:31:18  razeto
 * Added variables from raw data
 *
 * Revision 1.21  2006-08-21 11:08:01  razeto
 * Updated to new detector_interface
 *
 * Revision 1.20  2005/12/12 19:08:09  razeto
 * Split position and energy reco results (auth from maintainers)
 *
 * Revision 1.19  2005/12/03 15:15:54  razeto
 * Added position assigment from mctruth frame (auth from Davide)
 *
 * Revision 1.18  2005/11/18 16:42:54  razeto
 * Added new integer charge variable to decoded hit and cluster (auth from davide)
 *
 * Revision 1.17  2005/04/29 13:26:37  razeto
 * Fixed a bug, laben times need to be double precision (commit authorized by mantainer)
 *
 * Revision 1.16  2005/03/14 13:44:13  ddangelo
 * corrected to be 1-based
 *
 * Revision 1.15  2005/03/14 13:07:26  ddangelo
 * added "order_in_channel" variable in 3 levels: raw, decoded, clastered.
 * Added computation of the first one in raw event constructor
 *
 * Revision 1.14  2004/11/26 15:25:10  razeto
 * Thanks to Davide who so kindly underlined a syntactical error, I fixed the "Maintainer" word
 *
 * Revision 1.13  2004/11/26 14:06:22  razeto
 * Added Mantainer field
 *
 * Revision 1.12  2004/09/22 11:24:46  ddangelo
 * fixed 2 mispelling
 *
 * Revision 1.11  2004/09/22 11:15:26  ddangelo
 * added checks on sub-detector enabling status in all raw constructors.
 *
 * Revision 1.10  2004/08/31 13:29:03  ddangelo
 * added charge to laben hit
 *
 * Revision 1.9  2004/07/24 16:24:21  ddangelo
 * debugging (patch by Alessandro)
 *
 * Revision 1.8  2004/07/07 12:56:18  ddangelo
 * laben cluster and relative hit classes filled with meaningful variables (by Daniela)
 *
 * Revision 1.7  2004/06/07 17:17:09  ddangelo
 * b_is_valid removed from decoded hit.
 * time renamed to raw_time. getter renamed accordingly.
 *
 * Revision 1.6  2004/05/31 16:48:01  ddangelo
 * non-english name 'clusterized' replaced by 'clustered'
 *
 * Revision 1.5  2004/05/31 14:16:04  ddangelo
 * added clusterized level classes (still empty)
 *
 * Revision 1.4  2004/05/31 13:39:16  ddangelo
 * applied patch with 2 bx_laben_decoded_hit::operator=() (by Alessandro)
 *
 * Revision 1.3  2004/05/20 12:55:46  ddangelo
 * applied patch by Alessandro.
 * classes bx_laben_decoded_hit and bx_laben_decoded_event introduced.
 *
 * Revision 1.2  2004/04/27 16:58:57  ddangelo
 * minor updates
 *
 * Revision 1.1  2004/04/27 09:41:09  ddangelo
 * Event structure redisegned.
 * bx_reco_event now handles event header and has object for different detector portions.
 * Inheritance from different reconstruction status is kept within each detector portion obj.
 * Hit lists are linked by const ptr (maybe a temporary solution).
 *
 *
 */
