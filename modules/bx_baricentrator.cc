/* BOREXINO Reconstruction program
 *
 * Author: Alessandro Razeto <Alessandro.Razeto@ge.infn.it>
 * Maintainer: Alessandro Razeto <Alessandro.Razeto@ge.infn.it>
 *
 * $Id: bx_baricentrator.cc,v 1.26 2011/02/18 14:12:39 ddangelo Exp $
 *
 * Implementation of bx_baricentrator
 *
 */
#include "bx_baricentrator.hh"
#include "messenger.hh"
#include "db_channel.hh"
#include "bx_echidna_event.hh"
#include "barn_interface.hh"

#include <math.h>

// ctor
bx_baricentrator::bx_baricentrator (): bx_base_module("bx_baricentrator", bx_base_module::main_loop), mean_algorithm_name_map("mean_algorithm_name_map") {
  require_event_stage (bx_detector::laben, bx_base_event::clustered);

  mean_algorithm_name_map["Q"] = Q;
  mean_algorithm_name_map["T"] = T;
  mean_algorithm_name_map["QT"] = QT;
}

// module interface
void bx_baricentrator::begin () {
  f8_first_hits_gate = get_parameter ("first_hits_gate").get_float ();
  algo = mean_algorithm_name_map[get_parameter ("mean_algorithm").get_string ()];
  f8_Q_normalization = get_parameter ("Q_normalization").get_float ();
  f8_T_normalization = get_parameter ("T_normalization").get_float ();
  f8_QT_normalization = get_parameter ("QT_normalization").get_float ();
  f8_T_filter_time = get_parameter ("T_filter_time").get_float ();
  f8_QT_filter_time = get_parameter ("QT_filter_time").get_float ();
  float ref_index = get_parameter ("refidx").get_float ();
  if (ref_index > 0) path.set_refraction_index (ref_index);

  if (algo == T) b_use_charge = false;
  else b_use_charge = true;

  baricenter_best_radius = new TH1F ("baricenter_radius", "Best error radius of the charge baricenter", 100, 0, 10);
  barn_interface::get ()->store (barn_interface::file, baricenter_best_radius, this);

}

bx_echidna_event* bx_baricentrator::doit (bx_echidna_event *ev) {
    // Loop on every cluster
  Int_t size_clusters = ev->get_laben ().get_nclusters ();
  Int_t size_mclusters = ev->get_laben ().get_nclusters_muons (); 
  for (int i = 0; i < size_clusters + size_mclusters; i++) {
    bx_laben_cluster& cluster = (i < size_clusters ) ? ev->get_laben().get_cluster(i) : ev->get_laben().get_cluster_muon(i-size_clusters);
    bx_baricenter& b = cluster.get_baricenter();

      // 1) -------------------------------------------------------------------------
      // Calculate the mean values for baricenter position and hit time distribution
      // Mean values
    int count = 0;
    double mean_x = 0., mean_y = 0., mean_z = 0., mean_t = 0.;

      // Loop on every hit
    for (int j = 0; j < cluster.get_clustered_nhits (); j++) {
        // Get clustered hit reference
      const bx_laben_clustered_hit& hit = cluster.get_clustered_hit (j);
      
        // Discard hits too late
      if (hit.get_time () > f8_first_hits_gate) continue;
      
        // Get db_channel pointer
      const db_channel_laben* ch_info = hit.get_decoded_hit ().get_db_channel ();

        // Define the multeplicity of this hit
      int pe_count = get_hit_charge (hit);

        // Accumulate values
      double t_filtered;
      switch (algo) {
	case Q:
      	  mean_x += ch_info->pmt_x () * pe_count * f8_Q_normalization;
          mean_y += ch_info->pmt_y () * pe_count * f8_Q_normalization;
          mean_z += ch_info->pmt_z () * pe_count * f8_Q_normalization;
	  break;
	case T:
	  t_filtered = hit.get_time () + f8_T_filter_time;
      	  mean_x += ch_info->pmt_x () / t_filtered * f8_T_normalization;
          mean_y += ch_info->pmt_y () / t_filtered * f8_T_normalization;
          mean_z += ch_info->pmt_z () / t_filtered * f8_T_normalization;
	  break;
	case QT:
	  t_filtered = hit.get_time () + f8_QT_filter_time;
      	  mean_x += ch_info->pmt_x () / t_filtered * pe_count * f8_QT_normalization;
          mean_y += ch_info->pmt_y () / t_filtered * pe_count * f8_QT_normalization;
          mean_z += ch_info->pmt_z () / t_filtered * pe_count * f8_QT_normalization;
	  break;
      }
      mean_t += hit.get_time () * pe_count;	// mean_t does not depend on the algorithm
      count += pe_count;
    } 
      // Ignore empty clusters (even if not a regular condition)
    if (!count) continue;
    
      // Extract the values
    mean_x /= count; mean_y /= count; mean_z /= count; mean_t /= count;

      // 2) -------------------------------------------------------------------------
      // Calculate the mean optical time of flight from the baricenter position
    double mean_tof = 0;
    
    for (int j = 0; j < cluster.get_clustered_nhits (); j++) {
      const bx_laben_clustered_hit& hit = cluster.get_clustered_hit (j);
      if (hit.get_time () > f8_first_hits_gate) continue;
    
        // Initialize path
      path.init (mean_x, mean_y, mean_z, hit.get_decoded_hit ().get_db_channel ()); 

        // Accumulate tof
      mean_tof += path.get_time () * get_hit_charge (hit);
    }

      // Extract the value (count does not changes during loops since depends only on the hit 
      // distribution and on the cuts)
    mean_tof /= count;
    

      // 3) -------------------------------------------------------------------------
      // Define start_time as the difference between mean_tof - mean_t
    double start_time = mean_tof - mean_t;
    if (start_time < -10) get_message (bx_message::log) << "negative start time " << start_time << "ns for event " << ev->get_event_number () << " cluster " << i << dispatch;
    
      // 4) -------------------------------------------------------------------------
      // Calculate the sigma of the difference from the optical tof and the measured tof for each hit
    double sigma = 0;
    
    for (int j = 0; j < cluster.get_clustered_nhits (); j++) {
      const bx_laben_clustered_hit& hit = cluster.get_clustered_hit (j);
      if (hit.get_time () > f8_first_hits_gate) continue;
      
      path.init (mean_x, mean_y, mean_z, hit.get_decoded_hit ().get_db_channel ()); 
      
      double d = path.get_time () - (hit.get_time () + start_time);
      
      sigma += d * d * get_hit_charge (hit);
    }

      // Extract the value
    sigma = ::sqrtf (sigma / count) * path.get_c_medium ();
    baricenter_best_radius->Fill (sigma);


      // 5) -------------------------------------------------------------------------
      // Fill the event data
    b.f4_t = start_time; b.f4_x = mean_x; b.f4_y = mean_y; b.f4_z = mean_z;
    b.f4_best_error_radius = sigma;
  }

  ev->get_laben ().mark_stage (bx_base_event::baricentered);
  return ev;
}

int bx_baricentrator::get_hit_charge (const bx_laben_clustered_hit& hit) {
  if (!b_use_charge) return 1;
  return hit.get_decoded_hit ().get_charge_npe ();
}

void bx_baricentrator::end () {
}

/*
 * $Log: bx_baricentrator.cc,v $
 * Revision 1.26  2011/02/18 14:12:39  ddangelo
 * added running on parallel neutron cluster in muon events.
 * to be tested.
 *
 * Revision 1.25  2008-10-19 10:05:12  razeto
 * quieter
 *
 * Revision 1.24  2007-07-08 09:53:31  razeto
 * Depend only on clustered event but neutrino trigger
 *
 * Revision 1.23  2006-08-21 11:20:23  razeto
 * Updated to new barn_interface
 *
 * Revision 1.22  2006/05/08 17:31:31  razeto
 * Added db_channel patch for fadc (sent to the mailing list)
 *
 * Revision 1.21  2006/01/02 21:23:46  razeto
 * Changed test target to file target for db_barn
 *
 * Revision 1.20  2005/11/18 16:45:44  razeto
 * Now using integer charge from decoded hit
 *
 * Revision 1.19  2005/10/30 11:18:43  razeto
 * Added refraction index overriding in flight_path
 *
 * Revision 1.18  2005/06/18 15:28:46  razeto
 * Removed Minuit minimization, now this module only calculate the baricenter
 * position with some algorithms, see docs.
 *
 * Revision 1.17  2005/03/18 16:35:08  razeto
 * Updated to use flight_path
 *
 * Revision 1.16  2005/03/18 12:05:56  razeto
 * Removed some useless histos
 *
 * Revision 1.15  2005/03/01 15:17:15  razeto
 * Merged with cycle_2
 *
 * Revision 1.14  2004/12/15 11:08:05  razeto
 * Adapted to the new event structure
 *
 * Revision 1.13.2.1  2004/12/14 16:38:44  razeto
 * Suppressed minuit warns but added echidna log messages on fail
 *
 * Revision 1.13  2004/12/03 14:44:15  razeto
 * Restored mean_radius old formula (after some discussions with Sasha)
 *
 * Revision 1.12  2004/12/03 14:11:04  razeto
 * Fixed an error and fixed mean radius RMS formula
 *
 * Revision 1.11  2004/12/03 13:36:09  razeto
 * Upgraded according to Sasha suggestions
 *
 * Revision 1.10  2004/11/26 15:25:11  razeto
 * Thanks to Davide who so kindly underlined a syntactical error, I fixed the "Maintainer" word
 *
 * Revision 1.9  2004/11/26 14:01:20  razeto
 * Added Mantainer field
 *
 * Revision 1.8  2004/10/19 16:26:52  razeto
 * Updated to use vdt::get_bool
 * Added mark_stage (still need friend declaration in laben_event).
 *
 * Revision 1.7  2004/09/27 13:31:46  razeto
 * Added require trigger type
 *
 * Revision 1.6  2004/09/22 13:27:06  razeto
 * Changed bx_reco_event to bx_echidna_event
 *
 * Revision 1.5  2004/09/22 12:21:12  razeto
 * Updated to follow new sub event getters in bx_base_even
 *
 * Revision 1.4  2004/09/22 10:49:16  razeto
 * Added an histogram; added support for results writing in the cluster
 *
 * Revision 1.3  2004/09/22 10:38:33  razeto
 * Updated to follow sub_detector enum in bx_detector
 *
 * Revision 1.2  2004/09/21 09:11:13  razeto
 * Changed roundf to floorf since roundf is unknown to the old gcc
 *
 * Revision 1.1  2004/09/17 13:35:37  razeto
 * Added bx_baricentrator
 *
 */
