/* BOREXINO Reconstruction program
 *
 * Author: Alessandro Razeto <Alessandro.Razeto@ge.infn.it>
 * Maintainer: Alessandro Razeto <Alessandro.Razeto@ge.infn.it>
 *
 * $Id: bx_position_reco.cc,v 1.14 2009/10/23 14:00:04 koshio Exp $
 *
 * Implementation of bx_position_reco
 *
 */
#include "bx_position_reco.hh"
#include "messenger.hh"
#include "db_channel.hh"
#include "bx_echidna_event.hh"
#include "bx_mctruth_event.hh"
#include "barn_interface.hh"
#include <algorithm>

#include <math.h>

// ctor
bx_position_reco::bx_position_reco (): bx_base_module("bx_position_reco", bx_base_module::main_loop), position_algorithm_name_map("position_algorithm_name_map") {
  require_event_stage (bx_detector::laben, bx_base_event::clustered);

  position_algorithm_name_map["baricenter"] = baricenter;
  position_algorithm_name_map["milano"] = milano;
  position_algorithm_name_map["lngs"] = lngs;
  position_algorithm_name_map["moscow"] = moscow;
  position_algorithm_name_map["dubna"] = dubna;
  position_algorithm_name_map["mach4"] = mach4;
  position_algorithm_name_map["mctruth"] = mctruth;
  position_algorithm_name_map["fixed"] = fixed;
}

// module interface
void bx_position_reco::begin () {
  algo = position_algorithm_name_map[get_parameter ("position_algorithm").get_string ()];
  float ref_index = get_parameter ("refidx").get_float ();
  if (ref_index > 0) path.set_refraction_index (ref_index);
  for (int i = 0; i < 3 && algo == fixed; i++)
    fixed_positions[i] = get_parameter("fixed_position").get_vector ()[i].get_float ();

}

bx_echidna_event* bx_position_reco::doit (bx_echidna_event *ev) {
    // Loop on every cluster
  for (int i = 0; i < ev->get_laben ().get_nclusters (); i++) {
      // Get cluster reference
    const bx_laben_cluster& cluster = ev->get_laben ().get_cluster (i);
    bx_laben_rec_event &laben_rec = dynamic_cast<bx_laben_rec_event &>(ev->get_laben ());
    laben_rec.rec_clusters.push_back (cluster);
    bx_laben_rec_cluster &rec_cluster = laben_rec.get_rec_cluster(i);
    bx_position &position = rec_cluster.position;

    switch (algo) {
      case baricenter:
	if (ev->get_laben ().check_stage (bx_base_event::baricentered)) position = cluster.get_baricenter ();
	else 
	  get_message (bx_message::error) << "baricenter position not available for event " << ev->get_event_number () << " cluster " << i << dispatch;
	break;
      case milano:
	if (ev->get_laben ().check_stage (bx_base_event::reconstructed_mi)) position = cluster.get_position_mi ();
	else get_message (bx_message::error) << "milano position not available for event " << ev->get_event_number () << " cluster " << i << dispatch;
	break;
      case lngs:
	if (ev->get_laben ().check_stage (bx_base_event::reconstructed_lngs)) position = cluster.get_position_lngs ();
	else get_message (bx_message::error) << "lngs position not available for event " << ev->get_event_number () << " cluster " << i << dispatch;
	break;
      case moscow:
	if (ev->get_laben ().check_stage (bx_base_event::reconstructed_msk)) position = cluster.get_position_msk ();
	else get_message (bx_message::error) << "moscow position not available for event " << ev->get_event_number () << " cluster " << i << dispatch;
	break;
      case dubna:
	if (ev->get_laben ().check_stage (bx_base_event::reconstructed_dbn)) position = cluster.get_position_dbn ();
	else get_message (bx_message::error) << "dubna position not available for event " << ev->get_event_number () << " cluster " << i << dispatch;
	break;
      case mach4:
	if (ev->get_laben ().check_stage (bx_base_event::reconstructed_mach4)) position = cluster.get_position_mach4 ();
	else get_message (bx_message::error) << "mach4 position not available for event " << ev->get_event_number () << " cluster " << i << dispatch;
	break;
      case mctruth:
	if (ev->is_mctruth_enabled ()) position = ev->get_mctruth ().get_frame(i); 
	else get_message (bx_message::error) << "mctruth data not available for event " << ev->get_event_number () << " cluster " << i << dispatch;
	break;
      case fixed:
	position = bx_base_position(fixed_positions[0], fixed_positions[1], fixed_positions[2]);
	break;
    }; 
    tof_hits (rec_cluster);
  }

  ev->get_laben ().mark_stage (bx_base_event::reconstructed);
  return ev;
}

void bx_position_reco::tof_hits (bx_laben_rec_cluster &rec_cluster) {
  const bx_laben_cluster& cluster = *rec_cluster.p_cluster;

  for (int j = 0; j < cluster.get_clustered_nhits (); j++) {
    const bx_laben_clustered_hit& hit = cluster.get_clustered_hit (j);
    
      // Initialize path
    path.init (rec_cluster.get_position (), hit.get_decoded_hit ().get_db_channel ()); 

      // Calculate tof
    float tof = path.get_time ();
    
      // Create new_hit
    rec_cluster.rec_hits.push_back (hit);  
    rec_cluster.rec_hits[j].f8_time = hit.get_time () - tof;
  }
  std::sort (rec_cluster.rec_hits.begin (), rec_cluster.rec_hits.end ());

  double t0 = rec_cluster.rec_hits[0].f8_time;
  for (int j = 0; j < rec_cluster.get_rec_nhits (); j++) rec_cluster.rec_hits[j].f8_time -= t0;
}
  
void bx_position_reco::end () {
}

/*
 * $Log: bx_position_reco.cc,v $
 * Revision 1.14  2009/10/23 14:00:04  koshio
 * Add the lngs postion reconstruction
 *
 * Revision 1.13  2009-07-16 15:53:00  razeto
 * Added mach4 position reco
 *
 * Revision 1.12  2009-07-16 10:54:01  ddangelo
 * matched a modified variable name in laben event
 *
 * Revision 1.11  2008-12-15 11:46:55  razeto
 * Added fixed position (for sources)
 *
 * Revision 1.10  2007-07-08 09:53:31  razeto
 * Depend only on clustered event but neutrino trigger
 *
 * Revision 1.9  2007-05-10 16:40:42  razeto
 * T0 subtraction fixed
 *
 * Revision 1.8  2007-05-04 16:52:40  razeto
 * Rec hits start from zero
 *
 * Revision 1.7  2007-03-15 19:53:58  ddangelo
 * casting to match last event modifications
 *
 * Revision 1.6  2006/08/21 11:17:25  razeto
 * Updated to new barn_interface
 *
 * Revision 1.5  2006/01/25 13:12:00  misiaszek
 * Added a missing include (auth from maintainer)
 *
 * Revision 1.4  2005/12/30 09:37:29  razeto
 * Added generation of hits with time of flight correction
 *
 * Revision 1.3  2005/12/03 15:18:57  razeto
 * Added mctruth algorith (which just copy mc data to bx_position)
 *
 * Revision 1.2  2005/10/30 11:18:43  razeto
 * Added refraction index overriding in flight_path
 *
 * Revision 1.1  2005/06/20 14:17:26  razeto
 * Added position reco
 *
 *
 */
