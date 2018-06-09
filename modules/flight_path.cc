/* BOREXINO Reconstruction program
 *
 * Author: Alessandro Razeto <Alessandro.Razeto@ge.infn.it>
 * Maintainer: Alessandro Razeto <Alessandro.Razeto@ge.infn.it>
 *
 * $Id: flight_path.cc,v 1.19 2008/10/27 09:32:43 razeto Exp $
 *
 * Implementation of flight_path
 *
 */
#include "flight_path.hh"
#include "messenger.hh"
#include "db_channel.hh"
#include "bx_dbi.hh"
#include "db_profile.hh"
#include "db_calib.hh"
#include "light_guide.hh"
#include "barn_interface.hh"
#include "bx_laben_event.hh"
#include "constants.hh"
#include "TH2.h"
#include <math.h>

int flight_path::create_histo = -1;
TH2F *flight_path::path = 0;


flight_path::flight_path (): bx_named ("flight_path"), inited(false) {
  if (create_histo < 0) {
    create_histo = get_parameter ("create_histo").get_int ();
    if (create_histo < 0) get_message (bx_message::critic) << "create_histo parameter can not be negative" << dispatch;
    if (create_histo) {
      path = new TH2F ("path", "Theta vs distance distribution", 240, 0, 12, 90, 0, 90);
      barn_interface::get ()->store (barn_interface::file, path, this);
    }
  }
}
void flight_path::set_refraction_index (float ref_index) {
  _init ();
  f4_ref_index = ref_index;
  f4_c_medium_m_ns = constants::physic::c_m_ns / f4_ref_index;
}
  
void flight_path::_init () {
  if (!inited) {  // init internal static data (singleton like)
    inited = true;
    float ref_index = get_parameter ("refidx").get_float ();
    if (ref_index <= 0) ref_index = bx_dbi::get ()->get_calib ().get_refraction_index_data ();
    set_refraction_index (ref_index);
  }
  f4_distance = f4_straight_distance = f4_theta = -1;
}

void flight_path::init (double v[3], const db_channel_laben* ch_info) {
  if (!ch_info->is_ordinary ()) get_message (bx_message::error) << "trying to evaluate a distance to a non-ordinary channel " << ch_info->get_lg () << dispatch;
    // Copy arguments locally
  b_has_cone = ch_info->pmt_has_cone ();
  f4_pmt[0] = ch_info->pmt_x (); f4_pmt[1] = ch_info->pmt_y (); f4_pmt[2] = ch_info->pmt_z ();
  f4_v[0] = v[0]; f4_v[1] = v[1]; f4_v[2] = v[2];
  
    // Init static cached variables
  _init ();
}

void flight_path::init (double x, double y, double z, const db_channel_laben* ch_info) {
  if (!ch_info->is_ordinary ()) get_message (bx_message::error) << "trying to evaluate a distance to a non-ordinary channel " << ch_info->get_lg () << dispatch;
    // Copy arguments locally
  b_has_cone = ch_info->pmt_has_cone ();
  f4_pmt[0] = ch_info->pmt_x (); f4_pmt[1] = ch_info->pmt_y (); f4_pmt[2] = ch_info->pmt_z ();
  f4_v[0] = x; f4_v[1] = y; f4_v[2] = z;
  
    // Init static cached variables
  _init ();
}

void flight_path::init (const bx_base_position& position, const db_channel_laben* ch_info) {
  if (!ch_info->is_ordinary ()) get_message (bx_message::error) << "trying to evaluate a distance to a non-ordinary channel " << ch_info->get_lg () << dispatch;
    // Copy arguments locally
  b_has_cone = ch_info->pmt_has_cone ();
  f4_pmt[0] = ch_info->pmt_x (); f4_pmt[1] = ch_info->pmt_y (); f4_pmt[2] = ch_info->pmt_z ();
  f4_v[0] = position.get_x (); f4_v[1] = position.get_y (); f4_v[2] = position.get_z ();
  
    // Init static cached variables
  _init ();
}

void flight_path::straight_distance () {
  if (f4_straight_distance < 0) { // calculate just once after init
    float dx = f4_v[0] - f4_pmt[0];
    float dy = f4_v[1] - f4_pmt[1];
    float dz = f4_v[2] - f4_pmt[2];
    f4_straight_distance = ::sqrtf (dx*dx + dy*dy + dz*dz);
  }
}

void flight_path::theta () {
  if (f4_theta < 0) {  // calculate just once after init
      // Assure straight_distance is calculated
    straight_distance ();
    if (f4_straight_distance == 0.) { f4_theta = 0; return; }

      // Check if X is outside the pmt sphere
    float x_radius2 = f4_v[0] * f4_v[0] + f4_v[1] * f4_v[1] + f4_v[2] * f4_v[2];
    float collector_radius2 = f4_pmt[0] * f4_pmt[0] + f4_pmt[1] * f4_pmt[1] + f4_pmt[2] * f4_pmt[2];
    if (x_radius2 >= collector_radius2) { f4_theta = -90; return; }
  
      // Calculate sintheta
    float ratio = f4_v[0] * f4_pmt[0] + f4_v[1] * f4_pmt[1] + f4_v[2] * f4_pmt[2];
    ratio /= collector_radius2;
  
    float d2 = 0, t;
    t = f4_v[0] - f4_pmt[0] * ratio;
    d2 += t * t;
    t = f4_v[1] - f4_pmt[1] * ratio;
    d2 += t * t;
    t = f4_v[2] - f4_pmt[2] * ratio;
    d2 += t * t;
    if (d2 == 0) { f4_theta = 0; return; }
  
    float sintheta = ::sqrtf (d2) / f4_straight_distance;
    if (sintheta > 1.) {
      get_message (bx_message::warn) << "simple_theta failed X |" << f4_v[0] << ", " << f4_v[1] << ", " << f4_v[2] << "| = " << ::sqrtf (x_radius2) << ", PMT |" << f4_pmt[0] << ", " << f4_pmt[1] << ", " << f4_pmt[2] << "| = " << ::sqrtf (collector_radius2) << ", sintheta = " << sintheta << dispatch;
      sintheta = 1;
    }
  
    f4_theta = constants::number::rad_to_deg (::asinf (sintheta));
  }
}

void flight_path::distance () { 
  if (f4_distance < 0) { // calculate just once after init

      // Assure straight_distance and theta are calculated
    straight_distance ();
    
    f4_distance = f4_straight_distance;

    if (b_has_cone) { 
      theta (); 
      f4_distance += light_guide::get ()->get_path (f4_theta);
    } else f4_distance += light_guide::get ()->get_path (f4_theta);
    
    if (create_histo) path->Fill (f4_distance, f4_theta);
  }
}

void flight_path::get_distance_derivative (float v[3]) {
    // Assure straight_distance is calculated
  straight_distance ();  // Do derivative only of straight distance (suppose light guide has low effect)

  float dx = f4_v[0] - f4_pmt[0];
  float dy = f4_v[1] - f4_pmt[1];
  float dz = f4_v[2] - f4_pmt[2];

  v[0] = dx / f4_straight_distance;
  v[1] = dy / f4_straight_distance;
  v[2] = dz / f4_straight_distance;
}
/*
 * $Log: flight_path.cc,v $
 * Revision 1.19  2008/10/27 09:32:43  razeto
 * do not error for sintheta>1
 *
 * Revision 1.18  2007-06-13 11:08:56  razeto
 * Do not rely on the sphere collector radius from db, instead calculate it pmt by pmt
 *
 * Revision 1.17  2007-05-20 12:48:11  razeto
 * Fixed initialization bug (fixes refidx problem)
 *
 * Revision 1.16  2006-09-05 13:12:04  razeto
 * Added refidx option to flight path too
 *
 * Revision 1.15  2006/08/21 11:22:06  razeto
 * Updated to new barn_interface
 *
 * Revision 1.14  2006/05/08 17:31:31  razeto
 * Added db_channel patch for fadc (sent to the mailing list)
 *
 * Revision 1.13  2006/01/02 21:23:46  razeto
 * Changed test target to file target for db_barn
 *
 * Revision 1.12  2005/10/30 11:17:50  razeto
 * Now the refraction index can be overridden from user(module) but not directly from its own parameters
 *
 * Revision 1.11  2005/10/11 15:33:19  razeto
 * Fixed a bug in derivative calculation
 *
 * Revision 1.10  2005/10/11 14:21:22  razeto
 * Upgraded to cache step results and to calculate derivative
 *
 * Revision 1.9  2005/10/06 21:27:56  razeto
 * Added a small optimization (cached pmt position)
 *
 * Revision 1.8  2005/09/27 16:07:32  razeto
 * Fixed a bug in light_guide time calculation
 *
 * Revision 1.7  2005/09/06 16:29:24  razeto
 * Added a simple control
 *
 * Revision 1.6  2005/09/06 11:58:45  razeto
 * Added refraction index parameter
 *
 * Revision 1.5  2005/06/20 14:09:39  razeto
 * Updated
 *
 * Revision 1.4  2005/05/16 12:40:00  razeto
 * Added an initializer and a simple method for getting c
 *
 * Revision 1.3  2005/03/20 09:02:54  razeto
 * Upgraded to cache some variables and to fill on demand an histogram
 *
 * Revision 1.2  2005/03/18 16:30:09  razeto
 * Upgraded to working version
 *
 * Revision 1.1  2005/03/18 11:22:49  razeto
 * Added initial version
 *
 */
