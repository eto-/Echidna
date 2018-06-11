/* BOREXINO Reconstruction program
 *
 * Author: Davide D'Angelo <davide.dangelo@lngs.infn.it>
 * Maintainer: Davide D'Angelo <davide.dangelo@lngs.infn.it>
 *
 * $Id: bx_mctruth_event.cc,v 1.20 2009/07/17 17:05:41 ddangelo Exp $
 *
 * Implementation of bx_mctruth_event
 *
 */

#include "bx_mctruth_event.hh"
#include "bx_event_disk_format.h"
#include "bx_detector.hh"
#include <stdio.h>
#include <algorithm>

namespace {
  struct sort_time_operation {
    sort_time_operation () {}
    bool operator() (const bx_mctruth_hit& a, const bx_mctruth_hit& b) { return a.get_time () < b.get_time (); }
  };
};


bx_mctruth_hit::bx_mctruth_hit (uint16_t lg, float time) : u2_lg(lg), f4_time(time) {
}

bx_mctruth_daughter::bx_mctruth_daughter (mctruth_daughter_disk_format *raw) : i4_id     (raw->id), 
									       i4_pdg    (raw->pdg), 
									       f8_time   (raw->time), 
									       f4_energy (raw->energy) {
  for (int i=0; i<3; i++) {
    f4_position_v [i] = raw->position [i]/1000.;
    f4_direction_v[i] = raw->direction[i];
  }
}

bx_mctruth_deposit::bx_mctruth_deposit (mctruth_deposit_disk_format *raw) : i4_pdg_parent(raw->pdg_parent), f4_energy(raw->energy) {
  for (int i=0; i<3; i++)
    f4_position_v[i] = raw->position[i]/1000.;
}

bx_mctruth_user::bx_mctruth_user (mctruth_user_disk_format *raw) : i4_int1(raw->i1), 
								   i4_int2(raw->i2), 
								   f4_float1(raw->f1), 
								   f4_float2(raw->f2), 
								   f8_double(raw->d) {
}

bx_mctruth_frame::bx_mctruth_frame (const char *disk_event) {
  mctruth_frame_disk_format *frame = (mctruth_frame_disk_format *)disk_event;
  const char* end_of_frame = disk_event + frame->length;
  u2_file_id           = frame->file_id;
  f8_elec_event_time   = frame->elec_event_time;
  i4_event_id          = frame->event_id;
  i4_n_sequence        = frame->n_sequence;
  i4_isotope_coinc     = frame->isotope_coinc;
  i4_pdg               = frame->pdg;
  f8_time              = frame->time;
  f4_energy            = frame->energy;
  f4_visible_energy    = frame->visible_energy;
  for (int i = 0 ; i<3; i++) {
    f4_position_v  [i] = frame->position  [i];
    f4_baricenter_v[i] = frame->baricenter[i]/1000.;
    f4_direction_v [i] = frame->direction [i]/1000.;
  }
  i4_id_npe            = frame->id_npe;
  i4_od_npe            = frame->od_npe;
  i4_n_daughters       = frame->n_daughters;
  i4_n_deposits        = frame->n_deposits;
  i4_n_users           = frame->n_users;
  i4_n_id_photons      = frame->n_id_photons;
  i4_n_od_photons      = frame->n_od_photons;

  disk_event += sizeof(mctruth_frame_disk_format);

  for (int i = 0; i < i4_n_daughters; i++) {
    daughters_v.push_back( bx_mctruth_daughter((mctruth_daughter_disk_format *)disk_event));
    disk_event += sizeof(mctruth_daughter_disk_format);
  }

  for (int i = 0 ; i < i4_n_deposits; i++ ) {
    deposits_v.push_back( bx_mctruth_deposit((mctruth_deposit_disk_format *)disk_event));    
    disk_event += sizeof(mctruth_deposit_disk_format);
  }

  for (int i = 0 ; i < i4_n_users; i++ ) {
    users_v.push_back( bx_mctruth_user((mctruth_user_disk_format *)disk_event)); 
    disk_event += sizeof(mctruth_user_disk_format);
  }

  bx_message msg(bx_message::critic, "bx_mctruth_event: ");



  while (disk_event < end_of_frame) {
    mctruth_subframe_disk_format *subframe = (mctruth_subframe_disk_format *)disk_event;
    if (subframe->type == pe_id_list) {
      uint16_t n = (subframe->length - 4) / 6;
      for (int i = 0 ; i < n; i++ ) 	
        hit_id_v.push_back( bx_mctruth_hit(*(uint16_t *)(disk_event + 4 + i*2), *(float *)(disk_event + 4 + n*2 + i*4)));    
      std::sort (hit_id_v.begin (), hit_id_v.end (), sort_time_operation ());
    } else if (subframe->type == pe_od_list) {
      uint16_t n = (subframe->length - 4) / 6;
      for (int i = 0 ; i < n; i++ ) 	
        hit_od_v.push_back( bx_mctruth_hit(*(uint16_t *)(disk_event + 4 + i*2), *(float *)(disk_event + 4 + n*2 + i*4)));    
      std::sort (hit_od_v.begin (), hit_od_v.end (), sort_time_operation ());
    } else if (subframe->type == scattered_info) {
      ;
    } else msg << "unknown subframe type " << subframe->type << dispatch;
    if (!subframe->length) msg << "subframe size reported to be zero" << dispatch;
    disk_event += subframe->length;
  }
}

bx_mctruth_event::bx_mctruth_event (const char *disk_event) {
  if (!detector_interface::get ()->is_mctruth_enabled ()) return;

  mctruth_header_disk_format *header = (mctruth_header_disk_format *)disk_event;
  b_is_data = (header->length == 0xbe7bcbb0);
  if (header->length <= 0 || b_is_data) return; // the second check is against the arbitrary value used by builder for real data

  i2_trigger_jitter = header->trigger_jitter;

  disk_event += sizeof(mctruth_header_disk_format);
  for (int i = 0; i < header->frames; i++) {
    mctruth_frame_disk_format *frame = (mctruth_frame_disk_format *)disk_event;
    frame_v.push_back(disk_event);
    disk_event += frame->length;
  }
}

/*
 * $Log: bx_mctruth_event.cc,v $
 * Revision 1.20  2009/07/17 17:05:41  ddangelo
 * added a flag for real data
 *
 * Revision 1.19  2008-12-03 14:13:55  ddangelo
 * fixed a bug: file_id filled
 *
 * Revision 1.18  2008-02-18 14:24:40  ddangelo
 * changed units from mm to m for 4 variables
 *
 * Revision 1.17  2007-11-07 13:51:05  ddangelo
 * debugging
 *
 * Revision 1.16  2007-11-06 14:45:11  ddangelo
 * debugging mctruth
 *
 * Revision 1.15  2007-10-31 15:42:14  ddangelo
 * added quenched energy (mctruth)
 * disk format, internal and root event
 *
 * Revision 1.14  2007-10-29 16:28:49  ddangelo
 * reading new format.
 * subframe types for daughters, deposits and users abandoned.
 * Hierarchy assumed instead.
 * file_id added.
 * hits reintroduced, sorting added.
 *
 * Revision 1.13  2007-10-11 11:23:53  ddangelo
 * new mctruth format (internal AND root event)
 *
 * Revision 1.12  2006-08-21 11:08:16  razeto
 * Updated to new detector_interface
 *
 * Revision 1.11  2005/12/03 15:12:54  razeto
 * Added cm to meter conversion (Auth from davide)
 *
 * Revision 1.10  2005/09/23 16:50:38  razeto
 * Added mctruth hits sorting in time (auth from davide)
 *
 * Revision 1.9  2005/07/27 16:48:47  ddangelo
 * added a check on subframe length (by Ale)
 *
 * Revision 1.8  2005/05/11 14:00:29  ddangelo
 * fixed a bug
 *
 * Revision 1.7  2005/01/19 13:15:44  ddangelo
 * added elec_event_time
 *
 * Revision 1.6  2004/11/26 15:25:10  razeto
 * Thanks to Davide who so kindly underlined a syntactical error, I fixed the "Maintainer" word
 *
 * Revision 1.5  2004/11/26 14:06:22  razeto
 * Added Mantainer field
 *
 * Revision 1.4  2004/10/28 11:22:11  ddangelo
 * fixed a bug in prt arithmetics
 *
 * Revision 1.3  2004/09/22 11:24:26  ddangelo
 * added check for sub detector enabling status
 *
 * Revision 1.2  2004/07/13 13:35:47  ddangelo
 * added check to skip real data.
 *
 * Revision 1.1  2004/07/12 17:03:56  ddangelo
 * added mctruth.
 * Still undebugged, tmp commit to allow debug.
 * 
 *
 */
