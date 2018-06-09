/* BOREXINO Reconstruction program
 *
 * Author: Davide D'Angelo <davide.dangelo@lngs.infn.it>
 * based on work by Razeto&Pallas 
 * Maintainer: Davide D'Angelo <davide.dangelo@lngs.infn.it>
 *
 * $Id: bx_echidna_event.cc,v 1.6 2011/02/18 17:10:05 ddangelo Exp $
 *
 */

#include "bx_echidna_event.hh"
#include "bx_event_disk_format.h"
#include "messenger.hh"

bx_echidna_event::bx_echidna_event (const char *disk_event): 
  trigger(disk_event + sizeof (event_header_disk_format)), 
  laben  (disk_event + ((event_header_disk_format *)disk_event)->laben_offset),
  muon   (disk_event + ((event_header_disk_format *)disk_event)->muon_offset),
  mctruth(disk_event + ((event_header_disk_format *)disk_event)->mctruth_offset) {
	
  event_header_disk_format *header = (event_header_disk_format *)disk_event;
	
  u4_event_size_bytes         = header->event_size_bytes;
  u4_run_number               = header->run_number;
  u4_event_number             = header->event_number;
  u4_enabled_crates           = header->enabled_crates;
  u4_builder_time_seconds     = header->builder_time_seconds;
  u4_builder_time_nanoseconds = header->builder_time_nanoseconds;

//  if (!(u4_event_number%100)) m_dump_crates_status();

}

// debugging tool, currently unused. It will be possibly removed in future
void bx_echidna_event::m_dump_crates_status () {
  bx_message msg(bx_message::info, "bx_echidna_event ");
  msg << " Enabled crates (event " << u4_event_number << "; word " << u4_enabled_crates << " ): " << std::endl;
  msg << "laben crates daq: ";
  for (int i = 1; i <= constants::laben::ncrates; i++)
    msg << i << (is_crate_enabled_daq(i) ? ":ON " : ":OFF "); 
  msg << std::endl;
  msg << "laben crates runtime: ";
  for (int i = 1; i <= constants::laben::ncrates; i++)
    msg << i << (is_crate_enabled_runtime(i) ? ":ON " : ":OFF "); 
  msg << std::endl;
  msg << "muon crate daq: " << (is_crate_enabled_daq(15) ? "ON " : "OFF ") 
      << " runtime: " << ( is_crate_enabled_runtime(15) ? "ON " : "OFF ")  << std::endl;
  msg << dispatch;
}

/*
 * $Log: bx_echidna_event.cc,v $
 * Revision 1.6  2011/02/18 17:10:05  ddangelo
 * major code cleanup: removed fadc throughout the program
 *
 * Revision 1.5  2007-05-30 16:02:48  ddangelo
 * muon has_cluster variable converted to int. -1 for mcr disabled condition.
 * runtime enabled condition of crates recovered from builder logic
 *
 * Revision 1.4  2005-07-11 17:11:45  ddangelo
 * removed global event
 * added pid event
 *
 * Revision 1.3  2004/11/26 15:25:10  razeto
 * Thanks to Davide who so kindly underlined a syntactical error, I fixed the "Maintainer" word
 *
 * Revision 1.2  2004/11/26 14:06:22  razeto
 * Added Mantainer field
 *
 * Revision 1.1  2004/09/24 09:11:17  razeto
 * Updated bx_reco_event into bx_echidna_event
 *
 * Revision 1.7  2004/07/12 17:03:56  ddangelo
 * added mctruth.
 * Still undebugged, tmp commit to allow debug.
 *
 * Revision 1.6  2004/06/01 11:38:23  ddangelo
 * updated to new constants syntax
 *
 * Revision 1.5  2004/05/27 14:47:45  ddangelo
 * added support for fadc event
 *
 * Revision 1.4  2004/05/18 16:42:39  ddangelo
 * fixed a bug in initialization. A debug dump private method added.
 *
 * Revision 1.3  2004/04/28 11:03:38  ddangelo
 * implemented 2 (const/nonconst) getters for sub-objects, for reading/writing within modules
 *
 * Revision 1.2  2004/04/27 12:13:08  ddangelo
 * added getters for enabled_crates and builder times
 *
 * Revision 1.1  2004/04/27 09:41:09  ddangelo
 * Event structure redisegned.
 * bx_echidna_event now handles event header and has object for different detector portions.
 * Inheritance from different reconstruction status is kept within each detector portion obj.
 * Hit lists are linked by const ptr (maybe a temporary solution).
 *
 *
 */
