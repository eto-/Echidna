/* BOREXINO Reconstruction program
 *
 * Author: Davide D'Angelo <davide.dangelo@lngs.infn.it>
 * Maintainer: Davide D'Angelo <davide.dangelo@lngs.infn.it>
 *
 * $Id: bx_muon_event.cc,v 1.29 2012/01/16 13:04:55 ddangelo Exp $
 *
 * Implementation of bx_muon_event.hh
 *
 */

#include "bx_muon_event.hh"
#include "bx_event_disk_format.h"
#include "messenger.hh"
#include "bx_detector.hh"

#define NEW_TDC_CODE 42 // the answer

// old TDC hit ctr
bx_muon_raw_hit::bx_muon_raw_hit (unsigned short muon_channel, const bx_muon_edge& first_edge, const bx_muon_edge& second_edge) : u2_muon_channel(muon_channel), u2_lead_time(first_edge.time), u2_trail_time(second_edge.time) {}

// new TDC hit ctr
bx_muon_raw_hit::bx_muon_raw_hit (unsigned short muon_channel, unsigned short lead_time, unsigned short trail_time) : u2_muon_channel(muon_channel), u2_lead_time(lead_time), u2_trail_time(trail_time) {}


// Performs preliminary checks on the edge couple. Old TDC.
bx_muon_raw_event::check_t bx_muon_raw_event::m_check_validity (const bx_muon_edge& el, const bx_muon_edge& et) const {
  if (el.slope  == et.slope || el.tdc_chip != et.tdc_chip || el.tdc_channel != et.tdc_channel )
    return missing_edge;
  if (el.time  < et.time )
    return inverted;
  return fine;
}

  
bx_muon_raw_event::bx_muon_raw_event (const char *disk_event) {
  if (!detector_interface::get ()->is_muon_enabled ()) return;

  muon_header_disk_format *head = (muon_header_disk_format *)disk_event;
  if (head->nwords == 0) return;

  u4_nedges = head->nwords - sizeof (muon_header_disk_format) / sizeof (unsigned long);
  u4_trgid = head->trgid;
  u4_error_flag = head->error_flag;

  raw_hits.clear ();

  if (u4_error_flag != NEW_TDC_CODE) { // we use the unused error word as a flag to understand which TDC are mounted.
  // old TDC unpacking
	  unsigned char chip = 0;
	  for (unsigned long i = 0; i < u4_nedges; i++) {
		  unsigned short muon_channel = chip * constants::muon::tdc::channels_per_chip; // still without edge channel 
		  bx_muon_edge *el = (bx_muon_edge*)(disk_event + sizeof (muon_header_disk_format) + i * sizeof (muon_edge_disk_format));
		  if ( el->is_last ) chip++;
		  if ( el->is_event_number || el->is_invalid ) continue;

		  bx_muon_edge *et = (bx_muon_edge*)(disk_event + sizeof (muon_header_disk_format) + (++i) * sizeof (muon_edge_disk_format));
		  if ( et->is_last ) chip++;
		  if ( et->is_event_number || et->is_invalid ) continue;
		  switch (m_check_validity(*el, *et)) {
			  case fine:
			  	muon_channel += el->tdc_channel; // now add edge channel
				if (muon_channel > 255) {
					  bx_message msg(bx_message::error, "bx_muon_raw_event::bx_muon_raw_event() ");
					  msg << "Invalid muon channel " << muon_channel << dispatch;
				}
				else raw_hits.push_back (bx_muon_raw_hit(muon_channel, *el, *et));
				break;
			  case missing_edge:
				  /*bx_message msg(bx_message::debug, "bx_muon_raw_event: ");
				    msg << "event " << u4_trgid << "; mch " << (board * constants::muon::tdc::chips_per_board + el->tdc_chip) * constants::muon::tdc::channels_per_chip + el->tdc_channel 
				    << "; missing edge" << dispatch;*/
				  i--;
				  if ( et->is_last ) chip--; // to correct double counting    
				  break;
			  case inverted:
				  /*bx_message msg(bx_message::debug, "bx_muon_raw_event: ");
				    msg << "event " << u4_trgid << "; mch " << (board * constants::muon::tdc::chips_per_board + el->tdc_chip) * constants::muon::tdc::channels_per_chip + el->tdc_channel 
				    << "; inverted time ordering; lead " << el->time << "; trail " << et->time << ";\n" << dispatch;*/
				  break;
		  } // end of switch    
	  } // end of loop on edge pairs
	  if (chip != (constants::muon::tdc::max_boards * constants::muon::tdc::chips_per_board)) {
		  bx_message msg(bx_message::error, "bx_muon_raw_event: ");
		  msg << "Trgid " << u4_trgid << ". Inconsitent muon data. Found tdc chips #" << int(chip) << dispatch;
	  }
  } // end of is old TDC
  else {
    // new TDC unpacking
	std::vector<unsigned short> edge_times         [constants::muon::channels];
	std::vector<bool>           edge_trailing_flags[constants::muon::channels];
	for (unsigned long i = 0; i < u4_nedges; i++) {
		unsigned long *edge_ptr= (unsigned long*)(disk_event + sizeof (muon_header_disk_format) + i * sizeof (muon_edge_disk_format));
 		unsigned short time         =  *edge_ptr      & 0xffff; //16 bits
		unsigned char  tdc_channel  = (*edge_ptr>>19) & 0x7f  ; // 7 bits
		bool           is_trailing  = (*edge_ptr>>26) & 1     ;
		bool           is_board_1   = (*edge_ptr>>28) & 1     ;
		bool           is_board_2   = (*edge_ptr>>29) & 1     ;
		bool           is_invalid   = (*edge_ptr>>31) & 1     ;
		// AAA: TDC boards have been inverted during replacement on 13/01/2012. To correct this we cross them here. 
		int 	       mch = is_board_1 * constants::muon::tdc::channels_per_board + tdc_channel; 

		if (is_invalid) {
			bx_message msg(bx_message::warn, "bx_muon_raw_event: ");
		  	msg << "Trgid=" << u4_trgid << " OD TDC edge=" << i << " on mch=" << mch << " is invalid." << dispatch;	
			continue;
		}

		// 2) Check edge info autoconsistency
		if ( (is_board_1 == is_board_2) ) { 
			bx_message msg(bx_message::critic, "bx_muon_raw_event: ");
			msg << "event " << u4_trgid << " OD TDC edge=" << i << " on mch=" << mch << " corrupted data format. B1=" << is_board_1 << " B2=" << is_board_2 << dispatch;
		}

		edge_times         [mch].push_back(time);
		edge_trailing_flags[mch].push_back(is_trailing);
		
	} // end of loop on edges

        for (int iCh = 0; iCh < constants::muon::channels; iCh++) {
		unsigned short starting_ied = 0;
		unsigned short nedges = edge_trailing_flags[iCh].size();
		if (!nedges) continue;
		if (edge_trailing_flags[iCh][0]) { // first edge is trailing: caught half hit staring before the gate, skip it.
//			bx_message msg(bx_message::debug, "bx_muon_raw_event: ");
//			msg << "event " << u4_trgid << " mch=" << iCh << " starting with trailing edge. skipping it." << dispatch;
			starting_ied=1;	
		}
/*		if (!edge_trailing_flags[iCh][nedges-1]) { // last edge is leading: caught half hit ending after the gate, will be skipped automatically
			bx_message msg(bx_message::debug, "bx_muon_raw_event: ");
			msg << "event " << u4_trgid << " mch=" << iCh << " ending with leading edge. will be skipped." << dispatch;	
		}*/
		unsigned short prev_time = 0;
		for (unsigned short iEd = starting_ied; iEd < edge_times[iCh].size(); iEd++) {
			unsigned short time = edge_times[iCh][iEd];
			bool is_trailing = edge_trailing_flags[iCh][iEd];
			if (time < prev_time) {
				bx_message msg(bx_message::error, "bx_muon_raw_event: ");
				msg << "event " << u4_trgid << " mch=" << iCh << " invalid time ordering. prev=" << prev_time << " now=" << time << dispatch;		
				prev_time = time;
				continue;
			}
			if (is_trailing) {
				unsigned short length = time - prev_time;
				raw_hits.push_back(bx_muon_raw_hit(iCh, prev_time, time));
				if (length > 3000) { // QTC should not deliver signals longer then ~2us
					bx_message msg(bx_message::warn, "bx_muon_raw_event: ");
					msg << "event " << u4_trgid << " mch=" << iCh << " iEd=" << iEd << " extra long hit " << length << " ticks" << dispatch;		
				}
			} // end of is_trailing
			prev_time = time;
		} // end of loop in edges in the same channel
	} // end of loop on channels
  } // end of is new TDC
}

bx_muon_decoded_event::bx_muon_decoded_event () : i4_decoded_npmts(0), 
						  f4_decoded_charge(0.),
						  b_is_aligned(false) { 
  nhits_per_channel_dec = new int[constants::muon::channels]; 
  std::fill_n(nhits_per_channel_dec, constants::muon::channels, 0); 
}

bx_muon_clustered_event::bx_muon_clustered_event (): b_has_cluster_sss(false),
						     b_has_cluster_floor(false),
						     i4_npmts(0), 
						     i4_nhits_sss(0), 
						     i4_nhits_floor(0), 
						     f4_start_time_sss(0.), 
						     f4_start_time_floor(0.), 
						     f4_charge_sss(0.),
						     f4_charge_floor(0.) {
  nhits_per_channel = new int[constants::muon::channels]; 
  std::fill_n(nhits_per_channel, constants::muon::channels, 0); 
}

/*
 * $Log: bx_muon_event.cc,v $
 * Revision 1.29  2012/01/16 13:04:55  ddangelo
 * cleaned up and commented
 *
 * Revision 1.28  2012-01-16 09:48:38  ddangelo
 * decoding upgraded to handle new tdc format
 *
 * Revision 1.27  2012-01-14 17:08:20  ddangelo
 * muon raw event construction for new TDC format
 *
 * Revision 1.26  2011-10-13 17:40:16  ddangelo
 * added patch to handle new OD TDCs
 *
 * Revision 1.25  2009-10-26 11:23:16  ddangelo
 * removed to debug level messages
 *
 * Revision 1.24  2008-12-11 17:42:34  ddangelo
 * is_aligned flag initialized to false
 *
 * Revision 1.23  2008-10-22 10:29:51  ddangelo
 * filling raw hits vector only if channel is reasonable. discards junk data
 *
 * Revision 1.22  2008-08-26 13:39:39  ddangelo
 * added is_aligned flag
 *
 * Revision 1.21  2008-02-21 15:58:31  ddangelo
 * added muon track as a separate class.
 * Currently used twice: laben and muon (later also globally fit).
 * both inner and roor event
 * added tracked stage in base event enum
 * bx_muon_rec_event renamed as bx_muon_tracked_event
 *
 * Revision 1.20  2007-12-20 18:41:29  ddangelo
 * handling of individual nhits for sss and floor
 * filling some varaibles left over before.
 * some more debugging
 *
 * Revision 1.19  2007-11-29 15:00:51  ddangelo
 * debugging
 *
 * Revision 1.18  2007-11-26 14:07:37  ddangelo
 * clustered and rec level modified
 *
 * Revision 1.17  2007-11-14 17:07:55  ddangelo
 * new muon clustering variables.
 * indipendent up/floor clustered hits vectors
 * (internal and root event)
 * filling and inizialization. tested.
 *
 * Revision 1.16  2007-10-12 16:32:13  ddangelo
 * added rec stage
 *
 * Revision 1.15  2007-10-11 11:23:53  ddangelo
 * new mctruth format (internal AND root event)
 *
 * Revision 1.14  2007-07-11 13:22:50  ddangelo
 * lowering the error level for failed tdc #chips in data from critic to error.
 *
 * Revision 1.13  2007-05-30 16:02:49  ddangelo
 * muon has_cluster variable converted to int. -1 for mcr disabled condition.
 * runtime enabled condition of crates recovered from builder logic
 *
 * Revision 1.12  2007-04-02 10:47:32  ddangelo
 * debugging
 *
 * Revision 1.11  2006-08-21 11:08:01  razeto
 * Updated to new detector_interface
 *
 * Revision 1.10  2006/07/01 14:06:08  ddangelo
 * test on tdc chip number is now skipped in case error word is set.
 *
 * Revision 1.9  2004/11/26 15:25:10  razeto
 * Thanks to Davide who so kindly underlined a syntactical error, I fixed the "Maintainer" word
 *
 * Revision 1.8  2004/11/26 14:06:22  razeto
 * Added Mantainer field
 *
 * Revision 1.7  2004/09/22 11:24:46  ddangelo
 * fixed 2 mispelling
 *
 * Revision 1.6  2004/09/22 11:15:26  ddangelo
 * added checks on sub-detector enabling status in all raw constructors.
 *
 * Revision 1.5  2004/06/07 12:47:34  ddangelo
 * fixed some messages
 *
 * Revision 1.4  2004/06/01 11:38:23  ddangelo
 * updated to new constants syntax
 *
 * Revision 1.3  2004/05/20 18:17:26  ddangelo
 * low level data containers simplified.
 * class bx_muon_edge changed to a bitfield.
 * removed static masks and shifts for old data parsing.
 * raw hit construction directly into hit vector.
 * check_validity moved from hit to event class.
 *
 * Revision 1.1  2004/04/27 09:41:09  ddangelo
 * Event structure redisegned.
 * bx_reco_event now handles event header and has object for different detector portions.
 * Inheritance from different reconstruction status is kept within each detector portion obj.
 * Hit lists are linked by const ptr (maybe a temporary solution).
 *
 *
 */
