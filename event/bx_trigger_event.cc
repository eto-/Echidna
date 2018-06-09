/* BOREXINO Reconstruction program
 *
 * Author: Marco Pallavicini <pallas@ge.infn.it>
 * Maintainer: Davide D'Angelo <davide.dangelo@lngs.infn.it>
 *
 * $Id: bx_trigger_event.cc,v 1.10 2007/05/25 14:37:35 ddangelo Exp $
 *
 * Implementation of bx_trigger_event
 *
 */

#include "bx_trigger_event.hh"
#include "bx_event_disk_format.h"
#include <iostream>
#include <stdio.h>

bx_trigger_raw_event::bx_trigger_raw_event (const char *disk_event) {

  trigger_disk_format *head = (trigger_disk_format *)disk_event;

  u2_length = head->length;
  u2_error = head->error;
  u2_blocks = head->blocks;
  u2_version = head->version;
  u4_dsp_registers = head->dsp_registers;
  u2_btb_threshold = head->btb_threshold & 0xff; // upper 8 bits to be ignored
  u2_btb_firmware = head->btb_firmware;
  u4_trg_time = head->trg_time;
  u4_evid = head->evid;
  u1_trgtype = head->trgtype;
  u1_btb_inputs = head->btb_inputs; // casting into the bitfield through ptr casting. Can it be done better?
  u2_tab_sum = head->tab_sum && 0xff; // upper 8 bits to be ignored
  u4_gps1 = head->gps1;
  u4_gps2 = head->gps2;
  u4_gps3 = head->gps3;

}

/*
 * $Log: bx_trigger_event.cc,v $
 * Revision 1.10  2007/05/25 14:37:35  ddangelo
 * added management of btb flags.
 * Redesigned internal low level access.
 *
 * Revision 1.9  2007-02-19 17:07:17  ddangelo
 * debugging
 *
 * Revision 1.8  2006/12/14 15:13:22  ddangelo
 * added 8 bit masking for btb threshold and tab sum.
 *
 * Revision 1.7  2004/11/26 15:25:10  razeto
 * Thanks to Davide who so kindly underlined a syntactical error, I fixed the "Maintainer" word
 *
 * Revision 1.6  2004/11/26 14:10:38  razeto
 * Added Mantainer field
 *
 * Revision 1.5  2004/05/25 21:46:10  ddangelo
 * fixed a name duplication bug
 *
 * Revision 1.4  2004/05/20 15:29:20  ddangelo
 * removed 2 unnecessary includes
 *
 * Revision 1.3  2004/05/20 15:26:38  ddangelo
 * fixed a byte ordering problem
 *
 * Revision 1.2  2004/05/18 17:29:33  ddangelo
 * trigger disk event format expanded. class bx_trigger_raw_event expanded accordingly.
 * Getters expanded also. A bitfield added for btb inputs (8 bit only, should be portable).
 *
 * Revision 1.1  2004/04/27 09:41:09  ddangelo
 * Event structure redisegned.
 * bx_reco_event now handles event header and has object for different detector portions.
 * Inheritance from different reconstruction status is kept within each detector portion obj.
 * Hit lists are linked by const ptr (maybe a temporary solution).
 * 
 *
 */
