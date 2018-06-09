/* BOREXINO Reconstruction program
 * 
 * Author: Davide D'Angelo <davide.dangelo@lngs.infn.it>
 *
 * $Id: constants.cc,v 1.23 2012/01/16 09:48:02 ddangelo Exp $
 *
 * Implementation for constants.hh
 *
 */

#include "constants.hh"

const float constants::number::pi = 3.1415927;
const float constants::number::f4_rad_to_deg_factor = 180. / constants::number::pi;

const float constants::physic::c = 3E8;
const float constants::physic::c_m_ns = 0.3;

const unsigned char  constants::fiber::number_of_bundles = 29;
const unsigned char constants::fiber::number_of_fibers_in_bundle[] = { 0, 80, 80, 80, 80, 80, 80, 80, 80, 80, 80, 80, 80, 80, 80, 80, 80, 80, 80, 80, 80, 80, 80, 80, 80, 80, 80, 70, 70, 20 };

const unsigned char  constants::trigger::crate_number = 0;

const unsigned short constants::laben::pmts = 2212;
const unsigned short constants::laben::channels = 2240;
const unsigned char  constants::laben::ncrates = 14;
const unsigned short  constants::laben::nboards = 280;
const unsigned char  constants::laben::board_per_rack = 20;
const unsigned char  constants::laben::channels_per_rack = 160;
const unsigned char  constants::laben::channels_per_board = 8;
const unsigned char  constants::laben::frontend::board_per_rack = 14;

const unsigned char  constants::muon::tdc::max_boards = 2;
const unsigned char  constants::muon::tdc::chips_per_board = 4;
const unsigned char  constants::muon::tdc::channels_per_chip = 32;
const unsigned char  constants::muon::tdc::channels_per_board = 128;
const float          constants::muon::tdc::ns_per_clock = 1.0416;
const float          constants::muon::tdc::new_ns_per_clock = 0.8;

const unsigned short constants::muon::channels = 256;
const unsigned short constants::muon::channel_offset = 3000;
const unsigned long  constants::muon::max_hits = 600;
/* exact value for QTC gain of 1 channel QTC 6 ch 1 (no PMT taken into account: assumed 1.6 pC/pe gain) */
const float          constants::muon::pe_per_clock = 0.07295;
const unsigned char  constants::muon::crate_number = 15;

const unsigned short constants::precalibration_nevents = 1000;

/*
 * $Log: constants.cc,v $
 * Revision 1.23  2012/01/16 09:48:02  ddangelo
 * added gate width for new tdc
 *
 * Revision 1.22  2011-10-13 17:40:01  ddangelo
 * added patch to handle new OD TDCs
 *
 * Revision 1.21  2011-02-18 17:10:04  ddangelo
 * major code cleanup: removed fadc throughout the program
 *
 * Revision 1.20  2007-12-07 17:33:32  ddangelo
 * removed fixed muon tdc clock range
 *
 * Revision 1.19  2007-11-14 19:05:18  ddangelo
 * restored original gate width.
 * it should be removed and read from db
 *
 * Revision 1.18  2007-10-30 17:37:12  razeto
 * Added nboards for laben namespace
 *
 * Revision 1.17  2007-10-26 09:25:16  ddangelo
 * widening the muon tdc time buffer to 32 us to match enlarged laben trigger gate
 *
 * Revision 1.16  2007-05-04 09:09:12  ddangelo
 * added laben::pmts
 *
 * Revision 1.15  2006-09-08 10:30:43  ddangelo
 * added getters to check channel nature
 *
 * Revision 1.14  2005/09/20 16:01:09  razeto
 * Added some constants
 *
 * Revision 1.13  2005/03/18 16:35:57  razeto
 * Added some values
 *
 * Revision 1.12  2004/11/30 15:48:15  sukhotin
 * fadc::average constant was wrong. Thanks to Davide, who fixed this mistake.
 *
 * Revision 1.11  2004/11/29 16:38:24  ddangelo
 * some fadc constants added (by Sergei)
 *
 * Revision 1.8  2004/10/01 10:24:06  razeto
 * Added some laser constants
 *
 * Revision 1.7  2004/06/01 11:21:55  ddangelo
 * reorganized with nested classes
 *
 * Revision 1.6  2004/05/18 14:48:38  ddangelo
 * fixed muon tdc time buffer
 *
 * Revision 1.5  2004/05/18 14:21:59  razeto
 * Fixed a number
 *
 * Revision 1.4  2004/04/27 13:35:23  ddangelo
 * trigger and fadc sub-obj added. crate numbers added.
 *
 * Revision 1.3  2004/04/12 16:11:38  razeto
 * Added a laben constant
 *
 * Revision 1.2  2004/04/02 11:24:26  ddangelo
 * some constants added
 *
 * Revision 1.1  2004/04/01 10:50:33  ddangelo
 * added static class for wide spread constants.
 * muon raw event classes redisigned.
 *
 *
 */
