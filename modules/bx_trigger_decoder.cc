/* BOREXINO Reconstruction program
 *
 * Author: Marco Pallavicini <pallas@ge.infn.it>
 * Maintainer: Chiara Ghiano <chiara.ghiano@lngs.infn.it>
 *
 * $Id: bx_trigger_decoder.cc,v 1.22 2015/07/20 12:54:12 lukyanch Exp $
 *
 * Implementation for bx_trigger_decoder.hh
 * 
*/

#include <cmath>
#include <ctime> 
#include <iostream>
#include "bx_trigger_decoder.hh"
#include "bx_echidna_event.hh"
#include "bx_trigger_event.hh"
#include "bx_dbi.hh"
#include "db_profile.hh"
#include "barn_interface.hh"
#include <TH1F.h>

bx_trigger_decoder::bx_trigger_decoder() : bx_base_module("bx_trigger_decoder", bx_base_module::main_loop) {
}

void bx_trigger_decoder::begin () {
  u4_prev_trgid = 0;
  b_was_a_muon = false;
  i4_missing_events = 0;
  i4_missing_neutron_triggers = 0;
  i4_prev_year = 0;
  i4_prev_mon = 0;
  i4_prev_day = 0;
  i4_prev_hour = 0;
  i4_prev_min = 0;
  i4_prev_sec = 0;   
  prev_gps_times[0] = prev_gps_times[1] = 0;
  build_dt = new TH1F ("build_dt", "Time difference between gps and build time", 1000, 0, 3000);
  ppc0_dt = new TH1F ("ppc0_dt", "Time difference between gps and ppc0 time", 1000, 0, 3000);
  barn_interface::get ()->store (barn_interface::file, build_dt, this);
  barn_interface::get ()->store (barn_interface::file, ppc0_dt, this);
}

bx_echidna_event* bx_trigger_decoder::doit (bx_echidna_event *ev) {
  const bx_trigger_event& er = ev->get_trigger();
  // cast to a trigger event
  bx_trigger_decoded_event &ew = dynamic_cast<bx_trigger_decoded_event&>(ev->get_trigger());

  // dump btb threshold
  if(!ev->get_event_number())
    get_message(bx_message::log) << "BTB threshold = " << int32_t(er.get_btb_threshold()) << dispatch;

  m_decode_trgtype(er, ew);

  // call method to decode GPS data
  gpsclock_decode (er.get_gps1 (), er.get_gps2 (), er.get_gps3 (),
                  ev->get_event_number (),
	          ew.i4_day, ew.i4_mon, ew.i4_year, 
	          ew.i4_hour, ew.i4_min, ew.i4_sec, 
		  ew.i4_millisec, ew.i4_microsec,
		  ew.i4_time1, ew.i4_time2, ew.i4_time_t);

  if (! (ev->get_event_number () % 100)) { 
    bx_message &msg = get_message(bx_message::debug);
    msg << "Data: " << er.get_day () << "-" << er.get_mon () << "-" << er.get_year ();
    msg << "  Ora: " << er.get_hour () << ":" << er.get_min () << ":" << er.get_sec () << dispatch;
  }

  // check if trigger id is consecutive
  if (ev->get_event_number() > (u4_prev_trgid+1))
    i4_missing_events += ev->get_event_number() - u4_prev_trgid - 1;
  
  if (b_was_a_muon && (ev->get_trigger().get_trgtype() != bx_trigger_event::neutron))
    i4_missing_neutron_triggers +=1;

  // save evnum for check on trg 128 presence
  u4_prev_trgid = ev->get_event_number();
  b_was_a_muon = (ev->get_trigger().get_trgtype() == bx_trigger_event::neutrino) && 
 		 (ev->get_trigger().has_btb_flag(bx_trigger_event::mtb_flag));

  // Calculate dt  
  int32_t dt_gps_ppc0 = er.get_time_t() - ev-> get_trigger().get_trg_time();
  int32_t dt_build_gps = ev->get_builder_time_seconds()- er.get_time_t();

  // Selecting events with time difference between ppc0 and gps > 200 sec or <-200 sec
  if ( (ev-> get_event_number () != 0) && (dt_gps_ppc0 > 200 || dt_gps_ppc0 < -200 )){
   get_message(bx_message::error) << " for event= " << ev << " wrong gps time, fixing did not work " << er.get_time_t() << " " << ev->get_trigger().get_trg_time() << dispatch;    
  }
  
 // Fill histo  
     build_dt->Fill (dt_build_gps);
     ppc0_dt->Fill (dt_gps_ppc0 ); 
    
 // Evaluating gps time difference and  additional checks for gps dt
  uint32_t gps_times[2];                                           
  ev->get_trigger ().get_gps_time (gps_times[0], gps_times[1]);         
  int32_t gps_dt_s = gps_times[0] - prev_gps_times[0];              
  int32_t gps_dt_ns = gps_times[1] - prev_gps_times[1];                 
  double dt_gps = double(gps_dt_s) * 1e9 + double(gps_dt_ns);
  double dt_gps_sec = dt_gps * 1e-9;
	  	          
  if ((ev-> get_event_number () != 0) && (dt_gps < 0)) get_message (bx_message::error) << " negative gps_dt (" << dt_gps_sec << ") sec for event " << ev->get_event_number () <<  dispatch;
  if ((ev-> get_event_number () != 0) && (dt_gps >2e9)){ 
    get_message (bx_message::warn) << "long dt between event time and previous event time  (" << dt_gps_sec << ") sec for event " << ev->get_event_number () <<  dispatch;
  }	 
   // Saving infos for next event
   i4_prev_year = ev-> get_trigger().get_year();
   i4_prev_mon = ev-> get_trigger().get_mon();
   i4_prev_day = ev-> get_trigger().get_day();
   i4_prev_hour = ev-> get_trigger().get_hour();
   i4_prev_min = ev-> get_trigger().get_min();
   i4_prev_sec = ev-> get_trigger().get_sec(); 
   ev->get_trigger ().get_gps_time (prev_gps_times[0], prev_gps_times[1]);
   
  return ev;  
}

void bx_trigger_decoder::end () {
    get_message(bx_message::log) << "Missing events " << i4_missing_events 
    				   << " missing neutron triggers " << i4_missing_neutron_triggers << dispatch;
}

void bx_trigger_decoder::m_decode_trgtype(const bx_trigger_event& er, bx_trigger_decoded_event& ew) {
  const db_profile& profile_info = bx_dbi::get()->get_profile();
  if (er.get_trgtype() == profile_info.neutrino_trigger_tag()) ew.i4_trigger_type = bx_trigger_event::neutrino;
  else if (er.get_trgtype() == profile_info.muon_trigger_tag()) ew.i4_trigger_type = bx_trigger_event::muon;
  else if (er.get_trgtype() == profile_info.neutron_trigger_tag()) ew.i4_trigger_type = bx_trigger_event::neutron;
  else if (er.get_trgtype() == profile_info.laser266_trigger_tag()) ew.i4_trigger_type = bx_trigger_event::laser266;
  else if (er.get_trgtype() == profile_info.laser355_trigger_tag()) ew.i4_trigger_type = bx_trigger_event::laser355;
  else if (er.get_trgtype() == profile_info.laser394_trigger_tag()) ew.i4_trigger_type = bx_trigger_event::laser394;
  else if (er.get_trgtype() == profile_info.pulser_trigger_tag()) ew.i4_trigger_type = bx_trigger_event::pulser;
  else if (er.get_trgtype() == profile_info.random_trigger_tag()) ew.i4_trigger_type = bx_trigger_event::random;
}

// member function to decode GPS clock raw data info
// accepts as input 5 raw data words (coming from
// the clock )
// returns date and time in 2 formats
// standard one
// number of seconds and nano-seconds since 01-01-2000 at 0:0:0.0000000000
// (the latter is useful to compute time differences among events )
// D. Manuzio 12-11-2001
//            23-03-2004 slight changes to be class member
void bx_trigger_decoder::gpsclock_decode (
    uint32_t g1, uint32_t g2, uint32_t g3, // GPS clock raw data (input)
    int32_t evnum,     
    int32_t& mday, int32_t& mon, int32_t& year,  // day, month, year (out)
    int32_t& hour, int32_t& min, int32_t& sec,   // hour, minute, second (out)
    int32_t& msec, int32_t& microsec,        // milli-sec and micro-sec (out)
    uint32_t& sec2k, uint32_t& nsec2k,  // secs and nano-sec since 01-01-2000 at 0:0:0
    time_t& timet) {

// YEAR (*year)    current year	
// MDAY (*mday)      1 -> 31	
// MON (*mon)        1 -> 12
// HOUR (*hour)      0 -> 23
// MIN (*min)	     0 -> 59	
// SEC (*sec)	     0 -> 59	
	
  uint16_t wrd;
  int32_t first,second,third,fourth,five;
  int32_t leg;
  uint32_t nsec;
    
  // compute micro-seconds and nano-seconds
  wrd = (g1>>16)&0xFFFF;
  first = (wrd & 15);
  wrd >>= 4;
  second = (wrd & 15);
  wrd >>= 4;
  third= (wrd & 15);
  wrd >>= 4;
  fourth = (wrd & 15);
  microsec = second*1+third*10+fourth*100;
  nsec = first*100;

  // compute seconds and milli-secs
  wrd = (g1&0xFFFF);
  first = (wrd & 15);
  wrd >>= 4;
  second = (wrd & 15);
  wrd >>= 4;
  third= (wrd & 15);
  
  wrd = ((g2>>16)&0xFFFF);
  fourth = (wrd & 15);
  wrd >>= 4;
  five = (wrd & 7);
  msec = first+second*10+third*100;
  sec =fourth*1+five*10;
  if (sec<0 || sec>59)
    get_message(bx_message::warn) << "Gps Clock: Error in decoding seconds!  " << sec << dispatch;

  // compute minutes 
  wrd = ((g2>>16)&0xFFFF);
  wrd >>= 8;
  first = (wrd & 15);
  wrd >>= 4;
  second = (wrd & 7);
  min = first*1 + second*10;
  if (min<0 || min>59) 
    get_message(bx_message::warn) << "Gps Clock: Error in decoding minutes!  " << min << dispatch;
  
  // hour and solar/legal time 
  wrd = (g2&0xFFFF);
  first = (wrd & 15);
  wrd >>= 4;
  second = (wrd & 3);
  wrd >>= 2;
  third = (wrd & 1);
  hour = first*1+second*10;
  if (hour<0 || hour>23) 
    get_message(bx_message::warn) << "Gps Clock: Error in decoding hours!  " << hour << dispatch;
  leg = third;

  // month day 
  wrd = (g2&0xFFFF);
  wrd >>= 8;
  first = (wrd & 15);
  wrd >>= 4;
  second = (wrd & 3);
  mday = first*1+second*10;
  if (mday<1 || mday>31) 
    get_message(bx_message::warn) << "Gps Clock: Error in decoding days!  " << mday << dispatch;

  // month 
  wrd = g3&0xFFFF;
  first = (wrd & 15);
  wrd >>= 4;
  second = (wrd & 0x1);
  mon = first + second*10;
  if (mon<1 || mon>12) 
    get_message(bx_message::warn) << "Gps Clock: Error in decoding months!  " << mon << dispatch;
    
  // year 
  wrd = g3&0xFFFF;
  wrd >>= 8;
  first = (wrd & 15);
  wrd >>= 4;
  second = (wrd & 15);
  year = first*1+second*10;
  if (year<1 || ( year>20 && year < 92 ) ) 
    get_message(bx_message::warn) << "Gps Clock: Error in decoding year!  " << year << dispatch;
  
  // If gps3 or gps2 is wrong, pathc with previus event infos
  if ((year == 0) || ( mday == 0)) {
    year = i4_prev_year;
    mon = i4_prev_mon;
    mday = i4_prev_day;
    hour = i4_prev_hour;
    get_message (bx_message::warn) << "fixing gps time for event " << evnum << dispatch; 	
  }
  
/*  // compute number of seconds since 01-01-2000  */
  struct tm working_date;
  working_date.tm_sec = sec;
  working_date.tm_min = min;
  working_date.tm_hour = hour;
  working_date.tm_mday = mday;
  working_date.tm_mon = mon-1;
  working_date.tm_year = (year > 90) ? year : (year+100);
  working_date.tm_isdst = 0;
  
  timet = mktime(&working_date);

  // patch for gps overflow error occurred at run 16789 (2011-10-09)
  if (year > 90) {
    timet += 1024*7*86400;
    time_t tmp_timet = timet;
    if (hour == 0) tmp_timet += 3600;
    struct tm *fixed_date = gmtime(&tmp_timet);
    mday = fixed_date->tm_mday;
    mon  = fixed_date->tm_mon+1;
    year = fixed_date->tm_year;
  }

  sec2k = timet - 946684800L; // going from 1970 to 2000 based number

  if (timet > 1435708800L) sec2k += 4;      //1 Jul 2015 4 leap seconds
  else if (timet > 1341100800L) sec2k += 3; //1 Jul 2012 3 leap seconds
  else if (timet > 1230768000L) sec2k += 2; //31 Dic 2008 2 leap seconds
  else if (timet > 1136073600L) sec2k += 1; //1 Jan 2005 1 leap second

  nsec2k = (uint32_t)(msec*1000000 + microsec*1000 + nsec);
}

/*
 * $Log: bx_trigger_decoder.cc,v $
 * Revision 1.22  2015/07/20 12:54:12  lukyanch
 * 30 Jun 2015 leap second added
 *
 * Revision 1.21  2014/12/31 01:48:05  misiaszek
 * Leap second for 1st Jul 2012 added
 *
 * Revision 1.20  2011/10/14 09:33:44  ddangelo
 * patched for fixing gps overflow error
 *
 * Revision 1.19  2011-03-22 11:13:11  razeto
 * Fixed leap second start
 *
 * Revision 1.18  2011-02-22 08:11:27  razeto
 * GPS start from 2000 UTC + added leap seconds + timet
 *
 * Revision 1.17  2009-10-26 11:19:37  ddangelo
 * trgtype 'muon_tct' (or 'muon_totc') renamed as 'neutron' after real use
 * trgtype 'muon_mtb' renamed as 'muon' as duplicity with the above one disappeared and 'mtb' is associated to btb_inputs flag
 *
 * Revision 1.16  2009-10-23 13:34:31  ghiano
 * debugging
 *
 * Revision 1.15  2009-10-19 10:30:47  ghiano
 * gps decoding fixed for bad data
 *
 * Revision 1.14  2009-10-15 15:49:55  ghiano
 * gps decoding fixed for bad data
 *
 * Revision 1.13  2009-07-17 14:52:42  ddangelo
 * added histo of dt (ppc0 - gps)
 *
 * Revision 1.12  2008-05-08 19:07:24  razeto
 * Added time checking
 *
 * Revision 1.11  2008-04-14 17:04:20  ddangelo
 * added checks for trgid consecutiveness and for trg128 presence after muon
 *
 * Revision 1.10  2008-04-10 13:29:47  ddangelo
 * fixed the bug of the bisestile year
 *
 * Revision 1.9  2006-11-24 14:19:06  ddangelo
 * added a log msg with btb threshold
 *
 * Revision 1.8  2006/08/21 11:21:40  razeto
 * Updated to event format
 *
 * Revision 1.7  2004/11/29 13:21:23  razeto
 * Added Mantainer field
 *
 * Revision 1.6  2004/09/22 14:00:07  ddangelo
 * updated bx_reco_event into bx_echidna_event (on behalf of Marco)
 *
 * Revision 1.5  2004/09/22 12:26:31  ddangelo
 * updated a getter name (on behalf of Marco)
 *
 * Revision 1.4  2004/05/20 12:22:34  ddangelo
 * added a method to decode trgtype.
 * Fixed include problem.
 * commit with delegation from module coordinator
 *
 * Revision 1.3  2004/04/28 11:04:49  ddangelo
 * Modules adapted to call the right const/nonconst event portion getter
 *
 * Revision 1.2  2004/04/27 14:27:35  ddangelo
 * modifications to match new event structure.
 * (A getter is now required in event reading)
 * Cnd before dynamic cast in event writing).
 *
 * Revision 1.1  2004/04/27 10:02:48  ddangelo
 * module bx_trigger_decode_module deleted and re-added as bx_trigger_decoder
 *
 *
 */

