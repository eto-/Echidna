/* BOREXINO Reconstruction program
 *
 * Author: Marco Pallavicini <pallas@ge.infn.it>
 * Maintainer: Davide D'Angelo <davide.dangelo@lngs.infn.it>
 *
 * $Id: bx_trigger_event.hh,v 1.20 2011/02/22 08:11:27 razeto Exp $
 *
 * The trigger event object
 * 
 */
#ifndef _BX_TRIGGER_EVENT_HH
#define _BX_TRIGGER_EVENT_HH

#include <time.h>
#include "bx_rec_general.hh"
#include "bx_base_event.hh"
#include "constants.hh"



class bx_trigger_raw_event {
  public:
    bx_trigger_raw_event (const char *disk_trigger_event);

    enum btb_flag {
      mtb_flag     = 4,  // BTB input #3 muon MTB
      neutron_flag = 8,  // BTB input #4 thought for muon analog, then used for neutron in old trg setup 
      l355_flag    = 16, // BTB input #5 laser355 
      l266_flag    = 32, // BTB input #6 laser266 
      lcr_flag     = 64  // BTB input #7 common to laser394, calibration and random
    };

    uint16_t get_length          () const { return u2_length; }
    uint16_t get_error           () const { return u2_error; }    
    uint16_t get_blocks          () const { return u2_blocks; }
    uint16_t get_version         () const { return u2_version; } 
    uint32_t  get_dsp_registers   () const { return u4_dsp_registers; }
    uint16_t get_btb_threshold   () const { return u2_btb_threshold; }
    uint16_t get_btb_firmware    () const { return u2_btb_firmware; }
    uint32_t  get_trg_time        () const { return u4_trg_time; }
    uint32_t  get_evid            () const { return u4_evid; }
    uint16_t get_tab_sum         () const { return u2_tab_sum; }
    uint8_t  get_btb_inputs      () const { return u1_btb_inputs; }
    bool           has_btb_flag        ( btb_flag flag ) const { return u1_btb_inputs & flag; }
    uint8_t  get_trgtype         () const { return u1_trgtype; }
    uint32_t  get_gps1            () const { return u4_gps1; }
    uint32_t  get_gps2            () const { return u4_gps2; }
    uint32_t  get_gps3            () const { return u4_gps3; }

  private:
    uint16_t u2_length;
    uint16_t u2_error;
    uint16_t u2_blocks;
    uint16_t u2_version;
    uint32_t  u4_dsp_registers;
    uint16_t u2_btb_firmware;
    uint16_t u2_btb_threshold;
    uint32_t  u4_trg_time;
    uint32_t  u4_evid;
    uint8_t  u1_trgtype;
    uint8_t  u1_btb_inputs;
    uint16_t u2_tab_sum;
    uint32_t  u4_trgw;
    uint32_t  u4_gps1;
    uint32_t  u4_gps2;
    uint32_t  u4_gps3;
};

class bx_trigger_decoded_event {
  public:
    bx_trigger_decoded_event () {}

    enum trigger_type {
      neutrino,
      muon,
      neutron,
      laser266,
      laser355,
      laser394,
      pulser,
      random,
    };

    trigger_type get_trigger_type () const { return i4_trigger_type; }
    bool is_neutrino () const { return i4_trigger_type == neutrino; }
    bool is_muon     () const { return i4_trigger_type == muon;     }
    bool is_neutron  () const { return i4_trigger_type == neutron;  }
    bool is_laser266 () const { return i4_trigger_type == laser266; }
    bool is_laser355 () const { return i4_trigger_type == laser355; }
    bool is_laser394 () const { return i4_trigger_type == laser394; }
    bool is_laser    () const { return i4_trigger_type == laser266 || i4_trigger_type == laser355 || i4_trigger_type == laser394; } 
    bool is_pulser   () const { return i4_trigger_type == pulser;   }
    bool is_random   () const { return i4_trigger_type == random;   }
    bool is_service  () const { return is_laser() || is_pulser() || is_random();  }

    // getters
    int get_day      () const { return i4_day; }
    int get_mon      () const { return i4_mon; }
    int get_year     () const { return i4_year; }
    int get_hour     () const { return i4_hour; }
    int get_min      () const { return i4_min; }
    int get_sec      () const { return i4_sec; }
    int get_millisec () const { return i4_millisec; }
    int get_microsec () const { return i4_microsec; }
    time_t get_time_t () const { return i4_time_t; }

    // return event time in seconds, and microseconds since Jan 01, 2000 0:0:0;
    void get_gps_time    (uint32_t &secs, uint32_t &nano) const { secs=i4_time1; nano=i4_time2; }

  private:
    int i4_day, i4_mon, i4_year;
    int i4_hour, i4_min, i4_sec, i4_millisec, i4_microsec;
    uint32_t i4_time1, i4_time2;  // absolute secs and microsecs since 01/01/2000 UTC
    time_t i4_time_t;

    trigger_type i4_trigger_type;

  friend class bx_trigger_decoder;
};

class bx_trigger_event: 
  public bx_base_event,
  public bx_trigger_raw_event, 
  public bx_trigger_decoded_event {
  public:
    bx_trigger_event (const char *disk_event): bx_trigger_raw_event(disk_event) {}
    virtual ~bx_trigger_event () {}

  private:
    void mark_stage(event_stage es, bool v = true) { stages[es] = v; }

  friend class bx_trigger_decoder;
};

#endif
/*
 * $Log: bx_trigger_event.hh,v $
 * Revision 1.20  2011/02/22 08:11:27  razeto
 * GPS start from 2000 UTC + added leap seconds + timet
 *
 * Revision 1.19  2009-10-26 11:19:44  ddangelo
 * trgtype 'muon_tct' (or 'muon_totc') renamed as 'neutron' after real use
 * trgtype 'muon_mtb' renamed as 'muon' as duplicity with the above one disappeared and 'mtb' is associated to btb_inputs flag
 *
 * Revision 1.18  2008-06-20 16:22:38  razeto
 * Added an include to compile with gcc 4.3
 *
 * Revision 1.17  2008-05-08 19:06:36  razeto
 * Added get_time_t and GetTimeT from gps time
 *
 * Revision 1.16  2007-05-25 14:37:35  ddangelo
 * added management of btb flags.
 * Redesigned internal low level access.
 *
 * Revision 1.15  2006-12-14 15:13:22  ddangelo
 * added 8 bit masking for btb threshold and tab sum.
 *
 * Revision 1.14  2006/11/24 14:17:14  ddangelo
 * added a masking in btb thresh getter
 *
 * Revision 1.13  2006/08/21 11:09:25  razeto
 * Upgrades
 *
 * Revision 1.12  2006/07/17 18:16:56  ddangelo
 * added a getter
 *
 * Revision 1.11  2005/12/18 12:02:20  razeto
 * Fixed a bad name (auth from davide)
 *
 * Revision 1.10  2004/11/26 15:25:10  razeto
 * Thanks to Davide who so kindly underlined a syntactical error, I fixed the "Maintainer" word
 *
 * Revision 1.9  2004/11/26 14:10:38  razeto
 * Added Mantainer field
 *
 * Revision 1.8  2004/09/22 11:11:40  ddangelo
 * added enum with trgtypes (by Ale)
 *
 * Revision 1.7  2004/05/25 21:46:10  ddangelo
 * fixed a name duplication bug
 *
 * Revision 1.6  2004/05/25 16:37:05  ddangelo
 * added support for event stage flagging.
 * Inheritance from bx_base_event introduced
 *
 * Revision 1.5  2004/05/20 15:26:38  ddangelo
 * fixed a byte ordering problem
 *
 * Revision 1.4  2004/05/20 10:25:00  ddangelo
 * bool locals for trigger types introduced in decoded event, along with relative getters.
 *
 * Revision 1.3  2004/05/18 17:29:33  ddangelo
 * trigger disk event format expanded. class bx_trigger_raw_event expanded accordingly.
 * Getters expanded also. A bitfield added for btb inputs (8 bit only, should be portable).
 *
 * Revision 1.2  2004/04/27 10:45:23  ddangelo
 * fixed module name in friend declaration
 *
 * Revision 1.1  2004/04/27 09:41:09  ddangelo
 * Event structure redisegned.
 * bx_reco_event now handles event header and has object for different detector portions.
 * Inheritance from different reconstruction status is kept within each detector portion obj.
 * Hit lists are linked by const ptr (maybe a temporary solution).
 *
 *
 */
