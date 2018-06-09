/* BOREXINO Reconstruction program
 *
 * Author: Alessandro Razeto <Alessandro.Razeto@ge.infn.it>
 * Maintainer: Alessandro Razeto <Alessandro.Razeto@ge.infn.it>
 *
 * $Id: db_channel.hh,v 1.36 2015/01/09 15:03:08 misiaszek Exp $
 *
 * The database interface for channel informations. READONLY
 *
 */
#ifndef _DB_CHANNEL_HH
#define _DB_CHANNEL_HH
#ifndef CYCLE_NUMBER
#define CYCLE_NUMBER 18
#endif
#include "db_profile.hh"

#include <TObject.h>

#include <stdlib.h>
#include <cmath>
#include <string>
#include <map>
#include <vector>

class db_run;
class bx_dbi;

class db_channel: public TObject {
  public:
    int get_lg 			() const { return i4_lg; }

    db_profile::channel_description_type get_channel_description () const { return channel_description; }

    bool is_ordinary 		() const { return channel_description == db_profile::ordinary; }
    bool is_trigger 		() const { return channel_description == db_profile::trigger; }
    bool is_laser 		() const { return channel_description == db_profile::laser; }
    bool is_empty 		() const { return channel_description == db_profile::empty; }


    db_channel (): TObject() {};

  protected:
    db_channel (int lg, const db_run& run, const db_profile& profile);

  private:
    int i4_lg;
    db_profile::channel_description_type channel_description;
    
  ClassDef(db_channel,CYCLE_NUMBER)
};

class db_channel_laben: public db_channel {
  public:
      // Profile interface
    bool is_cngs_id 		() const { return get_channel_description () == db_profile::cngs_id; }
    bool is_cngs_od 		() const { return get_channel_description () == db_profile::cngs_od; }
    bool is_cngs_trg 		() const { return get_channel_description () == db_profile::cngs_trg; }

      // Geometry interface
    int pmt_hole_id 		() const { return i4_pmt_hole_id; }
    bool pmt_has_cone 		() const { return b_pmt_has_cone; }
    float pmt_x			() const { return f_x; }
    float pmt_y			() const { return f_y; }
    float pmt_z			() const { return f_z; }
    float pmt_theta		() const { return f_theta; }
    float pmt_phi		() const { return f_phi; }
    int pmt_fiber_bundle	() const { return i4_pmt_fiber_bundle; }
    bool is_pmt_disconnected	() const { return b_disconnected; }

      // Channel property
    bool is_precalib_bad	() const { return b_precalib_bad; }
    bool is_precalib_off	() const { return b_precalib_off; }
    float time_offset		() const { return f_time_offset; }
    float time_sigma		() const { return f_time_sigma; }
    float charge_peak		() const { return f_charge_peak; }
    float charge_sigma		() const { return f_charge_sigma; }
    float dark_noise		() const { return f_dark_noise; }
    float dark_sigma		() const { return f_dark_sigma; }
    const std::string& pmt_status () const { return s_pmt_status; }
    const std::vector<std::string>& charge_base_status () const { return v_charge_base_status; }
    const std::vector<std::string>& charge_peak_status () const { return v_charge_peak_status; }
    const std::vector<std::string>& charge_status      () const { return v_charge_status; }
    const std::vector<std::string>& timing_status      () const { return v_timing_status; }
    const std::vector<std::string>& multiplicity       () const { return v_multiplicity; }

    db_channel_laben (): db_channel () {};
  private:
    db_channel_laben (int lg, const db_run& run, const db_profile& profile);
    bool b_pmt_has_cone;
    int i4_pmt_hole_id, i4_pmt_fiber_bundle;
    float f_x, f_y, f_z, f_theta, f_phi;
    bool b_precalib_bad, b_precalib_off;
    float f_time_offset, f_time_sigma, f_charge_peak, f_charge_sigma, f_dark_noise, f_dark_sigma;
    std::string s_pmt_status;
    std::vector<std::string> v_charge_base_status, v_charge_peak_status, v_charge_status, v_timing_status, v_multiplicity;
    bool b_disconnected;

  friend class db_channel_builder;
  ClassDef(db_channel_laben,CYCLE_NUMBER)
};

class db_channel_muon: public db_channel {
  public:
    db_channel_muon (): db_channel () {};
    int   get_hole_id      () const { return i4_hole_id; }
    float get_x	           () const { return f4_x; }
    float get_y	           () const { return f4_y; }
    float get_z	           () const { return f4_z; }
    float get_radius       () const { return ::sqrt(f4_x*f4_x+f4_y*f4_y+f4_z*f4_z); }
    float get_rc           () const { return ::sqrt(f4_x*f4_x+f4_y*f4_y); }
    float get_theta        () const { return ::acos(f4_z/get_radius()); }
    float get_phi          () const { return (f4_y > 0) ? ::acos(f4_x/get_rc()) : 2*M_PI-::acos(f4_x/get_rc()); }
    bool  is_up            () const { return get_row() > 0; }
    bool  is_down          () const { return get_row() < 0; }
    bool  is_sss           () const { return get_row() >= -2; }
    bool  is_floor         () const { return get_row() < -2; }
    int   get_row          () const { return i4_hole_id/100; }
    int   get_col          () const { return ::abs(i4_hole_id%100); }
    bool  is_disconnected  () const { return b_disconnected;  }
    float get_time_offset  () const { return f4_time_offset;  }
    float get_time_sigma   () const { return f4_time_sigma;   }
    float get_charge_peak  () const { return f4_charge_peak;  }
    float get_charge_sigma () const { return f4_charge_sigma; }
    float get_dark_noise   () const { return f4_dark_noise;   }
    float get_dark_sigma   () const { return f4_dark_sigma;   }
    const std::string& pmt_status () const { return s_pmt_status; }
    const std::vector<std::string>& multiplicity () const { return v_multiplicity; }

  private:
    int i4_hole_id;
    float f4_x, f4_y, f4_z;
    float f4_time_offset, f4_time_sigma, f4_charge_peak, f4_charge_sigma, f4_dark_noise, f4_dark_sigma;
    db_channel_muon (int lg, const db_run& run, const db_profile& profile);
    std::string s_pmt_status;
    std::vector<std::string> v_multiplicity;
    bool b_disconnected;

  friend class db_channel_builder;
  ClassDef(db_channel_muon,CYCLE_NUMBER)
};

class db_channel_builder {
  private:
    static db_channel* build (int lg, const db_run& run, const db_profile& profile);
  friend class bx_dbi;
};

#endif
/*  
 *  $Log: db_channel.hh,v $
 *  Revision 1.36  2015/01/09 15:03:08  misiaszek
 *  cycle_18 new unstable
 *
 *  Revision 1.35  2013/06/18 18:56:37  razeto
 *  cycle_17 new unstable
 *
 *  Revision 1.34  2013-02-02 09:01:49  razeto
 *  Incremented to cycle_16 (cycle 15 was lost)
 *
 *  Revision 1.33  2012-03-22 19:08:36  razeto
 *  Added cngs reference channels
 *
 *  Revision 1.32  2011-04-19 05:54:58  razeto
 *  Moved to cycle 15 unstable
 *
 *  Revision 1.31  2011-02-18 17:10:05  ddangelo
 *  major code cleanup: removed fadc throughout the program
 *
 *  Revision 1.30  2010-08-06 17:20:16  razeto
 *  Moving to cycle 14
 *
 *  Revision 1.29  2009-11-26 13:42:52  razeto
 *  Moved to cycle_13_unstable
 *
 *  Revision 1.28  2008-12-15 17:13:55  razeto
 *  New cycle (12)
 *
 *  Revision 1.27  2008-10-17 13:41:12  razeto
 *  new development cycle (11)
 *
 *  Revision 1.26  2008-08-20 14:03:42  ddangelo
 *  muon db channel completed
 *
 *  Revision 1.25  2008-08-19 17:51:31  ddangelo
 *  added muon multiplicity and pmt disconnected
 *
 *  Revision 1.24  2008-06-20 16:22:38  razeto
 *  Added an include to compile with gcc 4.3
 *
 *  Revision 1.23  2008-02-27 20:46:13  razeto
 *  new development cycle (10)
 *
 *  Revision 1.22  2008-02-27 20:26:30  razeto
 *  New clasdef(9) version and new cycle version
 *
 *  Revision 1.21  2007-11-26 14:57:09  ddangelo
 *  added getters (muon part)
 *
 *  Revision 1.20  2007-11-09 18:49:38  razeto
 *  Disconnected added to laben_channel
 *
 *  Revision 1.19  2007-11-06 14:31:02  ddangelo
 *  added a few getters
 *
 *  Revision 1.18  2007-10-11 10:49:54  razeto
 *  Cycle 8 deployed
 *
 *  Revision 1.17  2007-06-22 15:15:26  razeto
 *  Moved to cycle 7
 *
 *  Revision 1.16  2007-05-07 15:47:30  razeto
 *  Cycle number in root classdef
 *
 *  Revision 1.15  2007-03-24 15:59:12  ddangelo
 *  muon db channel improved
 *
 *  Revision 1.14  2007-03-22 19:36:58  ddangelo
 *  implemented muon db channel
 *
 *  Revision 1.13  2006/11/05 10:27:30  razeto
 *  Remove some useless include
 *
 *  Revision 1.12  2006/09/08 10:11:47  razeto
 *  Upgraded db_channels and added muon db_channel
 *
 *  Revision 1.11  2006/05/08 17:31:33  razeto
 *  Added db_channel patch for fadc (sent to the mailing list)
 *
 *  Revision 1.10  2005/11/04 00:43:29  misiaszek
 *  Inheritance from TObject added
 *
 *  Revision 1.9  2004/11/26 15:25:10  razeto
 *  Thanks to Davide who so kindly underlined a syntactical error, I fixed the "Maintainer" word
 *
 *  Revision 1.8  2004/11/26 14:01:20  razeto
 *  Added Mantainer field
 *
 *  Revision 1.7  2004/11/24 09:46:41  razeto
 *  Moved some prototyping in db_acl from db_*
 *
 *  Revision 1.6  2004/09/27 13:31:06  razeto
 *  Fibre bundle, theta and phi read from database
 *
 *  Revision 1.5  2004/07/13 10:39:21  razeto
 *  Added a test
 *
 *  Revision 1.4  2004/06/01 10:42:19  razeto
 *  Added a getter
 *
 *  Revision 1.3  2004/05/27 10:53:49  razeto
 *  Added some Davide's is_XXX getters
 *
 *  Revision 1.2  2004/05/19 11:06:27  razeto
 *  Updated to a working version of db_channel
 *
 *  Revision 1.1  2004/05/19 10:25:50  razeto
 *  Added db_channel as discussed at last cycle meeting
 *
 */
