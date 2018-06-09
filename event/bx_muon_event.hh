/* BOREXINO Reconstruction program
 *
 * Author: Davide D'Angelo <davide.dangelo@lngs.infn.it>
 * Maintainer: Davide D'Angelo <davide.dangelo@lngs.infn.it>
 *
 * $Id: bx_muon_event.hh,v 1.37 2013/10/03 14:55:07 re Exp $
 *
 * The muon event object
 * 
 */
#ifndef _BX_MUON_EVENT_HH
#define _BX_MUON_EVENT_HH

#include <cstdlib>
#include <vector>

#include "bx_rec_general.hh"
#include "bx_base_event.hh"
#include "constants.hh"
#include "bx_track.hh"

// ****** raw *************/
// bitfield to parse tdc word (32 bits)
struct bx_muon_edge {
  unsigned short time : 16;
  unsigned char slope : 1;
  bool is_overflow : 1;
  unsigned char tdc_channel : 5; 
  bool is_event_number : 1;
  unsigned char tdc_chip : 2;
  int : 4; // this bits are reserved by CAEN in the tdc word
  bool is_last : 1;
  bool is_invalid : 1;
};

class bx_muon_raw_hit {
  public:

    bx_muon_raw_hit (unsigned short muon_channel, const bx_muon_edge& lead, const bx_muon_edge& trail);
    bx_muon_raw_hit (unsigned short muon_channel, unsigned short lead_time, unsigned short trail_time);

    unsigned short get_muon_channel    () const { return u2_muon_channel; } // 0 based
    unsigned short get_logical_channel () const { return get_muon_channel () + constants::muon::channel_offset + 1; } // 1 based
    unsigned short get_lead_time       () const { return u2_lead_time; }
    unsigned short get_trail_time      () const { return u2_trail_time; }
    int            get_time_diff       () const { return abs(u2_lead_time - u2_trail_time); }

  private:
    unsigned short u2_muon_channel;
    unsigned short u2_lead_time;
    unsigned short u2_trail_time;
};
   

class bx_muon_raw_event {
  public:
    typedef std::vector<bx_muon_raw_hit> bx_muon_raw_hit_vector;

    bx_muon_raw_event (const char *disk_event);
    virtual ~bx_muon_raw_event () {}
    const bx_muon_raw_hit& get_raw_hit (int i) const { return raw_hits[i]; }
    int get_raw_nhits              ()      const {return raw_hits.size();}
    const bx_muon_raw_hit_vector&  get_raw_hits    ()      const {return raw_hits;}
    unsigned long get_nedges       ()      const {return u4_nedges;}
    unsigned long get_trgid        ()      const {return u4_trgid;}
    unsigned long get_error_flag   ()      const {return u4_error_flag;}

  private:
    enum check_t {fine, missing_edge, inverted};
    check_t m_check_validity  (const bx_muon_edge&, const bx_muon_edge&) const;

    unsigned long u4_nedges;
    unsigned long u4_trgid;
    unsigned long u4_error_flag;
    bx_muon_raw_hit_vector raw_hits;
};

/********** decoded *************/
class db_channel_muon;
class bx_muon_decoded_hit {
  public:
    bx_muon_decoded_hit ( float time, float charge, const bx_muon_raw_hit* raw_hit_ptr) : f4_time(time), f4_charge(charge), raw_hit(raw_hit_ptr) {}
    const bx_muon_raw_hit& get_raw_hit () const {return *raw_hit;}

    float get_time          () const { return f4_time; }
    float get_charge        () const { return f4_charge; }
    const db_channel_muon* get_db_channel () const { return p_db_channel; }
  
  private:
    float f4_time, f4_charge;
    const bx_muon_raw_hit* raw_hit;
    const db_channel_muon* p_db_channel;

  friend class bx_muon_decoder;
}; 

class bx_muon_decoded_event {
  public:
    bx_muon_decoded_event ();
    virtual ~bx_muon_decoded_event () { delete nhits_per_channel_dec; }

    const bx_muon_decoded_hit& get_decoded_hit (int i) const { return decoded_hits[i];    }
    int   get_decoded_nhits  ()            const { return decoded_hits.size();            }
    int   get_decoded_nhits  (int channel) const { return nhits_per_channel_dec[channel]; }
    float get_decoded_charge ()            const { return f4_decoded_charge;              }	    
    int   get_decoded_npmts  ()            const { return i4_decoded_npmts;               }
    bool  is_aligned         ()            const { return b_is_aligned;                   }

  private:
    int* nhits_per_channel_dec;
    int i4_decoded_npmts;
    float f4_decoded_charge;
    bool b_is_aligned;
    typedef std::vector<bx_muon_decoded_hit> bx_muon_decoded_hit_vector;
    bx_muon_decoded_hit_vector decoded_hits;

  friend class bx_muon_decoder;
};

/************ clustered ***************/
class bx_muon_clustered_hit {
  public:
    bx_muon_clustered_hit ( float time, float charge, const bx_muon_decoded_hit* decoded_hit_ptr) :  i4_affiliation(0), f4_time(time), f4_charge(charge), decoded_hit(decoded_hit_ptr) {}
    const bx_muon_decoded_hit& get_decoded_hit () const {return *decoded_hit;}

    int   get_affiliation   () const { return i4_affiliation; }    
    float get_time          () const { return f4_time       ; }
    float get_charge        () const { return f4_charge     ; }

  private:
    int   i4_affiliation;
    float f4_time, f4_charge;
    const bx_muon_decoded_hit* decoded_hit;

  public:
    typedef std::vector<bx_muon_clustered_hit> bx_muon_clustered_hit_vector;

  friend class bx_muon_findcluster;
}; 

class bx_muon_cluster {
  public:
    bx_muon_cluster (int id, float x, float y, float z, float c, float t) : i4_id(id), f4_x(x), f4_y(y), f4_z(z), f4_charge(c), f4_start_time(t) {}
    int   get_id        () const { return i4_id; }
    float get_charge    () const { return f4_charge; }
    float get_x         () const { return f4_x; }
    float get_y         () const { return f4_y; }
    float get_z         () const { return f4_z; }
    float get_radius    () const { return sqrt(f4_x*f4_x+f4_y*f4_y+f4_z*f4_z); }
    float get_rc        () const { return sqrt(f4_x*f4_x+f4_y*f4_y); }
    float get_theta     () const { return ::acos(f4_z/get_radius()); }
    float get_phi       () const { return (f4_y > 0) ? ::acos(f4_x/get_rc()) : 2*constants::number::pi-::acos(f4_x/get_rc()); }
    float get_start_time() const { return f4_start_time; }
    bool  is_up         () const { return f4_z > 0.; }
    bool  is_down       () const { return (f4_z < 0.) && (f4_z > -5.); }
    bool  is_sss        () const { return f4_z > -5.; }
    bool  is_floor      () const { return f4_z <= -5.; }

  private:
    int   i4_id;
    float f4_x,f4_y,f4_z;
    float f4_charge;
    float f4_start_time; 
  public:
    typedef std::vector<bx_muon_cluster> bx_muon_cluster_vector;

  friend class bx_muon_findcluster;
};

class bx_muon_clustered_event {
  public:
    bx_muon_clustered_event ();
    virtual ~bx_muon_clustered_event () { delete nhits_per_channel; }

    const bx_muon_cluster& get_cluster (int i) const { return clusters[i];   }
    int   get_nclusters                  () const { return clusters.size(); }
    const bx_muon_clustered_hit& get_clustered_hit (int i) const { return clustered_hits[i];   }
    int   get_clustered_nhits            () const { return clustered_hits.size(); }
    int   get_clustered_nhits (int channel) const { return nhits_per_channel[channel]; }
    bool  has_cluster_sss                () const { return b_has_cluster_sss; }
    bool  has_cluster_floor              () const { return b_has_cluster_floor; }
    bool  has_cluster                    () const { return b_has_cluster_sss || b_has_cluster_floor; }
    float get_start_time_sss             () const { return f4_start_time_sss; }
    float get_start_time_floor           () const { return f4_start_time_floor; }
    float get_start_time                 () const { return (f4_start_time_sss < f4_start_time_floor) ? 
    								f4_start_time_sss : f4_start_time_floor; }
    float get_charge_sss                 () const { return f4_charge_sss;    }	    
    float get_charge_floor               () const { return f4_charge_floor; }    
    float get_charge                     () const { return f4_charge_sss+f4_charge_floor; }	    
    int   get_npmts                      () const { return i4_npmts; }
    int   get_clustered_nhits_sss        () const { return i4_nhits_sss; }
    int   get_clustered_nhits_floor      () const { return i4_nhits_floor; }

  private:
    int*  nhits_per_channel;
    bool  b_has_cluster_sss;
    bool  b_has_cluster_floor;
    int   i4_npmts;
    int   i4_nhits_sss;
    int   i4_nhits_floor;
    float f4_start_time_sss;
    float f4_start_time_floor;
    float f4_charge_sss, f4_charge_floor;
    bx_muon_cluster::bx_muon_cluster_vector clusters;
    bx_muon_clustered_hit::bx_muon_clustered_hit_vector clustered_hits;

  friend class bx_muon_findcluster;
};

class bx_muon_tracked_event {
  public:
    bx_muon_tracked_event () : b_is_tracked(false) {}
    virtual ~bx_muon_tracked_event () { }

    bool is_tracked () const { return b_is_tracked; }

    const bx_track& get_track () const { return track; }
//    const bx_track_by_points& get_track () const { return track; }
    bx_track_by_points& get_track () { return track; }

  private:
    bool b_is_tracked;
    bx_track_by_points track; 

  friend class bx_muon_tracker;
};

class bx_muon_event: 
  public bx_base_event,
  public bx_muon_raw_event, 
  public bx_muon_decoded_event, 
  public bx_muon_clustered_event,
  public bx_muon_tracked_event {
  public:
    bx_muon_event (const char *disk_event): bx_muon_raw_event(disk_event) { if (get_raw_nhits ()) mark_stage (bx_base_event::raw); }
    virtual ~bx_muon_event () {}

  private:
    void mark_stage(event_stage es, bool v = true) { stages[es] = v; }

  friend class bx_muon_decoder;
  friend class bx_muon_findcluster;
  friend class bx_muon_tracker;
};

#endif
/*
 * $Log: bx_muon_event.hh,v $
 * Revision 1.37  2013/10/03 14:55:07  re
 * Add #include <cstdlib> so that Echidna compiles in SL6
 *
 * Revision 1.36  2012-01-16 09:48:38  ddangelo
 * decoding upgraded to handle new tdc format
 *
 * Revision 1.35  2012-01-14 17:08:20  ddangelo
 * muon raw event construction for new TDC format
 *
 * Revision 1.34  2011-10-13 17:40:16  ddangelo
 * added patch to handle new OD TDCs
 *
 * Revision 1.33  2009-07-30 15:01:49  ddangelo
 * b_is_tracked initialized
 *
 * Revision 1.32  2008-08-26 13:39:39  ddangelo
 * added is_aligned flag
 *
 * Revision 1.31  2008-05-09 15:46:23  ddangelo
 * muon track implemented as "by_points" and "fitted" classes inheriting from common interface (pure virtual class bx_track)
 * laben and muon host an istance "by_points"
 * echidna event hosts an istance "fitted"
 * structure duplicated at root level and filled
 * GetDistance duplicated in BxPosition class
 *
 * Revision 1.30  2008-02-21 15:58:31  ddangelo
 * added muon track as a separate class.
 * Currently used twice: laben and muon (later also globally fit).
 * both inner and roor event
 * added tracked stage in base event enum
 * bx_muon_rec_event renamed as bx_muon_tracked_event
 *
 * Revision 1.29  2007-12-20 18:41:29  ddangelo
 * handling of individual nhits for sss and floor
 * filling some varaibles left over before.
 * some more debugging
 *
 * Revision 1.28  2007-12-06 16:47:12  ddangelo
 * better math
 *
 * Revision 1.27  2007-11-29 15:00:51  ddangelo
 * debugging
 *
 * Revision 1.26  2007-11-26 15:32:43  ddangelo
 * better use of constants
 *
 * Revision 1.25  2007-11-26 14:07:37  ddangelo
 * clustered and rec level modified
 *
 * Revision 1.24  2007-11-14 19:00:57  ddangelo
 * added plain getters for entry/exit point tracked positions
 *
 * Revision 1.23  2007-11-14 17:07:55  ddangelo
 * new muon clustering variables.
 * indipendent up/floor clustered hits vectors
 * (internal and root event)
 * filling and inizialization. tested.
 *
 * Revision 1.22  2007-10-12 16:32:13  ddangelo
 * added rec stage
 *
 * Revision 1.21  2007-05-30 16:02:49  ddangelo
 * muon has_cluster variable converted to int. -1 for mcr disabled condition.
 * runtime enabled condition of crates recovered from builder logic
 *
 * Revision 1.20  2007-04-02 10:47:32  ddangelo
 * debugging
 *
 * Revision 1.19  2007-03-27 15:18:20  ddangelo
 * variables npe_conc, charge_conc, nhits_conc added to laben cluster
 * f4_pe renamed as f4_charge in bx_laben_event.hh
 * decoded_charge and decoded_npmts added tu muon event
 *
 * Revision 1.18  2007-03-22 16:08:43  ddangelo
 * added some stuff for clustered level
 *
 * Revision 1.17  2007-02-21 18:49:44  ddangelo
 * fixed a few types
 *
 * Revision 1.16  2007/02/21 16:02:27  ddangelo
 * added reconstructed stage
 *
 * Revision 1.15  2006/09/08 11:11:47  ddangelo
 * modified to match db_channel specialization
 *
 * Revision 1.14  2006/08/24 18:03:40  ddangelo
 * changed a getter to return a ref rather then ptr.
 *
 * Revision 1.13  2004/11/26 15:25:10  razeto
 * Thanks to Davide who so kindly underlined a syntactical error, I fixed the "Maintainer" word
 *
 * Revision 1.12  2004/11/26 14:06:22  razeto
 * Added Mantainer field
 *
 * Revision 1.11  2004/07/25 13:31:14  ddangelo
 * added mark stage for raw event
 *
 * Revision 1.10  2004/06/01 11:38:23  ddangelo
 * updated to new constants syntax
 *
 * Revision 1.9  2004/05/27 14:54:14  ddangelo
 *  removed useless friend declarations
 *
 * Revision 1.8  2004/05/27 14:52:30  ddangelo
 * p_db_channel changed to const
 *
 * Revision 1.7  2004/05/25 16:37:05  ddangelo
 * added support for event stage flagging.
 * Inheritance from bx_base_event introduced
 *
 * Revision 1.6  2004/05/21 11:11:51  ddangelo
 * some names updated for consistency
 *
 * Revision 1.5  2004/05/20 18:32:26  ddangelo
 * added support for db_channel
 *
 * Revision 1.4  2004/05/20 18:17:26  ddangelo
 * low level data containers simplified.
 * class bx_muon_edge changed to a bitfield.
 * removed static masks and shifts for old data parsing.
 * raw hit construction directly into hit vector.
 * check_validity moved from hit to event class.
 *
 * Revision 1.2  2004/05/20 11:12:15  ddangelo
 * added a getter for the hit vector (required by algorythm use)
 *
 * Revision 1.1  2004/04/27 09:41:09  ddangelo
 * Event structure redisegned.
 * bx_reco_event now handles event header and has object for different detector portions.
 * Inheritance from different reconstruction status is kept within each detector portion obj.
 * Hit lists are linked by const ptr (maybe a temporary solution).
 *
 *
 */
