/* BOREXINO Reconstruction program
 *
 * Author: Davide D'Angelo <davide.dangelo@lngs.infn.it>
 * based on work by Razeto&Pallas 
 * Maintainer: Davide D'Angelo <davide.dangelo@lngs.infn.it>
 *
 * $Id: bx_echidna_event.hh,v 1.15 2011/02/18 17:10:05 ddangelo Exp $
 *
 * This is the generic event object
 *
 */

#ifndef _BX_ECHIDNA_EVENT_HH
#define _BX_ECHIDNA_EVENT_HH

#include "bx_rec_general.hh"

#include "bx_trigger_event.hh"
#include "bx_laben_event.hh"
#include "bx_muon_event.hh"
#include "bx_mctruth_event.hh"
#include "bx_track.hh"
#include "bx_neutron_event.hh"

class bx_echidna_event {
  public:
    bx_echidna_event (const char *disk_event);
    virtual ~bx_echidna_event () {}

    uint32_t get_event_size_bytes  () const { return u4_event_size_bytes; }
    uint32_t get_run_number        () const { return u4_run_number;       }
    uint32_t get_event_number      () const { return u4_event_number;     }
    uint16_t get_enabled_crates   () const { return (uint16_t)(~u4_enabled_crates >> 16); }

    bool          is_laben_enabled (int c) const { return (c>0 && c<=constants::laben::ncrates) ? is_crate_enabled (c) : false; }
    bool          is_laben_enabled      () const { bool ret_val = false; 
                                                   for (int i=1; i<=constants::laben::ncrates; i++) 
						     ret_val |= is_crate_enabled (i); 
                                                   return ret_val;}
    bool          is_muon_enabled       () const { return is_crate_enabled (constants::muon::crate_number); }
    bool          is_mctruth_enabled    () const { return mctruth.get_nframes() > 0; }
    bool          is_neutron_enabled    () const { return neutron.is_associated(); }

    uint32_t get_builder_time_seconds     () const { return u4_builder_time_seconds;     }
    uint32_t get_builder_time_nanoseconds () const { return u4_builder_time_nanoseconds; }
    bool          is_tracked_global     () const { return track_global.is_valid(); }
    bool          is_tracked_cmt        () const { return !track_cmt.get_error(); }

    // getters for different sub-objects
    // Two functions are supplied for each sub-obj:
    // 1) The const getter should be used for reading purposes.
    // 2) The nonconst getter should be used in the module dynamic cast for writing purposes.
    const bx_trigger_event& get_trigger () const { return trigger; }
    const bx_laben_event&   get_laben   () const { return laben;   }
    const bx_muon_event&    get_muon    () const { return muon;    }
    const bx_mctruth_event& get_mctruth () const { return mctruth; }
    const bx_neutron_event& get_neutron () const { return neutron; }
    const bx_track&         get_track_global () const { return track_global;   }
    const bx_track&         get_track_cmt    () const { return track_cmt;   }
    bx_trigger_event& get_trigger () { return trigger; }
    bx_laben_event&   get_laben   () { return laben;   }
    bx_muon_event&    get_muon    () { return muon;    } 
    bx_mctruth_event& get_mctruth () { return mctruth; }
    bx_neutron_event& get_neutron () { return neutron; }
    bx_track_fitted&  get_track_global () { return track_global;   }
    bx_track_by_points&  get_track_cmt () { return track_cmt;   }

  private:
    bool is_crate_enabled (uint8_t c) const { return is_crate_enabled_daq (c) && is_crate_enabled_runtime (c); }
    bool is_crate_enabled_daq (uint8_t c) const { if (c==0) return true; 
                                                        if (c <=16) return u4_enabled_crates & (1 << (c-1)); 
						        else return false; }
    bool is_crate_enabled_runtime (uint8_t c) const { if (c==0) return true; 
                                                            if (c <=16) return ! (u4_enabled_crates & (1 << (15+c))); 
						            else return false; }
    void  m_dump_crates_status();

    uint32_t u4_event_size_bytes;
    uint32_t u4_run_number;
    uint32_t u4_event_number;
    uint32_t u4_enabled_crates; //0-15 bits: daq (0=OFF, 1=ON); 16-31 bits: runtime (0=ON, 1=OFF); AAA: INVERTED LOGIC!!!
    uint32_t u4_builder_time_seconds;
    uint32_t u4_builder_time_nanoseconds;

    bx_trigger_event trigger;
    bx_laben_event   laben;
    bx_muon_event    muon;
    bx_mctruth_event mctruth;
    bx_neutron_event neutron;
    bx_track_fitted  track_global;
    bx_track_by_points track_cmt;
};

#endif

/*
 * $Log: bx_echidna_event.hh,v $
 * Revision 1.15  2011/02/18 17:10:05  ddangelo
 * major code cleanup: removed fadc throughout the program
 *
 * Revision 1.14  2010-05-21 15:55:14  ddangelo
 * different things on muon tracks
 *
 * 1.a) old laben_track renamed as laben_track_energy
 * new laben_track_tof added
 *
 * 1.b) (global) track renemed as track_global at base event level
 * track_cmt added at base event level (track by points)
 *
 * 1) all getters updated/integrated
 * is_tracked variable updated/integrated accordingly. inizialization.
 * job ported to root event as well. copy done.
 * friendship with old/new module updated
 *
 * 2) bxtrack_by_points class:
 * - theta, phi and impact added as variables.
 * - errors added on all of the above.
 * - error code variable requested by cmt tracker added
 *
 * Revision 1.13  2009-07-31 15:39:50  ddangelo
 * debugging the work of the [previous commit
 *
 * Revision 1.12  2009-07-30 16:20:44  ddangelo
 * - removed IsTracked() getters from BxTrack classes (root event)
 * + added is_tracked variable and getter in BxEvent and BxLaben classes (root event)
 * + fixed the copy of is_tracked variable in BxMuon class (root event)
 * + added is_tracked variable and getter for bx_echidna event class (internal event). to be filled by global tracker
 *
 * Revision 1.11  2008-08-26 15:30:30  ddangelo
 * added neutron enabled and association flag, muon aligned flags
 *
 * Revision 1.10  2008-07-11 17:06:28  ddangelo
 * added classes for neutron system (code by S. Davini)
 *
 * Revision 1.9  2008-05-09 15:46:23  ddangelo
 * muon track implemented as "by_points" and "fitted" classes inheriting from common interface (pure virtual class bx_track)
 * laben and muon host an istance "by_points"
 * echidna event hosts an istance "fitted"
 * structure duplicated at root level and filled
 * GetDistance duplicated in BxPosition class
 *
 * Revision 1.8  2007-06-01 15:57:00  ddangelo
 * added the enabled crates word to root event.
 * One plain getter in echidna event
 *
 * Revision 1.7  2007-05-30 16:02:48  ddangelo
 * muon has_cluster variable converted to int. -1 for mcr disabled condition.
 * runtime enabled condition of crates recovered from builder logic
 *
 * Revision 1.6  2007-03-15 19:17:19  ddangelo
 * pid event removed.
 * laben event upgraded with classes: bx_laben_shaped_cluster and bx_laben_ab_cluster
 * bx_laben_rec_cluster is now a parent class for the 2 new ones.
 * BxEvent modified accordingly: BxLabenRecHit and BxLabenRecCluster added.
 * BxPidEvent removed.
 *
 * Revision 1.5  2005/07/11 17:11:45  ddangelo
 * removed global event
 * added pid event
 *
 * Revision 1.4  2005/03/14 18:13:11  ddangelo
 * added a getter for mctruth existence
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
 * Revision 1.16  2004/09/22 12:13:15  ddangelo
 * removed '_nonconst' from sub event getters
 *
 * Revision 1.15  2004/07/12 17:03:56  ddangelo
 * added mctruth.
 * Still undebugged, tmp commit to allow debug.
 *
 * Revision 1.14  2004/07/12 10:28:14  ddangelo
 * added bx_global_event
 *
 * Revision 1.13  2004/06/01 11:38:23  ddangelo
 * updated to new constants syntax
 *
 * Revision 1.12  2004/05/30 12:30:59  ddangelo
 * debugging
 *
 * Revision 1.11  2004/05/27 14:47:45  ddangelo
 * added support for fadc event
 *
 * Revision 1.10  2004/05/18 16:45:44  ddangelo
 * some low level getters pushed private.
 * getters for crate enabled status implemented with daq rather than runtime (empty).
 *
 * Revision 1.9  2004/04/28 11:03:38  ddangelo
 * implemented 2 (const/nonconst) getters for sub-objects, for reading/writing within modules
 *
 * Revision 1.8  2004/04/27 13:35:48  ddangelo
 * debugging
 *
 * Revision 1.7  2004/04/27 13:24:06  ddangelo
 * a few more fancy getters
 *
 * Revision 1.6  2004/04/27 12:13:08  ddangelo
 * added getters for enabled_crates and builder times
 *
 * Revision 1.5  2004/04/27 09:41:09  ddangelo
 * Event structure redisegned.
 * bx_echidna_event now handles event header and has object for different detector portions.
 * Inheritance from different reconstruction status is kept within each detector portion obj.
 * Hit lists are linked by const ptr (maybe a temporary solution).
 *
 * Revision 1.4  2004/03/23 15:31:57  razeto
 * Added bx_decoded_event inheritance to bx_echidna_event
 *
 * Revision 1.3  2004/03/21 18:27:19  razeto
 * Fixed a typo in the log tag
 *
 *
 */
