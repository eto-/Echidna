/* BOREXINO Reconstruction program
 *
 * Author: Alessandro Razeto <Alessandro.Razeto@ge.infn.it>
 * Maintainer: Alessandro Razeto <Alessandro.Razeto@ge.infn.it>
 *
 * $Id: bx_detector.hh,v 1.36 2015/01/09 15:03:07 misiaszek Exp $
 *
 * A class to describe some detector status
 *
 */
#ifndef _BX_DETECTOR_HIT
#define _BX_DETECTOR_HIT
#define CYCLE_NUMBER 18

#include <map>
#include <string>
#include <vector>
#include <TObject.h>
#include <RVersion.h>

class bx_detector;
class bx_base_module;

#ifndef __CINT__
#include "bx_named.hh"
class bx_echidna_event;
class bx_reader;
class detector_interface: public bx_named {
  public:
    static void init (const bx_reader&) { me = new detector_interface; }
    static bx_detector* get () { return me->detector; }
    static void post_init (const bx_echidna_event &e, const bx_reader&) { me->set_event_detector_status (e); }
  private:
    detector_interface ();
    void set_event_detector_status (const bx_echidna_event &e);
    static detector_interface* me;
    bx_detector *detector;
};
#endif

class bx_detector: public TObject {
#ifndef __CINT__
  private:	// This way db_detector can not be istantiated in Echidna
  friend class detector_interface;
#else
  public:
#endif
    bx_detector (): TObject() {}
  public:
    enum sub_detector {
      trigger,
      laben,
      muon,
      mctruth,
    };

    bool is_laben_enabled	() const { return b_laben_enabled; }
    bool is_muon_enabled	() const { return b_muon_enabled; }
    bool is_mctruth_enabled	() const { return b_mctruth_enabled; }

    const std::vector<int>& get_disabled_channels () const { return disabled_channels_v; }
    const std::vector<int>& get_disabled_charge () const { return disabled_charge_v; }
    const std::vector<int>& get_user_disabled_channels () const { return user_disabled_channels_v; }
    const std::vector<std::string>& get_disabled_laben_channel_properties () const { return disabled_laben_channel_properties_v; }
    const std::vector<std::string>& get_disabled_muon_channel_properties  () const { return disabled_muon_channel_properties_v; }

    const std::map<std::string, int>& get_skipped_events (int trg_type) { return skipped_events[trg_type]; }
    int get_total_skipped_events () const { return total_skipped_events; }
#if !defined (__CINT__) && !defined (_ECHIDNA_ROOTLIB_)
    void skip_event (bx_base_module *module, int trg_type);
    void add_disabled_channel (int lg, int evnum, int type, const bx_named*);
    void read_disabled_channels (int evnum);
#endif
  private:
    bool b_laben_enabled, b_muon_enabled, b_mctruth_enabled;
    int total_skipped_events;
    std::map<int, std::map<std::string, int> > skipped_events; // index are trigger type and module_name
    std::vector<int> disabled_channels_v, disabled_charge_v;
    std::vector<int> user_disabled_channels_v;
    std::vector<std::string> disabled_laben_channel_properties_v;
    std::vector<std::string> disabled_muon_channel_properties_v;

#if !defined (__CINT__) && !defined (_ECHIDNA_ROOTLIB_)
    bool check_laben_property  (const std::string& p);
    bool check_muon_property   (const std::string& p);
    bool check_vector_property (const std::string& p, const std::vector<std::string> &v);
    bool check_laben_property  (const std::string& p, const std::vector<std::string> &v);
    bool check_muon_property   (const std::string& p, const std::vector<std::string> &v);
    void compute_disabled_channels (int run_number);
    void disable_lg (int, const std::string&);
    void disable_charge (int, const std::string&);
    int last_evnum_read, next_evnum_to_read;
    float scaling_factor;
    bool use_db;
#endif

  ClassDef(bx_detector,CYCLE_NUMBER)
};
#endif
/*
 * $Log: bx_detector.hh,v $
 * Revision 1.36  2015/01/09 15:03:07  misiaszek
 * cycle_18 new unstable
 *
 * Revision 1.35  2013/06/18 18:56:37  razeto
 * cycle_17 new unstable
 *
 * Revision 1.34  2013-02-02 09:01:49  razeto
 * Incremented to cycle_16 (cycle 15 was lost)
 *
 * Revision 1.33  2011-04-19 05:54:58  razeto
 * Moved to cycle 15 unstable
 *
 * Revision 1.32  2011-02-28 18:29:06  razeto
 * Added interface disable_charge
 *
 * Revision 1.31  2011-02-18 17:10:05  ddangelo
 * major code cleanup: removed fadc throughout the program
 *
 * Revision 1.30  2010-08-06 17:20:16  razeto
 * Moving to cycle 14
 *
 * Revision 1.29  2009-11-26 13:42:51  razeto
 * Moved to cycle_13_unstable
 *
 * Revision 1.28  2008-12-15 17:13:55  razeto
 * New cycle (12)
 *
 * Revision 1.27  2008-10-17 13:41:12  razeto
 * new development cycle (11)
 *
 * Revision 1.26  2008-09-25 13:12:35  razeto
 * Scaling factor added for disabling channels (intended for mc)
 *
 * Revision 1.25  2008-09-25 11:52:17  razeto
 * Added db property for bx_detector::disabled_laben_channel_properties
 *
 * Revision 1.24  2008-09-24 17:14:16  razeto
 * Added reading runtime disabled channels from db
 *
 * Revision 1.23  2008-08-20 12:45:10  ddangelo
 * fixing old and new bugs...
 *
 * Revision 1.22  2008-08-19 17:53:00  ddangelo
 * added muon channel disabling
 *
 * Revision 1.21  2008-02-27 20:46:13  razeto
 * new development cycle (10)
 *
 * Revision 1.20  2008-02-27 20:26:30  razeto
 * New clasdef(9) version and new cycle version
 *
 * Revision 1.19  2007-12-18 09:51:26  ludhova
 * declaration of  check_vector_property
 *
 * Revision 1.18  2007-12-07 14:26:22  ddangelo
 * added vector for bad charge channels, getter, initialization and add() method (auth by maintainer)
 *
 * Revision 1.17  2007-11-27 15:08:30  razeto
 * Log the disable reason
 *
 * Revision 1.16  2007-10-11 10:49:54  razeto
 * Cycle 8 deployed
 *
 * Revision 1.15  2007-06-22 15:15:26  razeto
 * Moved to cycle 7
 *
 * Revision 1.14  2007-05-07 13:40:26  ddangelo
 * applying patch to flag TObjects with cycle numbers
 *
 * Revision 1.13  2007-05-02 15:58:12  razeto
 * Upgraded bx_detector to support variable disabled channel list. Few names changed
 *
 * Revision 1.12  2007-03-10 15:18:22  ddangelo
 * Now reader init the bx_detector (even with the event detector status).
 *
 * Revision 1.11  2006-11-27 10:59:01  razeto
 * Upgraded bx_detector to have more detailed check on bad channels multiplicity:
 * now bad_multiplicity do not exist any more but there are other fields.
 * Upgraded laben_decoder: the fix_bad_channel flag do not exists anymore. To
 * avoid bad channels skipping, set detector_interface.bad_laben_channel_properties
 * to {}.
 * Configurations updated.
 *
 * Revision 1.10  2006/09/10 14:41:22  razeto
 * Removed inline since creates some problems with gcc 2.95
 *
 * Revision 1.9  2006/09/09 18:59:30  razeto
 * Upgraded to handle configuration for bad channel
 *
 * Revision 1.8  2006/09/09 14:06:56  razeto
 * Upgraded to have a draft calculation of bad channels
 *
 * Revision 1.7  2006/08/21 14:08:14  razeto
 * Fixed a typo
 *
 * Revision 1.6  2006/08/21 13:51:07  razeto
 * Fixed an other root bug with template of template
 *
 * Revision 1.5  2006/08/21 11:11:59  razeto
 * Improved bx_detector and created detector_interface
 *
 * Revision 1.4  2004/11/26 15:25:10  razeto
 * Thanks to Davide who so kindly underlined a syntactical error, I fixed the "Maintainer" word
 *
 * Revision 1.3  2004/11/26 14:01:19  razeto
 * Added Mantainer field
 *
 * Revision 1.2  2004/09/22 11:18:52  razeto
 * Added mctruth to bx_detector::sub_detector
 *
 * Revision 1.1  2004/09/22 10:17:16  razeto
 * Added
 *
 */
