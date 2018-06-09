/* BOREXINO Reconstruction program
 *
 * Author: Alessandro Razeto <Alessandro.Razeto@ge.infn.it>
 * 	   Marco Pallavicini <pallas@ge.infn.it>
 * Maintainer: Alessandro Razeto <Alessandro.Razeto@ge.infn.it>
 *
 * $Id: bx_base_module.hh,v 1.30 2011/02/18 17:10:05 ddangelo Exp $
 *
 * This is the pure base class for any module
 * A module is a piece of code performing any action on the data
 * See class bx_reco_framework for details
 * Each module has a pure virtual interface made of three methods:
 *    begin()    the operations the module should do BEFORE the event loop
 *    doit()     operation on each event
 *    end()      operations performed after the end of the event loop
 * Even if templates are possible for [sg]et_param_* they are not used to avoid
 * complicated code.
 *
 */
#ifndef _BX_BASE_MODULE_H
#define _BX_BASE_MODULE_H

#include "bx_rec_general.hh"
#include "vdt.hh"
#include "messenger.hh"
#include "cmap.hh"
#include "bx_base_event.hh"
#include "bx_trigger_event.hh"
#include "bx_named.hh"
#include "bx_detector.hh"

#include <string>
#include <vector>
#include <iostream>
#include <iomanip>

class bx_echidna_event;

class bx_base_module: public bx_named {
  public:
    enum module_role {
      reader,
      precalib_cycle1,
      precalib_cycle2,
      precalib_cycle3,
      precalib_cycle4,
      main_loop,
      writer,
      all,
    };
    bx_base_module (const std::string& myname, module_role role);
    virtual ~bx_base_module () {}
    
      // Get the module role, status, priority
    module_role get_role () const { return i_role; }
    bool is_enabled () 	const { return b_enabled; }
    int get_priority () const { return get_parameter ("priority").get_int (); } 

      // Override bx_named::set_parameter to keep track of b_enabled
    virtual void set_parameter (const std::string& name, const vdt& value);
    
      // Some operators
    bool operator< (const bx_base_module& b) const { 
      if (get_role () == b.get_role ()) return get_priority () <= b.get_priority ();
      else return get_role () < b.get_role ();
    }
    bool operator== (const std::string& name) const { return name == get_name (); }
    bool operator== (module_role role) const { return is_enabled () ? role == get_role () : false; }

      // Operations for the framework (needed to make some job in the 
      // bx_base_module infrastructure before/after the actual call to 
      // the real module method {begin|doit|end})
      // The role argument is the role this operation is required to, to
      // be compared to the current module role.
    void framework_begin (module_role role);
    bx_echidna_event* framework_doit (bx_echidna_event *ev, module_role role);
    void framework_end (module_role role);
    
      // Operation exported by the modules, internally called by the framework_ variant
    virtual void begin () = 0;
    virtual bx_echidna_event* doit (bx_echidna_event *ev) = 0;
    virtual void end () = 0;

    typedef std::vector<bx_base_module *> bx_base_module_vector;

  protected:
    bool b_has_data; 	// A variable intended to be used by modules which do operation at the end:
    			// bx_base_module initialize it to false; usage is left to inherited modules

    void require_event_stage (bx_detector::sub_detector d, bx_base_event::event_stage stage);
    bool check_required_event_stages (const bx_echidna_event *ev);

    void require_trigger_type (bx_trigger_event::trigger_type trg_type) { trg_type_requirements.push_back (trg_type); }
    bool check_trigger_type (const bx_echidna_event *ev);
    virtual void free_internal_buffers () {}
  private:
    module_role i_role;
    bool b_enabled;
    long int i4_priority;

    std::vector<bx_base_event::event_stage> trigger_requirements, laben_requirements, muon_requirements;
    std::vector<bx_trigger_event::trigger_type> trg_type_requirements;

    bx_base_module (const bx_base_module& t): bx_named (t.get_name ()) {} 
  public:
      // A list of module operation follow; a module operation can be used for
      // algorithms (see framework)
    struct begin_operation {
      begin_operation (module_role role): i_role(role) {}
      void operator() (bx_base_module *t) { t->framework_begin (i_role); }
      module_role i_role;
    };
    struct end_operation {
      end_operation (module_role role): i_role(role) {}
      void operator() (bx_base_module *t) { t->framework_end (i_role); }
      module_role i_role;
    };
    struct doit_operation {
      doit_operation (bx_echidna_event *ev, module_role role): ev_saved(ev), i_role(role) {}
      void operator() (bx_base_module *t) { t->framework_doit (ev_saved, i_role); }
      bx_echidna_event *ev_saved;
      module_role i_role;
    };
    struct compare_name_operation {
      compare_name_operation (const std::string& module_name): s_name(module_name) {}
      bool operator() (const bx_base_module *t) { return *t == s_name; }
      std::string s_name;
    };
    struct compare_role_operation {
      compare_role_operation (bx_base_module::module_role role): i_role(role) {}
      bool operator() (const bx_base_module *t) { return *t == i_role; }
      bx_base_module::module_role i_role;
    };
    struct sort_operation {
      sort_operation () {}
      bool operator() (const bx_base_module *a, const bx_base_module *b) { return *a < *b; }
    };
    struct delete_operation {
      delete_operation () {}
      void operator() (bx_base_module *t) { delete t; }
    };
    struct skip_event { 
      skip_event (const std::string& who): name(who) {}
      std::string name;
    };
};

inline std::ostream& operator<< (std::ostream &s, const bx_base_module* module) {
#if __GNUC__ >= 3
  s << "module: " << std::setw(30) << std::left << module->get_name ();
#else
  s << "module: " << std::setw(30) << module->get_name ();
#endif
  s << " role " << module->get_role ();
  s << " priority " << module->get_priority ();
  s << " enabled " << module->is_enabled ();
  return s;
}

#endif
/*
 * $Log: bx_base_module.hh,v $
 * Revision 1.30  2011/02/18 17:10:05  ddangelo
 * major code cleanup: removed fadc throughout the program
 *
 * Revision 1.29  2006-08-21 10:59:37  razeto
 * Added event skipping support
 *
 * Revision 1.28  2004/11/26 15:25:11  razeto
 * Thanks to Davide who so kindly underlined a syntactical error, I fixed the "Maintainer" word
 *
 * Revision 1.27  2004/11/26 14:01:20  razeto
 * Added Mantainer field
 *
 * Revision 1.26  2004/09/22 13:24:45  razeto
 * Changed bx_reco_event to bx_echidna_event
 *
 * Revision 1.25  2004/09/22 10:32:29  razeto
 * Moved sub_detector to bx_detector; added require_trigger_type
 *
 * Revision 1.24  2004/07/22 10:13:01  razeto
 * Jumbo patch for the configuration manager applied, read the documentation for usage
 *
 * Revision 1.23  2004/06/30 10:11:30  razeto
 * Added a level of indirection to the {begin|doit|end} call by the framework;
 * this is needed to mantain the bx_base_module infrastructure and to have
 * some check before calling the real module operation. These changes are
 * transparent for modules.
 *
 * Revision 1.22  2004/06/10 12:38:19  razeto
 * Added bx_named and parameter_broker interface
 *
 * Revision 1.21  2004/06/07 09:38:12  razeto
 * Added some operators for internal usage
 *
 * Revision 1.20  2004/05/31 14:57:37  razeto
 * Introduced new methods the modules can use to filter the events they receive:
 * using the require* methods a module can ask to receive in doit only the events
 * matching the require statements. Usefull for checking event status
 * before doing stuff.
 * Changed enabled method to is_enabled.
 *
 * Revision 1.19  2004/05/31 10:30:32  razeto
 * Upgraded to support named vdt upgrade
 *
 * Revision 1.18  2004/05/21 08:39:42  razeto
 * Added has data field
 *
 * Revision 1.17  2004/05/18 16:48:56  ddangelo
 * an #include changed to class declaration to avoid recursion
 *
 * Revision 1.16  2004/04/24 17:36:45  razeto
 * Added a check if the parameter already exist in setting (added an init method too)
 *
 * Revision 1.15  2004/04/18 10:30:14  razeto
 * Fixed to compile wit g++ 2.95 present on the cluster
 *
 * Revision 1.14  2004/04/12 16:06:21  razeto
 * Fixed a bug. Changed syntax of begin and and operator
 *
 * Revision 1.13  2004/04/06 12:42:17  razeto
 * Added module role to base_module, and the code to handle everywhere. Some other minor fixes
 *
 * Revision 1.12  2004/04/05 13:07:45  razeto
 * Added messanger for exception handling. Updated a comment
 *
 * Revision 1.11  2004/04/02 14:04:42  razeto
 * Added some comments. Fixed a typo
 *
 * Revision 1.10  2004/04/02 13:39:36  razeto
 * Added messenger. Changed from map to cmap
 *
 * Revision 1.9  2004/04/01 12:08:53  razeto
 * Moved the base_module_operation operators in the bx_base_module class
 * scope; moved the bx_base_module_vector typedef in the bx_base_module
 * class.
 * Added a event delete to the frame.
 *
 * Revision 1.8  2004/03/29 13:15:40  razeto
 * Updated bx_base_module to the vtd usage
 *
 * Revision 1.7  2004/03/26 16:33:59  razeto
 * Added check_parameter to see if a parameter is defined without generating an exception
 *
 * Revision 1.6  2004/03/23 13:01:08  razeto
 * Moved base_module_operatos to bx_base_module.hh.
 * Added compare_operation.
 * Added get_module by name to framework.
 *
 * Revision 1.5  2004/03/22 17:26:10  razeto
 * Added module parameter interface
 *
 * Revision 1.4  2004/03/22 14:30:19  razeto
 * Moved bx_base_module_vector typedef in bx_base_module.hh
 *
 * Revision 1.3  2004/03/21 18:52:27  razeto
 * Removed "using std::string". Some cosmetic changes
 *
 * Revision 1.2  2004/03/20 17:39:45  pallas
 * First release
 *
 * Revision 1.1  2004/03/20 14:16:47  razeto
 * Added sources
 *
 */
