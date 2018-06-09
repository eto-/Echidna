/* BOREXINO Reconstruction program
 *
 * Author: Alessandro Razeto <Alessandro.Razeto@ge.infn.it>
 * 	   Marco Pallavicini <pallas@ge.infn.it>
 * Maintainer: Alessandro Razeto <Alessandro.Razeto@ge.infn.it>
 *
 * $Id: bx_base_module.cc,v 1.28 2011/02/18 17:10:05 ddangelo Exp $
 *
 * Implementation of bx_base_module
 *
 */
#include "bx_base_module.hh"
#include "bx_echidna_event.hh"
#include "db_profile.hh"
#include "bx_dbi.hh"


// ctor
bx_base_module::bx_base_module (const std::string& myname, bx_base_module::module_role role): bx_named(myname), i_role(role) {
  b_has_data = false;

  if (role == bx_base_module::all)
    get_message(bx_message::critic) << "can not istantiate a module with role = all" << dispatch;

  b_enabled = (check_parameter ("enable") && get_parameter ("enable").get_int ()) ? true : false;
}

void bx_base_module::set_parameter (const std::string& name, const vdt& value) {
  bx_named::set_parameter (name, value);
  if (name == "enable") {
    b_enabled = get_parameter ("enable").get_int () ? true : false;
  }
}

void bx_base_module::require_event_stage (bx_detector::sub_detector d, bx_base_event::event_stage stage) {
  switch (d) {
    case bx_detector::trigger:
      trigger_requirements.push_back (stage);
      break;
    case bx_detector::laben:
      laben_requirements.push_back (stage);
      break;
    case bx_detector::muon:
      muon_requirements.push_back (stage);
      break;
    case bx_detector::mctruth:
      get_message(bx_message::error) << "require_event_stage not supported for mctruth yet" << dispatch;
      break;
  }
}
bool bx_base_module::check_trigger_type (const bx_echidna_event *ev) {
  if (!trg_type_requirements.size ()) return true;
  
  for (unsigned i = 0; i < trg_type_requirements.size (); i++) {
    const db_profile& db = bx_dbi::get ()->get_profile ();
    const bx_trigger_event& trg = ev->get_trigger ();
    switch (trg_type_requirements[i]) {
      case bx_trigger_event::neutrino:
	if (trg.get_trgtype () == db.neutrino_trigger_tag ()) return true;
	break;
      case bx_trigger_event::muon:
	if (trg.get_trgtype () == db.muon_trigger_tag ()) return true;
	break;
      case bx_trigger_event::neutron:
	if (trg.get_trgtype () == db.neutron_trigger_tag ()) return true;
	break;
      case bx_trigger_event::laser266:
	if (trg.get_trgtype () == db.laser266_trigger_tag ()) return true;
	break;
      case bx_trigger_event::laser355:
	if (trg.get_trgtype () == db.laser355_trigger_tag ()) return true;
	break;
      case bx_trigger_event::laser394:
	if (trg.get_trgtype () == db.laser394_trigger_tag ()) return true;
	break;
      case bx_trigger_event::pulser:
	if (trg.get_trgtype () == db.pulser_trigger_tag ()) return true;
	break;
      case bx_trigger_event::random:
	if (trg.get_trgtype () == db.random_trigger_tag ()) return true;
	break;
    }
  }
  return false;
}

bool bx_base_module::check_required_event_stages (const bx_echidna_event *ev) {
  for (unsigned i = 0; i < trigger_requirements.size (); i++) 
    if (!ev->get_trigger ().check_stage (trigger_requirements[i])) return false;

  for (unsigned i = 0; i < laben_requirements.size (); i++) 
    if (!ev->get_laben ().check_stage (laben_requirements[i])) return false;

  for (unsigned i = 0; i < muon_requirements.size (); i++) 
    if (!ev->get_muon ().check_stage (muon_requirements[i])) return false;

  return true;
}

void bx_base_module::framework_begin (module_role role) { 
  if (is_enabled () && (role == get_role () || role == all)) 
    begin (); 
}

bx_echidna_event* bx_base_module::framework_doit (bx_echidna_event *ev, module_role role) {
  if (is_enabled () && (role == get_role () || role == all) && check_required_event_stages (ev) && check_trigger_type (ev)) 
    if (!doit (ev)) {
      detector_interface::get ()->skip_event (this, ev->get_trigger ().get_trgtype ());
      throw skip_event(get_name ());
    }

  return ev;
}

void bx_base_module::framework_end (module_role role) {
  if (is_enabled () && (role == get_role () || role == all)) 
    end ();
}

/*
 * $Log: bx_base_module.cc,v $
 * Revision 1.28  2011/02/18 17:10:05  ddangelo
 * major code cleanup: removed fadc throughout the program
 *
 * Revision 1.27  2009-10-26 11:19:37  ddangelo
 * trgtype 'muon_tct' (or 'muon_totc') renamed as 'neutron' after real use
 * trgtype 'muon_mtb' renamed as 'muon' as duplicity with the above one disappeared and 'mtb' is associated to btb_inputs flag
 *
 * Revision 1.26  2006-08-21 10:59:37  razeto
 * Added event skipping support
 *
 * Revision 1.25  2006/05/10 12:14:39  razeto
 * Removed some useless include
 *
 * Revision 1.24  2004/11/26 15:25:11  razeto
 * Thanks to Davide who so kindly underlined a syntactical error, I fixed the "Maintainer" word
 *
 * Revision 1.23  2004/11/26 14:01:20  razeto
 * Added Mantainer field
 *
 * Revision 1.22  2004/09/27 17:06:24  razeto
 * Fixed a bug
 *
 * Revision 1.21  2004/09/22 13:24:45  razeto
 * Changed bx_reco_event to bx_echidna_event
 *
 * Revision 1.20  2004/09/22 11:22:15  razeto
 * Added a warn
 *
 * Revision 1.19  2004/09/22 11:19:28  razeto
 * Added a null case to avoid a warning
 *
 * Revision 1.18  2004/09/22 10:32:41  razeto
 * Moved sub_detector to bx_detector; added require_trigger_type
 *
 * Revision 1.17  2004/07/22 10:13:01  razeto
 * Jumbo patch for the configuration manager applied, read the documentation for usage
 *
 * Revision 1.16  2004/07/07 13:28:23  razeto
 * Fixed a mistype
 *
 * Revision 1.15  2004/06/30 10:11:30  razeto
 * Added a level of indirection to the {begin|doit|end} call by the framework;
 * this is needed to mantain the bx_base_module infrastructure and to have
 * some check before calling the real module operation. These changes are
 * transparent for modules.
 *
 * Revision 1.14  2004/06/10 12:38:19  razeto
 * Added bx_named and parameter_broker interface
 *
 * Revision 1.13  2004/05/31 14:57:37  razeto
 * Introduced new methods the modules can use to filter the events they receive:
 * using the require* methods a module can ask to receive in doit only the events
 * matching the require statements. Usefull for checking event status
 * before doing stuff.
 * Changed enabled method to is_enabled.
 *
 * Revision 1.12  2004/05/31 10:30:32  razeto
 * Upgraded to support named vdt upgrade
 *
 * Revision 1.11  2004/05/21 08:39:42  razeto
 * Added has data field
 *
 * Revision 1.10  2004/04/24 17:36:45  razeto
 * Added a check if the parameter already exist in setting (added an init method too)
 *
 * Revision 1.9  2004/04/06 12:42:17  razeto
 * Added module role to base_module, and the code to handle everywhere. Some other minor fixes
 *
 * Revision 1.8  2004/04/03 09:21:16  razeto
 * changed the message syntax
 *
 * Revision 1.7  2004/04/02 13:39:36  razeto
 * Added messenger. Changed from map to cmap
 *
 * Revision 1.6  2004/03/29 13:15:40  razeto
 * Updated bx_base_module to the vtd usage
 *
 * Revision 1.5  2004/03/26 16:33:59  razeto
 * Added check_parameter to see if a parameter is defined without generating an exception
 *
 * Revision 1.4  2004/03/22 17:26:10  razeto
 * Added module parameter interface
 *
 * Revision 1.3  2004/03/21 19:00:07  razeto
 * Some cosmetic changes
 *
 * Revision 1.2  2004/03/21 18:52:27  razeto
 * Removed "using std::string". Some cosmetic changes
 *
 * Revision 1.1  2004/03/20 17:39:45  pallas
 * First release
 *
 * Revision 1.1  2004/03/20 14:16:47  razeto
 * Added sources
 *
 */
