/* BOREXINO Reconstruction program
 *
 * Author: Alessandro Razeto <Alessandro.Razeto@ge.infn.it>
 * Maintainer: Alessandro Razeto <Alessandro.Razeto@ge.infn.it>
 *
 * $Id: bx_filter_trigger.cc,v 1.2 2009/10/26 11:19:37 ddangelo Exp $
 *
 */

#include "bx_filter_trigger.hh"
#include "bx_echidna_event.hh"

// ctor
bx_filter_trigger::bx_filter_trigger (): bx_base_module("bx_filter_trigger", bx_base_module::main_loop) {
  name_map_translation["neutrino"] = bx_trigger_event::neutrino;
  name_map_translation["muon"    ] = bx_trigger_event::muon;
  name_map_translation["neutron" ] = bx_trigger_event::neutron;
  name_map_translation["laser266"] = bx_trigger_event::laser266;
  name_map_translation["laser355"] = bx_trigger_event::laser355;
  name_map_translation["laser394"] = bx_trigger_event::laser394;
  name_map_translation["pulser"  ] = bx_trigger_event::pulser;
  name_map_translation["random"  ] = bx_trigger_event::random;
}

void bx_filter_trigger::begin () {
  const vdt::vdt_vector& v = get_parameter("skip_trigger_types").get_vector ();
  for (unsigned i = 0; i < v.size (); i++) {
    if(!name_map_translation.check (v[i].get_string ())) get_message (bx_message::error) << "unknown trigger type " << v[i] << dispatch;
    else skip_trigger_types.push_back (name_map_translation[v[i].get_string ()]);
  }
}


bx_echidna_event* bx_filter_trigger::doit (bx_echidna_event *ev) {
  for (unsigned int i = 0; i < skip_trigger_types.size (); i++)
    if (ev->get_trigger ().get_trigger_type () == skip_trigger_types[i]) return 0;

  return ev;
}

void bx_filter_trigger::end () {
}

/*
 * $Log: bx_filter_trigger.cc,v $
 * Revision 1.2  2009/10/26 11:19:37  ddangelo
 * trgtype 'muon_tct' (or 'muon_totc') renamed as 'neutron' after real use
 * trgtype 'muon_mtb' renamed as 'muon' as duplicity with the above one disappeared and 'mtb' is associated to btb_inputs flag
 *
 * Revision 1.1  2006-08-21 11:01:48  razeto
 * New filter module
 *
 */
