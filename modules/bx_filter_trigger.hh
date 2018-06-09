/* BOREXINO Reconstruction program
 *
 * Author: Alessandro Razeto <Alessandro.Razeto@ge.infn.it>
 * Maintainer: Alessandro Razeto <Alessandro.Razeto@ge.infn.it>
 *
 * $Id: bx_filter_trigger.hh,v 1.1 2006/08/21 11:01:48 razeto Exp $
 *
 * Simple module to skip all event with trigger type different from the
 * specified ones
 * 
 */

#ifndef _BX_FILTER_TRIGGER_HH
#define _BX_FILTER_TRIGGER_HH

#include "bx_base_module.hh"
#include "bx_trigger_event.hh"

class bx_filter_trigger: public bx_base_module {
  public:
    bx_filter_trigger ();
    virtual ~bx_filter_trigger () {}

    virtual void begin ();
    virtual bx_echidna_event* doit (bx_echidna_event *ev);
    virtual void end ();
  private:
    std::vector<bx_trigger_event::trigger_type> skip_trigger_types;
    std::cmap<std::string, bx_trigger_event::trigger_type> name_map_translation;
};

#endif
/*
 * $Log: bx_filter_trigger.hh,v $
 * Revision 1.1  2006/08/21 11:01:48  razeto
 * New filter module
 *
 */
