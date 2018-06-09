/* BOREXINO Reconstruction program
 *
 * Author: Alberto Di Cienzo <Alberto.dicienzo@lngs.infn.it>
 * Maintainer: Alberto Di Cienzo <Alberto.dicienzo@lngs.infn.it>
 *
 * $Id: bx_snews.hh,v 1.1 2008/10/23 09:13:46 dicienzo Exp $
 *
 * 
 */

#ifndef _BX_SNEWS_HH
#define _BX_SNEWS_HH

#include "bx_base_module.hh"
#include "bx_rec_general.hh"
#include <list>

class bx_snews: public bx_base_module {
  public:
  // this section if fixed; do not edit
    bx_snews ();
    virtual ~bx_snews () {}

    virtual void begin ();
    virtual bx_echidna_event* doit (bx_echidna_event *ev);
    virtual void end ();
    
    struct EventData{
      double 		Time;
      int		NHits;
      unsigned long 	ID;
    };
    
  private:
    std::list <EventData> Train;
    bool registered;
    double last_linger;
    int linger;
    int allarmreturn;
    bool connect_brain;
};
#endif
/*
 * $Log: bx_snews.hh,v $
 * Revision 1.1  2008/10/23 09:13:46  dicienzo
 * Added snews monitor
 *
 */
