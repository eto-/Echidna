/* BOREXINO Reconstruction program
 *
 * Author: Alessandro Razeto <Alessandro.Razeto@ge.infn.it>
 * 	   Marco Pallavicini <pallas@ge.infn.it>
 * Maintainer: Davide D'Angelo <davide.dangelo@lngs.infn.it>
 *
 * $Id: bx_base_event.hh,v 1.17 2009/10/23 14:00:03 koshio Exp $
 *
 * This is the pure base class for the event hierarchy
 * 
 */
#ifndef _BX_BASE_EVENT_HH
#define _BX_BASE_EVENT_HH

#include <vector>

#include "bx_rec_general.hh"

class bx_base_event {
  public:
    enum event_stage {
      raw = 0,
      decoded,
      clustered,
      baricentered,
      reconstructed_mi,
      reconstructed_lngs,
      reconstructed_msk,
      reconstructed_dbn,
      reconstructed_mach4,
      reconstructed,
      tracked,
      split,
      pid_ab,
      max // to be pushed higher when needed
    };
    const bool check_stage (event_stage v) const { return stages[v]; }
                                                                                                            
    bx_base_event () { stages.insert (stages.end (), max + 1, false); }
    virtual ~bx_base_event () {}

  protected:
    std::vector<bool> stages; 
};

#endif
/*
 * $Log: bx_base_event.hh,v $
 * Revision 1.17  2009/10/23 14:00:03  koshio
 * Add the lngs postion reconstruction
 *
 * Revision 1.16  2009-07-16 10:17:42  ddangelo
 * infrastructure for m4 position reco (patch by steve&ale)
 *
 * Revision 1.15  2008-02-22 11:48:22  ddangelo
 * event staging slitghly reorganized
 *
 * Revision 1.14  2008-02-21 15:58:31  ddangelo
 * added muon track as a separate class.
 * Currently used twice: laben and muon (later also globally fit).
 * both inner and roor event
 * added tracked stage in base event enum
 * bx_muon_rec_event renamed as bx_muon_tracked_event
 *
 * Revision 1.13  2005-07-13 14:15:26  ddangelo
 * removed rough clustered stage
 *
 * Revision 1.12  2005/07/11 17:11:45  ddangelo
 * removed global event
 * added pid event
 *
 * Revision 1.11  2005/06/20 16:45:23  ddangelo
 * added bx_position_reco_dbn
 *
 * Revision 1.10  2004/12/14 20:02:03  ddangelo
 * 2 stages added.
 *
 * Revision 1.9  2004/11/26 15:25:10  razeto
 * Thanks to Davide who so kindly underlined a syntactical error, I fixed the "Maintainer" word
 *
 * Revision 1.8  2004/11/26 14:06:22  razeto
 * Added Mantainer field
 *
 * Revision 1.7  2004/09/22 11:07:59  ddangelo
 * added baricentered event stage.
 * Indentation fixed to match Alessandro's "manias"
 *
 * Revision 1.6  2004/07/09 13:38:15  ddangelo
 * enum with event stages enlarged
 *
 * Revision 1.5  2004/05/25 16:37:05  ddangelo
 * added support for event stage flagging.
 * Inheritance from bx_base_event introduced
 *
 * Revision 1.4  2004/04/16 15:02:22  pallas
 * Start developing of root event structure
 * Not defined yet, but preprared infrastructure
 *
 * Revision 1.3  2004/03/20 14:13:11  razeto
 * Added main ifndef directive
 *
 * Revision 1.2  2004/03/20 12:27:18  razeto
 * Added bx_rec_general.hh include, fixed some tabs and removed a * line from logs
 *
 * Revision 1.1.1.1  2004/03/19 18:23:49  razeto
 * Imported sources
 *
 */
