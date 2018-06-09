/* BOREXINO Reconstruction program
 *
 * Author: Alessandro Razeto <Alessandro.Razeto@ge.infn.it>
 * Maintainer: Alessandro Razeto <Alessandro.Razeto@ge.infn.it>
 *
 * $Id: bx_position_reco.hh,v 1.6 2009/10/23 14:00:04 koshio Exp $
 *
 */
#ifndef _BX_POSITION_RECO_H
#define _BX_POSITION_RECO_H

#include "bx_base_module.hh"
#include "bx_rec_general.hh"
#include "flight_path.hh"
#include "cmap.hh"
#include <string>

class bx_laben_rec_cluster;
class bx_laben_clustered_hit;
class bx_position_reco: public bx_base_module {
  public:
    bx_position_reco ();
    virtual ~bx_position_reco () {}

    virtual void begin ();
    virtual bx_echidna_event* doit (bx_echidna_event *ev);
    virtual void end ();

  private:
    enum position_algorithm {
      baricenter,
      milano,
      lngs,
      moscow,
      dubna,
      mach4,
      mctruth,
      fixed,
    } algo;
    std::cmap<std::string, position_algorithm> position_algorithm_name_map;
    float fixed_positions[3];
    
    int get_hit_charge (const bx_laben_clustered_hit& hit);
    void tof_hits (bx_laben_rec_cluster &rec_cluster);
    
    flight_path path;
};

#endif
/*
 * $Log: bx_position_reco.hh,v $
 * Revision 1.6  2009/10/23 14:00:04  koshio
 * Add the lngs postion reconstruction
 *
 * Revision 1.5  2009-07-16 15:53:00  razeto
 * Added mach4 position reco
 *
 * Revision 1.4  2008-12-15 11:46:55  razeto
 * Added fixed position (for sources)
 *
 * Revision 1.3  2005-12-03 15:18:57  razeto
 * Added mctruth algorith (which just copy mc data to bx_position)
 *
 * Revision 1.2  2005/10/11 14:13:55  razeto
 * Added a declaration
 *
 * Revision 1.1  2005/06/20 14:17:26  razeto
 * Added position reco
 *
 *
 */
