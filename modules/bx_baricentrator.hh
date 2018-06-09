/* BOREXINO Reconstruction program
 *
 * Author: Alessandro Razeto <Alessandro.Razeto@ge.infn.it>
 * Maintainer: Alessandro Razeto <Alessandro.Razeto@ge.infn.it>
 *
 * $Id: bx_baricentrator.hh,v 1.11 2005/11/18 16:45:44 razeto Exp $
 *
 * Search charge baricenter: using the first hit of a cluster (defined
 * by a gate) calculate the charge baricenter; the uncertainty of the 
 * baricenter is measured calculating the mean quadratic difference 
 * between the time of flight of the photons (from baricenter to pmts) 
 * and the measured time of the hits.
 * 
 */
#ifndef _BX_BARICENTRATOR_H
#define _BX_BARICENTRATOR_H

#include "bx_base_module.hh"
#include "bx_rec_general.hh"
#include "TH1F.h"
#include "cmap.hh"
#include "flight_path.hh"
#include <string>

class bx_laben_clustered_hit;

class bx_baricentrator: public bx_base_module {
  public:
    bx_baricentrator ();
    virtual ~bx_baricentrator () {}

    virtual void begin ();
    virtual bx_echidna_event* doit (bx_echidna_event *ev);
    virtual void end ();

  private:
    double f8_first_hits_gate;
    double f8_Q_normalization;
    double f8_T_normalization, f8_T_filter_time;
    double f8_QT_normalization, f8_QT_filter_time;
    enum mean_algorithm {
      Q,
      T,
      QT,
    };
    mean_algorithm algo;
    bool b_use_charge;
    std::cmap<std::string, mean_algorithm> mean_algorithm_name_map;
    
    int get_hit_charge (const bx_laben_clustered_hit& hit);
    
    TH1F *baricenter_best_radius;
    
    flight_path path;
};

#endif
/*
 * $Log: bx_baricentrator.hh,v $
 * Revision 1.11  2005/11/18 16:45:44  razeto
 * Now using integer charge from decoded hit
 *
 * Revision 1.10  2005/06/18 15:28:46  razeto
 * Removed Minuit minimization, now this module only calculate the baricenter
 * position with some algorithms, see docs.
 *
 * Revision 1.9  2005/03/18 16:35:08  razeto
 * Updated to use flight_path
 *
 * Revision 1.8  2005/03/18 12:05:56  razeto
 * Removed some useless histos
 *
 * Revision 1.7  2005/03/01 15:17:32  razeto
 * Merged with cycle_2
 *
 * Revision 1.6.2.1  2004/12/14 16:38:44  razeto
 * Suppressed minuit warns but added echidna log messages on fail
 *
 * Revision 1.6  2004/12/03 13:36:09  razeto
 * Upgraded according to Sasha suggestions
 *
 * Revision 1.5  2004/11/26 15:25:11  razeto
 * Thanks to Davide who so kindly underlined a syntactical error, I fixed the "Maintainer" word
 *
 * Revision 1.4  2004/11/26 14:01:20  razeto
 * Added Mantainer field
 *
 * Revision 1.3  2004/09/22 13:27:06  razeto
 * Changed bx_reco_event to bx_echidna_event
 *
 * Revision 1.2  2004/09/22 10:44:50  razeto
 * Added an histogram
 *
 * Revision 1.1  2004/09/17 13:35:37  razeto
 * Added bx_baricentrator
 *
 */
