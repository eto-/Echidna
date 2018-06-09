/* BOREXINO Reconstruction program
 *
 * Author: Livia Ludhova <livia.ludhova@mi.infn.it> and Alesandro Razeto <Alessandro.Razeto@ge.infn.it>
 * Maintainer: Livia Ludhova <livia.ludhova@mi.infn.it>
 *
 * 1) in: list of lg
 * 2) out: crate, fe, laben occupancy
 */

#ifndef _LABEN_EVENT_SHAPE_H
#define _LABEN_EVENT_SHAPE_H

#include "bx_named.hh"
#include <vector>


class laben_event_shape: public bx_named {
  public:
    laben_event_shape(const std::vector<int>&);

    const std::vector<float>& get_crate_occupancy (); 
    const std::vector<float>& get_hvb_occupancy ();   
    const std::vector<float>& get_feb_occupancy ();   
    const std::vector<float>& get_lbnb_occupancy (); 
    const std::vector<float>& get_channel_occupancy (); 

       //getters for occupancy not normalized to the number of hits in the event
    const std::vector<float>& get_nhits_in_crate ();
    const std::vector<float>& get_nhits_in_hvb ();
    const std::vector<float>& get_nhits_in_feb ();
    const std::vector<float>& get_nhits_in_lbnb ();
    const std::vector<float>& get_nhits_in_channel ();

  private:
    int nhits;
    int nlg;
    int n_good_lg;
    int count_bad_channels;
    bool cr_occupancy_done;
    bool hvb_occupancy_done;
    bool feb_occupancy_done;
    bool lbnb_occupancy_done;
    bool channel_occupancy_done;
    bool nhits_in_crate_done;
    bool nhits_in_hvb_done;
    bool nhits_in_lbnb_done;
    bool nhits_in_feb_done;
    bool nhits_in_channel_done;

    std::vector<int>   lg_map;
    std::vector<float> crate_occupancy;
    std::vector<float> hvb_occupancy;
    std::vector<float> feb_occupancy;
    std::vector<float> lbnb_occupancy;
    std::vector<float> channel_occupancy;
    std::vector<float> nhits_in_crate;
    std::vector<float> nhits_in_hvb;
    std::vector<float> nhits_in_feb;
    std::vector<float> nhits_in_lbnb;
    std::vector<float> nhits_in_channel;
    std::vector<bool>   lg_weight; 
    void init_lg_weight ();  

};
#endif






