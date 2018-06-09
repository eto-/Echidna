/* BOREXINO Reconstruction program
 *
 * Author: Livia Ludhova <livia.ludhova@mi.infn.it> and Alesandro Razeto <Alessandro.Razeto@ge.infn.it>
 * Maintainer: Livia Ludhova <livia.ludhova@mi.infn.it>
 *
 * 1) in: list of lg
 * 2) out: crate, fe, laben occupancy
 */

#include "laben_event_shape.hh"
#include "bx_dbi.hh"
#include "db_run.hh"
#include "db_profile.hh"
#include "constants.hh"
#include "messenger.hh"
#include "bx_detector.hh"
#include "db_channel.hh"

#include <math.h>

laben_event_shape::laben_event_shape(const std::vector<int>& lg_list): bx_named("laben_event_shape"),
                                                    lg_map(constants::laben::channels, 0), 
						    crate_occupancy(constants::laben::ncrates, 0), 
						    hvb_occupancy(constants::laben::ncrates*constants::laben::frontend::board_per_rack/2, 0),
                                                    feb_occupancy(constants::laben::ncrates*constants::laben::frontend::board_per_rack, 0), 
						    lbnb_occupancy(constants::laben::ncrates*constants::laben::board_per_rack, 0),
						    channel_occupancy(constants::laben::channels, 0),
						    nhits_in_crate(constants::laben::ncrates, 0),
						    nhits_in_hvb(constants::laben::ncrates*constants::laben::frontend::board_per_rack/2, 0),
						    nhits_in_feb(constants::laben::ncrates*constants::laben::frontend::board_per_rack, 0), 
   					            nhits_in_lbnb(constants::laben::ncrates*constants::laben::board_per_rack, 0),
						    nhits_in_channel(constants::laben::channels, 0),
						    lg_weight(constants::laben::channels, true)
{

  for (unsigned int i=0; i < lg_list.size (); i++) lg_map[lg_list[i] - 1]++;
  count_bad_channels = (get_parameter_other ("detector_interface", "disabled_laben_channel_properties").get_vector ().size () > 0);
  nhits = lg_list.size ();
  nlg = constants::laben::channels;
  init_lg_weight ();
  cr_occupancy_done = false;
  hvb_occupancy_done = false;
  feb_occupancy_done = false;
  lbnb_occupancy_done = false;
  channel_occupancy_done = false;
  nhits_in_crate_done = false;
  nhits_in_hvb_done = false;
  nhits_in_feb_done = false;
  nhits_in_lbnb_done = false;
  nhits_in_channel_done = false;

}
 

void laben_event_shape::init_lg_weight () {
  lg_weight.resize (constants::laben::channels, 0);
  const std::vector<int> &bad_channels_list = detector_interface::get ()->get_disabled_channels ();
   
    //not consider bad channels if decoder fix_bad_channels parameter is on
  for (unsigned int i=0; i < bad_channels_list.size () && count_bad_channels; i++) {
    if (bad_channels_list[i] <= constants::laben::channels) 
      lg_weight[bad_channels_list[i] - 1] = false;
  }
  //not consider non-ordinary channels 
  for (int ilg = 1;  ilg <= nlg; ilg++) {
     if(!(bx_dbi::get ()->get_channel (ilg).is_ordinary ()) ) lg_weight[ilg - 1] = false; 
  }	 
}				       

  //crate occupancy = Nhits_in_crate_in_considered_lg / Nconsidered_lg_in_crate ( e.g. crate_weight) / Nhits_in_considered_lg * 168 (lg in crate)
const std::vector<float>& laben_event_shape::get_crate_occupancy (){
  if (cr_occupancy_done) return crate_occupancy;
  cr_occupancy_done = true;

  int n_crates = constants::laben::ncrates; //14
  int nhits_in_good_lg = 0;
    //crate weight = number of lg in crate which are considered
  std::vector<double> crate_weight(n_crates, 0.000001); //start at something !=0 to avoid division by 0 in normalisation
  for (int ilg = 1;  ilg <= nlg; ilg++) {
    int i_crate = constants::crate (ilg);
    if(lg_weight[ilg - 1]){
      crate_weight[i_crate - 1] ++;
      crate_occupancy[i_crate - 1] += lg_map[ilg - 1];
      nhits_in_good_lg += lg_map[ilg - 1];
    }
  } 

    //normalisation
  for (int cr = 0; cr < n_crates; cr ++) crate_occupancy[cr] = crate_occupancy[cr] / nhits_in_good_lg / crate_weight[cr] * 160;
  return crate_occupancy;
}

 //nhits_in_crate = Nhits_in_crate_in_considered_lg / Nconsidered_lg_in_crate ( e.g. crate_weight) * 168 (lg in crate)
const std::vector<float>& laben_event_shape::get_nhits_in_crate (){
  if (nhits_in_crate_done) return nhits_in_crate;
  nhits_in_crate_done = true;

  int n_crates = constants::laben::ncrates; //14
  int nhits_in_good_lg = 0;
    //crate weight = number of lg in crate which are considered
  std::vector<double> crate_weight(n_crates, 0.000001); //start at something !=0 to avoid division by 0 in normalisation
  for (int ilg = 1;  ilg <= nlg; ilg++) {
    int i_crate = constants::crate (ilg);
    if(lg_weight[ilg - 1]){
      crate_weight[i_crate - 1] ++;
      nhits_in_crate[i_crate - 1] += lg_map[ilg - 1];
      nhits_in_good_lg += lg_map[ilg - 1];
    }
  } 

    //normalisation
  for (int cr = 0; cr < n_crates; cr ++) nhits_in_crate[cr] = nhits_in_crate[cr] / crate_weight[cr] * 160;
  return nhits_in_crate;
}


 //hvb occupancy = Nhits_in_hvb_in_considered_lg / Nconsidered_lg_in_hvb ( e.g. feb_weight) / Nhits_in_considered_lg * 24 (lg in hvb)
const std::vector<float>& laben_event_shape::get_hvb_occupancy (){
  if (hvb_occupancy_done) return hvb_occupancy;
  hvb_occupancy_done = true;
  int n_crates = constants::laben::ncrates; //14
  int n_hvb = n_crates *  constants::laben::frontend::board_per_rack / 2; // 14 * 14/2 = 196/2 = 98
  int nhits_in_good_lg = 0;
    //hvb weight = number of lg in hvb which are considered
  std::vector<double> hvb_weight(n_hvb, 0.000001); //start at something !=0 to avoid division by 0 in normalisation
  for (int ilg = 1;  ilg <= nlg; ilg++) {
    int i_crate = constants::crate (ilg);
    int i_hv_board= (int) ((float) constants::laben::frontend::board_in_rack (ilg)/2. + 0.6);
    int i_hvb = (int) ((float) i_hv_board + (float) (i_crate - 1) * constants::laben::frontend::board_per_rack/2. + 0.6); 
    
    if(lg_weight[ilg - 1]){    
      hvb_weight[i_hvb - 1] ++;
      hvb_occupancy[i_hvb - 1] += lg_map[ilg - 1];
      nhits_in_good_lg += lg_map[ilg - 1];
    }
  } 
    //normalisation
  for (int hvb = 0; hvb < n_hvb; hvb ++) hvb_occupancy[hvb] = hvb_occupancy[hvb] / nhits_in_good_lg / hvb_weight[hvb] * 24;
    
  return hvb_occupancy;
}
 
 //nhits_in_hvb = Nhits_in_hvb_in_considered_lg / Nconsidered_lg_in_hvb ( e.g. feb_weight) * 24 (lg in hvb)
const std::vector<float>& laben_event_shape::get_nhits_in_hvb (){
  if (nhits_in_hvb_done) return nhits_in_hvb;
  nhits_in_hvb_done = true;
  int n_crates = constants::laben::ncrates; //14
  int n_hvb = n_crates *  constants::laben::frontend::board_per_rack / 2; // 14 * 14/2 = 196/2 = 98
  int nhits_in_good_lg = 0;
    //hvb weight = number of lg in hvb which are considered
  std::vector<double> hvb_weight(n_hvb, 0.000001); //start at something !=0 to avoid division by 0 in normalisation
  for (int ilg = 1;  ilg <= nlg; ilg++) {
    int i_crate = constants::crate (ilg);
    int i_hv_board= (int) ((float) constants::laben::frontend::board_in_rack (ilg)/2. + 0.6);
    int i_hvb = (int) ((float) i_hv_board + (float) (i_crate - 1) * constants::laben::frontend::board_per_rack/2. + 0.6); 
    
    if(lg_weight[ilg - 1]){    
      hvb_weight[i_hvb - 1] ++;
      nhits_in_hvb[i_hvb - 1] += lg_map[ilg - 1];
      nhits_in_good_lg += lg_map[ilg - 1];
    }
  } 
    //normalisation
  for (int hvb = 0; hvb < n_hvb; hvb ++) nhits_in_hvb[hvb] =  nhits_in_hvb[hvb] / hvb_weight[hvb] * 24;
    
  return  nhits_in_hvb;
}
 


 //feb occupancy = Nhits_in_feb_in_considered_lg / Nconsidered_lg_in_feb ( e.g. feb_weight) / Nhits_in_considered_lg * 12 (lg in feb)
const std::vector<float>& laben_event_shape::get_feb_occupancy (){
  if (feb_occupancy_done) return feb_occupancy;
  feb_occupancy_done = true;
  int n_crates = constants::laben::ncrates; //14
  int n_feb = n_crates *  constants::laben::frontend::board_per_rack; // 14 * 14 = 196
  int nhits_in_good_lg = 0;
    //feb weight = number of lg in feb which are considered
  std::vector<double> feb_weight(n_feb, 0.000001); //start at something !=0 to avoid division by 0 in normalisation
  for (int ilg = 1;  ilg <= nlg; ilg++) {
    int i_crate = constants::crate (ilg);
    int i_fe_board= constants::laben::frontend::board_in_rack (ilg);
    int i_feb = i_fe_board + (i_crate - 1) * constants::laben::frontend::board_per_rack;
    
    if(lg_weight[ilg - 1]){    
      feb_weight[i_feb - 1] ++;
      feb_occupancy[i_feb - 1] += lg_map[ilg - 1];
      nhits_in_good_lg += lg_map[ilg - 1];
    }
  } 
    //normalisation
  for (int feb = 0; feb < n_feb; feb ++) feb_occupancy[feb] = feb_occupancy[feb] / nhits_in_good_lg / feb_weight[feb] * 12;
    
  return feb_occupancy;
}
 
//nhits_in_feb = Nhits_in_feb_in_considered_lg / Nconsidered_lg_in_feb ( e.g. feb_weight) * 12 (lg in feb)
const std::vector<float>& laben_event_shape::get_nhits_in_feb (){
  if (nhits_in_feb_done) return nhits_in_feb;
  nhits_in_feb_done = true;
  int n_crates = constants::laben::ncrates; //14
  int n_feb = n_crates *  constants::laben::frontend::board_per_rack; // 14 * 14 = 196
  int nhits_in_good_lg = 0;
    //feb weight = number of lg in feb which are considered
  std::vector<double> feb_weight(n_feb, 0.000001); //start at something !=0 to avoid division by 0 in normalisation
  for (int ilg = 1;  ilg <= nlg; ilg++) {
    int i_crate = constants::crate (ilg);
    int i_fe_board= constants::laben::frontend::board_in_rack (ilg);
    int i_feb = i_fe_board + (i_crate - 1) * constants::laben::frontend::board_per_rack;
    
    if(lg_weight[ilg - 1]){    
      feb_weight[i_feb - 1] ++;
      nhits_in_feb[i_feb - 1] += lg_map[ilg - 1];
      nhits_in_good_lg += lg_map[ilg - 1];
    }
  } 
    //normalisation
  for (int feb = 0; feb < n_feb; feb ++) nhits_in_feb[feb] = nhits_in_feb[feb] / nhits_in_good_lg / feb_weight[feb] * 12;
    
  return nhits_in_feb;
}

 //laben occupancy = Nhits_in_laben_board_in_considered_lg / Nconsidered_lg_in_laben_board ( e.g. laben_weight) / Nhits_in_considered_lg * 8 (lg in laben_board)
const std::vector<float>& laben_event_shape::get_lbnb_occupancy (){
  if (lbnb_occupancy_done) return lbnb_occupancy;
  lbnb_occupancy_done = true;
  int n_crates = constants::laben::ncrates; //14
  int n_laben = n_crates * constants::laben::board_per_rack; // 14 * 20 = 280
  int nhits_in_good_lg = 0;
    //laben weight = number of lg in laben board  which are considered
  std::vector<double> laben_weight(n_laben, 0.000001); //start at something !=0 to avoid division by 0 in normalisation
  for (int ilg = 1;  ilg <= nlg; ilg++) {
    int i_lab = (ilg - 1) / constants::laben::channels_per_board + 1;
    if(lg_weight[ilg - 1]){    
      laben_weight[i_lab - 1] ++;
      lbnb_occupancy[i_lab - 1] += lg_map[ilg - 1];
      nhits_in_good_lg += lg_map[ilg - 1];
    }
  } 
    //normalization
  for (int laben = 0; laben < n_laben;  laben ++) lbnb_occupancy[laben] = lbnb_occupancy[laben] / nhits_in_good_lg / laben_weight[laben] * 8;
  
  return lbnb_occupancy;
}
 
 //nhits_in_laben = Nhits_in_laben_board_in_considered_lg / Nconsidered_lg_in_laben_board ( e.g. laben_weight) * 8 (lg in laben_board)
const std::vector<float>& laben_event_shape::get_nhits_in_lbnb (){
  if (nhits_in_lbnb_done) return nhits_in_lbnb;
  nhits_in_lbnb_done = true;
  int n_crates = constants::laben::ncrates; //14
  int n_laben = n_crates * constants::laben::board_per_rack; // 14 * 20 = 280
  int nhits_in_good_lg = 0;
    //laben weight = number of lg in laben board  which are considered
  std::vector<double> laben_weight(n_laben, 0.000001); //start at something !=0 to avoid division by 0 in normalisation
  for (int ilg = 1;  ilg <= nlg; ilg++) {
    int i_lab = (ilg - 1) / constants::laben::channels_per_board + 1;
    if(lg_weight[ilg - 1]){    
      laben_weight[i_lab - 1] ++;
      nhits_in_lbnb[i_lab - 1] += lg_map[ilg - 1];
      nhits_in_good_lg += lg_map[ilg - 1];
    }
  } 
    //normalization
  for (int laben = 0; laben < n_laben;  laben ++) nhits_in_lbnb[laben] = nhits_in_lbnb[laben] / laben_weight[laben] * 8;
  
  return nhits_in_lbnb;
}
 
 //channel occupancy = Nhits_in_considered_lg / Nhits_in_considered_lg 
const std::vector<float>& laben_event_shape::get_channel_occupancy (){
  if (channel_occupancy_done) return channel_occupancy;
  channel_occupancy_done = true;
  int nhits_in_good_lg = 0;
  for (int ilg = 1;  ilg <= nlg; ilg++) {
    if(lg_weight[ilg - 1]){    
      channel_occupancy[ilg - 1] += lg_map[ilg - 1];
      nhits_in_good_lg += lg_map[ilg - 1];
    }
  } 
    //normalisation
  for (int ilg = 0; ilg < nlg;  ilg ++) channel_occupancy[ilg] = channel_occupancy[ilg] / nhits_in_good_lg ;
  
  return channel_occupancy;
}

//hits_in_channel  = Nhits_in_considered_lg 
const std::vector<float>& laben_event_shape::get_nhits_in_channel (){
  if (nhits_in_channel_done) return nhits_in_channel;
  nhits_in_channel_done = true;
  int nhits_in_good_lg = 0;
  for (int ilg = 1;  ilg <= nlg; ilg++) {
    if(lg_weight[ilg - 1]){    
      nhits_in_channel[ilg - 1] += lg_map[ilg - 1];
      nhits_in_good_lg += lg_map[ilg - 1];
    }
  } 
  
  return nhits_in_channel;
}

