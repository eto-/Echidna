/* BOREXINO Reconstruction program
 *
 * Author: Davide D'Angelo <davide.dangelo@mi.infn.it> (starting from bx_calib_muon_ectronics)
 * Maintainer: Davide D'Angelo <davide.dangelo@mi.infn.it> 
 *
 * $Id: bx_calib_muon_electronics.cc,v 1.9 2009/10/26 11:19:37 ddangelo Exp $
 * 
 * Implemenentation of bx_calib_muon_electronics
 *
 */

#include <vector>
#include "bx_calib_muon_electronics.hh"
#include "bx_echidna_event.hh"
#include "bx_dbi.hh"
#include "db_channel.hh"
#include "db_run.hh"
#include "barn_interface.hh"

bx_calib_muon_electronics::bx_calib_muon_electronics() : bx_base_module("bx_calib_muon_electronics", bx_base_module::main_loop) {
  require_event_stage (bx_detector::muon, bx_base_event::decoded);
/*  require_trigger_type (bx_trigger_event::neutrino);
  require_trigger_type (bx_trigger_event::pulser);
  require_trigger_type (bx_trigger_event::laser394);
  require_trigger_type (bx_trigger_event::muon);*/

  multiplicity_translation_map["dead_in_pulser"          ] = dead_in_pulser          ;
  multiplicity_translation_map["dead_in_laser"           ] = dead_in_laser           ;
  multiplicity_translation_map["dead_in_neutrino"        ] = dead_in_neutrino        ;
  multiplicity_translation_map["dead_in_muon"            ] = dead_in_muon            ;
  multiplicity_translation_map["low_eff_in_pulser"       ] = low_eff_in_pulser       ;	      
  multiplicity_translation_map["low_eff_in_laser"        ] = low_eff_in_laser        ;	     
  multiplicity_translation_map["low_eff_in_neutrino"     ] = low_eff_in_neutrino     ;	     
  multiplicity_translation_map["low_eff_in_muon"         ] = low_eff_in_muon         ;	     
  multiplicity_translation_map["hot_in_pulser"           ] = hot_in_pulser           ;	     
  multiplicity_translation_map["hot_in_laser"            ] = hot_in_laser            ;	     
  multiplicity_translation_map["hot_in_neutrino"         ] = hot_in_neutrino         ;	    
  multiplicity_translation_map["hot_in_muon"             ] = hot_in_muon             ;	    
  multiplicity_translation_map["retriggering_in_pulser"  ] = retriggering_in_pulser  ;  
  multiplicity_translation_map["retriggering_in_laser"   ] = retriggering_in_laser   ;   
  multiplicity_translation_map["retriggering_in_neutrino"] = retriggering_in_neutrino;
  multiplicity_translation_map["retriggering_in_muon"    ] = retriggering_in_muon    ;

  trg_names[pulser]   = "pulser"  ;
  trg_names[laser]    = "laser"   ;
  trg_names[neutrino] = "neutrino"; 
  trg_names[muon]     = "muon"    ; 
}
      
//BEGIN
void bx_calib_muon_electronics::begin () {

  // Configuration parameters
  f_dead         = get_parameter ("dead_channel_thresh").get_float ();
  f_low_eff[0]   = get_parameter ("low_eff_thresh"     ).get_float ();
  f_hot    [0]   = get_parameter ("hot_thresh"         ).get_float ();
  f_low_eff[1]   = get_parameter ("low_eff_thresh"     ).get_float ();
  f_hot    [1]   = get_parameter ("hot_thresh"         ).get_float ();
  f_low_eff[2]   = get_parameter ("low_eff_thresh"     ).get_float ();
  f_hot    [2]   = get_parameter ("hot_thresh"         ).get_float ();
  f_low_eff[3]   = get_parameter ("low_eff_thresh"     ).get_float ();
  f_hot    [3]   = get_parameter ("hot_thresh"         ).get_float ();
  f_retriggering = get_parameter ("retriggering_thresh").get_float ();

  // Define histograms
  h_nhits_vs_mch[pulser]   = new TH1F ("h_nhits_vs_mch_pulser"  ,"decoded_hits per channel, pulser_trg"  , constants::muon::channels, 0, constants::muon::channels);  
  h_nhits_vs_mch[laser]    = new TH1F ("h_nhits_vs_mch_laser"   ,"decoded_hits per channel, laser_trg"   , constants::muon::channels, 0, constants::muon::channels);  
  h_nhits_vs_mch[neutrino] = new TH1F ("h_nhits_vs_mch_neutrino","decoded_hits per channel, neutrino_trg", constants::muon::channels, 0, constants::muon::channels);
  h_nhits_vs_mch[muon]     = new TH1F ("h_nhits_vs_mch_muon"    ,"decoded_hits per channel, muon_trg"    , constants::muon::channels, 0, constants::muon::channels);
  h_dt_vs_mch   [pulser]   = new TH2F ("h_dt_vs_mch_pulser"     ,"time 2nd_hit - 1st_hit, pulser_trg"    , constants::muon::channels, 0, constants::muon::channels, 2000, 0., 2000.);
  h_dt_vs_mch   [laser]    = new TH2F ("h_dt_vs_mch_laser"      ,"time 2nd_hit - 1st_hit, laser_trg"     , constants::muon::channels, 0, constants::muon::channels, 2000, 0., 2000.);
  h_dt_vs_mch   [neutrino] = new TH2F ("h_dt_vs_mch_neutrino"   ,"time 2nd_hit - 1st_hit, neutrino_trg"  , constants::muon::channels, 0, constants::muon::channels, 2000, 0., 2000.);
  h_dt_vs_mch   [muon]     = new TH2F ("h_dt_vs_mch_muon"       ,"time 2nd_hit - 1st_hit, muon_trg"      , constants::muon::channels, 0, constants::muon::channels, 2000, 0., 2000.);

  for(int tt = 0; tt < 4; tt++) {
    barn_interface::get ()->store (barn_interface::file, h_nhits_vs_mch[tt], this);
    barn_interface::get ()->store (barn_interface::file, h_dt_vs_mch   [tt], this);

    nevents_per_tt[tt] = 0;

      // Internal vectors and arrays
    nhits_per_tt_mch[tt] = new int[constants::muon::channels];
    std::fill_n (nhits_per_tt_mch[tt], constants::muon::channels, 0); 
  }

  get_message(bx_message::debug) << "begin" << dispatch;
}


//DOIT
bx_echidna_event* bx_calib_muon_electronics::doit (bx_echidna_event *ev) {

  int tt;
  if      (ev->get_trigger ().is_pulser   ()) tt = pulser;
  else if (ev->get_trigger ().is_laser394 ()) tt = laser;
  else if (ev->get_trigger ().is_muon     () || ev->get_trigger ().has_btb_flag(bx_trigger_raw_event::mtb_flag)) tt = muon;
  else if (ev->get_trigger ().is_neutrino () || ev->get_trigger ().is_random ()) tt = neutrino;
  else return ev;

//  get_message(bx_message::info) << "doit " << ev->get_event_number() << " tt " << tt << dispatch;

  nevents_per_tt[tt]++;
 
    //vectors for the calculation of retrigger dt
  std::vector<int>   nhits_per_mch_singev(constants::muon::channels, 0);
  std::vector<float> time_per_mch_singev (constants::muon::channels, 0);

    //loop on decoded hits 
  for (int i = 0; i < ev->get_muon ().get_decoded_nhits (); i++) {
    const bx_muon_decoded_hit& dhit = ev->get_muon ().get_decoded_hit (i);
    int   mch        = dhit.get_raw_hit().get_muon_channel();
    float hit_time   = dhit.get_time   ();

    nhits_per_tt_mch[tt][mch]++;
    nhits_per_mch_singev[mch]++;
    h_nhits_vs_mch[tt]->Fill(mch);
    
      //retrigger
    if (nhits_per_mch_singev[mch] > 1) {
      float dt = hit_time - time_per_mch_singev[mch];
      h_dt_vs_mch[tt]->Fill (mch, dt);
    }
    time_per_mch_singev[mch] = hit_time;
        
  } // end of loop on decoded hits
  
  return ev;  
}


//END
void bx_calib_muon_electronics::end () {

    //vectors for later setting of visitors 
  std::vector<std::vector<multiplicity> > multiplicity_pl_vec    (constants::muon::channels); 
  std::vector<std::vector<multiplicity> > multiplicity_n_vec     (constants::muon::channels);
  std::vector<std::vector<multiplicity> > prev_multiplicity_n_vec(constants::muon::channels);
  std::vector<std::vector<multiplicity> > multiplicity_m_vec     (constants::muon::channels);
  std::vector<std::vector<multiplicity> > prev_multiplicity_m_vec(constants::muon::channels);

  //inizialize vectors
  for(int ich = 0; ich < constants::muon::channels; ich++) {
    multiplicity_pl_vec    .push_back (std::vector<multiplicity>());
    multiplicity_n_vec     .push_back (std::vector<multiplicity>());
    prev_multiplicity_n_vec.push_back (std::vector<multiplicity>());
    multiplicity_m_vec     .push_back (std::vector<multiplicity>());
    prev_multiplicity_m_vec.push_back (std::vector<multiplicity>());
  }

    //find nch_x from the previous run 
  int prev_nch_dead        [4] = {};
  int prev_nch_low_eff     [4] = {};
  int prev_nch_hot         [4] = {};
  int prev_nch_retriggering[4] = {};
  
  int module_says_DB_write = 1;
             
    //vectors of mch, input 0 = does not have this characteristics, input 1 = has this characteristics, 
  std::vector<int> prev_mch_dead_p        (constants::muon::channels, 0);
  std::vector<int> prev_mch_dead_l        (constants::muon::channels, 0);
  std::vector<int> prev_mch_dead_n        (constants::muon::channels, 0);
  std::vector<int> prev_mch_dead_m        (constants::muon::channels, 0);
  std::vector<int> prev_mch_low_eff_p     (constants::muon::channels, 0);
  std::vector<int> prev_mch_low_eff_l     (constants::muon::channels, 0);
  std::vector<int> prev_mch_low_eff_n     (constants::muon::channels, 0);
  std::vector<int> prev_mch_low_eff_m     (constants::muon::channels, 0);
  std::vector<int> prev_mch_hot_p         (constants::muon::channels, 0);
  std::vector<int> prev_mch_hot_l         (constants::muon::channels, 0);
  std::vector<int> prev_mch_hot_n         (constants::muon::channels, 0);
  std::vector<int> prev_mch_hot_m         (constants::muon::channels, 0);
  std::vector<int> prev_mch_retriggering_p(constants::muon::channels, 0);
  std::vector<int> prev_mch_retriggering_l(constants::muon::channels, 0);
  std::vector<int> prev_mch_retriggering_n(constants::muon::channels, 0);
  std::vector<int> prev_mch_retriggering_m(constants::muon::channels, 0);
  
  db_run& run_info_prev = bx_dbi::get ()->get_run ();

  // loop on channels to extract prev run info
  for(int ich = 0; ich < constants::muon::channels; ich++) {
     
    const std::vector<std::string>& muon_multiplicity_status_v = run_info_prev.get_muon_multiplicity  (ich+constants::muon::channel_offset+1);
    for(unsigned i = 0; i < muon_multiplicity_status_v.size (); i ++ ){
      switch (multiplicity_translation_map[muon_multiplicity_status_v[i]]) {
        case dead_in_pulser:            prev_nch_dead        [0]++; prev_mch_dead_p        [ich] = 1; break;
        case dead_in_laser:             prev_nch_dead        [1]++; prev_mch_dead_l        [ich] = 1; break;
        case dead_in_neutrino:          prev_nch_dead        [2]++; prev_mch_dead_n        [ich] = 1; prev_multiplicity_n_vec[ich].push_back(dead_in_neutrino)        ; break;
        case dead_in_muon:              prev_nch_dead        [3]++; prev_mch_dead_m        [ich] = 1; prev_multiplicity_m_vec[ich].push_back(dead_in_muon)            ; break;
        case low_eff_in_pulser:         prev_nch_low_eff     [0]++; prev_mch_low_eff_p     [ich] = 1; break;
        case low_eff_in_laser:          prev_nch_low_eff     [1]++; prev_mch_low_eff_l     [ich] = 1; break;
        case low_eff_in_neutrino:       prev_nch_low_eff     [2]++; prev_mch_low_eff_n     [ich] = 1; prev_multiplicity_n_vec[ich].push_back(low_eff_in_neutrino)     ; break;
        case low_eff_in_muon:           prev_nch_low_eff     [3]++; prev_mch_low_eff_m     [ich] = 1; prev_multiplicity_m_vec[ich].push_back(low_eff_in_muon)         ; break;
        case hot_in_pulser:             prev_nch_hot         [0]++; prev_mch_hot_p         [ich] = 1; break;
        case hot_in_laser:              prev_nch_hot         [1]++; prev_mch_hot_l         [ich] = 1; break;
        case hot_in_neutrino:           prev_nch_hot         [2]++; prev_mch_hot_n         [ich] = 1; prev_multiplicity_n_vec[ich].push_back(hot_in_neutrino)         ; break;
        case hot_in_muon:               prev_nch_hot         [3]++; prev_mch_hot_m         [ich] = 1; prev_multiplicity_m_vec[ich].push_back(hot_in_muon)             ; break;
        case retriggering_in_pulser:    prev_nch_retriggering[0]++; prev_mch_retriggering_p[ich] = 1; break;
        case retriggering_in_laser:     prev_nch_retriggering[1]++; prev_mch_retriggering_l[ich] = 1; break; 
        case retriggering_in_neutrino:  prev_nch_retriggering[2]++; prev_mch_retriggering_n[ich] = 1; prev_multiplicity_n_vec[ich].push_back(retriggering_in_neutrino); break;
        case retriggering_in_muon:      prev_nch_retriggering[3]++; prev_mch_retriggering_m[ich] = 1; prev_multiplicity_m_vec[ich].push_back(retriggering_in_muon)    ; break;
      } 
    } // end of loop on flags vector elements
    
  } // end of loop on channels

  // dump prev run situation
  for (int tt = 0; tt < 4; tt++) get_message(bx_message::log) << "prev_nch_dead_in_"         << trg_names[tt] << " " << prev_nch_dead        [tt] << dispatch;
  for (int tt = 0; tt < 4; tt++) get_message(bx_message::log) << "prev_nch_low_eff_in_"      << trg_names[tt] << " " << prev_nch_low_eff     [tt] << dispatch;
  for (int tt = 0; tt < 4; tt++) get_message(bx_message::log) << "prev_nch_hot_in_"          << trg_names[tt] << " " << prev_nch_hot         [tt] << dispatch;
//  for (int tt = 0; tt < 4; tt++) get_message(bx_message::log) << "prev_nch_retriggering_in_" << trg_names[tt] << " " << prev_nch_retriggering[tt] << dispatch;
								      
    // current run
  int nch_dead        [4] = {};
  int nch_low_eff     [4] = {};
  int nch_hot         [4] = {};
  //int nch_retriggering[4] = {};
  
    //last two bins of the following vectors are used to store
    //fore-last: number of mch which did loose the characteristics
    //     last: number of mch which did "gain" the characteristics
  std::vector<int> mch_dead_p        (constants::muon::channels + 2, 0);
  std::vector<int> mch_dead_l        (constants::muon::channels + 2, 0);
  std::vector<int> mch_dead_n        (constants::muon::channels + 2, 0);
  std::vector<int> mch_dead_m        (constants::muon::channels + 2, 0);
  std::vector<int> mch_low_eff_p     (constants::muon::channels + 2, 0);
  std::vector<int> mch_low_eff_l     (constants::muon::channels + 2, 0);
  std::vector<int> mch_low_eff_n     (constants::muon::channels + 2, 0);
  std::vector<int> mch_low_eff_m     (constants::muon::channels + 2, 0);
  std::vector<int> mch_hot_p         (constants::muon::channels + 2, 0);
  std::vector<int> mch_hot_l         (constants::muon::channels + 2, 0);
  std::vector<int> mch_hot_n         (constants::muon::channels + 2, 0);
  std::vector<int> mch_hot_m         (constants::muon::channels + 2, 0);
  std::vector<int> mch_retriggering_p(constants::muon::channels + 2, 0);
  std::vector<int> mch_retriggering_l(constants::muon::channels + 2, 0);
  std::vector<int> mch_retriggering_n(constants::muon::channels + 2, 0);
  std::vector<int> mch_retriggering_m(constants::muon::channels + 2, 0);
  
  double mean_nhits[4] = {};
  double rms_nhits [4] = {};

    // calculate mean and rms for h_nhits_pulser, laser, neutrino, excluding dead in pulser mch 
    // (This amchorithm is due to Knuth,[1] who cites Welford.[2] from Wikipedia
  for(int tt = 0; tt < 4; tt++) {
    int nch_to_check = 0; 
    float mean = 0;     
    float S = 0;  
    if (!nevents_per_tt[tt]) continue;
    for (int ich = 0; ich < constants::muon::channels; ich++) {
      if (!bx_dbi::get ()->get_channel (ich+constants::muon::channel_offset+1).is_ordinary()) continue; //ref channels NOT used to calculate mean
      float eff_pulser = nhits_per_tt_mch[0][ich]/ (float)nevents_per_tt[0] ;
      // not dead_in_pulser
      if ( eff_pulser > f_dead ) {
	nch_to_check++;
	float delta = nhits_per_tt_mch[tt][ich] - mean;
	mean += delta / nch_to_check;
	S += delta * (nhits_per_tt_mch[tt][ich] - mean);
      }
    }
    mean_nhits[tt] = mean;

    if (!nch_to_check) {
      get_message(bx_message::error) << "All mch appear dead in pulser" << dispatch;  
      module_says_DB_write = 0;
      continue;
    } 

    //calculate rms
    rms_nhits[tt] = std::sqrt(S / nch_to_check);
    get_message(bx_message::log) << "Mean nhits per ch, tt " << trg_names[tt] << ":  " << mean_nhits[tt] << " +- " << rms_nhits[tt] << " in " << nevents_per_tt[tt] <<  " triggers " << dispatch;
    //check if the pulser mean is compatible with number of pulser triggers
    if(tt == 0 && (mean_nhits[0] > 1.2 * nevents_per_tt[0] || mean_nhits[0] < 0.8 * nevents_per_tt[0]) ) 
      get_message(bx_message::info) << "Mean nhits per mch "<< mean_nhits[0] << " +- " << rms_nhits[0] << " in pulser differs too much from the number of pulser triggers: " << nevents_per_tt[0] << dispatch;
    //check if the laser mean is compatible with number of laser triggers, considering led efficiency
    if(tt == 1 && mean_nhits[1] < 0.02 * nevents_per_tt[1]) 
      get_message(bx_message::warn) << "Mean nhits per mch "<< mean_nhits[1] << " +- " << rms_nhits[1] << " in laser is too low for the number of laser triggers: " << nevents_per_tt[1] << dispatch;
    if(tt == 1 && mean_nhits[1] > 0.1 * nevents_per_tt[1])  
      get_message(bx_message::warn) << "Mean nhits per mch "<< mean_nhits[1] << " +- " << rms_nhits[1] << " in laser is too high for the number of laser triggers: " << nevents_per_tt[1] << dispatch;
       
    //if rms and sqrt(mean) too different, calculate truncated mean and rms
    if(rms_nhits[tt] > 5 * sqrt(mean_nhits[tt])){
      get_message(bx_message::log) << "nhits distribution in " <<  trg_names[tt] << " fluctuates, rms/sqrt(mean) = " << int(rms_nhits[tt] / sqrt(mean_nhits[tt])) 
	       << " trying truncated mean/rms..." << dispatch;
      nch_to_check = 0; 
      mean = 0;     
      S = 0; 
      for(int ich = 0; ich < constants::muon::channels ; ich++) {
	if (!bx_dbi::get ()->get_channel (ich+constants::muon::channel_offset+1).is_ordinary()) continue;  //ref channels NOT used to calculate mean
	float eff_pulser = nhits_per_tt_mch[0][ich]/ (float)nevents_per_tt[0] ;
	//not dead_in_pulser and nhits is 1 rms around mean value, use the value to calculate truncated mean and rms
	if (eff_pulser > f_dead && std::fabs(nhits_per_tt_mch[tt][ich] - mean_nhits[tt]) < rms_nhits[tt]) {
          nch_to_check++;
	  float delta = nhits_per_tt_mch[tt][ich] - mean;
	  mean += delta / nch_to_check;
	  S += delta * (nhits_per_tt_mch[tt][ich] - mean);
	}
      }
      mean_nhits[tt] = mean;
      rms_nhits [tt] = std::sqrt(S / nch_to_check);
      get_message(bx_message::log) << "Mean nhits per ch, tt " << trg_names[tt] << ":  " << mean_nhits[tt] << " +- " << rms_nhits[tt] << " in " << nevents_per_tt[tt] <<  " triggers (truncated)" << dispatch;
    } // end of truncated mean/rms calculation
  } // end of loop on tt

  // Status of the single electronic channel (dead, hot, low_eff)
  for(int tt = 0; tt < 4; tt++) {
    for(int ich = 0; ich < constants::muon::channels; ich++) {

      if(!bx_dbi::get ()->get_channel (ich+constants::muon::channel_offset+1).is_ordinary()) continue;
      
      if(!nevents_per_tt[tt] || !rms_nhits[tt]) continue;

      // is dead  ?
      if (nhits_per_tt_mch[tt][ich] <= f_dead * mean_nhits[tt]) {
	nch_dead[tt]++;
	if(tt == 0) mch_dead_p[ich] = 1;
	if(tt == 1) mch_dead_l[ich] = 1;
	if(tt == 2) mch_dead_n[ich] = 1;
	if(tt == 3) mch_dead_m[ich] = 1;
	get_message(bx_message::log) << "Mch " << ich << ": Dead_in_" << trg_names[tt] << "  (" << nhits_per_tt_mch[tt][ich] << " in " << nevents_per_tt[tt] << "  triggers )" << dispatch;
	if(tt == 0) multiplicity_pl_vec[ich].push_back(dead_in_pulser  );
	if(tt == 1) multiplicity_pl_vec[ich].push_back(dead_in_laser   );
	if(tt == 2) multiplicity_n_vec [ich].push_back(dead_in_neutrino);
	if(tt == 3) multiplicity_m_vec [ich].push_back(dead_in_muon    );
	continue; // if dead no other check makes sense
      }
      // is hot?
      if (nhits_per_tt_mch[tt][ich] > (mean_nhits[tt] + f_hot[tt] * rms_nhits[tt]) ){
	nch_hot[tt]++;
	if(tt == 0) mch_hot_p[ich] = 1;
	if(tt == 1) mch_hot_l[ich] = 1;
	if(tt == 2) mch_hot_n[ich] = 1;
	if(tt == 3) mch_hot_m[ich] = 1;
	get_message(bx_message::log) << "Mch " << ich << ": Hot_in_" << trg_names[tt] << " " << nhits_per_tt_mch[tt][ich] << " hits in " << nevents_per_tt[tt] << " triggers" << dispatch;
	if(tt == 0) multiplicity_pl_vec[ich].push_back(hot_in_pulser  );
	if(tt == 1) multiplicity_pl_vec[ich].push_back(hot_in_laser   );
	if(tt == 2) multiplicity_n_vec [ich].push_back(hot_in_neutrino);
	if(tt == 3) multiplicity_m_vec [ich].push_back(hot_in_muon    );
      }
      // has low_eff ?
      if(nhits_per_tt_mch[tt][ich] < f_low_eff[tt] * mean_nhits[tt]){ 
	nch_low_eff[tt]++;
	if(tt == 0) mch_low_eff_p[ich] = 1;
	if(tt == 1) mch_low_eff_l[ich] = 1;
	if(tt == 2) mch_low_eff_n[ich] = 1;
	if(tt == 3) mch_low_eff_m[ich] = 1;
	get_message(bx_message::log) << "Mch " << ich << ": Low_eff_in_" <<  trg_names[tt] << " " << nhits_per_tt_mch[tt][ich] << " hits in " << nevents_per_tt[tt] << " triggers" << dispatch;
	if(tt == 0) multiplicity_pl_vec[ich].push_back(low_eff_in_pulser  );
	if(tt == 1) multiplicity_pl_vec[ich].push_back(low_eff_in_laser   );
	if(tt == 2) multiplicity_n_vec [ich].push_back(low_eff_in_neutrino);
	if(tt == 3) multiplicity_m_vec [ich].push_back(low_eff_in_muon    );
      }		

      // retriggering?
      /*TH1D* h_dt_vs_mch_singch   = h_dt_vs_mch[tt]->ProjectionY ("h_dt_vs_mch_singch", ich+1, ich+1);
      if((h_dt_vs_mch_singch->Integral () / nhits_per_tt_mch[tt][ich]) > f_retriggering) { 
	nch_retriggering[tt]++;
	if(tt == 0) mch_retriggering_p[ich] = 1;
	if(tt == 1) mch_retriggering_l[ich] = 1;
	if(tt == 2) mch_retriggering_n[ich] = 1;
	if(tt == 3) mch_retriggering_m[ich] = 1;
	get_message(bx_message::log) << "Mch " << ich << ": Retriggering_in_" << trg_names[tt] << " (" << h_dt_vs_mch_singch->Integral () << " second hits from the total of " 
		<< nhits_per_tt_mch[tt][ich] << " in " << nevents_per_tt[tt] << " triggers )" << dispatch;
	if(tt == 0) multiplicity_pl_vec[ich].push_back(retriggering_in_pulser  );
	if(tt == 1) multiplicity_pl_vec[ich].push_back(retriggering_in_laser   );
	if(tt == 2) multiplicity_n_vec [ich].push_back(retriggering_in_neutrino); 
	if(tt == 3) multiplicity_m_vec [ich].push_back(retriggering_in_muon    ); 
      }*/
    } // end of loop on mch
  } // end of loop on tt

  get_message(bx_message::debug) << "end" << dispatch;
  
  for(int tt = 0; tt < 4; tt ++)  get_message(bx_message::log) << "nch_dead_in_"         <<  trg_names[tt] << " " << nch_dead        [tt] << dispatch;
  for(int tt = 0; tt < 4; tt ++)  get_message(bx_message::log) << "nch_low_eff_in_"      <<  trg_names[tt] << " " << nch_low_eff     [tt] << dispatch;
  for(int tt = 0; tt < 4; tt ++)  get_message(bx_message::log) << "nch_hot_in_"          <<  trg_names[tt] << " " << nch_hot         [tt] << dispatch;
//  for(int tt = 0; tt < 4; tt ++)  get_message(bx_message::log) << "nch_retriggering_in_" <<  trg_names[tt] << " " << nch_retriggering[tt] << dispatch;

    // to check how many channels did change status 
  for(int ich = 0; ich < constants::muon::channels; ich ++){
    if (( mch_dead_p        [ich] - prev_mch_dead_p        [ich]) == -1)  mch_dead_p        [constants::muon::channels    ]++;
    if (( mch_dead_p        [ich] - prev_mch_dead_p        [ich]) ==  1)  mch_dead_p        [constants::muon::channels + 1]++;
    if (( mch_dead_l        [ich] - prev_mch_dead_l        [ich]) == -1)  mch_dead_l        [constants::muon::channels    ]++;
    if (( mch_dead_l        [ich] - prev_mch_dead_l        [ich]) ==  1)  mch_dead_l        [constants::muon::channels + 1]++;
    if (( mch_dead_n        [ich] - prev_mch_dead_n        [ich]) == -1)  mch_dead_n        [constants::muon::channels    ]++;
    if (( mch_dead_n        [ich] - prev_mch_dead_n        [ich]) ==  1)  mch_dead_n        [constants::muon::channels + 1]++;
    if (( mch_dead_m        [ich] - prev_mch_dead_m        [ich]) == -1)  mch_dead_m        [constants::muon::channels    ]++;
    if (( mch_dead_m        [ich] - prev_mch_dead_m        [ich]) ==  1)  mch_dead_m        [constants::muon::channels + 1]++;
    if (( mch_low_eff_p     [ich] - prev_mch_low_eff_p     [ich]) == -1)  mch_low_eff_p     [constants::muon::channels    ]++;
    if (( mch_low_eff_p     [ich] - prev_mch_low_eff_p     [ich]) ==  1)  mch_low_eff_p     [constants::muon::channels + 1]++;
    if (( mch_low_eff_l     [ich] - prev_mch_low_eff_l     [ich]) == -1)  mch_low_eff_l     [constants::muon::channels    ]++;
    if (( mch_low_eff_l     [ich] - prev_mch_low_eff_l     [ich]) ==  1)  mch_low_eff_l     [constants::muon::channels + 1]++;
    if (( mch_low_eff_n     [ich] - prev_mch_low_eff_n     [ich]) == -1)  mch_low_eff_n     [constants::muon::channels    ]++;
    if (( mch_low_eff_n     [ich] - prev_mch_low_eff_n     [ich]) ==  1)  mch_low_eff_n     [constants::muon::channels + 1]++;
    if (( mch_low_eff_m     [ich] - prev_mch_low_eff_m     [ich]) == -1)  mch_low_eff_m     [constants::muon::channels    ]++;
    if (( mch_low_eff_m     [ich] - prev_mch_low_eff_m     [ich]) ==  1)  mch_low_eff_m     [constants::muon::channels + 1]++;
    if (( mch_hot_p         [ich] - prev_mch_hot_p         [ich]) == -1)  mch_hot_p         [constants::muon::channels    ]++;
    if (( mch_hot_p         [ich] - prev_mch_hot_p         [ich]) ==  1)  mch_hot_p         [constants::muon::channels + 1]++;
    if (( mch_hot_l         [ich] - prev_mch_hot_l         [ich]) == -1)  mch_hot_l         [constants::muon::channels    ]++;
    if (( mch_hot_l         [ich] - prev_mch_hot_l         [ich]) ==  1)  mch_hot_l         [constants::muon::channels + 1]++;
    if (( mch_hot_n         [ich] - prev_mch_hot_n         [ich]) == -1)  mch_hot_n         [constants::muon::channels    ]++;
    if (( mch_hot_n         [ich] - prev_mch_hot_n         [ich]) ==  1)  mch_hot_n         [constants::muon::channels + 1]++;
    if (( mch_hot_m         [ich] - prev_mch_hot_m         [ich]) == -1)  mch_hot_m         [constants::muon::channels    ]++;
    if (( mch_hot_m         [ich] - prev_mch_hot_m         [ich]) ==  1)  mch_hot_m         [constants::muon::channels + 1]++;
    if (( mch_retriggering_p[ich] - prev_mch_retriggering_p[ich]) == -1)  mch_retriggering_p[constants::muon::channels    ]++;
    if (( mch_retriggering_p[ich] - prev_mch_retriggering_p[ich]) ==  1)  mch_retriggering_p[constants::muon::channels + 1]++;
    if (( mch_retriggering_l[ich] - prev_mch_retriggering_l[ich]) == -1)  mch_retriggering_l[constants::muon::channels    ]++;
    if (( mch_retriggering_l[ich] - prev_mch_retriggering_l[ich]) ==  1)  mch_retriggering_l[constants::muon::channels + 1]++;
    if (( mch_retriggering_n[ich] - prev_mch_retriggering_n[ich]) == -1)  mch_retriggering_n[constants::muon::channels    ]++;
    if (( mch_retriggering_n[ich] - prev_mch_retriggering_n[ich]) ==  1)  mch_retriggering_n[constants::muon::channels + 1]++;
    if (( mch_retriggering_m[ich] - prev_mch_retriggering_m[ich]) == -1)  mch_retriggering_m[constants::muon::channels    ]++;
    if (( mch_retriggering_m[ich] - prev_mch_retriggering_m[ich]) ==  1)  mch_retriggering_m[constants::muon::channels + 1]++;
  }

  get_message(bx_message::log) << "mch_dead_p        : " <<  mch_dead_p        [constants::muon::channels] << " mch lost and " <<  mch_dead_p        [constants::muon::channels + 1] << " mch gained it" << dispatch;
  get_message(bx_message::log) << "mch_dead_l        : " <<  mch_dead_l        [constants::muon::channels] << " mch lost and " <<  mch_dead_l        [constants::muon::channels + 1] << " mch gained it" << dispatch;
  get_message(bx_message::log) << "mch_dead_n        : " <<  mch_dead_n        [constants::muon::channels] << " mch lost and " <<  mch_dead_n        [constants::muon::channels + 1] << " mch gained it" << dispatch;
  get_message(bx_message::log) << "mch_dead_m        : " <<  mch_dead_m        [constants::muon::channels] << " mch lost and " <<  mch_dead_m        [constants::muon::channels + 1] << " mch gained it" << dispatch;
  get_message(bx_message::log) << "mch_low_eff_p     : " <<  mch_low_eff_p     [constants::muon::channels] << " mch lost and " <<  mch_low_eff_p     [constants::muon::channels + 1] << " mch gained it" << dispatch;
  get_message(bx_message::log) << "mch_low_eff_l     : " <<  mch_low_eff_l     [constants::muon::channels] << " mch lost and " <<  mch_low_eff_l     [constants::muon::channels + 1] << " mch gained it" << dispatch;
  get_message(bx_message::log) << "mch_low_eff_n     : " <<  mch_low_eff_n     [constants::muon::channels] << " mch lost and " <<  mch_low_eff_n     [constants::muon::channels + 1] << " mch gained it" << dispatch;
  get_message(bx_message::log) << "mch_low_eff_m     : " <<  mch_low_eff_m     [constants::muon::channels] << " mch lost and " <<  mch_low_eff_m     [constants::muon::channels + 1] << " mch gained it" << dispatch;
  get_message(bx_message::log) << "mch_hot_p         : " <<  mch_hot_p         [constants::muon::channels] << " mch lost and " <<  mch_hot_p         [constants::muon::channels + 1] << " mch gained it" << dispatch;
  get_message(bx_message::log) << "mch_hot_l         : " <<  mch_hot_l         [constants::muon::channels] << " mch lost and " <<  mch_hot_l         [constants::muon::channels + 1] << " mch gained it" << dispatch;
  get_message(bx_message::log) << "mch_hot_n         : " <<  mch_hot_n         [constants::muon::channels] << " mch lost and " <<  mch_hot_n         [constants::muon::channels + 1] << " mch gained it" << dispatch;
  get_message(bx_message::log) << "mch_hot_m         : " <<  mch_hot_m         [constants::muon::channels] << " mch lost and " <<  mch_hot_m         [constants::muon::channels + 1] << " mch gained it" << dispatch;
/*  get_message(bx_message::log) << "mch_retriggering_p: " <<  mch_retriggering_p[constants::muon::channels] << " mch lost and " <<  mch_retriggering_p[constants::muon::channels + 1] << " mch gained it" << dispatch;
  get_message(bx_message::log) << "mch_retriggering_l: " <<  mch_retriggering_l[constants::muon::channels] << " mch lost and " <<  mch_retriggering_l[constants::muon::channels + 1] << " mch gained it" << dispatch;
  get_message(bx_message::log) << "mch_retriggering_n: " <<  mch_retriggering_n[constants::muon::channels] << " mch lost and " <<  mch_retriggering_n[constants::muon::channels + 1] << " mch gained it" << dispatch;
  get_message(bx_message::log) << "mch_retriggering_m: " <<  mch_retriggering_m[constants::muon::channels] << " mch lost and " <<  mch_retriggering_m[constants::muon::channels + 1] << " mch gained it" << dispatch;*/

  // Delete vectors
  for (int i = 0; i < 4; i++) delete [] nhits_per_tt_mch[i];
   
    //check if the situation did not change with respect to the previous run
  int i_nch_status_change   = get_parameter ("nch_status_change").get_int ();

  for(int tt = 0; tt < 4; tt ++){
    if( abs(nch_dead[tt] - prev_nch_dead[tt]) >  i_nch_status_change)
      get_message(bx_message::log) << "nch_dead_in_" << trg_names[tt] << " changed by " << (nch_dead[tt] - prev_nch_dead[tt]) << " from the previous run "  << dispatch;
  }
  for(int tt = 0; tt < 4; tt ++){
    if( abs(nch_low_eff[tt] - prev_nch_low_eff[tt]) >  i_nch_status_change)
      get_message(bx_message::log) << "nch_low_eff_in_" << trg_names[tt] << " changed by " << (nch_low_eff[tt] - prev_nch_low_eff[tt]) << " from the previous run "  << dispatch;
  } 
  for(int tt = 0; tt < 4; tt ++){
    if( abs(nch_hot[tt] - prev_nch_hot[tt]) >  i_nch_status_change)
      get_message(bx_message::log) << "nch_hot_in_" << trg_names[tt] << " changed by " << (nch_hot[tt] - prev_nch_hot[tt]) << " from the previous run "  << dispatch;
  }
/*  for(int tt = 0; tt < 4; tt ++){
    if( abs(nch_retriggering[tt] - prev_nch_retriggering[tt]) >  i_nch_status_change)
      get_message(bx_message::log) << "nch_retriggering_in_" << trg_names[tt] << " changed by " << (nch_retriggering[tt] - prev_nch_retriggering[tt]) << " from the previous run "  << dispatch;
  }*/

    //set visitors and write to the DB
  if (get_parameter ("db_write").get_bool () &&  nch_dead[2] < 40 && mean_nhits[1] > 50 &&  mean_nhits[0] > 100 && module_says_DB_write == 1){ // AAA: ID values, to be tuned for OD
    db_run& run_info = bx_dbi::get ()->get_run ();
    std::vector<std::string> vstr;
    for(int ich = 0; ich < constants::muon::channels; ich++) {
      
      vstr.clear ();
      for (unsigned i = 0; i < multiplicity_pl_vec[ich]      .size (); i++) vstr.push_back (multiplicity_translation_map.rfind (multiplicity_pl_vec    [ich][i])->first);
      if (mean_nhits[2] > 100) //enough neutrino triggers using current run 
	for (unsigned i = 0; i < multiplicity_n_vec[ich]     .size (); i++) vstr.push_back (multiplicity_translation_map.rfind (multiplicity_n_vec     [ich][i])->first);
      else //not enough neutrino triggers using previous run  
	for (unsigned i = 0; i < prev_multiplicity_n_vec[ich].size (); i++) vstr.push_back (multiplicity_translation_map.rfind (prev_multiplicity_n_vec[ich][i])->first);
      if (mean_nhits[3] > 100) //enough muon triggers using current run 
	for (unsigned i = 0; i < multiplicity_m_vec[ich]     .size (); i++) vstr.push_back (multiplicity_translation_map.rfind (multiplicity_m_vec     [ich][i])->first);
      else //not enough neutrino triggers using previous run  
	for (unsigned i = 0; i < prev_multiplicity_m_vec[ich].size (); i++) vstr.push_back (multiplicity_translation_map.rfind (prev_multiplicity_m_vec[ich][i])->first);
      run_info.set_muon_multiplicity  (ich+constants::muon::channel_offset+1, vstr, this);
    }

     //write to the DB
    run_info.write_muon_electronic_channel (true, this);
    get_message(bx_message::log) << "Writing to DB" << dispatch; 
  }
  
}
/*
 * $Log: bx_calib_muon_electronics.cc,v $
 * Revision 1.9  2009/10/26 11:19:37  ddangelo
 * trgtype 'muon_tct' (or 'muon_totc') renamed as 'neutron' after real use
 * trgtype 'muon_mtb' renamed as 'muon' as duplicity with the above one disappeared and 'mtb' is associated to btb_inputs flag
 *
 * Revision 1.8  2008-09-29 17:52:02  ddangelo
 * debugging
 *
 * Revision 1.7  2008-09-23 16:41:11  ddangelo
 * debugging
 *
 * Revision 1.6  2008-08-19 17:53:48  ddangelo
 * debugging
 * some parameter tuning
 *
 * Revision 1.5  2008-08-13 12:25:48  ddangelo
 * added "muon" tt checks
 *
 * Revision 1.4  2008-08-12 17:13:53  ddangelo
 * starting to converge
 *
 * Revision 1.3  2008-08-11 16:51:36  ddangelo
 * cleaned up a bit
 *
 */
