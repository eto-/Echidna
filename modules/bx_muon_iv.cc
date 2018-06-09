/* BOREXINO Reconstruction program
 *
 * Author: Livia Ludhova<Livia.Ludhova@mi.infn.it>
 * 
 * Maintainer: Livia Ludhova<Livia.Ludhova@mi.infn.it>
 *
 */

#include "bx_muon_iv.hh"
#include "bx_echidna_event.hh"
#include "barn_interface.hh"
#include "db_run.hh"
#include "db_profile.hh"
#include "db_channel.hh"
#include "bx_dbi.hh"
#include "TH1F.h"
#include "TH2F.h"

// ctor
bx_muon_iv::bx_muon_iv (): bx_base_module("bx_muon_iv", bx_base_module::main_loop) {
   require_event_stage (bx_detector::laben, bx_base_event::clustered);
   require_trigger_type (bx_trigger_event::neutrino);
}


  // BEGIN
void bx_muon_iv::begin () {
  npe_satur_thresh = get_parameter ("npe_satur_thresh").get_int ();
  hits_satur_thresh = get_parameter ("hits_satur_thresh").get_float ();
  pmt_hitted_thresh = get_parameter ("pmt_hitted_thresh").get_float ();
  pmt_satur_thresh = get_parameter ("pmt_satur_thresh").get_float ();
  mean_time_thresh = get_parameter ("mean_time_thresh").get_float ();
  notcone_cone_thresh = get_parameter ("notcone_cone_thresh").get_float ();

  hit_lg_junk = new TH1F ("hit_lg_junk","hit_lg_junk", constants::laben::channels, 1,constants::laben::channels + 1);
  satur_hit_lg_junk = new TH1F ("satur_hit_lg_junk","satur_hit_lg_junk", constants::laben::channels, 1,constants::laben::channels +1);

  notcone_cone = new TH1F ("notcone_cone","notcone_cone", 150,0.,1.5);


  barn_interface::get ()->store (barn_interface::junk, hit_lg_junk, this);
  barn_interface::get ()->store (barn_interface::junk, satur_hit_lg_junk, this);
  barn_interface::get ()->store (barn_interface::file, notcone_cone, this);
  
  N_ordinary_lg = 0. ; 
  for(int ch = 1; ch < (constants::laben::channels + 1); ch ++) {
    if ( bx_dbi::get ()->get_channel (ch).is_ordinary () ) N_ordinary_lg ++;  
  }
}


  //DOIT
bx_echidna_event* bx_muon_iv::doit (bx_echidna_event *ev) {
 
  const bx_laben_event& er = ev->get_laben ();
  double ev_number = ev->get_event_number ();
  

  // Check if data are present, otherwise return
  if (!er.get_decoded_nhits()) 
    return ev;
  
  if (er.get_nclusters () > 0){ // we do not care 0-cluster events here
      //we loop through all clusters
    for (int cluster = 0; cluster < er.get_nclusters (); cluster ++) {
      int muon_suspect = 0;  
      double N_satur_hits = 0.;
      double N_hits = er.get_cluster (cluster).get_clustered_nhits ();
      double mean_time = er.get_cluster (cluster).get_mean_time ();
      double N_cone_hits = 0;
      double N_notcone_hits = 0;

        //we loop through all hits in the cluster
      for (int hit = 0; hit < N_hits;hit++){
	double npe = er.get_cluster (cluster).get_clustered_hit (hit).get_decoded_hit ().get_charge_npe ();
	double lg = er.get_cluster (cluster).get_clustered_hit (hit).get_decoded_hit ().get_raw_hit ().get_logical_channel ();
	hit_lg_junk-> Fill(lg);
	if(npe > npe_satur_thresh){
	  satur_hit_lg_junk-> Fill(lg);
	  N_satur_hits ++;
	}

	const db_channel_laben &ch_info = dynamic_cast<const db_channel_laben& >(bx_dbi::get ()->get_channel (int(lg)));
	if ( ch_info.pmt_has_cone () ) N_cone_hits ++;          
	if ( !(ch_info.pmt_has_cone ())) N_notcone_hits ++;          

	//if ( bx_dbi::get ()->get_channel (lg).pmt_has_cone () ) N_cone_hits ++;          
	//if ( !(bx_dbi::get ()->get_channel (lg).pmt_has_cone ()) ) N_notcone_hits ++;          
	
      }
      
      //check notcone/cone ratio
      notcone_cone->Fill (N_notcone_hits / N_cone_hits);
      if ((N_notcone_hits / N_cone_hits) > notcone_cone_thresh){
	std::cout << "high notcone/cone  " <<  N_notcone_hits / N_cone_hits <<  ", EvNum " <<  ev_number << std::endl;
        muon_suspect ++;
      }

        //check if the fraction of saturated hits is too high
      if (N_satur_hits/N_hits >  hits_satur_thresh) {
	get_message (bx_message::log) << "high fraction of saturated hits " << N_satur_hits/N_hits << ", EvNum" <<  ev_number << dispatch;
	std::cout << "high fraction of saturated hits " << N_satur_hits/N_hits << ", Nhits " << N_hits << " , EvNum " <<  ev_number << std::endl;
        muon_suspect ++;
      }

        //check if too many lg  hit or saturated
      double N_hit_lg = 0. ; 
      double N_satur_lg = 0. ; 
      for(int ch = 1; ch < (constants::laben::channels + 1); ch ++) {
	if ( (hit_lg_junk->GetBinContent (ch)) > 0 )N_hit_lg ++;
	if ( (satur_hit_lg_junk->GetBinContent (ch)) > 0 ) N_satur_lg ++;
      }
      
      if (N_hit_lg/N_ordinary_lg > pmt_hitted_thresh){
	get_message (bx_message::log) << "high fraction of hit PMT: " << N_hit_lg/N_ordinary_lg << " ,Ev_Num " << ev_number << dispatch;
	std::cout << "high fraction of hit PMT: " << N_hit_lg/N_ordinary_lg << "Ev_Num " << ev_number << std::endl;
	muon_suspect ++;
      }
      if (N_satur_lg/N_ordinary_lg > pmt_satur_thresh){
	get_message (bx_message::log) << "high fraction of satur: " << N_satur_lg/N_ordinary_lg << ", Ev_Num " << ev_number << dispatch;
	std::cout << "high fraction of satur PMT:" << N_satur_lg/N_ordinary_lg << ", Ev_Num " << ev_number << std::endl;
	muon_suspect ++;
      }
 
      //check if the cluster mean time is above the threshold
      if (mean_time > mean_time_thresh){
	get_message (bx_message::log) << "high cluster mean time: " << mean_time << " ns, Ev_Num " << ev_number << dispatch;
	std::cout << "high cluster mean time: " << mean_time << " ns, Ev_Num " << ev_number << std::endl;
	muon_suspect ++;
      }

      //      if(muon_suspect > 0){
      //      if(  N_notcone_hits / N_cone_hits > 1.0 ) {
      	std::cout << "Ev_Num " << ev_number << " muon suspected" << std::endl;
      	std::cout << "cluster mean time: " << mean_time << " ns" << std::endl; 
      	std::cout << "fraction of hit PMT: " << N_hit_lg/N_ordinary_lg << std::endl;
      	std::cout << "fraction of satur PMT:" << N_satur_lg/N_ordinary_lg << std::endl;
      	std::cout << "fraction of saturated hits " << N_satur_hits/N_hits << std::endl;
      	std::cout << "N_hits " << N_hits << std::endl;
              std::cout << "notcone/cone  " <<  N_notcone_hits / N_cone_hits << std::endl;
      	std::cout << "cluster " << cluster << " while event has " << er.get_nclusters () << " clusters" << std::endl;
      	std::cout << "N_ordinary_lg " << N_ordinary_lg << std::endl;
      	std::cout << "N_hit_lg " << N_hit_lg << std::endl;
      	std::cout << "N_satur_lg " << N_satur_lg << std::endl;
      	std::cout << "*********************************" << std::endl;
      
	//
      satur_hit_lg_junk->Reset ();
      hit_lg_junk->Reset ();

   } //end of the loop in clusters  
  }
  return ev;
}

 //END
void bx_muon_iv::end () {
}

