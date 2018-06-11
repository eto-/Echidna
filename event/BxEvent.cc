/* BOREXINO Reconstruction program
 *
 * Author: Davide D'Angelo <davide.dangelo@lngs.infn.it>
 * based on original design by Marco Pallavicini
 * Maintainer: Davide D'Angelo <davide.dangelo@lngs.infn.it>
 *
 * $Id: BxEvent.cc,v 1.168 2015/08/26 11:22:08 misiaszek Exp $
 *
 * Implementation of BxEvent 
 *
 */

#include "TVector3.h"
#include "BxEvent.hh"
#include "bx_barn.hh"
#include "constants.hh"
#include <iostream>
#include <string.h>

ClassImp(BxTrigger);
ClassImp(BxLabenRawHit);
ClassImp(BxLabenDecodedHit);
ClassImp(BxLabenClusteredHit);
ClassImp(BxLabenCluster);
ClassImp(BxLabenRecHit);
ClassImp(BxLabenRecCluster);
ClassImp(BxLaben);
ClassImp(BxPosition);
ClassImp(BxEnergy);
ClassImp(BxMuonRawHit);
ClassImp(BxMuonDecodedHit);
ClassImp(BxMuonCluster);
ClassImp(BxMuonClusteredHit);
ClassImp(BxMuon);
ClassImp(BxMcTruthHit);
ClassImp(BxMcTruthFrame);
ClassImp(BxMcTruth);
ClassImp(BxEvent);
ClassImp(BxTrack);
ClassImp(BxTrackByPoints);
ClassImp(BxTrackFitted);
ClassImp(BxDistance);
ClassImp(BxNeutron);
ClassImp(BxNeutronPulse);

// These constructors depend on echidna internal events.
// The dependency is removed in the ROOT shared library with this #ifndef
#ifndef _ECHIDNA_ROOTLIB_
#include "bx_echidna_event.hh"
#include "cmap.hh"
#include "bx_dbi.hh"
#include "db_run.hh"

/*********** BxTrigger *************/
void BxTrigger::operator=(const bx_trigger_event& e) {
  trgtype = e.get_trgtype(); 
  trgtime = e.get_trg_time();
  btb_threshold = e.get_btb_threshold();
  btb_inputs = (uint8_t) e.get_btb_inputs ();
  e.get_gps_time(gpstimes[0], gpstimes[1]);
  timet = e.get_time_t();
}


/*********** BxLaben *************/
BxLabenRawHit::BxLabenRawHit(const bx_laben_raw_hit& h) : lg(h.get_logical_channel()),
							  time1(h.get_time_1()), 
							  time2(h.get_time_2()),
							  gray(h.get_gray_counter()), 
							  base(h.get_base()),
							  peak(h.get_peak()), 
							  order(h.get_order_in_channel()),
							  flags_board(h.get_flags_board ()),
							  flags_ch(h.get_flags_ch ()) {
}

BxLabenDecodedHit::BxLabenDecodedHit(const bx_laben_decoded_hit& h, uint16_t raw_index) : lg(h.get_raw_hit().get_logical_channel()),
												raw_time(h.get_raw_time()),
						//						time_error(h.get_time_error()),
												flag(h.get_flag()),
												order(h.get_order_in_channel()),
						//						d80(h.get_d80()),
												raw_charge(h.get_uncorrected_charge_bin()),
												charge(h.get_charge_pe()),
												charge_mean(h.get_charge_mean_pe()),
//												npe(h.get_charge_npe()),
												raw_index(raw_index),
												num_cluster(0),
												rec_time(0.),
												short_cluster(false) {
}

BxLabenClusteredHit::BxLabenClusteredHit(const int cluster, const bx_laben_clustered_hit& h, uint16_t decoded_index) : charge(h.get_decoded_hit().get_charge_pe()) {
}

BxLabenCluster::BxLabenCluster(const bx_laben_cluster& c, uint16_t decoded_index) : 
											  npmts                ( c.get_npmts                     ()),
											  npmts_conc           ( c.get_npmts_conc                ()),
											  npmts_short          ( c.get_npmts_short               ()),
											  npmts_dt1            ( c.get_npmts_dt1                 ()),
											  npmts_dt2            ( c.get_npmts_dt2                 ()),
											  npmts_pos            ( c.get_npmts_pos                 ()),
											  nhits                ( c.get_clustered_nhits           ()), 
											  nhits_conc           ( c.get_clustered_nhits_conc      ()),
											  nhits_short          ( c.get_clustered_nhits_short     ()),
											  nhits_pos            ( c.get_clustered_nhits_pos       ()),
//											  npe                  ( c.get_npe                       ()),
//								 			  npe_conc             ( c.get_npe_conc                  ()),
											  charge               ( c.get_charge                    ()),
											  charge_conc          ( c.get_charge_conc               ()),
											  charge_short         ( c.get_charge_short              ()),
											  charge_dt1           ( c.get_charge_dt1                ()),
											  charge_dt2           ( c.get_charge_dt2                ()),
											  charge_pos           ( c.get_charge_pos                ()),
											  charge_npmts         ( c.get_charge_npmts              ()),
											  charge_noavg_dt1     ( c.get_charge_noavg_dt1          ()),
											  charge_noavg_dt2     ( c.get_charge_noavg_dt2          ()),
											  charge_noavg         ( c.get_charge_noavg              ()),
											  charge_noavg_short   ( c.get_charge_noavg_short        ()),
											  start_time           ( c.get_start_time                ()),
//											  rough_time           ( c.get_rough_time                ()),
											  mean_time            ( c.get_mean_time                 ()),
											  mean_time_short      ( c.get_mean_time_short           ()),
											  rms_time             ( c.get_rms_time                  ()),
											  rms_time_short       ( c.get_rms_time_short            ()),
											  duration             ( c.get_duration                  ()),
											  duration_short       ( c.get_duration_short            ()),
											  flag                 ( c.get_flag                      ()),
											  is_neutron           ( c.is_neutron                    ()),
											  decoded_index        ( decoded_index                     ),
											  baricenter           ( c.get_baricenter                ()),
											  position_lngs        ( c.get_position_lngs             ()),
                                                                                          position_noavg       ( c.get_position_noavg            ()),
                                                                                          npmt_geo_weight      ( c.get_npmt_geo_weight           ()),
											  npmt_QE_weight       ( c.get_npmt_QE_weight            ()),
                                                                                          npmt_geo_QE_weight   ( c.get_npmt_geo_QE_weight        ()),
                                                                                          charge_geo_weight    ( c.get_charge_geo_weight         ()),
											  charge_QE_weight     ( c.get_charge_QE_weight          ()),						
                                                                                          charge_geo_QE_weight ( c.get_charge_geo_QE_weight      ()){
  peak_times.clear();
  for (int i = 0 ; i < c.get_split_npeaks(); i++) {
    peak_times.push_back(c.get_split_peak(i).get_start_time());
    peak_charges.push_back(c.get_split_peak(i).get_nhits());
  } 
  npeaks = c.get_split_npeaks();
}

BxLabenRecHit::BxLabenRecHit(const int cluster, const bx_laben_rec_hit& h) {
}

//BxLabenRecCluster::BxLabenRecCluster(const bx_laben_ab_mach4_cluster& c) : BxPosition(c.get_position()), 
BxLabenRecCluster::BxLabenRecCluster(const bx_laben_positron_cluster& c) : BxPosition(c.get_position()), 
//								     BxEnergy(c.get_energy()),
								     ns_asymmetry(c.get_ns_asymmetry()),
								     sphere_chi2(c.get_sphere_chi2()), 
								     sphere_lkl(c.get_sphere_lkl()),
								     sphere_rel_var(c.get_sphere_rel_var()),
								     plane_cos(c.get_plane_cos()),
								     plane_chi2(c.get_plane_chi2()),
								     h_plane_chi2(c.get_h_plane_chi2()),
								     quality_flags(c.get_quality_flags()),
								     gatti(c.get_gatti()),
								     lkl(c.get_lkl()),
								     gattic(c.get_gattic()),
								     lklc(c.get_lklc()),
								     rise_time(c.get_rise_time()),
								     rms(c.get_rms()),
								     rms_c11(c.get_rms_c11()),
								     kurtosis(c.get_kurtosis()),
								     kurtosis_c11(c.get_kurtosis_c11()),
								     mlp_ab(c.get_mlp_ab()),
//								     peak_mach4(c.get_peak_mach4()),
//								     mean_mach4(c.get_mean_mach4()),
//								     rms_mach4 (c.get_rms_mach4 ()),
//								     skew_mach4(c.get_skew_mach4()),
//								     kurt_mach4(c.get_kurt_mach4()),
								     gatti_ops_beta(c.get_gatti_ops_beta()),
								     gatti_c11_beta(c.get_gatti_c11_beta()),
								     gatti_ops_nops(c.get_gatti_ops_nops()){
   std::copy(c.get_tailtot  (), c.get_tailtot  () + 10, tailtot   );
   std::copy(c.get_tailtot_ab_mlp  (), c.get_tailtot_ab_mlp  () + 10, tailtot_ab_mlp   );
   std::copy(c.get_tailtot_c11_mva  (), c.get_tailtot_c11_mva  () + 10, tailtot_c11_mva   );
//   std::copy(c.get_tailtot_mach4(), c.get_tailtot_mach4() + 17, tailtot_mach4);
//   std::copy(c.get_gatti_mach4()  , c.get_gatti_mach4()   +  4, gatti_mach4  );
   std::copy(c.get_sh_power (), c.get_sh_power () +  4, sh_power  );
}


BxPosition::BxPosition(const bx_base_position& p) : time (p.get_t()), 
						    x (p.get_x()),
						    y (p.get_y()),
						    z (p.get_z()),
						    dt (p.get_dt()),
						    dx (p.get_dx()),
						    dy (p.get_dy()),
						    dz (p.get_dz()),
						    user (p.get_user()),
						    converged (p.is_converged()),
						    matrix (p.get_matrix ()) {
}
BxEnergy::BxEnergy(const bx_base_energy& e) : nhits    (e.get_nhits()),
					      npe      (e.get_npe()),
					      charge   (e.get_charge()) {
}

BxLaben::BxLaben(const bx_write_opts& opts) : has_raw(opts.laben.raw), 
					      has_decoded(opts.laben.decoded),
					      has_clustered(opts.laben.clustered), 
					      has_rec(opts.laben.rec) {
}

void BxLaben::operator=(const bx_laben_event& e) {
  empty_boards         = e.get_empty_boards   ();
  trigger_time         = e.get_trigger_rawt   ();
  laser_time           = e.get_laser_rawt     ();
//  npe                  = e.get_npe            ();
  npmts                = e.get_npmts          ();
  charge               = e.get_charge         ();
  n_live_pmts          = e.get_n_live_pmts    ();
  n_live_charge        = e.get_n_live_charge  ();
  n_hits_on_empty      = e.get_nhits_on_empty ();
//  cluster_window_limit = e.get_window_limit   ();
  track_energy         = e.get_track_energy   ();
  is_tracked_energy    = e.is_tracked_energy  ();
  track_tof            = e.get_track_tof      ();
  is_tracked_tof       = e.is_tracked_tof     ();
  npmts_win1           = e.get_npmts_win1();
  npmts_win2           = e.get_npmts_win2();
  charge_win1          = e.get_charge_win1();
  charge_win2          = e.get_charge_win2();
  raw_hits      .clear();
  decoded_hits  .clear();
  clusters      .clear();
  clusters_muons.clear();
  clustered_hits.clear();
  rec_clusters  .clear();
  rec_hits      .clear();

  n_raw_hits = e.get_raw_nhits();
  n_raw_hits_fw = e.get_raw_nhits_fw();
  for (int i = 0; i < 8; i++) n_raw_hits_flags[i] = e.get_raw_nhits_flag (bx_laben_raw_hit::flags(i));
  n_invalid_pmts = e.get_invalid_pmts();
  n_invalid_charge = e.get_invalid_charge();
  n_decoded_hits = e.get_decoded_nhits();
  n_clusters = e.get_nclusters();
  n_clusters_muons   = e.get_nclusters_muons();
  n_clusters_found   = e.get_nclusters_found();
  n_clusters_old     = e.get_nclusters_old();
  n_clusters_neutron = e.get_nclusters_neutron();
  n_clustered_hits = 0;
  for (int i = 0 ; i < e.get_nclusters(); i++)
    n_clustered_hits += e.get_cluster(i).get_clustered_nhits();

  std::map<const bx_laben_raw_hit*, uint16_t> raw_map;
  std::map<const bx_laben_decoded_hit*, uint16_t> decoded_map;

  if(has_raw) {
    for (int i = 0 ; i < e.get_raw_nhits(); i++) {
      raw_hits.push_back(e.get_raw_hit(i));
      raw_map[&(e.get_raw_hit(i))] = i;
    }
  }
  
  if(has_decoded) {
    for (int i = 0 ; i < e.get_decoded_nhits(); i++) {
      decoded_hits.push_back(BxLabenDecodedHit(e.get_decoded_hit(i), 
					       raw_map[&(e.get_decoded_hit(i).get_raw_hit())] ));
      decoded_map[&(e.get_decoded_hit(i))] = i;
    }
    for (int i = 0 ; i < e.get_nclusters(); i++) {
      for (int j = 0 ; j < e.get_cluster(i).get_clustered_nhits(); j++) {
          const bx_laben_decoded_hit& dec_hit = e.get_cluster(i).get_clustered_hit(j).get_decoded_hit();
	  uint16_t k = decoded_map[&dec_hit]; 
          decoded_hits[k].num_cluster = i+1;
	  decoded_hits[k].short_cluster = e.get_cluster(i).get_clustered_hit(j).is_short_cluster();
      }
    }
    for (int i = 0 ; i < e.get_nrec_clusters(); i++) {
      for (int j = 0 ; j < e.get_rec_cluster(i).get_rec_nhits(); j++) {
          const bx_laben_decoded_hit& dec_hit = e.get_rec_cluster(i).get_rec_hit(j).get_clustered_hit().get_decoded_hit();
	  uint16_t k = decoded_map[&dec_hit]; 
	  decoded_hits[k].rec_time = e.get_rec_cluster(i).get_rec_hit(j).get_time();
      }
    }
  }

  if(has_clustered > 0) {
    for (int i = 0 ; i < e.get_nclusters(); i++) {
      clusters.push_back(BxLabenCluster(e.get_cluster(i), 
					decoded_map[&(e.get_cluster(i).get_clustered_hit(0).get_decoded_hit())] ));
      if (has_clustered > 1) {
	for (int j = 0 ; j < e.get_cluster(i).get_clustered_nhits(); j++)
	  clustered_hits.push_back(BxLabenClusteredHit(i+1, e.get_cluster(i).get_clustered_hit(j), decoded_map[&(e.get_cluster(i).get_clustered_hit(j).get_decoded_hit())] ));
      }
    } 
    for (int i = 0 ; i < e.get_nclusters_muons(); i++)
      clusters_muons.push_back(BxLabenCluster(e.get_cluster_muon(i), 
					decoded_map[&(e.get_cluster_muon(i).get_clustered_hit(0).get_decoded_hit())] ));
  }

  if(has_rec > 0) {
    for (int i = 0 ; i < e.get_nrec_clusters(); i++) {
      rec_clusters.push_back(BxLabenRecCluster(e.get_positron_cluster(i)));
      if (has_rec > 1) {
	for (int j = 0 ; j < e.get_rec_cluster(i).get_rec_nhits(); j++)	{
	  rec_hits.push_back(BxLabenRecHit(i+1, e.get_rec_cluster(i).get_rec_hit(j)));
	}
      }
    } 
  }
  
}



/*********** BxMuon *************/
BxMuonRawHit::BxMuonRawHit(const bx_muon_raw_hit& h) : mch(h.get_muon_channel()), 
						       lead_time(h.get_lead_time()),
						       trail_time(h.get_trail_time()) {
}

BxMuonDecodedHit::BxMuonDecodedHit(const bx_muon_decoded_hit& h) : run(bx_dbi::get()->get_run().get_number()), 
								   mch(h.get_raw_hit().get_muon_channel()), 
								   time(h.get_time()),
								   charge(h.get_charge()) {
}

BxMuonCluster::BxMuonCluster(const bx_muon_cluster& c) : id(c.get_id()),
							 x(c.get_x()), y(c.get_y()), z(c.get_z()),
						         charge(c.get_charge()),
						         start_time(c.get_start_time()) {
}

BxMuonClusteredHit::BxMuonClusteredHit(const bx_muon_clustered_hit& h) : run(bx_dbi::get()->get_run().get_number()), 
								   	 mch(h.get_decoded_hit().get_raw_hit().get_muon_channel()), 
								         time(h.get_time()),
								         charge(h.get_charge()) {
}

BxMuon::BxMuon(const bx_write_opts& opts) : has_raw(opts.muon.raw), 
					    has_decoded(opts.muon.decoded),
					    has_clustered(opts.muon.clustered) {
}

void BxMuon::operator=(const bx_muon_event& e) {
  raw_hits.clear();
  decoded_hits.clear();
  clusters.clear();
  is_aligned             = e.is_aligned();
  n_raw_hits             = e.get_raw_nhits();
  n_decoded_hits         = e.get_decoded_nhits();
  n_clustered_hits_sss   = e.get_clustered_nhits_sss();
  n_clustered_hits_floor = e.get_clustered_nhits_floor();
  n_clusters             = e.get_nclusters();
  decoded_charge         = e.get_decoded_charge();
  decoded_npmts          = e.get_decoded_npmts();
  has_cluster_sss        = e.has_cluster_sss();
  has_cluster_floor      = e.has_cluster_floor();
  charge_sss             = e.get_charge_sss();
  charge_floor           = e.get_charge_floor();
  npmts                  = e.get_npmts();
  start_time_sss         = e.get_start_time_sss();
  start_time_floor       = e.get_start_time_floor();
  is_tracked             = e.is_tracked ();
  track                  = e.get_track ();

  if(has_raw) {
    for (int i = 0 ; i < e.get_raw_nhits(); i++)
      raw_hits.push_back(e.get_raw_hit(i));
  }

  if(has_decoded) {
    for (int i = 0 ; i < e.get_decoded_nhits(); i++)
      decoded_hits.push_back(e.get_decoded_hit(i));
  }
 
  if(has_clustered > 0) {
    for (int i = 0 ; i < e.get_nclusters(); i++)
      clusters.push_back(e.get_cluster(i));
    /*if (has_clustered > 1) {
      for (int i = 0 ; i < e.get_clustered_nhits(); i++)
        clustered_hits.push_back(e.get_clustered_hit(i));
    } */
 }
}

/********** BxTrack ***************/
void BxTrackByPoints::operator=(const bx_track& tt ) {
  const bx_track_by_points& t = dynamic_cast<const bx_track_by_points&>(tt);
  t1 = t.get_t1();
  t2 = t.get_t2();
  x1 = t.get_x1();
  y1 = t.get_y1();
  z1 = t.get_z1();
  x2 = t.get_x2();
  y2 = t.get_y2();
  z2 = t.get_z2();
  dx1 = t.get_dx1();
  dy1 = t.get_dy1();
  dz1 = t.get_dz1();
  dx2 = t.get_dx2();
  dy2 = t.get_dy2();
  dz2 = t.get_dz2();
  theta   = t.get_theta  ();
  phi     = t.get_phi    ();
  dtheta  = t.get_dtheta ();
  dphi    = t.get_dphi   ();
  impact  = t.get_impact ();
  dimpact = t.get_dimpact();
  idnhits = t.get_laben_normhits();
  error   = t.get_error  ();
  downward = t.is_downward();
}

void BxTrackFitted::operator=( const bx_track& tt ) {
  const bx_track_fitted& t = dynamic_cast<const bx_track_fitted&>(tt);
  alpha       = t.get_alpha       ();
  beta        = t.get_beta        ();
  gamma       = t.get_gamma       ();
  delta       = t.get_delta       ();
  alpha_error = t.get_alpha_error ();
  beta_error  = t.get_beta_error  ();
  gamma_error = t.get_gamma_error ();
  delta_error = t.get_delta_error ();
  chi2        = t.get_chi2        ();
  points      = t.get_points      ();
  theta       = t.get_theta  ();
  phi         = t.get_phi    ();
  dtheta      = t.get_dtheta ();
  dphi        = t.get_dphi   ();
  impact      = t.get_impact ();
  idnhits     = t.get_laben_normhits();
  dimpact     = t.get_dimpact();
  downward    = t.is_downward();
}

/*********** BxMcTruth *************/
BxMcTruthHit::BxMcTruthHit(const Int_t frame, const bx_mctruth_hit& h) : num_frame(frame), 
									 lg(h.get_lg()), 
									 time(h.get_time()) {
}

BxMcTruthDaughter::BxMcTruthDaughter(const Int_t frame, const bx_mctruth_daughter& d) :  num_frame(frame), 
		 									 id(d.get_id()), 
											 pdg(d.get_pdg()),
											 time(d.get_time()),
											 energy(d.get_energy()) {
  for (int i = 0; i < 3; i++ ) {
    position[i]  = d.get_position(i);
    direction[i] = d.get_direction(i);
  }
}

BxMcTruthDeposit::BxMcTruthDeposit(const Int_t frame, const bx_mctruth_deposit& d) : num_frame(frame), 
										 pdg_parent(d.get_pdg_parent()), 
	 									 energy(d.get_energy()) {
  for (int i = 0; i < 3; i++ ) 
    position[i]  = d.get_position(i);
}

BxMcTruthUser::BxMcTruthUser(const Int_t frame, const bx_mctruth_user& u) : num_frame(frame), 
									 int1(u.get_int1()), 
									 int2(u.get_int2()),
									 float1(u.get_float1()),
									 float2(u.get_float2()),
									 double1(u.get_double()) {
}

BxMcTruthFrame::BxMcTruthFrame(const bx_mctruth_frame& f) : file_id        (f.get_file_id        ()),
							    elec_event_time(f.get_elec_event_time()),
                                                            event_id       (f.get_event_id       ()), 
							    n_sequence     (f.get_n_sequence     ()),
							    isotope_coinc  (f.get_isotope_coinc  ()),
							    pdg            (f.get_pdg            ()),
							    time           (f.get_time           ()),
							    energy         (f.get_energy         ()),
							    visible_energy (f.get_visible_energy ()),
							    id_npe         (f.get_id_npe         ()),
							    od_npe         (f.get_od_npe         ()),
							    n_daughters    (f.get_n_daughters    ()), 
							    n_deposits     (f.get_n_deposits     ()),
							    n_users        (f.get_n_users        ()),
  							    n_id_photons   (f.get_n_id_photons   ()),
							    n_od_photons   (f.get_n_od_photons   ()) {
  for (int i = 0; i < 3; i++ ) { 
    position  [i] = f.get_position  (i);
    baricenter[i] = f.get_baricenter(i);
    direction [i] = f.get_direction (i);
  }
}

BxMcTruth::BxMcTruth(const bx_write_opts& opts) : write_flag(opts.mctruth) {
}

void BxMcTruth::operator=(const bx_mctruth_event& e) {
  frames.clear();
  hits_id.clear();
  hits_od.clear();
  daughters.clear();
  deposits.clear();
  users.clear();

//std::cout << " in mctruth operator=\n";

  n_hits_id = 0, n_hits_od = 0, n_daughters = 0, n_deposits = 0, n_users = 0;
  for (uint32_t i = 0 ; i < e.get_nframes(); i++) {
    n_hits_id   += e.get_frame(i).get_id_npe();
    n_hits_od   += e.get_frame(i).get_od_npe();
    n_daughters += e.get_frame(i).get_n_daughters();
    n_deposits  += e.get_frame(i).get_n_deposits(); 
    n_users     += e.get_frame(i).get_n_users();
  }
  n_frames = e.get_nframes();

  if (write_flag > 0) {
    for (uint32_t i = 0 ; i < e.get_nframes(); i++) {
      frames.push_back(e.get_frame(i));
      if (write_flag > 1) {
        for (int j = 0; j < e.get_frame(i).get_id_npe(); j++)  
	  hits_id.push_back(BxMcTruthHit(i+1, e.get_frame(i).get_hit_id(j)));
        for (int j = 0; j < e.get_frame(i).get_od_npe(); j++)  
	  hits_od.push_back(BxMcTruthHit(i+1, e.get_frame(i).get_hit_od(j)));
        for (int j = 0; j < e.get_frame(i).get_n_daughters(); j++)  
	  daughters.push_back(BxMcTruthDaughter(i+1, e.get_frame(i).get_daughter(j)));
        for (int j = 0; j < e.get_frame(i).get_n_deposits(); j++)  
	  deposits.push_back(BxMcTruthDeposit(i+1, e.get_frame(i).get_deposit(j)));
        for (int j = 0; j < e.get_frame(i).get_n_users(); j++)  
	  users.push_back(BxMcTruthUser(i+1, e.get_frame(i).get_user(j)));
      }
    }
  }
}

/********** Neutron ***************/
BxNeutronPulse::BxNeutronPulse(const bx_neutron_pulse& p) : charge(p.get_charge()),
                                                            amplitude(p.get_amplitude()),
                                                            peak_time(p.get_peak_time()),
                                                            rise_time(p.get_rise_time()),
                                                            fall_time(p.get_fall_time()),
                                                            x(p.get_x()),y(p.get_y()),z(p.get_z()),
                                                            dx(p.get_dx()),dy(p.get_dy()),dz(p.get_dz()){
}

void BxNeutron::operator=(const bx_neutron_event& e) {
  is_enabled = e.is_enabled();
  is_associated = e.is_associated();
  n_neutrons=e.get_n_neutrons();
  pulses.clear();
  for (unsigned u=0; u<e.get_neutron_pulses().size(); u++)
    pulses.push_back(e.get_neutron_pulses()[u]);
}



/*********** BxEvent *************/
BxEvent::BxEvent(const bx_write_opts& opts) : trigger(),
					      laben(opts),
					      muon(opts),
					      mctruth(opts),n_m4(0) {
}

void BxEvent::operator=(const bx_echidna_event& ev) { 
  run = ev.get_run_number();
  evnum = ev.get_event_number(); 
  enabled_crates = ev.get_enabled_crates();
  is_tracked_global = ev.is_tracked_global();
  //is_tracked_cmt = ev.is_tracked_cmt();
  trigger = ev.get_trigger();
  laben = ev.get_laben();
  muon = ev.get_muon();
  mctruth = ev.get_mctruth();
  //neutron = ev.get_neutron();
  track_global = ev.get_track_global();
  //track_cmt = ev.get_track_cmt();
}

#endif // of #ifndef _ECHIDNA_ROOTLIB_

Float_t BxTrackByPoints::GetDedx() const {
  Float_t i = GetImpact();
  if (i>6.85) return 0;
  else if (i>4.25) return -idnhits/500./GetPathSSS();
  else return +(idnhits-320.*GetPathBuffer())/500./GetPathIV();
}

Float_t BxTrackFitted::GetDedx() const {
  Float_t i = GetImpact();
  if (i>6.85) return 0;
  else if (i>4.25) return -idnhits/500./GetPathSSS();
  else return +(idnhits-320.*GetPathBuffer())/500./GetPathIV();
}

BxPosition BxTrackByPoints::GetProjection(const BxPosition& p) const {
  TVector3 v_a(x1, y1, z1);
  TVector3 v_ab(x2-x1, y2-y1, z2-z1);
  float l_ab = v_ab.Mag();
  TVector3 v_ac(p.GetX()-x1, p.GetY()-y1, p.GetZ()-z1);
  TVector3 v_pout = v_a + v_ab.Dot(v_ac)/(l_ab*l_ab) * v_ab;

  BxPosition pout(v_pout(0),v_pout(1),v_pout(2));
  return pout;
}

BxPosition BxTrackFitted::GetProjection(const BxPosition& p) const {
  TVector3 v_a(0, alpha, gamma);
  TVector3 v_ab(1, beta, delta);
  float l_ab = v_ab.Mag();
  TVector3 v_ac(p.GetX()-v_a(0), p.GetY()-v_a(1), p.GetZ()-v_a(2));
  TVector3 v_pout = v_a + v_ab.Dot(v_ac)/(l_ab*l_ab) * v_ab;

  BxPosition pout(v_pout(0),v_pout(1),v_pout(2));
  return pout;
}

Float_t BxTrackByPoints::GetDistance(const BxPosition& p) const {
  TVector3 v_in (x1, y1, z1);
  TVector3 v_out(x2, y2, z2);
  TVector3 v_pos(p.GetX(), p.GetY(), p.GetZ());

  return (sin((v_out-v_in).Angle(v_in-v_pos))*(v_in-v_pos).Mag());
}

Float_t BxTrackFitted::GetDistance(const BxPosition& p) const {
  TVector3 v_in (0, alpha, gamma);
  TVector3 v_out(1, alpha+beta, gamma+delta);
  TVector3 v_pos(p.GetX(), p.GetY(), p.GetZ());

  return (sin((v_out-v_in).Angle(v_in-v_pos))*(v_in-v_pos).Mag());
}


BxPhysTags::BxPhysTags (): tags(0), tags_user(0), tags_solar(0), tags_cno(0), coincidence_id(0), coincidence_index(0) {
}


BxFwfdCluster::BxFwfdCluster(Int_t v1, Double_t v2, Float_t v3, Float_t v4, Float_t v5, Float_t v6, Float_t v7, Int_t v8, Float_t v9, Float_t v10, Float_t v11) {
  peak_pos = v1; time_prev = v2; peak_ampl = v3; asum_charge = v4; atten_asum = v5; dsum_charge = v6; dsum_charge_corr = v7; 
  num_ech_cluster = v8; x = v9; y = v10; z = v11; muon_tag = -100; noise_tag = -100; user1 = 0; user2 = 0;
}

BxFwfdCluster::BxFwfdCluster() {peak_pos = 0; time_prev = 0.; peak_ampl = 0.; asum_charge = 0.; atten_asum = 0.; dsum_charge = 0.; dsum_charge_corr = 0.; num_ech_cluster = 0; x = 0.; y = 0.; z = 0.;}

BxFwfd::BxFwfd(Bool_t v1, Bool_t v2, Int_t v3, Int_t v4, uint32_t v5, uint32_t v6, Int_t v7, Int_t v8, Int_t v9, uint8_t v10, Double_t v11, Int_t v12, UInt_t v13, uint8_t* v14) : clusters() {

  is_present = v1; is_odsum_present = v2; n_fwfd_evs = v3; unix_time = v4; gpstimes[0] = v5; gpstimes[1] = v6; run = v7; evnum = v8; evnum_bx = v9; error = v10; raw_time = v11; trgtype = v12; n_clusters = v13;
  for (int i = 0; i < 16; i++)  dcode[i] = v14[i];

}

BxFwfd::BxFwfd() : clusters() {  is_present = false; is_odsum_present = false; n_fwfd_evs = 0; unix_time = 0; gpstimes[0] = 0; gpstimes[1] = 0; run = 0; evnum = 0; evnum_bx = 0; error = 0; raw_time = 0; trgtype = 0; n_clusters = 0;
  for (int i = 0; i < 16; i++)  dcode[i] = 0;}

void BxFwfd::SetCluster(UInt_t n_clusters, UInt_t num_cluster, Int_t v1, Double_t v2, Float_t v3, Float_t v4, Float_t v5, Float_t v6, Float_t v7, Int_t v8, Float_t v9, Float_t v10, Float_t v11) {
      if (clusters.size() != n_clusters)  clusters.resize(n_clusters);
      clusters.at(num_cluster) = BxFwfdCluster(v1,v2,v3,v4,v5,v6,v7,v8,v9,v10,v11);
     }

void BxFwfd::ClearClusters() {
      if (!clusters.empty())  clusters.clear();
}

void BxFwfd::ClearWForms() {
      if (!wform_asum.empty())  wform_asum.clear();
      if (!wform_dsum.empty())  wform_dsum.clear();
      if (!wform_odsum.empty())  wform_odsum.clear();
}


double BxEvent::GetTimeDifference (const uint32_t* prev_gps_times, Double_t prev_laben_trigger_time) const {
  static const long int gray_window = (1 << 16) * 50;
  
  long int gps_dt_s = GetTrigger ().GetGpsTimeSec () - prev_gps_times[0];
  long int gps_dt_ns = GetTrigger ().GetGpsTimeNs () - prev_gps_times[1];
  double dt_gps_us = double(gps_dt_s) * 1e6 + double(gps_dt_ns) * 1e-3;

  if (::fabs (dt_gps_us) > (gray_window * 1e-3 - 20)) { // 20 is to have some margin when dt is to close to gray_window
    return dt_gps_us;
  } else {
    double curr_laben_trigger_time = GetLaben ().GetTriggerTime ();
    bool positive_dt = curr_laben_trigger_time > prev_laben_trigger_time;
    if (::fabs (dt_gps_us) < 1) dt_gps_us = positive_dt ? 1 : -1;
    if (prev_laben_trigger_time > gray_window) prev_laben_trigger_time -= gray_window;
    if (curr_laben_trigger_time > gray_window) curr_laben_trigger_time -= gray_window;

    double dt_laben_simple = curr_laben_trigger_time - prev_laben_trigger_time;
    if (dt_gps_us > 0 && dt_laben_simple < 0) dt_laben_simple += gray_window;
    else if (dt_gps_us < 0 && dt_laben_simple > 0) dt_laben_simple -= gray_window;
    
    if (::fabs (dt_laben_simple) > gray_window) std::cout << "unexpected high value of dt_laben " << dt_laben_simple << " ns" << std::endl;
    return dt_laben_simple * 1e-3;
  }
}

double BxEvent::GetTimeDifference (int current_cluster, const BxEvent& prev_event, int prev_cluster) const {
  if (prev_cluster < 0 || current_cluster < 0) return GetTimeDifference (prev_event);
  return GetTimeDifference (current_cluster, prev_event.GetTrigger ().GetGpsTimes (), prev_event.GetLaben ().GetCluster (prev_cluster).GetStartTime ());
}

double BxEvent::GetTimeDifference (int current_cluster, const uint32_t* prev_gps_times, Double_t prev_laben_cluster_time) const {
  static const long int gray_window = (1 << 16) * 50;
  
  long int gps_dt_s = GetTrigger ().GetGpsTimeSec () - prev_gps_times[0];
  long int gps_dt_ns = GetTrigger ().GetGpsTimeNs () - prev_gps_times[1];
  double dt_gps_us = double(gps_dt_s) * 1e6 + double(gps_dt_ns) * 1e-3;

  if (::fabs (dt_gps_us) > (gray_window * 1e-3 - 20)) { // 20 is to have some margin when dt is to close to gray_window
    return dt_gps_us;
  } else {
    double curr_laben_cluster_time = GetLaben ().GetCluster (current_cluster).GetStartTime ();
    bool positive_dt = curr_laben_cluster_time > prev_laben_cluster_time;
    if (::fabs (dt_gps_us) < 0.1) dt_gps_us = positive_dt ? 1 : -1;
    if (prev_laben_cluster_time > gray_window) prev_laben_cluster_time -= gray_window;
    if (curr_laben_cluster_time > gray_window) curr_laben_cluster_time -= gray_window;

    double dt_laben_simple = curr_laben_cluster_time - prev_laben_cluster_time;
    if (dt_gps_us > 0 && dt_laben_simple < 0) dt_laben_simple += gray_window;
    else if (dt_gps_us < 0 && dt_laben_simple > 0) dt_laben_simple -= gray_window;
    
    if (::fabs (dt_laben_simple) > gray_window) std::cout << "unexpected high value of dt_laben " << dt_laben_simple << " ns" << std::endl;
    return dt_laben_simple * 1e-3;
  }
}


Float_t BxLabenRecCluster::GetTailTot(Int_t tail) const {
  if ( tail < 40 || tail > 130 )
    return -1.;
  int index = (tail - 40)/10;
  return tailtot[index];
}

Float_t BxLabenRecCluster::GetTailTotAbMlp(Int_t tail_index) const {
  if ( tail_index < 0 || tail_index > 9 )
    return -1.;
  return tailtot_ab_mlp[tail_index];
}

Float_t BxLabenRecCluster::GetTailTotC11Mva(Int_t tail_index) const {
  if ( tail_index < 0 || tail_index > 9 )
    return -1.;
  return tailtot_c11_mva[tail_index];
}

/*Float_t BxLabenRecCluster::GetTailTotMach4(Int_t tail) const {
  if ( tail < 30 || tail > 110 )
    return -1.;
  int index = (tail - 30)/5;
  return tailtot_mach4[index];
}*/


static double jd (time_t t) {
  struct tm *s = gmtime (&t);
  double year = s->tm_year + 1900;
  double month = s->tm_mon + 1;
  double day = s->tm_mday;
  double hour = s->tm_hour + s->tm_min / 60.;
  return 367 * year - floor (7. / 4 * (floor ((month + 9) / 12) + year)) + floor (275 * month / 9) + day - 730531.5  + hour / 24; 
}

static const double pi=3.1415926535897932384626433832795029L;
static double _2r (double grad) { return grad / 180 * pi; }
static double _2g (double rad) { return rad * 180 / pi; }

float BxTrigger::GetSunAltitude (time_t t, float& azimuth) {
  const double latitude = 42.421;    // LNGS aerial is at N 42° 25' 15.6"
  const double longitude = 13.51466; // LNGS aerial is at E 13° 30' 52.8"
  double jj = jd (t);
  double T = jj / 36525; /* Julian Century */

  double M = fmod (357.5291 + 35999.0503 * T - 0.0001559 * T * T - 0.00000045 * T * T * T, 360);
  double Lo = fmod (280.46645 + 36000.76983 * T + 0.0003032 * T * T,  360);
  double DL = (1.9146 - 0.004817 * T - 0.000014 * T * T) * sin (_2r (M)) + (0.019993 - 0.000101 * T)* sin (_2r (2 * M)) + 0.00029 * sin (_2r (3 * M));
  double L = Lo + DL;

  double eps = 23.43999 - 0.013 * T;
  double delta=  _2g (asin (sin (_2r (L)) * sin (_2r (eps))));
  double RA = fmod (_2g (atan2 (cos (_2r (eps)) * sin (_2r (L)), cos (_2r (L)))) + 360, 360);

  double GMST = fmod (280.46061837 + 360.98564736629 * jj + 0.000387933 * T * T - T * T * T /38710000 + 360,  360);
  double LMST = GMST + longitude;
  double H = LMST - RA;
  azimuth = fmod (_2g (atan2 (-sin (_2r (H)), cos (_2r (latitude)) * tan (_2r (delta))- sin (_2r (latitude)) * cos (_2r (H)))) + 360,  360);
  float altitude = _2g (asin (sin (_2r (latitude)) * sin (_2r (delta)) + cos (_2r (latitude)) * cos (_2r (delta))*cos(_2r (H))));
  return altitude;
}

time_t BxTrigger::GetSunRise (time_t t) {
  time_t day = t - t % 86400;
  time_t sunrise = day;

  int sign = 1;
  float unused;
  for (int step = 86400 / 2; step > 30; step /= 2) {
    sunrise += sign * step;
    float altitude = GetSunAltitude (sunrise, unused);
    if (::fabs(altitude) < 0.1) break;
    sign = altitude > 0 ? -1 : 1;
  }
  return sunrise;
}

time_t BxTrigger::GetSunSet (time_t t) {
  time_t day = t - t % 86400;
  time_t sunfall = day;

  int sign = 1;
  float unused;
  for (int step = 86400 / 2; step > 30; step /= 2) {
    sunfall += sign * step;
    float altitude = GetSunAltitude (sunfall, unused);
    if (::fabs(altitude) < 0.1) break;
    sign = altitude > 0 ? 1 : -1;
  }
  return sunfall;
}

time_t BxTrigger::GetMidday (time_t t) {
  time_t day = t - t % 86400;
  time_t midday = day;

  int sign = 1;
  float azimuth;
  for (int step = 86400 / 2; step > 30; step /= 2) {
    midday += sign * step;
    GetSunAltitude (midday, azimuth);
    if (::fabs (azimuth - 180) < 0.1) break;
    sign = azimuth > 180 ? -1 : 1;
  }
  return midday;
}

/*
 * $Log: BxEvent.cc,v $
 * Revision 1.168  2015/08/26 11:22:08  misiaszek
 * new rms/kurtosis for c11 psd
 *
 * Revision 1.167  2015/08/25 12:28:14  misiaszek
 * mlp_ab variable added
 *
 * Revision 1.166  2015/08/06 12:25:25  koun
 * added Normalize_geo/QE_Pmts/Charge variables
 *
 * Revision 1.165  2015/08/03 15:42:15  misiaszek
 * cno tags added
 *
 * Revision 1.164  2015/07/29 09:23:12  ilia.drachnev
 * added position_reco_noavg variable
 *
 * Revision 1.163  2015/07/28 09:14:54  misiaszek
 * tailtot for c11/b discrimination added
 *
 * Revision 1.162  2015/07/14 13:55:00  lukyanch
 * Remove inheritance from TObject for BxFwfd and add two user variables to BxFwfdCluster
 *
 * Revision 1.161  2015/07/14 09:54:46  misiaszek
 * charge_noavg_dt1,  charge_noavg_dt2,  charge_noavg , charge_noavg_short added
 *
 * Revision 1.160  2015/07/13 17:38:53  misiaszek
 * var_user variables & solar tags added
 *
 * Revision 1.159  2015/07/13 15:06:48  misiaszek
 * tailtot_ab_mlp for a/b discrimination with MLP algorithm added
 *
 * Revision 1.158  2015/07/02 15:35:25  lukyanch
 * Add muon and noise tags to FADC event
 *
 * Revision 1.157  2015/01/02 02:09:52  misiaszek
 * npmts, nhits & charge postition corrected energy variables added to BxLabenCluster
 *
 * Revision 1.156  2014/12/15 12:10:13  dekor
 * Added charge_win1 charge_win2 stuff
 *
 * Revision 1.155  2014/12/12 20:11:46  misiaszek
 * track_cmt removed
 *
 * Revision 1.154  2014/12/11 21:27:12  wurm
 * added getters and variables to BxEvent that concern C-11 likelihood tagging
 * * for tracks: new distance-to- and projection-on-track getters, dE/dx calculation
 * * counting of neutron-like clusters in tt128, new flag of neutron clusters
 *
 * Revision 1.153  2014/12/09 16:50:20  misiaszek
 * Not used energy variables removed
 *
 * Revision 1.152  2014/12/05 17:25:02  misiaszek
 * position_mi/_dbn/_msk/_mach4/_mach4_fixed and energy_mc/_lik/_msk/_dbn removed
 *
 * Revision 1.151  2014/12/03 17:07:56  misiaszek
 * charge_dt1 and charge_dt2 to clusters added
 *
 * Revision 1.150  2014/11/07 12:54:44  litvinov
 * update in FADC stuff: (1) num_ech_cluster and (2) corrected dsum_charge added in clusters, (3) wform of the OD analog sum and (4) whether it is present or not added in base BxFwfd
 *
 * Revision 1.149  2013/01/29 12:32:29  litvinov
 * reconstructed x,y,z for FADC (spatial reco is being implemented in fadc reader)
 *
 * Revision 1.148  2013-01-22 09:22:17  mosteiro
 * dt1,2 vars calculated for each cluster
 *
 * Revision 1.147  2012-10-30 16:21:48  ddangelo
 * added vectors for npmts calculation in tt64
 * internal and root event with getters and copy.
 *
 * Revision 1.146  2012-10-30 15:39:17  ddangelo
 * fixed typos
 *
 * Revision 1.145  2012-10-22 15:56:03  ddangelo
 * added npmts_dt1, npmts_dt2 to laben (clustered) event.
 * internal and root, with getters and cpy operator.
 *
 * Revision 1.144  2012-10-18 14:15:16  litvinov
 * FADC raw_time is being double, ulong is not enough at 32bits (bxmaster)
 *
 * Revision 1.143  2012-10-17 18:24:57  litvinov
 * FADC raw_time now uint32_t (was UInt_t)
 *
 * Revision 1.142  2011-04-13 10:37:35  ddangelo
 * added variable nclusters_old
 *
 * Revision 1.141  2011-03-23 16:04:48  ddangelo
 * added charge_clean, zeroing, copy. internal and root event
 *
 * Revision 1.140  2011-03-23 15:28:42  ddangelo
 * fwfd ctrs zero variables
 *
 * Revision 1.139  2011-02-22 08:11:27  razeto
 * GPS start from 2000 UTC + added leap seconds + timet
 *
 * Revision 1.138  2011-02-21 10:03:39  davini
 * debug: copy of bx_pid_positron
 *
 * Revision 1.137  2011-02-19 13:29:11  davini
 * bx_pid_positron stuffs
 *
 * Revision 1.136  2011-02-18 18:22:01  ddangelo
 * event support to mach4 a/b discrimination commented out. variable writing with the module commented out.
 * added event support for gatti ops. to be improved.
 *
 * Revision 1.135  2011-02-18 17:10:05  ddangelo
 * major code cleanup: removed fadc throughout the program
 *
 * Revision 1.134  2011-02-18 16:04:22  davini
 * added npmts_400 nhits_400 charge_400 on Oleg's request; added charge_npmts based on Alessandro, Stefano and Livia idea;
 *
 * Revision 1.133  2011-02-11 16:28:46  litvinov
 * FWFW raw_time unsigned instead of signed int
 *
 * Revision 1.132  2010-12-10 10:32:43  litvinov
 * added FADC raw_time. Waveforms are now vectors
 *
 * Revision 1.131  2010-11-29 16:03:01  litvinov
 * added run to BxFwfd; n_clusters was uint8_t, now integer
 *
 * Revision 1.130  2010-08-04 08:46:15  wurm
 * copy phi, theta, imp for fitted track
 *
 * Revision 1.129  2010-08-03 15:59:03  wurm
 * removed commented lines
 *
 * Revision 1.128  2010-08-03 15:58:05  wurm
 * introduced theta, phi, impact variables for fitted tracks
 *
 * Revision 1.127  2010-07-01 18:21:35  razeto
 * Added counter hit flag for new laben fw and nhits_fw from laben boards
 *
 * Revision 1.126  2010-06-30 10:58:02  litvinov
 * added variable n_fwfd_evs to BxFwfd
 *
 * Revision 1.125  2010-06-24 16:38:33  litvinov
 * added constructor with parameters for BxFwfd
 *
 * Revision 1.124  2010-06-17 14:44:22  litvinov
 * more FWFD stuff: default constructor for BxFwfdCluster, dtors, ClearClusters()
 *
 * Revision 1.123  2010-06-16 12:54:20  litvinov
 * added FWFD stuff. Commit authorized by ddangelo
 *
 * Revision 1.122  2010-05-27 12:49:35  ddangelo
 * debugging
 *
 * Revision 1.121  2010-05-26 15:12:35  ddangelo
 * last commit on c12 brought to c13:
 * bug fix in BxLaben::operator= rec_clusters not always aligned to clusters.
 *
 * Revision 1.120  2010-05-21 15:55:14  ddangelo
 * different things on muon tracks
 *
 * 1.a) old laben_track renamed as laben_track_energy
 * new laben_track_tof added
 *
 * 1.b) (global) track renemed as track_global at base event level
 * track_cmt added at base event level (track by points)
 *
 * 1) all getters updated/integrated
 * is_tracked variable updated/integrated accordingly. inizialization.
 * job ported to root event as well. copy done.
 * friendship with old/new module updated
 *
 * 2) bxtrack_by_points class:
 * - theta, phi and impact added as variables.
 * - errors added on all of the above.
 * - error code variable requested by cmt tracker added
 *
 * Revision 1.119  2010-03-10 21:33:53  razeto
 * Added n_m4
 *
 * Revision 1.118  2009-11-18 11:43:30  ddangelo
 * added rms_time, duration and npmts short version for back compatibility.
 * npmt variables renamed to npmts also in internal event.
 * mach4_n1700 renamed as mach4_fixed throughout the event classes.
 *
 * Revision 1.117  2009-11-10 11:37:40  razeto
 * Fixed rec_time assignement bug
 *
 * Revision 1.116  2009-10-26 19:17:22  ddangelo
 * added bx_position_mach4_n1700 class (internal and root event, copy and getters)
 * in (base)postion class variable likelihood renamed as user. parallel getters to use it as likelihood or refraction index
 *
 * Revision 1.115  2009-10-23 16:42:22  ddangelo
 * cluster_neutrons renamed as clusters_muons
 *
 * Revision 1.114  2009-10-23 14:00:03  koshio
 * Add the lngs postion reconstruction
 *
 * Revision 1.113  2009-10-23 09:07:15  ddangelo
 * added a laben cluster vector for parallel neutron clustering in the muon gate. empty for non-muon events.
 *
 * Revision 1.112  2009-10-22 16:06:45  ddangelo
 * in bx_laben_cluster added nhits_bkg for evaluate of bkg contribution in case of piled-up events.
 * internal and root event, copy and getters
 *
 * Revision 1.111  2009-10-08 15:45:59  ddangelo
 * implemented pid event variables removal and addition.
 * Internal and root, inizialization, copy and getters.
 *
 * Revision 1.110  2009-10-06 17:49:58  ddangelo
 * 1) laben decoded event: added n_invalid_charge (internal and root event, copy and getters)
 * 2) laben raw hit (root): getter IsInvalid now reads dedicated bit 0x80;
 * 3) NormalizePmts() and NormalizeCharge() methods now correct for invalid channels as well.
 *
 * Revision 1.109  2009-10-06 13:36:17  ddangelo
 * laben decoded event: added variable n_invalid_pmts to account for FE||FF condition only on ordinary && !disabled channels
 * internal and root event, getters and copy
 *
 * Revision 1.108  2009-10-05 15:18:16  ddangelo
 * added npmts at laben decoded event level (inner and outer, copy and getters);
 *
 * Revision 1.107  2009-10-05 14:43:07  ddangelo
 * added invalid flag getter in bx laben raw hit
 * corrected for correct array size raw_flags in raw event
 * added 7 getters for individual flag counters
 *
 * Revision 1.106  2009-09-18 22:20:12  ddangelo
 * added sphere_rel_var to laben shaped cluster. Internal and root event, inizialization and copy.
 *
 * Revision 1.105  2009-09-17 14:28:51  ddangelo
 * in laben clustered/decoded (internal/root) hit added the flag short_cluster to say if hit belonged to the old (c11) cluster
 *
 * Revision 1.104  2009-09-17 13:57:46  ddangelo
 * in laben cluster added nhits, charge and mean_time '_short' for old clustering values
 *
 * Revision 1.103  2009-09-16 15:44:28  ddangelo
 * npe (and npe_conc) removed from BxLaben, BxLabenCluster and BxLabenDecodedHit classes (root event only)
 *
 * Revision 1.102  2009-09-16 15:36:52  ddangelo
 * BxLabenCluster::end_time renamed as duration
 * neutron, tags and raw_index branches disabled
 *
 * Revision 1.101  2009-09-16 12:36:51  ddangelo
 * in BxLabenDecodedHit num_cluster and rec_time inizialized to 0 (it stays for non clustered hits)
 * num_cluster numbering moved to 1-based
 *
 * Revision 1.100  2009-08-31 10:53:56  ddangelo
 * laben decoded hits: added num cluster and rec time
 *
 * Revision 1.99  2009-07-31 15:39:49  ddangelo
 * debugging the work of the [previous commit
 *
 * Revision 1.98  2009-07-30 16:20:44  ddangelo
 * - removed IsTracked() getters from BxTrack classes (root event)
 * + added is_tracked variable and getter in BxEvent and BxLaben classes (root event)
 * + fixed the copy of is_tracked variable in BxMuon class (root event)
 * + added is_tracked variable and getter for bx_echidna event class (internal event). to be filled by global tracker
 *
 * Revision 1.97  2009-07-22 10:41:22  ddangelo
 * in position classes, both internal and root event:
 * - n_iterations removed
 * + matrix and converged variables added with getters.
 *
 * Revision 1.96  2009-07-17 17:46:17  ddangelo
 * In root event only:
 * charge (and charge_mean for now) copied also to clustered hits.
 * d80 and time_error commented out from decoded hits
 *
 * Revision 1.95  2009-07-17 15:39:51  ddangelo
 * laben cluster:
 * + added npmts_thresh, nhits_thresh and charge_thresh, to be computed with >0.2pe hits only (internal and root event).
 * - removed cluster rough time (internal and root).
 * - removed npe and npe_conc (root only)
 *
 * Revision 1.94  2009-07-17 14:44:08  ddangelo
 * ppc0 cpu time reintroduced in trg event (for cmp with GPS)
 *
 * Revision 1.93  2009-07-17 13:08:34  ddangelo
 * alpha beta varaibles reorganized, rise_time added
 *
 * Revision 1.92  2009-07-16 15:50:31  ddangelo
 * copyng data implemented
 *
 * Revision 1.91  2009-07-16 15:17:57  ddangelo
 * mach a/b ported to root event
 * debugging tailtot getter for m4
 * other debugging
 *
 * Revision 1.90  2009-07-16 10:17:42  ddangelo
 * infrastructure for m4 position reco (patch by steve&ale)
 *
 * Revision 1.89  2009-06-04 09:53:12  razeto
 * Fixed a small bug in DT calculation
 *
 * Revision 1.88  2009-04-20 13:50:40  ddangelo
 * added errors on rec coordinates and time in class bx_position and BxPosition. initialization and copy included.
 *
 * Revision 1.87  2009-04-15 17:12:38  ddangelo
 * n_point changed into a bitfield, internal and root event (BxTrackFitted class)
 *
 * Revision 1.86  2009-03-31 17:11:06  ddangelo
 * debugging
 *
 * Revision 1.85  2008-12-15 12:12:44  razeto
 * Added rms_time to cluster (to improve PID)
 *
 * Revision 1.84  2008-12-10 11:40:24  razeto
 * Added mean charge and peak from tt1
 *
 * Revision 1.83  2008-12-02 13:53:32  ddangelo
 * added n_points to fitted track class and getters. internal and root event
 *
 * Revision 1.82  2008-11-11 17:21:45  razeto
 * Reorganized tags
 *
 * Revision 1.81  2008-10-01 16:17:04  ddangelo
 * removed is_pointlike variables and related stuff (bx_filters will perform this task)
 *
 * Revision 1.80  2008-08-26 15:30:30  ddangelo
 * added neutron enabled and association flag, muon aligned flags
 *
 * Revision 1.79  2008-08-21 21:47:27  razeto
 * added midday
 *
 * Revision 1.78  2008-08-21 15:56:55  razeto
 * Reorganization and upgrades in sun code
 *
 * Revision 1.77  2008-08-18 12:17:03  razeto
 * New coordinates for GPS aerial
 *
 * Revision 1.76  2008-08-08 12:48:03  razeto
 * New tag format redesign completed
 *
 * Revision 1.75  2008-08-05 16:30:52  ddangelo
 * added fall time, some variables renamed
 *
 * Revision 1.74  2008-08-05 15:30:50  ddangelo
 * added chi2
 *
 * Revision 1.73  2008-08-05 14:16:59  razeto
 * Init tags accordling to new tags specs
 *
 * Revision 1.72  2008-07-29 18:28:19  razeto
 * Azimut and latitude experimental calculation added
 *
 * Revision 1.71  2008-07-17 15:55:00  ddangelo
 * added track direction, removed unused getter
 *
 * Revision 1.70  2008-07-11 17:06:27  ddangelo
 * added classes for neutron system (code by S. Davini)
 *
 * Revision 1.69  2008-07-11 14:30:18  ddangelo
 * fixed math in BxDistance class
 *
 * Revision 1.68  2008-05-10 23:34:05  ddangelo
 * BxDistance class defined and implemented
 * GetDistance and GetDisatnceError implemented in BxTrackFitted
 * some more work on track and position classes
 *
 * Revision 1.67  2008-05-09 15:46:23  ddangelo
 * muon track implemented as "by_points" and "fitted" classes inheriting from common interface (pure virtual class bx_track)
 * laben and muon host an istance "by_points"
 * echidna event hosts an istance "fitted"
 * structure duplicated at root level and filled
 * GetDistance duplicated in BxPosition class
 *
 * Revision 1.66  2008-02-28 15:48:21  ddangelo
 * debugging
 *
 * Revision 1.65  2008-02-27 11:48:30  ddangelo
 * debugging
 *
 * Revision 1.64  2008-02-26 18:29:25  ddangelo
 * added is_pointlike/tracklike variable and getters in rec (shaped) cluster
 * both inner and root event
 *
 * Revision 1.63  2008-02-26 17:27:44  ddangelo
 * added
 * BxTrack::GetDistance(BxPosition&)
 * BxPosition::GetDistance(BxTrack&)
 * P.S.: thanks to yura
 *
 * Revision 1.62  2008-02-21 15:58:31  ddangelo
 * added muon track as a separate class.
 * Currently used twice: laben and muon (later also globally fit).
 * both inner and roor event
 * added tracked stage in base event enum
 * bx_muon_rec_event renamed as bx_muon_tracked_event
 *
 * Revision 1.61  2007-12-20 18:41:29  ddangelo
 * handling of individual nhits for sss and floor
 * filling some varaibles left over before.
 * some more debugging
 *
 * Revision 1.60  2007-12-10 15:45:02  ddangelo
 * debugging
 *
 * Revision 1.59  2007-12-07 14:09:47  ddangelo
 * added n_live_charge (internal and root), getters and normalize method
 * renamed a variable
 *
 * Revision 1.58  2007-12-06 16:48:09  ddangelo
 * upgraded to new muon clustering.
 * muon tracking side debugged and improved.
 *
 * Revision 1.57  2007-11-26 14:08:11  ddangelo
 * muon event writing commented out (tmp)
 *
 * Revision 1.56  2007-11-14 19:27:03  ddangelo
 * quiter
 *
 * Revision 1.55  2007-11-14 19:02:11  ddangelo
 * added tracking variables to root event
 *
 * Revision 1.54  2007-11-14 17:07:55  ddangelo
 * new muon clustering variables.
 * indipendent up/floor clustered hits vectors
 * (internal and root event)
 * filling and inizialization. tested.
 *
 * Revision 1.53  2007-11-12 15:31:12  razeto
 * Added cluster flags and removed broad variable (davide auth)
 *
 * Revision 1.52  2007-11-06 14:45:10  ddangelo
 * debugging mctruth
 *
 * Revision 1.51  2007-11-05 23:39:02  razeto
 * Added new hits_on_empty variable to laben data
 *
 * Revision 1.50  2007-10-31 17:12:53  razeto
 * Added a new variable for clustering
 *
 * Revision 1.49  2007-10-31 15:42:13  ddangelo
 * added quenched energy (mctruth)
 * disk format, internal and root event
 *
 * Revision 1.48  2007-10-30 18:40:41  ddangelo
 * added multi condition flag to laben decoded hit
 * getter, smart getters, copy.
 *
 * Revision 1.47  2007-10-30 18:12:27  ddangelo
 * added empty_boards to root event too.
 * getters
 *
 * Revision 1.46  2007-10-30 17:06:04  ddangelo
 * warning removal
 *
 * Revision 1.45  2007-10-30 15:45:02  ddangelo
 * added # of laben cluster found by algorythm (different from saved one for high multiplicity events)
 * added end hit time of a cluster
 * internal and root event. getters, copy, initialization, etc...
 *
 * Revision 1.44  2007-10-29 17:20:44  ddangelo
 * adding nphotons for msk energy reco
 *
 * Revision 1.43  2007-10-29 16:30:54  ddangelo
 * added file_id in mctruth
 * added energy in msk reco
 *
 * Revision 1.42  2007-10-26 12:04:51  ddangelo
 * supporting the splitting of msk position and energy reco in 2 modules.
 * Both internal and root event
 *
 * Revision 1.41  2007-10-26 09:22:01  ddangelo
 * supporting separation of msk position and energy reco in 2 modules
 *
 * Revision 1.40  2007-10-25 15:46:03  ddangelo
 * added a bit field to shaped event, meaning to be assigned. smart getters to be added at that time. internal and external event.
 * all variables zeroed in ctor for shaped event.
 *
 * Revision 1.39  2007-10-25 14:43:22  ddangelo
 * added a variable to rec clsuter
 *
 * Revision 1.38  2007-10-11 11:23:52  ddangelo
 * new mctruth format (internal AND root event)
 *
 * Revision 1.37  2007-09-28 17:07:55  razeto
 * Fixed behaviour for dt close to the 3.2ms edge
 *
 * Revision 1.36  2007-09-13 06:01:02  razeto
 * const attribute for GetTimeDiff family
 *
 * Revision 1.35  2007-08-08 21:04:17  razeto
 * GetTimeDiff behaves good even if there are no clusters
 *
 * Revision 1.34  2007-06-06 10:51:38  razeto
 * Init tags to zero
 *
 * Revision 1.33  2007-06-01 15:57:00  ddangelo
 * added the enabled crates word to root event.
 * One plain getter in echidna event
 *
 * Revision 1.32  2007-06-01 14:56:47  ddangelo
 * removed a debug printout
 *
 * Revision 1.31  2007-05-25 15:56:49  ddangelo
 * added npe and charge for laben decode event. Internal and root.
 *
 * Revision 1.30  2007-05-25 14:37:34  ddangelo
 * added management of btb flags.
 * Redesigned internal low level access.
 *
 * Revision 1.29  2007-05-09 12:25:02  razeto
 * BxLabenRecCluster::GetTailTot moved outside ifdef (non maintanier urgent commit)
 *
 * Revision 1.28  2007-05-07 16:41:42  ddangelo
 * n_live_channels renamed as n_live_pmts as requested by ale.
 * getters and normalize improved
 *
 * Revision 1.27  2007-05-04 16:33:24  pallas
 * Changed variables for alpha beta
 * Now tailtot is an array of 10 values
 *
 * Revision 1.26  2007-05-03 17:25:39  ddangelo
 * added n_live_channels for run-time information. added relative getter.
 * added functions Normalize() for common data types.
 * Both internal and root event
 *
 * Revision 1.25  2007-04-27 14:23:29  pallas
 * Small changes to alpha beta variables
 *
 * Revision 1.24  2007-03-30 12:42:41  razeto
 * rec clusters can be 0, if bx_position_reco is disabled
 *
 * Revision 1.23  2007-03-28 17:40:23  ddangelo
 * introduced dbn "energy" class.
 * introduced muon decoded charge/npmts.
 * completed the job on Laben Rec Cluster.
 *
 * Revision 1.22  2007-03-27 18:01:29  ddangelo
 * energy class upgraded (inner event only for now, tmp)
 *
 * Revision 1.21  2007-03-27 15:18:20  ddangelo
 * variables npe_conc, charge_conc, nhits_conc added to laben cluster
 * f4_pe renamed as f4_charge in bx_laben_event.hh
 * decoded_charge and decoded_npmts added tu muon event
 *
 * Revision 1.20  2007-03-23 19:47:52  ddangelo
 * debugging
 *
 * Revision 1.19  2007-03-22 19:36:57  ddangelo
 * implemented muon db channel
 *
 * Revision 1.18  2007-03-22 16:09:35  ddangelo
 * added muon clustered level
 *
 * Revision 1.17  2007-03-15 19:17:19  ddangelo
 * pid event removed.
 * laben event upgraded with classes: bx_laben_shaped_cluster and bx_laben_ab_cluster
 * bx_laben_rec_cluster is now a parent class for the 2 new ones.
 * BxEvent modified accordingly: BxLabenRecHit and BxLabenRecCluster added.
 * BxPidEvent removed.
 *
 * Revision 1.16  2007-03-15 10:43:58  razeto
 * Calculate time difference between clusters
 *
 * Revision 1.15  2007-02-21 18:49:44  ddangelo
 * fixed a few types
 *
 * Revision 1.14  2007/02/21 15:30:30  ddangelo
 * fixed the size of muon decoded hits
 *
 * Revision 1.13  2006/12/14 15:01:02  ddangelo
 * added BTB trhreshold in root file
 *
 * Revision 1.12  2006/10/23 15:34:34  ddangelo
 * applied ale's patch to inlcude laben integrity flags
 *
 * Revision 1.11  2006/10/13 15:31:18  razeto
 * Added variables from raw data
 *
 * Revision 1.10  2006-08-24 18:03:40  ddangelo
 * changed a getter to return a ref rather then ptr.
 *
 * Revision 1.9  2006/08/21 16:52:53  ddangelo
 * added decoded index to clustered hit and plain getter.
 * Computation from map in BxLaben::operator=.
 * More fancy getters to come.
 *
 * Revision 1.8  2006/08/21 11:12:56  razeto
 * Updated to use new bx_barn
 *
 * Revision 1.7  2006/05/15 09:50:03  ddangelo
 * added 5 variables for improved fadc daq
 *
 * Revision 1.6  2006/05/15 09:23:06  ddangelo
 * added class BxPhysTags to be filled by physics tools
 * Removed IgnoreObjectStreamer() calls from all constructors.
 *
 * Revision 1.5  2006/05/08 17:31:33  razeto
 * Added db_channel patch for fadc (sent to the mailing list)
 *
 * Revision 1.4  2006/02/02 09:59:21  razeto
 * Added iostream include
 *
 * Revision 1.3  2006/01/09 16:18:49  razeto
 * Moved db_barn getters to .cc, to avoid db_barn.hh inclusion in BxEven.hh
 *
 * Revision 1.2  2005/12/18 14:35:54  razeto
 * Implememted GetTimeDifferece method as asked in the echidna meeting
 *
 * Revision 1.1  2005/12/12 19:37:57  razeto
 * Moved BxEvent from root to event
 *
 * Revision 1.44  2005/12/12 19:08:09  razeto
 * Split position and energy reco results (auth from maintainers)
 *
 * Revision 1.43  2005/12/06 16:15:04  ddangelo
 * added GPS time
 *
 * Revision 1.42  2005/12/02 18:13:45  misiaszek
 *
 * initialization for lg and run variables added for decoded/clustered laben hits
 *
 * Revision 1.41  2005/11/18 16:42:55  razeto
 * Added new integer charge variable to decoded hit and cluster (auth from davide)
 *
 * Revision 1.40  2005/10/04 20:20:36  razeto
 * Added iterations value during minimization
 *
 * Revision 1.39  2005/09/26 15:56:20  razeto
 * Added rough time from clustering (auth from davide)
 *
 * Revision 1.38  2005/09/26 12:37:57  razeto
 * Fixed a bug and introduced likelihood in BxPosition
 *
 * Revision 1.37  2005/09/20 17:17:12  razeto
 * Fixed pid shape event (auth from mantainer)
 *
 * Revision 1.36  2005/09/19 13:03:45  ddangelo
 * added peaks from splitting.
 * A tmp size is also added to account for ROOT mishandling of vector<float>.
 *
 * Revision 1.35  2005/08/22 11:26:47  ddangelo
 * added Pid classes
 *
 * Revision 1.34  2005/07/27 16:58:57  ddangelo
 * added flag for "broad" (i.e. spaparanzeted) clusters
 *
 * Revision 1.33  2005/07/11 17:16:59  ddangelo
 * removed global event
 * added pid event (partially)
 * untested
 *
 * Revision 1.32  2005/07/06 16:36:41  ddangelo
 * 'Bool_t is_out_of_gate' re-introduced in class BxLabenDecodedHit
 *
 * Revision 1.31  2005/06/20 16:46:46  ddangelo
 * added bx_position_reco_dbn
 *
 * Revision 1.30  2005/04/21 16:02:28  ddangelo
 * - removed out_of_gate variable from laben decode hit since no longer used.
 * - fixed a small bug in filling the n_decoded_hits and the the n_clusters variables in laben event.
 *
 * Revision 1.29  2005/03/14 19:34:54  ddangelo
 * added indexes for backtracing laben hits
 *
 * Revision 1.28  2005/03/14 13:45:22  ddangelo
 * added "order" in laben raw, decoded and clustered stages
 *
 * Revision 1.27  2005/03/01 15:19:10  razeto
 * Merged with cycle_2
 *
 * Revision 1.26  2005/02/10 16:59:36  ddangelo
 * added a variable with size for every vector in every class and relative getters.
 * added some comments
 * removed compilation warnings
 *
 * Revision 1.25  2004/12/22 17:02:48  ddangelo
 * minor improvements
 *
 * Revision 1.24  2004/12/15 18:18:18  ddangelo
 * added position reco results. First draft, to be checked.
 *
 * Revision 1.23.2.3  2005/02/04 10:06:16  ddangelo
 * added elec event time to mctruth
 *
 * Revision 1.23.2.2  2005/02/01 16:51:37  ddangelo
 * nhits in BxLabenCluster is now filled
 *
 * Revision 1.23.2.1  2004/12/06 11:45:44  ddangelo
 * some warnings removed
 *
 * Revision 1.23  2004/12/03 11:59:47  ddangelo
 * added dsum waveform
 * added raw channel waveform (very haevy, disabled by default)
 * Both temporarily use fixed length C-style arrays. To be improved sometime.
 *
 * Revision 1.22  2004/12/02 16:44:09  ddangelo
 * added digital sum ampl+charge+peak
 * waveform still missing
 *
 * Revision 1.21  2004/12/01 15:13:41  ddangelo
 * added classes for fadc event.
 * added a few vairiables.
 * Work in progress, more stuff to come.
 *
 * Revision 1.20  2004/11/26 15:25:11  razeto
 * Thanks to Davide who so kindly underlined a syntactical error, I fixed the "Maintainer" word
 *
 * Revision 1.19  2004/11/26 14:06:23  razeto
 * Added Mantainer field
 *
 * Revision 1.18  2004/11/24 21:08:21  ddangelo
 * added lists of hits outside clusters/fragments
 * Use of std::vector introduced instead of TClonesArray
 * Required ROOT v4.00/??
 *
 * Revision 1.17  2004/10/06 13:49:35  ddangelo
 * fixed a memory leak in BxLabenCluster
 *
 * Revision 1.16  2004/09/23 16:20:39  ddangelo
 * added TClonesArray with hits in cluster class.
 * Splitting doesn't work yet, but file is written correctly.
 *
 * Revision 1.15  2004/09/22 13:25:39  ddangelo
 * updated bx_reco_event to bx_echidna_event
 *
 * Revision 1.14  2004/09/22 11:33:21  ddangelo
 * fixed a typo.
 * added a delete.
 *
 * Revision 1.13  2004/08/31 13:28:41  ddangelo
 * added charge to decoded hit
 *
 * Revision 1.12  2004/07/26 17:36:33  ddangelo
 * fBits and fUniqueID suppressed in ROOT file.
 *
 * Revision 1.11  2004/07/25 16:41:43  ddangelo
 * added logical channel to laben decoded hits.
 *
 * Revision 1.10  2004/07/25 16:29:31  ddangelo
 * some development, not much active at the moment
 *
 * Revision 1.9  2004/07/13 14:50:48  ddangelo
 * added BxClusteredHit. Currently used due to the TClonesArray problem.
 *
 * Revision 1.8  2004/07/13 13:37:10  ddangelo
 * added McTruth and McTruthFrame to root event.
 * McTruthHit commented out for the moment, due to ROOT problems.
 * To be debugged.
 *
 * Revision 1.7  2004/07/07 15:45:26  ddangelo
 * added BxLabenCluster.
 * Some minor debugging.
 *
 * Revision 1.6  2004/06/25 14:46:39  ddangelo
 * fixed a bug.
 * removed a comment.
 *
 * Revision 1.5  2004/06/22 13:02:52  ddangelo
 * added laben laser and trigger time in BxLaben obj
 *
 * Revision 1.4  2004/06/07 18:59:02  ddangelo
 * fixed a mistyping
 *
 * Revision 1.3  2004/06/07 17:14:03  ddangelo
 * conditional compile macros introduced.
 * Laben raw and decoded hit introduced.
 * Muon raw and decoded hit implemented with 2 different classes.
 *
 * Revision 1.2  2004/06/03 15:00:46  ddangelo
 * (Copy) constructor replaced by operator=() for every class.
 * Still in a tmp state.
 *
 * Revision 1.1  2004/05/30 11:54:48  ddangelo
 * A first working version of root file classes.
 * Not many physical variables yet;
 * Global Event still commented.
 * names match ROOT standards not echidna ones.
 * Makefile updated (file names).
 *
 *
 */
