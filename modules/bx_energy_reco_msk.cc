/* BOREXINO Reconstruction program
 *
 * Author: Evgeny Litvinovich <litvinov@lngs.infn.it>
 * Maintainer: Evgeny Litvinovich <litvinov@lngs.infn.it>
 *
 * $Id: 
 *
 * Implementation for bx_energy_reco_msk.hh
*/

#include <algorithm>
#include "bx_energy_reco_msk.hh"
#include "bx_echidna_event.hh"
#include "db_channel.hh"
#include "db_calib.hh"
#include "bx_dbi.hh"
#include "constants.hh"
#include "TMinuit.h"
#include "TString.h"
#include "TMath.h"
#include "barn_interface.hh"
#include "bx_detector.hh"

#include <math.h>
#define pi 3.141592653589793238
#define speed_of_light_m_ns 0.3

namespace { 
  bx_energy_reco_msk *current_module;
  void fcn_e(int32_t &npar, double *gin, double &f, double *x, int32_t iflag) 
            { f = current_module->msk_fcn_e(x); }
};

bx_energy_reco_msk::bx_energy_reco_msk() : bx_base_module("bx_energy_reco_msk", bx_base_module::main_loop) {
  require_event_stage (bx_detector::laben, bx_base_event::reconstructed);
  current_module = this;
}

void bx_energy_reco_msk::begin() {

  msk_nphotons = new TH1F("msk_nphotons", "Nphotons/4#pi", 5000, 0, 200000);
  barn_interface::get()->store (barn_interface::file, msk_nphotons, this);

//  msk_charge = new TH1F("msk_charge", "reconstructed charge", 2000, 0, 10000);
//  barn_interface::get()->store (barn_interface::file, msk_charge, this);

  p_minuit = new TMinuit(4);
  p_minuit->SetPrintLevel(-1);
  p_minuit->Command("SET NOW");  // no warnings

  i4_n_minuit_warns = 0;
  f8_attenuation_length = get_parameter("attenuation_length").get_float();
//  f4_ref_index = get_parameter("refidx").get_float();
  f4_cathode_efficiency = get_parameter("cathode_efficiency").get_float();
//  f8_coef_to_return_npe_cone = get_parameter("coef_to_return_npe_cone").get_float();
//  f8_coef_to_return_npe_nocone = get_parameter("coef_to_return_npe_nocone").get_float();

  f4_pmt_positions = new std::vector<TVector3>(constants::laben::channels);
  b_pmt_cone = new std::vector<bool>(constants::laben::channels);

  const std::vector<int32_t>& disabled_channels = detector_interface::get()->get_disabled_channels();
  for (int32_t k = 0; k < (int32_t)disabled_channels.size(); ++k)  i4_disabled_channels.push_back(disabled_channels[k]);

  f4_ref_index = bx_dbi::get()->get_calib().get_refraction_index_data();
  f4_lg_entry_radius = bx_dbi::get()->get_profile().light_guide_entry_aperture();
  f4_cathode_radius = bx_dbi::get()->get_profile().pmt_chathode_radius();

  for (int32_t i = 0; i < constants::laben::channels; ++i)  {
       const db_channel_laben& channel_info = dynamic_cast<const db_channel_laben&>(bx_dbi::get()->get_channel(i+1));
       (*b_pmt_cone)[i] = channel_info.pmt_has_cone(); 
       if ((*b_pmt_cone)[i])  {
           ((*f4_pmt_positions)[i])(0) = channel_info.pmt_x();
           ((*f4_pmt_positions)[i])(1) = channel_info.pmt_y();
           ((*f4_pmt_positions)[i])(2) = channel_info.pmt_z();
	  }
       else  {  // no cone
           ((*f4_pmt_positions)[i])(0) = 1.03634 * channel_info.pmt_x();
           ((*f4_pmt_positions)[i])(1) = 1.03634 * channel_info.pmt_y();
           ((*f4_pmt_positions)[i])(2) = 1.03634 * channel_info.pmt_z();
          }
      }
}

bx_echidna_event* bx_energy_reco_msk::doit(bx_echidna_event *ev)  {

  double f8_valueX, f8_valueY, f8_valueZ, f8_valueE, f8_error;

  int32_t n_clusters = ev->get_laben().get_nclusters();
  if (n_clusters <= 0)  return ev;

  for (int32_t i = 0; i < n_clusters; ++i)  {

       if (ev->get_laben().get_cluster(i).get_clustered_nhits() < 75)  continue;
//       if (ev->get_laben().get_cluster(i).get_charge() < 100)  continue;

  f4_collected_charge.assign(constants::laben::channels, 0);

  const bx_position& pos = ev->get_laben().get_rec_cluster(i).get_position();
  bx_energy_msk& ene = ev->get_laben().get_cluster(i).get_energy_msk();

  f8_valueX = pos.get_x();
  f8_valueY = pos.get_y();
  f8_valueZ = pos.get_z();

  p_fit_ev = ev;
  i4_fit_cluster = i;

// Start reconstructing energy
  p_minuit->SetFCN(fcn_e);
  p_minuit->DefineParameter(0, "E", 2000., 0.5, 5.0, 60000.0);  // Nph/steradian
  p_minuit->DefineParameter(1, "X", f8_valueX, 0.1, -7.0, 7.0);
  p_minuit->DefineParameter(2, "Y", f8_valueY, 0.1, -7.0, 7.0);
  p_minuit->DefineParameter(3, "Z", f8_valueZ, 0.1, -7.0, 7.0);

  p_minuit->FixParameter(1);
  p_minuit->FixParameter(2);
  p_minuit->FixParameter(3);

  TVector3 minuit_position(f8_valueX, f8_valueY, f8_valueZ);
  TVector3 distance;

/*  double r = minuit_position.Mag();
  if (r > 6)  continue;
*/
  f8_omega.assign(constants::laben::channels, 0);
  f8_attenuation_factor.assign(constants::laben::channels, 0);

  for (int32_t k = 0; k < constants::laben::channels; ++k)  {
       const db_channel& ch_info = bx_dbi::get()->get_channel(k+1);
       if (!ch_info.is_ordinary())  continue;

       int32_t index = ch_info.get_lg()-1;

       if (std::find(i4_disabled_channels.begin(),i4_disabled_channels.end(),index+1) != i4_disabled_channels.end())
       continue;

       distance = (*f4_pmt_positions)[index] - minuit_position;
       double path = distance.Mag();

       if ((*b_pmt_cone)[index])  {
            double costheta = (*f4_pmt_positions)[index] * distance * (1./(6.33*path));
            (f8_omega)[index] = pi*f4_lg_entry_radius*f4_lg_entry_radius*pow(costheta,1.9)/(path*path);
            (f8_attenuation_factor)[index] = exp(-1.*path/f8_attenuation_length);
          }
       else  {  // no cone
            double costheta = (*f4_pmt_positions)[index] * distance * (1./(6.56*path));
            (f8_omega)[index] = pi*f4_cathode_radius*f4_cathode_radius*costheta/(path*path);
            (f8_attenuation_factor)[index] = exp(-1.*path/f8_attenuation_length);
	  }
      }

  p_minuit->Command("MIGRAD");

  if (p_minuit->fCstatu != "CONVERGED ")  {
//      get_message(bx_message::warn) << "event " << ev->get_event_number() << ", Minuit returned status " << p_minuit->fCstatu << dispatch;
      ++i4_n_minuit_warns;
     }

  p_minuit->GetParameter(0, f8_valueE, f8_error);

  p_minuit->Release(1); // Fixed parameters
  p_minuit->Release(2); // need to be released back.
  p_minuit->Release(3); // Otherwise they remain fixed next Minuit call

// I'm ready to return reconstructed energy:
//double charge = 0.;
  for (int32_t kk = 0; kk < constants::laben::channels; ++kk)  {
       const db_channel& ch_info = bx_dbi::get()->get_channel(kk+1);
       if (!ch_info.is_ordinary())  continue;
       int32_t index = ch_info.get_lg()-1;
       if (find(i4_disabled_channels.begin(),i4_disabled_channels.end(),index+1) != i4_disabled_channels.end())
       continue;

//       if ((*b_pmt_cone)[index])  { charge += f8_valueE * f8_coef_to_return_npe_cone; }
//                            else  { charge += f8_valueE * f8_coef_to_return_npe_nocone; }
      }

  ene.f4_nphotons = float(f8_valueE * 12.6); //into 4pi
//  ene.f4_charge = (float)charge;

// Since we decided to have also f4_charge, I'm obliged to do a trick.
// I shall just correct i4_npe with the coeff., followed from the ratio charge/npe,
// __measured__ by the detector within current event:
//  ene.f4_charge = npe *
//  (ev->get_laben().get_cluster(i).get_charge()/ev->get_laben().get_cluster(i).get_npe());

// f4_likelihood; i4_n_iterations

  msk_nphotons->Fill(ene.get_nphotons());
//  msk_charge->Fill(ene.get_charge());

 }

//  ev->get_laben().mark_stage(bx_base_event::reconstructed_msk);

  return ev;
}

void bx_energy_reco_msk::end() {

  if (i4_n_minuit_warns) get_message (bx_message::info) << i4_n_minuit_warns << " events processed with the warns from MINUIT" << dispatch;

  delete f4_pmt_positions;
  delete b_pmt_cone;
  delete p_minuit;
}


double bx_energy_reco_msk::msk_fcn_e (double x[0]) {

  double f8_result = 0.;

  const bx_laben_cluster& cluster = p_fit_ev->get_laben().get_cluster(i4_fit_cluster);

  for (int32_t i = 0; i < cluster.get_clustered_nhits(); ++i)  {
       const bx_laben_clustered_hit& hit = cluster.get_clustered_hit(i);
       const bx_laben_decoded_hit& dhit = hit.get_decoded_hit();
       const db_channel* ch_info = dhit.get_db_channel();

       int32_t index = ch_info->get_lg()-1;

       f4_collected_charge[index] = dhit.get_charge_pe();
      }

  for (int32_t k = 0; k < constants::laben::channels; ++k)  {
       const db_channel& ch_info = bx_dbi::get()->get_channel(k+1);
       if (!ch_info.is_ordinary())  continue;

       int32_t index = ch_info.get_lg()-1;

       if (find(i4_disabled_channels.begin(),i4_disabled_channels.end(),index+1) != i4_disabled_channels.end())
       continue;

       double expected_q = x[0]*f4_cathode_efficiency*(f8_omega)[index]*(f8_attenuation_factor)[index];

       double f8_prob = TMath::Poisson(f4_collected_charge[index], expected_q);
//       double f8_prob = pow(f4_collected_charge[index], expected_q) / exp(expected_q) / lgamma(collected_charge[index] + 1);

//       double lnpoisson = f4_collected_charge[index] * log(expected_q) - expected_q - lgamma(f4_collected_charge[index] + 1.);
//       double f8_prob = exp(lnpoisson);
       f8_result += -2. * log(f8_prob);
      }

  return f8_result;
}

/*
 * $Log: bx_energy_reco_msk.cc,v $
 * Revision 1.9  2011/03/02 13:43:57  litvinov
 * "reconstructed" event_stage is required
 *
 * Revision 1.8  2011-03-02 05:39:21  litvinov
 * using rec_cluster for position
 *
 * Revision 1.7  2011-02-17 13:13:03  litvinov
 * threshold to perform reco: 75 instead of 60 nhits
 *
 * Revision 1.6  2011-02-15 14:52:05  litvinov
 * using position_lngs instead of mi
 *
 * Revision 1.5  2008-11-21 16:54:36  litvinov
 * do not reconstruct events with nhits<60
 *
 * Revision 1.4  2008-10-17 10:07:01  litvinov
 * use position_mi as input instead of dbn
 *
 * Revision 1.3  2008-06-20 16:21:55  razeto
 * Added an include to compile with gcc 4.3
 *
 * Revision 1.2  2007-11-07 19:32:14  litvinov
 * using "dubna" reconstructed position; using "charge" as input instead of "npe"
 *
 */
