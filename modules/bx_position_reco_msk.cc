/* BOREXINO Reconstruction program
 *
 * Author: Sergey Sukhotin <sukhotin@in2p3.fr>
 *         Evgeny Litvinovich <litvinov@lngs.infn.it>
 * Maintainer: Evgeny Litvinovich <litvinov@lngs.infn.it>
 *
 * $Id: bx_position_reco_msk.cc,v 1.23 2007/11/11 13:59:47 litvinov Exp $
 *
 * Implementation for bx_position_reco_msk.hh
*/

#include "bx_position_reco_msk.hh"
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
  bx_position_reco_msk *current_module;
  void fcn_t(int32_t &npar, double *gin, double &f, double *x, int32_t iflag) 
            { f = current_module->msk_fcn_t(x); }
};

bx_position_reco_msk::bx_position_reco_msk() : bx_base_module("bx_position_reco_msk", bx_base_module::main_loop) {
  require_event_stage (bx_detector::laben, bx_base_event::clustered);
  current_module = this;
}

void bx_position_reco_msk::begin() {

  msk_x = new TH1F("msk_x", "X coordinate of the event", 2000, -10, 10);
  msk_y = new TH1F("msk_y", "Y coordinate of the event", 2000, -10, 10);
  msk_z = new TH1F("msk_z", "Z coordinate of the event", 2000, -10, 10);
  msk_t = new TH1F("msk_t", "time-of-flight by Minuit", 200, 0, 100);

  msk_r = new TH1F("msk_r", "reconstructed R", 1400, 0, 14);
  
  barn_interface::get()->store(barn_interface::file, msk_x, this);
  barn_interface::get()->store(barn_interface::file, msk_y, this);
  barn_interface::get()->store(barn_interface::file, msk_z, this);
  barn_interface::get()->store(barn_interface::file, msk_t, this);

  barn_interface::get()->store(barn_interface::file, msk_r, this);

  p_minuit = new TMinuit(4);
  p_minuit->SetPrintLevel(-1);
  p_minuit->Command("SET NOW");  // no warnings

  i4_n_minuit_warns = 0;
  f8_landau_mpv = get_parameter("landau_mpv").get_float();
  f8_landau_sigma = get_parameter("landau_sigma").get_float();
  f8_time_cut = get_parameter("time_cut").get_float();
//  f4_ref_index = get_parameter("refidx").get_float();

  f4_pmt_positions = new std::vector<TVector3>(constants::laben::channels);
  b_pmt_cone = new std::vector<bool>(constants::laben::channels);

  const std::vector<int32_t>& disabled_channels = detector_interface::get()->get_disabled_channels();
  for (int32_t k = 0; k < (int32_t)disabled_channels.size(); ++k)  i4_disabled_channels.push_back(disabled_channels[k]);

  f4_ref_index = bx_dbi::get()->get_calib().get_refraction_index_data();
  
//  f4_ref_index = 1.53;

  for (int32_t i = 0; i < constants::laben::channels; ++i)  {
       const db_channel_laben& channel_info = dynamic_cast<const db_channel_laben&>(bx_dbi::get()->get_channel(i+1));
       (*b_pmt_cone)[i] = channel_info.pmt_has_cone(); 
       if ((*b_pmt_cone)[i])  {
           ((*f4_pmt_positions)[i])(0) = channel_info.pmt_x();
           ((*f4_pmt_positions)[i])(1) = channel_info.pmt_y();
           ((*f4_pmt_positions)[i])(2) = channel_info.pmt_z();
	  }
       else  {  // no cone
           ((*f4_pmt_positions)[i])(0) = 1.03634 * channel_info.pmt_x(); //1.03634
           ((*f4_pmt_positions)[i])(1) = 1.03634 * channel_info.pmt_y();
           ((*f4_pmt_positions)[i])(2) = 1.03634 * channel_info.pmt_z();
          }
      }
}

bx_echidna_event* bx_position_reco_msk::doit(bx_echidna_event *ev)  {

  double f8_valueX, f8_valueY, f8_valueZ, f8_valueT, f8_error;

  for (int32_t i = 0; i < ev->get_laben().get_nclusters(); ++i)  {

//    const bx_baricenter& b = ev->get_laben().get_cluster(i).get_baricenter();
//    float r = b.get_r();
//    float t0 = (6.85 - r)*f4_ref_index/speed_of_light_m_ns;

  bx_position_msk& pos = ev->get_laben().get_cluster(i).get_position_msk();

  p_fit_ev = ev;
  i4_fit_cluster = i;

// Start reconstructing position
  p_minuit->SetFCN(fcn_t);
  p_minuit->DefineParameter(0, "X", 0, 0.1, -10.0, 10.0);
  p_minuit->DefineParameter(1, "Y", 0, 0.1, -10.0, 10.0);
  p_minuit->DefineParameter(2, "Z", 0, 0.1, -10.0, 10.0);
  p_minuit->DefineParameter(3, "T", 35, 0.2, 1.0, 100.0);
  p_minuit->Command("MIGRAD");

  if (p_minuit->fCstatu != "CONVERGED ")  {
//      get_message(bx_message::warn) << "event " << ev->get_event_number() << ", Minuit returned status " << p_minuit->fCstatu << dispatch;
      ++i4_n_minuit_warns;
     }

  p_minuit->GetParameter(0, f8_valueX, f8_error);
  pos.f4_x = f8_valueX;
  p_minuit->GetParameter(1, f8_valueY, f8_error);
  pos.f4_y = f8_valueY;
  p_minuit->GetParameter(2, f8_valueZ, f8_error);
  pos.f4_z = f8_valueZ;
  p_minuit->GetParameter(3, f8_valueT, f8_error);
  pos.f4_t = f8_valueT;

// f4_likelihood; i4_n_iterations

  msk_x->Fill(f8_valueX); msk_y->Fill(f8_valueY); msk_z->Fill(f8_valueZ);
  msk_t->Fill(f8_valueT);
  
  msk_r->Fill(pow(f8_valueX*f8_valueX+f8_valueY*f8_valueY+f8_valueZ*f8_valueZ,0.5));

 }

  ev->get_laben().mark_stage(bx_base_event::reconstructed_msk);

  return ev;
}

void bx_position_reco_msk::end() {

  if (i4_n_minuit_warns) get_message (bx_message::info) << i4_n_minuit_warns << " events processed with the warnings from MINUIT" << dispatch;

  delete f4_pmt_positions;
  delete b_pmt_cone;
  delete p_minuit;
}

double bx_position_reco_msk::msk_fcn_t (double *x) {

  TVector3 minuit_position(x);
  TVector3 distance;
  double f8_result = 0.;

  const bx_laben_cluster& cluster = p_fit_ev->get_laben().get_cluster(i4_fit_cluster);

  for (int32_t i = 0; i < cluster.get_clustered_nhits(); ++i)  {
       const bx_laben_clustered_hit& hit = cluster.get_clustered_hit(i);
       const bx_laben_decoded_hit& dhit = hit.get_decoded_hit();
       const db_channel* ch_info = dhit.get_db_channel();

       int32_t index = ch_info->get_lg()-1;

       double time = hit.get_time();

       if ( time < f8_time_cut )  {
            distance = (*f4_pmt_positions)[index] - minuit_position;
            double path = distance.Mag();

            if ((*b_pmt_cone)[index])  {
                double costheta = (*f4_pmt_positions)[index] * distance * (1./(6.33*path));

                // G4 fit for the cone (Igor)
                if (costheta < 0.92)  { path += 0.2898; } // +=0.2898
                else  { path += 3.134 - costheta*5.621 + costheta*costheta*2.748; }

                double delay = time + x[3] - ((f4_ref_index*path)/speed_of_light_m_ns);
                double pdf_value = TMath::Landau(delay, f8_landau_mpv, f8_landau_sigma);
                f8_result += -2. * log(pdf_value);
               }
            else  {  // no cone

                double delay = time + x[3] - ((f4_ref_index*path)/speed_of_light_m_ns);
                double pdf_value = TMath::Landau(delay, f8_landau_mpv, f8_landau_sigma);
                f8_result += -2. * log(pdf_value);
               }
          }
      }

  return f8_result;
}


/*
 * $Log: bx_position_reco_msk.cc,v $
 * Revision 1.23  2007/11/11 13:59:47  litvinov
 * larger minuit limits
 *
 * Revision 1.22  2007-10-25 16:20:01  litvinov
 * decoupling "Moscow" energy and position reconstruction onto 2 independent modules
 *
 * Revision 1.21  2007-06-15 11:48:28  litvinov
 * reading refractive index from DB instead of echidna.cfg
 *
 * Revision 1.20  2007-05-25 16:59:59  litvinov
 * minuit warnings suppressed (Ale's request)
 *
 * Revision 1.19  2007-05-24 09:20:14  litvinov
 * getting status of Migrad convergency and putting it into log if problems
 *
 * Revision 1.18  2007-05-18 17:13:16  litvinov
 * extending minuit limits for the parameters
 *
 * Revision 1.17  2007-05-05 16:27:42  litvinov
 * Now the module returns reconstructed energy in terms of i4_npe and f4_charge.
 * i4_nhits still missed as IMHO useless
 *
 * Revision 1.16  2007-05-02 16:28:26  litvinov
 * change getter name from get_bad_channels to get_disabled_channels
 * according to last modification in bx_detector
 *
 * Revision 1.15  2007-05-02 11:41:44  litvinov
 * optimization for faster calculations
 *
 * Revision 1.14  2007-04-29 15:10:19  litvinov
 * moving parameters used by module into echidna.cfg
 *
 * Revision 1.13  2007-04-28 17:42:28  litvinov
 * Fixing the params x,y,z after time-of-flight minimization
 * converged and energy reco started
 *
 * Revision 1.12  2007-04-27 18:17:30  litvinov
 * two fcn-functions instead of one;
 * coming back po Poissonian statistics;
 * preparation to return npe/nhits/charge
 *
 * Revision 1.11  2007-04-23 11:16:28  litvinov
 * check for bad channels via detector_interface
 *
 * Revision 1.10  2007-03-27 18:02:07  ddangelo
 * tmp. writing to event commented out.
 *
 * Revision 1.9  2007-03-08 16:10:58  razeto
 * Add include to allow compiling with newest root (v5.15)
 *
 * Revision 1.8  2006-11-17 17:51:29  litvinov
 * Module has been fully rewritten. No more groups of PMTs again.
 * Many new features. Work is in progress.
 *
 * Revision 1.7  2006-08-21 11:17:25  razeto
 * Updated to new barn_interface
 *
 * Revision 1.6  2006/05/08 17:31:31  razeto
 * Added db_channel patch for fadc (sent to the mailing list)
 *
 * Revision 1.5  2006/01/11 12:01:15  litvinov
 * Implementation of geometrical sharing of the PMTs to finite number of groups.
 * This approx. 6-7 times reduces machine time.
 *
 * Revision 1.4  2006/01/09 16:21:03  razeto
 * Updated to the new root_barn target
 *
 * Revision 1.3  2005/10/16 16:25:34  litvinov
 * updated to call the Poisson predefined in ROOT
 *
 * Revision 1.2  2004/12/15 16:34:41  litvinov
 * updated according to Davide's changes of the event
 *
 * Revision 1.1  2004/12/13 14:04:39  litvinov
 * Moscow event's position & energy reconstruction module is implemented
 *
 *
 */
