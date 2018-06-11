/* BOREXINO Reconstruction program
 *
 * Author: Alessandro Razeto <Alessandro.Razeto@ge.infn.it>
 * Maintainer: Alessandro Razeto <Alessandro.Razeto@ge.infn.it>
 *
 * $Id: bx_position_reco_dbn.cc,v 1.19 2009/10/26 19:17:55 ddangelo Exp $
 *
 * Implementation of bx_baricentrator
 *
 */
#include "bx_position_reco_dbn.hh"
#include "messenger.hh"
#include "bx_dbi.hh"
#include "db_profile.hh"
#include "db_channel.hh"
#include "bx_echidna_event.hh"
#include "barn_interface.hh"
#include "TMinuit.h"

#include <math.h>

namespace { // MINUIT ROOT workaround
  bx_position_reco_dbn *current_module;
  void fcn (int32_t &npar, double *gin, double &f, double *x, int32_t iflag) { f = current_module->my_fcn (npar, x, gin, iflag); }
};

const uint16_t bx_position_reco_dbn::n_shell;
const float bx_position_reco_dbn::shell_external_radius[n_shell] = { 0.7, 1.5, 2.5, 3.5, 100 };  // last shell have a dummy external radius

// ctor
bx_position_reco_dbn::bx_position_reco_dbn (): bx_base_module("bx_position_reco_dbn", bx_base_module::main_loop), i4_minuit_error_count(0) {
  require_event_stage (bx_detector::laben, bx_base_event::clustered);

  current_module = this; // see minuit root workaround

}

// module interface
void bx_position_reco_dbn::begin () {
  f8_first_hits_gate = get_parameter ("first_hits_gate").get_float ();
  b_use_grad = get_parameter ("use_grad").get_bool ();
  const vdt::vdt_vector& gauss_sigma_v = get_parameter ("gauss_sigma").get_vector ();
  const vdt::vdt_vector& t0_shift_v = get_parameter ("t0_shift").get_vector ();
  if (gauss_sigma_v.size () != n_shell || t0_shift_v.size () != n_shell)
    get_message (bx_message::critic) << "gauss_sigma and t0_shift vectors must have size = " << n_shell << dispatch;
  for (int32_t i = 0; i < n_shell; i++) {
    f4_gauss_sigma[i] = gauss_sigma_v[i].get_float ();
    f4_t0_shift[i] = t0_shift_v[i].get_float ();
  }
  float ref_index = get_parameter ("refidx").get_float ();
  if (ref_index > 0) path.set_refraction_index (ref_index);

  f4_sphere_radius = bx_dbi::get ()->get_profile ().light_collector_radius ();
  
    // Initialize minuit
  p_minuit = new TMinuit (4);
  p_minuit->SetPrintLevel (-1);
  p_minuit->Command ("SET NOW");
  p_minuit->SetFCN (fcn);
  if (b_use_grad) p_minuit->Command ("SET GRAD 1");

  p_iconv = new TH2F ("dbn_iconv", "Iconv vs Ierr flags for dbn position reco", 4, 0, 4, 13, 0, 13);
  barn_interface::get ()->store (barn_interface::file, p_iconv, this);
}

bx_echidna_event* bx_position_reco_dbn::doit (bx_echidna_event *ev) {
    // Loop on every cluster
  for (int32_t i = 0; i < ev->get_laben ().get_nclusters (); i++) {
      // Get cluster reference
    const bx_baricenter& b = ev->get_laben ().get_cluster (i).get_baricenter ();

      // to be used in my_fcn (see minuit root workaround)
    p_fit_ev = ev;
    i4_fit_cluster = i;
    //i4_n_iterations = 0;
    
      // Search for the radius shell to use and set pdf sigma
    float r = b.get_r ();
    int32_t shell = n_shell - 1; // skip most external shell
    for (; shell > 0; shell--) if (shell_external_radius[shell - 1] < r) break;
    f4_current_gauss_sigma_square = f4_gauss_sigma[shell] * f4_gauss_sigma[shell];
    float t0 = (f4_sphere_radius - r) / path.get_c_medium () + f4_t0_shift[shell];
    
      // Define parameters
    p_minuit->DefineParameter (0, "Event Time", t0, 1, -50, 50);
    p_minuit->DefineParameter (1, "X", b.get_x (), 0.10, -6.0, 6.0);
    p_minuit->DefineParameter (2, "Y", b.get_y (), 0.10, -6.0, 6.0);
    p_minuit->DefineParameter (3, "Z", b.get_z (), 0.10, -6.0, 6.0);

      // MINIMIZE
    int32_t ierr = p_minuit->Command ("MIGRAD");
    ierr = p_minuit->Command ("HESSE");
    if (!ierr) {
      //get_message (bx_message::log) << "minimization failed on event " << ev->get_event_number () << dispatch;
      i4_minuit_error_count ++;
    }

      // Get results
    double par[4], junk_error;
    p_minuit->GetParameter (0, par[0], junk_error);
    p_minuit->GetParameter (1, par[1], junk_error);
    p_minuit->GetParameter (2, par[2], junk_error);
    p_minuit->GetParameter (3, par[3], junk_error);
    double likelihood = my_fcn (4, par, 0, 0);

    bx_position_dbn& p = ev->get_laben ().get_cluster (i).get_position_dbn ();
    p.f4_t = par[0]; p.f4_x = par[1]; p.f4_y = par[2]; p.f4_z = par[3]; p.f4_user = likelihood; //p.i4_n_iterations = i4_n_iterations;

      // Get ICONV
    int32_t nvpar,nparx,iconv;
    double amin,edm,errdef;
    p_minuit->mnstat(amin,edm,errdef,nvpar,nparx,iconv);

    p_iconv->Fill (iconv, ierr);
  }

  ev->get_laben ().mark_stage (bx_base_event::reconstructed_dbn);
  return ev;
}

void bx_position_reco_dbn::end () {
  if (i4_minuit_error_count) get_message (bx_message::log) << "minimization failed " << i4_minuit_error_count << " times" << dispatch;
  delete p_minuit;
}


double bx_position_reco_dbn::my_fcn (int32_t npar, double *x, double *grad, int32_t iflag) {
  // x[0] is t
  
  const bx_laben_cluster& cluster = p_fit_ev->get_laben().get_cluster (i4_fit_cluster);
  bool store_gradients = (iflag == 2);
  double likelihood = 0;
  float X_derivative[3];

    // Set to 0 the gradient since after the += operator will be used
  if (store_gradients) std::fill_n (grad, 4, 0.);

    // Loop on the cluster hits
  for (int32_t i = 0; i < cluster.get_clustered_nhits (); i++) {
    const bx_laben_clustered_hit& hit = cluster.get_clustered_hit (i);
    if (hit.get_time () > f8_first_hits_gate) continue;

      // Calculate delta as time_of_hit - time_of_flight
    path.init (x + 1, hit.get_decoded_hit ().get_db_channel ());
    double hit_time = hit.get_time () + x[0];
    double time_of_flight = path.get_time ();
    double delta = hit_time - time_of_flight;

      // Apply PDF 
    if (delta > 0) {
      likelihood += delta;
      if (store_gradients) {
	  // Calculate and assign t and X derivative
	path.get_time_derivative (X_derivative);
	grad[0] += 1;
	for (int32_t j = 0; j < 3; j++) grad [j + 1] -= X_derivative[j];
      }
    } else {  
      likelihood += delta * delta / f4_current_gauss_sigma_square;
      if (store_gradients) {
	  // Calculate and assign t and X derivative
	path.get_time_derivative (X_derivative);
	grad[0] += 2 * delta / f4_current_gauss_sigma_square;
	for (int32_t j = 0; j < 3; j++) grad [j + 1] -= X_derivative[j] * 2 * delta / f4_current_gauss_sigma_square;
      }
    }
  } 

  //i4_n_iterations ++;
  return likelihood;
}

/*
 * $Log: bx_position_reco_dbn.cc,v $
 * Revision 1.19  2009/10/26 19:17:55  ddangelo
 * in (base)postion class variable likelihood renamed as user. modules modified accordingly
 *
 * Revision 1.18  2009-07-22 10:41:53  ddangelo
 * n_iterations use commented to allow compilation
 *
 * Revision 1.17  2007-07-08 09:53:31  razeto
 * Depend only on clustered event but neutrino trigger
 *
 * Revision 1.16  2007-05-26 19:54:43  razeto
 * Be quiter
 *
 * Revision 1.15  2006/08/21 11:17:25  razeto
 * Updated to new barn_interface
 *
 * Revision 1.14  2006/01/02 21:24:39  razeto
 * Changed test target to file target for db_barn
 *
 * Revision 1.13  2005/10/30 11:51:23  razeto
 * Upgrade to have a better starting t
 *
 * Revision 1.12  2005/10/30 11:25:53  razeto
 * Added refraction index overriding in flight_path
 *
 * Revision 1.11  2005/10/13 13:42:07  razeto
 * Added shell dependent pdf
 *
 * Revision 1.10  2005/10/12 14:44:26  razeto
 * Added derivative calculation, still in testing
 *
 * Revision 1.9  2005/10/06 21:30:16  razeto
 * Working toward analitic derivative
 *
 * Revision 1.8  2005/10/04 20:24:12  razeto
 * Added iterations value during minimization
 *
 * Revision 1.7  2005/09/26 12:38:46  razeto
 * Store likelihood value in event
 *
 * Revision 1.6  2005/09/22 11:54:21  razeto
 * Now the pdf can be slightly modified with 2 parameters
 *
 * Revision 1.5  2005/07/07 15:14:09  razeto
 * Reduced the minimization steps and added a debugging histogram
 *
 * Revision 1.4  2005/06/30 14:42:02  razeto
 * Added position writing in event
 *
 * Revision 1.3  2005/06/27 15:44:29  razeto
 * Fixed a typo
 *
 * Revision 1.2  2005/06/22 13:27:21  razeto
 * Fixed a typo
 *
 * Revision 1.1  2005/06/21 12:06:56  razeto
 * Added to repository
 *
 *
 */
