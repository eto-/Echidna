/* BOREXINO Reconstruction program
 *
 *
 * Event's shape variable for Particle Identification and Background rejection
 */
#include <algorithm>
#include "bx_pid_shape.hh"
#include "bx_echidna_event.hh"
#include "db_channel.hh"
#include "db_calib.hh"
#include "bx_dbi.hh"
#include "constants.hh"
#include "barn_interface.hh"
#include "db_profile.hh"
#include "constants.hh"
#include <assert.h>
#include <cmath> 
#include <complex> 

Double_t plane (Double_t* x, Double_t* p){
  return p[0];
}

const float &pi = constants::number::pi;

// ctor
bx_pid_shape::bx_pid_shape() : bx_base_module("bx_pid_shape", bx_base_module::main_loop) {
  require_event_stage (	bx_detector::laben, bx_base_event::reconstructed); 
 //require_trigger_type (bx_trigger_event::neutrino);

} 
// begin
void bx_pid_shape::begin () {
  i4_min_nhits = get_parameter ("min_nhits").get_int ();
  i4_nbins = get_parameter ("nbins").get_int ();
  event_debug = get_parameter ("event_debug").get_int ();


  n_crate = constants::laben::ncrates; //14
  n_feb = n_crate *  constants::laben::frontend::board_per_rack; // 14 * 14 = 196
  n_laben = n_crate * constants::laben::board_per_rack; // 14 * 20 = 280
  nlg = constants::laben::channels;

    // front end vs laben board
  n_good_lg = 0;
  n_on_crate = 0;
  n_on_feb = 0;
  n_on_laben = 0;

   //crate weight
  percentage_of_nhits_crate = new TH2F ("percentage_of_nhits_crate","percentage_of_nhits_crate", n_crate, 1., 1. + n_crate, 120, 0., 60.);
  percentage_of_nhits_feb = new TH2F ("percentage_of_nhits_feb","percentage_of_nhits_feb", n_feb, 1., 1. + n_feb, 40, 0., 20.);
  percentage_of_nhits_laben = new TH2F ("percentage_of_nhits_laben","percentage_of_nhits_laben", n_laben, 1., 1. + n_laben, 40, 0., 20.);
  percentage_of_nhits_lg = new TH2F ("percentage_of_nhits_lg","percentage_of_nhits_lg", nlg, 1., 1. + nlg, 60, 0., 6.);

  nhitted_crates_vs_nhits = new TH2F("nhitted_crates_vs_nhits","nhitted_crates_vs_nhits",300, 1.,600., n_crate + 1, 0., n_crate + 2.);
  nhitted_feb_vs_nhits = new TH2F("nhitted_feb_vs_nhits","nhitted_feb_vs_nhits",500, 1.,1000., n_feb + 1, 0., n_feb + 2.);
  nhitted_laben_vs_nhits = new TH2F("nhitted_laben_vs_nhits","nhitted_laben_vs_nhits",600, 1.,1200., n_laben + 1, 0., n_laben + 2.);
  nhitted_lg_vs_nhits = new TH2F("nhitted_lg_vs_nhits","nhitted_lg_vs_nhits",1000, 1.,2000., nlg + 1, 0., nlg + 2.);

  h_crate_weight = new TH1F ("h_crate_weight", "weight crate" , n_crate, 1., 1. + n_crate);
  h_feb_weight =  new TH1F ("h_feb_weight", "weight feb", n_feb , 1., 1. + n_feb);
  h_laben_weight = new TH1F ("h_laben_weight", "weight laben" , n_laben, 1.,  1. + n_laben);
  h_lg_weight = new TH1F ("h_lg_weight", "weight lg" , nlg, 1.,  1. + nlg);

  //for all events
  crate_all = new TH1F ("crate_all", "not weighted" , n_crate, 1., 1. + n_crate);
  feb_all =  new TH1F ("feb_all", "not weighted", n_feb , 1., 1. + n_feb);
  laben_all = new TH1F ("laben_all", "not weighted" , n_laben, 1.,  1. + n_laben);

   //for all events weighted
  crate_all_weighted = new TH1F ("crate_all_weighted", "per good lg in crate, weighted" , n_crate, 1., 1. + n_crate);
  feb_all_weighted =  new TH1F ("feb_all_weighted", "per good lg in feb, weighted", n_feb , 1., 1. + n_feb);
  laben_all_weighted = new TH1F ("laben_all_weighted", "per good lg in laben, weighted" , n_laben, 1.,  1. + n_laben);

    //for individual events
  crate_1cl = new TH1F ("crate_1cl", "crate", n_crate, 1., 1. + n_crate);
  feb_1cl =  new TH1F ("feb_1cl","feb", n_feb, 1., 1. + n_feb);
  laben_1cl = new TH1F ("laben_1cl", "laben", n_laben, 1., 1. + n_laben);
  lg_1cl = new TH1F ("lg_1cl", "lg", nlg, 1., 1. + nlg);
  
    //found noisy for single events
  bad_crate = new TH1F ("bad_crates", "noisy",  n_crate, 1., 1. + n_crate);
  bad_feb = new TH1F ("bad_feb", "noisy", n_feb, 1., 1. + n_feb);
  bad_laben = new TH1F ("bad_laben","noisy",  n_laben, 1., 1. + n_laben);
  bad_lg = new TH1F ("bad_lg","noisy",  nlg, 1., 1. + nlg);

    //barn_interface
  barn_interface::get ()->store (barn_interface::file, h_crate_weight, this);
  barn_interface::get ()->store (barn_interface::file, h_feb_weight, this);
  barn_interface::get ()->store (barn_interface::file, h_laben_weight, this);
  barn_interface::get ()->store (barn_interface::file, h_lg_weight, this);

  barn_interface::get ()->store (barn_interface::file, percentage_of_nhits_crate , this);
  barn_interface::get ()->store (barn_interface::file, percentage_of_nhits_feb , this);
  barn_interface::get ()->store (barn_interface::file, percentage_of_nhits_laben , this);
  barn_interface::get ()->store (barn_interface::file, percentage_of_nhits_lg , this);  

  barn_interface::get ()->store (barn_interface::file, nhitted_crates_vs_nhits, this);
  barn_interface::get ()->store (barn_interface::file, nhitted_feb_vs_nhits, this);
  barn_interface::get ()->store (barn_interface::file, nhitted_laben_vs_nhits, this);
  barn_interface::get ()->store (barn_interface::file, nhitted_lg_vs_nhits, this);

  barn_interface::get ()->store (barn_interface::file, crate_all, this);
  barn_interface::get ()->store (barn_interface::file, feb_all, this);
  barn_interface::get ()->store (barn_interface::file, laben_all, this);

  barn_interface::get ()->store (barn_interface::file, crate_all_weighted, this);
  barn_interface::get ()->store (barn_interface::file, feb_all_weighted, this);
  barn_interface::get ()->store (barn_interface::file, laben_all_weighted, this);

  barn_interface::get ()->store (barn_interface::junk, crate_1cl,this );
  barn_interface::get ()->store (barn_interface::junk, feb_1cl, this);
  barn_interface::get ()->store (barn_interface::junk, laben_1cl, this);
  barn_interface::get ()->store (barn_interface::junk, lg_1cl, this);

  barn_interface::get ()->store (barn_interface::file, bad_crate, this);
  barn_interface::get ()->store (barn_interface::file, bad_feb, this);
  barn_interface::get ()->store (barn_interface::file, bad_laben, this);
  barn_interface::get ()->store (barn_interface::file, bad_lg, this);
  
  
    // vectors initialization
  crate_weight.resize (n_crate, 0);
  feb_weight.resize (n_feb, 0);
  laben_weight.resize (n_laben, 0);
  lg_weight.resize (nlg, 0);
	
  std::fill_n (crate_weight.begin (), n_crate, 0);
  std::fill_n (feb_weight.begin (), n_feb, 0);
  std::fill_n (laben_weight.begin (), n_laben, 0);
  std::fill_n (lg_weight.begin (), nlg, 0);

  // map of disabled channels from bx_detector
  p_disabled_lg = new bool[constants::laben::channels + 1];
  std::fill_n (p_disabled_lg, constants::laben::channels + 1, false); // Init vector at 0

  //vector of disabled lg
  const std::vector<int32_t>& v = detector_interface::get ()->get_disabled_channels ();
  //do a map from v
  for (int32_t i = 0; i < (int32_t) v.size (); i++) if (v[i] <= constants::laben::channels) {
    int32_t disabled_lg = v[i];
    p_disabled_lg[disabled_lg] = true;
  }

    // count good ordinary lg for crates, feb and laben boards
  for (int32_t ilg = 1;  ilg < (1 + nlg); ilg++) {
    if((bx_dbi::get ()->get_channel (ilg).is_ordinary ()) ) {
      if(!p_disabled_lg[ilg]) {
	n_good_lg ++;
	
	int32_t crate = constants::crate (ilg);
	int32_t fe_board= constants::laben::frontend::board_in_rack (ilg);
	int32_t feb = fe_board + (crate - 1) * constants::laben::frontend::board_per_rack;
	int32_t lab = (ilg - 1) / constants::laben::channels_per_board + 1;
	
	crate_weight[crate - 1] += 1;
	feb_weight[feb - 1] += 1;
	laben_weight[lab - 1] += 1;
	lg_weight[ilg - 1] += 1;       
      }
    }
  }
    // some histograms for test and monitoring purpose
    // phi- theta histogram 
    // the number of bins MUST be an even number
  theta_vs_phi_PMT = new TH2F ("theta_vs_phi_PMT", "phi versus cos(theta) working PMT" , i4_nbins, 0., 2*pi, i4_nbins, -1., 1.);
  theta_vs_phi_PMT_ev_alive = new TH2F ("theta_vs_phi_PMT_ev_alive", "phi versus cos(theta) working PMT" , i4_nbins, 0., 2*pi, i4_nbins, -1., 1.);
  theta_vs_phi_PMT_ev_tot = new TH2F ("theta_vs_phi_PMT_ev_tot", "phi versus cos(theta) all PMT" , i4_nbins, 0., 2*pi, i4_nbins, -1., 1.);
  theta_vs_phi_all = new TH2F ("theta_vs_phi_all", "phi versus cos(theta) all clusters" , i4_nbins, 0., 2*pi, i4_nbins, -1., 1.);
  theta_vs_phi = new TH2F ("theta_vs_phi", "phi versus cos(theta)" , i4_nbins, 0., 2.*pi, i4_nbins, -1., 1.);

  virt_pmt_rel_var = new TH1F ("virt_pmt_rel_var", "realtive variance of virtual PMTs" , 1000, 0., 20.);

  ns_asym = new TH1F ("ns_asym", "north hits in %", 100 , 0., 100.);
  plane_cos_chi2 = new TH2F ("plane_cos_chi2", "Ctheta_vs_phi plane vs chi2", 100 , 0.5, 1.,1000,0,10);
  h_plane_chi2 = new TH1F ("h_plane_chi2", "horisontal plane fit chi2", 1000,0,10);

  sphere_chi2 = new TH1F ("sphere_chi2", "chi square plot", 500 , 0., 50.);
  sphere_lkl = new TH1F ("sphere_lkl", "ln(likelihood) plot", 1000 ,0., 100.);

      // spherical harmonics coefficients
  aYp0 = new TH1F ("aYp0", "spherical harmonics: 0th order coeff. ", 5000 , 0., 5.);
  aYp1 = new TH1F ("aYp1", "spherical harmonics: 1st order coeff. ", 5000 , 0., 5.);
  aYp2 = new TH1F ("aYp2", "spherical harmonics: 2nd order coeff. ", 5000 , 0., 5.);
  aYp3 = new TH1F ("aYp3", "spherical harmonics: 3rd order coeff. ", 5000 , 0., 5.);


  //fill weight histos
  for(int32_t crate = 1; crate <= n_crate; crate ++) h_crate_weight->SetBinContent (crate, crate_weight[crate - 1]);
  for(int32_t feb = 1; feb <= n_feb; feb ++) h_feb_weight->SetBinContent (feb, feb_weight[feb - 1]);
  for(int32_t lbnb = 1; lbnb <= n_laben; lbnb ++) h_laben_weight->SetBinContent (lbnb, laben_weight[lbnb - 1]);
  for(int32_t lg = 1; lg <= nlg; lg ++) {
    h_lg_weight->SetBinContent (lg, lg_weight[lg - 1]);
    if(lg_weight[lg - 1])  {
      const db_channel_laben& info_channel = dynamic_cast<const db_channel_laben&>(bx_dbi::get()->get_channel(lg));
      const double theta = info_channel.pmt_theta ();
      const double phi = info_channel.pmt_phi ();

      const double theta_rad = (90. - theta) * pi/180.;
      const double phi_rad = phi * pi/180.;

      theta_vs_phi_PMT->Fill (phi_rad,::cosf(theta_rad));
    }
 }

  //count number of on crate/feb/laben
  for(int32_t t = 1; t < (n_crate + 1);t++){
    if(crate_weight[t-1] > 0) n_on_crate++;
  }
  
  for(int32_t t = 1; t < (n_feb + 1);t++){
    if(feb_weight[t-1] > 0) n_on_feb++;
  }
  
  for(int32_t t = 1; t < (n_laben + 1);t++){
    if(laben_weight[t-1] > 0) n_on_laben++;
  }
  
  get_message(bx_message::log) << "n_good_lg " << n_good_lg << dispatch;
  get_message(bx_message::log) << "n_on_crate " << n_on_crate << dispatch;
  get_message(bx_message::log) << "n_on_feb " << n_on_feb << dispatch;
  get_message(bx_message::log) << "n_on_laben " << n_on_laben << dispatch;


    // storage of  histograms in the root barn 
  barn_interface::get ()->store (barn_interface::file, theta_vs_phi_PMT, this);
  barn_interface::get ()->store (barn_interface::junk, theta_vs_phi_PMT_ev_tot, this);
  barn_interface::get ()->store (barn_interface::junk, theta_vs_phi_PMT_ev_alive, this);
  barn_interface::get ()->store (barn_interface::file, theta_vs_phi_all, this);
  barn_interface::get ()->store (barn_interface::junk, theta_vs_phi, this);

  barn_interface::get ()->store (barn_interface::file, virt_pmt_rel_var, this);

  barn_interface::get ()->store (barn_interface::file, ns_asym, this);
  barn_interface::get ()->store (barn_interface::file, plane_cos_chi2, this);
  barn_interface::get ()->store (barn_interface::file, h_plane_chi2, this);

  barn_interface::get ()->store (barn_interface::file, sphere_chi2, this);
  barn_interface::get ()->store (barn_interface::file, sphere_lkl, this);
  barn_interface::get ()->store (barn_interface::file, aYp0, this);
  barn_interface::get ()->store (barn_interface::file, aYp1, this);
  barn_interface::get ()->store (barn_interface::file, aYp2, this);
  barn_interface::get ()->store (barn_interface::file, aYp3, this);

}


//DO IT

bx_echidna_event* bx_pid_shape::doit (bx_echidna_event *ev) {


  int32_t evnum = ev->get_event_number ();

  int32_t n_crate = constants::laben::ncrates; //14
  int32_t n_feb = n_crate *  constants::laben::frontend::board_per_rack; // 14 * 12 = 196
  int32_t n_laben = n_crate * constants::laben::board_per_rack; // 14 * 20 = 280
    
  int32_t index_hot_electronics = 0;
  int32_t index_low_electronics = 0;
  
  for (int32_t cluster = 0; cluster < ev->get_laben ().get_nclusters (); cluster++) {

    //reset histos 
    feb_1cl->Reset ();
    laben_1cl->Reset ();
    crate_1cl->Reset ();
    lg_1cl->Reset ();

    float north_counts = 0;
    float total_counts = 0;
    theta_vs_phi->Reset ();
    theta_vs_phi_PMT_ev_tot->Reset ();
    theta_vs_phi_PMT_ev_alive->Reset ();

        //number of clustered hits
    int32_t nhits =  ev->get_laben ().get_cluster (cluster).get_clustered_nhits (); 
    //double ev_charge =  ev->get_laben ().get_cluster (cluster).get_charge (); 
    
       //get shaped_cluster
    bx_laben_shaped_cluster& cluster_ref = ev->get_laben().get_shaped_cluster (cluster);


    const double mi_x = ev->get_laben ().get_rec_cluster (cluster).get_position ().get_x();
    const double mi_y = ev->get_laben ().get_rec_cluster (cluster).get_position ().get_y();
    const double mi_z = ev->get_laben ().get_rec_cluster (cluster).get_position ().get_z();
      
    const TVector3 mi_pos (mi_x, mi_y, mi_z);

    //calculate live PMT weight transhormed in the eve. position
    for(int32_t lg = 1; lg <= nlg; lg ++) {
      if(!bx_dbi::get ()->get_channel (lg).is_empty ()) {
	const db_channel_laben& info_channel = dynamic_cast<const db_channel_laben&>(bx_dbi::get()->get_channel(lg));
	const TVector3 PMT_pos (info_channel.pmt_x (), info_channel.pmt_y (), info_channel.pmt_z ());
        
	const TVector3 distance = PMT_pos - mi_pos;
  
	const double cos_theta = distance.CosTheta ();
	const double phi_rad = distance.Phi () + pi;
	theta_vs_phi_PMT_ev_tot->Fill (phi_rad,cos_theta);
	if(lg_weight[lg - 1]) theta_vs_phi_PMT_ev_alive->Fill (phi_rad,cos_theta);
      }
    }
  

  
        // spherical harmonics coefficients
    
    const std::complex <float> im_i(0., 1.); // imaginary unit
    const std::complex <float> re_1 = 1.; // real unit
   
    std::complex <float> aY00 = 0.;
    std::complex <float> aY10 = 0;
    std::complex <float> aY1p1 = 0;
    std::complex <float> aY1n1 = 0.;
    std::complex <float> aY20 = 0.;
    std::complex <float> aY2p1 = 0.;
    std::complex <float> aY2n1 = 0.;
    std::complex <float> aY2p2 = 0.;
    std::complex <float> aY2n2 = 0.;
    std::complex <float> aY30 = 0.;
    std::complex <float> aY3n1 = 0.;
    std::complex <float> aY3p1 = 0.;
    std::complex <float> aY3n2 = 0.;
    std::complex <float> aY3p2 = 0.;
    std::complex <float> aY3n3 = 0.;
    std::complex <float> aY3p3 = 0.;

    if (nhits < i4_min_nhits) {
      cluster_ref.f4_ns_asymmetry = -1.;
      cluster_ref.f4_sphere_chi2 = -1.;
      cluster_ref.f4_sphere_lkl = -10.;
      cluster_ref.f4_sphere_rel_var = -1;

      cluster_ref.f4_plane_cos = 2.;
      cluster_ref.f4_plane_chi2 = -1.;
      cluster_ref.f4_h_plane_chi2 = -1.;


      cluster_ref.v_sh_power[0] = -10.; 
      cluster_ref.v_sh_power[1] = -10.;
      cluster_ref.v_sh_power[2] = -10.;
      cluster_ref.v_sh_power[3] = -10.;	

      continue;
    }

       	
    // Loop on every clustered hit
    for (int32_t hit = 0; hit < nhits; hit++) {

        // Get db_channel pointer for hit
      const db_channel_laben* ch_info =  ev->get_laben ().get_cluster (cluster). get_clustered_hit (hit).get_decoded_hit ().get_db_channel ();
      //        int32_t npe =  ev->get_laben ().get_cluster (cluster). get_clustered_hit (hit).get_decoded_hit ().get_charge_npe ();
	
	//  *******************************
	// crates, fe and laben boards
	//  *******************************
	
	  // get crate, fe and laben boards corresponding to each clustered hit 
      int32_t  lg = ch_info->get_lg ();

      int32_t crate = constants::crate (lg);
      int32_t fe_board= constants::laben::frontend::board_in_rack (lg);
      int32_t feb = fe_board + (crate - 1) * constants::laben::frontend::board_per_rack;
      int32_t lab = (lg - 1) / constants::laben::channels_per_board + 1;
      
      // fill histos
      //lg
      lg_1cl->Fill (lg);
      
      // crates
      if(crate_weight[crate - 1]){
	crate_all->Fill (crate);
	crate_all_weighted->Fill (crate, 1 / crate_weight[crate - 1]);
	crate_1cl->Fill (crate);
	//crate_1cl->Fill (crate, 160 / crate_weight[crate - 1]);
      }
      else if(crate_weight[crate - 1] ==  0)
	get_message(bx_message::warn) << "Crate " << crate << " has 0 good_lg, but lg " << lg << " gives clustered hit "<<dispatch;	  
      
      //feb
      if(feb_weight[feb - 1]){
	feb_all->Fill (feb);
	feb_all_weighted->Fill (feb, 1 / feb_weight[feb - 1]);
	feb_1cl->Fill (feb);
	//feb_1cl->Fill (feb, 12 / feb_weight[feb - 1]);
      }
      else if(feb_weight[feb - 1] == 0)
	get_message(bx_message::warn) << "Feb " << feb << " has 0 good_lg, but lg " << lg << " gives clustered hit "<<dispatch;
      
      //laben
      if(laben_weight[lab - 1]){
	laben_all->Fill (lab);
	laben_all_weighted->Fill (lab, 1 / laben_weight[lab - 1]);
	laben_1cl->Fill (lab);
	//laben_1cl->Fill (lab, 8 / laben_weight[lab - 1]);
      }
      else if(laben_weight[lab - 1] == 0)
	get_message(bx_message::warn) << "Laben board " << lab << " has 0 good_lg, but lg " << lg << " gives clustered hit "<<dispatch;
      
	
      //  *****************
      // sphericity etc..
      //  *****************
      
      // get pmt position of this hit
      const TVector3 hit_pmt_pos (ch_info->pmt_x (), ch_info->pmt_y (), ch_info->pmt_z ());
      
      const TVector3 distance2 = hit_pmt_pos - mi_pos;      
      
      const double cos_theta = distance2.CosTheta ();
      const double phi_rad =  distance2.Phi () + pi;
      
      double hit_charge =  ev->get_laben ().get_cluster (cluster). get_clustered_hit (hit).get_decoded_hit ().get_charge_pe ();
      if(ev->get_trigger ().is_neutron()  && hit_charge<0.1) hit_charge=1;//noavg correction for charge in non-neutrino trigger 
      if(hit_charge) theta_vs_phi->Fill (phi_rad,cos_theta,hit_charge);
						
	
      total_counts ++;

      if(cos_theta > 0.) north_counts ++;


      // check for event sfericity with spherical harmonics	
      const float ct = cos_theta;
      const float st = ::sqrtf (1 - ct * ct);
      const  float norm = ::sqrtf (pi) *2 / nhits;
      
      const std::complex <float> i_phi (0,phi_rad); 
      const std::complex <float> i2_phi (0,2*phi_rad); 
      const std::complex <float> i3_phi (0,3*phi_rad); 
      
      aY00  += norm * ::sqrtf (1. / 4. / pi);
		
      aY10  += norm * ::sqrtf (3. / 4. / pi) * ct;;
      aY1p1 -= norm * ::sqrtf (3. / 8. / pi) * st * exp(i_phi);
      aY1n1 += norm * ::sqrtf (3. / 8. / pi) * st * exp(-i_phi);
      
      aY20  += norm * ::sqrtf (5. / 16 / pi) * (3 * ct * ct - 1);
      aY2p1 -= norm * ::sqrtf (15. / 8 / pi) * ct * st * exp(i_phi);
      aY2n1 += norm * ::sqrtf (15. / 8 / pi) * ct * st * exp(-i_phi);
      aY2p2 += norm * ::sqrtf (15. / 32 / pi) * st * st * exp(i2_phi);
      aY2n2 += norm * ::sqrtf (15. / 32 / pi) * st * st * exp(-i2_phi);
      
      aY30  += norm * ::sqrtf (7. / 16 / pi) * ct * (5 * ct * ct - 3);
      aY3p1 -= norm * ::sqrtf (21. / 64 / pi) * st * (5 * ct * ct - 1) * exp(i_phi);
      aY3n1 += norm * ::sqrtf (21. / 64 / pi) * st * (5 * ct * ct - 1) * exp(-i_phi);
      aY3p2 += norm * ::sqrtf (105. / 32 / pi) * ct * st * st * exp(i2_phi);
      aY3n2 += norm * ::sqrtf (105. / 32 / pi) * ct * st * st * exp(-i2_phi);
      aY3p3 -= norm * ::sqrtf (35. / 64 / pi) * st * st * st * exp(i3_phi);
      aY3n3 += norm * ::sqrtf (35. / 64 / pi) * st * st * st * exp(-i3_phi);
      
    } // end of loop on clustered hits
    


    //  *************************************************
    // OUT of hits loop: crates, fe and laben boards
    //  *************************************************

    // double hits_per_lg_expected = (double) nhits / (double) n_good_lg; 

    double hot_hit_fraction_crate = 0.40;
    double hot_hit_fraction_feb = 0.12;
    double hot_hit_fraction_laben = 0.10;
    double hot_hit_fraction_lg = 0.045;


    /*
    //mean and rms
    double mean_hits_in_crate = (double) nhits / (double) n_on_crate;
    //    double mean_hits_in_crate =  hits_per_lg_expected * 160;
    if(mean_hits_in_crate < 1.)  mean_hits_in_crate = 1; 
    double rms_hits_in_crate = sqrt(mean_hits_in_crate);

    double mean_hits_in_feb = (double) nhits / (double) n_on_feb;
    //double mean_hits_in_feb =  hits_per_lg_expected * 12;
    if(mean_hits_in_feb < 1.)  mean_hits_in_feb = 1; 
    double rms_hits_in_feb = sqrt(mean_hits_in_feb);

    double mean_hits_in_laben = (double) nhits / (double) n_on_feb;
    //double mean_hits_in_laben = hits_per_lg_expected * 8;
    if(mean_hits_in_laben < 1.)  mean_hits_in_laben = 1; 
    double rms_hits_in_laben = sqrt(mean_hits_in_laben);
    */

    
   // crate distribution for this event
    int32_t nhitted_crates = 0;
    for(int32_t cr = 1; cr < (1 + n_crate); cr++){
      double bin_content =  crate_1cl->GetBinContent (cr);
      if (bin_content > 0.) nhitted_crates ++;
      percentage_of_nhits_crate->Fill(cr, bin_content/nhits*100.);
      if(bin_content >= hot_hit_fraction_crate * nhits) {
	bad_crate->Fill (cr);
	index_hot_electronics = 1;
	if(event_debug) get_message(bx_message::warn) << " noisy crate " << cr << " Event number " << evnum << dispatch;
      }
    }
    nhitted_crates_vs_nhits->Fill(nhits, nhitted_crates);
    if(nhits < 20 && nhitted_crates < 4) index_low_electronics = 1;
    if(nhits >=20 && nhits < 100 && nhitted_crates < 0.1 * nhits) index_low_electronics = 1;
    if(nhits >= 100 && nhits < 350 &&  nhitted_crates < (0.016 * nhits + 8.4) ) index_low_electronics = 1;
    if(nhits >= 350 &&  nhitted_crates < 14 ) index_low_electronics = 1;
    
    
        // feb distribution for this event
    int32_t nhitted_feb = 0;
    for(int32_t fe = 1; fe < (1 + n_feb); fe++){
      double bin_content =  feb_1cl->GetBinContent (fe);
      if (bin_content > 0.) nhitted_feb ++;
      percentage_of_nhits_feb->Fill(fe,bin_content/nhits*100.);
      if(bin_content >= hot_hit_fraction_feb * nhits) {
	index_hot_electronics = 1;
	bad_feb->Fill (fe);
	if(event_debug) get_message(bx_message::warn) << " noisy feb " << fe << " Event number " << evnum << dispatch;
      }
    }
    nhitted_feb_vs_nhits->Fill(nhits, nhitted_feb);
    if(nhits < 40 && nhitted_feb < 10) index_low_electronics = 1;
    if(nhits >= 40  && nhits < 300 &&  nhitted_feb < n_on_feb/196 * (0.346 * nhits - 3.8)) index_low_electronics = 1;
    if(nhits >= 300 && nhits < 800 &&  nhitted_feb < n_on_feb/196 * (0.12 * nhits + 64) ) index_low_electronics = 1;
    if(nhits >= 800 && nhitted_feb < n_on_feb/196 * 160) index_low_electronics = 1;
   
    
      // laben distribution for this event
    int32_t nhitted_laben = 0;
    for(int32_t lab = 1; lab < (1 + n_laben); lab++){
      double bin_content =  laben_1cl->GetBinContent (lab);
      if (bin_content > 0.) nhitted_laben ++;
      percentage_of_nhits_laben->Fill(lab, bin_content/nhits*100.);
      if(bin_content >= hot_hit_fraction_laben * nhits) {
	bad_laben->Fill (lab);
	index_hot_electronics = 1;
	if(event_debug) get_message(bx_message::warn) << " noisy laben " << lab << " Event number " << evnum << dispatch;
      }
    }
    nhitted_laben_vs_nhits->Fill(nhits, nhitted_laben);
    if(nhits < 40   && nhitted_laben < 10) index_low_electronics = 1;
    if(nhits >= 40 && nhits < 400 &&  nhitted_laben < n_on_laben/278 * (0.417 * nhits - 6.68) ) index_low_electronics = 1;
    if(nhits >= 400 && nhits < 900 &&  nhitted_laben < n_on_laben/278 * (0.18 * nhits + 88.) ) index_low_electronics = 1;
    if(nhits >= 900 && nhitted_laben < n_on_laben/278 * 250) index_low_electronics = 1;
    
   
    
      // lg distribution for this event
    int32_t nhitted_lg = 0;
    for(int32_t l = 1; l < (1 + constants::laben::channels); l++){
      double bin_content =  lg_1cl->GetBinContent (l);
      if (bin_content > 0.) nhitted_lg ++;
      percentage_of_nhits_lg->Fill(l, bin_content/nhits*100.);
      if(bin_content >= hot_hit_fraction_lg * nhits) {
	bad_lg->Fill (l);
	index_hot_electronics = 1;
	if(event_debug) get_message(bx_message::warn) << " noisy lg " << l << " Event number " << evnum << dispatch;
      }
    }
    nhitted_lg_vs_nhits->Fill(nhits, nhitted_lg);
    if(nhits < 50 && nhitted_lg < 0.5 * nhits) index_low_electronics = 1;
    if(nhits >= 50 && nhits < 1800 && nhitted_lg <  n_good_lg/1911. * (0.9 * nhits - 20.)) index_low_electronics = 1;
    if(nhits >= 1800 &&  nhitted_lg < n_good_lg/1911. * 1600.) index_low_electronics = 1;
     


      //store strange event in the vector
    if(index_low_electronics || index_hot_electronics ) ev_el.push_back (evnum);
      //for hot events setting bit 0
    if(index_hot_electronics) cluster_ref.i1_quality_flags |= 1 << 0;
      // for low occupancy events setting bit 1
    if(index_low_electronics) cluster_ref.i1_quality_flags |= 1 << 1;

    
    //  *******************************
    // OUT of hits loop: sphericity etc
    //  *******************************

    // 1) check for north/south symmetry
    double ns_asymm = north_counts / total_counts * 100.;
    ns_asym->Fill(ns_asymm);



    // 2) check for event sfericity: angle of simple plane z = ax + by + d with xy plane z= 0  (fitted in theta - phi distribution )


    //normalize theta_vs_phi with alived PMTs in each virtual PMT (taken from theta_vs_phi_PMT) 
    double sum_xi = 0.; //sum of theta_vs_phi histo
    double sum_xi2 = 0.;

    for(int32_t i = 1; i <= i4_nbins; i++) {
      for(int32_t j = 1; j <= i4_nbins;j++) {
	double norm = 1;
        const double alive = theta_vs_phi_PMT_ev_alive->GetBinContent(i, j);
	const double all = theta_vs_phi_PMT_ev_tot->GetBinContent(i, j);
        if(all) norm = alive/ all; //per n PMT in 1 virtual PMT
	if(norm) {
	  theta_vs_phi->SetBinContent(i,j,theta_vs_phi->GetBinContent(i,j)/norm);
	  //summing in _all histo
	  theta_vs_phi_all->SetBinContent(i,j, theta_vs_phi_all->GetBinContent(i,j) + theta_vs_phi->GetBinContent(i,j));	  

	  double bin_content = theta_vs_phi->GetBinContent(i,j);
	  sum_xi  +=  bin_content;
	  sum_xi2 +=  bin_content*bin_content;
	}
      }
    }      
    
    const double f8_virt_pmt_mean = sum_xi /(i4_nbins * i4_nbins);
    const double f8_virt_pmt_var = std::sqrt(sum_xi2 /(i4_nbins * i4_nbins) - f8_virt_pmt_mean*f8_virt_pmt_mean);
    const double f8_virt_var = (f8_virt_pmt_var/f8_virt_pmt_mean);
    virt_pmt_rel_var->Fill(f8_virt_var);
  
    const int32_t n_bins = i4_nbins * i4_nbins;
    const double mean = sum_xi /(i4_nbins * i4_nbins);

    //set fit function: plane with angle
    TF2 f2("f2", "[0]*x + [1]*y + [2]", 0, 2*pi, -1, 1); 
    f2.SetParameters(0.0, 0.0, mean);
    theta_vs_phi->Fit("f2","QRN");
    double par0 = f2.GetParameter (0);
    double par1 = f2.GetParameter (1);
    double f8_plane_chi2 = f2.GetChisquare ()/f2.GetNDF ();
    //calculate angle
    double cos_a = fabs ( 1 / (::sqrt(par0*par0 + par1*par1 + 1)));
    plane_cos_chi2->Fill (cos_a, f8_plane_chi2);


    TF2 f3("f3", plane, 0, 2*pi, -1, 1, 1); 
    f3.SetParameter(0, mean);
    theta_vs_phi->Fit("f3","RQN");
    double f8_h_plane_chi2 = f3.GetChisquare ()/f3.GetNDF ();
    h_plane_chi2->Fill (f8_h_plane_chi2);


    // compute lkl and chi2
    double sum_log_fact = 0.; // sum of the logaritms of ni!
    double a = 0., b = 0.;

    for(int32_t i = 1; i <= i4_nbins; i++) {
      for(int32_t j = 1; j <= i4_nbins; j++) {
	double ni =  theta_vs_phi->GetBinContent (i, j);//= bin content
	a += ni;
	b += ni*ni;
	
	for(int32_t k = 1; k <= (int32_t) (ni + 0.5); k++) //+0.5 to have (int32_t) 3.9 = 4!!
	  sum_log_fact += ::log((double) k); //sum = log (ni!)
	
      }
    }
    
      //log of Poissonian likelihood
    double ln_lkl = -(sum_xi * std::log(mean) - mean * n_bins - sum_log_fact)/sum_xi;//normalized per sum_xi and set as positive
    //    double ln_lkl = (nhits * std::log(mean) - mean * n_bins - sum);

      //Pearsons chi2 test
    double chi2_value = mean + b/sum_xi - 2*a/n_bins; 
    sphere_chi2->Fill (chi2_value);
    sphere_lkl->Fill (ln_lkl);

    // compute spherical harmonics powers
    
    const float p0 = norm(aY00);
    const float p1 = norm (aY10) + norm (aY1p1) + norm(aY1n1);
    const float p2 = norm (aY20) + norm (aY2p1) + norm(aY2n1) + norm (aY2p2) + norm(aY2n2);
    const float p3 = norm (aY30) + norm (aY3p1) + norm(aY3n1) + norm (aY3p2) + norm(aY3n2) + norm (aY3p3) + norm(aY3n3);
    

    // Filling histograms in the root barn 
    aYp0->Fill (p0);
    aYp1->Fill (p1);
    aYp2->Fill (p2);
    aYp3->Fill (p3);
    
    //Filling variables in the Echidna bx_cluster_refent:
    cluster_ref.f4_ns_asymmetry = ns_asymm;
    cluster_ref.f4_sphere_chi2 = chi2_value;
    cluster_ref.f4_sphere_lkl = ln_lkl;
    cluster_ref.f4_sphere_rel_var = f8_virt_var;

    cluster_ref.f4_plane_cos = cos_a;
    cluster_ref.f4_plane_chi2 = f8_plane_chi2;
    cluster_ref.f4_h_plane_chi2 = f8_h_plane_chi2;

    
    cluster_ref.v_sh_power[0] = p0;
    cluster_ref.v_sh_power[1] = p1;
    cluster_ref.v_sh_power[2] = p2;
    cluster_ref.v_sh_power[3] = p3;	
    
  } // end of loop on clusters
  return ev;
}

void bx_pid_shape::end () {
  delete crate_1cl;
  delete feb_1cl;
  delete laben_1cl;
  delete lg_1cl;
  delete theta_vs_phi;
  delete theta_vs_phi_PMT_ev_tot;
  delete theta_vs_phi_PMT_ev_alive;
  
  delete [] p_disabled_lg;

  //evnum of suspicious_events
  int32_t number_of_bad_events = ev_el.size ();
  get_message(bx_message::log) << "Number of electronics bad events is " << number_of_bad_events << dispatch;
  suspicious_events = new TH1F ("suspicious_events", "suspicious_events", number_of_bad_events , 1., number_of_bad_events + 1.);
  for (int32_t i=1;i <= number_of_bad_events; i++) suspicious_events->SetBinContent(i, ev_el[i-1]);
  
  barn_interface::get ()->store (barn_interface::file, suspicious_events, this);
  
}

/* 
 * $Log: bx_pid_shape.cc,v $
 * Revision 1.34  2015/07/28 14:04:03  ilia.drachnev
 * set production in all triggers; noavg corrected charge in non-neutrino trigger
 *
 * Revision 1.33  2011/03/02 10:45:25  ludhova
 * required reconstructed
 *
 * Revision 1.32  2011-03-02 08:40:00  ludhova
 * rec cluster position instead of milano position
 *
 * Revision 1.31  2009-10-26 17:25:27  ludhova
 * removed debug messages
 *
 * Revision 1.30  2009-10-23 07:53:01  ludhova
 * removed log message
 *
 * Revision 1.29  2009-10-08 17:35:25  ludhova
 * new variables
 *
 * Revision 1.28.2.1  2009-09-28 11:56:06  ludhova
 * distribution in cosTheta-Phi, charge and position corrections
 *
 * Revision 1.28  2009-01-24 12:03:43  razeto
 * Amin commit: call fit with N to avoid fit drawing
 *
 * Revision 1.27  2008-10-01 16:17:52  ddangelo
 * removed is_pointlike and related stuff (bx_filter will perform this task now)
 *
 * Revision 1.26  2008-08-21 16:32:38  ddangelo
 * patched to consider only laben channels during disabling
 *
 * Revision 1.25  2008-04-23 13:27:33  ludhova
 * redefining is_pointlike based on the new IDF
 *
 * Revision 1.24  2008-02-26 18:31:56  ddangelo
 * added filling of is_pointlike variable based on meantime and peaktime values
 * added paramaters with limits and their readout from config
 * Events without any peak are considered pointlike.
 *
 * Commit done while maintainer was out of her mind.
 *
 * Revision 1.23  2007-10-30 13:55:04  ludhova
 * few changes in spherical harmonics and limits on hot and low el.
 *
 * Revision 1.22  2007-10-29 09:31:35  ludhova
 * few changes in spherical harmonics
 *
 * Revision 1.21  2007-10-26 12:17:56  ludhova
 * setting  cluster_ref.i1_quality_flags
 *
 * Revision 1.20  2007-10-26 10:46:01  ludhova
 * bug removal
 *
 * Revision 1.19  2007-10-26 09:55:29  ludhova
 * new binning for spherical harmonics histos
 *
 * Revision 1.18  2007-10-26 09:42:14  ludhova
 * setting of some rec_clus variables for too low E clusters to default values  +
 *
 * Revision 1.17  2007-10-25 14:43:22  ludhova
 * filling new variable sphere_cos
 *
 * Revision 1.16  2007-10-25 10:38:49  ludhova
 * some cosmetics
 *
 * Revision 1.15  2007-10-24 13:59:18  ludhova
 * some deletes
 *
 * Revision 1.14  2007-10-24 13:53:32  ludhova
 * new search for events with strange electronics distribution
 *
 * Revision 1.13.4.1  2007-10-20 15:40:50  ludhova
 * deletes of histos sent to barn removed
 *
 * Revision 1.13  2007-04-04 12:20:10  ddangelo
 * correct resource dealllocation for TF and TH objects
 * cycle on clusters reintroduced
 * writing to new event structure introduced
 * changes agreed with maintainer
 *
 * Revision 1.12  2007-03-15 19:52:48  ddangelo
 * writing to (old) pid event commented out.
 *
 * Revision 1.11  2006/11/02 13:57:42  ludhova
 * nhits replaced by npe for symmetry calculations
 *
 * Revision 1.10  2006-08-30 16:01:43  ludhova
 * removed verbose mode of fitting
 *
 * Revision 1.9  2006/08/30 15:51:05  ludhova
 * identificattion of bad crates and boards change *
 * Revision 1.7  2006/05/08 17:31:31  razeto
 * Added db_channel patch for fadc (sent to the mailing list)
 *
 * Revision 1.6  2006/01/09 16:21:03  razeto
 * Updated to the new rootrn target
 *
 * Revision 1.5  2005/09/23 14:14:25  zavatare
 * Several fixes
 *
 * Revision 1.4  2005/09/21 14:00:43  razeto
 * Added binning in spheric harmonic calculation
 *
 * Revision 1.3  2005/09/21 12:49:33  razeto
 * Fixed normalization coefficient for spherical harmonic
 
 * Revision 1.2  2005/09/20 17:19:46  razeto
 * Updated to a just working code
 */

