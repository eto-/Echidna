/* BOREXINO Reconstruction program
 *
 * Author: Michael Wurm <michael.wurm@ph.tum.de>
 * Maintainer: Davide D'Angelo <Davide.Dangelo@lngs.infn.it>
 *
 * $Id: bx_laben_tof_tracker.cc,v 1.19 2014/12/11 21:27:12 wurm Exp $
 *
 * Implementation of bx_laben_tof_tracker
 * If you cut and paste from this file, 
 * remember to remove comments
 * 
 * A few lines of sample code are commented out
 * to allow clean compilation.
 * C-style comments are used in this cases, 
 * while C++style are used for explanations
 */

#include "bx_laben_tof_tracker.hh"
#include "bx_echidna_event.hh"
#include "barn_interface.hh"
#include "db_channel.hh"
#include "db_run.hh"
#include "db_profile.hh"
#include "bx_dbi.hh"
#include "TH1F.h"
#include "TH2F.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "math.h"

// constructor
bx_laben_tof_tracker::bx_laben_tof_tracker (): bx_base_module("bx_laben_tof_tracker", bx_base_module::main_loop) {
  require_event_stage (bx_detector::laben, bx_base_event::tracked);
  require_trigger_type (bx_trigger_event::neutrino);
}  // event has to be reconstructed by energy tracker first, as variables of the entry point are re-used


// module Interface
void bx_laben_tof_tracker::begin () {
	//get_message (bx_message::debug) << "test begin" << dispatch;

  // resource initialization
	i4_laben_events = 0;
	i4_tracked_events = 0;
	i4_low_nhits = 1000; 
	
	ep_time_cut = 50.;
	pi = constants::number::pi;
	ep_max_theta_low  = 0.3;
	ep_max_theta_high = 0.1;
	ep_max_phi_low  = 0.3;
	ep_max_phi_high = 0.2;	
	ep_diff_hits = 15000;
	ep_par_threshold = 7.5;
	ep_theta_pos_threshold = 5.;
	
	xp_plains = 12;
	xp_theta_diff = 0.75;
	xp_nhits_diff = 4500;
	xp_theta_bins = 60;
}


bx_echidna_event* bx_laben_tof_tracker::doit (bx_echidna_event *ev) {

	// Get const references to database visitors with 1 or more keys (run, profile, calib_profile)
	//const db_run& run_info = bx_dbi::get ()->get_run ();
	//const db_profile& profile_info = bx_dbi::get ()->get_profile ();

	//get_message(bx_message::log) << "event " << ev->get_event_number() << dispatch;

	i4_laben_events ++;

	const bx_laben_event& er = ev->get_laben();
	const bx_track_by_points& et = (bx_track_by_points&) er.get_track_energy();
	bx_laben_tracked_event &ew = dynamic_cast<bx_laben_tracked_event&>(ev->get_laben());
	ew.b_is_tracked_tof = false;
	
	
	// only events with at least one cluster
	if (!er.get_nclusters()) {
		//get_message (bx_message::debug) << "no cluster present" << dispatch;
		return ev;
	}

	// exclude low-hit muons from xp reconstruction
	bool lowhits = false;
	if (er.get_decoded_nhits()<i4_low_nhits) {
		//get_message (bx_message::log) << "low hits, modified reconstruction" << dispatch;
		lowhits = true;
		return ev;
	}
	
	// exclude muons w/o entry point  from laben energy tracking
	if (et.get_dx1()==0) {
		//get_message (bx_message::log) << "no entry point available, skip" << dispatch;
		return ev;
	}
	
	// ENTRY POINT -------------------------------------------------
	
	// call old entry point
	coo oep;
	oep.x = et.get_x1();
	oep.y = et.get_y1();
	oep.z = et.get_z1();
	oep.phi = atan2(oep.y, oep.x);
	oep.rad = m_get_radius(oep);
	oep.theta = acos(oep.z/oep.rad);

	//get_message (bx_message::log) << "old ep: x=" << oep.x << ", y=" << oep.y << ", z=" << oep.z << ", phi=" << oep.phi << ", theta=" << oep.theta << ", rad=" << oep.rad << dispatch;
	

	// read clustered hits for new entry point fit
	
	// define maximum number of read hits, block writing of second hits in channel
	std::vector<coo> ep_pmts;
	
	UInt_t lasthit = constants::laben::channels;
	UInt_t  firstclusterhits = er.get_cluster(0).get_npmts();//get_clustered_nhits();
	if (firstclusterhits<lasthit) lasthit = firstclusterhits;
	Float_t mean_rad = 0;
	//get_message (bx_message::log) << lasthit << dispatch;
	
	
	// loop over hits, convert them to ep coordinates, save them to ep_pmts
	for (UInt_t ihit=0; ihit<lasthit; ihit++) {
		
		const bx_laben_clustered_hit& chit = er.get_cluster(0).get_clustered_hit(ihit);
		const db_channel_laben* c = chit.get_decoded_hit().get_db_channel();
		if (c->get_lg() > constants::laben::channels) continue;
		//	Int_t lg = c->pmt_lg();		
	//	if (lg > constants::laben::channels) continue;
	//	if (occupied[lg]) continue;
	//	occupied[lg] = true;
		
		coo pmt;
		pmt.time = chit.get_time();
		/*if ( pmt.time > ep_time_cut ) {
			get_message (bx_message::log) << pmt.time << " ns ---> time cut effective" << dispatch;
			continue;
		}*/
		pmt.x = c->pmt_x();
		pmt.y = c->pmt_y();
		pmt.z = c->pmt_z();
		pmt.rad = m_get_radius(pmt);
		mean_rad += pmt.rad;
		
		coo pmt_trans = m_rotate_y(m_rotate_z(pmt,pi-oep.phi),pi/2-oep.theta);
		pmt_trans.phi = m_limit_angle(atan2(pmt_trans.y, pmt_trans.x)+pi,-pi,pi);
		pmt_trans.theta = m_limit_angle(acos(pmt_trans.z/pmt_trans.rad)-pi/2,-pi/2,pi/2);
		//get_message (bx_message::log) << "pmt_trans.z/rad: " << pmt_trans.z << " " << pmt_trans.rad << dispatch;
		//get_message (bx_message::log) << ihit << " : pmt_trans.phi: " << pmt_trans.phi << ", pmt_trans.theta: " << pmt_trans.theta << dispatch;

		ep_pmts.push_back(pmt_trans);
		
	}
	mean_rad /= ep_pmts.size();
	
	
	// produce graph for phi and theta fits;
	
	Float_t graph_phi[constants::laben::channels], graph_time[constants::laben::channels];
	Int_t graph_n = 0;
		
	for (UInt_t ihit=0; ihit<ep_pmts.size() && graph_n<constants::laben::channels; ihit++) {
		if (fabs(ep_pmts[ihit].theta) > ep_max_theta_low) continue;
		if (er.get_decoded_nhits() > ep_diff_hits && fabs(ep_pmts[ihit].theta) > ep_max_theta_high) continue;
		graph_phi[graph_n] = ep_pmts[ihit].phi;
		graph_time[graph_n] = ep_pmts[ihit].time;
		//get_message (bx_message::log) << graph_n << " : " << ihit << " : " << graph_phi[graph_n] << " : " << graph_time[graph_n] << dispatch;
		graph_n++;
	}
	if (!graph_n) return ev; // AAA: Line added by DD to avoid crash (21/08/10)
	//get_message (bx_message::warn) << "1" << dispatch;
			
	TGraph ep_graph_phi(graph_n, graph_phi, graph_time);
	//TF1 ep_fit_phi("ep_fit_phi","[1]*(-1)*(-1+TMath::Sign(1,x-[0]))*(sqrt(abs(x-[0])/[2]))+[3]*(1+TMath::Sign(1,x-[0]))*(sqrt((x-[0])/[4]))+[5]",-pi,pi);
	TF1 ep_fit_phi("ep_fit_phi","[1]*(-1)*(-1+TMath::Sign(1,x-[0]))*(-1*sin((x-[0])/[2]))+[3]*(1+TMath::Sign(1,x-[0]))*(sin((x-[0])/[4]))+[5]",-pi,pi);
	ep_fit_phi.SetParameters(0,4.,2.,4.,2.,0);	
	ep_fit_phi.SetParLimits(1,0,1e6);
	ep_fit_phi.SetParLimits(0,-0.5,0.5);
	ep_fit_phi.SetParLimits(5,-10,10);
	ep_graph_phi.Fit(&ep_fit_phi,"RQ0","",-1.,1.);
	
	if (isnan(ep_fit_phi.GetParameter(0)) || isnan(ep_fit_phi.GetParError(0))) {
		//get_message (bx_message::warn) << "ep phi fit did not converge" << dispatch;
		return ev;
	}
	//get_message (bx_message::warn) << "2" << dispatch;
				
	Float_t graph_theta[constants::laben::channels];
	graph_n = 0;
	
	for (UInt_t ihit=0; ihit<ep_pmts.size() && graph_n<constants::laben::channels; ihit++) {
		if (fabs(ep_pmts[ihit].phi) > ep_max_phi_low) continue;
		if (er.get_decoded_nhits() > ep_diff_hits && fabs(ep_pmts[ihit].phi) > ep_max_phi_high) continue;
		//get_message (bx_message::warn) << ep_pmts[ihit].theta << dispatch;
			
		graph_theta[graph_n] = ep_pmts[ihit].theta;
		graph_time[graph_n] = ep_pmts[ihit].time;
//		get_message (bx_message::warn) << graph_n << " : " << ihit << " => " << graph_theta[graph_n] << " : " << graph_time[graph_n] << dispatch;
		graph_n++;
	}
	if (!graph_n) return ev; // AAA: Line added by DD to avoid crash (21/08/10)
	//get_message (bx_message::info) << "theta: ep_pmts.size(): " << ep_pmts.size() << " graph_n: " << graph_n << dispatch;
			
	TGraph ep_graph_theta(graph_n, graph_theta, graph_time);
	TF1 ep_fit_theta("ep_fit_theta","[1]*(-1)*(-1+TMath::Sign(1,x-[0]))*(-1*sin((x-[0])/[2]))+[3]*(1+TMath::Sign(1,x-[0]))*(sin((x-[0])/[4]))+[5]",-pi,pi);
	ep_fit_theta.SetParameters(0,4.,2.,4.,2.,0);	
	ep_fit_theta.SetParLimits(1,0,1e6);
	ep_fit_theta.SetParLimits(0,-0.5,0.5);
	ep_fit_theta.SetParLimits(5,-10,10);
	ep_graph_theta.Fit(&ep_fit_theta,"RQ0","",-1.,1.);
	
	if (isnan(ep_fit_theta.GetParameter(0)) || isnan(ep_fit_theta.GetParError(0))) {
		//get_message (bx_message::warn) << "ep theta fit did not converge" << dispatch;
		return ev;
	}
	//get_message (bx_message::log) << "ep trans angles: phi " << ep_fit_phi.GetParameter(0) << ", theta " << ep_fit_theta.GetParameter(0) << dispatch; 

	// check for bad fits and double muon events
	if (ep_fit_phi.GetParameter(2)==0 || ep_fit_theta.GetParameter(2)==0 || ep_fit_phi.GetParameter(4)==0 || ep_fit_theta.GetParameter(4)==0) {
		//get_message (bx_message::warn) << "strange ep fit parameters, possibly double muon" << dispatch;
		return ev;		
	}
	if ( !(ep_fit_phi.GetParameter(1)/ep_fit_phi.GetParameter(2) > ep_par_threshold) || 
		 !(ep_fit_phi.GetParameter(3)/ep_fit_phi.GetParameter(4) > ep_par_threshold) ||
		 !(ep_fit_theta.GetParameter(1)/ep_fit_theta.GetParameter(2) > ep_par_threshold) ||
		 !(ep_fit_theta.GetParameter(3)/ep_fit_theta.GetParameter(4) > ep_theta_pos_threshold) ) {
		//get_message (bx_message::warn) << "strange ep fit parameters, possibly double muon" << dispatch;
		//get_message (bx_message::warn) << ep_fit_phi.GetParameter(1)/ep_fit_phi.GetParameter(2) << " " << ep_fit_phi.GetParameter(3)/ep_fit_phi.GetParameter(4) << " " << ep_fit_theta.GetParameter(1)/ep_fit_theta.GetParameter(2) << " " << ep_fit_theta.GetParameter(3)/ep_fit_theta.GetParameter(4) << dispatch;
		return ev;					
	}

	coo nep_trans = m_construct_xyz(mean_rad, Float_t(ep_fit_phi.GetParameter(0)-pi), Float_t(ep_fit_theta.GetParameter(0)+pi/2), Float_t(ep_fit_phi.GetParError(0)), Float_t(ep_fit_theta.GetParError(0)));

	//get_message (bx_message::log) << "new ep_trans: phi=" << ep_fit_phi.GetParameter(0) << "(" << ep_fit_phi.GetParError(0) << "), theta=" << ep_fit_theta.GetParameter(0) << "(" << ep_fit_theta.GetParError(0) << ")" << dispatch;




	coo nep = m_rotate_z(m_rotate_y(nep_trans, -(pi/2-oep.theta)), -(pi-oep.phi));

	//get_message (bx_message::log) << "old ep: x=" << oep.x << ", y=" << oep.y << ", z=" << oep.z << dispatch;
	//get_message (bx_message::log) << "new ep: x=" << nep.x << "(" << nep.dx << "), y=" << nep.y<< "(" << nep.dy << "), z=" << nep.z << "(" << nep.dz << ")" << dispatch;

	bx_track_by_points& t = ew.get_track_tof(); 
	t.f4_t1 = er.get_cluster(0).get_start_time() - er.get_trigger_rawt();
	t.f4_x1 = nep.x;
	t.f4_y1 = nep.y;
	t.f4_z1 = nep.z;
	t.f4_dx1 = fabs(nep.dx); 
	t.f4_dy1 = fabs(nep.dy);
	t.f4_dz1 = fabs(nep.dz);
	t.f4_t2 = 0;
	t.f4_x2 = 0;
	t.f4_y2 = 0;
	t.f4_z2 = 0;
	t.f4_dx2 = 0;
	t.f4_dy2 = 0;
	t.f4_dz2 = 0;
	t.f4_phi = 0;
	t.f4_theta = 0;
	t.f4_impact = 0;
	t.f4_labennormhits = -1;
        t.b_downward = true;
	ew.b_is_tracked_tof = false;
	
	if (lowhits) return ev;

	// EXIT POINT --------------------------------------------------
	
	coo ep_temp = m_normalize_vector(nep);
	coo ep_norm = m_construct_rpt(ep_temp.x, ep_temp.y, ep_temp.z);
	//get_message (bx_message::log) << "ep: phi " << ep_norm.phi << ", theta " << ep_norm.theta << dispatch;
	
	// DETERMINE XP PHI

	Float_t plain_imp[xp_plains], plain_dimp[xp_plains], plain_phi[xp_plains], plain_dphi[xp_plains], plain_c2[xp_plains], dev_c2[xp_plains];
	Int_t plain_counter=0;

	Float_t old_3 = pi;
	
	for (Int_t imp=0; imp<xp_plains; imp++) {
		
		Float_t plain_impact = (xp_plains-1)/2.-imp;
		// reuse old vectors for phi plains
		graph_n = 0;
		mean_rad = 0;
		//get_message (bx_message::log) << "imp " << imp << dispatch;
		
		// loop over hits
		for (UInt_t ihit = 0; ihit<lasthit && graph_n<constants::laben::channels; ihit++) {
			
			const bx_laben_clustered_hit& chit = er.get_cluster(0).get_clustered_hit(ihit);
			
			coo pmt;
			pmt.time = chit.get_time();
			//if (chit.get_order_in_channel() > 1) continue;  //AAA 0-based or 1-based ??? to be checked
			const db_channel_laben* c =chit.get_decoded_hit().get_db_channel();
			pmt.x = c->pmt_x();
			pmt.y = c->pmt_y();
			pmt.z = c->pmt_z();
			
			// tof cut
			Float_t tof = sqrt(pow(pmt.x-nep.x,2)+pow(pmt.y-nep.y,2)+pow(pmt.z-nep.z,2))/0.3;
			if (pmt.time<tof-5.) continue;
			if (pmt.time>tof*1.66+20.) continue;
			
			// plain cut
			Float_t plain_distance = pmt.x*ep_norm.x+pmt.y*ep_norm.y+pmt.z*ep_norm.z - plain_impact;
			if (fabs(plain_distance) > 0.5) continue; 
			
			pmt.rad = m_get_radius(pmt);
						
			pmt.x -= ep_norm.x*plain_impact;
			pmt.y -= ep_norm.y*plain_impact;
			pmt.z -= ep_norm.z*plain_impact;
			
			coo pmt_trans = m_rotate_y(m_rotate_z(pmt,-ep_norm.phi),ep_norm.theta);
			pmt_trans.time = pmt.time;
		
			pmt_trans.phi = m_limit_angle(atan2(pmt_trans.y, pmt_trans.x),0,2*pi);
			//pmt_trans.phi = atan2(pmt_trans.y, pmt_trans.x);
			//if (pmt_trans.phi<0) pmt_trans.phi += 2*pi;
			pmt_trans.rad = m_get_radius(pmt_trans);
			//get_message (bx_message::log) << imp << " => " << ihit << " : " << pmt_trans.phi << "--->" << pmt_trans.time << dispatch;	
			mean_rad += pmt_trans.rad;
			graph_phi[graph_n] = pmt_trans.phi;
			graph_time[graph_n] = pmt_trans.time;
			graph_n++;	
		}
		
		if (graph_n==0) continue;
		
		// define ring graph and fit
		mean_rad /= graph_n;		
		TGraph xp_graph_phi(graph_n, graph_phi, graph_time);
		TF1 xp_fit_phi("xp_fit_phi","0.3/[0]*sqrt([1]*[1]+[2]*[2]-2*[1]*[2]*(cos(x)*cos([3])+sin(x)*sin([3])))+[4]",0,2*pi);
		xp_fit_phi.SetParameters(1, mean_rad, mean_rad/2., old_3, 0);
		xp_fit_phi.SetParLimits(0,0,10);
		xp_fit_phi.FixParameter(1,mean_rad);
		xp_fit_phi.SetParLimits(2,0.1,mean_rad);
		xp_fit_phi.SetParLimits(3,-pi,3*pi);
		xp_graph_phi.Fit(&xp_fit_phi,"RQ0");
		Float_t xp_fit_phi_error = xp_fit_phi.GetParError(3);
		if (xp_fit_phi_error < 0.015) xp_fit_phi_error = 0.3;
		//if (xp_fit_phi.GetParameter(2) < 0.1) xp_fit_phi_error = 1e6;

		// fill in the graph of phis
		plain_imp[plain_counter]  = imp;
		plain_dimp[plain_counter] = 0;
		plain_phi[plain_counter]  = xp_fit_phi.GetParameter(3);	  	
		if (xp_fit_phi.GetNDF()>0) plain_c2[plain_counter] = xp_fit_phi.GetChisquare()/xp_fit_phi.GetNDF();
		else plain_c2[plain_counter] = 1e6;
		plain_dphi[plain_counter] = xp_fit_phi_error;
		//old_3 = xp_fit_phi.GetParameter(3);
		//get_message (bx_message::log) << imp << ":" << graph_n << " hits, " << plain_counter << " => plain phi : " << plain_phi[plain_counter] << "+/-" << plain_dphi[plain_counter] << " ; reduced c2 : " << plain_c2[plain_counter]<< dispatch;
		plain_counter++;
	}
	
	// check whether some of the chi2-values are substantially larger than others;
	// in this case, give phi of this plain a very large error bar
 	Float_t mean_c2 = 0;
  	for (Int_t p=0; p<plain_counter; p++) mean_c2 += plain_c2[p];
  	mean_c2 /= plain_counter;
  	Float_t sigma_c2=0;
  	for (Int_t p=0; p<plain_counter; p++) {
		dev_c2[p] = plain_c2[p]-mean_c2;
		sigma_c2 += dev_c2[p]*dev_c2[p];
	}
	sigma_c2 = sqrt(sigma_c2/plain_counter);
	for (Int_t p=0; p<plain_counter; p++) {
		if (dev_c2[p]>sigma_c2) {
			plain_dphi[p] = 1e6;
			//get_message(bx_message::log) << "imp " << p << " has dev " << dev_c2[p] << ">" << sigma_c2 << " => dphi set to 1e6" << dispatch;
		}
	}
	//get_message(bx_message::log) << "aa" << dispatch;
			
	// determine if many of the phis are larger than pi;
	// gives orientation for the later fit to xp theta;
	Float_t frac_larger_than_pi = 0;
	Float_t frac_counter = 0;
	for (Int_t p=0; p<plain_counter; p++) {
		//if ( plain_dphi[p] == 1e6 ) continue;
		frac_counter ++;
		if ( (plain_phi[p] > pi && plain_phi[p] < 2.*pi) || plain_phi[p]<0) frac_larger_than_pi++;
		plain_phi[p] = m_limit_angle(plain_phi[p],0,pi);
		//get_message (bx_message::log) << p << "=>" << plain_phi[p] << dispatch;
	}
	frac_larger_than_pi /= frac_counter;
	
	//get_message(bx_message::log) << plain_counter  << dispatch;
	
	if (plain_counter>xp_plains) plain_counter=xp_plains;
	if (plain_counter==0) return ev;
	// create and fit graph of phis to get xp phi
	TGraphErrors xp_graph_phis(plain_counter, plain_imp, plain_phi, plain_dimp, plain_dphi); 
	TF1 xp_fit_phis("xp_fit_phis","[0]",0,xp_plains);
	xp_fit_phis.SetParameter(0,pi);
	xp_graph_phis.Fit(&xp_fit_phis,"RQ0");
	//get_message(bx_message::log) << "first xp phi fit: " << xp_fit_phis.GetParameter(0) << dispatch;
	coo xp_trans;
	xp_trans.phi = xp_fit_phis.GetParameter(0);
	xp_trans.dphi = xp_fit_phis.GetParError(0);
	Float_t xp_fit_phis_chi2 = xp_fit_phis.GetChisquare();
	//if (drightangle > 2*TMath::Pi()) drightangle = 2*TMath::Pi(); ?
		
	// try again with shifted phi values
	for (Int_t p=0; p<plain_counter; p++) {
		if (plain_phi[p] < 0.9*pi) plain_phi[p] += pi;
		//get_message (bx_message::log) << p << "=>" << plain_phi[p] << dispatch;
	}  	
	TGraphErrors xp_graph_phis_shift(plain_counter, plain_imp, plain_phi, plain_dimp, plain_dphi); 
	xp_fit_phis.SetParameter(0,pi);
	xp_graph_phis_shift.Fit(&xp_fit_phis,"RQ0");
	//get_message(bx_message::log) << "second xp phi fit: " << xp_fit_phis.GetParameter(0) << dispatch;
	if (xp_fit_phis.GetChisquare() < xp_fit_phis_chi2) {
		xp_trans.phi = xp_fit_phis.GetParameter(0);
		xp_trans.dphi = xp_fit_phis.GetParError(0);
	}	
		
	// choose the right orientation for theta fit
	if (frac_larger_than_pi>0.5 && xp_trans.phi < pi ) xp_trans.phi += pi;
	if (frac_larger_than_pi<0.5 && xp_trans.phi > pi ) xp_trans.phi -= pi;

	if (isnan(xp_trans.phi)) {
		get_message (bx_message::warn) << "xp phi fit did not converge!" << dispatch;
		return ev;
	}
		
	//get_message (bx_message::log) << "xp phi " << xp_trans.phi << "(" << xp_trans.dphi << ")" << dispatch;	

	// DETERMINE XP THETA
	
	graph_n = 0;
	mean_rad = 0;

	// get normalized vector to define right phi plain
	coo xp_phi_norm;
	xp_phi_norm.x = -sin(xp_trans.phi);
	xp_phi_norm.y = cos(xp_trans.phi);
	xp_phi_norm.z = 0.;
		
	// loop over hits
	for (UInt_t ihit = 0; ihit<lasthit && graph_n<constants::laben::channels; ihit++) {
		
		const bx_laben_clustered_hit& chit = er.get_cluster(0).get_clustered_hit(ihit);
		
		coo pmt;
		pmt.time = chit.get_time();
		//if (chit.get_order_in_channel() > 1) continue;  //AAA 0-based or 1-based ??? to be checked
		const db_channel_laben* c = chit.get_decoded_hit().get_db_channel();
		pmt.x = c->pmt_x();
		pmt.y = c->pmt_y();
		pmt.z = c->pmt_z();
		pmt.rad = m_get_radius(pmt);

		// tof cut
		Float_t tof = sqrt(pow(pmt.x-nep.x,2)+pow(pmt.y-nep.y,2)+pow(pmt.z-nep.z,2))/0.3;
		if (pmt.time<tof-5.) continue;
		if (pmt.time>tof*1.66+10.) continue;
		
		// transform coordinates to ep zenith
		coo pmt_trans = m_rotate_y(m_rotate_z(pmt,-ep_norm.phi),ep_norm.theta);
	
		// cut on phi plain
		Float_t plain_distance = pmt_trans.x*xp_phi_norm.x + pmt_trans.y*xp_phi_norm.y + pmt_trans.z*xp_phi_norm.z;
		if (fabs(plain_distance)>xp_theta_diff) continue;
		
		if (pmt_trans.x*cos(xp_trans.phi)+pmt_trans.y*sin(xp_trans.phi)>0) graph_theta[graph_n] = atan2( sqrt(pmt_trans.x*pmt_trans.x+pmt_trans.y*pmt_trans.y), pmt_trans.z);
		else graph_theta[graph_n] = atan2( -sqrt(pmt_trans.x*pmt_trans.x+pmt_trans.y*pmt_trans.y), pmt_trans.z ) + 2*pi;

		graph_time[graph_n] = pmt.time;
		graph_n++;
		mean_rad += pmt.rad;
	}
	if (!graph_n) return ev; // AAA: Line added by DD to avoid crash (21/08/10)	
	mean_rad /= graph_n;

	TGraph xp_graph_theta(graph_n,graph_theta,graph_time);

	// find good starting parameters for the fit
	// by comparing nhits with impact parameters from old id tracking
	
	Float_t impact_angle = 0.05;
	if (et.get_impact() < 6.5) impact_angle = acos(et.get_impact()/6.5); // maybe a parameter to play with
	xp_trans.theta = 2.*impact_angle;
	
	// set fit range
	Float_t cone_angle_low = xp_trans.theta - 1.68;
	if (cone_angle_low < 0.) cone_angle_low = 0;
	Float_t cone_angle_high = xp_trans.theta + 1.68;
	if (cone_angle_high > 2*pi) cone_angle_high = 2*pi;

	if (er.get_decoded_nhits() > xp_nhits_diff) { // use realistic curved function
		TF1 xp_fit_theta("xp_fit_theta","[0]*sqrt(2.-2.*cos(x))/(0.3 * sqrt (1. - 1./[2]**2)) * ( sin(acos(1./[2]) - abs(x/2. - [1]/2.) )  + [2]*sin( abs(x/2. - [1]/2.))) + [3]",0.,1.05*pi);
		xp_fit_theta.SetParameters(mean_rad,xp_trans.theta,1.5,10.);
		xp_fit_theta.FixParameter(0,mean_rad);  // impact
		xp_fit_theta.SetParLimits(1,0,1.25*pi);    // exit angle
		xp_fit_theta.SetParLimits(2,1.2,1.8);      // refractive index
		xp_fit_theta.SetRange(cone_angle_low,cone_angle_high);
		xp_graph_theta.Fit(&xp_fit_theta,"RQ0");
		xp_trans.theta = xp_fit_theta.GetParameter(1); // first iteration

		// new fit range
		cone_angle_low = xp_trans.theta - 1.68;
		if (cone_angle_low < 0.) cone_angle_low = 0;
		cone_angle_high = xp_trans.theta + 1.68;
		if (cone_angle_high > 2*pi) cone_angle_high = 2*pi;
		xp_fit_theta.SetRange(cone_angle_low,cone_angle_high);
		xp_graph_theta.Fit(&xp_fit_theta,"RQ0");
		xp_trans.theta = xp_fit_theta.GetParameter(1); // second iteration

		// new fit range
		cone_angle_low = xp_trans.theta - 1.68;
		if (cone_angle_low < 0.) cone_angle_low = 0;
		cone_angle_high = xp_trans.theta + 1.68;
		if (cone_angle_high > 2*pi) cone_angle_high = 2*pi;
		xp_fit_theta.SetRange(cone_angle_low,cone_angle_high);
		xp_graph_theta.Fit(&xp_fit_theta,"RQ0");
		xp_trans.theta = xp_fit_theta.GetParameter(1); // third iteration
		xp_trans.dtheta = xp_fit_theta.GetParError(1);
	}
	else { // use effective pointed function
		
		// produce theta binned second time minimum graph
		Float_t binned_theta[xp_theta_bins], binned_minimum[xp_theta_bins], binned_2nd_minimum[xp_theta_bins];
		for (Int_t ib=0; ib<xp_theta_bins; ib++) {
			binned_theta[ib] = (ib+0.5)/xp_theta_bins*(2*pi);
			binned_minimum[ib] = binned_2nd_minimum[ib] = 100.;
		}
		for (Int_t ih=0; ih<graph_n; ih++) {
			Int_t ib = Int_t(graph_theta[ih]/(2*pi)*xp_theta_bins);
			if (graph_time[ih] < binned_minimum[ib]) {
				binned_2nd_minimum[ib] = binned_minimum[ib];
				binned_minimum[ib]     = graph_time[ih];
				continue;
			}
			if (graph_time[ih] < binned_2nd_minimum[ib]) binned_2nd_minimum[ib] = graph_time[ih]; 
		}
		// if second hit is very late (because there is low hits in this bin), use the first one
		for (Int_t ib=0; ib<xp_theta_bins; ib++) {
			if (binned_2nd_minimum[ib]-15.>binned_minimum[ib] && binned_minimum[ib]>0) binned_2nd_minimum[ib] = binned_minimum[ib];
		}
		// if a bin is empty, cut it from the graph
		Int_t xp_real_bins = xp_theta_bins;
		for (Int_t ib=0; ib<xp_real_bins; ib++) {
			if (binned_2nd_minimum[ib]==100) {
				xp_real_bins--;
				for (Int_t jb=ib; jb<xp_real_bins; jb++) {
					binned_2nd_minimum[jb] = binned_2nd_minimum[jb+1];
					binned_theta[jb] = binned_theta[jb+1];
				}
				ib--;
			}
		}

		if (!xp_real_bins) return ev;

		// fit
		TGraph xp_graph_binned(xp_real_bins, binned_theta, binned_2nd_minimum);
		TF1 xp_fit_binned("xp_fit_binned","[0]*( sqrt(sin([3]*(x-[1]))) )+[2]",0,3);
		xp_fit_binned.SetParameters(1,1,0,0.5);
		xp_fit_binned.SetParLimits(0,0,100);
		xp_fit_binned.SetParLimits(1,0,100);
		xp_fit_binned.SetParLimits(3,0,100);
		xp_graph_binned.Fit(&xp_fit_binned,"RQ0");

		xp_trans.theta = xp_fit_binned.GetParameter(1);
		xp_trans.dtheta = xp_fit_binned.GetParError(1);
	}

	if (isnan(xp_trans.theta)) {
		get_message (bx_message::warn) << "xp theta fit did not converge!" << dispatch;
		return ev;	
	}

	//get_message (bx_message::log) << "xp theta " << xp_trans.theta << "(" << xp_trans.dtheta << ")" << dispatch;

	// build xp and backtransform
	xp_trans = m_construct_xyz(mean_rad, xp_trans.phi, xp_trans.theta, xp_trans.dphi, xp_trans.dtheta);
	coo xp = m_rotate_z(m_rotate_y(xp_trans,-ep_norm.theta),ep_norm.phi);
	
	// refit phi for xp
	
	graph_n = 0;

	// get normalized vector to define right theta plain
	coo xp_theta_norm = m_construct_xyz(1., xp_trans.phi, xp_trans.theta-pi/2, 0, 0);
			
	// loop over hits
	// select hits on a circle perpendicular to the one defined by ep and xp
	
	for (UInt_t ihit = 0; ihit<lasthit && graph_n<constants::laben::channels; ihit++) {
		
		const bx_laben_clustered_hit& chit = er.get_cluster(0).get_clustered_hit(ihit);
		
		coo pmt;
		pmt.time = chit.get_time();
		//if (chit.get_order_in_channel() > 1) continue;  //AAA 0-based or 1-based ??? to be checked
		const db_channel_laben* c = chit.get_decoded_hit().get_db_channel();
		pmt.x = c->pmt_x();
		pmt.y = c->pmt_y();
		pmt.z = c->pmt_z();
		pmt.rad = m_get_radius(pmt);

		// tof cut
		Float_t tof = sqrt(pow(pmt.x-nep.x,2)+pow(pmt.y-nep.y,2)+pow(pmt.z-nep.z,2))/0.3;
		if (pmt.time<tof-5.) continue;
		if (pmt.time>tof*1.66+10.) continue;
		
		// transform coordinates to ep zenith
		coo pmt_trans = m_rotate_y(m_rotate_z(pmt,-ep_norm.phi),ep_norm.theta);
	
		// cut on phi plain
		Float_t plain_distance = pmt_trans.x*xp_theta_norm.x + pmt_trans.y*xp_theta_norm.y + pmt_trans.z*xp_theta_norm.z;
		if (fabs(plain_distance)>xp_theta_diff) continue;
		
		coo pmt_2nd_trans = m_rotate_y(m_rotate_z(pmt_trans,-xp_trans.phi),xp_trans.theta-pi/2);
		graph_phi[graph_n] = m_limit_angle(atan2(pmt_2nd_trans.y, pmt_2nd_trans.x),-pi,pi);
		graph_time[graph_n] = pmt.time-tof;
		graph_n++;
	}
	

	if (!graph_n) return ev; // AAA: Line added by DD to avoid crash (21/08/10)
	TGraph xp_graph_2nd_phi(graph_n,graph_phi,graph_time);
	TF1 xp_fit_2nd_phi("xp_fit_2nd_phi","[1]*(-1)*(-1+TMath::Sign(1,x-[0]))*(-1*sin((x-[0])/[2]))+[3]*(1+TMath::Sign(1,x-[0]))*(sin((x-[0])/[4]))+[5]",-pi,pi);
	xp_fit_2nd_phi.SetParameters(0,4,2.,4.,2.,0);
	xp_fit_2nd_phi.SetParLimits(1,0,1e6);
	xp_fit_2nd_phi.SetParLimits(0,-0.5,0.5);
	xp_fit_2nd_phi.SetParLimits(2,0,1e6);
   	xp_fit_2nd_phi.SetParLimits(3,0,1e6);
        xp_fit_2nd_phi.SetParLimits(5,-10,10);
	//xp_fit_2nd_phi.FixParameter(5,0);
	xp_graph_2nd_phi.Fit(&xp_fit_2nd_phi,"RQ0","",-pi,pi);
	
	if (isnan(xp_fit_2nd_phi.GetParameter(0))) {
		get_message (bx_message::warn) << "2nd xp phi fit did not converge!" << dispatch;
		return ev;	
	}

	Float_t xp_phi_trans_trans = xp_fit_2nd_phi.GetParameter(0); 

	//get_message (bx_message::log) << "xp 2nd phi " << xp_phi_trans_trans << "(" << xp_fit_2nd_phi.GetParError(0) << ")" << dispatch;

	coo xp_trans_trans = m_construct_xyz(mean_rad, xp_fit_2nd_phi.GetParameter(0), pi/2, xp_fit_2nd_phi.GetParError(0), xp_trans.dtheta);

	xp_trans = m_rotate_z(m_rotate_y(xp_trans_trans,-xp_trans.theta+pi/2),xp_trans.phi);

	if (fabs(xp_phi_trans_trans)<0.5) {
		xp = m_rotate_z(m_rotate_y(xp_trans,-ep_norm.theta),ep_norm.phi);
		//get_message (bx_message::log) << "new xp: x=" << xp.x << "(" << xp.dx << "), y=" << xp.y << "(" << xp.dy << "), z=" << xp.z << "(" << xp.dz << ")" << dispatch;
	}

	// construct auxiliary parameters
	Float_t dex = sqrt(pow(nep.x-xp.x,2)+pow(nep.y-xp.y,2)+pow(nep.z-xp.z,2));
	Float_t cosalpha = (-nep.x*(xp.x-nep.x)-nep.y*(xp.y-nep.y)-nep.z*(xp.z-nep.z))/(dex*nep.rad);
	coo pedal;
	pedal.x = nep.x + (xp.x-nep.x)*nep.rad*cosalpha/dex;
	pedal.y = nep.y + (xp.y-nep.y)*nep.rad*cosalpha/dex;
	pedal.z = nep.z + (xp.z-nep.z)*nep.rad*cosalpha/dex;
	pedal.rad = m_get_radius(pedal);


	t.f4_t2 = dex/0.3;
	t.f4_x2 = xp.x;
	t.f4_y2 = xp.y;
	t.f4_z2 = xp.z;
	t.f4_dx2 = fabs(xp.dx);
	t.f4_dy2 = fabs(xp.dy);
	t.f4_dz2 = fabs(xp.dz);
	t.f4_phi = m_limit_angle(atan2(nep.y-xp.y, nep.x-xp.x),0,2*pi);
	t.f4_theta = acos((nep.z-xp.z)/dex);
	t.f4_impact = pedal.rad;
	t.b_downward = (nep.z>xp.z);
	if (ev->get_laben().get_decoded_nhits() > 2250) ew.b_is_tracked_tof = true;

        // write normalized number of decoded hits into track variable
        if ( ev->get_laben().get_nclusters()>0 )
           t.f4_labennormhits = ev->get_laben().get_decoded_nhits() / (ev->get_laben().get_n_live_pmts() - ev->get_laben().get_invalid_pmts()) * 2000.;
//	get_message (bx_message::log) << "old xp: x=" << et.get_x2() << ", y=" << et.get_y2() << ", z=" << et.get_z2() << dispatch;

	i4_tracked_events++;
	return ev;
}

void bx_laben_tof_tracker::end () {
	get_message (bx_message::info) << i4_tracked_events << " of " << i4_laben_events << " were tracked." << dispatch;
}

// Private method(s)

coo bx_laben_tof_tracker::m_rotate_x(coo input, Float_t phi) {
	coo output;
	output.x = input.x;
	output.y = input.y*cos(phi)-input.z*sin(phi);
	output.z = input.y*sin(phi)+input.z*cos(phi);
	output.dx = input.dx;
	output.dy = input.dy*cos(phi)-input.dz*sin(phi);
	output.dz = input.dy*sin(phi)+input.dz*cos(phi);	
	output.rad = input.rad;
	output.time = input.time;	
	return output;	
}

coo bx_laben_tof_tracker::m_rotate_y(coo input, Float_t phi) {
	coo output;
	output.x = input.x*cos(phi)-input.z*sin(phi);
	output.y = input.y;
	output.z = input.x*sin(phi)+input.z*cos(phi);			
	output.dx = input.dx*cos(phi)-input.dz*sin(phi);
	output.dy = input.dy;
	output.dz = input.dx*sin(phi)+input.dz*cos(phi);		
	output.rad = input.rad;
	output.time = input.time;	
	return output;	
}

coo bx_laben_tof_tracker::m_rotate_z(coo input, Float_t phi) {
	coo output;
	output.x = input.x*cos(phi)-input.y*sin(phi);
	output.y = input.x*sin(phi)+input.y*cos(phi);
	output.z = input.z;	
	output.dx = input.dx*cos(phi)-input.dy*sin(phi);
	output.dy = input.dx*sin(phi)+input.dy*cos(phi);
	output.dz = input.dz;	
	output.rad = input.rad;
	output.time = input.time;	
	return output;	
}
							   
coo bx_laben_tof_tracker::m_construct_xyz(Float_t rad, Float_t phi, Float_t theta, Float_t dphi, Float_t dtheta) {
	coo output;
	output.x = rad*cos(phi)*sin(theta);
	output.y = rad*sin(phi)*sin(theta);
	output.z = rad*cos(theta);
	output.rad = rad;
	output.dx = rad*sqrt(pow(sin(phi)*sin(theta)*dphi,2)+pow(cos(phi)*cos(theta)*dtheta,2));
	output.dy = rad*sqrt(pow(cos(phi)*sin(theta)*dphi,2)+pow(sin(phi)*cos(theta)*dtheta,2));
	output.dz = rad*sin(theta)*dtheta;
	output.phi = phi;
	output.theta = theta;
	output.dphi = dphi;
	output.dtheta = dtheta;
	output.time = 0;
	return output;
}

coo bx_laben_tof_tracker::m_construct_rpt(Float_t x, Float_t y, Float_t z) {
	coo output;
	output.x = x;
	output.y = y;
	output.z = z;
	output.rad = sqrt(x*x+y*y+z*z);
	output.dx = 0;
	output.dy = 0;
	output.dz = 0;
	output.phi = atan2(y, x);
	output.theta = acos(z/output.rad);
	output.time = 0;
	return output;
}

Float_t bx_laben_tof_tracker::m_limit_angle(Float_t phi, Float_t min_angle, Float_t max_angle) {
	Float_t corrector = max_angle-min_angle;
	while (phi<min_angle) phi += corrector;
	while (phi>max_angle) phi -= corrector;
	return phi;
}

Float_t bx_laben_tof_tracker::m_get_radius(coo input) {
	return sqrt(input.x*input.x+input.y*input.y+input.z*input.z);
}

Float_t bx_laben_tof_tracker::m_get_radius(Float_t x, Float_t y, Float_t z=0) {
	return sqrt(x*x+y*y+z*z);
}


coo bx_laben_tof_tracker::m_normalize_vector(coo input) {
	Float_t factor = 1./input.rad;
	coo output;
	output.x = input.x*factor;
	output.y = input.y*factor;
	output.z = input.z*factor;
	output.rad = input.rad*factor;
	output.dx = input.dx*factor;
	output.dy = input.dy*factor;
	output.dz = input.dz*factor;
	output.theta = input.theta;
	output.phi = input.phi;
	output.time = input.time;
	return output;	
}

//bool bx_laben_tof_tracker::m_check_this_and_that(const bx_laben_event& er) {
  // ...
  //return true;
//}

/*
 * $Log: bx_laben_tof_tracker.cc,v $
 * Revision 1.19  2014/12/11 21:27:12  wurm
 * added getters and variables to BxEvent that concern C-11 likelihood tagging
 * * for tracks: new distance-to- and projection-on-track getters, dE/dx calculation
 * * counting of neutron-like clusters in tt128, new flag of neutron clusters
 *
 * Revision 1.18  2011/03/09 17:05:56  wurm
 * save definition of graphs
 *
 * Revision 1.17  2011-03-09 10:24:26  wurm
 * put some safety to initializing TGraphs
 *
 * Revision 1.16  2010-08-24 14:46:27  wurm
 * introduced 1000 hits reconstruction threshold, fixed counter
 *
 * Revision 1.14.2.1  2010-08-24 13:46:09  wurm
 * introduced a 1000 hits threshold for reconstruction, fixed tracked counter
 *
 * Revision 1.14  2010-08-21 02:52:42  ddangelo
 * debugging: added a few checks on number of elements > 0 before creating TGraph objs
 *
 * Revision 1.13  2010-08-10 16:18:24  wurm
 * adjusted the xp phi fit for virginia
 *
 * Revision 1.11  2010-08-06 15:17:19  wurm
 * set parameters for phi fit
 *
 * Revision 1.10  2010-08-06 11:10:53  wurm
 * increased verbosity
 *
 * Revision 1.9  2010-08-05 15:23:08  wurm
 * changed the conditions for laben.is_tracked_tof, adjusted global tracker
 *
 * Revision 1.8  2010-08-04 13:40:19  wurm
 * muted warning messages
 *
 * Revision 1.7  2010-08-04 08:29:35  wurm
 * removed unneccesary variables, set uncertainties for global tracker to 0.5 m
 *
 * Revision 1.6  2010-08-03 15:57:34  wurm
 * 1) Introduced theta, phi, impact variables for fitted tracks
 * 2) set fixed uncertainties for global tracking
 * 3) condition for execution of laben tracking changed to MTF OR MCF
 * 4) lowered verbosity of all modules
 *
 * Revision 1.5  2010-07-30 17:06:59  wurm
 * require EP and XP for is_tracked variable to be set true
 *
 * Revision 1.4  2010-07-01 09:52:52  wurm
 * removed small bugs from xp theta/phi fits
 *
 * Revision 1.3  2010-06-28 10:44:03  wurm
 * new laben tof tracking, modifications to muon tracking (Introducing phi,theta,impact to the event) and modification of global tracker to prefer laben tof over laben energy (and conditions)
 *
 * Revision 1.2  2010-05-21 13:17:37  ddangelo
 * adding laben_tof_tracker and cmt_tracker
 *
 * Revision 1.1  2010-05-21 12:34:01  ddangelo
 * bx_laben_tracker renamed as bx_laben_energy_tracker
 * added bx_laben_tof_tracker and bx_cmt_tracker
 *
 * Revision 1.15  2006-08-21 11:19:02  razeto
 * Updated to new barn_interface
 *
 * Revision 1.14  2006/01/09 16:21:03  razeto
 * Updated to the new root_barn target
 *
 * Revision 1.13  2005/02/10 17:37:39  ddangelo
 * removed afew compilation warning
 *
 * Revision 1.12  2005/02/03 19:00:18  ddangelo
 * maintainer changed (Razeto->D'Angelo)
 * Module really implemented with many possible examples for developers.
 * It features all aspects described in "Module's programmers guide" in the docs.
 *
 * Revision 1.11  2004/11/26 15:25:11  razeto
 * Thanks to Davide who so kindly underlined a syntactical error, I fixed the "Maintainer" word
 *
 * Revision 1.10  2004/11/26 14:01:21  razeto
 * Added Mantainer field
 *
 * Revision 1.9  2004/09/22 13:28:37  razeto
 * Changed bx_reco_event to bx_echidna_event
 *
 * Revision 1.8  2004/04/06 12:42:17  razeto
 * Added module role to base_module, and the code to handle everywhere. Some other minor fixes
 *
 * Revision 1.7  2004/04/03 09:20:35  razeto
 * Added messenger
 *
 * Revision 1.6  2004/04/02 13:40:01  razeto
 * Added messenger
 *
 * Revision 1.5  2004/03/21 18:55:05  razeto
 * Some cosmetic changes
 *
 * Revision 1.4  2004/03/20 19:10:03  pallas
 * Debugging
 *
 * Revision 1.3  2004/03/20 19:05:17  pallas
 * Debugging
 *
 * Revision 1.2  2004/03/20 18:56:46  pallas
 * Change counter
 *
 * Revision 1.1  2004/03/20 17:39:45  pallas
 * First release
 *
 * Revision 1.1  2004/03/20 14:16:47  razeto
 * Added sources
 *
 */
