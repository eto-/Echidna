/*
Author: Timo Lewke <timo.lewke@ph.tum.de>, Michael Wurm <mwurm@ph.tum.de>
Maintainer: Timo lewke <timo.lewke@ph.tum.de>

Module for muon veto pulser, to distinguish good and bad channels
*/

#include <cmath>
#include <iostream>
#include <vector>
#include <algorithm>
#include "bx_calib_muon_pulser.hh"
#include "bx_echidna_event.hh"
#include "bx_dbi.hh"
#include "db_profile.hh"
#include "db_channel.hh"
#include "db_run.hh"
#include "barn_interface.hh"
#include "TH1F.h"
#include "TF1.h"

bx_calib_muon_pulser::bx_calib_muon_pulser() : bx_base_module("bx_calib_muon_pulser", bx_base_module::main_loop) {
  require_event_stage (bx_detector::muon, bx_base_event::decoded);
  require_trigger_type (bx_trigger_event::pulser);
  
}

  
//BEGIN
void bx_calib_muon_pulser::begin () {
	
//	int32_t lg = constants::muon::channels;

	// histograms for barn
	pulser_hits_vs_lg = new TH1F("pulser_hits_vs_mch","decoded hits", constants::muon::channels,0, constants::muon::channels);
	working_channels = new TH1F("dead_channels","Broken QTC-Channels", constants::muon::channels,0, constants::muon::channels);
	
	pulser_hits_vs_lg->SetXTitle("Muon Channel");
	pulser_hits_vs_lg->SetYTitle("Hits");
	working_channels->SetXTitle("Muon Channel");
	working_channels->SetYTitle("Dead");	
	
	barn_interface::get ()->store (barn_interface::file, pulser_hits_vs_lg, this);
	barn_interface::get ()->store (barn_interface::file, working_channels , this);
		
}


// DOIT
bx_echidna_event* bx_calib_muon_pulser::doit (bx_echidna_event *ev) {
	
	const bx_muon_event& er = ev->get_muon (); 
	 
	//loop on hits
	for(int32_t i = 0; i < er.get_decoded_nhits (); i++) {
	float lg = er.get_decoded_hit (i).get_raw_hit().get_muon_channel ();
	pulser_hits_vs_lg -> Fill (lg);
	}
	
	for(int32_t lg = 0; lg < constants::muon::channels; lg++){
	double counter = pulser_hits_vs_lg -> GetBinContent (lg);
	int32_t status; 
	if (counter == 1000) status = 1; 
	else status = 0;
	working_channels -> SetBinContent (lg,status);	
	}
	
	for (int32_t lg = 0; lg < constants::muon::channels; lg++ ){
	a[lg]=0;
	double aus = working_channels -> GetBinContent (lg);
	if(aus == 1) a[lg]++; 
	}
		
	return ev;

}

// end
void bx_calib_muon_pulser::end () {
for (int32_t lg = 0; lg < constants::muon::channels; lg++ ){
	if (a[lg]!=0) get_message(bx_message::warn) << " Dead QTC-Channel in mch: " << lg-1 << dispatch;	
	}
}

