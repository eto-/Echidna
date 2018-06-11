/*
 * Author: Juergen Winter <juergen.winter@ph.tum.de>
 * Maintainer: Juergen Winter <juergen.winter@ph.tum.de>
 *
 * Module for tagging misalignment between muon data and mcr
*/

#include <cmath>
#include <iostream>
#include <vector>
#include <algorithm>
#include "bx_calib_muon_alignment.hh"
#include "messenger.hh" 
#include "bx_echidna_event.hh"
#include "bx_dbi.hh"
#include "db_profile.hh"
#include "db_channel.hh"
#include "db_run.hh"
#include "barn_interface.hh"
#include "TH1F.h"
#include "TF1.h"

bx_calib_muon_alignment::bx_calib_muon_alignment() : bx_base_module("bx_calib_muon_alignment", bx_base_module::main_loop) {
  //require_event_stage (bx_detector::muon, bx_base_event::decoded);
//  require_trigger_type (bx_trigger_event::pulser);  
}

  
//BEGIN
void bx_calib_muon_alignment::begin () {


	//initialisiation
	nmis_thresh  = get_parameter("nmis_thresh").get_int ();   //cut on number of subsequent bad pulser events indicating a misalignment
	nhits_thresh = get_parameter("nhits_thresh").get_int ();  //cut on pulser events with too few decoded nhits  

	b_check_last = true;	//false== misalignment, true==no misalignment
	b_check_mcr  = true;    // true== mcr up
	
	num =0;

	u4_evnum_up=0;
	u4_evnum_aligned=0;
	u4_evnum=0;

	time_up=0;
	time_aligned=0;
	time=0;

}


// DOIT
bx_echidna_event* bx_calib_muon_alignment::doit (bx_echidna_event *ev) {
	
	const bx_muon_event& er = ev->get_muon (); 

	u4_evnum = ev->get_event_number();
	time  = ev->get_trigger().get_time_t();

	//mcr down?
	if(!(ev->is_muon_enabled ()) && b_check_mcr){
	        b_check_mcr=false;
	        get_message (bx_message::log) << "mcr down "<< dispatch;
	}

	if(b_check_mcr){
		u4_evnum_up = u4_evnum;
		time_up = time;
	}

        if(!ev->get_trigger().is_pulser()) return ev;

	// subsequent misaligned events?
        if(b_check_last && b_check_mcr){
       		if(er.get_decoded_nhits () >  nhits_thresh ){
			u4_evnum_aligned = u4_evnum;
			time_aligned = time;
			num=0;
		}
		else{ 	num++;
			get_message (bx_message::log) << "misaligned event, ev no "<< u4_evnum<< "     num "<< num <<dispatch;
		}		

		if(num==nmis_thresh)  b_check_last=false;
	}

return ev;
}

// end
void bx_calib_muon_alignment::end () {

// writing to visitors

        db_run& run_info = bx_dbi::get()->get_run ();
        time_t		time_tot     = time;
        uint32_t 	u4_evnum_tot = u4_evnum;

	if(b_check_last && b_check_mcr){
		u4_evnum_up = u4_evnum_tot;
		time_up     = time_tot;
		u4_evnum_aligned = u4_evnum_tot;
		time_aligned     = time_tot;
		get_message(bx_message::log) << "no misalignment "  << dispatch;
	}


	else{	get_message(bx_message::warn) << "misalignment "  << dispatch;
		
		//no mcr crash
		if(b_check_mcr){
			u4_evnum_up=u4_evnum_tot;
		        time_up = time_tot;
		}
		
		//mcr crash without misalignment before
		if(!b_check_mcr && b_check_last) {
			u4_evnum_aligned=u4_evnum_up;
			time_aligned = time_up;
		}
           }


// setting visitors
	run_info.set_muon_alignment (u4_evnum_tot,time_tot,u4_evnum_up,time_up,u4_evnum_aligned,time_aligned,this); 

// write to DB
	if(get_parameter ("db_write").get_bool () ){
	  run_info.write_muon_alignment (true, this);
	}

}
