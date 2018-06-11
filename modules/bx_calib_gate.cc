/* BOREXINO Reconstruction program
 *
 * Author: Barbara Caccianiga <Barbara.Caccianiga@mi.infn.it>
 * Maintainer: Barbara Caccianiga <Barbara.Caccianiga@mi.infn.it>
 * 
 *
 *
 */
#include "bx_calib_gate.hh"
#include "bx_echidna_event.hh"
#include "bx_trigger_event.hh"
#include "barn_interface.hh"
#include "db_profile.hh"
#include "db_calib.hh"
#include "db_channel.hh"
#include "db_run.hh"
#include "bx_dbi.hh"
#include "bx_trigger_event.hh"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "constants.hh"
#include <math.h>
#include <sstream>
#include <iostream>


namespace { 
  double myfitfunction(double *x, double *par) {
    double function ;
    if(x[0] < par[0]) {
      function = par[2] ;  
    } else if (x[0] >= par[0] && x[0] <= par[1] ) {
      function = par[3] ;
    } else {
      function = par[4] ;
    }
    return function ;
  }
};


// ctor
bx_calib_gate::bx_calib_gate():bx_base_module("bx_calib_gate", bx_base_module::main_loop) {
  i4_times = 0;
  require_event_stage(bx_detector::laben, bx_base_event::decoded);
   
  i4_bin = 1000;
  up_limit = 15000.;
  low_limit = -15000.;   
}

// module interface
void bx_calib_gate::begin() {
  get_message(bx_message::debug) << "bx_calib_gate " << dispatch;

  p_tempo_random   = new TH1F("p_tempo_random",   "time random",   i4_bin,low_limit,up_limit);
  p_tempo_laser    = new TH1F("p_tempo_laser",    "time laser",    i4_bin,low_limit,up_limit);
  p_tempo_pulser   = new TH1F("p_tempo_pulser",   "time pulser",   i4_bin,low_limit,up_limit);
  p_tempo_neutrino = new TH1F("p_tempo_neutrino", "time neutrino", i4_bin,low_limit,up_limit);

  barn_interface::get()->store(barn_interface::file,p_tempo_random,this);
  barn_interface::get()->store(barn_interface::file,p_tempo_laser,this);
  barn_interface::get()->store(barn_interface::file,p_tempo_pulser,this);
  barn_interface::get()->store(barn_interface::file,p_tempo_neutrino,this);
}

bx_echidna_event* bx_calib_gate::doit (bx_echidna_event *ev) {
  const bx_laben_event& er = ev->get_laben();
  i4_times = i4_times+1;  
  int32_t nhits = er.get_decoded_nhits();
  for(int32_t i = 0; i < nhits; i++){
     const bx_laben_decoded_hit& hit = er.get_decoded_hit(i);
     const bx_laben_raw_hit& rhit = hit.get_raw_hit();
     double t = hit.get_raw_time() - er.get_trigger_rawt();
     int32_t lg = rhit.get_logical_channel (); 
     if(bx_dbi::get ()->get_channel (lg).is_ordinary ()) {
       if(ev->get_trigger().is_random())
	 p_tempo_random->Fill(t);
       if(ev->get_trigger().is_laser394()) 
	 p_tempo_laser->Fill(t);   
       if(ev->get_trigger().is_pulser()) 
	 p_tempo_pulser->Fill(t);  
       if(ev->get_trigger().is_neutrino()) 
	 p_tempo_neutrino->Fill(t);
     }
  }
  
  return ev;
}

void bx_calib_gate::end () {

  if(p_tempo_neutrino->GetEntries() && p_tempo_laser->GetEntries() && p_tempo_pulser->GetEntries() && p_tempo_random->GetEntries() >
  2000 ) {
    

    float maxbin_laser    = p_tempo_laser->GetMaximumBin()*(up_limit - low_limit)/i4_bin + low_limit;
    float maxbin_pulser   = p_tempo_pulser->GetMaximumBin()*(up_limit - low_limit)/i4_bin + low_limit;
    float maxbin_neutrino = p_tempo_neutrino->GetMaximumBin()*(up_limit - low_limit)/i4_bin + low_limit;

    TF1 *func = new TF1("myfitfunction",myfitfunction,low_limit,up_limit,5);
    func->SetParameter(0,-6600);
    func->SetParameter(1,500);    
    func->SetParameter(2,p_tempo_random->Integral()/i4_bin/10.);
    func->SetParameter(3,p_tempo_random->Integral()/i4_bin*2.5);
    func->SetParameter(4,0);
    p_tempo_random->Fit(func,"QN");
    p_tempo_random->Fit(func,"QN");
    p_tempo_random->Fit(func,"QN");
    
    db_run& run_info = bx_dbi::get ()->get_run ();

    float start_of_gate   = func->GetParameter(0) ;
    float end_of_gate     = func->GetParameter(1) ;
    float width_of_gate   = end_of_gate - start_of_gate; 
    
    float pre_width   = run_info.get_laben_gate_width();	     
    float pre_start   = run_info.get_laben_gate_start();
    float pre_laser   = run_info.get_laben_laser_offset();
    float pre_pulser  = run_info.get_laben_pulser_offset();
    float pre_cluster = run_info.get_laben_cluster_offset();
    int32_t write_on_db = 1;
       
    //difference old and new values in ns  
    float c_width   =    fabs(pre_width - width_of_gate);
    float c_start   =    fabs(pre_start - start_of_gate);
    float c_laser   =    fabs(pre_laser - maxbin_laser);
    float c_pulser  =    fabs(pre_pulser - maxbin_pulser);
    float c_cluster =    fabs(pre_cluster - maxbin_neutrino);
    

    if(c_width > 250 || c_start > 150 || c_laser > 250 || c_pulser  > 250 || c_cluster > 300 ) {
      write_on_db = 0;
      get_message(bx_message::error)<< "Trigger properties are too different from the last calibration; the database will be not updated unless you set db_write = 2" << dispatch;
      if(c_width > 250)  get_message(bx_message::warn)<< "Gate width differs by " << c_width << " ns " << dispatch;
      if(c_start > 150)  get_message(bx_message::warn)<< "Gate start differs by " << c_start << " ns " << dispatch;
      if(c_laser > 250)  get_message(bx_message::warn)<< "Laser position differs by " << c_laser << " ns " << dispatch;
      if(c_pulser > 250) get_message(bx_message::warn)<< "Pulser position differs by " << c_pulser << " ns " << dispatch;
      if(c_cluster > 300)get_message(bx_message::warn)<< "Cluster position differs by " << c_cluster << " ns " << dispatch;

    }
    
    if(get_parameter ("db_write").get_int () == 2) write_on_db = 1 ;
    if(!run_info.is_trigger_parameters_present () && get_parameter ("db_write").get_int () && write_on_db){ 


      run_info.set_laben_gate_width	   (func->GetParameter(1) - func->GetParameter(0), this);
      run_info.set_laben_gate_start	   (func->GetParameter(0), this);
      run_info.set_laben_laser_offset      (maxbin_laser, this);
      run_info.set_laben_pulser_offset     (maxbin_pulser, this);
      run_info.set_laben_cluster_offset    (maxbin_neutrino, this);
   
      run_info.write_trigger_parameters(true,this);
      get_message(bx_message::info)<< "Writing trigger parameters to DB" << dispatch;
 
    }   
    get_message(bx_message::log)<<"Start of the gate         =  "   << func->GetParameter(0) << "  ns  "  << dispatch;
    get_message(bx_message::log)<<"End of the gate           =  "   << func->GetParameter(1) << "  ns  "  << dispatch;
    get_message(bx_message::log)<<"Width of the gate         =  "   << func->GetParameter(1) - func->GetParameter(0)<< "  ns  "  << dispatch;
    get_message(bx_message::log)<<"Cluster Position =  "   << maxbin_neutrino << "  ns  "  << dispatch;
    get_message(bx_message::log)<<"Laser Position    =  "   << maxbin_laser  << "  ns  "     << dispatch;
    get_message(bx_message::log)<<"Pulser Position   =  "   << maxbin_pulser << "  ns  "   << dispatch;
  

  } else {
    get_message(bx_message::warn)<<"Not enough statistics for evaluating gate properties"   << dispatch;  
  }
}

