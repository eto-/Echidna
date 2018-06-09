/* BOREXINO Reconstruction program
 *
 * Author: Elena Guardincerri <elena.guardincerri@ge.infn.it>
 *
 * Maintainer: Elena Guardincerri <elena.guardincerri@ge.infn.it>
 *
 * $Id: bx_calib_monitor.cc,v 1.24 2011/03/02 18:51:57 guardi Exp $
 *
 */

#include "bx_calib_monitor.hh"
#include "bx_echidna_event.hh"
#include "barn_interface.hh"
#include "db_run.hh"
#include "db_profile.hh"
#include "bx_dbi.hh"
#include "TH1F.h"
#include "TH2F.h"

// ctor
bx_calib_monitor::bx_calib_monitor (): bx_base_module("bx_calib_monitor", bx_base_module::main_loop) {

  require_event_stage (bx_detector::laben, bx_base_event::reconstructed);
  // Option: you can require to run only on events with a given trigger type
  require_trigger_type (bx_trigger_event::neutrino);
}

bx_calib_monitor::~bx_calib_monitor (){
  m_delete_histos();
  delete neutrino_z_t_14C;
  delete neutrino_z_t_Bi;
  delete neutrino_x2y2_t_14C;
  delete neutrino_x2y2_t_Bi;
  delete neutrino_x_t_14C;
  delete neutrino_x_t_Bi;
  delete neutrino_y_t_14C;
  delete neutrino_y_t_Bi;
}

void bx_calib_monitor::begin () {

  get_message (bx_message::debug) << "bx_calib_monitor begin" << dispatch;

  f_buffer_time =  get_parameter("time_interval").get_int();

  // Create 1-D histograms
  neutrino_z_t_14C = new TH1F("_z_t_C"," z versus t for 14C events", 1000, -6.0,6.0);  
  neutrino_z_t_14C->GetXaxis()->SetTitle("time (s)");
  neutrino_z_t_14C->GetYaxis()->SetTitle("z (m)");
  neutrino_z_t_Bi = new TH1F("_z_t_Bi"," z versus t for BiPo events", 1000, -6.0,6.0);  
  neutrino_z_t_Bi->GetXaxis()->SetTitle("time (s)");
  neutrino_z_t_Bi->GetYaxis()->SetTitle("z (m)");
  neutrino_x2y2_t_14C = new TH1F("_x2y2_t_C","sqrt_x^2+y^2_  versus t for 14C events", 600, 0,6.0); 
  neutrino_x2y2_t_14C->GetXaxis()->SetTitle("time (s)");
  neutrino_x2y2_t_14C->GetYaxis()->SetTitle("sqrt(x^2 + y^2) (m)");
  neutrino_x2y2_t_Bi = new TH1F("_x2y2_t_Bi"," sqrt_x^2+y^2_ versus t for BiPo events", 600, 0,6.0); 
  neutrino_x2y2_t_Bi->GetXaxis()->SetTitle("time (s)");
  neutrino_x2y2_t_Bi->GetYaxis()->SetTitle("sqrt(x^2 + y^2) (m)");
  neutrino_x_t_14C = new TH1F("_x_t_C"," x versus t for 14C events", 1000, -6.0, 6.0); 
  neutrino_x_t_14C->GetXaxis()->SetTitle("time (s)");
  neutrino_x_t_14C->GetYaxis()->SetTitle("x (m)");
  neutrino_x_t_Bi = new TH1F("_x_t_Bi"," x versus t for BiPo events", 1000, -6.0, 6.0); 
  neutrino_x_t_Bi->GetXaxis()->SetTitle("time (s)");
  neutrino_x_t_Bi->GetYaxis()->SetTitle("x (m)");
  neutrino_y_t_14C = new TH1F("_y_t_C"," y versus t for 14C events", 1000, -6.0, 6.0); 
  neutrino_y_t_14C->GetXaxis()->SetTitle("time (s)");
  neutrino_y_t_14C->GetYaxis()->SetTitle("y (m)");
  neutrino_y_t_Bi = new TH1F("_y_t_Bi"," y versus t for BiPo events", 1000, -6.0, 6.0); 
  neutrino_y_t_Bi->GetXaxis()->SetTitle("time (s)");
  neutrino_y_t_Bi->GetYaxis()->SetTitle("y (m)");

  neutrino_x_14C_c = new TH1F("_x_C_c","cumulative x  for 14C events", 1000, -6.0, 6.0);
  neutrino_x_14C_c->GetXaxis()->SetTitle("x(m)");
  neutrino_x_14C_c->GetYaxis()->SetTitle("counts");
  neutrino_y_14C_c = new TH1F("_y_C_c"," cumulative y  for 14C events", 1000, -6.0, 6.0);
  neutrino_y_14C_c->GetXaxis()->SetTitle("y(m)");
  neutrino_y_14C_c->GetYaxis()->SetTitle("counts");
  neutrino_z_14C_c = new TH1F("_z_C_c"," cumulative z  for 14C events", 1000, -6.0, 6.0);
  neutrino_z_14C_c->GetXaxis()->SetTitle("z(m)");
  neutrino_z_14C_c ->GetYaxis()->SetTitle("counts");

  neutrino_x_Bi_c = new TH1F("_x_Bi_c","cumulative x  for Bi events", 1000, -6.0, 6.0);
  neutrino_x_Bi_c->GetXaxis()->SetTitle("x(m)");
  neutrino_x_Bi_c->GetYaxis()->SetTitle("counts");
  neutrino_y_Bi_c = new TH1F("_y_Bi_c","cumulative y  for Bi events", 1000, -6.0, 6.0);
  neutrino_y_Bi_c->GetXaxis()->SetTitle("y(m)");
  neutrino_y_Bi_c->GetYaxis()->SetTitle("counts");
  neutrino_z_Bi_c = new TH1F("_z_Bi_c","cumulative z  for Bi events", 1000, -6.0, 6.0);
  neutrino_z_Bi_c->GetXaxis()->SetTitle("z(m)");
  neutrino_z_Bi_c->GetYaxis()->SetTitle("counts");

  neutrino_x_HE_c = new TH1F("_x_HE_c","cumulative x  for HE events", 1000, -6.0, 6.0);
  neutrino_x_HE_c->GetXaxis()->SetTitle("x(m)");
  neutrino_x_HE_c->GetYaxis()->SetTitle("counts");
  neutrino_y_HE_c = new TH1F("_y_HE_c","cumulative y  for HE events", 1000, -6.0, 6.0);
  neutrino_y_HE_c->GetXaxis()->SetTitle("y(m)");
  neutrino_y_HE_c->GetYaxis()->SetTitle("counts");
  neutrino_z_HE_c = new TH1F("_z_HE_c","cumulative z  for HE events", 1000, -6.0, 6.0);
  neutrino_z_HE_c->GetXaxis()->SetTitle("z(m)");
  neutrino_z_HE_c->GetYaxis()->SetTitle("counts");

  neutrino_x_VHE_c = new TH1F("_x_VHE_c","cumulative x  for VHE events", 1000, -6.0, 6.0);
  neutrino_x_VHE_c->GetXaxis()->SetTitle("x(m)");
  neutrino_x_VHE_c->GetYaxis()->SetTitle("counts");
  neutrino_y_VHE_c = new TH1F("_y_VHE_c","cumulative y  for VHE events", 1000, -6.0, 6.0);
  neutrino_y_VHE_c->GetXaxis()->SetTitle("y(m)");
  neutrino_y_VHE_c->GetYaxis()->SetTitle("counts");
  neutrino_z_VHE_c = new TH1F("_z_VHE_c","cumulative z  for VHE events", 1000, -6.0, 6.0);
  neutrino_z_VHE_c->GetXaxis()->SetTitle("z(m)");
  neutrino_z_VHE_c->GetYaxis()->SetTitle("counts");

  energy_14C = new TH1F("energy_14C", "charge of 14C  (p.e.)", 100,0,100);
  energy_14C->GetXaxis()->SetTitle("charge (p.e.)");
  energy_14C->GetYaxis()->SetTitle("counts");
  energy_Bi = new TH1F("energy_Bi", "charge of Bi events in BiPo chains (p.e.)", 1000,0,2000);
  energy_Bi->GetXaxis()->SetTitle("charge (p.e.)");
  energy_Bi->GetYaxis()->SetTitle("counts");
  energy_HE = new TH1F("energy_HE", "charge of high energy ( > 100 p. e.) events", 2500,0,5000);
  energy_HE->GetXaxis()->SetTitle("charge (p.e.)");
  energy_HE->GetYaxis()->SetTitle("counts");

  neutrino_E_spectrum_c = new TH1F("energy_spectrum_c", "charge spectrum (p.e.)", 1000,0,2000);
  neutrino_E_spectrum_c->GetXaxis()->SetTitle("charge (p.e.)");
  neutrino_E_spectrum_c->GetYaxis()->SetTitle("counts");

  hits_t_trigger_t = new TH1F("hits_t_trigger_t", "decoded hits time - trigger time", 13000,-25000,1000);
  hits_t_trigger_t->GetXaxis()->SetTitle("time diff");
  hits_t_trigger_t->GetYaxis()->SetTitle("counts");

  start_t_trigger_t = new TH1F("start_t_trigger_t", "cluster start time - trigger time", 13000,-25000,1000);
  start_t_trigger_t->GetXaxis()->SetTitle("time diff");
  start_t_trigger_t->GetYaxis()->SetTitle("counts");


  //create cumulative 2D histos
  neutrino_z_x2y2_14C_c  = new TH2F("_z_x2y2_C_c"," z versus sqrt_x^2+y^2_ for 14C events (cumulative)", 1000, 0,6.0,600, -6.0,6.0 );
  neutrino_z_x2y2_14C_c->GetXaxis()->SetTitle("sqrt(x^2 + y^2) (m)");
  neutrino_z_x2y2_14C_c->GetYaxis()->SetTitle("z (m)");
  neutrino_z_x2y2_Bi_c  = new TH2F("_z_x2y2_Bi_c"," z versus sqrt_x^2+y^2_ for BiPo events (cumulative)", 1000,0,6.0,600, -6.0,6.0 );
  neutrino_z_x2y2_Bi_c->GetXaxis()->SetTitle("sqrt(x^2 + y^2) (m)");
  neutrino_z_x2y2_Bi_c->GetYaxis()->SetTitle("z (m)");
  neutrino_z_x2y2_HE_c  = new TH2F("_z_x2y2_HE_c"," z versus sqrt_x^2+y^2_ for high energy ( > 100 p.e.) events (cumulative) ", 1000,0,6.0,600, -6.0,6.0 );
  neutrino_z_x2y2_HE_c->GetXaxis()->SetTitle("sqrt(x^2 + y^2) (m)");
  neutrino_z_x2y2_HE_c->GetYaxis()->SetTitle("z (m)");
  neutrino_z_x2y2_VHE_c  = new TH2F("_z_x2y2_VHE_c"," z versus sqrt_x^2+y^2_ for very high energy ( > 500 p.e.) events (cumulative) ", 1000,0,6.0,600, -6.0,6.0 );
  neutrino_z_x2y2_VHE_c->GetXaxis()->SetTitle("sqrt(x^2 + y^2) (m)");
  neutrino_z_x2y2_VHE_c->GetYaxis()->SetTitle("z (m)");
  neutrino_x_y_14C_c = new  TH2F("_x_y_C_c"," x versus y for 14C events (cumulative) ", 1000, -6.0,6.0,1000, -6.0,6.0 );
  neutrino_x_y_14C_c->GetXaxis()->SetTitle("x (m)");
  neutrino_x_y_14C_c->GetYaxis()->SetTitle("y (m)");
  neutrino_x_y_Bi_c = new  TH2F("_x_y_Bi_c"," x versus y for BiPo events (cumulative)", 1000, -6.0,6.0,1000, -6.0,6.0 );
  neutrino_x_y_Bi_c->GetXaxis()->SetTitle("x (m)");
  neutrino_x_y_Bi_c->GetYaxis()->SetTitle("y (m)");
  neutrino_x_y_HE_c = new  TH2F("_x_y_HE_c"," x versus y for high energy ( > 100 p.e.) events (cumulative)", 1000, -6.0,6.0,1000, -6.0,6.0 );
  neutrino_x_y_HE_c->GetXaxis()->SetTitle("x (m)");
  neutrino_x_y_HE_c->GetYaxis()->SetTitle("y (m)");
  neutrino_x_y_VHE_c = new  TH2F("_x_y_VHE_c"," x versus y for high energy ( > 500 p.e.) events (cumulative)", 1000, -6.0,6.0,1000, -6.0,6.0 );
  neutrino_x_y_VHE_c->GetXaxis()->SetTitle("x (m)");
  neutrino_x_y_VHE_c->GetYaxis()->SetTitle("y (m)");

  event_rate = new TH1F("ev_rate","echidna event rate before prescaling" ,1000,0,25000);
  event_rate ->  GetXaxis()->SetTitle("time (s)");
  event_rate ->  GetYaxis()->SetTitle("rate (Hz)");
  rate_buffer_up = new TH1F("rate_buffer_up","event rate in the buffer, z > 3" ,1000,0,25000);
  rate_buffer_up -> GetXaxis()->SetTitle( "time (s)");
  rate_buffer_up -> GetYaxis()->SetTitle("rate (Hz)");
  rate_buffer_equator = new TH1F("rate_buffer_equator","event rate in the buffer, |z| < 3" ,1000,0,25000);
  rate_buffer_equator -> GetXaxis()->SetTitle( "time (s)");
  rate_buffer_equator -> GetYaxis()->SetTitle("rate (Hz)");
  rate_buffer_down = new TH1F("rate_buffer_down","event rate in the buffer, z < -3" ,1000,0,25000);
  rate_buffer_down -> GetXaxis()->SetTitle( "time (s)");
  rate_buffer_down -> GetYaxis()->SetTitle("rate (Hz)");

  barn_interface::get ()->store (barn_interface::junk, neutrino_z_t_14C, this);
  barn_interface::get ()->store (barn_interface::junk, neutrino_z_t_Bi, this);
  barn_interface::get ()->store (barn_interface::junk,neutrino_x2y2_t_14C, this);
  barn_interface::get ()->store (barn_interface::junk,neutrino_x2y2_t_Bi, this);
  barn_interface::get ()->store (barn_interface::junk,neutrino_x_t_14C, this);
  barn_interface::get ()->store (barn_interface::junk,neutrino_x_t_Bi, this);
  barn_interface::get ()->store (barn_interface::junk,neutrino_y_t_14C, this);
  barn_interface::get ()->store (barn_interface::junk,neutrino_y_t_Bi, this);
  barn_interface::get ()->store (barn_interface::junk, neutrino_z_x2y2_14C_c, this);
  barn_interface::get ()->store (barn_interface::junk, neutrino_z_x2y2_Bi_c, this);
  barn_interface::get ()->store (barn_interface::junk, neutrino_z_x2y2_HE_c, this);
  barn_interface::get ()->store (barn_interface::junk, neutrino_z_x2y2_VHE_c, this);
  barn_interface::get ()->store (barn_interface::junk,neutrino_x_y_14C_c, this);
  barn_interface::get ()->store (barn_interface::junk,neutrino_x_y_Bi_c, this);
  barn_interface::get ()->store (barn_interface::junk,neutrino_x_y_HE_c, this);
  barn_interface::get ()->store (barn_interface::junk,neutrino_x_y_VHE_c, this);
  barn_interface::get ()->store (barn_interface::junk, energy_14C, this);
  barn_interface::get ()->store (barn_interface::junk, energy_Bi, this);
  barn_interface::get ()->store (barn_interface::junk, energy_HE, this);
  barn_interface::get ()->store (barn_interface::junk, neutrino_x_14C_c, this);
  barn_interface::get ()->store (barn_interface::junk, neutrino_y_14C_c, this);
  barn_interface::get ()->store (barn_interface::junk, neutrino_z_14C_c, this);
  barn_interface::get ()->store (barn_interface::junk, neutrino_x_Bi_c, this);
  barn_interface::get ()->store (barn_interface::junk, neutrino_y_Bi_c, this);
  barn_interface::get ()->store (barn_interface::junk, neutrino_z_Bi_c, this);
  barn_interface::get ()->store (barn_interface::junk, neutrino_x_HE_c, this);
  barn_interface::get ()->store (barn_interface::junk, neutrino_y_HE_c, this);
  barn_interface::get ()->store (barn_interface::junk, neutrino_z_HE_c, this);
  barn_interface::get ()->store (barn_interface::junk, neutrino_x_VHE_c, this);
  barn_interface::get ()->store (barn_interface::junk, neutrino_y_VHE_c, this);
  barn_interface::get ()->store (barn_interface::junk, neutrino_z_VHE_c, this);
  barn_interface::get ()->store (barn_interface::junk, neutrino_E_spectrum_c, this);
  barn_interface::get ()->store (barn_interface::junk, event_rate, this);
  barn_interface::get ()->store (barn_interface::junk, hits_t_trigger_t, this);
  barn_interface::get ()->store (barn_interface::junk, start_t_trigger_t, this);
  barn_interface::get ()->store (barn_interface::junk, rate_buffer_up, this);
  barn_interface::get ()->store (barn_interface::junk, rate_buffer_equator, this);
  barn_interface::get ()->store (barn_interface::junk, rate_buffer_down, this);

  f_run_histos_list = new TList();

  f_run_histos_list->Add(energy_14C);
  f_run_histos_list->Add(energy_Bi);
  f_run_histos_list->Add(energy_HE);
  f_run_histos_list->Add(neutrino_z_x2y2_14C_c);
  f_run_histos_list->Add(neutrino_z_x2y2_Bi_c);
  f_run_histos_list->Add(neutrino_z_x2y2_HE_c);
  f_run_histos_list->Add(neutrino_z_x2y2_VHE_c);
  f_run_histos_list->Add(neutrino_x_y_14C_c);
  f_run_histos_list->Add(neutrino_x_y_Bi_c);
  f_run_histos_list->Add(neutrino_x_y_HE_c);
  f_run_histos_list->Add(neutrino_x_y_VHE_c);


  f_run_histos_list->Add(neutrino_x_14C_c);
  f_run_histos_list->Add(neutrino_y_14C_c);
  f_run_histos_list->Add(neutrino_z_14C_c);
  f_run_histos_list->Add(neutrino_x_Bi_c);
  f_run_histos_list->Add(neutrino_y_Bi_c);
  f_run_histos_list->Add(neutrino_z_Bi_c);
  f_run_histos_list->Add(neutrino_x_HE_c);
  f_run_histos_list->Add(neutrino_y_HE_c);
  f_run_histos_list->Add(neutrino_z_HE_c);
  f_run_histos_list->Add(neutrino_x_VHE_c);
  f_run_histos_list->Add(neutrino_y_VHE_c);
  f_run_histos_list->Add(neutrino_z_VHE_c);
  f_run_histos_list->Add(neutrino_E_spectrum_c);
  f_run_histos_list->Add(hits_t_trigger_t);
  f_run_histos_list->Add(start_t_trigger_t);
  f_run_histos_list->Add(rate_buffer_up);
  f_run_histos_list->Add(rate_buffer_equator);
  f_run_histos_list->Add(rate_buffer_down);

  f_run_histos_list->Add(event_rate);

  f_run_histos_list->SetOwner(kTRUE);

  f_temp_histos_list = new TList();

  f_temp_histos_list->SetOwner(kTRUE);

  f_prev_bunch_time = 0;
  f_prev_bi_event_time = 0;
  f_Bi_E = 0;
  f_14C_E = 0;
  f_HE_E = 0;
  f_start_time = 0;
  f_n_events = 0;
  f_n_events_buffer_up = 0;
  f_n_events_buffer_equator = 0;
  f_n_events_buffer_down = 0;
  f_total_live_time = 0;
  f_begin_new_bunch = 0;
  f_end_new_bunch = 0;

}


bx_echidna_event* bx_calib_monitor::doit (bx_echidna_event *ev) {

  f_n_events++;

  bool is_bipo = false;
  bool is_bi = false;
  bool is_carbon = false;
  bool is_HE = false;
  bool is_VHE = false;
  bool fill_event_rate = false;

  const bx_laben_event& er = ev->get_laben();
  const bx_trigger_event& et = ev->get_trigger();
  bool trg_type = et.get_trgtype();
  int btb_inputs = et.get_btb_inputs();
  int n_clusters = er.get_nclusters();
  double n_decoded = er.get_decoded_nhits();
  double trigger_time = er.get_trigger_rawt();
  double mean_time = 0;
  double start_time = 0;
  if(n_clusters > 0) {
    mean_time = er.get_cluster(0).get_mean_time();
    start_time = er.get_cluster(0).get_start_time();
  }
  if(trg_type == 1 && btb_inputs == 0 && n_clusters == 1 ){
    const bx_base_position& pos = er.get_rec_cluster(0).get_position();
    double charge = er.get_cluster(0).get_charge();
    neutrino_E_spectrum_c->Fill(charge);
    start_t_trigger_t->Fill(start_time -trigger_time);


    double r = pos.get_r();
    double z = pos.get_z();
    double x = pos.get_x();
    double y = pos.get_y();
    double x_y = ::sqrtf(x*x + y*y); 

    unsigned long event_time_nsec;
    unsigned long event_time_sec;
    et.get_gps_time ( event_time_sec, event_time_nsec);
    double time_this_event =  event_time_sec + (double) event_time_nsec * 1E-9 ; 

    if(f_start_time == 0){
      f_start_time = time_this_event;
      f_begin_new_bunch = time_this_event;
      f_end_new_bunch = time_this_event;
    }
    if(time_this_event - f_start_time < 200E-3){
      f_end_new_bunch = time_this_event;
      fill_event_rate = false;
    }
    else{
      fill_event_rate = true;
      f_total_live_time += f_end_new_bunch - f_begin_new_bunch;
      f_begin_new_bunch = time_this_event;
      f_end_new_bunch = time_this_event;
    }
    if(charge > 200 && charge < 500) {   //Po candidate
      if(f_prev_bi_event_time > 0){
	double delta_time = time_this_event - f_prev_bi_event_time;
	if(delta_time > 13E-6 && delta_time < 500E-6){
	  is_bipo = true;
	}
      } 
    }
    
    if (charge < 90) {
      is_carbon = true;
      f_14C_E = charge;
    }
    else if (charge > 100){
      is_HE = true;
      f_HE_E = charge;
    }
    if(charge > 500){
      is_VHE = true;

    }
    if(charge > 90 && charge < 1600) {  //Bi candidate
      f_prev_bi_event_time = time_this_event; 
      is_bi = true;
      f_Bi_E = charge;
    }
    unsigned long elapsed_time_sec = event_time_sec - f_prev_bunch_time;
    if(elapsed_time_sec > f_buffer_time){
      if(f_prev_bunch_time > 0){
	m_send_and_reset_histos();
      }
      f_prev_bunch_time = event_time_sec;
    }

    //Fill histograms
    if(is_carbon && mean_time < 1000){
      neutrino_z_x2y2_14C_c->Fill(x_y,z);
      neutrino_x_y_14C_c->Fill(x,y);
      energy_14C->Fill(f_14C_E);
    }
    if(is_bipo && mean_time < 1000){
      neutrino_z_x2y2_Bi_c->Fill(f_old_x_y, f_old_z);
      neutrino_x_y_Bi_c->Fill(f_old_x,f_old_y);
      energy_Bi->Fill(f_Bi_E);
    }
    if(is_HE && mean_time < 1000){
      neutrino_z_x2y2_HE_c->Fill(x_y, z);
      neutrino_x_y_HE_c->Fill(x,y);
      energy_HE->Fill(f_HE_E);
    }
    if(is_VHE &&  mean_time > 1000){
      neutrino_z_x2y2_VHE_c->Fill(x_y, z);
      neutrino_x_y_VHE_c->Fill(x,y);
    }

    if(is_bi == true){
      f_old_x = x;
      f_old_y = y;
      f_old_z = z;
      f_old_r = r;
      f_old_x_y = x_y;
    }
    if( r > 4.25){  //event reconstructer in the buffer
      if( z < -3 ) f_n_events_buffer_down++;
      else if( z > 3 ) f_n_events_buffer_up++;
      else  f_n_events_buffer_equator++;
    }
    if(fill_event_rate == true){
      event_rate->Fill(time_this_event - f_start_time,(double)f_n_events/f_total_live_time);
      rate_buffer_up->Fill(time_this_event - f_start_time,(double)f_n_events_buffer_up/f_total_live_time);
      rate_buffer_equator->Fill(time_this_event - f_start_time,(double)f_n_events_buffer_equator/f_total_live_time);
      rate_buffer_down->Fill(time_this_event - f_start_time,(double)f_n_events_buffer_down/f_total_live_time);
    }

  }
  else if (trg_type == 1 && btb_inputs == 0 ){
    for(int i = 0; i < n_decoded; i++){
      double hits_time = er.get_decoded_hit(i).get_raw_time();
      hits_t_trigger_t->Fill(hits_time - trigger_time);
    }
  }

  return ev;
}

void bx_calib_monitor::end () {
  get_message (bx_message::debug) << "bx_calib_monitor end" << dispatch;

}

void  bx_calib_monitor::m_send_and_reset_histos(){
  TIter iterator(f_temp_histos_list);
  TH2F* histo;
  while ((histo = (TH2F* )iterator.Next())){
    double n_entries = histo->GetEntries();
    if(n_entries < 500) histo->SetMarkerStyle(2);
    barn_interface::get ()->network_send (histo, this);
    histo->Reset();
  } 
  TIter iterator_run(f_run_histos_list);
  while ((histo = (TH2F* )iterator_run.Next())){
    double  n_entries = histo->GetEntries();
    if(n_entries < 500) histo->SetMarkerStyle(2);
    barn_interface::get ()->network_send (histo, this);
  }
}

void bx_calib_monitor::m_delete_histos(){
  delete f_run_histos_list;
  delete f_temp_histos_list;
}

double bx_calib_monitor::m_fit2D(double &x, double &y, TH2F * histo){
  return 0;

}
/*
 * $Log: bx_calib_monitor.cc,v $
 * Revision 1.24  2011/03/02 18:51:57  guardi
 * dependence on the reconstruction algorithm removed
 *
 * Revision 1.23  2011-03-02 14:13:44  guardi
 * milano position reconstruction replaced with lngs one
 *
 * Revision 1.22  2009-10-23 09:22:15  guardi
 * debug messages removed
 *
 * Revision 1.21  2009-03-31 15:21:14  guardi
 * name of an histogram corrected
 *
 * Revision 1.20  2009-03-30 08:49:51  guardi
 * plots for event rate in the buffer (upper part, equator, bottom ) added and more non cumulative plots removed
 *
 * Revision 1.19  2009-03-30 08:17:37  guardi
 * non cumulative histograms removed
 *
 * Revision 1.18  2009-03-23 10:26:19  guardi
 *  histograms for hits_time - trigger_time and start_time - trigger_time added
 *
 * Revision 1.17  2008-10-07 17:57:53  guardi
 * fix
 *
 * Revision 1.16  2008-10-07 17:40:30  guardi
 * added event rate histogram
 *
 * Revision 1.15  2008-10-07 16:20:54  guardi
 * debug
 *
 * Revision 1.14  2008-10-07 13:48:16  guardi
 * cut on mean time different for VHE plots
 *
 * Revision 1.13  2008-10-07 13:14:28  guardi
 * added check on mean time
 *
 * Revision 1.12  2008-10-06 15:49:21  guardi
 * debug
 *
 * Revision 1.11  2008-10-04 11:07:12  guardi
 * very high energy events threshold lowered to 500 p.e.
 *
 * Revision 1.10  2008-10-04 10:44:45  guardi
 * histograms renamed and coincidence time window for BiPo search reduced
 *
 * Revision 1.9  2008-10-04 10:27:10  guardi
 * added full spectrum and plots for the position of >1000 p.e. events
 *
 * Revision 1.8  2008-10-04 08:58:09  guardi
 *  FSR on axis set to 6 m
 *
 * Revision 1.7  2008-10-04 08:50:55  guardi
 *  debug
 *
 * Revision 1.6  2008-10-04 08:20:04  guardi
 * added cumulative 1D plots and debug
 *
 * Revision 1.5  2008-10-04 07:33:45  guardi
 * added 1D position distributions
 *
 * Revision 1.4  2008-10-03 08:18:51  razeto
 * Upgraded (by elena)
 *
 * Revision 1.3  2008-10-02 16:21:23  razeto
 * Require only position_mi
 *
 * Revision 1.2  2008-10-01 16:53:28  guardi
 * charge spectra and axis titles added
 *
 * Revision 1.1  2008-10-01 15:55:42  guardi
 * module to plot the position of the events online added. Purpose: determine the position of the source during calibration
 *
 *
 */
