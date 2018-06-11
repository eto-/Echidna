/* BOREXINO Reconstruction program
 *
 * Author: Alessandro Fontana <alessandro.fontana@mi.infn.it>
 * Maintainer: Alessandro Fontana <alessandro.fontana@mi.infn.it>
 *
 *module to get some information about radial and oblique laser runs
 */
#include "bx_calib_laser_transparency.hh"
#include "bx_echidna_event.hh"
#include "barn_interface.hh"
#include "bx_dbi.hh"
#include "db_calib.hh"
#include "db_profile.hh"
#include "db_channel.hh"
#include "bx_trigger_event.hh"
#include "constants.hh"
#include <math.h>
#include <sstream>
#include <iostream>
#define speed_light 0.3

// ctor
bx_calib_laser_transparency::bx_calib_laser_transparency():bx_base_module("bx_calib_laser_transparency", bx_base_module::main_loop) {
  x_bar = 0;
  y_bar = 0;
  N_bar = 0;
  Nbins = 200;
  distance = 0;
  dir_peak_tot = 0;
  abs_range = 10;
  A = Nbins/2/abs_range;  //meters-bins conversion (range is [-abs_range, abs_range] meters) 
  B = Nbins/2;
  counter_ev = 0;
  processed_events = 0;
  R = 6.33;
  zero_charge_dir_peak = 0;
  require_event_stage(bx_detector::laben, bx_base_event::decoded);
  require_trigger_type(bx_trigger_event::laser394);  
}

const float bx_calib_laser_transparency::buffer_source_theta[19] = {26.4,26.4,26.4,26.4,26.4,26.4,7.,7.,7.,7.,7.,-26.4,-26.4,-26.4,-26.4,-26.4,-26.4,-7.,-7.};
const float bx_calib_laser_transparency::buffer_source_phi[19] = {335.,50.,5.,125.,215.,260.,10.,90.,170.,190.,270.,185.,155.,140.,35.,305.,320.,270.,90.}; 
const int32_t bx_calib_laser_transparency::buffer_target[19] = {942,1248,314,448,1427,475,1934,1604,147,319,1779,922,1113,1606,1615,1932,1777,1779,1116};
const float bx_calib_laser_transparency::radial_source_theta[12] = {67.15,50.7,26.4,26.4,26.4,7.,-26.4,-26.4,-7.,-50.7,-67.15,-26.4};
const float bx_calib_laser_transparency::radial_source_phi[12] = {190.,250.,65.,185.,305.,210.,125.,245.,30.,70.,10.,5.};
// module interface
void bx_calib_laser_transparency::begin() {
  get_message(bx_message::debug) << " Begin " << dispatch;

  time_vs_angle = new TH2F("time_vs_angle","time vs angle",200, 0, 200,180 ,0 , 180);
  time_dis = new TH1F("time_dis","time distribution",200,0.,200.); 
  time_dis_calib = new TH1F("time_dis_calib","time distribution for calibration",200,0.,200.);
  frontpro = new TH2F("frontpro", "frontpro",Nbins,-abs_range ,abs_range ,Nbins , -abs_range, abs_range);
  backpro = new TH2F("backpro", "backpro",Nbins,-abs_range ,abs_range ,Nbins , -abs_range, abs_range);
  top = new TH2F("top", "top",Nbins,-abs_range ,abs_range ,Nbins , -abs_range, abs_range);
  bottom = new TH2F("bottom", "bottom",Nbins,-abs_range ,abs_range ,Nbins , -abs_range, abs_range);
  charge_dis = new TH1F("charge_dis","charge for direct target",150,0.,15);
  raw_charge_dis = new TH1F("raw_charge_dis","raw charge for direct target",250,0.,250);
  lg_dis = new TH1F ("lg_dis","lg_dis",2240,1,2241);
  
  barn_interface::get()->store(barn_interface::file,time_dis_calib,this); 
  barn_interface::get()->store(barn_interface::file,time_dis,this);
  barn_interface::get()->store(barn_interface::file,time_vs_angle,this);
  barn_interface::get()->store(barn_interface::file,frontpro,this);
  barn_interface::get()->store(barn_interface::file,backpro,this);
  barn_interface::get()->store(barn_interface::file,top,this);
  barn_interface::get()->store(barn_interface::file,bottom,this);
  barn_interface::get()->store(barn_interface::file,charge_dis,this);
  barn_interface::get()->store(barn_interface::file,raw_charge_dis,this);
  barn_interface::get()->store(barn_interface::file,lg_dis,this);
  if (get_parameter ("refidx").get_float () > 0) refidx = get_parameter ("refidx").get_float ();
  else refidx = bx_dbi::get()->get_calib().get_refraction_index_data();
  if (get_parameter ("time_calib_events").get_int () > 0) time_calib_events = get_parameter ("time_calib_events").get_int ();
  else time_calib_events = 1000;
  if (get_parameter ("threshold").get_float () > 0) threshold = get_parameter ("threshold").get_float ();
  else threshold = 0.02;
  if (get_parameter ("peak_width").get_float () > 0) peak_width = get_parameter ("peak_width").get_float ();
  else peak_width = 3;
  if (get_parameter ("resolution").get_int () > 0) resolution = get_parameter ("resolution").get_int ();
  else resolution = 3;
  type = get_parameter ("type").get_int();
  patch_panel = get_parameter ("patch_panel").get_int();
  delay = get_parameter ("delay").get_float();
  get_message (bx_message::log) << " Refraction Index = " << refidx <<dispatch;
  
  if (type == 1) {
    get_message (bx_message::log) << " Laser type = Radial" <<dispatch;
    theta_src = radial_source_theta[patch_panel - 1]*pi/180;
    phi_src = radial_source_phi[patch_panel - 1]*pi/180;
    x_s = cos(phi_src)*cos(theta_src);
    y_s = sin(phi_src)*cos(theta_src);
    z_s = sin(theta_src);
    STe = sin( -theta_src + pi/2);
    CTe = cos( -theta_src + pi/2);
    SPe = sin(phi_src + pi/2);
    CPe = cos(phi_src + pi/2);
  }
   if (type == 2) {
    get_message (bx_message::log) << " Laser type = Oblique" <<dispatch;
    theta_src = buffer_source_theta[patch_panel - 1]*pi/180;
    phi_src = buffer_source_phi[patch_panel - 1]*pi/180;
    nominal_target = buffer_target[patch_panel -1];
  }
  
  get_message (bx_message::log) << " Patch panel = " << patch_panel <<dispatch; 
  
}

bx_echidna_event* bx_calib_laser_transparency::doit (bx_echidna_event *ev) {
  const bx_laben_event& er = ev->get_laben();
  counter_ev++;
  int32_t nhits = er.get_decoded_nhits();
  
  //Filling calibration histo
  if (counter_ev <= time_calib_events){
    for(int32_t i = 0; i < nhits; i++){
      const bx_laben_decoded_hit& hit = er.get_decoded_hit(i);
      double itime = hit.get_raw_time() - er.get_laser_rawt() - delay;
      int32_t npe = hit.get_charge_npe();
      time_dis_calib->SetBinContent(int32_t(itime), time_dis_calib->GetBinContent(int32_t(itime)) + npe);
    }
  }
  //Searching peaks
  if (counter_ev == time_calib_events){
    dir_peak_time = time_dis_calib->GetMaximumBin();
    time_inf = dir_peak_time - 5.;
    time_sup = dir_peak_time + 5.;
    TSpectrum * time_calib = new TSpectrum(10); 
    time_calib->Search(time_dis_calib ,peak_width ,"new", threshold);
    peak_time = time_calib->GetPositionX();
    get_message(bx_message::log) << " Number of peaks in time distribution: " << time_calib->GetNPeaks() <<dispatch;
    for(int32_t k = 0; k < time_calib->GetNPeaks(); k++){
      get_message(bx_message::log) << " Time position for peak number " << k + 1 << " " << peak_time[k] << " in meters " << peak_time[k]*0.3/refidx <<dispatch;
    }
    get_message(bx_message::log) << " Direct peak flight time: " << dir_peak_time << " in meters " << (dir_peak_time)*0.3/refidx <<dispatch;
    get_message(bx_message::log) << " Time window direct peak: " << time_inf << " - " << time_sup << " ns " <<dispatch;
  }
  
  //Search for target PMT in direct peak time gate
  if (counter_ev > time_calib_events && counter_ev < time_calib_events + 1000){
   for(int32_t i = 0; i < nhits; i++){
      const bx_laben_decoded_hit& hit = er.get_decoded_hit(i);
      const db_channel* ch_info = hit.get_db_channel();
      ich = ch_info->get_lg();
      double itime = hit.get_raw_time();
      itime -= er.get_laser_rawt() + delay;
      if(itime>time_inf && itime < time_sup){
        int32_t npe = hit.get_charge_npe();
        ich = ch_info->get_lg();
        lg_dis->SetBinContent(ich, lg_dis->GetBinContent(ich) + npe);
      }
    }
  
  }
  
  
  if (counter_ev == time_calib_events + 1000){
    target = lg_dis->GetMaximumBin();
    get_message(bx_message::log) << " Logical channel for target = " << target <<dispatch;
    if(type == 2 && nominal_target != target) get_message(bx_message::warn) <<" Mismatch between nominal and effective target " << dispatch;
  }
  
  
  if(counter_ev > time_calib_events + 1000){
    processed_events++;
    counter_charge_dir_peak = 0;
    dir_peak_ev = 0;
    for(int32_t j = 0; j<nhits; j++){
      fbflag = 0;
      const bx_laben_decoded_hit& hit = er.get_decoded_hit(j);
      const db_channel* ch_info = hit.get_db_channel();
      ich = ch_info->get_lg();
      double itime = hit.get_raw_time();
      itime -= er.get_laser_rawt() + delay;
      int32_t npe = hit.get_charge_npe();
      raw_charge = hit.get_charge_bin();
      time_dis->SetBinContent(int32_t(itime), time_dis->GetBinContent(int32_t(itime)) + npe);
      if(itime > time_inf && itime < time_sup){
      dir_peak_ev++;
        lg_dis->Fill(ich);
        counter_charge_dir_peak +=npe;
        if(ich==target){    
          charge_dis->Fill(npe);
          raw_charge_dis->Fill(raw_charge);
        }
      }
   
      const db_channel_laben& channel_info = dynamic_cast<const db_channel_laben&>(bx_dbi::get()->get_channel(ich));
      //posizione pmt colpito
      x_pmt = channel_info.pmt_x(); 
      y_pmt = channel_info.pmt_y();
      z_pmt = channel_info.pmt_z();
      omega1 = acos(fabs(z_pmt)/R);
      mod1 = R*omega1; 
      norma_tb = ::sqrtf(x_pmt*x_pmt + y_pmt*y_pmt);
      x_tb = mod1*x_pmt/norma_tb;
      y_tb = mod1*y_pmt/norma_tb;  
      norma_pmt = ::sqrtf(x_pmt*x_pmt + y_pmt*y_pmt + z_pmt*z_pmt);
      if(type == 1 /*radial*/){
        //rotazione (asse z lungo la congiungente sorgente centro)
        x_r = x_pmt*CPe + SPe*y_pmt;
        y_r = +x_pmt*CTe*SPe - y_pmt*CTe*CPe - z_pmt*STe;
        z_r = -x_pmt*STe*SPe + y_pmt*STe*CPe - z_pmt*CTe;
        x_pmt /= norma_pmt;
        y_pmt /= norma_pmt;
        z_pmt /= norma_pmt;
        //angolo al centro tra sorgente e pmt in esame
        alpha = acos((x_s*x_pmt)+(y_s*y_pmt)+(z_s*z_pmt))*180/pi;
        time_vs_angle->SetBinContent(int32_t(itime),int32_t(alpha),time_vs_angle->GetBinContent(int32_t(itime),int32_t(alpha))+npe);
        if(z_pmt > 0.) top->SetBinContent(int32_t(x_tb*A + Nbins/2), int32_t(y_tb*A + B), top->GetBinContent(int32_t(x_tb*A + B),int32_t(y_tb*A + B))+1);
        else bottom->SetBinContent(int32_t(x_tb*A + Nbins/2), int32_t(y_tb*A + B) , bottom->GetBinContent(int32_t(x_tb*A + B),int32_t(y_tb*A + B))+1);
        //calcolo per le variabili proiettate (srotolamento delle due semisfere)
        omega2 = acos(z_r/R);
        mod2 = R*omega2;
        norma_r = sqrt(x_r*x_r + y_r*y_r);
        x_p = mod2*x_r/norma_r;
        y_p = mod2*y_r/norma_r;
        if(omega2>pi/2.) {omega2 = pi - omega2; fbflag = 1;}
        mod2 = R*omega2;
        norma_r = sqrt(x_r*x_r + y_r*y_r);
        if(norma_r != 0.){
          x_p = mod2*x_r/norma_r;
          y_p = mod2*y_r/norma_r;
        }
        //riempimento scatterplot illustrativi
        if(itime > time_inf && itime < time_sup){
        if( fbflag == 0){
          frontpro->SetBinContent(int32_t(x_p*A + B), int32_t(y_p*A + B) , frontpro->GetBinContent(int32_t(x_p*A + B),int32_t(y_p*A + B)) + npe);
          x_bar += x_p*npe;
          y_bar += y_p*npe;
          N_bar +=npe;}
        else
        backpro->SetBinContent(int32_t(x_p*A + B), int32_t(y_p*A + B) , backpro->GetBinContent(int32_t(x_p*A + B),int32_t(y_p*A + B)) + npe);
        }
      }
     
   }
   if(counter_charge_dir_peak == 0) zero_charge_dir_peak++;

  }
  return 0;
}






void bx_calib_laser_transparency::end () {
  get_message(bx_message::debug) << " End " <<dispatch;
  if (type == 1){
    x_bar /= N_bar;
    y_bar /= N_bar;
    get_message(bx_message::log) << " This is a Radial Run " <<dispatch;
    get_message(bx_message::log) << " Total charge collected (in NPE excluded calibration event) = " << time_dis->Integral() <<dispatch;
    get_message(bx_message::log) << " center of light position  X =  " << x_bar << " Y = " << y_bar << dispatch;
    for(int32_t i=1; i <= Nbins; i++){
      for (int32_t j=1; j <= Nbins; j++){
        x_p = (i - B)/A;
	y_p = (j - B)/A;
	int32_t npe = int32_t(frontpro->GetBinContent(i,j));
        distance += sqrt(pow((x_p - x_bar),2) + pow((y_p - y_bar),2))*npe;
      }
    }
    distance /= N_bar;
    get_message(bx_message::log) << " Radius of spot of light in front of source = " << distance <<dispatch;
    get_message(bx_message::log) << " Back charge in semisphere containing source = " << time_vs_angle->Integral(1,160,1,90)/time_dis->Integral()*100 << " % " <<dispatch;
    get_message(bx_message::log) << " Charge in spot R=1 m = " << frontpro->Integral( int32_t((-1+x_bar)*A+B) ,int32_t((1+x_bar)*A+B), int32_t((-1+y_bar)*A+B),int32_t((1+y_bar)*A+B))/time_dis->Integral()*100 << " % " <<dispatch;
    get_message(bx_message::log) << " Ratio Spot R=1 charge/back charge " <<  frontpro->Integral( int32_t((-1+x_bar)*A+B) ,int32_t((1+x_bar)*A+B),int32_t((-1+y_bar)*A+B),int32_t((1+y_bar)*A+B))/time_vs_angle->Integral(1,160,1,90) <<dispatch;
    get_message(bx_message::log) << " Zero charge in direct peak (% on total event)= " << zero_charge_dir_peak/processed_events*100 <<dispatch;
  }
  if(type == 2){
    double saturation = (raw_charge_dis->Integral(210,250)/raw_charge_dis->Integral())*100;
    get_message(bx_message::log) << " This is an oblique Run " <<dispatch;
    get_message(bx_message::log) << " Fraction of saturated events = " << saturation <<dispatch;
    get_message(bx_message::log) << " Mean charge on target PMT = " << charge_dis->GetMean() <<dispatch;
    get_message(bx_message::log) << " Zero charge events on target PMT (%) = " << ( - charge_dis->GetEntries() + processed_events)/(processed_events)*100 <<dispatch;
    get_message(bx_message::log) << " Mean zero charge events included = " << charge_dis->GetMean()*(1 + (charge_dis->GetEntries() - processed_events)/processed_events) <<dispatch;
  } 
}


/*
 *
 */
