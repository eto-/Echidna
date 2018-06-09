/* BOREXINO Reconstruction program
 *
 * Author: Barbara Caccianiga <Barbara.Caccianiga@mi.infn.it>
 * Maintainer: Barbara Caccianiga <Barbara.Caccianiga@mi.infn.it>
 *
 * $Id: bx_calib_fiber_bundle.cc,v 1.8 2007/02/06 11:24:50 bcaccian Exp $
 *
 * Implementation of bx_calib_fiber_bundle
 *
 */
#include "bx_calib_fiber_bundle.hh"
#include "bx_echidna_event.hh"
#include "barn_interface.hh"
#include "bx_dbi.hh"
#include "db_profile.hh"
#include "db_channel.hh"
#include "bx_trigger_event.hh"
#include "TH1F.h"
#include "TH2F.h"
#include "constants.hh"
#include <math.h>
#include <sstream>
#include <iostream>

// ctor
bx_calib_fiber_bundle::bx_calib_fiber_bundle():bx_base_module("bx_calib_fiber_bundle", bx_base_module::main_loop) {
  i4_times = 0;
  require_event_stage(bx_detector::laben, bx_base_event::decoded);  
  require_trigger_type(bx_trigger_event::laser394);
}

// module interface
void bx_calib_fiber_bundle::begin() {
  get_message(bx_message::debug) << "bx_calib_fiber_bundle " << dispatch;
  n_bundles = constants::fiber::number_of_bundles;
  f_laser_low = get_parameter("laser_low").get_float();
  f_laser_high = get_parameter("laser_high").get_float();
  f_dark_low = get_parameter("dark_low").get_float();
  f_dark_high = get_parameter("dark_high").get_float();
  get_message(bx_message::info) << "number of fiber bundles " << n_bundles <<dispatch;
  get_message(bx_message::info) << "laser events between  "<< f_laser_low << "  and  " << f_laser_high <<  dispatch;
  get_message(bx_message::info) << "dark noise events between  "<< f_dark_low << "  and  " << f_dark_high <<  dispatch;
  p_timetot = new TH1F("time", "time" , 20000,-10000.,10000.);
  p_laser = new TH2F ("laser", "channel vs bundle for laser event", n_bundles,1.,(n_bundles+1),2240,1.,2241.);
  p_dark_noise = new TH2F ("darknoise", "channel vs bundle for dark noise", n_bundles,1.,(n_bundles+1),2240,1.,2241.);
  barn_interface::get()->store(barn_interface::file,p_timetot,this);
  barn_interface::get()->store(barn_interface::file,p_laser,this);
  barn_interface::get()->store(barn_interface::file,p_dark_noise,this);
}

bx_echidna_event* bx_calib_fiber_bundle::doit (bx_echidna_event *ev) {
  const bx_laben_event& er = ev->get_laben();
  i4_times = i4_times+1;  
  int nhits = er.get_decoded_nhits();
  for(int i = 0; i < nhits; i++){
    const bx_laben_decoded_hit& hit = er.get_decoded_hit(i);
    const bx_laben_raw_hit& rawhit = hit.get_raw_hit();
    unsigned short ch = rawhit.get_logical_channel();
    int ibundle = hit.get_db_channel()->pmt_fiber_bundle();
    double t = hit.get_raw_time();
    t = t-er.get_laser_rawt();
    p_timetot->Fill(t);
    if(t > f_laser_low && t < f_laser_high){
      p_laser->Fill(ibundle,ch);
    }
    if(t > f_dark_low && t < f_dark_high){
      p_dark_noise->Fill(ibundle,ch);
    }
  }
  return ev;
}

void bx_calib_fiber_bundle::end () {

  TH1D* bundles_laser = p_laser->ProjectionX("bundles_laser");
  TH1D* channels_laser = p_laser->ProjectionY("channels_laser");
  TH1D* bundles_dark = p_dark_noise->ProjectionX("bundles_dark");
  TH1D* channels_dark = p_dark_noise->ProjectionY("channels_dark");
  TH1D* bundle_v[n_bundles+1];
  TH1D* bundle_dark_v[n_bundles+1];
  for(int i = 1; i < n_bundles+1; i++){
    const int n_fibers = constants::fiber::number_of_fibers_in_bundle[i];
    std::ostringstream name;
    std::ostringstream name_dark;
    name <<"bundle_"<<i;
    name_dark <<"bundle_dark"<<i;
    TH1D* bundle_i = p_laser->ProjectionY(name.str().c_str(), i, i);
    TH1D* bundle_dark_i = p_dark_noise->ProjectionY(name_dark.str().c_str(), i, i);
    bundle_v[i] = bundle_i;
    bundle_dark_v[i] = bundle_dark_i;
    barn_interface::get()->store(barn_interface::file,bundle_i,this);
    barn_interface::get()->store(barn_interface::file,bundle_dark_i,this);
    float Integral = bundle_i->Integral();
    float Integral_dark = bundle_dark_i->Integral();
    int n_broken = 0;
    int n_broken_dark = 0;
    float rms = 0.;
    float rms_dark = 0.;
    float Mean = Integral / n_fibers;
    float Mean_dark = Integral_dark / n_fibers;
    int fiber[81] = {0,};
    int pmt[81] = {0,};
    for(int j = 1; j < 2241; j++){
      double binc = p_laser->GetBinContent(i,j);
      double binc_dark = p_dark_noise->GetBinContent(i,j);
      int ib = dynamic_cast<const db_channel_laben&>(bx_dbi::get()->get_channel(j)).pmt_fiber_bundle();
      if(ib == i){
        if( binc != 0.) rms = rms + (binc-Mean)*(binc-Mean);
	else if( binc == 0.){
          n_broken = n_broken + 1;
	  fiber[n_broken] = j;
        }
        if( binc_dark != 0.) rms_dark = rms_dark + (binc_dark-Mean_dark)*(binc_dark-Mean_dark);
	else if( binc_dark == 0.){
          n_broken_dark = n_broken_dark + 1;
	  pmt[n_broken_dark] = j;
        }	
      }
    }
    Mean = Mean * 100.;
    rms = rms / n_fibers ;
    rms = ::sqrtf(rms) * 100.;
    Mean_dark = Mean_dark * 100. *3. / 1000.; //normalization to take into acocunt different time windows 
    rms_dark = rms_dark / n_fibers ;
    rms_dark = ::sqrtf(rms_dark) * 100. *3. / 1000.;
    
//
    bx_message &msg = get_message(bx_message::info);
    msg << "BUNDLE =  "<<i <<std::endl;
    msg <<"Mean (%)= " <<Mean/i4_times<<" rms=  "<<rms/i4_times <<"     Dark rate= "<<Mean_dark/i4_times<< " rms (dark rate)= " <<rms_dark/i4_times <<std::endl;
    msg <<"Total number of fibers in bundle: " << n_fibers <<"  Non working fibers: " << n_broken <<"  Corresponding to logical channels: "<< std::endl;
    for (int i = 1; i < n_broken + 1; i++){
      msg << fiber[i]<<", " ;
    }
    msg <<std::endl;
    msg <<"Non working PMT: " << n_broken_dark <<"  Corresponding to logical channels: "<< std::endl;
    for (int i = 1; i < n_broken_dark + 1; i++){
      msg << pmt[i]<< ", ";
    }
      msg <<""<<std::endl;
      msg <<""<<std::endl;
      msg <<dispatch;
  }
  barn_interface::get()->store(barn_interface::file,bundles_laser,this);
  barn_interface::get()->store(barn_interface::file,channels_laser,this);
  barn_interface::get()->store(barn_interface::file,bundles_dark,this);
  barn_interface::get()->store(barn_interface::file,channels_dark,this);
  
  get_message(bx_message::info)<<"Number of times in bx_calib_fiber_bundle= " <<i4_times <<dispatch;
}


/*
 * $Log: bx_calib_fiber_bundle.cc,v $
 * Revision 1.8  2007/02/06 11:24:50  bcaccian
 *
 * Corrected a bug: "return ev" instead of "return 0" at the end of doit
 *
 * Revision 1.7  2006-08-21 11:20:08  razeto
 * Updated to new barn_interface
 *
 * Revision 1.6  2006/05/08 17:31:31  razeto
 * Added db_channel patch for fadc (sent to the mailing list)
 *
 * Revision 1.5  2006/01/09 16:21:25  razeto
 * Updated to the new root_barn target
 *
 * Revision 1.4  2004/11/26 15:25:11  razeto
 * Thanks to Davide who so kindly underlined a syntactical error, I fixed the "Maintainer" word
 *
 * Revision 1.3  2004/11/26 14:16:54  razeto
 * Added Mantainer field
 *
 * Revision 1.2  2004/11/04 11:52:00  bcaccian
 * CVSChanged a histogram name conflicting with an other module
 *
 * Revision 1.1  2004/10/01 10:10:16  bcaccian
 * Added
 *
 */
