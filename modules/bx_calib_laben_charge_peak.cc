/* BOREXINO Reconstruction program
 *
 * Author: Maria Elena Monzani <monzani@mi.infn.it>
 * Maintainer: Oleg Smirnov, Livia Ludhova
 *
 * $Id: bx_calib_laben_charge_peak.cc,v 1.20 2008/11/26 14:22:38 ludhova Exp $ 
 *
 * Implementation of bx_calib_laben_charge_peak
 *
 */
#include "bx_calib_laben_charge_peak.hh"
#include "messenger.hh"
#include "bx_dbi.hh"
#include "db_profile.hh"
#include "db_run.hh"
#include "db_channel.hh"
#include "barn_interface.hh"
#include <fstream>

// ctor
bx_calib_laben_charge_peak::bx_calib_laben_charge_peak (): bx_base_module("bx_calib_laben_charge_peak", bx_base_module::main_loop) {
  require_event_stage (bx_detector::laben, bx_base_event::decoded);
  require_trigger_type (bx_trigger_event::laser394);
}

void bx_calib_laben_charge_peak::begin () {
  get_message(bx_message::debug) << "begin" << dispatch;
  // creation of the 2-dim histogram 
  int32_t nch = constants::laben::channels;
  bx_adc_charge_calib = new TH2S ("bx_adc_charge_calib", "ADC Charge Calibration Histogram", nch, 1, nch + 1, 511, -255.5, 255.5);
  barn_interface::get ()->store (barn_interface::file, bx_adc_charge_calib, this);

  bx_q_charge_calib = new TH2S ("bx_q_charge_calib", "Charge Calibration Histogram", nch, 1, nch + 1, 120, -1., 5.);
  barn_interface::get ()->store (barn_interface::file, bx_q_charge_calib, this);

  p0_vs_lg = new TH1F ("p0_vs_lg", "p0_vs_lg", nch, 1, nch + 1);
  barn_interface::get ()->store (barn_interface::file, p0_vs_lg, this);

  mu_vs_lg = new TH1F ("mu_vs_lg", "mu_vs_lg", nch, 1, nch + 1);
  barn_interface::get ()->store (barn_interface::file, mu_vs_lg, this);

  fit1_vs_lg = new TH1F ("fit1_vs_lg", "fit1_vs_lg", nch, 1, nch + 1);
  barn_interface::get ()->store (barn_interface::file, fit1_vs_lg, this);

  mean_vs_lg = new TH1F ("mean_vs_lg", "mean_vs_lg", nch, 1, nch + 1);
  barn_interface::get ()->store (barn_interface::file, mean_vs_lg, this);

  rms_vs_lg = new TH1F ("rms_vs_lg", "rms_vs_lg", nch, 1, nch + 1);
  barn_interface::get ()->store (barn_interface::file, rms_vs_lg, this);

  h_mean = new TH1F ("h_mean", "h_mean", 200, 0., 100.);
  barn_interface::get ()->store (barn_interface::file, h_mean, this);

  h_rms = new TH1F ("h_rms", "h_rms", 200, 0., 100.);
  barn_interface::get ()->store (barn_interface::file, h_rms, this);

  h_fit1 = new TH1F ("h_fit1", "h_fit1", 200, 0., 100.);
  barn_interface::get ()->store (barn_interface::file, h_fit1, this);

  
  // creation of the result vectors
  f_charge_peak = new float[nch]; 
  f_charge_sigma  = new float[nch]; 
  
  nTriggers = 0;
}

//DO IT
// filling of the 2-dim histogram
bx_echidna_event* bx_calib_laben_charge_peak::doit (bx_echidna_event *ev) {
  const bx_laben_event& el = ev->get_laben();
  // Update has data field
  if (el.get_raw_nhits ()) b_has_data = true;
  int32_t nhits = el.get_decoded_nhits();
  for(int32_t i = 0; i < nhits; i++){
    const bx_laben_decoded_hit& hit = el.get_decoded_hit(i);
    const bx_laben_raw_hit& rawhit = hit.get_raw_hit();
    uint16_t ch = rawhit.get_logical_channel();
    float chbin = hit.get_charge_bin();
    float q = hit.get_charge_pe();
    //TIME
    double raw_time =  hit.get_raw_time () - ev->get_laben ().get_laser_rawt ();
    if (raw_time > 610 && raw_time < 690) {
      bx_adc_charge_calib->Fill (ch,chbin);
      bx_q_charge_calib->Fill (ch,q);
    }
  }
  nTriggers++;
  return ev;
}

//END
void bx_calib_laben_charge_peak::end () {

  get_message(bx_message::log) << "Run has " << nTriggers << " laser triggers " << dispatch;
  
  //if no laser data...
  if (!b_has_data) {
    get_message(bx_message::error) << "No laser data are present in Run: charge calibration not performed." << dispatch;
    // deallocate the vectors
    delete [] f_charge_peak;
    delete [] f_charge_sigma;
    return;
  } 
  db_run& run_info = bx_dbi::get()->get_run ();

  //initialize vectors which will set the visitors at the end to the values from the previous run
  for (int32_t i = 0; i < constants::laben::channels; i++) {
    f_charge_peak[i] = run_info.get_laben_charge_peak (i + 1);
    f_charge_sigma[i]  = run_info.get_laben_charge_sigma (i + 1);
  }

  // mean value of the laser peak for the projection of all channels (mean of the run)
  TH1D *mean_proj = bx_adc_charge_calib->ProjectionY("mean_proj");
  mean_proj->SetAxisRange(5,250);
  float peak_run = mean_proj->GetBinCenter(mean_proj->GetMaximumBin());
  int32_t nent = (Int_t) mean_proj->Integral();
  float mean_entries = (float) nent/(constants::laben::channels);
  mean_proj->Delete();
  get_message(bx_message::log) << "Mean value of the Laser Peak in this run: " << peak_run << " ADC bin" << dispatch;
  get_message(bx_message::log) << "Mean number of entries per Channel: " << mean_entries << dispatch;
  
  // 1 Gauss: fit function for the single channel histos
  TF1* g0 = new TF1("g0","gaus",0.,255.);
  Double_t tpar[3];
 
  int32_t not_connected = 0;
  int32_t no_light = 0;
  int32_t not_conv = 0;
  int32_t large_drift = 0;

    //looop on single channels
  for (int32_t i = 0; i < constants::laben::channels; i++) {

      // discard non Pmt channels (empty)  and PMTs which are disconnected
    const db_channel_laben &ch_info = dynamic_cast <const db_channel_laben&> (bx_dbi::get ()->get_channel (i + 1));
    if (!ch_info.is_ordinary() || ch_info.is_pmt_disconnected()){
      get_message(bx_message::log) << "Channel " << i+1 << " not connected to any Pmt" << dispatch;
       not_connected++;
       //check the value which was written to DB before for peak postion, 
       //if the value from the previous run was having rms=0 (contains just mean value of that run), we set peak=peak_run
       //and set rms to 0
       if (f_charge_peak[i] < 5. || f_charge_peak[i] > 220. || f_charge_sigma[i] == 0. ) f_charge_peak[i] = peak_run;
       f_charge_sigma[i]  = 0;
       continue;
    }    

     //single channels projections
    TH1D *proj = bx_adc_charge_calib->ProjectionY("proj", i+1, i+1);
    proj->SetAxisRange(5, 220); //overflow peak starts from 220 in some channels
     //single channel mean and rms
    double mean = proj->GetMean();
    double rms = proj->GetRMS();
    nent = (int32_t) proj->Integral();
    

    // discard channels with no direct light
    if (nent < (int32_t) mean_entries/100 || nent < 100) {
      get_message(bx_message::log) << "Connected Pmt in ordinary lg " << i+1 << ": low statistics; entries = " << nent << dispatch;
      proj->Delete();
      no_light++;  
      //check the value which was written to DB before for peak postion, 
      //if the value from the previous run was having rms=0 (contains just mean value of that run), we set peak=peak_run
      //and set rms to 0
      if (f_charge_peak[i] < 5. || f_charge_peak[i] > 220. || f_charge_sigma[i] == 0. ) f_charge_peak[i] = peak_run;
      f_charge_sigma[i]  = 0;
      continue;
    }

         //fitting 1 Gaussian
    double rmin = mean - rms;
    if (rmin < 5) rmin = 5;
    double rmax = mean + 2*rms;
    if (rmax >= 220) rmax = 220;
    if (rmax <= rmin) rmax = rmin + 50;
    g0->SetRange(rmin, rmax);
    tpar[0] = proj->GetBinContent(proj->GetMaximumBin());
    tpar[1] = proj->GetBinCenter(proj->GetMaximumBin());
    tpar[2] = rms;
    g0->SetParameters(&tpar[0]);
    int32_t FitResult = proj->Fit ("g0","QRL0"); //if fit does not fail FitResult = 0
    g0->GetParameters(tpar);
 

    //calculate mu
    double P0 = (double) nent / (double) nTriggers;
    double Mu = 0;
    if (P0 < 1.) Mu = -log(1.-P0);
    p0_vs_lg->SetBinContent (i+1, P0);
    mu_vs_lg->SetBinContent (i+1, Mu);
   
     //logs about the fit and mean for each channel
    get_message(bx_message::log) << "Pmt in lg " << i+1 << ": entries = " << nent << ", mean = " << mean << ", rms = "<<rms<<dispatch;
    get_message(bx_message::log) << "Pmt in lg " << i+1 << ": FitResult ok? " << FitResult << " peak from fit = " << tpar[1] << ", sig = " << tpar[2] << ", chisq = " << g0->GetChisquare()<<"/" << rmax-rmin <<dispatch;
    get_message(bx_message::log) << "Pmt in lg " << i+1 <<  ", Mu = " << Mu <<", <q1> = " << mean*(1.-Mu/2.) << dispatch;
    fit1_vs_lg->SetBinContent (i+1, tpar[1]);
    mean_vs_lg->SetBinContent (i+1, mean);
    rms_vs_lg->SetBinContent (i+1, rms);
    h_fit1->Fill (tpar[1]);
    h_mean->Fill (mean);
    h_rms->Fill (rms);
   

    //check if fit converged 
    if (FitResult != 0) {
      not_conv++;  
      get_message(bx_message::log)<<" For LG = " << i+1 << " fit failed"<< dispatch;
    }     

    //warning for potentially bad PMTs (reset the vector for visitors to mean values, not the fit results!)
    if (FitResult!=0 || mean < 10 || tpar[1] < 10  || fabs(tpar[1] - mean) > 15.0) {
    //for these channels , the mean value is set (if ok)
    //if mean value bad, values from the previous run kept (if they were ok) 
    //if also previous run values bad, mean of this run is set
      get_message(bx_message::log) <<"(Potentially) Bad Pmt (Lg = " << i+1 << dispatch;
      if (mean > 5. && rms > 1.) {
	f_charge_peak[i] = mean;
	f_charge_sigma[i]  = rms;
      }
      else if ( run_info.get_laben_charge_peak (i+1) < 5. ||  run_info.get_laben_charge_peak (i+1) > 220. ||  run_info.get_laben_charge_sigma (i+1) == 0. ) {//checking values from previous run
	f_charge_peak[i] = peak_run;
	f_charge_sigma[i]  = 0;
	continue ; //to skip big drift warnings
      }
    }
    else {
      f_charge_peak[i] = tpar[1];
      f_charge_sigma[i]  = tpar[2];
    }

	
    // warn for channels with large peak (=fit or mean) drift compared to previous run
    if (fabs ( run_info.get_laben_charge_peak (i+1) - f_charge_peak[i] ) > 15 && run_info.get_laben_charge_sigma (i+1) != 0 && f_charge_sigma[i] != 0) {
      get_message(bx_message::log) << "Pmt " << i+1 << ": large peak drift; new = " <<  f_charge_peak[i] << " , rms " << f_charge_sigma[i]
				   << ", old = " << run_info.get_laben_charge_peak(i+1) << "; sigma = " <<  run_info.get_laben_charge_sigma (i+1) << dispatch;
      large_drift++;
      // if requested upon configuration, discard these data, use previous run of this channel (if was ok), or this run mean
      if (get_parameter ("large_drift").get_bool ()) {
	if ((run_info.get_laben_charge_peak (i+1) > 5. && run_info.get_laben_charge_peak (i+1) < 220.) || run_info.get_laben_charge_sigma (i+1) != 0. ) {
	  f_charge_peak[i] = run_info.get_laben_charge_peak (i+1);
	  f_charge_sigma[i]  = run_info.get_laben_charge_sigma (i+1);
	}
	else {
	  f_charge_peak[i] = peak_run;
	  f_charge_sigma[i] = 0;
	}
	
      }
    }
    proj->Delete();
  } //end of loop in channels	 
  
  g0->Delete();
  
  int32_t good_ch = constants::laben::channels - not_connected - no_light - not_conv;
      
  get_message(bx_message::info) << not_connected << " channels not connected to any Pmt" << dispatch;
  get_message(bx_message::info) << no_light << " channels with no direct light" << dispatch;
  get_message(bx_message::info) << not_conv << " channels with bad calibration fit" << dispatch;
  get_message(bx_message::info) << large_drift << " channels have large peak drift" << dispatch;
  get_message(bx_message::info) << good_ch << " channels correctly calibrated" << dispatch;
  
    //fill the visitors to db_run
  if (get_parameter ("write_calib").get_bool () && (float)good_ch/(constants::laben::channels - not_connected) > 0.8){
    //set the visitors
    for (int32_t i=0; i<constants::laben::channels; i++) {
      run_info.set_laben_charge_peak (i + 1, f_charge_peak[i], this);
      run_info.set_laben_charge_sigma (i + 1, f_charge_sigma[i], this);
      get_message(bx_message::log) << "Values in DB: lG: "<< i + 1 << " peak " <<  f_charge_peak[i] << " rms " <<  f_charge_sigma[i]  << dispatch;
    }
    run_info.write_laben_laser_charge (true, this);
  }
  else if (get_parameter ("write_calib").get_bool ())
    get_message(bx_message::warn) << "Too may channels have bad charge calibration: DB tables will not be updated!" << dispatch;

  
  get_message(bx_message::debug) << "end" << dispatch;
  
  // deallocate the vectors
  delete [] f_charge_peak;
  delete [] f_charge_sigma;
 
}
/*
 * $Log: bx_calib_laben_charge_peak.cc,v $
 * Revision 1.20  2008/11/26 14:22:38  ludhova
 * debug
 *
 * Revision 1.19  2008-11-19 17:04:04  ludhova
 * debug
 *
 * Revision 1.18  2008-11-19 15:50:38  ludhova
 * some new histos + new binning
 *
 * Revision 1.17  2008-09-24 12:14:33  razeto
 * Added parentheses for bool operations
 *
 * Revision 1.16  2007-12-15 10:17:57  ludhova
 * charge histos only for laser hits in the direct peak (time cuts)
 *
 * Revision 1.15  2007-11-26 15:25:57  ludhova
 * removed warn for laser validation
 *
 * Revision 1.14  2007-11-22 11:54:41  ludhova
 * cosmetics
 *
 * Revision 1.13  2007-11-21 13:28:07  ludhova
 * more corrections
 *
 * Revision 1.12  2007-11-21 12:57:59  ludhova
 * more corrections
 *
 * Revision 1.11  2007-11-20 16:31:02  ludhova
 * sorry... forgot to uncomment writing to DB
 *
 * Revision 1.10  2007-11-20 16:30:11  ludhova
 * different logics in settings of visitors
 *
 * Revision 1.9  2007-11-13 16:00:19  smirnov
 * Fit range reduced to the max of 250 (from 255)
 *
 * Revision 1.8  2007-11-13 15:31:22  smirnov
 * Redundant comments removed
 *
 * Revision 1.7  2007-11-13 14:12:51  smirnov
 * Gaussian fit tuned. Default values in the case of the fit failure will be set to the adc mean value
 *
 * Revision 1.6  2007-04-13 13:53:45  razeto
 * Added real charge histogram (for displaying only)
 *
 * Revision 1.5  2006-11-17 17:40:08  razeto
 * Never return 0
 *
 * Revision 1.4  2006/08/21 11:20:08  razeto
 * Updated to new barn_interface
 *
 * Revision 1.3  2006/01/02 21:27:18  razeto
 * Changed test target to file target for db_barn (missing maintainer)
 *
 * Revision 1.2  2005/05/05 17:13:40  monzani
 * S/Getters names updated according to the interface modification. Debugging.
 *
 * Revision 1.1  2005/05/05 10:03:48  monzani
 * Laser charge calibrations added (some debug still needed).
 *
 */
