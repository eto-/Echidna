/* BOREXINO Reconstruction program
 *
 * Author: Livia Ludhova and Gemma Testera
 */
#include "bx_calib_laben_charge_tt1.hh"
#include "messenger.hh"
#include "bx_dbi.hh"
#include "db_profile.hh"
#include "db_run.hh"
#include "db_channel.hh"
#include "barn_interface.hh"
#include <fstream>

// ctor
bx_calib_laben_charge_tt1::bx_calib_laben_charge_tt1 (): bx_base_module("bx_calib_laben_charge_tt1", bx_base_module::main_loop) {
  require_event_stage (bx_detector::laben, bx_base_event::clustered);
  require_trigger_type (bx_trigger_event::neutrino);
}

void bx_calib_laben_charge_tt1::begin () {
  get_message(bx_message::debug) << "begin" << dispatch;
  // creation of the 2-dim histogram 
  int nch = constants::laben::channels;
  bx_adc_charge_tt1 = new TH2S ("bx_adc_charge_tt1", "ADC Charge Calibration Histogram", nch, 1, nch + 1, 511, -255.5, 255.5);
  barn_interface::get ()->store (barn_interface::file, bx_adc_charge_tt1, this);
   

  Amp_vs_lg = new TH1F ("Amp_vs_lg", "Amp_vs_lg", nch, 1, nch + 1);
  barn_interface::get ()->store (barn_interface::file, Amp_vs_lg, this);
  C1_vs_lg = new TH1F ("C1_vs_lg", "C1_vs_lg", nch, 1, nch + 1);
  barn_interface::get ()->store (barn_interface::file, C1_vs_lg, this);
  Sig1_vs_lg = new TH1F ("Sig1_vs_lg", "Sig1_vs_lg", nch, 1, nch + 1);
  barn_interface::get ()->store (barn_interface::file, Sig1_vs_lg, this);
  R12_vs_lg = new TH1F ("R12_vs_lg", "R12_vs_lg",  nch, 1, nch + 1);
  barn_interface::get ()->store (barn_interface::file, R12_vs_lg, this);
  R13_vs_lg = new TH1F ("R13_vs_lg", "R13_vs_lg",nch, 1, nch + 1);
  barn_interface::get ()->store (barn_interface::file, R13_vs_lg, this);
  Chi2_vs_lg = new TH1F ("Chi2_vs_lg", "Chi2_vs_lg",  nch, 1, nch + 1);
  barn_interface::get ()->store (barn_interface::file, Chi2_vs_lg, this);
  Mean_vs_lg = new TH1F ("Mean_vs_lg", "Mean_vs_lg",  nch, 1, nch + 1);
  barn_interface::get ()->store (barn_interface::file, Mean_vs_lg, this);
  Rms_vs_lg = new TH1F ("Rms_vs_lg", "Rms_vs_lg",  nch, 1, nch + 1);
  barn_interface::get ()->store (barn_interface::file, Rms_vs_lg, this);
  P0_vs_lg = new TH1F ("P0_vs_lg", "P0_vs_lg",  nch, 1, nch + 1);
  barn_interface::get ()->store (barn_interface::file, P0_vs_lg, this);


  hC1 = new TH1F ("hC1", "hC1", 600, 0., 150.);
  barn_interface::get ()->store (barn_interface::file, hC1, this);
  hSig1 = new TH1F ("hSig1", "hSig1", 150, 0., 50.);
  barn_interface::get ()->store (barn_interface::file, hSig1, this);
  hR12 = new TH1F ("hR12", "hR12",250, 0., 0.5);
  barn_interface::get ()->store (barn_interface::file, hR12, this);
  hR13 = new TH1F ("hR13", "hR13", 200 , 0., 0.2);
  barn_interface::get ()->store (barn_interface::file, hR13, this);
  hChi2 = new TH1F ("hChi2", "hChi2",150, 0., 3.);
  barn_interface::get ()->store (barn_interface::file, hChi2, this);
  hMean = new TH1F ("hMean", "hMean",600, 0., 150.);
  barn_interface::get ()->store (barn_interface::file, hMean, this);
  hRms = new TH1F ("hRms", "hRms",150, 0., 50.);
  barn_interface::get ()->store (barn_interface::file, hRms, this);
  hP0 = new TH1F ("hP0", "hP0",200, 0.9, 1.);
  barn_interface::get ()->store (barn_interface::file, hP0, this);


  // creation of the result vectors
   f_charge_peak = new float[nch]; 
   f_charge_sigma  = new float[nch]; 
   f_charge_mean  = new float[nch]; 
   f_charge_rms  = new float[nch]; 
   f_p0  = new float[nch]; 

   nTriggers = 0;
   current_run = -10;
 }

 //DO IT
 // filling of the 2-dim histogram
 bx_echidna_event* bx_calib_laben_charge_tt1::doit (bx_echidna_event *ev) {

   if(current_run == -10) current_run = ev->get_run_number ();
   //select 1 cluster events and btb inputs
   if (ev->get_laben ().get_nclusters () == 1 && ev->get_trigger().get_btb_inputs() == 0) {
     //select low energy events an dposition R < 4.2
     int nhits =  ev->get_laben ().get_cluster (0).get_clustered_nhits (); 
     if(nhits < 100) {
       double radius = ev->get_laben().get_cluster(0).get_baricenter().get_r ();
       if (radius < 4.2){
	   nTriggers++;
	 for(int i = 0; i < nhits; i++){
	   int lg = ev->get_laben ().get_cluster (0).get_clustered_hit(i).get_decoded_hit ().get_raw_hit ().get_logical_channel ();
	   float adc = ev->get_laben ().get_cluster (0).get_clustered_hit(i).get_decoded_hit ().get_charge_bin();
	   bx_adc_charge_tt1->Fill (lg, adc);
	}
      }
    }
  }
  return ev;
}

 double G3(double *x, double *par) {
    double xx =x[0];
    double rad= 1./sqrt(3.1415926*2);
    double m1      = par[1];
    double m2      = 2*par[1];
    double m3      = 3*par[1];
    double sigma1  = par[2];
    double sigma2  = sqrt(2.)*par[2];
    double sigma3  = sqrt(3.)* par[2];
    double arg1= pow((xx-m1),2)/(2*pow(sigma1,2));
    double arg2=pow((xx-m2),2)/(2*pow(sigma2,2));
    double arg3=pow((xx-m3),2)/(2*pow(sigma3,2));
    //double fun= par[0]*( exp(-arg1) + par[3]*exp(-arg2) + par[4]*exp(-arg3));
    double fun= rad * par[0] * (exp(-arg1)/sigma1 + par[3]*exp(-arg2)/sigma2 + par[4]*exp(-arg3)/sigma3);
    return fun;
  }



//END
void bx_calib_laben_charge_tt1::end () {

   db_run& run_info = bx_dbi::get()->get_run ();
  
  //initialize vectors which will set the visitors at the end to the values from the previous run
  for (int i = 0; i < constants::laben::channels; i++) {
    f_charge_peak[i] = run_info.get_laben_charge_tt1_peak (i + 1);
    f_charge_sigma[i]  = run_info.get_laben_charge_tt1_sigma (i + 1);
    f_charge_mean[i]  = run_info.get_laben_charge_tt1_mean (i + 1);
    f_charge_rms[i]  = run_info.get_laben_charge_tt1_rms (i + 1);
    f_p0[i]  = run_info.get_laben_charge_tt1_p0 (i + 1);
  }
  
  
  // mean value of the peak for the projection of all channels (mean of the run)
  TH1D *mean_proj = bx_adc_charge_tt1->ProjectionY("mean_proj");
  mean_proj->SetAxisRange(5,220);
  float peak_run = mean_proj->GetBinCenter(mean_proj->GetMaximumBin());
  float mean_run = mean_proj->GetMean ();
  int entries = (Int_t) mean_proj->Integral();
  float mean_entries = (float) entries/(constants::laben::channels);
  mean_proj->Delete();
  get_message(bx_message::log) << "Laser Peak (projection for all lg) in this run: " << peak_run << " ADC bin" << dispatch;
  get_message(bx_message::log) << "Mean value (projection for all lg) in this run: " << mean_run << " ADC bins" << dispatch;
  get_message(bx_message::log) << "Mean number of entries per Channel: " << mean_entries << dispatch;
  
  // fit function for the single channel histos

  int not_connected = 0;
  int no_light = 0;
  int too_low = 0;
  int not_conv = 0;
  int large_drift = 0;

    //some parameters for the fit 
  int MaxBin;
  double maxBinCont; 
  double mean;
  double rms;
  double P0 = 0;
  double Mu = 0;

    //fit paramteres
  double full_thresh = 15;
  int maxFullBin = 0; //

    //number of bins with negative charge
  int ofs = 255; 

  double Amp1, erAmp1;
  double C1, erC1;    
  double Sig1, erSig1;    
  double r21, er21;
  double r31, er31; 
  double Chi2 = -1;
  double NDF = 0;
  double Chi2_NDF = -1 ;

    //looop on single channels
  for (int i = 0; i < constants::laben::channels; i++) {

      // discard non Pmt channels (empty)  and PMTs which are disconnected
    const db_channel_laben &ch_info = dynamic_cast <const db_channel_laben&> (bx_dbi::get ()->get_channel (i + 1));
    if (!ch_info.is_ordinary() || ch_info.is_pmt_disconnected()){
      get_message(bx_message::log) << "Channel " << i+1 << " not connected to any Pmt" << dispatch;
       not_connected++;
       //check the value which was written to DB before for peak postion, 
       //if the value from the previous run was having rms=0 (contains just mean value of that run), we set peak=peak_run
       //and set rms to 0
       if (f_charge_peak[i] < 5. || f_charge_peak[i] > 220. || f_charge_sigma[i] == 0. ) {
	 f_charge_peak[i] = peak_run;
	 f_charge_mean[i] = mean_run;
       }
       f_charge_sigma[i]  = 0;
       f_charge_rms[i]  = 0;
       f_p0[i]  = 0;
       continue;
    }    

     //single channels projections
    TH1D *proj = bx_adc_charge_tt1->ProjectionY("proj", i+1, i+1);
    proj->SetAxisRange(5, 220); //overflow peak starts from 220 in some channels
     //single channel mean and rms
    mean = proj->GetMean();
    rms = proj->GetRMS();
    entries = (int) proj->Integral();
    P0 = 1. - ((double) entries / (double) nTriggers);
    if (P0 < 1.) Mu = -log(P0);
    
    get_message(bx_message::log) << "Pmt in lg " << i+1 << ": entries = " << entries << ", mean = " << mean << ", rms = "<<rms<<dispatch;
    get_message(bx_message::log) << "Pmt in lg " << i+1 <<  ", P0 = " << P0 << ", Mu = " << Mu <<", <q1> = " << mean*(1.-Mu/2.) << dispatch;


    // discard channels with no direct light
    if (entries < (int) mean_entries/100 || entries < 100) {
      get_message(bx_message::log) << "Connected Pmt in ordinary lg " << i+1 << ": low statistics; entries = " << entries << dispatch;
      proj->Delete();
      no_light++;  
      //check the value which was written to DB before for peak postion, 
      //if the value from the previous run was having rms=0 (contains just mean value of that run), we set peak=peak_run
      //and set rms to 0
      if (f_charge_peak[i] < 5. || f_charge_peak[i] > 220. || f_charge_sigma[i] == 0. ) {
	f_charge_peak[i] = peak_run;
	f_charge_mean[i] = mean_run;
      }
      f_charge_sigma[i]  = 0;
      f_charge_rms[i]  = 0;
      f_p0[i]  = 0;
      continue;
    }

 
    //fitting 3 Gaussians
    
    //search for maximum full bin
    int minFullBin = -255;
    for(int ibin = 1 + ofs; ibin < (255 +ofs); ibin++) {
      if(proj->GetBinContent(ibin) > full_thresh && minFullBin == -255)  minFullBin = ibin -ofs;
      if(proj->GetBinContent(ibin) > full_thresh)  maxFullBin = ibin -ofs;
    }

    int FitResult;

    if (maxFullBin < 10) {    //if only low energy below 10 adc bins are filled, not fitting
      f_charge_peak[i] = proj->GetMaximumBin () - ofs;
      f_charge_sigma[i] = rms;
      f_charge_mean[i] = mean;
      f_charge_rms[i] = rms;
      f_p0[i] = P0;
      get_message(bx_message::log) << "Connected Pmt in ordinary lg " << i+1 << " hits only below 10 ADC channels, P0: " << P0 << ", mean = " << mean << ", rms = " << rms << ", entries = " << entries  << dispatch;
      too_low ++;
      continue;
    }
    else { //we fit

      proj->SetAxisRange(-255, 255);
      TF1 *fr = new TF1("myfun",G3,-255.,255.,5);
    
        //position of maximum (normally cca. ofs + 25 = 280)
      MaxBin = proj->GetMaximumBin ();
      if (MaxBin > 10+ofs) {//check if max is above 10
	maxBinCont = proj->GetBinContent (MaxBin); 
      }
      else {   //maximum is below 10
	maxBinCont = proj->GetBinContent(int(mean)+ofs);
	MaxBin = int(mean)+ofs;
      }
    
      // set and limit fit paramteres
      fr->SetParameter(0,maxBinCont*rms); 
      fr->SetParName(0,"A"); 
      
      fr->SetParameter(1,MaxBin - ofs);  
      fr->SetParLimits(1,8,260);
      fr->SetParName(1,"Center1");
      
      fr->SetParameter(2,(MaxBin-ofs)/3.);	 
      fr->SetParLimits(2,5,80);
      fr->SetParName(2,"sigma");
      
      fr->SetParameter(3,0.02); 
      if(rms > 15) fr->SetParLimits(3,0.,0.5); 
      else fr->SetParLimits(3,0.,0.2); 
      fr->SetParName(3,"r21");
      
      fr->SetParameter(4,0.005); 
      if(rms > 15) fr->SetParLimits(4,0.,0.2); 
      else fr->SetParLimits(4,0.,0.1); 
      fr->SetParName(4,"r31");
      
      proj->GetXaxis()->SetRangeUser(-10,100);
      
      double low_fit_region = (MaxBin-ofs)/2.;
      if (minFullBin > low_fit_region) low_fit_region = minFullBin;
      double high_fit_region = 3.5 * (MaxBin-ofs);
      if (maxFullBin < high_fit_region) high_fit_region = maxFullBin;
      
      FitResult = proj->Fit("myfun","LQN","",low_fit_region, high_fit_region);
      
      Amp1 = fr->GetParameter(0);
      erAmp1 = fr->GetParError(0);
      
      C1   = fr->GetParameter(1);	
      erC1 = fr->GetParError(1);
      
      Sig1 = fr->GetParameter(2);
      erSig1 = fr->GetParError(2);
      
      r21 = fr->GetParameter(3);
      er21 = fr->GetParError(3);
      
      r31  = fr->GetParameter(4);
      er31 = fr->GetParError(4);
      
      Chi2=fr->GetChisquare();
      NDF=fr->GetNDF();
      
      if(NDF!=0)  Chi2_NDF = Chi2/NDF;
      if(NDF==0)  Chi2_NDF = -1.;
      
      Amp_vs_lg->SetBinContent(i+1, Amp1);
      Amp_vs_lg->SetBinError(i+1, erAmp1);
      
      hC1->Fill(C1);
      C1_vs_lg->SetBinContent(i+1, C1);
      C1_vs_lg->SetBinError(i+1, erC1);
      
      hSig1->Fill(Sig1);
      Sig1_vs_lg->SetBinContent(i+1, Sig1);
      Sig1_vs_lg->SetBinError(i+1, erSig1);
      
      hR12->Fill(r21);
      R12_vs_lg->SetBinContent(i+1, r21);
      R12_vs_lg->SetBinError(i+1, er21);
      
      hR13->Fill(r31);
      R13_vs_lg->SetBinContent(i+1, r31);
      R13_vs_lg->SetBinError(i+1, er31);
      
      hChi2->Fill(Chi2_NDF);
      Chi2_vs_lg->SetBinContent(i+1, Chi2_NDF); 
      
      hP0->Fill(P0);
      P0_vs_lg->SetBinContent(i+1, P0); 

      fr->Delete();

      //logs about the fit and mean for each channel
      get_message(bx_message::log) << "Pmt in lg " << i+1 << ": FitResult ok? " << FitResult << " peak from fit = " << C1 << ", Sigma = " << Sig1 << ", chi2/NDF = " << Chi2_NDF <<dispatch;
    }
    
    //warning for potentially bad PMTs (reset the vector for visitors to mean values, not the fit results!)
    if (Chi2_NDF > 3. && Chi2_NDF < 10.) 
      get_message(bx_message::log) <<"Increased Chi2 for PMT in lg  " << i+1 << dispatch;


    if (Chi2_NDF > 10.) {
      //    if ( FitResult!=0 || Chi2_NDF > 2) {
      not_conv++;  
      get_message(bx_message::log)<<" For LG = " << i+1 << " fit failed"<< dispatch;
      
        //for these channels , the mean value is set (if ok)
        //if mean value bad, values from the previous run kept (if they were ok) 
        //if also previous run values bad, mean of this run is set
      get_message(bx_message::log) <<"(Potentially) Bad Pmt (Lg = " << i+1 << dispatch;
      if (mean > 5. && rms > 1.) {
	f_charge_peak[i] = mean;
	f_charge_sigma[i]  = rms;
	f_charge_mean[i] = mean;
	f_charge_rms[i]  = rms;
	f_p0[i] = P0;
      }
      else if (run_info.get_laben_charge_tt1_peak (i+1) < 5. ||  run_info.get_laben_charge_tt1_peak (i+1) > 220. ||  run_info.get_laben_charge_tt1_sigma (i+1) == 0. ) {//checking values from previous run
	//we should never arrive here
	get_message(bx_message::log) << "Pmt " << i+1 << "STRANGE CHECK IT" << dispatch;
	f_charge_peak[i] = peak_run;
	f_charge_sigma[i]  = 0;
	f_charge_mean[i] = mean_run;
	f_charge_rms[i]  = 0;
	f_p0[i] = P0;
	continue ; //to skip big drift warnings
      }
    }
    else {
      f_charge_peak[i] = C1;
      f_charge_sigma[i]  = Sig1;
      f_charge_mean[i] = mean;
      f_charge_rms[i]  = rms;
      f_p0[i] = P0;
      if (C1 > 30.) get_message(bx_message::log) << "Pmt " << i+1 << "too high gain" << dispatch;
      if (C1 < 20.) get_message(bx_message::log) << "Pmt " << i+1 << "too low gain" << dispatch;
    }

    // warn for channels with large peak (=fit or mean) drift compared to previous run
    if (fabs (run_info.get_laben_charge_tt1_peak (i+1) - f_charge_peak[i] ) > 10 && run_info.get_laben_charge_tt1_sigma (i+1) != 0 && f_charge_sigma[i] != 0) {
      get_message(bx_message::log) << "Pmt " << i+1 << ": large peak drift; new = " <<  f_charge_peak[i] << " , rms " << f_charge_sigma[i]
				   << ", old = " << run_info.get_laben_charge_tt1_peak(i+1) << "; sigma = " <<  run_info.get_laben_charge_tt1_sigma (i+1) << dispatch;
      large_drift++;
    }
    
    proj->Delete();
  }//end of loop in channels	 
  
  
  int good_ch = constants::laben::channels - not_connected - no_light - too_low - not_conv;
      
  get_message(bx_message::info) << not_connected << " channels not connected to any Pmt" << dispatch;
  get_message(bx_message::info) << no_light << " channels with no direct light" << dispatch;
  get_message(bx_message::info) << not_conv << " channels with bad calibration fit" << dispatch;
  get_message(bx_message::info) << too_low << " channels with hits only below 10" << dispatch;
  get_message(bx_message::info) << large_drift << " channels have large peak drift" << dispatch;
  get_message(bx_message::info) << good_ch << " channels correctly calibrated" << dispatch;

  //find last laser and neutrino calib run
  //  std::ostringstream where_str;
  //bx_dbi *dbi = bx_dbi::get ();
  //int last_laser_pmt_calib_run = dbi->query (bx_dbi::bx_calib, "\"LaserPmtCalibration\"", where_str.str (), "MAX(\"RunNumber\")", "", -1)["max"][0].get_int ();
  //int last_neutrino_pmt_calib_run = dbi->query (bx_dbi::bx_calib, "\"NeutrinoPmtCalibration\"", where_str.str (), "MAX(\"RunNumber\")", "", -1)["max"][0].get_int ();
  //std::cout << "Laser calib is " << last_laser_pmt_calib_run << " and neutrino run " << last_neutrino_pmt_calib_run << " and current run " << current_run << std::endl;

    //fill the visitors to db_run
  int final_db_write = 0;
  if (get_parameter ("db_write").get_bool () == 2) final_db_write = 1;
  if (nTriggers > get_parameter ("min_trig").get_int () && get_parameter ("db_write").get_bool () && (float)good_ch/(constants::laben::channels - not_connected) > 0.8) final_db_write = 1;

  if (final_db_write) {
    //set the visitors
    for (int i=0; i<constants::laben::channels; i++) {
      run_info.set_laben_charge_tt1_peak (i + 1, f_charge_peak[i], this);
      run_info.set_laben_charge_tt1_sigma (i + 1, f_charge_sigma[i], this);
      run_info.set_laben_charge_tt1_mean (i + 1, f_charge_mean[i], this);
      run_info.set_laben_charge_tt1_rms (i + 1, f_charge_rms[i], this);
      run_info.set_laben_charge_tt1_p0 (i + 1, f_p0[i], this);
      get_message(bx_message::log) << "Values in DB: lG: "<< i + 1 << " peak " <<  f_charge_peak[i] << " sigma " <<  f_charge_sigma[i]  << " mean " << f_charge_mean[i] << " rms " << f_charge_rms[i] << " P0 " << f_p0[i]  << dispatch;
    }
    run_info.write_laben_tt1_charge (true, this);
  }
  else if (get_parameter ("db_write").get_bool ()){
    if (nTriggers < get_parameter ("min_trig").get_int ()) get_message(bx_message::warn) << "Not enough statistics: DB tables will not be updated!" << dispatch;
    if ((float)good_ch/(constants::laben::channels - not_connected) < 0.8) get_message(bx_message::warn) << "Too many channels badly calibrated: DB tables will not be updated!" << dispatch;
  }
  get_message(bx_message::debug) << "end" << dispatch;
  
  // deallocate the vectors
  delete [] f_charge_peak;
  delete [] f_charge_sigma;
  delete [] f_charge_mean;
  delete [] f_charge_rms;
  delete [] f_p0;
}
