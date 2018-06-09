/* BOREXINO Reconstruction program
 *
 * Author: Davide D'Angelo <davide.dangelo@mi.infn.it>
 * Based on code from Maria Elena Monzani
 * Maintainer: Michael Wurm <mwurm@ph.tum.de>
 *
 * $Id: bx_calib_muon_charge_peak.cc,v 1.13 2012/06/05 14:36:12 meindl Exp $
 *
 * Implementation of bx_calib_muon_channel
 *
 */
#include "bx_calib_muon_charge_peak.hh"
#include "messenger.hh"
#include "bx_dbi.hh"
#include "db_profile.hh"
#include "db_run.hh"
#include "barn_interface.hh"
#include "TF1.h"
#include <cmath>


// ctor
bx_calib_muon_charge_peak::bx_calib_muon_charge_peak (): bx_base_module("bx_calib_muon_charge_peak", bx_base_module::main_loop) {
  require_event_stage (bx_detector::muon, bx_base_event::decoded);
  require_trigger_type (bx_trigger_event::laser394);
}

void bx_calib_muon_charge_peak::begin () {
  get_message(bx_message::debug) << "begin" << dispatch;

  // creation of the 2-dim histogram 
  int nch = constants::muon::channels;
  muon_charge_calib = new TH2S ("muon_charge_calib", "Charge Calibration Histogram",nch,0,nch, 256 ,  0, 256);
  barn_interface::get ()->store (barn_interface::file, muon_charge_calib, this);

  b_has_data = false;
  i4_nevents = 0;

  //get time offset and range parameters
  f4_time_allowance = get_parameter ("time_allowance").get_float();
  f4_time_offset    = get_parameter ("time_offset")   .get_float();
  f4_chi2_max       = get_parameter ("chi2_max")      .get_float();
  f4_chi2_min       = get_parameter ("chi2_min")      .get_float();
 
  if ( get_parameter ("write_calib_charge").get_bool () && bx_dbi::get()->get_run ().is_muon_laser_present () ) 
    get_message (bx_message::warn) << "Laser calibrations already present for this run." << dispatch; 
}

bx_echidna_event* bx_calib_muon_charge_peak::doit (bx_echidna_event *ev) {
  const bx_muon_event& er = ev->get_muon();

  int nhits = er.get_decoded_nhits();

  // Update has data field
  if (!nhits) return ev;
  
  i4_nevents++;
  b_has_data |= true;

  for(int i = 0; i < nhits; i++){
    const bx_muon_decoded_hit& hit = er.get_decoded_hit(i);
    unsigned short mch = hit.get_raw_hit().get_muon_channel();
    float time   = hit.get_time  ();
    float charge = hit.get_charge();
    // filling of the 2-dim histogram, only positive charge values
    if (charge>0 && fabs(time-f4_time_offset)<f4_time_allowance) muon_charge_calib->Fill(mch,charge);
  }

  return ev;
}

void bx_calib_muon_charge_peak::end () {
  if (!b_has_data) {
    get_message(bx_message::error) << "No laser data are present in Run: charge calibration not performed." << dispatch;
    return;
  }  

  // retrieve parameters
  float minimum_efficiency = get_parameter ("minimum_efficiency").get_float () ;
  float max_charge_offset  = get_parameter ("max_charge_offset" ).get_float () ;
  float min_charge_rms     = get_parameter ("min_charge_rms"    ).get_float () ;
  float max_charge_rms     = get_parameter ("max_charge_rms"    ).get_float () ;
  float charge_drift       = get_parameter ("charge_drift"      ).get_float () ;

  const db_profile& profile_info = bx_dbi::get()->get_profile ();
            db_run&     run_info = bx_dbi::get()->get_run     ();

  // mean value of the laser peak
  TH1D *mean_proj_charge = muon_charge_calib->ProjectionY("mean_proj_charge");
  int avg_entries = int(mean_proj_charge->Integral()/208.);
  float avg_peak;
  avg_peak = (float) mean_proj_charge->GetMean();

  mean_proj_charge->Delete();
  get_message(bx_message::log) << "Average number of entries: " << avg_entries             << dispatch;
  get_message(bx_message::log) << "Average charge peak: "       << avg_peak    << " ticks" << dispatch;

  // creation of the result vectors
  float *charge_peak   = new float[constants::muon::channels]; 
  float *charge_sigma  = new float[constants::muon::channels]; 
  std::fill_n (charge_peak ,constants::muon::channels, avg_peak);
  std::fill_n (charge_sigma,constants::muon::channels,       0.);

  // some ctrs
  int not_connected = 0, no_light = 0;
  int not_conv_charge = 0;
//  int no_charge_mean = 0;
  int large_drift_charge = 0;
  int good_in_charge = 0;


  // main loop
  for (int i=0; i<constants::muon::channels; i++) {
    // discard non Pmt channels
    if (profile_info.logical_channel_description (i+constants::muon::channel_offset+1) != db_profile::ordinary){
      charge_peak[i]=0.;
      get_message(bx_message::log) << "Mch " << i << " is service channel." << dispatch;
      not_connected++;
      continue;
    }
    TH1D *proj_charge = muon_charge_calib->ProjectionY("proj_charge", i+1, i+1);
    int entries = int(proj_charge->Integral());
    // discard channels with no direct light
    if (float(entries)/float(i4_nevents) < minimum_efficiency) {
      get_message(bx_message::log) << "Mch " << i << ": low statistics; entries = " << entries << " ; required minimum entries = " << minimum_efficiency * float(i4_nevents) << dispatch;
      proj_charge->Delete();
      no_light++;  
      continue;
    }



    // *** CHARGE CALIBRATION ***
    //
    TF1* func_charge = new TF1("func_charge","gaus");
    double cpar[3];

    //Find Gauss-Peak of spe-spektrum
    int maximumbin = proj_charge->GetMaximumBin();

    //Average over peak
    int n=1;	//1 bin left and right of the maximumbin
    cpar[0] = 0;
    if (maximumbin >=2) for (int q=-n;q<=n;q++) {cpar[0] += proj_charge->GetBinContent(maximumbin + q);}
    else cpar[0] = proj_charge->GetBinContent(maximumbin);
    cpar[0] /= float(2*n+1);
    cpar[1] = proj_charge->GetBinCenter(maximumbin);
    cpar[2] = max_charge_rms/2.;

//    get_message(bx_message::log) << "Mch " << i << ": charge : Fit input parameters : cpar[0] = " << cpar[0] << "; cpar[1] = " << cpar[1] << "; cpar[2] = " << cpar[2] << dispatch;

    func_charge->SetParameters(cpar);

    func_charge->SetParLimits(0,0,entries/2.);	//height of gaussian
    func_charge->SetParLimits(1,0,208);		//mean f gaussian
    func_charge->SetParLimits(2,0,208);		//sigma of gaussian
    double charge_reduc_chi;


   //Define Fit-Range
   double leftthreshold = 2./3. * cpar[0];
   double rightthreshold = 1./4. * cpar[0];
   int leftrangebin  = proj_charge->GetXaxis()->GetFirst();
   int rightrangebin = proj_charge->GetXaxis()->GetLast();
   bool leftrange = false; bool rightrange = false;

   //Search for left and right bin of the fit range
   for (int bin = leftrangebin; bin <=  maximumbin; bin++)
   {
      //Leftside
      if (bin < maximumbin && proj_charge->GetBinContent(bin) > leftthreshold && leftrange == false) {leftrangebin = bin; leftrange = true; break;}

      if (bin == maximumbin && leftrange == false) 
      {
         if (leftrangebin < maximumbin - 2) leftrangebin = maximumbin - 2;
	 leftrange = true;
      }
   }

   for (int bin = rightrangebin; bin >=  maximumbin; bin--)
   {
      //Rightside
      if (bin > maximumbin && proj_charge->GetBinContent(bin) > rightthreshold && rightrange == false) {rightrangebin = bin; rightrange = true; break;}

      if (bin == maximumbin && rightrange == false)
      {
         if (rightrangebin > maximumbin + 2) rightrangebin = maximumbin + 2;
         rightrange = true;
      }
   }

//   if (leftrange && rightrange) get_message(bx_message::log) << "Left and right side of Gauss-Peak found. Left side : " << proj_charge->GetBinCenter(leftrangebin) << " ; Righ side : " <<  proj_charge->GetBinCenter(rightrangebin) << dispatch;

    proj_charge->Fit("func_charge","QL0", "", proj_charge->GetBinCenter(leftrangebin), proj_charge->GetBinCenter(rightrangebin));
    charge_reduc_chi = func_charge->GetChisquare()/func_charge->GetNDF();
    func_charge->GetParameters(cpar);

   get_message(bx_message::log) << "Mch " << i << ": charge : Fit result           : cpar[0] = " << cpar[0] << "; cpar[1] = " << cpar[1] << "; cpar[2] = " << cpar[2] << "; reduced chisquare = " << charge_reduc_chi << dispatch;
    delete func_charge;


    bool fitfail = false;
    if (fabs (cpar[1]-avg_peak) > max_charge_offset)
    {get_message(bx_message::log) << "Mch " << i << ": charge : bad charge fit due to peak position." << dispatch; fitfail = true;};
 
    if (cpar[2] <= min_charge_rms || cpar[2] > max_charge_rms)
    {get_message(bx_message::log) << "Mch " << i << ": charge : bad charge fit due to bad sigma." << dispatch; fitfail = true;};

    if (cpar[2] > cpar[1])
    {get_message(bx_message::log) << "Mch " << i << ": charge : bad charge fit due to too large sigma (sigma > peak position)." << dispatch; fitfail = true;};

    if (cpar[1] < 1.75)
    {get_message(bx_message::log) << "Mch " << i << ": charge : peak position too small." << dispatch; fitfail = true;};
	       
    if (charge_reduc_chi > f4_chi2_max || charge_reduc_chi < f4_chi2_min)
    {get_message(bx_message::log) << "Mch " << i << ": charge : bad charge fit due to bad chisquare." << dispatch; fitfail = true;};

    if (fitfail){
       not_conv_charge++;  
       cpar[1] = (Double_t) proj_charge->GetMean();
       cpar[2] = (Double_t) proj_charge->GetRMS();
    }

     proj_charge->Delete();
     // warn for channels with large peak drift compared to previous run
     if ( fabs (run_info.get_muon_charge_peak (i+constants::muon::channel_offset+1) - cpar[1] ) > charge_drift && run_info.get_muon_charge_sigma (i+constants::muon::channel_offset+1) != 0 ) {
       get_message(bx_message::log) << "Mch " << i << ": charge : large peak drift in charge; new peak = " << cpar[1] << ", new sigma = " << cpar[2] <<
 	"; old peak = " << run_info.get_muon_charge_peak(i+constants::muon::channel_offset+1) <<  dispatch;
       large_drift_charge++;
       // if requested upon configuration, discard these data
       if (get_parameter ("large_drift").get_bool ()) continue;
     }
     charge_peak [i] = (float) cpar[1];
     charge_sigma[i] = (float) cpar[2];
     good_in_charge++;
 
     get_message(bx_message::debug) << "Mch " << i << ": charge : calibrated in charge; entries = " << entries 
 				   << " ; peak = " << cpar[1] << " ; sigma = " << cpar[2] << dispatch;

/**/get_message(bx_message::log) << "Mch " << i << ": charge : calibrated in charge; entries = " << entries 
				   << " ; peak = " << cpar[1] << " ; sigma = " << cpar[2] << dispatch;

  }


//  int good_ch_time   = constants::muon::channels - not_connected - no_light - not_conv_time   - large_drift_time;
//  int good_ch_charge = constants::muon::channels - not_connected - no_light - not_conv_charge - large_drift_charge;
  
  get_message(bx_message::info) << not_connected      << " channels not connected to a Pmt"  << dispatch;
  get_message(bx_message::info) << no_light           << " channels with no direct light"    << dispatch;
  get_message(bx_message::info) << not_conv_charge    << " channels with bad charge fit"     << dispatch;
  get_message(bx_message::info) << large_drift_charge << " channels have large charge drift" << dispatch;
  get_message(bx_message::info) << good_in_charge     << " channels with charge calibration" << dispatch;
 
  // fill the visitors to db_run
  run_info.set_muon_charge_peak (constants::muon::channel_offset+1, avg_peak, this);
  for (int i=1; i<constants::muon::channels; i++) {
    int lg = i + 1 + constants::muon::channel_offset;
    if (profile_info.logical_channel_description (lg) != db_profile::ordinary) continue;
    run_info.set_muon_charge_peak  (lg , charge_peak [i], this);
    run_info.set_muon_charge_sigma (lg , charge_sigma[i], this);
  }

  // deallocate the vectors
  delete [] charge_peak;
  delete [] charge_sigma;

  // decide whether to update DB tables
  if ( get_parameter ("write_calib_charge").get_bool () ) {
    if ( float(good_in_charge) / (constants::muon::channels-not_connected) > get_parameter ("update_threshold").get_float() )
      run_info.write_muon_laser_charge   (true, this);
    else
      get_message(bx_message::warn) << "Too may channels have bad charge calibration: DB charge table will not be updated." << dispatch;
  }

  get_message(bx_message::debug) << "end" << dispatch;
}

/*
 * $Log: bx_calib_muon_charge_peak.cc,v $
 * Revision 1.13  2012/06/05 14:36:12  meindl
 * Modified log-message to give additional information on channel calibration.
 *
 * Revision 1.12  2012-05-16 15:14:30  meindl
 * Muted log-message.
 *
 * Revision 1.11  2012-05-16 14:56:41  meindl
 * Modified range selection criteria for the charge fit and acceptance of charge calibration results.
 *
 * Revision 1.10  2012-04-06 08:55:54  wurm
 * increased range of charge calib histogram
 *
 * Revision 1.9  2008-12-15 14:31:43  wurm
 * lowered verbosity, added parameters for chi2
 *
 * Revision 1.8  2008-12-11 16:28:08  meindl
 * some more parameter tuning
 *
 * Revision 1.7  2008-12-11 16:10:29  meindl
 * changed fit parameter limits
 *
 * Revision 1.6  2008-12-11 15:46:50  meindl
 * included gauss fit with dynamic fit range; removed obsolete parts of the code
 *
 * Revision 1.5  2008-10-16 19:48:15  ddangelo
 * channel efficiency check improved,
 * time cut
 * some parameters tuning
 *
 * Revision 1.4  2008-10-14 15:25:32  wurm
 * debugging
 *
 * Revision 1.3  2008-10-14 14:32:28  wurm
 * split module of calib_muon_channel
 *
 * Revision 1.6  2008-10-07 14:42:54  wurm
 * omit charge calibration if time calibration files
 *
 * Revision 1.5  2008-10-07 14:16:27  wurm
 * modified time fit
 *
 * Revision 1.4  2008-10-01 09:44:23  wurm
 * adjusted histogram range of charge
 *
 * Revision 1.3  2008-09-29 13:32:47  wurm
 * rms for bad charge fits, further minor adjustments
 *
 * Revision 1.2  2006-09-11 20:56:40  ddangelo
 * added drift check, params readout, charge fit.
 * First params tuning, first run tests.
 * More to be done
 *
 * Revision 1.1  2006/09/11 14:12:38  ddangelo
 * bx_calib_muon_charge_peak and bx_calib_muon_time_align
 * replaced by
 * bx_calib_muon_channel
 *
 *
 */
