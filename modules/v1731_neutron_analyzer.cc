/* BOREXINO Reconstruction program
 *
 * Author: Stefano Davini <stefano.davini@ge.infn.it>
 * Maintainer: Stefano Davini <stefano.davini@ge.infn.it>
 *
 * $Id:
 *
 * Declaration of namespace v1731_neutron_analyzer
 * 
 *
 * 
 */

#include "v1731_neutron_analyzer.hh"
using namespace v1731_neutron_analyzer;
using v1731_event::v1731event;

long v1731_neutron_analyzer::i4_over_thr(long i4_size, const float* f_in_samples, float f_thr){
  long i2_over=0;
  for (long i4_s=0; i4_s<i4_size; ++i4_s){
    if (f_in_samples[i4_s]>f_thr) ++i2_over;
  }
  return i2_over;
}

float v1731_neutron_analyzer::f_area_over_thr(long i4_size, const float* f_in_samples, float f_thr){
  float f_area=0.;
  for (long i4_s=0; i4_s<i4_size-1; ++i4_s){
    const float f_over = (0.5*(f_in_samples[i4_s] + f_in_samples[i4_s+1])) - f_thr;
    if (f_over>0) f_area+=f_over;
  }
  return f_area;
}


Double_t v1731_neutron_analyzer::analitic_convolution(Double_t* x, Double_t* par){

  const Double_t k = 1./par[4] - 1./par[5];
  const Double_t a = 1./par[4] + par[2]/(par[3]*par[3]);
  const Double_t b = 1./(2.*par[3]*par[3]);
  const Double_t alpha = sqrt(b);
  const Double_t beta  = a/(2.*sqrt(b));
  const Double_t F = (x[0] <= par[2] ? x[0] : par[2]);
  
  const Double_t gauss_conv = par[1]*( (x[0] < par[2] ? exp(-(x[0]-par[2])*(x[0]-par[2])/(2*par[3]*par[3])) : 0.) - 1./(alpha*par[4])*exp(-x[0]/par[4] - par[2]*par[2]/(2*par[3]*par[3]) + a*a/(4.*b)) * sqrt(M_PI)*0.5 * ( 1. + TMath::Erf(alpha*F - beta)) );
  
  const Double_t expo_conv = par[1] *( (x[0] >= par[2] ? exp(-(x[0]-par[2])/par[5]): 0.) - 1./(k*par[4]) * exp(par[2]/par[5] - x[0]/par[4]) * ( exp(k*x[0]) - exp(k*par[2]) ));
  
  return (par[0] + gauss_conv + (x[0] >= par[2] ? expo_conv : 0.));

}

/******************** v1731_neutron_analyzer *****************/


v1731_neutron_analyzer::neutron_analyzer::neutron_analyzer(const v1731event & v1731ev, bool b_always_muon){
  
  b_is_muon = false;
  analog_neutrons.clear();

  if (v1731ev.get_active(0,0)) {

    /* analog sum (0,0) analysis: is it a muon? */    
    if ((b_is_muon = is_muon(v1731ev, b_always_muon))){  // = not a mistake: intialization and evaluation of b_is_muon performed at same time
      
      /* analog sum analysis: is there any neutron? */
      if (v1731ev.get_zle_enabled()) {// zle enabled 
	zle_analog_sum_analyze(v1731ev);
      }
      else { // zle not enabled
	analog_sum_analyze(v1731ev);
      }
    }
  }
  
  if (analog_neutrons.get_number()>0){
    //analog_neutrons.show();

    /* neutron: single channel analyzer; */
    /* TO BE DONE */
  }
  
}


bool v1731_neutron_analyzer::neutron_analyzer::is_muon(const v1731event & v1731ev, bool b_always_muon){
#define SIZE_FIRST 22000  // 22000 --> 44 us (muon zone; muon is ~at 24 us)
#define SIZE_TO_MEAN 3000 // 3000  --> 6 us (analyzed as base in !zle mode)
  
  const long i4_total_size=v1731ev.get_event_length(0,0)/2;

  f4_muon_time = 24000;   // fake initilization 
  float f4_base_mean = 0.;// mean of first samples (base)
  int32_t i2_intersection= 0; // used to disriminate pulser/muon events
  
  if ((i4_total_size-SIZE_FIRST)>0){ // avoids array problem
    if(!v1731ev.get_zle_enabled()){
      for (int32_t i2m=0; i2m<SIZE_TO_MEAN; ++i2m){
	f4_base_mean += v1731ev.get_sample_at_bin(0, 0, i2m+1);
      }
      f4_base_mean /= SIZE_TO_MEAN;
    }
    else if (v1731ev.get_zle_enabled()){
      const int32_t i2_ngoodzones=v1731ev.get_number_of_good_zones(0,0);
      if (i2_ngoodzones>0) {
	const long i4_start_bin=v1731ev.get_begin_of_good_zone(0,0,0)/2;  // start bin of 1st good zone
	const long i4_length_1gz=v1731ev.get_length_of_good_zone(0,0,0)/2;// size of 1st good zone
	const int32_t i2_size=(i4_length_1gz < 100 ? i4_length_1gz : 100);    // size of base to compute mean
	for (int32_t i2m=0; i2m<i2_size; ++i2m){
	  f4_base_mean+=v1731ev.get_sample_at_bin(0,0, i4_start_bin+i2m+1);
	}
	f4_base_mean/=i2_size;
      }
    }
    
    bool is_over = false;
    bool first_edge = true; // used to acquire muon edge time
    const float f4_over = (f4_base_mean+125. > 255. ? 225. : f4_base_mean+125.);

    for (int32_t i2_s=SIZE_TO_MEAN; i2_s<SIZE_FIRST; ++i2_s){    
      const float f4_tmp_sample = v1731ev.get_sample_at_bin(0, 0, i2_s+1);
      
      if ((f4_tmp_sample>f4_over) && (!is_over)){
	++i2_intersection;
	is_over=true;
      }
      if (f4_tmp_sample<f4_over) is_over=false;
      if ((i2_intersection==1) && (is_over)){ // muon's peak
	if (first_edge){
	  f4_muon_time= static_cast <float> (2.*i2_s);
	  first_edge=false;
	}
      }
    }

    if (i2_intersection==1) return (true);      // 1 intersection -> muon with small overshoot
    else if (i2_intersection>4) return (true);  // 5+ intersections -> muon with big overshoot
    else return (false | b_always_muon);
  }

  return false;
}

void v1731_neutron_analyzer::neutron_analyzer::analog_sum_analyze(const v1731event & v1731ev){

  const long i4_total_size=v1731ev.get_event_length(0,0)/2;
  const long i4_other_size=i4_total_size-SIZE_FIRST; // size of "after muon" array
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,0,0)
  double* f_other_samples= new double [i4_other_size];
  double* f_other_dest=new double [i4_other_size];
#else
  float* f_other_samples= new float [i4_other_size];
  float* f_other_dest=new float [i4_other_size];
#endif
  v1731ev.get_sample_array(0, 0, f_other_samples, SIZE_FIRST, i4_other_size);
  
  TSpectrum S;         // primary peak serching:
  const Int_t i2_peakfounds=S.SearchHighRes(f_other_samples, f_other_dest, i4_other_size, 25, 60, true, 2, true, 3);
  float* f_peak_pos=S.GetPositionX();  // array of peak positions
  
  if (i2_peakfounds>0){
    
    float f_tmp_pk_val, f_tmp_max_val, f_tmp_FWHM, f_tmp_area, f_tmp_mean;
    const int32_t i2_near_size= 1000;  // 2us (right,left)
    
    for (int32_t i2_p=0; i2_p<i2_peakfounds; ++i2_p){
      const long i4_tmp_pk_pos = static_cast <long> (f_peak_pos[i2_p]);
      int32_t i2_tmp_near= i2_near_size;                // half length of near zone
      long i4_begin_near, i4_end_near;
      
      do {
	i4_begin_near= static_cast<long> (i4_tmp_pk_pos - i2_tmp_near);   // begin of near zone   
	i4_end_near= static_cast<long> (i4_tmp_pk_pos + i2_tmp_near);
	
	if ((i4_begin_near<0) || (i4_end_near>i4_other_size)){   //avoids filling/allocation errors
	  i2_tmp_near/=2;
	} 
      } while ((i4_begin_near<0) || (i4_end_near>i4_other_size));
           
      const float* f_samples_near= &f_other_samples[i4_begin_near]; //2us left and right the peak -> near zone
      
      f_tmp_pk_val= f_other_samples[i4_tmp_pk_pos];                  // peak value
      f_tmp_mean= TMath::Mean(2*i2_tmp_near, f_samples_near);        // mean of the zone 
      const float f_HM= (f_tmp_pk_val + f_tmp_mean)*0.5;             // Half Maximum of the peak 
      
      f_tmp_max_val= TMath::MaxElement(2*i2_tmp_near, f_samples_near);
      
      const int32_t i2_FWHM= v1731_neutron_analyzer::i4_over_thr(2*i2_tmp_near, f_samples_near, f_HM);// Full Width HM of the peak
      f_tmp_FWHM = static_cast<float> (i2_FWHM);        
      
      long i4_start_nearEST;
      int32_t i2_size_nearEST;
      if ((i4_tmp_pk_pos - 2*i2_FWHM)< 0) {
	i4_start_nearEST=i4_begin_near;
	i2_size_nearEST=2*i2_tmp_near;
      }
      else {
	i4_start_nearEST=i4_tmp_pk_pos-2*i2_FWHM;
	i2_size_nearEST=4*i2_FWHM;
      }
      
      const float* f_samples_nearEST= &f_other_samples[i4_start_nearEST]; // points to nearEST zone (base and peak)
      f_tmp_area= v1731_neutron_analyzer::f_area_over_thr(i2_size_nearEST, f_samples_nearEST, f_tmp_mean); //area of nearEST zone
      
      if ((f_tmp_pk_val-f_tmp_mean>30) || (f_tmp_max_val-f_tmp_mean>30)){  // threshold
	const float f_time=2.*(SIZE_FIRST + f_peak_pos[i2_p]) - f4_muon_time;
	v1731_neutron_analyzer::pulse tmp_pulse;
	tmp_pulse.f8_time   = f_time;
	tmp_pulse.f4_charge = 2.* f_tmp_area;
	tmp_pulse.f4_rise  = 2.* f_tmp_FWHM;
	tmp_pulse.f4_peak   = f_tmp_max_val - f_tmp_mean;      
	analog_neutrons.push(tmp_pulse);
      } 
    }
    
  }
  
  delete [] f_other_samples;
  delete [] f_other_dest;
}

void v1731_neutron_analyzer::neutron_analyzer::zle_analog_sum_analyze(const v1731event & v1731ev){
  const int32_t i2_num_good_zones = v1731ev.get_number_of_good_zones(0,0);
  if (i2_num_good_zones<2) return;

  const long i4_start_analysis_smp = static_cast<long>((f4_muon_time + 14000)/2); // muon overshoot expires after 10-20us

  for (int32_t i2_gz=0; i2_gz<i2_num_good_zones; ++i2_gz){
    /* analysis of single good zone */
    const long i4_smp_start = v1731ev.get_begin_of_good_zone(0,0,i2_gz)/2;
    const long i4_smp_size  = v1731ev.get_length_of_good_zone(0,0,i2_gz)/2;

    if ((i4_smp_size!=0) && (i4_smp_start>i4_start_analysis_smp)){
      float* f_zone_samples = new float [i4_smp_size];
      float* f_zone_dest    = new float [i4_smp_size];
      v1731ev.get_sample_array(0, 0, f_zone_samples, i4_smp_start, i4_smp_size);
      peak_zone pz(i4_smp_size, f_zone_samples, f_zone_dest,i4_smp_start, f4_muon_time);
      for (unsigned u4p=0; u4p<pz.get_npeak_found(); ++u4p){
	if (true){ // threshold pulse / neutron
	  analog_neutrons.push(pz,u4p);
	}
      }
    delete [] f_zone_samples;
    delete [] f_zone_dest;
    }
  } 
}


/***************** peak_zone **********************/


v1731_neutron_analyzer::peak_zone::peak_zone(long i4_in_size, float* f_samples, float* f_samples_clean, long i4_smp_start, float f_mu_time){  
  i4_size = i4_in_size;
  
  TSpectrum S;
  const float f_pk_sigma = 25.;
  const double f_thres= 50.;
  const int32_t i2_iter = 2;
  const int32_t i2_aver_win = 3;
    
  i2_npeak_found=S.SearchHighRes(f_samples, f_samples_clean, i4_size, f_pk_sigma, f_thres, kTRUE, i2_iter, kTRUE, i2_aver_win);

  if (i2_npeak_found>0){
    Float_t* f_peak_pos = S.GetPositionX();                   // peak positions
    Float_t* f_peak_pos_sort = new float [i2_npeak_found];    // will contain peak positions sorted temporaly
    int32_t* i2_index = new int32_t [i2_npeak_found];
    TMath::Sort(i2_npeak_found, f_peak_pos, i2_index, false); // see TMath::Sort for reference
    for (int32_t i2_s=0; i2_s<i2_npeak_found; ++i2_s){            // this proedure is used to sort f_peak_pos array
      f_peak_pos_sort[i2_s]=f_peak_pos[i2_index[i2_s]];
    }
    delete [] i2_index;

    const int32_t i2_nextrema = i2_npeak_found + 1;                    // numer of analysis zone extrema
    long* i4_extrema= new long[i2_nextrema];                       // position of analysis zone extrema
    extrema(i2_nextrema, i4_extrema, f_peak_pos_sort, f_samples);  // calculate extrema position
    
    for (int32_t i2_i=0; i2_i<i2_npeak_found; ++i2_i){

      const long i4_subzone_default_len = i4_extrema[i2_i+1] - i4_extrema[i2_i];
      const long i4_len_back  = 350;
      const long i4_len_front = 650;
      const long i4_len_max   = i4_len_back + i4_len_front; 

      /*const long i4_len_subzone=(i4_extrema[i2_i+1]-i4_extrema[i2_i] < i4_max_len ? i4_extrema[i2_i+1] -i4_extrema[i2_i] : (i4_max_len));
	const long i4_subzone_start = (i4_len_subzone < i4_max_len ? i4_extrema[i2_i] : (long) f_peak_pos_sort[i2_i] - 400);*/
      const long i4_subzone_len = ( i4_subzone_default_len <= i4_len_max ? i4_subzone_default_len : (static_cast <long>(f_peak_pos_sort[i2_i]) + i4_len_front < i4_extrema[i2_i+1] ? i4_len_max : i4_extrema[i2_i+1] - static_cast<long>(f_peak_pos_sort[i2_i]) + i4_len_back));
      const long i4_subzone_start = (i4_subzone_default_len <= i4_len_max ? i4_extrema[i2_i] : static_cast <long>((f_peak_pos_sort[i2_i]) - i4_len_back));
      const float* f_subzone = &f_samples[i4_subzone_start];
    
      const float f4_tmp_mean   = TMath::Mean(i4_subzone_len, f_subzone);
      const float f4_tmp_max_val= TMath::MaxElement(i4_subzone_len, f_subzone);
      const long  i4_tmp_max_pos= TMath::LocMax(i4_subzone_len, f_subzone);
      const float f4_tmp_HM     = (f4_tmp_max_val + f4_tmp_mean)*0.5;  // Half Maximum of the peak
      const float f4_tmp_FWHM   = static_cast<float> (v1731_neutron_analyzer::i4_over_thr(i4_subzone_len, f_subzone, f4_tmp_HM));// Full Width HM of the peak
      const float f4_tmp_rise   = (f4_tmp_FWHM < 80. ? (f4_tmp_FWHM/log(2.)) : 10.); // "tapullo"
     
            
      /* Peak fit procedure:
	 -get samples near peak 
	 -fit gauss/expo (A = (max_val - f_mean); mu = max_pos; sigma = fwhm/sqrt(2*log2) :10 ; tau = 40)
      */

      TH1F h_samples("peak_zone_histo", "Neutron Pulse Responce;time(ns);ADC bin;", i4_subzone_len, 0, 2*i4_subzone_len);
      for (long i4s=0; i4s<i4_subzone_len; ++i4s){
	h_samples.Fill(2*i4s, f_subzone[i4s]);
	h_samples.SetBinError(i4s+1, 1./sqrt(3.));
      }

      TF1 Responce("responce", v1731_neutron_analyzer::analitic_convolution, 600, 2*i4_subzone_len, 6);
      Responce.SetParameter(0, f4_tmp_mean);                  // base
      Responce.SetParameter(1, f4_tmp_max_val - f4_tmp_mean); // peak amplitude
      Responce.SetParameter(2, 2*i4_tmp_max_pos);             // peak time
      Responce.SetParameter(3, f4_tmp_rise);                  // peak rise time
      Responce.SetParameter(4, 400.);                         // RC coupling
      Responce.SetParameter(5, 40.);                          // peak fall time
      Responce.SetParLimits(0, 0., 255.);
      Responce.SetParLimits(1, 0., 300.);
      Responce.SetParLimits(3, 0., 80.);
      Responce.SetParLimits(4, 0., 10000.);
      Responce.SetParLimits(5, 0., 500.);
      h_samples.Fit("responce","LLQBMRN");

      const float f_tmp_area = Responce.GetParameter(1)*(sqrt(M_PI/2.)*TMath::Abs(Responce.GetParameter(3)) + TMath::Abs(Responce.GetParameter(5)));
      
      //const float f_tmp_area_l=0.5*Responce.GetParameter(1)*(sqrt(2.*M_PI)*TMath::Abs(Responce.GetParameter(3)));
      //const float f_tmp_area_r=0.5*Responce.GetParameter(1)*Responce.GetParameter(5);
      v1731_neutron_analyzer::pulse tmp_peak;
      tmp_peak.f8_time  = 2*(i4_smp_start + i4_subzone_start) + TMath::Abs(Responce.GetParameter(2)) - f_mu_time;
      tmp_peak.f4_peak  = Responce.GetParameter(1);
      tmp_peak.f4_rise  = TMath::Abs(Responce.GetParameter(3));
      tmp_peak.f4_fall  = TMath::Abs(Responce.GetParameter(5));
      tmp_peak.f4_charge= f_tmp_area;
      peaks.push_back(tmp_peak);
    }
    delete [] f_peak_pos_sort;
    delete [] i4_extrema;
  }
  
}


void v1731_neutron_analyzer::peak_zone::extrema(int32_t i2_nextrema, long* i4_extrema, Float_t* f_in_peak_pos, float* f_in_samples){
  i4_extrema[0] = 0;
  i4_extrema[i2_nextrema-1] = i4_size-1;
  for (int32_t i2_e=0; i2_e<i2_npeak_found-1; ++i2_e){
    const float f_x1 = f_in_peak_pos[i2_e];        // position of e-th peak
    const float f_x2 = f_in_peak_pos[i2_e+1];      // position of (e+1)-th peak
    const float f_w1 = f_in_samples[static_cast<long> (f_x1)];  // weight of e-th peak
    const float f_w2 = f_in_samples[static_cast<long> (f_x2)];  // weight of (e+1)-th peak
    i4_extrema[i2_e+1]= static_cast<long> ((f_w1*f_x2 + f_w2*f_x1)/(f_w1+f_w2)); //extrema is a "inverted" weighted mean 
  }
}



/************* peak ***************/


v1731_neutron_analyzer::peak::peak(): bx_named("v1731_neutron_peak"){
}

void v1731_neutron_analyzer::peak::clear(){
  peaks.clear();
}

void v1731_neutron_analyzer::peak::push(const v1731_neutron_analyzer::peak_zone& pz, unsigned u4_peak){
  peaks.push_back(pz.get_peak(u4_peak));
}

void v1731_neutron_analyzer::peak::push(const v1731_neutron_analyzer::pulse& in_pulses){
  peaks.push_back(in_pulses);
}

void v1731_neutron_analyzer::peak::show(){
  bx_message &msg = get_message(bx_message::info);
  for (std::vector<v1731_neutron_analyzer::pulse>::const_iterator cit=peaks.begin(); cit!=peaks.end(); ++cit){
    msg<<" Time: "<<(*cit).f8_time<<" ns\t";
    msg<<" Charge: "<<(*cit).f4_charge<<" ns*bin\t";
    msg<<" Ampl: "<<(*cit).f4_peak<<" bin\t";
    msg<<" Rise: "<<(*cit).f4_rise<<" ns\t";
    msg<<" Fall: "<<(*cit).f4_fall<<" ns \n";
  }
  msg<<dispatch;
}

///////////////////////////////////////////////////////


/*
 * $Log: v1731_neutron_analyzer.cc,v $
 * Revision 1.10  2009/11/20 10:52:26  davini
 * removed MINOS call
 *
 * Revision 1.9  2008-12-15 13:50:42  davini
 * quiet
 *
 * Revision 1.8  2008-12-15 13:35:38  davini
 * bug fized: f4_muon_time always initialized
 *
 * Revision 1.7  2008-09-14 14:45:17  davini
 * removed an error about charge computation (hoping it's not useless now...)
 *
 * Revision 1.6  2008-08-18 16:04:55  davini
 * removed warning (thanks again to A.R.)
 *
 * Revision 1.5  2008-08-18 15:57:35  davini
 * removed warnings (thanks to A.R.)
 *
 * Revision 1.4  2008-08-14 18:42:35  davini
 * neutron search begins 14 us after muon
 *
 * Revision 1.3  2008-08-06 15:00:36  davini
 * added fall time;
 *
 *
 *
 */
