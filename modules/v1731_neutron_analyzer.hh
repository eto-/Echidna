/* BOREXINO Reconstruction program
 *
 * Author: Stefano Davini <stefano.davini@ge.infn.it>
 * Maintainer: Stefano Davini <stefano.davini@ge.infn.it>
 *
 * $Id:
 *
 * Definiton of namespace v1731_neutron_analyzer, 
 * whose class neutron_analyzer is used in bx_v1731sys 
 * to analyze v1731event 
 *
 * 
 */

#ifndef _V1731_NEUTRON_ANALYZER_HH
#define _V1731_NEUTRON_ANALYZER_HH

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include <TH1F.h>
#include <TF1.h>
#include <TMath.h>
#include <TSpectrum.h>
#include "v1731_event.hh" 


namespace v1731_neutron_analyzer{
  long i4_over_thr(long i4_size, const float* f_in_samples, float f_thr);
  float f_area_over_thr(long i4_size, const float* f_in_samples, float f_thr); 
  Double_t analitic_convolution(Double_t* x, Double_t* par);

  struct pulse{
    float  f4_charge;  // charge of orginal signa (before CR) - in (ADC bins*ns)
    double f8_time;    // neutron peak time - muon edge time in (ns)
    float  f4_peak;    // peak amplitude of signal (before CR) in (ADC bins)
    float  f4_rise;    // rise time of signal (before CR) in (ns)
    float  f4_fall;    // fall time of signal (before CR) in (ns)
    float  f4_base;
  };

  class neutron_analyzer;
  class peak_zone;
  class peak;

  using v1731_event::v1731event;
};

/******************** peak ****************************/

class v1731_neutron_analyzer::peak: public bx_named{
public:
  peak();
  virtual ~peak() {}
  void clear();
  void push(const v1731_neutron_analyzer::peak_zone&, unsigned u4_peak);
  void push(const v1731_neutron_analyzer::pulse&);
  void show();
  unsigned get_number () const {return peaks.size();}
  float get_time(unsigned u4)  const {return (u4<peaks.size() ? peaks[u4].f8_time  : -2);}
  float get_peak(unsigned u4)  const {return (u4<peaks.size() ? peaks[u4].f4_peak : -2);}
  float get_rise(unsigned u4)  const {return (u4<peaks.size() ? peaks[u4].f4_rise : -2);}
  float get_charge(unsigned u4)const {return (u4<peaks.size() ? peaks[u4].f4_charge : -2);}
  float get_fall(unsigned u4)  const {return (u4<peaks.size() ? peaks[u4].f4_fall : -2);}

private:
  std::vector <v1731_neutron_analyzer::pulse> peaks;
};


/********** neutron_analyzer**********/

class v1731_neutron_analyzer::neutron_analyzer{
public:
  neutron_analyzer(const v1731event&, bool b_always_muon = false);
  ~neutron_analyzer() {}
  bool get_is_muon()     const {return b_is_muon;}
  float get_muon_time()  const {return f4_muon_time;}
  bool get_has_neutron()   const {return (analog_neutrons.get_number()>0);}
  unsigned get_n_neutrons()const {return analog_neutrons.get_number();}
  float get_time(unsigned u4)  const {return analog_neutrons.get_time(u4);}
  float get_peak(unsigned u4)  const {return analog_neutrons.get_peak(u4);}
  float get_charge(unsigned u4)const {return analog_neutrons.get_charge(u4);}
  float get_rise(unsigned u4) const {return analog_neutrons.get_rise(u4);}
  float get_fall(unsigned u4) const {return analog_neutrons.get_fall(u4);}
private: 
  v1731_neutron_analyzer::peak analog_neutrons;
  bool b_is_muon;
  float f4_muon_time;
  int32_t i2_n_pulses;
  bool is_muon(const v1731event&, bool b_always_muon = false);
  void analog_sum_analyze(const v1731event&);
  void zle_analog_sum_analyze(const v1731event&);
};


/************************** peak_zone *********************/

class v1731_neutron_analyzer::peak_zone{
public:
  peak_zone(long i4_in_size, float* f_in_samples, float* f_in_samples_clean, long i4_in_sample_start, float f_mu_time = 0);
  ~peak_zone() {}
  unsigned get_npeak_found() const {return peaks.size();}
  const v1731_neutron_analyzer::pulse& get_peak(unsigned u4) const {return peaks[u4];}

private:
  long i4_size;
  Int_t i2_npeak_found;  
  std::vector <v1731_neutron_analyzer::pulse> peaks;
  void extrema(int32_t i2_nextrema, long* i4_extrema, Float_t* f_peak_pos, float* f_in_samples);  // returns extrema of peak zones 



};

////////////////////////////////////

#endif

/*
 * $Log: v1731_neutron_analyzer.hh,v $
 * Revision 1.3  2008/08/06 15:00:36  davini
 * added fall time;
 *
 *
 */
