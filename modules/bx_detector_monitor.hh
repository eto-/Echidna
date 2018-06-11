/* BOREXINO Reconstruction program
 *
 * Author: Livia Ludhova <livia.ludhova@mi.infn.it>
 * Maintainer: Livia Ludhova <livia.ludhova@mi.infn.it>
 *
 * 1) online control crate, fe-board and laben-board occupancy distribution after each n_events (from echidna.cfg) 
 */

#ifndef _BX_DETECTOR_MONITOR_HH
#define _BX_DETECTOR_MONITOR_HH

#include <vector>
#include "bx_base_module.hh"
#include "TVector3.h"
#include "TF1.h"
#include "TF2.h"
#include "TH1F.h"
#include "TH2F.h"
#include <list>
#include "constants.hh"
#include "bx_detector.hh"

class bx_echidna_event;

class bx_detector_monitor: public bx_base_module {
  public:
    bx_detector_monitor ();
    virtual ~bx_detector_monitor () {}
    
    // Operations for the framework 
    virtual void begin ();
    virtual bx_echidna_event* doit (bx_echidna_event *ev);
    virtual void end ();

  private:
     
      // occupancy
    int32_t total_triggers;
    int32_t n_triggers [4]; 
  
    int32_t check_after_precalib_done;

    int32_t calibration_mode;
    int32_t disabling_lg;

    int32_t n[5];
    int32_t n_on[5];

    int32_t start_lg_off;
    int32_t start_lg_on;
    int32_t end_lg_off;
    int32_t end_lg_on;
    int32_t start_lg_off_charge;
    int32_t start_lg_on_charge;
    int32_t end_lg_off_charge;
    int32_t end_lg_on_charge;

    std::vector<int32_t> times_found_noisy[5][4];
    std::vector<int32_t> times_found_empty[5][4];

      //precalibration check
    std::vector<float> nhits_in_precalib;
      
    //cumulative_off (in any trigger)
    std::vector<int32_t> cummulative_off_lg;
    std::vector<int32_t> cummulative_charge_off_lg;

    std::vector<int32_t> latest_off_lg_p;
    std::vector<int32_t> latest_off_lg_HV;
    std::vector<double> cummulative_occupancy_p;
    std::vector<double> cummulative_lbnb_occupancy_p;
    std::vector<double> cummulative_occupancy_HV;
    std::vector<double> cummulative_hvb_occupancy_HV;

      //how many time the lg changed its status to off (on) in pulser and HV triggers
    std::vector<int32_t> lg_Nx_to_off_p;
    std::vector<int32_t> lg_Nx_to_on_p;
    std::vector<int32_t> lg_Nx_to_off_HV;
    std::vector<int32_t> lg_Nx_to_on_HV;

      //the event number of the LATEST change
    std::vector<int32_t> ev_lg_to_off_p;
    std::vector<int32_t> ev_lg_to_on_p;
    std::vector<int32_t> ev_lg_to_off_HV;
    std::vector<int32_t> ev_lg_to_on_HV;

      // for neutrino trigger rate and nhits rate
    uint32_t n_very_first_event_tsec;
    int32_t n_trgrate_cycle;
    int32_t n_events_in_this_cycle;
    double n_time_last_previous_cycle;
    int32_t n_nhits_in_cycle;

    std::list<double> tlist;
    std::list<int32_t> neutrino_nhits_list;
    
      // for random nhits rate
    uint32_t r_very_first_event_tsec;
    int32_t r_trgrate_cycle;
    int32_t r_events_in_this_cycle;
    double r_time_last_previous_cycle;
    int32_t r_nhits_in_cycle;

    std::list<double> random_tlist;
    std::list<int32_t> random_nhits_list;

     
    //pulser
    double nhits_pulser_sum;
               
      //histo for electronics
      
    // 0 = pulser, 1 = random+neutrino+laser 0= off 1= on
    TH1F* h_lg_Nx_changed[2][2];

    TH1F* dec_occupancy[5][4];
    TH1F* neutrino_trigger_rate;    

    TH2F* neutrino_charge_vs_lg;    
    TH2F* charge_pulser_vs_lg;    
  
    TH1F* neutrino_nhits_rate;      
    TH1F* random_nhits_rate;    
    TH1F* neutrino_nhits_per_event;      
    TH1F* random_nhits_per_event;    
    TH1F* pulser_nhits;    


    int32_t pulser_bins;
    int32_t max_events;
  
    enum module_type{
      crate = 0,
      hvb = 1,
      feb = 2,
      lbnb = 3,
      channel = 4,
    };

  enum trigger_type{
      pulser = 0,
      random = 1,
      neutrino = 2,
      laser = 3,
   };

  std::cmap <int32_t, std::string> module_type_map;
  std::cmap <int32_t, std::string> trigger_type_map;

};
#endif
