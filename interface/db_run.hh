/* BOREXINO Reconstruction program
 *
 * Author: Alessandro Razeto <Alessandro.Razeto@ge.infn.it>
 * Maintainers: Alessandro Razeto <Alessandro.Razeto@ge.infn.it>,
 *
 * $Id: db_run.hh,v 1.58 2015/07/24 12:37:24 ilia.drachnev Exp $
 *
 * The database interface for run informations
 *
 */
#ifndef _BD_RUN_H
#define _BD_RUN_H
#if defined(__ROOTCLING__) || defined(__CINT__)
#define _ROOT_CINT_
#endif
#ifndef CYCLE_NUMBER
#define CYCLE_NUMBER 18
#endif
#include "db_acl.hh"
#include "db_acl.hh"

#include <TObject.h>

#include <string>
#include <map>
#include <vector>
#include <time.h>

class db_run: public db_acl, public TObject {
  public:
    enum run_type {
      normal,
      calibration,
      laser,
      laser394,
      laser355,
      laser266,
      random,
      pulser,
      custom,
      no_file,
      test,
      source,
      water,
      nofile,
    };

    enum disabled_type {
      charge = 1,
      timing = 10,
    };

    db_run() : db_acl(), TObject() {};

      // Getters for the run daq data
    int get_number	 	() const { return i_run_number; }
    int get_profile_id	 	() const { return i_profile_id; }
    int get_calib_profile 	() const { return i_calib_profile; }
    run_type get_type	 	() const { return i_run_type; }
    int get_max_time	 	() const { return i_max_time; }
    int get_max_events 		() const { return i_max_events; }
    int get_stop_reason		() const { return i_stop_reason; }
    int get_number_of_files 	() const { return i_number_of_files; }
    time_t get_start_time	() const { return u4_start_time; }
    time_t get_stop_time	() const { return u4_stop_time; }
    int get_duration 		() const { return i_duration; }
    int get_events 		() const { return i_events; }

      // S/Getters for the precalib data: a check method is also supplied to
      // check if the data is available (ie the precalib for the required data
      // has already run or the data are already stored on db).
    bool is_precalib_present 		() const { return b_precalib_present; }
    bool is_precalib_quality_present 	() const { return b_precalib_quality_present; }
    
    int    get_laben_precalib_low_bin  	 (int lg) const { return map_get(lg, laben_precalib_low_bin_v,"laben_precalib_low_bin_v"); }
    int    get_laben_precalib_high_bin 	 (int lg) const { return map_get(lg, laben_precalib_high_bin_v,"laben_precalib_high_bin_v"); }
    bool   get_laben_precalib_rising_on_even (int lg) const { return map_get(lg, laben_precalib_rising_on_even_v,"laben_precalib_rising_on_even_v"); }
    int    get_laben_precalib_gray_shift (int lg) const { return map_get(lg, laben_precalib_gray_shift_v,"laben_precalib_gray_shift_v"); }
    float  get_laben_precalib_delta80 	 (int lg) const { return map_get(lg, laben_precalib_delta80_v,"laben_precalib_delta80_v"); }
    uint32_t get_muon_precalib_pulse_time () const { return u4_muon_precalib_pulse_time; }
    double get_muon_precalib_pedestal	 (int lg) const { return map_get(lg, muon_precalib_pedestal_v,"muon_precalib_pedestal_v"); } 
    double get_muon_precalib_pedsigma	 (int lg) const { return map_get(lg, muon_precalib_pedsigma_v,"muon_precalib_pedsigma_v"); }
    float  get_laben_precalib_mean_time  ()       const { return  f_laben_precalib_mean_time; }
    float  get_laben_precalib_sigma_time ()      const { return f_laben_precalib_sigma_time; }
    const std::vector<int>& get_laben_precalib_bad_channels  () const { return laben_precalib_bad_channels_v; }
    const std::vector<int>& get_laben_precalib_off_channels  () const { return laben_precalib_off_channels_v; }
    
    bool check_laben_precalib_low_bin	(int lg) const { return map_check(lg, laben_precalib_low_bin_v); }
    bool check_laben_precalib_delta80	(int lg) const { return map_check(lg, laben_precalib_delta80_v); }
    bool check_laben_precalib_high_bin	(int lg) const { return map_check(lg, laben_precalib_high_bin_v); }
    bool check_laben_precalib_rising_on_even (int lg) const { return map_check(lg, laben_precalib_rising_on_even_v); }
    bool check_laben_precalib_gray_shift(int lg) const { return map_check(lg, laben_precalib_gray_shift_v); }
    bool check_muon_precalib_pedestal	(int lg) const { return map_check(lg, muon_precalib_pedestal_v); } 
    bool check_muon_precalib_pedsigma	(int lg) const { return map_check(lg, muon_precalib_pedsigma_v); }

     //effective quantum efficiensies
    double get_laben_qe (int lg) const { return map_get(lg,laben_qe_v ,"effectiveQE"); }
    double get_laben_qe_nocorr (int lg) const { return map_get(lg,laben_qe_nocorr_v ,"effectiveQE_no_corr"); }
    double get_laben_qe_dark_noise (int lg) const { return map_get(lg,laben_qe_dark_noise_v ,"effectiveQE_dark_noise"); }



#if !defined(_ECHIDNA_ROOTLIB_) && !defined(_ROOT_CINT_)
    void set_laben_precalib_low_bin	(int lg, int count, const bx_named* obj) { if (check_acl ("set_laben_precalib_low_bin", obj)) laben_precalib_low_bin_v[lg] = count; }
    void set_laben_precalib_high_bin	(int lg, int count, const bx_named* obj) { if (check_acl ("set_laben_precalib_high_bin", obj)) laben_precalib_high_bin_v[lg] = count; }
    void set_laben_precalib_gray_shift	(int lg, int shift, const bx_named* obj) { if (check_acl ("set_laben_precalib_gray_shift", obj)) laben_precalib_gray_shift_v[lg] = shift; }
    void set_laben_precalib_rising_on_even	(int lg, bool t, const bx_named* obj)  { if (check_acl ("set_laben_precalib_rising_on_even", obj)) laben_precalib_rising_on_even_v[lg] = t; }
    void set_laben_precalib_delta80	(int lg, float dt, const bx_named* obj)  { if (check_acl ("set_laben_precalib_delta80", obj)) laben_precalib_delta80_v[lg] = dt; }
    void set_muon_precalib_pulse_time	(uint32_t time, const bx_named* obj) { if (check_acl ("set_muon_precalib_pulse_time", obj)) u4_muon_precalib_pulse_time = time; }
    void set_muon_precalib_pedestal	(int lg, double pedestal, const bx_named* obj) { if (check_acl ("set_muon_precalib_pedestal", obj)) muon_precalib_pedestal_v[lg] = pedestal; }
    void set_muon_precalib_pedsigma	(int lg, double pedsigma, const bx_named* obj) { if (check_acl ("set_muon_precalib_pedsigma", obj)) muon_precalib_pedsigma_v[lg] = pedsigma; }
    void set_laben_precalib_mean_time   (float laben_precalib_mean_time, const bx_named* obj)	   { if (check_acl ("set_laben_precalib_mean_time", obj)) f_laben_precalib_mean_time = laben_precalib_mean_time; }  		
    void set_laben_precalib_sigma_time  (float laben_precalib_sigma_time, const bx_named* obj)	   { if (check_acl ("set_laben_precalib_sigma_time", obj)) f_laben_precalib_sigma_time = laben_precalib_sigma_time; }		
    void set_laben_precalib_bad_channels  (std::vector<int>& laben_precalib_bad_channels, const bx_named* obj) { if (check_acl ("set_laben_precalib_bad_channels", obj)) laben_precalib_bad_channels_v = laben_precalib_bad_channels; } 
    void set_laben_precalib_off_channels  (std::vector<int>& laben_precalib_off_channels, const bx_named* obj) { if (check_acl ("set_laben_precalib_off_channels", obj)) laben_precalib_off_channels_v = laben_precalib_off_channels; } 
    void reset_precalibrations		(const bx_named* obj);
    
    void write_laben_precalib 		(bool y, const bx_named* obj) { if (check_acl ("write_laben_precalib", obj)) b_write_laben_precalib = y; }
    void write_laben_precalib_quality	(bool y, const bx_named* obj) { if (check_acl ("write_laben_precalib_quality", obj)) b_write_laben_precalib_quality = y; }
    void write_muon_precalib 		(bool y, const bx_named* obj) { if (check_acl ("write_muon_precalib", obj)) b_write_muon_precalib = y; }
    void write_disabled_channels 	(bool y, const bx_named* obj) { if (check_acl ("write_disabled_channels", obj)) b_write_disabled_channels = y; }
#endif

      // S/Getters for some calibration data
    bool is_laben_calib_tt1_present          () const { return b_laben_calib_tt1_present         ; }
    bool is_laben_laser_present              () const { return b_laben_laser_present             ; }
    bool is_muon_laser_present 	             () const { return b_muon_laser_present              ; }
    bool is_trigger_parameters_present       () const { return b_trigger_parameters_present      ; }
    bool is_laben_dark_rates_present	     () const { return b_laben_dark_rates_present        ; }
    bool is_muon_dark_rates_present	     () const { return b_muon_dark_rates_present         ; }
    bool is_laben_electronic_channel_present () const { return b_laben_electronic_channel_present; }
    bool is_muon_electronic_channel_present  () const { return b_muon_electronic_channel_present ; }
    bool is_muon_alignment_present 	     () const { return b_muon_alignment_present         ; }
    bool is_disabled_channels_present 	     () const { return b_disabled_channels_present         ; }
    bool is_laben_dark_rates_parametrisation_present   () const { return b_laben_dark_rates_parametrisation_present; }
    
    float get_laben_time_offset		(int lg) const { return map_get(lg, laben_time_offset_v, "laben_time_offset_v"); }
    float get_laben_time_sigma		(int lg) const { return map_get(lg, laben_time_sigma_v, "laben_time_sigma_v"); }
    float get_laben_charge_peak		(int lg) const { return map_get(lg, laben_charge_peak_v, "laben_charge_peak_v"); }
    float get_laben_charge_sigma	(int lg) const { return map_get(lg, laben_charge_sigma_v, "laben_charge_sigma_v"); }
    float get_laben_charge_tt1_peak	(int lg) const { return map_get(lg, laben_charge_tt1_peak_v, "laben_charge_tt1_peak_v"); }
    float get_laben_charge_tt1_sigma	(int lg) const { return map_get(lg, laben_charge_tt1_sigma_v, "laben_charge_tt1_sigma_v"); }
    float get_laben_charge_tt1_mean	(int lg) const { return map_get(lg, laben_charge_tt1_mean_v, "laben_charge_tt1_mean_v"); }
    float get_laben_charge_tt1_rms	(int lg) const { return map_get(lg, laben_charge_tt1_rms_v, "laben_charge_tt1_rms_v"); }
    float get_laben_charge_tt1_p0	(int lg) const { return map_get(lg, laben_charge_tt1_p0_v, "laben_charge_tt1_p0_v"); }

    float get_muon_time_offset		(int lg) const { return map_get(lg, muon_time_offset_v, "muon_time_offset_v"); }
    float get_muon_time_sigma		(int lg) const { return map_get(lg, muon_time_sigma_v, "muon_time_sigma_v"); }
    float get_muon_charge_peak		(int lg) const { return map_get(lg, muon_charge_peak_v, "muon_charge_peak_v"); }
    float get_muon_charge_sigma	        (int lg) const { return map_get(lg, muon_charge_sigma_v, "muon_charge_sigma_v"); }

    float get_laben_gate_width		() const { return f_laben_gate_width; }
    float get_laben_gate_start		() const { return f_laben_gate_start; }
    float get_laben_laser_offset        () const { return f_laben_laser_offset; }
    float get_laben_pulser_offset	() const { return f_laben_pulser_offset; }
    float get_laben_cluster_offset      () const { return f_laben_cluster_offset; }

    float get_laben_mean_dark_noise     () const { return f_laben_mean_dark_noise; }
    float get_laben_mean_dark_sigma     () const { return f_laben_mean_dark_sigma; }
    int   get_laben_dead_cone           () const { return i_laben_dead_cone; }
    int   get_laben_dead_no_cone        () const { return i_laben_dead_no_cone; }
    int   get_laben_hot_cone            () const { return i_laben_hot_cone; }
    int   get_laben_hot_no_cone         () const { return i_laben_hot_no_cone;}

    float get_laben_dark_noise          (int lg) const  { return  map_get(lg, laben_dark_noise_v, "laben_dark_noise_v");}
    float get_laben_dark_sigma          (int lg) const  { return  map_get(lg, laben_dark_sigma_v, "laben_dark_sigma_v");}
    const std::string& get_laben_pmt_status (int lg) const  { return  map_get(lg, laben_pmt_status_v, "laben_pmt_status_v");}

    const std::map<int, std::map<int, db_run::disabled_type> >& get_disabled_channels () const { return disabled_channels_v; }

    float get_muon_mean_dark_noise     () const { return f_muon_mean_dark_noise; }
    float get_muon_mean_dark_sigma     () const { return f_muon_mean_dark_sigma; }
    int   get_muon_dead                () const { return i_muon_dead; }
    int   get_muon_hot                 () const { return i_muon_hot; }

    float get_muon_dark_noise          (int lg) const  { return  map_get(lg, muon_dark_noise_v, "muon_dark_noise_v");}
    float get_muon_dark_sigma          (int lg) const  { return  map_get(lg, muon_dark_sigma_v, "muon_dark_sigma_v");}
    const std::string& get_muon_pmt_status (int lg) const  { return  map_get(lg, muon_pmt_status_v, "muon_pmt_status_v");}

    int    get_muon_nevents             () const { return i_muon_nevents;         }
    time_t get_muon_time                () const { return t_muon_time;            }
    int    get_muon_up_nevents          () const { return i_muon_up_nevents;      }
    time_t get_muon_up_time             () const { return t_muon_up_time;         }
    int    get_muon_aligned_nevents     () const { return i_muon_aligned_nevents; }
    time_t get_muon_aligned_time        () const { return t_muon_aligned_time;    }

    //InnerDetectorDarkRateParametrisation
    float get_mean_dark_rate_per_used_pmt        () const { return f_mean_dark_rate_per_used_pmt      ;}
    float get_error_mean_dark_rate_per_used_pmt  () const { return f_error_mean_dark_rate_per_used_pmt;}
    float get_mu_win1                            () const { return f_mu_win1                          ;}
    float get_mu_win2                            () const { return f_mu_win2                          ;}
    float get_pois_const_win1                    () const { return f_pois_const_win1                  ;}            
    float get_error_pois_const_win1              () const { return f_error_pois_const_win1            ;}            
    float get_mu_fit_win1                        () const { return f_mu_fit_win1                      ;}                
    float get_error_mu_fit_win1                  () const { return f_error_mu_fit_win1                ;}                
    float get_exp_const_win1                     () const { return f_exp_const_win1                   ;}
    float get_error_exp_const_win1               () const { return f_error_exp_const_win1             ;}
    float get_tau_win1                           () const { return f_tau_win1                         ;}
    float get_error_tau_win1                     () const { return f_error_tau_win1                   ;}
    float get_pois_const_win2                    () const { return f_pois_const_win2                  ;}            
    float get_error_pois_const_win2              () const { return f_error_pois_const_win2            ;}            
    float get_mu_fit_win2                        () const { return f_mu_fit_win2                      ;}                
    float get_error_mu_fit_win2                  () const { return f_error_mu_fit_win2                ;}                
    float get_exp_const_win2                     () const { return f_exp_const_win2                   ;}
    float get_error_exp_const_win2               () const { return f_error_exp_const_win2             ;}
    float get_tau_win2                           () const { return f_tau_win2                         ;}
    float get_error_tau_win2                     () const { return f_error_tau_win2                   ;}

    /*const std::vector<std::string>& get_laben_charge_base_status  (int lg) const  { return  map_get(lg, laben_charge_base_status_v);}
    const std::vector<std::string>& get_laben_charge_peak_status  (int lg) const  { return  map_get(lg, laben_charge_peak_status_v);}
    const std::vector<std::string>& get_laben_charge_status       (int lg) const  { return  map_get(lg, laben_charge_status_v);}
    const std::vector<std::string>& get_laben_timing_status       (int lg) const  { return  map_get(lg, laben_timing_status_v);}
    const std::vector<std::string>& get_laben_multiplicity        (int lg) const  { return  map_get(lg, laben_multiplicity_v);}*/
    const std::vector<std::string>& get_laben_charge_base_status  (int lg) const  { return  laben_charge_base_status_v[lg];}
    const std::vector<std::string>& get_laben_charge_peak_status  (int lg) const  { return  laben_charge_peak_status_v[lg];}
    const std::vector<std::string>& get_laben_charge_status       (int lg) const  { return  laben_charge_status_v[lg];}
    const std::vector<std::string>& get_laben_timing_status       (int lg) const  { return  laben_timing_status_v[lg];}
    const std::vector<std::string>& get_laben_multiplicity        (int lg) const  { return  laben_multiplicity_v[lg];}
    const std::vector<std::string>& get_muon_multiplicity         (int lg) const  { return  muon_multiplicity_v[lg-3001];} // offset to be removed if reverted back to map container
    
    bool is_pmt_disconnected	       (int lg) const  { return disconnected_pmts_v.find (lg) != disconnected_pmts_v.end (); }

#if !defined(_ECHIDNA_ROOTLIB_) && !defined(_ROOT_CINT_)
    void set_laben_time_offset		(int lg, float laben_time_offset, const bx_named* obj) { if (check_acl ("set_laben_time_offset", obj)) laben_time_offset_v[lg] = laben_time_offset; }
    void set_laben_time_sigma	        (int lg, float laben_time_sigma, const bx_named* obj) { if (check_acl ("set_laben_time_sigma", obj)) laben_time_sigma_v[lg] = laben_time_sigma; }
    void set_laben_charge_peak		(int lg, float laben_charge_peak, const bx_named* obj) { if (check_acl ("set_laben_charge_peak", obj)) laben_charge_peak_v[lg] = laben_charge_peak; }
    void set_laben_charge_sigma		(int lg, float laben_charge_sigma, const bx_named* obj) { if (check_acl ("set_laben_charge_sigma", obj)) laben_charge_sigma_v[lg] = laben_charge_sigma; }
    void set_laben_charge_tt1_peak	(int lg, float laben_charge_tt1_peak, const bx_named* obj) { if (check_acl ("set_laben_charge_tt1_peak", obj)) laben_charge_tt1_peak_v[lg] = laben_charge_tt1_peak; }
    void set_laben_charge_tt1_sigma	(int lg, float laben_charge_tt1_sigma, const bx_named* obj) { if (check_acl ("set_laben_charge_tt1_sigma", obj)) laben_charge_tt1_sigma_v[lg] = laben_charge_tt1_sigma; }
    void set_laben_charge_tt1_mean	(int lg, float laben_charge_tt1_mean, const bx_named* obj) { if (check_acl ("set_laben_charge_tt1_mean", obj)) laben_charge_tt1_mean_v[lg] = laben_charge_tt1_mean; }
    void set_laben_charge_tt1_rms	(int lg, float laben_charge_tt1_rms, const bx_named* obj) { if (check_acl ("set_laben_charge_tt1_rms", obj)) laben_charge_tt1_rms_v[lg] = laben_charge_tt1_rms; }
    void set_laben_charge_tt1_p0	(int lg, float laben_charge_tt1_p0, const bx_named* obj) { if (check_acl ("set_laben_charge_tt1_p0", obj)) laben_charge_tt1_p0_v[lg] = laben_charge_tt1_p0; }
    void set_muon_time_offset		(int lg, float muon_time_offset, const bx_named* obj) { if (check_acl ("set_muon_time_offset", obj)) muon_time_offset_v[lg] = muon_time_offset; }
    void set_muon_time_sigma	        (int lg, float muon_time_sigma, const bx_named* obj) { if (check_acl ("set_muon_time_sigma", obj)) muon_time_sigma_v[lg] = muon_time_sigma; }
    void set_muon_charge_peak		(int lg, float muon_charge_peak, const bx_named* obj) { if (check_acl ("set_muon_charge_peak", obj)) muon_charge_peak_v[lg] = muon_charge_peak; }
    void set_muon_charge_sigma		(int lg, float muon_charge_sigma, const bx_named* obj) { if (check_acl ("set_muon_charge_sigma", obj)) muon_charge_sigma_v[lg] = muon_charge_sigma; }
    void set_laben_gate_width           (float laben_gate_width, const bx_named* obj) { if (check_acl ("set_laben_gate_width", obj)) f_laben_gate_width = laben_gate_width; }
    void set_laben_gate_start           (float laben_gate_start, const bx_named* obj) { if (check_acl ("set_laben_gate_start", obj)) f_laben_gate_start = laben_gate_start; }
    void set_laben_laser_offset         (float laben_laser_offset, const bx_named* obj) { if (check_acl ("set_laben_laser_offset", obj)) f_laben_laser_offset =  laben_laser_offset; }
    void set_laben_pulser_offset        (float laben_pulser_offset, const bx_named* obj) { if (check_acl ("set_laben_pulser_offset", obj)) f_laben_pulser_offset = laben_pulser_offset; }
    void set_laben_cluster_offset       (float laben_cluster_offset, const bx_named* obj) { if (check_acl ("set_laben_cluster_offset", obj)) f_laben_cluster_offset = laben_cluster_offset; }
 
    void set_laben_mean_dark_noise      (float dark_noise, const bx_named* obj) { if (check_acl ("set_laben_mean_dark_noise", obj)) f_laben_mean_dark_noise = dark_noise; }
    void set_laben_mean_dark_sigma      (float dark_sigma, const bx_named* obj) { if (check_acl ("set_laben_mean_dark_sigma", obj)) f_laben_mean_dark_sigma = dark_sigma; }
    void set_laben_dead_cone            (int dead_cone, const bx_named* obj)    { if (check_acl ("set_laben_dead_cone", obj)) i_laben_dead_cone = dead_cone;}
    void set_laben_dead_no_cone         (int dead_no_cone, const bx_named* obj) { if (check_acl ("set_laben_dead_no_cone", obj)) i_laben_dead_no_cone = dead_no_cone; }
    void set_laben_hot_cone             (int hot_cone, const bx_named* obj)     { if (check_acl ("set_laben_hot_cone", obj)) i_laben_hot_cone = hot_cone; }
    void set_laben_hot_no_cone          (int hot_no_cone, const bx_named* obj)  { if (check_acl ("set_laben_hot_no_cone", obj)) i_laben_hot_no_cone = hot_no_cone; }
    void set_laben_dark_noise           (int lg, float dark_noise, const bx_named* obj)  { if (check_acl ("set_laben_dark_noise", obj)) laben_dark_noise_v[lg] = dark_noise;}
    void set_laben_dark_sigma           (int lg, float dark_sigma, const bx_named* obj)  { if (check_acl ("set_laben_dark_sigma", obj)) laben_dark_sigma_v[lg] = dark_sigma;}
    void set_laben_pmt_status           (int lg, const std::string& pmt_status, const bx_named* obj)    { if (check_acl ("set_laben_pmt_status", obj)) laben_pmt_status_v[lg] = pmt_status;}    

    void add_disabled_channel		(int lg, int evnum, disabled_type type, const bx_named* obj) { if (check_acl ("add_disabled_channel", obj) && disabled_channels_v[evnum][lg] != timing) disabled_channels_v[evnum][lg] = type; }

    void set_muon_mean_dark_noise      (float dark_noise, const bx_named* obj) { if (check_acl ("set_muon_mean_dark_noise", obj)) f_muon_mean_dark_noise = dark_noise; }
    void set_muon_mean_dark_sigma      (float dark_sigma, const bx_named* obj) { if (check_acl ("set_muon_mean_dark_sigma", obj)) f_muon_mean_dark_sigma = dark_sigma; }
    void set_muon_dead                 (int dead, const bx_named* obj)    { if (check_acl ("set_muon_dead", obj)) i_muon_dead = dead;}
    void set_muon_hot                  (int hot, const bx_named* obj)     { if (check_acl ("set_muon_hot", obj)) i_muon_hot = hot; }
    void set_muon_dark_noise           (int lg, float dark_noise, const bx_named* obj)  { if (check_acl ("set_muon_dark_noise", obj)) muon_dark_noise_v[lg] = dark_noise;}
    void set_muon_dark_sigma           (int lg, float dark_sigma, const bx_named* obj)  { if (check_acl ("set_muon_dark_sigma", obj)) muon_dark_sigma_v[lg] = dark_sigma;}
    void set_muon_pmt_status           (int lg, const std::string& pmt_status, const bx_named* obj)    { if (check_acl ("set_muon_pmt_status", obj)) muon_pmt_status_v[lg] = pmt_status;}    

    void set_laben_charge_base_status   (int lg, const std::vector<std::string>& charge_base_status, const bx_named* obj) { if (check_acl ("set_laben_charge_base_status", obj)) laben_charge_base_status_v[lg] = charge_base_status;} 
    void set_laben_charge_peak_status   (int lg, const std::vector<std::string>& charge_peak_status, const bx_named* obj) { if (check_acl ("set_laben_charge_peak_status", obj)) laben_charge_peak_status_v[lg] = charge_peak_status;} 
    void set_laben_charge_status        (int lg, const std::vector<std::string>& charge_status, const bx_named* obj) { if (check_acl ("set_laben_charge_status", obj)) laben_charge_status_v[lg] = charge_status;}     
    void set_laben_timing_status        (int lg, const std::vector<std::string>& timing_status, const bx_named* obj) { if (check_acl ("set_laben_timing_status", obj)) laben_timing_status_v[lg] = timing_status;} 
    void set_laben_multiplicity         (int lg, const std::vector<std::string>& multiplicity, const bx_named* obj)  { if (check_acl ("set_laben_multiplicity", obj)) laben_multiplicity_v[lg] = multiplicity;}  
    void set_muon_multiplicity          (int lg, const std::vector<std::string>& multiplicity, const bx_named* obj)  { if (check_acl ("set_muon_multiplicity", obj)) muon_multiplicity_v[lg-3001] = multiplicity;} // offset to be removed if reverted back to map container
    void set_muon_alignment             (int n, time_t t, int un, time_t ut, int an, time_t at, const bx_named* obj)  { if (check_acl ("set_muon_alignment", obj)) i_muon_nevents = n, t_muon_time = t, i_muon_up_nevents = un, t_muon_up_time = ut, i_muon_aligned_nevents = an, t_muon_aligned_time = at;}

    //InnerDetectorDarkRateParametrisation
    void set_mean_dark_rate_per_used_pmt        (float mean_dark_rate_per_used_pmt,       const bx_named* obj) { if (check_acl ("set_mean_dark_rate_per_used_pmt"       , obj)) f_mean_dark_rate_per_used_pmt       = mean_dark_rate_per_used_pmt      ;}
    void set_error_mean_dark_rate_per_used_pmt  (float error_mean_dark_rate_per_used_pmt, const bx_named* obj) { if (check_acl ("set_error_mean_dark_rate_per_used_pmt" , obj)) f_error_mean_dark_rate_per_used_pmt = error_mean_dark_rate_per_used_pmt;}
    void set_mu_win1                            (float mu_win1,                           const bx_named* obj) { if (check_acl ("set_mu_win1"                           , obj)) f_mu_win1                           = mu_win1                          ;}
    void set_mu_win2                            (float mu_win2,                           const bx_named* obj) { if (check_acl ("set_mu_win2"                           , obj)) f_mu_win2                           = mu_win2                          ;}
    void set_pois_const_win1                    (float pois_const_win1,                   const bx_named* obj) { if (check_acl ("set_pois_const_win1"                   , obj)) f_pois_const_win1                   = pois_const_win1                  ;}            
    void set_error_pois_const_win1              (float error_pois_const_win1,             const bx_named* obj) { if (check_acl ("set_error_pois_const_win1"             , obj)) f_error_pois_const_win1             = error_pois_const_win1            ;}            
    void set_mu_fit_win1                        (float mu_fit_win1,                       const bx_named* obj) { if (check_acl ("set_mu_fit_win1"                       , obj)) f_mu_fit_win1                       = mu_fit_win1                      ;}                
    void set_error_mu_fit_win1                  (float error_mu_fit_win1,                 const bx_named* obj) { if (check_acl ("set_error_mu_fit_win1"                 , obj)) f_error_mu_fit_win1                 = error_mu_fit_win1                ;}                
    void set_exp_const_win1                     (float exp_const_win1,                    const bx_named* obj) { if (check_acl ("set_exp_const_win1"                    , obj)) f_exp_const_win1                    = exp_const_win1                   ;}
    void set_error_exp_const_win1               (float error_exp_const_win1,              const bx_named* obj) { if (check_acl ("set_error_exp_const_win1"              , obj)) f_error_exp_const_win1              = error_exp_const_win1             ;}
    void set_tau_win1                           (float tau_win1,                          const bx_named* obj) { if (check_acl ("set_tau_win1"                          , obj)) f_tau_win1                          = tau_win1                         ;}
    void set_error_tau_win1                     (float error_tau_win1,                    const bx_named* obj) { if (check_acl ("set_error_tau_win1"                    , obj)) f_error_tau_win1                    = error_tau_win1                   ;}
    void set_pois_const_win2                    (float pois_const_win2,                   const bx_named* obj) { if (check_acl ("set_pois_const_win2"                   , obj)) f_pois_const_win2                   = pois_const_win2                  ;}            
    void set_error_pois_const_win2              (float error_pois_const_win2,             const bx_named* obj) { if (check_acl ("set_error_pois_const_win2"             , obj)) f_error_pois_const_win2             = error_pois_const_win2            ;}            
    void set_mu_fit_win2                        (float mu_fit_win2,                       const bx_named* obj) { if (check_acl ("set_mu_fit_win2"                       , obj)) f_mu_fit_win2                       = mu_fit_win2                      ;}                
    void set_error_mu_fit_win2                  (float error_mu_fit_win2,                 const bx_named* obj) { if (check_acl ("set_error_mu_fit_win2"                 , obj)) f_error_mu_fit_win2                 = error_mu_fit_win2                ;}                
    void set_exp_const_win2                     (float exp_const_win2,                    const bx_named* obj) { if (check_acl ("set_exp_const_win2"                    , obj)) f_exp_const_win2                    = exp_const_win2                   ;}
    void set_error_exp_const_win2               (float error_exp_const_win2,              const bx_named* obj) { if (check_acl ("set_error_exp_const_win2"              , obj)) f_error_exp_const_win2              = error_exp_const_win2             ;}
    void set_tau_win2                           (float tau_win2,                          const bx_named* obj) { if (check_acl ("set_tau_win2"                          , obj)) f_tau_win2                          = tau_win2                         ;}
    void set_error_tau_win2                     (float error_tau_win2,                    const bx_named* obj) { if (check_acl ("set_error_tau_win2"                    , obj)) f_error_tau_win2                    = error_tau_win2                   ;}
   
    void write_laben_laser_time		        (bool y, const bx_named* obj) { if (check_acl ("write_laben_laser_time"        , obj)) b_write_laben_laser_time         = y; }
    void write_laben_laser_charge	        (bool y, const bx_named* obj) { if (check_acl ("write_laben_laser_charge"      , obj)) b_write_laben_laser_charge       = y; }
    void write_laben_tt1_charge	                (bool y, const bx_named* obj) { if (check_acl ("write_laben_tt1_charge"        , obj)) b_write_laben_tt1_charge         = y; }
    void write_muon_laser_time		        (bool y, const bx_named* obj) { if (check_acl ("write_muon_laser_time"         , obj)) b_write_muon_laser_time          = y; }
    void write_muon_laser_charge	        (bool y, const bx_named* obj) { if (check_acl ("write_muon_laser_charge"       , obj)) b_write_muon_laser_charge        = y; }
    void write_laben_dark_rates 	        (bool y, const bx_named* obj) { if (check_acl ("write_laben_dark_rates"        , obj)) b_write_laben_dark_rates         = y; }
    void write_muon_dark_rates   	        (bool y, const bx_named* obj) { if (check_acl ("write_muon_dark_rates"         , obj)) b_write_muon_dark_rates          = y; }
    void write_trigger_parameters 	        (bool y, const bx_named* obj) { if (check_acl ("write_trigger_parameters"      , obj)) b_write_trigger_parameters       = y; }
    void write_laben_electronic_channel	        (bool y, const bx_named* obj) { if (check_acl ("write_laben_electronic_channel", obj)) b_write_laben_electronic_channel = y; }
    void write_muon_electronic_channel	        (bool y, const bx_named* obj) { if (check_acl ("write_muon_electronic_channel" , obj)) b_write_muon_electronic_channel  = y; }
    void write_muon_alignment	                (bool y, const bx_named* obj) { if (check_acl ("write_muon_alignment"          , obj)) b_write_muon_alignment           = y; }
    void write_laben_dark_rates_parametrisation	(bool y, const bx_named* obj) { if (check_acl ("write_laben_dark_rates_parametrisation", obj)) b_write_laben_dark_rates_parametrisation = y; }
                                                                                                                                  
#endif                                                                                                                            

  private:
      // Daq run data
    int      i_run_number;
    int      i_profile_id;
    int      i_calib_profile;
    run_type i_run_type;
    int      i_max_time;
    int      i_max_events;
    int      i_stop_reason;
    int      i_number_of_files;
    time_t   u4_start_time;
    time_t   u4_stop_time;
    int      i_duration;
    int      i_events;

      // Precalib data
    bool                  b_precalib_present, b_precalib_quality_present;
    bool                  b_write_laben_precalib, b_write_laben_precalib_quality;
    bool                  b_write_muon_precalib;
    std::map<int, float>  laben_precalib_delta80_v;
    std::map<int, int>    laben_precalib_low_bin_v, laben_precalib_high_bin_v, laben_precalib_gray_shift_v;
    std::map<int, bool>   laben_precalib_rising_on_even_v;
    uint32_t	          u4_muon_precalib_pulse_time;
    std::map<int, double> muon_precalib_pedestal_v, muon_precalib_pedsigma_v;
    float                 f_laben_precalib_mean_time, f_laben_precalib_sigma_time;
    std::vector<int>      laben_precalib_bad_channels_v, laben_precalib_off_channels_v;
    
      // Calib data (laser) 
    bool b_laben_laser_present;
    bool b_write_laben_laser_time;
    bool b_write_laben_laser_charge;
    bool b_laben_calib_tt1_present;
    bool b_write_laben_tt1_charge;
    bool b_muon_laser_present;
    bool b_write_muon_laser_time;
    bool b_write_muon_laser_charge;
    std::map<int, float> laben_time_offset_v, laben_time_sigma_v, laben_charge_peak_v, laben_charge_sigma_v;
    std::map<int, float> laben_charge_tt1_peak_v, laben_charge_tt1_sigma_v, laben_charge_tt1_mean_v, laben_charge_tt1_rms_v, laben_charge_tt1_p0_v;
    std::map<int, float> muon_time_offset_v, muon_time_sigma_v, muon_charge_peak_v, muon_charge_sigma_v;

      // Calib data (gate)
    float f_laben_gate_width, f_laben_gate_start;
    float f_laben_laser_offset, f_laben_pulser_offset, f_laben_cluster_offset;

      // Calib data (dark noise)
    float f_laben_mean_dark_noise, f_laben_mean_dark_sigma;
    int   i_laben_dead_cone, i_laben_dead_no_cone, i_laben_hot_cone, i_laben_hot_no_cone;
    std::map<int, float> laben_dark_noise_v, laben_dark_sigma_v;
    std::map<int, std::string> laben_pmt_status_v;
    float f_muon_mean_dark_noise, f_muon_mean_dark_sigma;
    int   i_muon_dead, i_muon_hot;
    std::map<int, float> muon_dark_noise_v, muon_dark_sigma_v;
    std::map<int, std::string> muon_pmt_status_v;

      // Calib data (electronics calibration)
   /*std::map<int, std::vector<std::string> > laben_charge_base_status_v;
    std::map<int, std::vector<std::string> > laben_charge_peak_status_v;
    std::map<int, std::vector<std::string> > laben_charge_status_v;
    std::map<int, std::vector<std::string> > laben_timing_status_v;
    std::map<int, std::vector<std::string> > laben_multiplicity_v;*/
    std::vector<std::string> laben_charge_base_status_v[2241];
    std::vector<std::string> laben_charge_peak_status_v[2241];
    std::vector<std::string> laben_charge_status_v[2241];
    std::vector<std::string> laben_timing_status_v[2241];
    std::vector<std::string> laben_multiplicity_v[2241];
    std::vector<std::string> muon_multiplicity_v[256];
    std::map<int, bool> disconnected_pmts_v;

    std::map<int, std::map<int, db_run::disabled_type> > disabled_channels_v;

      //effectiveQE
    std::map<int, float> laben_qe_v;//effective QE from DB normalized by mean
    std::map<int, float> laben_qe_nocorr_v;//effective QE from DB normalized by mean without concentrator corr.
    std::map<int, float> laben_qe_dark_noise_v;//mean dark rate in eff. QE computation period in Hz

      // Muon alignement
    int i_muon_nevents, i_muon_up_nevents, i_muon_aligned_nevents;
    time_t t_muon_time, t_muon_up_time, t_muon_aligned_time;

      // Inner Detector Dark Rate Parametrisation 
    float f_mean_dark_rate_per_used_pmt, f_error_mean_dark_rate_per_used_pmt, f_mu_win1, f_mu_win2, f_pois_const_win1, f_error_pois_const_win1;
    float f_mu_fit_win1, f_error_mu_fit_win1, f_exp_const_win1, f_error_exp_const_win1, f_tau_win1, f_error_tau_win1, f_pois_const_win2;
    float f_error_pois_const_win2, f_mu_fit_win2, f_error_mu_fit_win2, f_exp_const_win2, f_error_exp_const_win2, f_tau_win2, f_error_tau_win2;

      // Status flags
    bool b_trigger_parameters_present;
    bool b_write_trigger_parameters;
    bool b_laben_dark_rates_present;
    bool b_write_laben_dark_rates;
    bool b_muon_dark_rates_present;
    bool b_write_muon_dark_rates;
    bool b_laben_electronic_channel_present;
    bool b_write_laben_electronic_channel;
    bool b_muon_electronic_channel_present;
    bool b_write_muon_electronic_channel;
    bool b_muon_alignment_present;
    bool b_write_muon_alignment;
    bool b_disabled_channels_present;
    bool b_write_disabled_channels;
    bool b_laben_dark_rates_parametrisation_present;
    bool b_write_laben_dark_rates_parametrisation;

      // Internal map to decode run types
    std::map<std::string, db_run::run_type> run_type_map;

      // Internal handler
    void m_read_runinfo                  (int run, bx_dbi *dbi);
    void m_read_laben_precalib           (int run, bx_dbi *dbi);
    void m_read_laben_precalib_quality   (int run, bx_dbi *dbi);
    void m_read_muon_precalib            (int run, bx_dbi *dbi);
    void m_read_laben_laser_calib        (int run, bx_dbi *dbi);
    void m_read_laben_tt1_calib          (int run, bx_dbi *dbi);
    void m_read_muon_laser_calib         (int run, bx_dbi *dbi);
    void m_read_laben_dark_rates         (int run, bx_dbi *dbi);
    void m_read_muon_dark_rates          (int run, bx_dbi *dbi);
    void m_read_trigger_parameters       (int run, bx_dbi *dbi);
    void m_read_laben_electronic_channel (int run, bx_dbi *dbi);
    void m_read_muon_electronic_channel  (int run, bx_dbi *dbi);
    void m_read_muon_alignment           (int run, bx_dbi *dbi);
    void m_read_disconnected_pmts	 (int run, bx_dbi *dbi);
    void m_read_disabled_channels	 (int run, bx_dbi *dbi);
    void m_read_laben_qe                 (int run, bx_dbi *dbi);
    void m_read_laben_dark_rates_parametrisation(int run, bx_dbi *dbi);
    void m_write_laben_precalib                         ();
    void m_write_laben_precalib_quality                 ();
    void m_write_muon_precalib                          ();
    void m_write_laben_laser_calib                      ();
    void m_write_laben_tt1_calib                        ();
    void m_write_muon_laser_calib                       ();
    void m_write_trigger_parameters                     ();
    void m_write_laben_dark_rates                       ();
    void m_write_muon_dark_rates                        ();
    void m_write_laben_electronic_channel               ();
    void m_write_muon_electronic_channel                ();
    void m_write_muon_alignment                         ();
    void m_write_disabled_channels	                ();
    void m_write_laben_dark_rates_parametrisation       ();
   

      // Private ctors and dctor only available to bx_dbi
    db_run (int run_number);
    void flush (); // flush data to bx_dbi


      // Only bx_dbi can istantiate and destroy db_run objects
    friend class bx_dbi;
    ClassDef(db_run,CYCLE_NUMBER)
};
#endif
/*  
 *  $Log: db_run.hh,v $
 *  Revision 1.58  2015/07/24 12:37:24  ilia.drachnev
 *  added effective QE vectors
 *
 *  Revision 1.57  2015/01/09 15:03:08  misiaszek
 *  cycle_18 new unstable
 *
 *  Revision 1.56  2013/06/18 18:56:37  razeto
 *  cycle_17 new unstable
 *
 *  Revision 1.55  2013-02-02 09:01:50  razeto
 *  Incremented to cycle_16 (cycle 15 was lost)
 *
 *  Revision 1.54  2013-01-16 22:47:06  misiaszek
 *  Visitors & getters for InnerDetectorDarkRateParametrisation added
 *
 *  Revision 1.53  2011-04-19 05:54:58  razeto
 *  Moved to cycle 15 unstable
 *
 *  Revision 1.52  2010-08-06 17:20:16  razeto
 *  Moving to cycle 14
 *
 *  Revision 1.51  2009-12-03 16:13:53  misiaszek
 *  Changes for new streamer and rootcint in disabled_channels_v and run_type_map
 *
 *  Revision 1.50  2009-11-26 13:42:52  razeto
 *  Moved to cycle_13_unstable
 *
 *  Revision 1.49  2008-12-15 17:13:55  razeto
 *  New cycle (12)
 *
 *  Revision 1.48  2008-11-26 13:53:45  ludhova
 *  mean, rms and P0 columns in NeutrinoPmtCalibration table
 *
 *  Revision 1.47  2008-11-17 14:46:31  ludhova
 *  visitors for NeutrinoPmtCalibration DB table
 *
 *  Revision 1.46  2008-10-17 13:41:12  razeto
 *  new development cycle (11)
 *
 *  Revision 1.45  2008-09-25 11:45:37  ludhova
 *  visitors for DisabledChannels
 *
 *  Revision 1.44  2008-09-24 17:13:26  razeto
 *  Added disabled_channels runtime map
 *
 *  Revision 1.43  2008-09-04 15:32:26  ddangelo
 *  added a getter
 *  fixed a typo
 *
 *  Revision 1.42  2008-08-26 18:05:07  ddangelo
 *  muon alignement: times implemented as timestamps
 *
 *  Revision 1.41  2008-08-25 16:37:58  ddangelo
 *  muon alignment flagging implemented. varaibles, sgetters, read/write methods, flags.
 *  to be tested
 *
 *  Revision 1.40  2008-08-19 15:02:40  ddangelo
 *  debug
 *
 *  Revision 1.39  2008-08-11 12:41:47  ddangelo
 *  added visitors for muon electronics calibration (multiplicity only for the moment)
 *  data, s/getters, read/write methods, initializiation, control flags and flushing
 *  b_electronic_channel_present modfied to b_laben_electronic_channel_present to avoid conflicts
 *
 *  Revision 1.38  2008-02-27 20:46:13  razeto
 *  new development cycle (10)
 *
 *  Revision 1.37  2008-02-27 20:26:30  razeto
 *  New clasdef(9) version and new cycle version
 *
 *  Revision 1.36  2007-11-09 18:38:58  razeto
 *  Added disabled pmts (thanks to Marcin for SQL support)
 *
 *  Revision 1.35  2007-10-11 10:49:54  razeto
 *  Cycle 8 deployed
 *
 *  Revision 1.34  2007-06-22 15:15:27  razeto
 *  Moved to cycle 7
 *
 *  Revision 1.33  2007-05-07 15:47:28  razeto
 *  Cycle number in root classdef
 *
 *  Revision 1.32  2007-02-23 19:21:33  ddangelo
 *  visitor upgraded to includ muon dark rates data. both r/w implemented.
 *  some other variables renamed to avoid name clashes: interface unbroken.
 *
 *  Revision 1.31  2007/02/19 15:11:23  razeto
 *  Added new run types
 *
 *  Revision 1.30  2007/01/30 11:03:07  razeto
 *  When doing precalib alway start from scratch
 *
 *  Revision 1.29  2006/11/17 11:42:36  razeto
 *  Added water run type
 *
 *  Revision 1.28  2006/11/05 10:27:30  razeto
 *  Remove some useless include
 *
 *  Revision 1.27  2006/09/11 15:07:08  razeto
 *  Added map name to map_get
 *
 *  Revision 1.26  2006/09/07 14:11:41  ddangelo
 *  added muon laser calibration parameters, sgetters and acl.
 *
 *  Revision 1.25  2006/08/30 09:43:09  ludhova
 *  new visitors
 *
 *  Revision 1.24  2006/07/19 10:29:23  dfranco
 *  Added visitors for dark rates and detector efficiency
 *
 *  Revision 1.23  2006/07/18 14:48:20  razeto
 *  Added support for PrecalibDecodingQuality writing
 *
 *  Revision 1.22  2006/07/18 10:42:26  razeto
 *  Merged all reading code from Ludhova
 *
 *  Revision 1.21  2006/07/18 09:56:42  dfranco
 *  Reorganization of the methods
 *
 *  Revision 1.20  2006/06/21 12:36:43  ludhova
 *  added precalibration_quality (missing write)
 *
 *  Revision 1.19  2006/01/25 13:24:01  misiaszek
 *  Moved from cmap to simple map (to work with root)
 *
 *  Revision 1.18  2005/05/05 17:12:34  monzani
 *  Changed the names of the visitors for laser calibrations (Laben specified).
 *
 *  Revision 1.17  2005/05/05 10:06:01  monzani
 *  DB update from laser calibration modules added.
 *
 *  Revision 1.16  2004/11/26 17:42:56  razeto
 *  Added muon precalib handling
 *
 *  Revision 1.15  2004/11/26 15:25:10  razeto
 *  Thanks to Davide who so kindly underlined a syntactical error, I fixed the "Maintainer" word
 *
 *  Revision 1.14  2004/11/26 14:01:20  razeto
 *  Added Mantainer field
 *
 *  Revision 1.13  2004/11/24 13:05:37  razeto
 *  Moved some variable position to simplify the ctor
 *
 *  Revision 1.12  2004/11/24 09:46:41  razeto
 *  Moved some prototyping in db_acl from db_*
 *
 *  Revision 1.11  2004/10/19 18:41:58  razeto
 *  Integrated visitors writing in the framework
 *
 *  Revision 1.10  2004/10/19 16:30:41  razeto
 *  Added laben predecoding reading.
 *  Added predecoding writing to the database (still experimental).
 *
 *  Revision 1.9  2004/08/14 21:21:40  razeto
 *  Added some fields used from calibrations (and some improvements)
 *
 *  Revision 1.8  2004/06/10 10:06:48  razeto
 *  Added some run type, and changed an error condition
 *
 *  Revision 1.7  2004/05/26 09:23:04  razeto
 *  Upgraded know database predecoding status
 *
 *  Revision 1.6  2004/05/18 14:26:59  razeto
 *  Updated
 *
 *  Revision 1.5  2004/04/26 13:48:29  razeto
 *  Added db_acl to check set calls for calling module to have the right privileges.
 *  Fixed the names of set/get methods.
 *  Modifications decided at the software Paris meeting.
 *
 *  Revision 1.4  2004/04/24 17:41:11  razeto
 *  Added a void method. To be developed
 *
 *  Revision 1.3  2004/04/20 11:37:31  ddangelo
 *  added a few muon precalibration varaibles and relative sgetters.
 *
 *  Revision 1.2  2004/04/12 16:10:58  razeto
 *  Updated
 *
 *  Revision 1.1  2004/04/09 08:05:00  razeto
 *  Added
 *
 */
