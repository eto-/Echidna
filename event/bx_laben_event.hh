/* BOREXINO Reconstruction program
 *
 * Author: Davide D'Angelo <davide.dangelo@lngs.infn.it>
 * based on work by Razeto&Pallas 
 * Maintainer: Davide D'Angelo <davide.dangelo@lngs.infn.it>
 *
 * $Id: bx_laben_event.hh,v 1.135 2015/08/26 11:21:44 misiaszek Exp $
 *
 * The laben event object
 * 
 */
#ifndef _BX_LABEN_EVENT_HH
#define _BX_LABEN_EVENT_HH

#include <vector>

#include "bx_rec_general.hh"
#include "bx_base_event.hh"
#include "constants.hh"
#include "bx_track.hh"

// ------- RAW --------- //
class bx_laben_raw_hit {
  public:
    bx_laben_raw_hit (const char *disk_hit, uint8_t order_in_channel);
    uint16_t get_logical_channel  ()  const { return u2_channel; }
    uint8_t  get_time_1           ()  const { return u1_time_1; }
    uint8_t  get_time_2           ()  const { return u1_time_2; }
    uint16_t get_gray_counter     ()  const { return u2_gray_counter; }
    uint8_t  get_base             ()  const { return u1_base; }
    uint8_t  get_peak             ()  const { return u1_peak; }
    uint16_t get_flags_board      ()  const { return (u2_flags >> 8) & 0xff; }
    uint16_t get_flags_ch         ()  const { return u2_flags & 0xff; }
    enum flags {
      good = 0,
      fifo_full,
      fifo_empty,
      counter,
      trg_jump,
      trg_jump_large,
      trg_in_busy,
      invalid,
      __max__,
    };
    bool	   check_flag		(flags flag) const { bool f = u2_flags & flags_bits[flag]; return (flag != good) ? f : !f; } // there is no flag for good
    uint8_t  get_order_in_channel ()  const { return u1_order_in_channel; }
  private:
    uint16_t u2_channel;
    uint8_t  u1_time_1, u1_time_2;
    uint16_t u2_gray_counter;
    uint8_t  u1_base, u1_peak;
    uint16_t u2_flags;
    uint8_t  u1_order_in_channel; // as in daq fifo: all hits are counted
    static uint16_t flags_bits[__max__]; // KEEP alligned to flags enum !!!
  public:
    typedef std::vector<bx_laben_raw_hit> bx_laben_raw_hit_vector;
};

class bx_laben_raw_event {
  public:
    bx_laben_raw_event (const char *disk_event);
    virtual ~bx_laben_raw_event () {}
    const bx_laben_raw_hit& get_raw_hit (int i) const { return raw_hits[i]; }
    int get_raw_nhits () const { return raw_hits.size (); }
    int get_raw_nhits_fw () const { return i4_nhits_fw; }
    int get_raw_nhits_flag (bx_laben_raw_hit::flags flag_id) const { return nhits_flag[flag_id]; }
    int get_empty_boards () const { return i4_empty_boards; }
  private:
    uint32_t u4_errors;
    int nhits_flag[bx_laben_raw_hit::__max__]; // keep alligned to size of bx_laben_raw_hit::flags enum
    int i4_empty_boards;
    int i4_nhits_fw;
    bx_laben_raw_hit::bx_laben_raw_hit_vector raw_hits;
  friend class bx_laben_raw_validator;
};


// ------- DECODED --------- //
class db_channel_laben;
class laben_time_hit;
class laben_charge_hit;
class bx_laben_decoded_hit {
  public:
    bx_laben_decoded_hit (const bx_laben_raw_hit &hit, uint16_t index): u1_flag(0), raw_hit(&hit) {}
    const bx_laben_raw_hit& get_raw_hit () const { return *raw_hit; }

      // TDC data from the laben_hit class
    double get_raw_time   ()                  const { return f8_raw_time;   } // Relative to a gray counter crossing window
    float  get_d80        ()                  const { return f4_d80;        }
    float  get_time_error ()                  const { return f4_time_error; }
    bool   is_timing_good (float error = .45) const { return bool(f4_time_error <= error); }
    enum flag_type {
      out_of_gate = 1,
      reflection  = 2,
      reference   = 4,
      retrigger   = 8,
      disabled    = 16,
    };
    uint8_t get_flag ()		  const { return u1_flag;               }
    bool           is_good ()		  const { return !u1_flag;              } 
    bool    is_out_of_gate ()		  const { return u1_flag & out_of_gate; } 
    bool     is_reflection ()		  const { return u1_flag & reflection;  } 
    bool      is_reference ()		  const { return u1_flag & reference;   } 
    bool      is_retrigger ()		  const { return u1_flag & retrigger;   } 
    uint8_t get_order_in_channel () const { return u1_order_in_channel;   }
    bool operator< (const bx_laben_decoded_hit& hit) const { return get_raw_time () < hit.get_raw_time (); }
	
      // ADC data
    float get_charge_bin ()             const { return f4_charge_bin; }
    float get_uncorrected_charge_bin () const { return f4_uncorrected_charge_bin; }
    int   get_charge_npe ()		const { return i4_charge_npe; }
    float get_charge_pe ()              const { return f4_charge_pe; }
    float get_uncorrected_charge_pe ()  const { return f4_uncorrected_charge_pe; }
    float get_charge_mean_pe ()              const { return f4_charge_mean_pe; }

    const db_channel_laben* get_db_channel () const { return p_db_channel; }

  private:
      // From laben_time_hit
    double f8_raw_time;
    float f4_d80, f4_time_error;
      // From laben_charge_hit
    int i4_charge_npe;
    float f4_charge_bin, f4_uncorrected_charge_bin;
    float f4_charge_pe, f4_uncorrected_charge_pe, f4_charge_mean_pe;
      // From laben_decoder
    uint8_t  u1_flag;
    uint8_t  u1_order_in_channel; // only decoded hits are counted
    const bx_laben_raw_hit* raw_hit;
    const db_channel_laben *p_db_channel;

    const bx_laben_decoded_hit& operator= (laben_time_hit& t_hit);
    const bx_laben_decoded_hit& operator= (laben_charge_hit& c_hit);

  public:
    typedef std::vector<bx_laben_decoded_hit> bx_laben_decoded_hit_vector;
  friend class bx_laben_decoder;
};

class bx_laben_decoded_event {
  public:
    bx_laben_decoded_event ();
    virtual ~bx_laben_decoded_event () {}

    const bx_laben_decoded_hit& get_decoded_hit (int i) const { return decoded_hits[i]; }
    int get_decoded_nhits () const { return decoded_hits.size (); }
    int get_npmts   () const { return i4_npmts; }
    float get_npe() const { return f4_npe; }
    float get_charge() const { return f4_charge; }	

    double get_trigger_rawt () const { return f8_trigger_rawt; }
    int get_trigger_reference_decoded_nhits () const { return trigger_reference_decoded_hits.size (); }
    const bx_laben_decoded_hit& get_trigger_reference_decoded_hit (int i) const { return trigger_reference_decoded_hits[i]; }

    double get_laser_rawt () const { return f8_laser_rawt; }
    int get_laser_reference_decoded_nhits () const { return laser_reference_decoded_hits.size (); }
    const bx_laben_decoded_hit& get_laser_reference_decoded_hit (int i) const { return laser_reference_decoded_hits[i]; }

    int get_n_live_pmts   () const { return i4_n_live_pmts; }
    int get_n_live_charge () const { return i4_n_live_charge; }
    float normalize_charge (float arg) const { return i4_n_live_charge ? (arg / i4_n_live_charge * constants::laben::pmts) : 0. ;}
    int normalize_pmts (int arg) const { return i4_n_live_pmts ? int (::roundf(arg / i4_n_live_pmts * constants::laben::pmts)) : 0 ;}

    int get_invalid_pmts   () const { return i4_invalid_pmts; }
    int get_invalid_charge () const { return i4_invalid_charge; }
    int get_nhits_on_empty () const { return i4_nhits_on_empty; }

private:
    bx_laben_decoded_hit::bx_laben_decoded_hit_vector decoded_hits, trigger_reference_decoded_hits, laser_reference_decoded_hits;
    double f8_trigger_rawt, f8_laser_rawt;
    int i4_npmts;
    float f4_npe, f4_charge;
    int i4_n_live_pmts;
    int i4_n_live_charge;
    int i4_invalid_pmts;
    int i4_invalid_charge;
    int i4_nhits_on_empty;
  friend class bx_laben_decoder;
};


// ------- CLUSTERED --------- //
class bx_laben_clustered_hit {
  public:
    bx_laben_clustered_hit (const bx_laben_decoded_hit &hit, uint16_t index) : f8_time(0.), decoded_hit(&hit) {}
    const bx_laben_decoded_hit& get_decoded_hit () const { return *decoded_hit; }

    double        get_time             () const { return f8_time;             }
    uint8_t get_order_in_channel () const { return u1_order_in_channel; }
    bool          is_short_cluster     () const { return b_short_cluster     ; }
    bool operator< (const bx_laben_clustered_hit& hit) const { return get_time () < hit.get_time (); }

  private:
    double f8_time;
    uint8_t  u1_order_in_channel; // only hits in cluster are counted
    bool b_short_cluster;
    const bx_laben_decoded_hit* decoded_hit;
    
  public:
    typedef std::vector<bx_laben_clustered_hit> bx_laben_clustered_hit_vector;
  friend class bx_laben_findcluster;
  friend class bx_laben_findcluster_time;
};

class bx_base_position {
  public:
    bx_base_position () : f4_t(0.), f4_x(10.), f4_y(10.), f4_z(10.), f4_dt(0.), f4_dx(0.), f4_dy(0.), f4_dz(0.), f4_user(0.), b_converged(false), i4_matrix(0) {}
    bx_base_position (float x, float y, float z) : f4_t(0.), f4_x(x), f4_y(y), f4_z(z), f4_dt(0.), f4_dx(0.), f4_dy(0.), f4_dz(0.), b_converged(false), i4_matrix(0) {}
    virtual ~bx_base_position () {}
    
    float get_t  () const { return f4_t; }
    float get_x  () const { return f4_x; }
    float get_y  () const { return f4_y; }
    float get_z  () const { return f4_z; }
    float get_dt () const { return f4_dt; }
    float get_dx () const { return f4_dx; }
    float get_dy () const { return f4_dy; }
    float get_dz () const { return f4_dz; }
    float get_r () const { return ::sqrtf (f4_x * f4_x + f4_y * f4_y + f4_z * f4_z); }
    float get_likelihood () const { return f4_user; }
    float get_ref_index  () const { return f4_user; }
    float get_user       () const { return f4_user; }
    bool  is_converged   () const { return b_converged; }
    bool  is_pos_def     () const { return i4_matrix== 1; }
    bool  is_not_pos_def () const { return i4_matrix==-1; }
    bool  is_approximate () const { return i4_matrix== 0; }
    int   get_matrix     () const { return i4_matrix; }

  protected:
    float f4_t, f4_x, f4_y, f4_z;
    float f4_dt, f4_dx, f4_dy, f4_dz; 
    float f4_user;
    bool  b_converged;
    int   i4_matrix;
    
  friend class bx_baricentrator;
  friend class bx_position_reco_mi;
  friend class bx_position_reco_lngs;
  friend class bx_position_reco_noavg;
  friend class bx_position_reco_msk;
  friend class bx_position_reco_dbn;
  friend class bx_position_reco_mach4;
};

class bx_base_energy {
  public:
    bx_base_energy () : i4_nhits(0), i4_npe(0), f4_charge(0.) {}
    virtual ~bx_base_energy () {}

    int   get_nhits  () const { return i4_nhits; }
    int   get_npe    () const { return i4_npe; }
    float get_charge () const { return f4_charge; }

  protected:
    int   i4_nhits;
    int   i4_npe;
    float f4_charge;

  friend class bx_energy_reco_mc;
  friend class bx_energy_reco_lik;
  friend class bx_energy_reco_msk;
  friend class bx_energy_reco_dbn;
};

class bx_baricenter: public bx_base_position {
  public:
    bx_baricenter () : f4_best_error_radius(0.) {}
    virtual ~bx_baricenter () {}
    
    float get_best_error_radius () const { return f4_best_error_radius; }
 
  private:
    float f4_best_error_radius;

  friend class bx_baricentrator;
};

class bx_position_mi: public bx_base_position {
  public:
    bx_position_mi () {}
    virtual ~bx_position_mi () {}
    
  friend class bx_position_reco_mi;
};

class bx_position_lngs: public bx_base_position {
  public:
    bx_position_lngs () {}
    virtual ~bx_position_lngs () {}
    
  friend class bx_position_reco_lngs;
};

class bx_position_noavg: public bx_base_position {
  public:
    bx_position_noavg () {}
    virtual ~bx_position_noavg () {}

  friend class bx_position_reco_noavg;
};


class bx_position_dbn: public bx_base_position {
  public:
    bx_position_dbn () {}
    virtual ~bx_position_dbn () {}
  
  friend class bx_position_reco_dbn;
};

class bx_position_msk: public bx_base_position {
  public:
    bx_position_msk () {}
    virtual ~bx_position_msk () {}

  friend class bx_position_reco_msk;
};

class bx_position_mach4: public bx_base_position {
  public:
    bx_position_mach4 () {}
    virtual ~bx_position_mach4 () {}

  friend class bx_position_reco_mach4;
};

class bx_energy_mc: public bx_base_energy {
  public:
    bx_energy_mc () {}
    virtual ~bx_energy_mc () {}

  friend class bx_energy_reco_mc;
};

class bx_energy_lik: public bx_base_energy {
  public:
    bx_energy_lik () {}
    virtual ~bx_energy_lik () {}

  friend class bx_energy_reco_lik;
};

class bx_energy_dbn: public bx_base_energy {
  public:
    bx_energy_dbn () {}
    virtual ~bx_energy_dbn () {}

  friend class bx_energy_reco_dbn;
};

class bx_energy_msk: public bx_base_energy {
  public:
    bx_energy_msk () {}
    virtual ~bx_energy_msk () {}

    float get_nphotons () const { return f4_nphotons; }

  private:
    float f4_nphotons;

  friend class bx_energy_reco_msk;
};

class bx_split_peak {
  public:
    float get_start_time() const { return f4_start_time; } 
    float get_nhits     () const { return f4_nhits;      }

  private:
    bx_split_peak ()               : f4_start_time(0.), f4_nhits(0.) {}
    bx_split_peak (float t, float n) : f4_start_time(t) , f4_nhits(n)  {}

    float f4_start_time;
    float f4_nhits;

  friend class bx_splitting_filter;
};

class bx_laben_cluster {
  public:
    bx_laben_cluster ();
    virtual ~bx_laben_cluster () {}

    int    get_npe             () const { return i4_npe;             }
    int    get_npe_conc        () const { return i4_npe_conc;        }
    float  get_charge          () const { return f4_charge;          }
    float  get_charge_conc     () const { return f4_charge_conc;     }
    float  get_charge_thresh   () const { return f4_charge_thresh;   }
    float  get_charge_short    () const { return f4_charge_short;    }
    float  get_charge_400      () const { return f4_charge_400;      }
    float  get_charge_dt1      () const { return f4_charge_dt1;      }
    float  get_charge_dt2      () const { return f4_charge_dt2;      }
    float  get_charge_pos      () const { return f4_charge_pos;      }
    float  get_charge_npmts    () const { return f4_charge_npmts;    }
    float  get_charge_mean     () const { return f4_charge_mean;     }
    float  get_charge_clean    () const { return f4_charge_clean;    }
    float  get_charge_noavg_dt1() const { return f4_charge_noavg_dt1;}
    float  get_charge_noavg_dt2() const { return f4_charge_noavg_dt2;}
    float  get_charge_noavg    () const { return f4_charge_noavg;    }
    float  get_charge_noavg_short() const { return f4_charge_noavg_short;    }
    int    get_npmts           () const { return i4_npmts;           }
    int    get_npmts_conc      () const { return i4_npmts_conc;      }
    int    get_npmts_thresh    () const { return i4_npmts_thresh;    }
    int    get_npmts_short     () const { return i4_npmts_short;     }
    int    get_npmts_400       () const { return i4_npmts_400;       }
    int    get_npmts_dt1       () const { return i4_npmts_dt1;       } 
    int    get_npmts_dt2       () const { return i4_npmts_dt2;       }
    int    get_npmts_pos       () const { return i4_npmts_pos;       }
    double get_start_time      () const { return f8_start_time;      }
//    double get_rough_time    () const { return f8_rough_time;      }
    float  get_mean_time       () const { return f4_mean_time;       }
    float  get_mean_time_short () const { return f4_mean_time_short; }
    float  get_rms_time        () const { return f4_rms_time;        }
    float  get_rms_time_short  () const { return f4_rms_time_short;  }
    float  get_duration        () const { return clustered_hits[clustered_hits.size () - 1].get_time ();    }
    float  get_duration_short  () const { return f4_duration_short;  }
    float  get_npmt_geo_weight      () const { return npmt_geo_weight; }
    float  get_npmt_QE_weight       () const { return npmt_QE_weight; }
    float  get_npmt_geo_QE_weight   () const { return npmt_geo_QE_weight; }
    float  get_charge_geo_weight    () const { return charge_geo_weight; }
    float  get_charge_QE_weight     () const { return charge_QE_weight; }
    float  get_charge_geo_QE_weight () const { return charge_geo_QE_weight; }
    enum flag_type {
      out_of_gate = 1,
      broad	  = 2,
      trigger     = 4,
    };
    uint8_t get_flag ()		  const { return u1_flag;               }
    bool           is_good ()		  const { return !u1_flag;              } 
    bool    is_out_of_gate ()		  const { return u1_flag & out_of_gate; } 
    bool          is_broad ()		  const { return u1_flag & broad;  	} 
    bool        is_trigger ()		  const { return u1_flag & trigger;     } 
    bool        is_neutron ()             const { return b_is_neutron;          }

    const bx_laben_clustered_hit& get_clustered_hit (int i) const { return clustered_hits[i]; }
    int   get_clustered_nhits        ()   const { return clustered_hits.size ()   ; }
    int   get_clustered_nhits_conc   ()   const { return i4_clustered_nhits_conc  ; }
    int   get_clustered_nhits_thresh ()   const { return i4_clustered_nhits_thresh; }
    int   get_clustered_nhits_short  ()   const { return i4_clustered_nhits_short ; }
    int   get_clustered_nhits_400    ()   const { return i4_clustered_nhits_400   ; }
    int   get_clustered_nhits_pos    ()   const { return i4_clustered_nhits_pos   ; }
    float get_clustered_nhits_bkg    ()   const { return f4_clustered_nhits_bkg   ; }

    const bx_baricenter&     get_baricenter     () const { return baricenter    ; }
    const bx_position_mi&    get_position_mi    () const { return position_mi   ; }
    const bx_position_lngs&  get_position_lngs  () const { return position_lngs   ; }
    const bx_position_noavg&  get_position_noavg  () const { return position_noavg   ; }
    const bx_position_dbn&   get_position_dbn   () const { return position_dbn  ; }
    const bx_position_msk&   get_position_msk   () const { return position_msk  ; }
    const bx_position_mach4& get_position_mach4 () const { return position_mach4; }
    const bx_position_mach4& get_position_mach4_fixed () const { return position_mach4_fixed; }
    const bx_energy_mc&      get_energy_mc      () const { return energy_mc     ; }
    const bx_energy_lik&     get_energy_lik     () const { return energy_lik    ; }
    const bx_energy_dbn&     get_energy_dbn     () const { return energy_dbn    ; }
    const bx_energy_msk&     get_energy_msk     () const { return energy_msk    ; }
    bx_baricenter&     get_baricenter       () { return baricenter    ; }
    bx_position_mi&    get_position_mi      () { return position_mi   ; }
    bx_position_lngs&  get_position_lngs    () { return position_lngs ; }
    bx_position_noavg& get_position_noavg   () { return position_noavg; } 
    bx_position_dbn&   get_position_dbn     () { return position_dbn  ; }
    bx_position_msk&   get_position_msk     () { return position_msk  ; }
    bx_position_mach4& get_position_mach4   () { return position_mach4; }
    bx_position_mach4& get_position_mach4_fixed () { return position_mach4_fixed; }
    bx_energy_mc&      get_energy_mc        () { return energy_mc     ; }
    bx_energy_lik&     get_energy_lik       () { return energy_lik    ; }
    bx_energy_dbn&     get_energy_dbn       () { return energy_dbn    ; }
    bx_energy_msk&     get_energy_msk       () { return energy_msk    ; }
        
    const bx_split_peak& get_split_peak (int i) const { return split_peaks[i]; }
    int get_split_npeaks () const { return split_peaks.size (); }
    std::vector<bx_split_peak>& get_split_peaks () { return split_peaks; }

  private:
    int    i4_clustered_nhits_conc, i4_clustered_nhits_thresh, i4_clustered_nhits_short, i4_clustered_nhits_400, i4_clustered_nhits_pos;
    float  f4_clustered_nhits_bkg;
    float  f4_charge, f4_charge_conc, f4_charge_mean, f4_charge_thresh, f4_charge_short, f4_charge_400, f4_charge_dt1,f4_charge_dt2,  f4_charge_pos, f4_charge_npmts, f4_charge_clean;
    float  f4_charge_noavg_dt1,  f4_charge_noavg_dt2,  f4_charge_noavg , f4_charge_noavg_short;
    int    i4_npe,    i4_npe_conc;
    int    i4_npmts,   i4_npmts_conc, i4_npmts_thresh, i4_npmts_short, i4_npmts_400;
    int i4_npmts_dt1, i4_npmts_dt2, i4_npmts_pos;
    double f8_start_time;//, f8_rough_time;
    float  f4_mean_time, f4_mean_time_short, f4_rms_time, f4_rms_time_short;
    float  f4_duration_short;
    uint8_t u1_flag;
    bool b_is_neutron;
    bx_laben_clustered_hit::bx_laben_clustered_hit_vector clustered_hits;
    bx_baricenter   baricenter;
    bx_position_mi  position_mi;
    bx_position_lngs  position_lngs;
    bx_position_noavg position_noavg;
    bx_position_dbn position_dbn;
    bx_position_msk position_msk;
    bx_position_mach4 position_mach4;
    bx_position_mach4 position_mach4_fixed;
    bx_energy_mc    energy_mc;
    bx_energy_lik   energy_lik;
    bx_energy_dbn   energy_dbn;
    bx_energy_msk   energy_msk;
    std::vector<bx_split_peak> split_peaks;

    float npmt_geo_weight;
    float npmt_QE_weight;
    float npmt_geo_QE_weight;
    float charge_geo_weight;
    float charge_QE_weight;
    float charge_geo_QE_weight;
  friend class bx_new_charge_weight;

  public:
    typedef std::vector<bx_laben_cluster> bx_laben_cluster_vector;
  friend class bx_laben_findcluster;
};

class bx_laben_clustered_event {
  public:
    bx_laben_clustered_event ();
    virtual ~bx_laben_clustered_event () {}

    const bx_laben_cluster& get_cluster (int i) const { return clusters[i]; }
    bx_laben_cluster& get_cluster (int i) { return clusters[i]; }
    const bx_laben_cluster& get_cluster_muon (int i) const { return clusters_muons[i]; }
    bx_laben_cluster& get_cluster_muon (int i) { return clusters_muons[i]; }
    int get_nclusters () const { return clusters.size (); }
    int get_nclusters_muons () const { return clusters_muons.size (); }
    int get_nclusters_found () const { return i4_nclusters_found; }
    int get_nclusters_neutron () const { return i4_nclusters_neutron; }
    int get_nclusters_old   () const { return i4_nclusters_old; }
    const std::vector<int>& get_npmts_win1 () const { return v_npmts_win1; }
    const std::vector<int>& get_npmts_win2 () const { return v_npmts_win2; }
    const std::vector<float>& get_charge_win1 () const { return v_charge_win1; }
    const std::vector<float>& get_charge_win2 () const { return v_charge_win2; }
    double get_window_limit () const { return f4_window_limit; }
    

  private:
    bx_laben_cluster::bx_laben_cluster_vector clusters;
    bx_laben_cluster::bx_laben_cluster_vector clusters_muons;
    int i4_nclusters_found, i4_nclusters_old, i4_nclusters_neutron;
    std::vector<int> v_npmts_win1, v_npmts_win2;
    std::vector<float> v_charge_win1, v_charge_win2;
    float f4_window_limit;

  friend class bx_laben_findcluster;
};



// ------- REC --------- //
class bx_laben_rec_hit {
  public:
    bx_laben_rec_hit (const bx_laben_clustered_hit &hit): p_clustered_hit(&hit) {}
    const bx_laben_clustered_hit& get_clustered_hit () const { return *p_clustered_hit; }

    double get_time () const { return f8_time; }
    bool operator< (const bx_laben_rec_hit& hit) const { return get_time () < hit.get_time (); }

  private:
    float f8_time;
    const bx_laben_clustered_hit* p_clustered_hit;

  public:
    typedef std::vector<bx_laben_rec_hit> bx_laben_rec_hit_vector;
  friend class bx_position_reco;
};


class bx_mctruth_frame;
class bx_position: public bx_base_position {
  public:
    bx_position () {}
    virtual ~bx_position () {}

    const bx_base_position& operator= (const bx_base_position& p) { return *(bx_base_position *)this = p; }
    const bx_base_position& operator= (const bx_mctruth_frame&);  // to assign from mctruth

  friend class bx_position_reco;
};


class bx_laben_rec_cluster {
  public:
    bx_laben_rec_cluster (const bx_laben_cluster &cluster): p_cluster(&cluster) { rec_hits.clear(); } // ctor
    virtual ~bx_laben_rec_cluster () {} // dtor

    const bx_laben_cluster& get_cluster () const { return *p_cluster; } // leftward getter
    const bx_laben_rec_hit& get_rec_hit (int i) const { return rec_hits[i]; } // upward getter
    int get_rec_nhits () const { return rec_hits.size (); }

    const bx_position&  get_position  () const { return position;  }
  //    const bx_energy&    get_energy    () const { return energy;    }

  private:
    bx_position position;
  //    bx_energy energy;
    const bx_laben_cluster* p_cluster;
    bx_laben_rec_hit::bx_laben_rec_hit_vector rec_hits;

  friend class bx_position_reco;
};

class bx_laben_shaped_cluster; 
class bx_laben_ab_cluster;

//--- SHAPED ---- //
class bx_laben_shaped_cluster : public bx_laben_rec_cluster {
  public:  
    bx_laben_shaped_cluster(const bx_laben_cluster &cluster);
    virtual ~bx_laben_shaped_cluster() {}

    float get_ns_asymmetry   () const { return f4_ns_asymmetry;   }
    float get_sphere_chi2    () const { return f4_sphere_chi2;    }
    float get_sphere_lkl     () const { return f4_sphere_lkl;     }
    float get_sphere_rel_var () const { return f4_sphere_rel_var; }
    float get_plane_cos      () const { return f4_plane_cos;      }
    float get_plane_chi2     () const { return f4_plane_chi2;     }
    float get_h_plane_chi2   () const { return f4_h_plane_chi2;   }
    const float* get_sh_power() const { return v_sh_power; }
    float get_sh_power (int order) const { return v_sh_power[order];}
    char  get_quality_flags  () const { return i1_quality_flags; }
    bool  is_electronics_hot () const { return i1_quality_flags & 0x1; }
    bool  is_occupancy_low   () const { return i1_quality_flags & 0x2; }

  private:
    float f4_ns_asymmetry;      // north-south charge asymmetry
    float f4_sphere_chi2;    	// chi2 to sphericity
    float f4_sphere_lkl;        // likelihood to sphericity
    float f4_sphere_rel_var;	// relative variance = sigma/mean of the hits in the cos(theta) - phi plane centered in the event position and normalized per hit charge.	
    float f4_plane_cos;	        // cos of the angle made by the plane fitting the theta-phi parameters space.	
    float f4_plane_chi2;        // chi2/NDF of the fit of a plane in a cos(theta)-phi distribution
    float f4_h_plane_chi2;      // chi2/NDF of the fit of a horizontal plane in a cos(theta)-phi distribution
    float v_sh_power[4]; 	// power of spherical armonics		
    char  i1_quality_flags;        // bitfield, meaning to be assigned

  friend class bx_pid_shape;
};

// ---- A/B DISCRIMINATED ----- //
class bx_laben_ab_cluster : public bx_laben_shaped_cluster {
  public:
    bx_laben_ab_cluster (const bx_laben_cluster &cluster) : bx_laben_shaped_cluster(cluster) {}
    virtual ~bx_laben_ab_cluster() {}
		

    float get_tailtot(int tail) const;
    float get_tailtot_ab_mlp(int tail) const;
    float get_tailtot_c11_mva(int tail) const;
    const float* get_tailtot()  const { return v_tailtot; }
    const float* get_tailtot_ab_mlp()  const { return v_tailtot_ab_mlp; }
    const float* get_tailtot_c11_mva()  const { return v_tailtot_c11_mva; }
    float get_gatti      () const { return f4_gatti;      }
    float get_gattic     () const { return f4_gattic;     }
    float get_lkl        () const { return f4_lkl;        }
    float get_lklc       () const { return f4_lklc;       }
    float get_rise_time  () const { return f4_rise_time;  }
    float get_rms	 () const { return rms;           }
    float get_rms_c11	 () const { return rms_c11;       }
    float get_kurtosis   () const { return kurtosis;      }		
    float get_kurtosis_c11() const { return kurtosis_c11; }		
    float get_mlp_ab	 () const { return mlp_ab;        }		

  private:
    float v_tailtot[10];             // tail-to-tot ratio with 10 different tail values 40-130 step 10
    float v_tailtot_ab_mlp[10];      // tail-to-tot ratio with 10 different tail values for a/b MLP algorithm
    float v_tailtot_c11_mva[10];     // tail-to-tot ratio with 10 different tail values for C11/b- MVA algorithms
    float f4_gatti;                  // gatti optimal filter variable
    float f4_gattic;                 // gatti optimal filter variable (cumulative)
    float f4_lkl;                    // likelihood ratio
    float f4_lklc;                   // likelihood ratio (cumulative)
    float f4_rise_time;              // rise_time
    float rms;			     // rms
    float rms_c11;		     // rms for c11 psd
    float kurtosis;		     // kurtosis
    float kurtosis_c11;		     // kurtosis for c11 psd
    float mlp_ab;	             // mlp for ab

  friend class bx_pid_alphabeta;
};

// ---- A/B DISCRIMINATED MACH4 ----- //
/*class bx_laben_ab_mach4_cluster : public bx_laben_ab_cluster {
  public:
    bx_laben_ab_mach4_cluster (const bx_laben_cluster &cluster) : bx_laben_ab_cluster(cluster) {}
    virtual ~bx_laben_ab_mach4_cluster() {}
		
    float get_tailtot_mach4(int tail) const;
    const float* get_tailtot_mach4 () const { return v_tailtot_mach4  ; }
    float get_gatti_mach4     (int i) const { return v_gatti_mach4[i] ; }
    const float* get_gatti_mach4   () const { return v_gatti_mach4    ; }
    float get_peak_mach4           () const { return f4_peak_mach4    ; }
    float get_mean_mach4           () const { return f4_mean_mach4    ; }
    float get_rms_mach4            () const { return f4_rms_mach4     ; }
    float get_skew_mach4           () const { return f4_skew_mach4    ; }
    float get_kurt_mach4           () const { return f4_kurt_mach4    ; }
		
  private:
    float v_tailtot_mach4[17];          // tail-to-tot ratio with 17 different tail values: 30-110ns step 5ns
    float v_gatti_mach4[4]; 		// gatti parameter calculated with 4 reference shapes from 4 different alpha/beta samples
    float f4_peak_mach4;		// peak time of hits
    float f4_mean_mach4;		// mean time of hits
    float f4_rms_mach4;			// root-mean-square of hit time distributio
    float f4_skew_mach4;		// skewness of hit time distribution
    float f4_kurt_mach4;                // kurtosis of hit time distribution

  public:
    typedef std::vector<bx_laben_ab_mach4_cluster> bx_laben_rec_cluster_vector;
  friend class bx_pid_ab_mach4;
};*/

// ---- POSITRON/ELECTRON DISCRIMINATION  ----- //
class bx_laben_positron_cluster : public bx_laben_ab_cluster { // to be modified if mach4 ab reintroduced
  public:
    bx_laben_positron_cluster (const bx_laben_cluster &cluster);// : bx_laben_ab_cluster(cluster);
    virtual ~bx_laben_positron_cluster() {}
	
    float get_gatti_ops_beta     () const { return f4_gatti_ops_beta ; }
    float get_gatti_c11_beta     () const { return f4_gatti_c11_beta ; }
    float get_gatti_ops_nops     () const { return f4_gatti_ops_nops ; }
		
  private:
    float f4_gatti_ops_beta; 	// gatti parameter for otho-positronium/beta discrimination; ops shapes from bxmc2 (tau = 3.1 ns); beta shapes from DATA 214Bi
    float f4_gatti_c11_beta; 	// gatti parameter for c11/beta discrimination; c11 shapes from DATA TFC; beta shapes from 214Bi
    float f4_gatti_ops_nops; 	// gatti parameter for otho-positronium/no-otho-positronium c11 discrimination; ops shapes from bxmc2 (tau = 3.1 ns); nops shapes from bxmc2

  public:
    typedef std::vector<bx_laben_positron_cluster> bx_laben_rec_cluster_vector;
  friend class bx_pid_positron;
};

class bx_laben_rec_event {
  public:
    bx_laben_rec_event () { rec_clusters.clear(); } // ctor
    virtual ~bx_laben_rec_event () {} // dtor

    bx_laben_rec_cluster&      get_rec_cluster      (int i) { return dynamic_cast<bx_laben_rec_cluster&>    (rec_clusters[i]); } // upward getter
    bx_laben_shaped_cluster&   get_shaped_cluster   (int i) { return dynamic_cast<bx_laben_shaped_cluster&> (rec_clusters[i]); } // upward getter
    bx_laben_ab_cluster&       get_ab_cluster       (int i) { return dynamic_cast<bx_laben_ab_cluster&>     (rec_clusters[i]); } // upward getter
//    bx_laben_ab_mach4_cluster& get_ab_mach4_cluster (int i) { return dynamic_cast<bx_laben_ab_mach4&>     (rec_clusters[i]); } // upward getter
    bx_laben_positron_cluster&       get_positron_cluster (int i) { return rec_clusters[i]; } // upward getter
    const bx_laben_rec_cluster&      get_rec_cluster      (int i) const { return dynamic_cast<const bx_laben_rec_cluster&>    (rec_clusters[i]); } // upward getter
    const bx_laben_shaped_cluster&   get_shaped_cluster   (int i) const { return dynamic_cast<const bx_laben_shaped_cluster&> (rec_clusters[i]); } // upward getter
    const bx_laben_ab_cluster&       get_ab_cluster       (int i) const { return dynamic_cast<const bx_laben_ab_cluster&>     (rec_clusters[i]); } // upward getter
//    const bx_laben_ab_mach4_cluster& get_ab_mach4_cluster (int i) const { return dynamic_cast<const bx_laben_ab_mach4&>     (rec_clusters[i]); } // upward getter
    const bx_laben_positron_cluster& get_positron_cluster (int i) const { return rec_clusters[i]; } // upward getter

    int get_nrec_clusters () const { return rec_clusters.size (); }
  private:
    bx_laben_positron_cluster::bx_laben_rec_cluster_vector rec_clusters;

  friend class bx_position_reco;
};

class bx_laben_tracked_event {
  public:
    bx_laben_tracked_event () : b_is_tracked_energy (false), b_is_tracked_tof(false) {} // ctor
    virtual ~bx_laben_tracked_event () {} // dtor

    bool is_tracked_energy () const { return b_is_tracked_energy; }
    bool is_tracked_tof    () const { return b_is_tracked_tof   ; }

    const bx_track& get_track_energy () const { return track_energy; }
    const bx_track& get_track_tof    () const { return track_tof; }
    bx_track_by_points& get_track_energy () { return track_energy; }
    bx_track_by_points& get_track_tof    () { return track_tof; }

  private:
    bool b_is_tracked_energy;
    bool b_is_tracked_tof;
    bx_track_by_points track_energy; 
    bx_track_by_points track_tof; 

  friend class bx_laben_energy_tracker;
  friend class bx_laben_tof_tracker;
};

// ------- derived general laben class --------- //
class bx_laben_event : 
  public bx_base_event,
  public bx_laben_raw_event,
  public bx_laben_decoded_event,
  public bx_laben_clustered_event,
  public bx_laben_rec_event,
  public bx_laben_tracked_event {
public:
  bx_laben_event (const char *disk_event): bx_laben_raw_event(disk_event) { if (get_raw_nhits ()) mark_stage (bx_base_event::raw); }
  virtual ~bx_laben_event () {}
private:
    void mark_stage(event_stage es, bool v = true) { stages[es] = v; }

  friend class bx_laben_decoder;
  friend class bx_laben_findcluster;
  friend class bx_laben_findcluster_time;
  friend class bx_baricentrator;
  friend class bx_position_reco_mi;
  friend class bx_position_reco_lngs;
  friend class bx_position_reco_noavg;
  friend class bx_position_reco_msk;  
  friend class bx_position_reco_dbn;  
  friend class bx_position_reco_mach4;
  friend class bx_position_reco;  
  friend class bx_pid_shape;
  friend class bx_pid_alphabeta;
  friend class bx_pid_positron;
  friend class bx_laben_energy_tracker;
  friend class bx_laben_tof_tracker;
};

#endif
/*
 * $Log: bx_laben_event.hh,v $
 * Revision 1.135  2015/08/26 11:21:44  misiaszek
 * new rms/kurtosis for c11 psd
 *
 * Revision 1.134  2015/08/25 12:25:49  misiaszek
 * mlp_ab variable added to rec_cluster
 *
 * Revision 1.133  2015/08/06 12:33:17  koun
 * added Normalize_geo/QE_Pmts/Charge variables
 *
 * Revision 1.132  2015/07/29 09:23:12  ilia.drachnev
 * added position_reco_noavg variable
 *
 * Revision 1.131  2015/07/28 09:14:54  misiaszek
 * tailtot for c11/b discrimination added
 *
 * Revision 1.130  2015/07/14 09:57:06  misiaszek
 * laben cluster has new friend ;-) - Ilia module for charge_noavg
 *
 * Revision 1.129  2015/07/14 09:42:21  misiaszek
 * f4_charge_noavg_dt1,  f4_charge_noavg_dt2,  f4_charge_noavg , f4_charge_noavg_short added
 *
 * Revision 1.128  2015/07/13 15:07:11  misiaszek
 * tailtot_ab_mlp for a/b discrimination with MLP algorithm added
 *
 * Revision 1.127  2015/01/02 02:13:43  misiaszek
 * npmts, nhits & charge postition corrected energy variables added
 *
 * Revision 1.126  2014/12/15 12:09:42  dekor
 * Added charge_win1 charge_win2 stuff
 *
 * Revision 1.125  2014/12/11 21:27:12  wurm
 * added getters and variables to BxEvent that concern C-11 likelihood tagging
 * * for tracks: new distance-to- and projection-on-track getters, dE/dx calculation
 * * counting of neutron-like clusters in tt128, new flag of neutron clusters
 *
 * Revision 1.124  2014/12/03 17:08:06  misiaszek
 * charge_dt1 and charge_dt2 to clusters added
 *
 * Revision 1.123  2013/01/24 16:43:13  ludhova
 * from:Pablo npmts_dt1/2 moved as cluster variables
 *
 * Revision 1.122  2012-10-30 16:21:48  ddangelo
 * added vectors for npmts calculation in tt64
 * internal and root event with getters and copy.
 *
 * Revision 1.121  2012-10-22 15:56:04  ddangelo
 * added npmts_dt1, npmts_dt2 to laben (clustered) event.
 * internal and root, with getters and cpy operator.
 *
 * Revision 1.120  2011-04-13 10:37:35  ddangelo
 * added variable nclusters_old
 *
 * Revision 1.119  2011-03-23 16:04:48  ddangelo
 * added charge_clean, zeroing, copy. internal and root event
 *
 * Revision 1.118  2011-03-02 17:21:13  razeto
 * Skipping events in reco_mi, default position is 10,10,10
 *
 * Revision 1.117  2011-02-19 13:29:11  davini
 * bx_pid_positron stuffs
 *
 * Revision 1.116  2011-02-18 18:22:02  ddangelo
 * event support to mach4 a/b discrimination commented out. variable writing with the module commented out.
 * added event support for gatti ops. to be improved.
 *
 * Revision 1.115  2011-02-18 16:04:27  davini
 * added npmts_400 nhits_400 charge_400 on Oleg's request; added charge_npmts based on Alessandro, Stefano and Livia idea;
 *
 * Revision 1.114  2011-01-15 20:57:40  razeto
 * Use label for max
 *
 * Revision 1.113  2010-11-26 10:18:29  meindl
 * Fixed empty_boards miscalculation (i.e. changing "nhits_flag[7]" to "nhits_flag[8]")
 *
 * Revision 1.112  2010-07-01 18:21:36  razeto
 * Added counter hit flag for new laben fw and nhits_fw from laben boards
 *
 * Revision 1.111  2010-07-01 13:28:25  ddangelo
 * nhits_bkg changed to float to accomodate info for tt1
 *
 * Revision 1.110  2010-05-21 15:55:14  ddangelo
 * different things on muon tracks
 *
 * 1.a) old laben_track renamed as laben_track_energy
 * new laben_track_tof added
 *
 * 1.b) (global) track renemed as track_global at base event level
 * track_cmt added at base event level (track by points)
 *
 * 1) all getters updated/integrated
 * is_tracked variable updated/integrated accordingly. inizialization.
 * job ported to root event as well. copy done.
 * friendship with old/new module updated
 *
 * 2) bxtrack_by_points class:
 * - theta, phi and impact added as variables.
 * - errors added on all of the above.
 * - error code variable requested by cmt tracker added
 *
 * Revision 1.109  2010-05-21 13:17:15  ddangelo
 * adding laben_tof_tracker and cmt_tracker
 *
 * Revision 1.108  2009-11-18 11:43:31  ddangelo
 * added rms_time, duration and npmts short version for back compatibility.
 * npmt variables renamed to npmts also in internal event.
 * mach4_n1700 renamed as mach4_fixed throughout the event classes.
 *
 * Revision 1.107  2009-10-26 19:17:23  ddangelo
 * added bx_position_mach4_n1700 class (internal and root event, copy and getters)
 * in (base)postion class variable likelihood renamed as user. parallel getters to use it as likelihood or refraction index
 *
 * Revision 1.106  2009-10-23 16:42:22  ddangelo
 * cluster_neutrons renamed as clusters_muons
 *
 * Revision 1.105  2009-10-23 14:00:03  koshio
 * Add the lngs postion reconstruction
 *
 * Revision 1.104  2009-10-23 09:07:16  ddangelo
 * added a laben cluster vector for parallel neutron clustering in the muon gate. empty for non-muon events.
 *
 * Revision 1.103  2009-10-22 16:06:46  ddangelo
 * in bx_laben_cluster added nhits_bkg for evaluate of bkg contribution in case of piled-up events.
 * internal and root event, copy and getters
 *
 * Revision 1.102  2009-10-08 15:45:59  ddangelo
 * implemented pid event variables removal and addition.
 * Internal and root, inizialization, copy and getters.
 *
 * Revision 1.100  2009-10-06 17:05:58  razeto
 * Use unused bit for invalid (0x8000) and mask. Check_flag fixed
 *
 * Revision 1.99  2009-10-06 13:36:17  ddangelo
 * laben decoded event: added variable n_invalid_pmts to account for FE||FF condition only on ordinary && !disabled channels
 * internal and root event, getters and copy
 *
 * Revision 1.98  2009-10-05 15:18:16  ddangelo
 * added npmts at laben decoded event level (inner and outer, copy and getters);
 *
 * Revision 1.97  2009-09-18 22:20:12  ddangelo
 * added sphere_rel_var to laben shaped cluster. Internal and root event, inizialization and copy.
 *
 * Revision 1.96  2009-09-17 14:28:51  ddangelo
 * in laben clustered/decoded (internal/root) hit added the flag short_cluster to say if hit belonged to the old (c11) cluster
 *
 * Revision 1.95  2009-09-17 13:57:47  ddangelo
 * in laben cluster added nhits, charge and mean_time '_short' for old clustering values
 *
 * Revision 1.94  2009-09-16 15:36:52  ddangelo
 * BxLabenCluster::end_time renamed as duration
 * neutron, tags and raw_index branches disabled
 *
 * Revision 1.93  2009-07-31 15:39:50  ddangelo
 * debugging the work of the [previous commit
 *
 * Revision 1.92  2009-07-22 10:41:23  ddangelo
 * in position classes, both internal and root event:
 * - n_iterations removed
 * + matrix and converged variables added with getters.
 *
 * Revision 1.91  2009-07-17 15:39:51  ddangelo
 * laben cluster:
 * + added npmts_thresh, nhits_thresh and charge_thresh, to be computed with >0.2pe hits only (internal and root event).
 * - removed cluster rough time (internal and root).
 * - removed npe and npe_conc (root only)
 *
 * Revision 1.90  2009-07-17 13:08:34  ddangelo
 * alpha beta varaibles reorganized, rise_time added
 *
 * Revision 1.89  2009-07-16 15:50:32  ddangelo
 * copyng data implemented
 *
 * Revision 1.88  2009-07-16 15:17:57  ddangelo
 * mach a/b ported to root event
 * debugging tailtot getter for m4
 * other debugging
 *
 * Revision 1.87  2009-07-16 14:24:03  ddangelo
 * more on mach4 a/b
 *
 * Revision 1.86  2009-07-16 10:50:50  ddangelo
 * new class for mach4 a/b discrimination
 *
 * Revision 1.85  2009-07-16 10:17:42  ddangelo
 * infrastructure for m4 position reco (patch by steve&ale)
 *
 * Revision 1.84  2009-04-20 13:50:40  ddangelo
 * added errors on rec coordinates and time in class bx_position and BxPosition. initialization and copy included.
 *
 * Revision 1.83  2008-12-15 12:12:44  razeto
 * Added rms_time to cluster (to improve PID)
 *
 * Revision 1.82  2008-12-15 11:46:13  razeto
 * Added position ctor
 *
 * Revision 1.81  2008-12-11 17:12:25  razeto
 * Do not remove low level flags (not to break time_decondig)
 *
 * Revision 1.80  2008-12-10 11:40:24  razeto
 * Added mean charge and peak from tt1
 *
 * Revision 1.79  2008-10-20 12:01:31  razeto
 * Fixed a bug: nhits_flag has to be long 7
 *
 * Revision 1.78  2008-10-07 14:03:20  razeto
 * Added invalid flag to event (and removed checker from laben_time_hit)
 *
 * Revision 1.77  2008-10-01 16:17:04  ddangelo
 * removed is_pointlike variables and related stuff (bx_filters will perform this task)
 *
 * Revision 1.76  2008-05-09 15:46:23  ddangelo
 * muon track implemented as "by_points" and "fitted" classes inheriting from common interface (pure virtual class bx_track)
 * laben and muon host an istance "by_points"
 * echidna event hosts an istance "fitted"
 * structure duplicated at root level and filled
 * GetDistance duplicated in BxPosition class
 *
 * Revision 1.75  2008-02-26 18:29:25  ddangelo
 * added is_pointlike/tracklike variable and getters in rec (shaped) cluster
 * both inner and root event
 *
 * Revision 1.74  2008-02-21 15:58:31  ddangelo
 * added muon track as a separate class.
 * Currently used twice: laben and muon (later also globally fit).
 * both inner and roor event
 * added tracked stage in base event enum
 * bx_muon_rec_event renamed as bx_muon_tracked_event
 *
 * Revision 1.73  2007-12-07 14:09:47  ddangelo
 * added n_live_charge (internal and root), getters and normalize method
 * renamed a variable
 *
 * Revision 1.72  2007-11-12 15:19:54  razeto
 * Added cluster flags and removed broad variable (davide auth)
 *
 * Revision 1.71  2007-11-05 23:39:02  razeto
 * Added new hits_on_empty variable to laben data
 *
 * Revision 1.70  2007-10-31 17:12:53  razeto
 * Added a new variable for clustering
 *
 * Revision 1.69  2007-10-30 18:12:27  ddangelo
 * added empty_boards to root event too.
 * getters
 *
 * Revision 1.68  2007-10-30 17:31:21  ddangelo
 * raw hit: out_of_gate condition replaced by a more flexible flag.
 * raw event: empty boards added
 *
 * Revision 1.67  2007-10-30 15:45:03  ddangelo
 * added # of laben cluster found by algorythm (different from saved one for high multiplicity events)
 * added end hit time of a cluster
 * internal and root event. getters, copy, initialization, etc...
 *
 * Revision 1.66  2007-10-26 14:45:23  ddangelo
 * added nphotons for msk energy reco
 *
 * Revision 1.65  2007-10-26 12:04:53  ddangelo
 * supporting the splitting of msk position and energy reco in 2 modules.
 * Both internal and root event
 *
 * Revision 1.64  2007-10-26 09:22:01  ddangelo
 * supporting separation of msk position and energy reco in 2 modules
 *
 * Revision 1.63  2007-10-25 15:46:03  ddangelo
 * added a bit field to shaped event, meaning to be assigned. smart getters to be added at that time. internal and external event.
 * all variables zeroed in ctor for shaped event.
 *
 * Revision 1.61  2007-05-25 15:56:49  ddangelo
 * added npe and charge for laben decode event. Internal and root.
 *
 * Revision 1.60  2007-05-07 16:41:43  ddangelo
 * n_live_channels renamed as n_live_pmts as requested by ale.
 * getters and normalize improved
 *
 * Revision 1.59  2007-05-07 12:49:33  ddangelo
 * cleaned up useless methods
 *
 * Revision 1.58  2007-05-04 16:33:24  pallas
 * Changed variables for alpha beta
 * Now tailtot is an array of 10 values
 *
 * Revision 1.57  2007-05-03 17:25:39  ddangelo
 * added n_live_channels for run-time information. added relative getter.
 * added functions Normalize() for common data types.
 * Both internal and root event
 *
 * Revision 1.56  2007-04-27 14:23:29  pallas
 * Small changes to alpha beta variables
 *
 * Revision 1.55  2007-03-30 12:42:41  razeto
 * rec clusters can be 0, if bx_position_reco is disabled
 *
 * Revision 1.54  2007-03-27 18:19:57  ddangelo
 * completing
 *
 * Revision 1.53  2007-03-27 18:01:29  ddangelo
 * energy class upgraded (inner event only for now, tmp)
 *
 * Revision 1.52  2007-03-27 15:18:20  ddangelo
 * variables npe_conc, charge_conc, nhits_conc added to laben cluster
 * f4_pe renamed as f4_charge in bx_laben_event.hh
 * decoded_charge and decoded_npmts added tu muon event
 *
 * Revision 1.51  2007-03-15 19:59:28  ddangelo
 * completing
 *
 * Revision 1.50  2007-03-15 19:17:19  ddangelo
 * pid event removed.
 * laben event upgraded with classes: bx_laben_shaped_cluster and bx_laben_ab_cluster
 * bx_laben_rec_cluster is now a parent class for the 2 new ones.
 * BxEvent modified accordingly: BxLabenRecHit and BxLabenRecCluster added.
 * BxPidEvent removed.
 *
 * Revision 1.49  2006/10/23 15:34:34  ddangelo
 * applied ale's patch to inlcude laben integrity flags
 *
 * Revision 1.48  2006/10/13 15:31:18  razeto
 * Added variables from raw data
 *
 * Revision 1.47  2006-08-21 17:22:40  ddangelo
 * ndexes to lower level hits removed (undo of previous commit).
 * Indexes for root event are computed in BxEvent::operator=.
 *
 * Revision 1.46  2006/06/29 15:02:44  razeto
 * Added indexes to lower level hits
 *
 * Revision 1.45  2006/05/08 17:31:33  razeto
 * Added db_channel patch for fadc (sent to the mailing list)
 *
 * Revision 1.44  2005/12/30 12:16:43  razeto
 * Added a comparison operator (for sorting hits)
 *
 * Revision 1.43  2005/12/12 19:08:09  razeto
 * Split position and energy reco results (auth from maintainers)
 *
 * Revision 1.42  2005/12/07 21:35:21  misiaszek
 *
 * class bx_energy_reco_lik is friend of bx_base_position (auth from Davide D'Angelo)
 *
 * Revision 1.41  2005/12/03 15:15:54  razeto
 * Added position assigment from mctruth frame (auth from Davide)
 *
 * Revision 1.40  2005/12/01 12:45:46  misiaszek
 *
 * bx_energy_reco_mc is friend class of bx_base_position (auth from Davide D'Angelo)
 *
 * Revision 1.39  2005/11/18 16:42:54  razeto
 * Added new integer charge variable to decoded hit and cluster (auth from davide)
 *
 * Revision 1.38  2005/10/13 13:42:39  razeto
 * Added radius getter (auth from davide)
 *
 * Revision 1.37  2005/10/04 20:20:37  razeto
 * Added iterations value during minimization
 *
 * Revision 1.36  2005/09/26 12:38:26  razeto
 * Introduced likelihood in base position
 *
 * Revision 1.35  2005/09/20 17:17:12  razeto
 * Fixed pid shape event (auth from mantainer)
 *
 * Revision 1.34  2005/07/27 16:58:57  ddangelo
 * added flag for "broad" (i.e. spaparanzeted) clusters
 *
 * Revision 1.33  2005/06/30 15:16:21  ddangelo
 * added a friend for dbn reco in base position
 *
 * Revision 1.32  2005/06/20 16:45:23  ddangelo
 * added bx_position_reco_dbn
 *
 * Revision 1.31  2005/04/29 13:26:37  razeto
 * Fixed a bug, laben times need to be double precision (commit authorized by mantainer)
 *
 * Revision 1.30  2005/03/15 18:08:45  ddangelo
 * added rec stage to laben event.
 * untested
 *
 * Revision 1.29  2005/03/14 13:07:26  ddangelo
 * added "order_in_channel" variable in 3 levels: raw, decoded, clastered.
 * Added computation of the first one in raw event constructor
 *
 * Revision 1.28  2004/12/22 13:48:27  ddangelo
 * debugging
 *
 * Revision 1.27  2004/12/21 18:03:23  ddangelo
 * in bx_split_peak:
 *   nhits mdoified to use a float.
 *   both constructors pushed private
 *
 * Revision 1.26  2004/12/21 16:47:27  ddangelo
 * added bx_peak class for splitting output
 *
 * Revision 1.25  2004/12/15 18:17:09  ddangelo
 * energy moved into bx_base_position class
 *
 * Revision 1.24  2004/12/14 20:01:25  ddangelo
 * base class bx_base_postion added.
 * classes bx_baricenter, bx_position_mi, bx_position_msk added at bx_laben_cluster level.
 * baricenter stuff moved into proper class.
 * friend declarations added where needed for modules bx_position_reco_mi, bx_position_reco_msk.
 *
 * Revision 1.23  2004/11/26 15:25:10  razeto
 * Thanks to Davide who so kindly underlined a syntactical error, I fixed the "Maintainer" word
 *
 * Revision 1.22  2004/11/26 14:06:22  razeto
 * Added Mantainer field
 *
 * Revision 1.21  2004/10/21 11:24:49  ddangelo
 * added friend class bx_baricentrator (sigh!)
 *
 * Revision 1.20  2004/09/22 10:45:39  ddangelo
 * added bx_baricentrator support (by Ale)
 *
 * Revision 1.19  2004/09/13 12:13:06  pallas
 * Changed the time references hits handling (by Alessandro)
 *
 * Revision 1.18  2004/08/31 13:29:03  ddangelo
 * added charge to laben hit
 *
 * Revision 1.17  2004/07/24 16:24:21  ddangelo
 * debugging (patch by Alessandro)
 *
 * Revision 1.16  2004/07/09 13:53:53  ddangelo
 * added friend declaration in bx_laben_event to allow mark stage by clustering modules
 *
 * Revision 1.15  2004/07/07 12:56:18  ddangelo
 * laben cluster and relative hit classes filled with meaningful variables (by Daniela)
 *
 * Revision 1.14  2004/06/07 17:17:09  ddangelo
 * b_is_valid removed from decoded hit.
 * time renamed to raw_time. getter renamed accordingly.
 *
 * Revision 1.13  2004/06/07 11:06:11  ddangelo
 * added operator< for decoded hit (by Alessandro)
 *
 * Revision 1.12  2004/06/01 18:23:40  ddangelo
 * applied patch by Alessandro
 * in bx_laben_decoded_event:
 * vectors of new private sub class time_reference_elements
 * + other upgrades
 *
 * Revision 1.11  2004/05/31 16:48:01  ddangelo
 * non-english name 'clusterized' replaced by 'clustered'
 *
 * Revision 1.10  2004/05/31 14:16:04  ddangelo
 * added clusterized level classes (still empty)
 *
 * Revision 1.9  2004/05/31 13:39:16  ddangelo
 * applied patch with 2 bx_laben_decoded_hit::operator=() (by Alessandro)
 *
 * Revision 1.8  2004/05/27 14:54:14  ddangelo
 *  removed useless friend declarations
 *
 * Revision 1.7  2004/05/27 14:52:30  ddangelo
 * p_db_channel changed to const
 *
 * Revision 1.6  2004/05/25 16:37:05  ddangelo
 * minor updates
 *
 * Revision 1.1  2004/04/27 09:41:09  ddangelo
 * Event structure redisegned.
 * bx_reco_event now handles event header and has object for different detector portions.
 * Inheritance from different reconstruction status is kept within each detector portion obj.
 * Hit lists are linked by const ptr (maybe a temporary solution).
 *
 *
 */
