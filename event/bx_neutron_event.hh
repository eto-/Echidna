/* BOREXINO Reconstruction program
 *
 * Author: Stefano Davini <stefano.davini@ge.infn.it>
 * Maintainer: Davide D'Angelo <davide.dangelo@mi.infn.it>
 *
 * $Id: bx_neutron_event.hh,v 1.3 2008/08/26 15:30:30 ddangelo Exp $
 *
 * The neutron event object (data from the v1731 system dedicated to neutron detection)
 * 
 */

#ifndef _BX_NEUTRON_EVENT_HH
#define _BX_NEUTRON_EVENT_HH

#include <vector>

#include "bx_rec_general.hh"
#include "bx_base_event.hh"
#include "constants.hh"

class bx_neutron_pulse {
  public:
    bx_neutron_pulse();
    virtual ~bx_neutron_pulse() {}

    float get_charge    () const { return f4_charge; }
    float get_amplitude () const { return f4_amplitude  ; }
    float get_peak_time () const { return f4_peak_time  ; }
    float get_rise_time () const { return f4_rise_time ; }
    float get_fall_time () const { return f4_fall_time ; }
    float get_x      () const { return f4_x     ; }
    float get_y      () const { return f4_y     ; }
    float get_z      () const { return f4_z     ; }
    float get_dx     () const { return f4_dx    ; }
    float get_dy     () const { return f4_dy    ; }
    float get_dz     () const { return f4_dz    ; }
  
  private:
    float f4_charge;    // charge as pulse area
    float f4_amplitude; // charge as peak amplitude
    float f4_peak_time;      // peak time position in ns after muon time
    float f4_rise_time;     // peak standard deviation
    float f4_fall_time;     // from the exponential
    float f4_x;         // position of pulse x
    float f4_y;
    float f4_z; 
    float f4_dx;
    float f4_dy;
    float f4_dz;
  
  friend class bx_v1731sys;
};


class bx_neutron_event {
  public:
    bx_neutron_event ();
    virtual ~bx_neutron_event () {}
  
    // getters
    bool is_enabled      () const { return b_is_enabled;  }
    bool is_associated   () const { return b_is_associated;  }
    int  get_n_pulses    () const { return pulses.size();   }
    int  get_n_neutrons  () const { return i4_n_neutrons; }
    const std::vector<bx_neutron_pulse>& get_neutron_pulses  () const { return pulses; }
    const bx_neutron_pulse& get_neutron_pulse(int i) const {return pulses.at(i);}
  
  private:  
    bool b_is_enabled;                     // File found
    bool b_is_associated;                  // trgid/time associated
    int i4_n_neutrons;                     // number of neutron candidates
    std::vector<bx_neutron_pulse> pulses;  // list of candidate pulses in V1731 boards
  
  friend class bx_v1731sys;
};

#endif
/*
 * $Log: bx_neutron_event.hh,v $
 * Revision 1.3  2008/08/26 15:30:30  ddangelo
 * added neutron enabled and association flag, muon aligned flags
 *
 * Revision 1.2  2008-08-05 16:30:52  ddangelo
 * added fall time, some variables renamed
 *
 * Revision 1.1  2008-07-11 17:05:45  ddangelo
 * added
 * 
 *
 */
