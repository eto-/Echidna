/* BOREXINO Reconstruction program
 *
 * Author: Davide D'Angelo <davide.dangelo@lngs.infn.it>
 * Maintainer: Davide D'Angelo <davide.dangelo@lngs.infn.it>
 *
 * $Id: bx_track.hh,v 1.14 2014/12/11 21:27:12 wurm Exp $
 *
 * A track class to be used by laben, muon and global event
 * 
 */
#ifndef _BX_TRACK_HH
#define _BX_TRACK_HH

//#include <vector>

#include "bx_rec_general.hh"
#include "bx_base_event.hh"
#include "constants.hh"


class bx_track {
  public:
  //  bx_track () {}
    virtual ~bx_track () { }

    virtual float get_pedal_t () const = 0;
    virtual float get_pedal_x () const = 0;
    virtual float get_pedal_y () const = 0; 
    virtual float get_pedal_z () const = 0;

    virtual float get_theta   () const = 0;
    virtual float get_phi     () const = 0;
    virtual float get_impact  () const { return ::sqrtf ( pow(get_pedal_x(), 2) + pow(get_pedal_y(), 2) + pow(get_pedal_z(), 2) ); }
    
    virtual bool  is_upward           () const = 0;
    virtual bool  is_downward         () const = 0;

};

class bx_track_by_points : public bx_track {
  public:
    bx_track_by_points ();
    virtual ~bx_track_by_points () { }

    float get_x1 () const { return f4_x1 ; }
    float get_y1 () const { return f4_y1 ; }
    float get_z1 () const { return f4_z1 ; }
    float get_x2 () const { return f4_x2 ; }
    float get_y2 () const { return f4_y2 ; }
    float get_z2 () const { return f4_z2 ; }
    float get_t1 () const { return f4_t1 ; }
    float get_t2 () const { return f4_t2 ; }
    float get_dx1() const { return f4_dx1; }
    float get_dy1() const { return f4_dy1; }
    float get_dz1() const { return f4_dz1; }
    float get_dx2() const { return f4_dx2; }
    float get_dy2() const { return f4_dy2; }
    float get_dz2() const { return f4_dz2; }

    virtual float get_pedal_t () const { return f4_t1 + (f4_t2-f4_t1) * get_OA() * get_cos_alpha() / get_AB(); }
    virtual float get_pedal_x () const { return f4_x1 + (f4_x2-f4_x1) * get_OA() * get_cos_alpha() / get_AB(); }
    virtual float get_pedal_y () const { return f4_y1 + (f4_y2-f4_y1) * get_OA() * get_cos_alpha() / get_AB(); }
    virtual float get_pedal_z () const { return f4_z1 + (f4_z2-f4_z1) * get_OA() * get_cos_alpha() / get_AB(); }

    virtual float get_theta   () const { return f4_theta  ; }
    virtual float get_phi     () const { return f4_phi    ; }
    virtual float get_dtheta  () const { return f4_dtheta ; }
    virtual float get_dphi    () const { return f4_dphi   ; }
    virtual float get_impact  () const { return f4_impact ; }
    virtual float get_dimpact () const { return f4_dimpact; } 
    virtual float get_laben_normhits () const { return f4_labennormhits; }
    int get_error () const { return i4_error; }
//    virtual float get_theta   () const { return ::acos ( (f4_z1-f4_z2) / get_AB()); }
//    virtual float get_phi     () const { float phi = ::acos ( (f4_x2-f4_x1) / ::sqrtf( pow((f4_x2-f4_x1), 2) + pow((f4_y2-f4_y1), 2) ) );
//                                 if (f4_y2 < f4_y1) phi = 2*constants::number::pi - phi;
//				 if (is_downward()) phi = constants::number::pi + phi; 
//				 return (phi>2*constants::number::pi) ? phi - 2*constants::number::pi : phi;
//			       }
    virtual bool  is_upward           () const { return !b_downward; }
    virtual bool  is_downward         () const { return b_downward ; }

  private:
    float get_OA        () const { return ::sqrtf( f4_x1*f4_x1 + f4_y1*f4_y1 + f4_z1 * f4_z1 ); }
    float get_AB        () const { return ::sqrtf( pow((f4_x1-f4_x2), 2) + pow((f4_y1-f4_y2), 2) + pow((f4_z1-f4_z2), 2) ); }
    float get_cos_alpha () const { return (-f4_x1*(f4_x2-f4_x1)-f4_y1*(f4_y2-f4_y1)-f4_z1*(f4_z2-f4_z1))/(get_AB()*get_OA()); }

    float f4_t1, f4_t2;
    float f4_x1, f4_y1, f4_z1;
    float f4_x2, f4_y2, f4_z2;
    float f4_dx1, f4_dy1, f4_dz1;
    float f4_dx2, f4_dy2, f4_dz2;
    float f4_theta, f4_phi;
    float f4_dtheta, f4_dphi;
    float f4_impact, f4_dimpact;
    float f4_labennormhits;
    int   i4_error;
    bool  b_downward;

  friend class bx_muon_tracker;
  friend class bx_laben_energy_tracker;
  friend class bx_laben_tof_tracker;
  friend class bx_cmt_tracker;
};


class bx_track_fitted : public bx_track {
  public:
    bx_track_fitted ();
    virtual ~bx_track_fitted () { }

    double get_alpha       () const { return f8_alpha      ; }
    double get_beta        () const { return f8_beta       ; }
    double get_gamma       () const { return f8_gamma      ; }
    double get_delta       () const { return f8_delta      ; }
    double get_alpha_error () const { return f8_alpha_error; }
    double get_beta_error  () const { return f8_beta_error ; }
    double get_gamma_error () const { return f8_gamma_error; }
    double get_delta_error () const { return f8_delta_error; }
    float  get_chi2        () const { return f4_chi2       ; }
    int    get_points      () const { return u1_points     ; }

    virtual float get_pedal_t () const { return 0.; }
    virtual float get_pedal_x () const { return -(f8_alpha*f8_beta+f8_gamma*f8_delta)/(1.+f8_beta*f8_beta+f8_delta*f8_delta); }
    virtual float get_pedal_y () const { return f8_alpha+f8_beta*get_pedal_x(); }
    virtual float get_pedal_z () const { return f8_gamma+f8_delta*get_pedal_x(); }

    /*virtual float get_phi     () const {
      float phi = atan2(f8_beta, 1.);
      if ( (f8_delta<0. && b_downward) || (f8_delta>0. && !b_downward) ) phi += constants::number::pi;
      while ( phi<0. ) { phi += 2.*constants::number::pi; }
      while ( phi>2.*constants::number::pi ) { phi -= 2.*constants::number::pi; }
      return phi;
    }
    virtual float get_theta   () const {
      float theta = acos(f8_delta/sqrt(1.+f8_beta*f8_beta+f8_delta*f8_delta));
      if ( theta<0. ) theta = constants::number::pi+theta;
      if ( (f8_delta<0. && b_downward) || (f8_delta>0. && !b_downward) ) theta = constants::number::pi-theta;
      return theta;
    }*/
    virtual float get_theta   () const { return f4_theta  ; }
    virtual float get_phi     () const { return f4_phi    ; }
    virtual float get_dtheta  () const { return f4_dtheta ; }
    virtual float get_dphi    () const { return f4_dphi   ; }
    virtual float get_impact  () const { return f4_impact ; }
    virtual float get_dimpact () const { return f4_dimpact; } 
   
    virtual bool  is_upward           () const { return !b_downward; }
    virtual bool  is_downward         () const { return b_downward ; }
    virtual float get_laben_normhits  () const { return f4_labennormhits; }

    bool  is_valid                    () const { return b_is_valid ; }

  private:
    float get_OA        () const { return ::sqrtf( 100 + pow(f8_alpha + f8_beta*10,2) + pow(f8_gamma + f8_delta*10,2) ); }
    float get_AB        () const { return ::sqrtf( pow(20, 2) + pow(f8_beta*20, 2) + pow(f8_delta*20, 2) ); }
    float get_cos_omega () const { return ( 10*(-20) - (f8_alpha + f8_beta*10)*(-f8_beta*20) - (f8_gamma + f8_delta*10)*(-f8_delta*20) ) / (get_AB()*get_OA()); }

    double f8_alpha, f8_beta, f8_gamma, f8_delta;   
    double f8_alpha_error, f8_beta_error, f8_gamma_error, f8_delta_error;
    float  f4_chi2;
    float  f4_phi, f4_theta, f4_impact, f4_dphi, f4_dtheta, f4_dimpact;
    float  f4_labennormhits;
    bool   b_downward;
    char   u1_points;
    bool   b_is_valid;

  friend class bx_global_tracker;
};
#endif
/*
 * $Log: bx_track.hh,v $
 * Revision 1.14  2014/12/11 21:27:12  wurm
 * added getters and variables to BxEvent that concern C-11 likelihood tagging
 * * for tracks: new distance-to- and projection-on-track getters, dE/dx calculation
 * * counting of neutron-like clusters in tt128, new flag of neutron clusters
 *
 * Revision 1.13  2010/08/03 15:58:05  wurm
 * introduced theta, phi, impact variables for fitted tracks
 *
 * Revision 1.12  2010-05-21 15:55:14  ddangelo
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
 * Revision 1.11  2010-05-21 13:17:15  ddangelo
 * adding laben_tof_tracker and cmt_tracker
 *
 * Revision 1.10  2009-07-31 15:39:50  ddangelo
 * debugging the work of the [previous commit
 *
 * Revision 1.9  2009-04-15 17:12:38  ddangelo
 * n_point changed into a bitfield, internal and root event (BxTrackFitted class)
 *
 * Revision 1.8  2008-12-11 08:02:46  wurm
 * simplified get_theta for track by points
 *
 * Revision 1.7  2008-12-02 13:53:32  ddangelo
 * added n_points to fitted track class and getters. internal and root event
 *
 * Revision 1.6  2008-10-01 16:42:23  wurm
 * changed theta and phi getters for fitted track
 *
 * Revision 1.5  2008-08-05 15:30:51  ddangelo
 * added chi2
 *
 * Revision 1.4  2008-07-17 15:55:00  ddangelo
 * added track direction, removed unused getter
 *
 * Revision 1.3  2008-07-11 14:58:54  ddangelo
 * added a few functions (by werner)
 *
 * Revision 1.2  2008-05-09 15:46:23  ddangelo
 * muon track implemented as "by_points" and "fitted" classes inheriting from common interface (pure virtual class bx_track)
 * laben and muon host an istance "by_points"
 * echidna event hosts an istance "fitted"
 * structure duplicated at root level and filled
 * GetDistance duplicated in BxPosition class
 *
 * Revision 1.1  2008-02-21 15:58:31  ddangelo
 * added muon track as a separate class.
 * Currently used twice: laben and muon (later also globally fit).
 * both inner and roor event
 * added tracked stage in base event enum
 * bx_muon_rec_event renamed as bx_muon_tracked_event
 *
 */
