/* BOREXINO Reconstruction program
 *
 * Author: Alessandro Fontana <alessandro.fontana@mi.infn.it>
 * Maintainer: Alessandro Fontana <alessandro.fontana@mi.infn.it>
 *
 * $Id: bx_calib_laser_transparency.hh,v 1.3 2006/08/23 13:10:23 fontana Exp $
 *
 * 
 * 
 */
#ifndef _BX_CALIB_LASER_TRANSPARENCY_H
#define _BX_CALIB_LASER_TRANSPARENCY_H
#include "bx_base_module.hh"
#include "bx_rec_general.hh"
#include "bx_echidna_event.hh"
#include "TH1F.h"
#include "TH2F.h"
#include "TSpectrum.h"
#define pi 3.141592654


class bx_calib_laser_transparency: public bx_base_module {
  public:
    
    bx_calib_laser_transparency ();
    virtual ~bx_calib_laser_transparency () {}

    virtual void begin ();
    virtual bx_echidna_event* doit (bx_echidna_event *ev);
    virtual void end ();
 
  private:
   
    int32_t counter_ev;
    int32_t processed_events;
    int32_t time_calib_events;
    float refidx;
    float threshold;
    float peak_width;
    int32_t resolution;
    TH1F *time_dis;
    TH1F *time_dis_calib;
    TH1F *spot_distribution;
    TSpectrum *time_calib;
    TH2F* time_vs_angle;
    TH2F* frontpro;
    TH2F* backpro;
    TH2F* top;
    TH2F* bottom;
    TH1F* charge_dis;
    TH1F* raw_charge_dis;
    TH1F* lg_dis;
    double x_bar;
    double y_bar;
    double N_bar;
    double distance;
    int32_t ich, nhits, xx, yy;
    int32_t Nbins;
    double dir_peak_ev, dir_peak_tot;
    int32_t lg;
    //int32_t npe;
    float raw_charge;
    double charge;
    double abs_range;
    double A;  //meters-bins conversion (range is [-abs_range, abs_range] meters)
    double B;
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,0,0)
    double * peak_time;
#else
    float * peak_time;
#endif
    float zero_charge_dir_peak;
    float counter_charge_dir_peak;
    //Vectors of source and target positions. For target logical channel is used
    const static float buffer_source_phi[19];
    const static float buffer_source_theta[19];
    const static int32_t buffer_target[19];
    const static float radial_source_phi[12];
    const static float radial_source_theta[12];
    int32_t type, patch_panel, target, fbflag;
    double delay, alpha, omega1, omega2, norma_r, theta_src, phi_src, nominal_target;
    double x_pmt, y_pmt, z_pmt, x_r, y_r, z_r, x_p, y_p, z_p, x_tb, y_tb, mod1, mod2, norma_pmt, norma_tb, tmp2;
    double x_t, y_t, z_t, R;
    float x_s, y_s, z_s;
    double STe, CTe, SPe, CPe;
    double dir_peak_time, time_sup, time_inf;
    
  

};

#endif
/*
 *
 */
