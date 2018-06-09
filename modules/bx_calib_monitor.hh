/* BOREXINO Reconstruction program
 *
 * Author: Elena Guardincerri <elena.guardincerri@ge.infn.it>
 * Maintainer: Elena Guardincerri <elena.guardincerri@ge.infn.it>
 *
 * $Id: bx_calib_monitor.hh,v 1.10 2009/03/30 08:49:51 guardi Exp $
 *
 */

#ifndef _BX_CALIB_MONITOR_HH
#define _BX_CALIB_MONUOTR_HH

#include "bx_base_module.hh"
#include "bx_rec_general.hh"

class bx_laben_event;
class TH1F;
class TH2F;

class bx_calib_monitor: public bx_base_module {
  public:
  
    bx_calib_monitor ();
    virtual ~bx_calib_monitor ();

    virtual void begin ();
    virtual bx_echidna_event* doit (bx_echidna_event *ev);
    virtual void end ();

  private:

    void m_send_and_reset_histos();
    void m_delete_histos();
    double m_fit2D(double&, double&, TH2F *);

  // histograms
  //
  // these 1-D histos are filled from a fit to the TH2F below
    TH1F* neutrino_z_t_14C;
    TH1F* neutrino_z_t_Bi;
    TH1F* neutrino_x2y2_t_14C;
    TH1F* neutrino_x2y2_t_Bi;
    TH1F* neutrino_x_t_14C;
    TH1F* neutrino_x_t_Bi;
    TH1F* neutrino_y_t_14C;
    TH1F* neutrino_y_t_Bi;
    TH1F* neutrino_x_14C_c;
    TH1F* neutrino_x_Bi_c;
    TH1F* neutrino_x_HE_c;
    TH1F* neutrino_x_VHE_c;
    TH1F* neutrino_y_14C_c;
    TH1F* neutrino_y_Bi_c;
    TH1F* neutrino_y_HE_c;
    TH1F* neutrino_y_VHE_c;
    TH1F* neutrino_z_14C_c;
    TH1F* neutrino_z_Bi_c;
    TH1F* neutrino_z_HE_c;
    TH1F* neutrino_z_VHE_c;
    TH1F* neutrino_E_spectrum_c;
    TH1F* rate_buffer_up;
    TH1F* rate_buffer_equator;
    TH1F* rate_buffer_down;
    

    TH1F* energy_14C;
    TH1F* energy_Bi;
    TH1F* energy_HE;

    TH1F* event_rate;

    TH1F* hits_t_trigger_t;
    TH1F* start_t_trigger_t;

    // these 2-D histos are produced every buffer_time seconds
    TH2F* neutrino_z_x2y2_14C;
    TH2F* neutrino_z_x2y2_Bi;
    TH2F* neutrino_z_x2y2_HE;
    TH2F* neutrino_z_x2y2_VHE;
    TH2F* neutrino_x_y_14C;
    TH2F* neutrino_x_y_Bi;
    TH2F* neutrino_x_y_HE;
    TH2F* neutrino_x_y_VHE;

    //cumulative 2D histograms
    TH2F* neutrino_z_x2y2_14C_c;
    TH2F* neutrino_z_x2y2_Bi_c;
    TH2F* neutrino_z_x2y2_HE_c;
    TH2F* neutrino_z_x2y2_VHE_c;
    TH2F* neutrino_x_y_14C_c;
    TH2F* neutrino_x_y_Bi_c;
    TH2F* neutrino_x_y_HE_c;
    TH2F* neutrino_x_y_VHE_c;

    // lists for handling the TObjects in this class
    TList* f_temp_histos_list;
    TList* f_run_histos_list;

    //time after which the 2-D histograms are redrawn (in seconds)
    unsigned long f_buffer_time; 
    double f_prev_bi_event_time;
    unsigned long f_prev_bunch_time;
    double f_old_x;
    double f_old_y;
    double f_old_z;
    double f_old_r;
    double f_old_x_y;
    double f_Bi_E;
    double f_14C_E;
    double f_HE_E;
    double f_start_time;
    int f_n_events;
    int f_n_events_buffer_up;
    int f_n_events_buffer_equator;
    int f_n_events_buffer_down;
    double f_total_live_time;
    double f_begin_new_bunch;
    double f_end_new_bunch;


};

#endif
/*
 * $Log: bx_calib_monitor.hh,v $
 * Revision 1.10  2009/03/30 08:49:51  guardi
 * plots for event rate in the buffer (upper part, equator, bottom ) added and more non cumulative plots removed
 *
 * Revision 1.9  2009-03-30 08:17:37  guardi
 * non cumulative histograms removed
 *
 * Revision 1.8  2009-03-23 10:26:19  guardi
 *  histograms for hits_time - trigger_time and start_time - trigger_time added
 *
 * Revision 1.7  2008-10-07 17:40:30  guardi
 * added event rate histogram
 *
 * Revision 1.6  2008-10-04 10:27:10  guardi
 * added full spectrum and plots for the position of >1000 p.e. events
 *
 * Revision 1.5  2008-10-04 08:20:04  guardi
 * added cumulative 1D plots and debug
 *
 * Revision 1.4  2008-10-04 07:33:45  guardi
 * added 1D position distributions
 *
 * Revision 1.3  2008-10-03 08:18:51  razeto
 * Upgraded (by elena)
 *
 * Revision 1.2  2008-10-01 16:53:28  guardi
 * charge spectra and axis titles added
 *
 * Revision 1.1  2008-10-01 15:55:43  guardi
 * module to plot the position of the events online added. Purpose: determine the position of the source during calibration
 *
 *
 */
