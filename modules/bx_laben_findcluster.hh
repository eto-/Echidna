/* BOREXINO Reconstruction program
 *
 * Author: Daniela Manuzio <dmanuzio@ge.infn.it>
 * Maintainer: Alessandro Razeto <Alessandro.Razeto@ge.infn.it>
 *
 * $Id: bx_laben_findcluster.hh,v 1.31 2015/07/20 15:23:30 ilia.drachnev Exp $
 *
 * Definition of bx_laben_findcluster
 * Module that finds the candidate fragments in the laben hit 
 * collection. No real splitting is done in this module.
 * Input: the list of time ordered laben hits
 * Output: a list of bx_cluster objects
 * 
*/

#ifndef _BX_LABEN_FINDCLUSTER_HH
#define _BX_LABEN_FINDCLUSTER_HH

#include "bx_rec_general.hh"
#include "bx_base_module.hh"
#include <vector>
#include "TROOT.h"
#include "TH2F.h"
#include "TH1F.h"
#include "db_run.hh"

class bx_echidna_event;
class bx_laben_cluster;
class bx_laben_event;

class bx_laben_findcluster : public bx_base_module {
  public:
    bx_laben_findcluster ();
    virtual ~bx_laben_findcluster () {}
      
    // Operations for the framework 
    virtual void begin ();
    virtual bx_echidna_event* doit (bx_echidna_event *ev);
    virtual void end ();

  private:
    // parameter cache
    int i4_start_threshold;
    int i4_neutron_start_threshold;
    int i4_count_threshold;
    int i4_he_threshold;
    int i4_ripple_bins;
    int i4_he_ripple_bins;
    int i4_tail_bins;
    int i4_long_bins;
    int i4_enable_histos;
    bool b_force_one_large_cluster;
    bool b_strict_gate;
	bool b_mach4_cluster;
    int i4_gate_start, i4_gate_end, i4_cluster_offset;

    int i4_ripple_count;
    int i4_he_ripple_count;
	float f_dark_rate;
  static const int i4_dt1_len_=230;
  static const int i4_dt2_len_=400;
  int i4_NFrames_win1_, i4_NFrames_win2_;
  float f4_trigger_start_;
  /*int* pablo_npmts_dt1;
    int* pablo_npmts_dt2;*/

      // Histograms
    TH2F *clusters_raw_times;
    TH2F *clusters_mean_times;
    TH2F *clusters_times;
    TH2F *clusters_hits;
    TH2F *clusters_npe;
    TH2F *clusters_pmts;
    TH1F* clusters_number;
    TH2F* histo_raw;
    TH2F* histo_cluster;
    
      // Internal methods
    void findcluster (bx_echidna_event *ev);
    void clusterize (int start_bin, int end_bin, const bx_echidna_event *ev);
    void clusterize_mach4 (int start_bin, int end_bin, const bx_echidna_event *ev);
    void clusterize_neutrons (int start_bin, int end_bin, const bx_echidna_event *ev);
    void clusterize_neutrons_in_muongate (int start_bin, int end_bin, const bx_echidna_event *ev);
    void findcluster_time (bx_echidna_event *ev);
    int find_start_hit (const bx_laben_cluster& cluster_ref, float& time_separation, int intervals);
    int evaluate_ripple_count (float dark_rate, int ripple_bins, int count_threshold, const char* message);
    int get_gate_index_ (const int gate_length, const Float_t dt);

    // vector with 16 ns binned times 
    int i4_time_bins;
    short int *binned_times;
    unsigned long prev_gps_times[2];

    struct cluster_data {
      cluster_data (): i4_short_end_bin(-1) {}
      int i4_start_bin, i4_long_end_bin, i4_short_end_bin, i4_n_hits;
      float f4_n_hits_background;
      bool bool_muon_clustered, b_is_neutron;
    };
    std::vector<cluster_data> clusters;
    int* fired_channel;
//**noavg code::db_run instance**
    db_run* run_info;
//*******************************
};

#endif
/*
 * $Log: bx_laben_findcluster.hh,v $
 * Revision 1.31  2015/07/20 15:23:30  ilia.drachnev
 * noavg charges added
 *
 * Revision 1.30  2014/12/11 21:27:12  wurm
 * added getters and variables to BxEvent that concern C-11 likelihood tagging
 * * for tracks: new distance-to- and projection-on-track getters, dE/dx calculation
 * * counting of neutron-like clusters in tt128, new flag of neutron clusters
 *
 * Revision 1.29  2013/06/05 15:16:02  mosteiro
 * stylistic shift
 *
 * Revision 1.28  2013-01-22 09:22:17  mosteiro
 * dt1,2 vars calculated for each cluster
 *
 * Revision 1.27  2012-11-09 04:47:31  mosteiro
 * added private member variables and functions needed for filling npmts_win1 and npmts_win2
 *
 * Revision 1.26  2012-11-01 02:12:18  mosteiro
 * Added lengths of npmts_dt1 and npmts_dt2 variables
 *
 * Revision 1.25  2011-04-01 12:23:30  meindl
 * New function "clusterize_neutrons_in_muongate" :
 * Uses cycle13 neutron clustering for neutron detection inside the muon gate
 *
 * Revision 1.24  2010-07-01 16:26:43  meindl
 * And the corresponding header file.
 *
 * Revision 1.23  2009-10-23 16:49:42  meindl
 * And the corresponding header-file
 *
 * Revision 1.22  2009-10-22 19:01:12  meindl
 * modified structure "cluster_data" to include background-hits
 *
 * Revision 1.21  2009-09-16 13:20:06  razeto
 * Added short/long clustering marking, waiting for event variable commit
 *
 * Revision 1.20  2009-09-16 11:16:10  razeto
 * Long cluster (aka m4) is now a parameter (defaulting at 1.5us)
 *
 * Revision 1.19  2009-08-04 16:01:19  alvaro
 * Updated dark_rate value and fixed unsigned-signed comparison warning
 *
 * Revision 1.17  2008-10-06 17:02:51  razeto
 * reverted to previous status (no rate limit)
 *
 * Revision 1.16  2008-10-05 13:30:45  razeto
 * Added rate_limit for calibrations and online echidna
 *
 * Revision 1.15  2008-02-15 10:45:58  razeto
 * clusterize internal method created, for more modularity. Added clusterize_neutrons (empty)
 *
 * Revision 1.14  2007-11-12 15:20:30  razeto
 * Filling new flag variable
 *
 * Revision 1.13  2007-11-12 09:51:58  razeto
 * Use int as specified in name
 *
 * Revision 1.12  2007-11-05 23:42:28  razeto
 * Code more stable
 *
 * Revision 1.11  2007-10-31 17:12:39  razeto
 * Code debugged for normal events; still in beta phase for neutrons
 *
 * Revision 1.10  2007-10-30 18:06:47  razeto
 * Experimental code: work in gateless mode for events close to the previous.
 *
 * Revision 1.9  2007-04-15 16:53:31  razeto
 * Added high energy mode to handle 140ns dead time
 *
 * Revision 1.8  2007-04-15 14:19:11  razeto
 * New end of cluster
 *
 * Revision 1.7  2005/09/19 11:41:29  razeto
 * Added an option to teste the splitting
 *
 * Revision 1.6  2005/07/15 13:51:10  razeto
 * Updated to use a new algorithm
 *
 * Revision 1.5  2005/07/13 12:32:08  razeto
 * Merged the the 2 laben clustering modules: now find_cluster_time is
 * a subroutine of the bx_laben_findcluster module
 *
 * Revision 1.4  2005/05/29 16:48:59  razeto
 * Changed algorithm
 *
 * Revision 1.3  2004/11/29 13:21:23  razeto
 * Added Mantainer field
 *
 * Revision 1.2  2004/09/22 14:20:47  dmanuzio
 * Changed bx_reco_event to bx_echidna_event
 *
 * Revision 1.1  2004/07/09 13:59:08  dmanuzio
 * Added 2 modules for Laben clustering
 *
 *
 */
