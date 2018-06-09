/* BOREXINO Reconstruction program
 *
 * Authors: Davide D'Angelo <Davide D'Angelo@mi.infn.it>, Micheal Wurm <mwurm@ph.tum.de
 * Maintainer: Davide D'Angelo <Davide D'Angelo@mi.infn.it>
 *
 * $Id: bx_muon_findcluster.hh,v 1.11 2009/11/08 15:35:06 wurm Exp $
 *
 * Definition of bx_muon_findcluster
 * Module that finds the candidate clusters in the muon hit 
 * collection. 
 * Input: the list of time ordered muon hits
 * Output: a list of bx_muon_cluster and bx_muon_clustered_hits objects
 * 
*/

#ifndef _BX_MUON_FINDCLUSTER_HH
#define _BX_MUON_FINDCLUSTER_HH

#include "bx_rec_general.hh"
#include "bx_base_module.hh"
#include <vector>
#include "TH1F.h"

class bx_echidna_event;
class bx_muon_clustered_hit;

class bx_muon_findcluster : public bx_base_module {
  public:
    bx_muon_findcluster ();
    virtual ~bx_muon_findcluster () {}
      
    // Operations for the framework 
    virtual void begin ();
    virtual bx_echidna_event* doit (bx_echidna_event *ev);
    virtual void end ();

  private:
    int  m_affiliate (int offset, std::vector<bx_muon_clustered_hit>&, float maxdt, float maxdr);
    bool m_split (int id_to_be_split);

    // binned distribution params
    int   i4_bin_width;
    int   i4_clustered_events;
    int   i4_start_threshold;
    int   i4_count_threshold;
    int   i4_ripple_count;
    int   i4_enable_histos;

    // binned distribution derived quantities
    int   i4_n_bins;
    int   i4_ripple_bins;

    // affiliation params
    float f4_mincharge;
    float f4_maxdt;
    float f4_maxdr_sss;
    float f4_maxdr_floor;
		float f4_hit_charge_threshold;    // minimum charge a hit has to have to be clustered

    // clustering weighting params 
    float f4_tau;
    int   i4_max_hits_sss;
    int   i4_max_hits_floor;

    // splitting params
    int   i4_split_minimum_hits;
    float f4_wimp_fraction;
    float f4_split_separation;
    float f4_split_window;
    float f4_ds_splitting;
    float f4_dt_splitting;

    // Histograms
    TH1F *clusters_start_times;
    TH1F *clusters_nhits;
    TH1F *clusters_charge;
    TH1F *clusters_npmts;

    // vector with 16 ns binned times 
    std::vector<int> binned_times_sss;
    std::vector<int> binned_times_floor;

    std::vector<bx_muon_clustered_hit> chits_sss;
    std::vector<bx_muon_clustered_hit> chits_floor;

    int* fired_channels;
};

#endif
/*
 * $Log: bx_muon_findcluster.hh,v $
 * Revision 1.11  2009/11/08 15:35:06  wurm
 * fixed bug in cluster splitting that caused empty clusters
 *
 * Revision 1.10  2008-02-02 15:29:30  wurm
 *
 *
 * added private variable
 *
 * Revision 1.9  2007-11-26 16:16:07  ddangelo
 * cleaned up
 *
 * Revision 1.8  2007-11-26 14:06:26  ddangelo
 * completely redone
 *
 * Revision 1.7  2007-11-14 15:46:33  ddangelo
 * new clustering
 * detector segmented in 2/3 parts.
 * individual hits array and cluster definition
 * event writing to be completed
 *
 * Revision 1.6  2007-05-25 15:09:56  ddangelo
 * minor things
 *
 * Revision 1.5  2007-03-28 17:50:16  ddangelo
 * code restyled and cleaned up. no longer helper functions.
 * n_hits_per_channel is now filled. fired channels removed.
 * binned array resized correctly. no extra bins.
 * constants used instead of wired numbers.
 *
 * Revision 1.4  2007-03-22 16:13:02  ddangelo
 * starting to do something...
 *
 * Revision 1.3  2007-02-22 19:59:55  ddangelo
 * first run tests, still to go
 *
 * Revision 1.2  2007/02/21 18:50:09  ddangelo
 * first development. still junk
 *
 * Revision 1.1  2007/02/21 15:48:59  ddangelo
 * added. blank.
 *
 *
 */
