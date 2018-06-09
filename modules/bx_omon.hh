/* BOREXINO On-line MONitor (server)
 *
 * Author: Andrew Sabelnikov <sabelnik@lngs.infn.it>
 * Maintainer: Andrew Sabelnikov <sabelnik@lngs.infn.it>
 *
 * $Id: bx_omon.hh,v 1.10 2008/08/23 11:48:05 sabelnik Exp $
 *
 * 
 */

#ifndef _BX_OMON_HH
#define _BX_OMON_HH

#define MAXCHID 2241
#define MAXFWFD 100
#define MAXOD   256

#include "bx_base_module.hh"
#include "bx_rec_general.hh"
#include <stdlib.h>
#include <stdio.h>

class bx_laben_event;
class TH1F;
class TH2F;

class bx_omon: public bx_base_module {
  public:
  // this section if fixed; do not edit
    bx_omon ();
    virtual ~bx_omon () {}

    virtual void begin ();
    virtual bx_echidna_event* doit (bx_echidna_event *ev);
    virtual void end ();
  private:
  // this section is free; add here members, methods, histo pointers
  //  int i4_times;
  //  std::vector<double> *my_vector;
  //  bool m_check_this_and_that(const bx_laben_event& er);
  //  TH1F* my_histo_check;
  //  TH2F* my_histo_compute;
    int law, after_e, after_t, count;
    time_t last_condition_time;
    std::string path;
    void update_histos(bx_echidna_event *ev); // update histograms with current event
    int condition();      // calculate export condition
    int reset_condition(); //
    void export_ascii();  // if condition then export
    void MyDump1DHisto(TH1F*, FILE* );  //
    void ResetHistos(); // reset all histograms
    TH1F* test;

// Laben and general purpose histograms
    TH1F* TriggerRateHistory; // почему нельзя на одной срочке???
    TH1F* Trigger1RateHistory;
    TH1F* Trigger2RateHistory;
    TH1F* Trigger3RateHistory;
    TH1F* Trigger4RateHistory;
    TH1F* Trigger5RateHistory;
    TH1F* Trigger6RateHistory;
    TH1F* Trigger7RateHistory;
    TH1F* Trigger8RateHistory;
    TH1F* Trigger9RateHistory;
    TH1F* TriggerRateHistogram;
    TH1F* EventSizeHistory;
    TH1F* EventSizeDistribution;
    TH1F* HitsTotalHistory;
    TH1F* HitsTotalDistribution;
    TH1F* HitsWConeHistory;
    TH1F* HitsWConeDistribution;
    TH1F* HitsNConeHistory;
    TH1F* HitsNConeDistribution;
    TH1F* TriggerTypeDistribution;
    TH1F* EventsTimeDistribution;
    TH1F* GPS_ppc0_Distribution;
    TH1F* GPS_build_Distribution;
    TH1F* DistribAverTime;
    TH1F* DistribRefTime;
    TH1F* QuasyEnergy;
    TH1F* QuasyEnergyT01;
    TH1F* QuasyEnergyT08;
    TH1F* QuasyEnergyT32;
    TH1F* QuasyEnergyT64;
    TH1F* QuasyEnergyCutted;
    std::vector<TH1F*> 
      RawChargeLCh, 
      ChargeCutLCh, 
      RawBaseLCh,
      TimeLCh,
      HistoryLCh;
    int HisPmtCounts[MAXCHID];
    float HisPmtOccupancy[MAXCHID];
    int events;
    int events_mu;
    int run_number;
    int trigger_type;
    time_t last_update;

// FWFD histograms
    std::vector<TH1F*> 
      PedCh,
      AmplCh,
      PeakCh,
      ChargeCh;
    TH1F* PeakNSum;
    TH1F* AmplNSum;
    TH1F* ChargeNSum;
    TH1F* NWin;
    TH1F* EvFlag;
    TH1F* TTrig;
    TH1F* NPulse;
    TH1F* TotSize;

// Outer Detector histograms
    std::vector<TH1F*> 
      ODLeadTimeLCh,
      ODTrailTimeLCh,
      ODChargeLCh,
      ODHistoryLCh;
    TH1F* muTriggerRateHistory;
    TH1F* muTriggerRateHistogram;
    TH1F* muHitsTotalHistory;
    TH1F* muHitsTotalDistribution;
    TH1F* ODQuasyEnergy;
    int MuonPmtCounts[MAXOD];
    float MuonPmtOccupancy[MAXOD];

};

#endif
/*
 * $Log: bx_omon.hh,v $
 * Revision 1.10  2008/08/23 11:48:05  sabelnik
 * OD mapping errors fixed
 *
 * Revision 1.9  2008-06-30 09:12:48  sabelnik
 * OD map in OMON
 *
 * Revision 1.8  2008-03-31 15:20:32  sabelnik
 * trigger types history
 *
 * Revision 1.7  2008-01-04 09:26:52  pallas
 * A small change to allow trigger type selection in online monitor histograms
 * Also, update law based on time implemented
 *
 * Revision 1.6  2007-04-02 18:14:38  sabelnik
 * triger types
 *
 * Revision 1.5  2006-09-04 16:04:40  sabelnik
 * export all histos (not filled!)
 *
 * Revision 1.4  2006/08/25 16:45:10  sabelnik
 * complete set of histograms booked
 *
 * Revision 1.3  2006/08/18 15:22:46  sabelnik
 * charge histograms added
 *
 * Revision 1.2  2006/08/08 11:21:56  sabelnik
 * detector and electronics maps implemented.
 * OMON web pages are generated from now on
 * by means of Echidna-based histogram server.
 *
 * Revision 1.1  2006/07/13 16:59:53  sabelnik
 * first entry
 *
 * Revision 1.0  2006/05/23 19:00:18  sabelnik
 * derivated from bx_test_module
 *
 */
