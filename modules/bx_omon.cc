/* BOREXINO On-line MONitor (server side)
 *
 * Author: Andrew Sabelnikov <sabelnik@lngs.infn.it>
 * Maintainer: Andrew Sabelnikov <sabelnik@lngs.infn.it>
 *
 * $Id: bx_omon.cc,v 1.27 2009/09/16 16:17:23 sabelnik Exp $
 *
 * Implementation of bx_omon
 */

#include "bx_omon.hh"
#include "bx_echidna_event.hh"
#include "barn_interface.hh"
#include "db_run.hh"
#include "db_profile.hh"
#include "bx_dbi.hh"
#include "TH1F.h"
#include "TH2F.h"

bx_omon::bx_omon (): bx_base_module("bx_omon", bx_base_module::main_loop) {
}

#if ROOT_VERSION_CODE >= ROOT_VERSION(6,0,0)
#define SET_REBIN(h) h->SetCanExtend(TH1::kAllAxes);
#else
#define SET_REBIN(h) h->SetBit(TH1::kCanRebin);
#endif

//-----------------------------------------------------------------------------
// module interface
void bx_omon::begin () {
  get_message (bx_message::debug) << "OMON begin" << dispatch;
  // Get user parameters
  law = (uint32_t)(get_parameter("law").get_int());
  path = get_parameter("path").get_string();
  after_e = (get_parameter("after_e").get_int());
  after_t = (get_parameter("after_t").get_int());
  if (law == 1) {count = after_e;};
  last_condition_time = 0;
  if (law==1)
  	get_message (bx_message::info) << "OMON condition based on number of events:" << after_e << dispatch;
  else if (law==2)
  	get_message (bx_message::info) << "OMON condition based on number of seconds:" << after_t << dispatch;
  else {
  	get_message (bx_message::warn) << "OMON wrong condition. Setting time to 30 s" << dispatch;
	law = 2;
	after_t = 30;
  }

  //  create all nessesary histograms here (hm...)
  char name1[100];
  sprintf(name1,"omon-test-histo");
  char title[100];
  sprintf(title,"OMON test Distribution");

  test = new TH1F(name1,title,201,-100,100);
  barn_interface::get ()->store (barn_interface::file, test, this);
  test->Fill(1);

  for(int32_t i=0; i<MAXCHID; i++) {
    HisPmtCounts[i] = 0;
    HisPmtOccupancy[i] = 0;

    sprintf(name1,"RawCharge_LCh=%04d", i);
    sprintf(title,"RawCharge_LCh=%04d", i);
    RawChargeLCh.push_back(new TH1F(name1, title, 256, 0., 256.));

    sprintf(name1,"RawBase_LCh=%04d", i);
    sprintf(title,"RawBase_LCh=%04d", i);
    RawBaseLCh.push_back(new TH1F(name1, title, 256, 0., 256.));

    sprintf(name1,"ChargeCut_LCh=%04d", i);
    sprintf(title,"ChargeCut_LCh=%04d", i);
    ChargeCutLCh.push_back(new TH1F(name1, title, 256, 0., 16.));

    sprintf(name1,"Time_LCh=%04d", i);
    sprintf(title,"Time_LCh=%04d", i);
    TimeLCh.push_back(new TH1F(name1, title, 256, 0., 500.));
    SET_REBIN(TimeLCh[i]);

    sprintf(name1,"History_LCh=%04d", i);
    sprintf(title,"History_LCh=%04d", i);
    HistoryLCh.push_back(new TH1F(name1, title, 3, 0., 3.));
    SET_REBIN(HistoryLCh[i]);
  };

  sprintf(name1,"DistribAverTime");
  sprintf(title,"DistribAverTime");
  DistribAverTime = new TH1F(name1, title, 340, -100000., 3400000.);

  sprintf(name1,"DistribRefTime");
  sprintf(title,"DistribRefTime");
  DistribRefTime = new TH1F(name1, title, 340, -6500., 500.);

  sprintf(name1,"TriggerRateHistory");
  sprintf(title,"TriggerRateHistory");
  TriggerRateHistory = new TH1F(name1, title, 3, 0., 3.);
  SET_REBIN(TriggerRateHistory);

// --- trg types beg
  sprintf(name1,"Trigger1ateHistory");
  sprintf(title,"Trigger1ateHistory");
  Trigger1RateHistory = new TH1F(name1, title, 3, 0., 3.);
  SET_REBIN(Trigger1RateHistory);

  sprintf(name1,"Trigger2ateHistory");
  sprintf(title,"Trigger2ateHistory");
  Trigger2RateHistory = new TH1F(name1, title, 3, 0., 3.);
  SET_REBIN(Trigger2RateHistory);

  sprintf(name1,"Trigger3ateHistory");
  sprintf(title,"Trigger3ateHistory");
  Trigger3RateHistory = new TH1F(name1, title, 3, 0., 3.);
  SET_REBIN(Trigger3RateHistory);

  sprintf(name1,"Trigger4ateHistory");
  sprintf(title,"Trigger4ateHistory");
  Trigger4RateHistory = new TH1F(name1, title, 3, 0., 3.);
  SET_REBIN(Trigger4RateHistory);

  sprintf(name1,"Trigger5ateHistory");
  sprintf(title,"Trigger5ateHistory");
  Trigger5RateHistory = new TH1F(name1, title, 3, 0., 3.);
  SET_REBIN(Trigger5RateHistory);

  sprintf(name1,"Trigger6ateHistory");
  sprintf(title,"Trigger6ateHistory");
  Trigger6RateHistory = new TH1F(name1, title, 3, 0., 3.);
  SET_REBIN(Trigger6RateHistory);

  sprintf(name1,"Trigger7ateHistory");
  sprintf(title,"Trigger7ateHistory");
  Trigger7RateHistory = new TH1F(name1, title, 3, 0., 3.);
  SET_REBIN(Trigger7RateHistory);

  sprintf(name1,"Trigger8ateHistory");
  sprintf(title,"Trigger8ateHistory");
  Trigger8RateHistory = new TH1F(name1, title, 3, 0., 3.);
  SET_REBIN(Trigger8RateHistory);

  sprintf(name1,"Trigger9ateHistory");
  sprintf(title,"Trigger9ateHistory");
  Trigger9RateHistory = new TH1F(name1, title, 3, 0., 3.);
  SET_REBIN(Trigger9RateHistory);

// --- trg types end

  sprintf(name1,"EventSizeHistory");
  sprintf(title,"EventSizeHistory");
  EventSizeHistory = new TH1F(name1, title, 3, 0., 3.);
  SET_REBIN(EventSizeHistory);

  sprintf(name1,"EventSizeDistribution");
  sprintf(title,"EventSizeDistribution");
  EventSizeDistribution = new TH1F(name1, title, 100, 0., 3.);
  SET_REBIN(EventSizeDistribution);

  sprintf(name1,"HitsTotalHistory");
  sprintf(title,"HitsTotalHistory");
  HitsTotalHistory = new TH1F(name1, title, 3, 0., 3.);
  SET_REBIN(HitsTotalHistory);

  sprintf(name1,"HitsTotalDistribution");
  sprintf(title,"HitsTotalDistribution");
  HitsTotalDistribution = new TH1F(name1, title, 100, 0., 3.);
  SET_REBIN(HitsTotalDistribution);

  sprintf(name1,"HitsWConeHistory");
  sprintf(title,"HitsWConeHistory");
  HitsWConeHistory = new TH1F(name1, title, 3, 0., 3.);
  SET_REBIN(HitsWConeHistory);

  sprintf(name1,"HitsWConeDistribution");
  sprintf(title,"HitsWConeDistribution");
  HitsWConeDistribution = new TH1F(name1, title, 100, 0., 3.);
  SET_REBIN(HitsWConeDistribution);

  sprintf(name1,"HitsNConeHistory");
  sprintf(title,"HitsNConeHistory");
  HitsNConeHistory = new TH1F(name1, title, 3, 0., 3.);
  SET_REBIN(HitsNConeHistory);

  sprintf(name1,"HitsNConeDistribution");
  sprintf(title,"HitsNConeDistribution");
  HitsNConeDistribution = new TH1F(name1, title, 100, 0., 3.);
  SET_REBIN(HitsNConeDistribution);

  sprintf(name1,"TriggerTypeDistribution");
  sprintf(title,"TriggerTypeDistribution");
  TriggerTypeDistribution = new TH1F(name1, title, 3, 0., 3.);
  SET_REBIN(TriggerTypeDistribution);

  sprintf(name1,"EventsTimeDistribution");
  sprintf(title,"EventsTimeDistribution");
  EventsTimeDistribution = new TH1F(name1, title, 82, 0., 9.);
  EventsTimeDistribution->GetXaxis()->SetBinLabel( 1,"1_ns");
  EventsTimeDistribution->GetXaxis()->SetBinLabel(10,"10_ns");
  EventsTimeDistribution->GetXaxis()->SetBinLabel(19,"100_ns");
  EventsTimeDistribution->GetXaxis()->SetBinLabel(28,"1_mks");
  EventsTimeDistribution->GetXaxis()->SetBinLabel(37,"10_mks");
  EventsTimeDistribution->GetXaxis()->SetBinLabel(46,"100_mks");
  EventsTimeDistribution->GetXaxis()->SetBinLabel(55,"1_ms");
  EventsTimeDistribution->GetXaxis()->SetBinLabel(64,"10_ms");
  EventsTimeDistribution->GetXaxis()->SetBinLabel(73,"100_ms");
  EventsTimeDistribution->GetXaxis()->SetBinLabel(82,"1_s");
  for(int32_t i=1; i<=82; i++) {
    EventsTimeDistribution->GetXaxis()->SetBinLabel(i,"-");
  }

  sprintf(name1,"GPS_ppc0_Distribution");
  sprintf(title,"GPS_ppc0_Distribution");
  GPS_ppc0_Distribution = new TH1F(name1, title, 201, -100, 100);

  sprintf(name1,"GPS_build_Distribution");
  sprintf(title,"GPS_build_Distribution");
  GPS_build_Distribution = new TH1F(name1, title, 201, -100, 100);

  sprintf(name1,"QuasyEnergy");
  sprintf(title,"QuasyEnergy");
  QuasyEnergy = new TH1F(name1, title, 700, 0., 3.);
  SET_REBIN(QuasyEnergy);

  sprintf(name1,"QuasyEnergyT01");
  sprintf(title,"QuasyEnergyT01");
  QuasyEnergyT01 = new TH1F(name1, title, 700, 0., 3.);
  SET_REBIN(QuasyEnergyT01);

  sprintf(name1,"QuasyEnergyT08");
  sprintf(title,"QuasyEnergyT08");
  QuasyEnergyT08 = new TH1F(name1, title, 700, 0., 3.);
  SET_REBIN(QuasyEnergyT08);

  sprintf(name1,"QuasyEnergyT32");
  sprintf(title,"QuasyEnergyT32");
  QuasyEnergyT32 = new TH1F(name1, title, 700, 0., 3.);
  SET_REBIN(QuasyEnergyT32);

  sprintf(name1,"QuasyEnergyT64");
  sprintf(title,"QuasyEnergyT64");
  QuasyEnergyT64 = new TH1F(name1, title, 700, 0., 3.);
  SET_REBIN(QuasyEnergyT64);

  sprintf(name1,"QuasyEnergyCutted");
  sprintf(title,"QuasyEnergyCutted");
  QuasyEnergyCutted = new TH1F(name1, title, 700, 0., 3.);
  SET_REBIN(QuasyEnergyCutted);

  sprintf(name1,"TriggerRateHistogram");
  sprintf(title,"TriggerRateHistogram");
  TriggerRateHistogram = new TH1F(name1, title, 3, 0., 3.); //!!!
  SET_REBIN(TriggerRateHistogram);

  for(int32_t i=0; i<MAXOD; i++) {

    sprintf(name1,"ODLeadTimeLCh=%03d", i);
    sprintf(title,"ODLeadTimeLCh=%03d", i);
    ODLeadTimeLCh.push_back(new TH1F(name1, title, 256, 0., 100.));
    SET_REBIN(ODLeadTimeLCh[i]);

    sprintf(name1,"ODTrailTimeLCh=%03d", i);
    sprintf(title,"ODTrailTimeLCh=%03d", i);
    ODTrailTimeLCh.push_back(new TH1F(name1, title, 256, 0., 100.));
    SET_REBIN(ODTrailTimeLCh[i]);

    sprintf(name1,"ODChargeLCh=%03d", i);
    sprintf(title,"ODChargeLCh=%03d", i);
    ODChargeLCh.push_back(new TH1F(name1, title, 256, 0., 100.));
    SET_REBIN(ODChargeLCh[i]);

    sprintf(name1,"ODHistoryLCh=%03d", i);
    sprintf(title,"ODHistoryLCh=%03d", i);
    ODHistoryLCh.push_back(new TH1F(name1, title, 3, 0., 3.));
    SET_REBIN(ODHistoryLCh[i]);

    MuonPmtOccupancy[i] = 0;
    MuonPmtCounts[i] = 0;
  }

  sprintf(name1,"ODQuasyEnergy");
  sprintf(title,"ODQuasyEnergy");
  ODQuasyEnergy = new TH1F(name1, title, 500, 0., 3.);
  SET_REBIN(ODQuasyEnergy);

  sprintf(name1,"muTriggerRateHistory");
  sprintf(title,"muTriggerRateHistory");
  muTriggerRateHistory = new TH1F(name1, title, 3, 0., 3.);
  SET_REBIN(muTriggerRateHistory);

  sprintf(name1,"muTriggerRateHistogram"); // !!!
  sprintf(title,"muTriggerRateHistogram");
  muTriggerRateHistogram = new TH1F(name1, title, 3, 0., 3.);
  SET_REBIN(muTriggerRateHistogram);

  sprintf(name1,"muHitsTotalHistory");
  sprintf(title,"muHitsTotalHistory");
  muHitsTotalHistory = new TH1F(name1, title, 3, 0., 3.);
  SET_REBIN(muHitsTotalHistory);

  sprintf(name1,"muHitsTotalDistribution");
  sprintf(title,"muHitsTotalDistribution");
  muHitsTotalDistribution = new TH1F(name1, title, 100, 0., 3.);
  SET_REBIN(muHitsTotalDistribution);

  sprintf(name1,"EvFlag");
  sprintf(title,"EvFlag");
  EvFlag = new TH1F(name1, title, 128, 0., 128.);

  sprintf(name1,"TTrig"); // trigger time in ticks
  sprintf(title,"TTrig");
  TTrig = new TH1F(name1, title, 200, -200., 0.);

  for(int32_t i=0; i<MAXFWFD; i++) {

    sprintf(name1,"PedCh=%03d", i); // mean pedestal in channel
    sprintf(title,"PedCh=%03d", i);
    PedCh.push_back(new TH1F(name1, title, 256, 0., 256.));

    sprintf(name1,"PeakCh=%03d", i); // peak position in channel
    sprintf(title,"PeakCh=%03d", i);
    PeakCh.push_back(new TH1F(name1, title, 100, 0., 100.));
    SET_REBIN(PeakCh[i]);

    sprintf(name1,"AmplCh=%03d", i); // peak position in channel
    sprintf(title,"AmplCh=%03d", i);
    AmplCh.push_back(new TH1F(name1, title, 256, 0., 256.));
    SET_REBIN(AmplCh[i]);

    sprintf(name1,"ChargeCh=%03d", i); // peak position in channel
    sprintf(title,"ChargeCh=%03d", i);
    ChargeCh.push_back(new TH1F(name1, title, 256, 0., 256.));
    SET_REBIN(ChargeCh[i]);
  }

  sprintf(name1,"PeakNSum"); // Total Sum Peak Position
  sprintf(title,"PeakNSum");
  PeakNSum = new TH1F(name1, title, 100, 0., 100.);
  SET_REBIN(PeakNSum);

  sprintf(name1,"AmplNSum"); // Total Sum Amplitude
  sprintf(title,"AmplNSum");
  AmplNSum = new TH1F(name1, title, 256, 0., 256.);
  SET_REBIN(AmplNSum);

  sprintf(name1,"ChargeNSum"); // Total Sum Charge
  sprintf(title,"ChargeNSum");
  ChargeNSum = new TH1F(name1, title, 256, 0., 256.);
  SET_REBIN(ChargeNSum);

  sprintf(name1,"NWin"); // Total Sum Charge
  sprintf(title,"NWin");
  NWin = new TH1F(name1, title, 100, 0., 100.);
  SET_REBIN(NWin);

  sprintf(name1,"NPulse"); // Total Sum Charge
  sprintf(title,"NPulse");
  NPulse = new TH1F(name1, title, 256, 0., 256.);
  SET_REBIN(NPulse);

  sprintf(name1,"TotSize"); // Total Sum Charge
  sprintf(title,"TotSize");
  TotSize = new TH1F(name1, title, 256, 0., 256.);
  SET_REBIN(TotSize);

  events = 0;
  events_mu = 1;
//  for(int32_t i=0; i<MAXCHID; i++) HisPmtCounts[i]=0; // почему-то эта строчка не работала? фиг скобки?
  run_number = 0;
  trigger_type = 0xFFFF;
  last_update = 0;
}


//-----------------------------------------------------------------------------
bx_echidna_event* bx_omon::doit (bx_echidna_event *ev) {

  get_message (bx_message::debug) << "OMON event " << dispatch;

  // get time and read trigger type file every 10 s
  // in case of change, reset histograms
  // if file not found, ignore and use all trigger types
  time_t now = time(0);
  uint32_t diff = now - last_update;
  if ( diff > 10 ) {
  	FILE *fp=fopen("/bxwww/data/omon_trigger_type.txt","rt");
	if (!fp) {
		get_message (bx_message::warn) << "Trigger type file not found. All trigger type used" << dispatch << std::endl;
		trigger_type = 0xFFFF;
	} else {
		char dummy[200];
		int32_t new_type = 0;
		fscanf(fp,"%s\n",dummy);
		fclose(fp);

		if ( strcmp(dummy,"all") == 0 ) {
			new_type = 0xFFFF;
		}
		else if ( strcmp(dummy,"neutrino") == 0 ) {
			new_type = 1;
		}
		else if ( strcmp(dummy,"random") == 0 ) {
			new_type = 64;
		}
		else if ( strcmp(dummy,"pulser") == 0 ) {
			new_type = 32;
		}
		else if ( strcmp(dummy,"laser") == 0 ) {
			new_type = 8;
		}
		else if ( strcmp(dummy,"neutron") == 0 ) {
			new_type = 128;
		}
		else if ( strcmp(dummy,"muon") == 0 ) {
			new_type = 2;
		}
		if ( new_type != trigger_type ) {
			get_message (bx_message::info) << "Trigger type changed. Old = " << trigger_type 
			                               << "   New = " << new_type << dispatch;
			ResetHistos(); 
		}

		trigger_type = new_type;
	}
	last_update = now;
  }

  // update histograms with current event
  // filter on trigger type if required
  if ( trigger_type & ev->get_trigger().get_trgtype() ) {
  	update_histos(ev);
  	// calculate export condition
  	if (condition()) {
  		// if condition then export
    		reset_condition();
    		export_ascii();
  	}
  }

  return ev;
}

//-----------------------------------------------------------------------------
void bx_omon::end () {
  get_message (bx_message::debug) << "OMON end" << dispatch;

  for(int32_t i=0; i<MAXCHID; i++) {
    delete RawChargeLCh[i];  
  }
}

//-----------------------------------------------------------------------------
// Private method(s)

//-----------------------------------------------------------------------------
void bx_omon::update_histos(bx_echidna_event *ev) // update histograms with current event
{
  get_message (bx_message::debug) << "update_histos. events:" << events << " mu " << events_mu << dispatch;
  // list to fill
  events ++;
  // cycle on hits in event
  const bx_laben_event& er = ev->get_laben ();

// ----- event based
// trigger rate history
  int32_t hour = ev->get_trigger ().get_hour();
  int32_t min = ev->get_trigger ().get_min();
  char tmstmp[20];
  sprintf(tmstmp, "%02d:%02d", hour, min);
  TriggerRateHistory->Fill(tmstmp, 1./60.);
  if (ev->get_trigger().get_trgtype() ==   1) { Trigger1RateHistory->Fill(tmstmp, 1./60.);};
  if (ev->get_trigger().get_trgtype() ==   2) { Trigger2RateHistory->Fill(tmstmp, 1./60.);};
//  if (ev->get_trigger().get_trgtype() ==   4) { Trigger3RateHistory->Fill(tmstmp, 1./60.);};
  if (ev->get_trigger().get_trgtype() ==   8) { Trigger3RateHistory->Fill(tmstmp, 1./60.);};
//  if (ev->get_trigger().get_trgtype() ==  16) { Trigger5RateHistory->Fill(tmstmp, 1./60.);};
  if (ev->get_trigger().get_trgtype() ==  32) { Trigger4RateHistory->Fill(tmstmp, 1./60.);};
  if (ev->get_trigger().get_trgtype() ==  64) { Trigger5RateHistory->Fill(tmstmp, 1./60.);};
  if (ev->get_trigger().get_trgtype() == 128) { Trigger6RateHistory->Fill(tmstmp, 1./60.);};
//  Trigger9RateHistory->Fill(tmstmp, 1./60.); // dirty

  HitsTotalHistory->Fill(tmstmp, er.get_raw_nhits());

  char TrigType[20];
  sprintf(TrigType,"%02d", ev->get_trigger().get_trgtype());
//                           ev->get_trigger ().get_trgtype ();
  TriggerTypeDistribution->Fill(TrigType, 1);

  QuasyEnergy->Fill(er.get_decoded_nhits(),1);

//FIXME!!!  QuasyEnergyT01->SetBins(QuasyEnergy->GetNbinsX(), QuasyEnergy->GetNbinsX());

  if (ev->get_trigger().get_trgtype() == 1)
  {
    QuasyEnergyT01->Fill(er.get_decoded_nhits(),1);
  }

  if (ev->get_trigger().get_trgtype() == 8)
  {
    QuasyEnergyT08->Fill(er.get_decoded_nhits(),1);
  }

  if (ev->get_trigger().get_trgtype() == 32)
  {
    QuasyEnergyT32->Fill(er.get_decoded_nhits(),1);
  }

  if (ev->get_trigger().get_trgtype() == 64)
  {
    QuasyEnergyT64->Fill(er.get_decoded_nhits(),1);
  }

  if (er.get_decoded_nhits() < 500) {
    QuasyEnergyCutted->Fill(er.get_decoded_nhits(),1);
  }

// ----- hit-based histograms

  for (int32_t i=0; i<er.get_decoded_nhits(); i++)
  {
    const bx_laben_decoded_hit& h = er.get_decoded_hit(i);
    int32_t lg = h.get_raw_hit().get_logical_channel();
    HisPmtCounts[lg]++;
    HisPmtOccupancy[lg] = (float)HisPmtCounts[lg]/(float)events;

    RawChargeLCh[lg]->Fill(h.get_raw_hit().get_peak() - h.get_raw_hit().get_base() );
    RawBaseLCh[lg]->Fill(h.get_raw_hit().get_base());
    TimeLCh[lg]->Fill(h.get_raw_time());
    HistoryLCh[lg]->Fill(tmstmp, 1./60.);
  };

  run_number = ev->get_run_number(); // not much elegant, sorry

// ----- muon histograms
  if (ev->get_trigger().get_trgtype() ==   2) { 
          events_mu ++;
	  const bx_muon_event& em = ev->get_muon ();
	//  get_message (bx_message::warn) << "total muon hits: " << em.get_decoded_nhits() << "; logical channels: " << dispatch;
	// ----- cycle on muon hits
	  for (int32_t i=0; i<em.get_decoded_nhits(); i++) // raw or decoded?
	  {
	    const bx_muon_decoded_hit &decoded_hit = em.get_decoded_hit (i);
	    int32_t ml = decoded_hit.get_raw_hit().get_logical_channel();
	//    get_message (bx_message::warn) << ml << ", " << dispatch;
	    MuonPmtCounts[ml - 3000]++;
	    MuonPmtOccupancy[ml - 3000] = (float)MuonPmtCounts[ml - 3000]/(float)events_mu;
        ODChargeLCh[ml - 3000]->Fill(decoded_hit.get_raw_hit().get_time_diff());
        ODLeadTimeLCh[ml - 3000]->Fill(decoded_hit.get_raw_hit().get_lead_time());
        ODTrailTimeLCh[ml - 3000]->Fill(decoded_hit.get_raw_hit().get_trail_time());
        ODHistoryLCh[ml - 3000]->Fill(tmstmp, 1./60.);
	  }
	  muHitsTotalHistory->Fill(tmstmp, em.get_raw_nhits());

  }
//  get_message (bx_message::warn)  << dispatch;

  test->Fill(2);
  return;
}

//-----------------------------------------------------------------------------
int32_t bx_omon::condition() // update histograms with current event
{
  get_message (bx_message::debug) << "condition" << count << dispatch;
  // pick up rule
  if (law == 1) {
  	count--;
  	return !count; //1=true; 0=false
  }
  if (law==2) {
  	time_t now = time(0);
	int32_t diff = now - last_condition_time;
	if ( diff > after_t )
		return 1;
  }
  return 0;
}

//-----------------------------------------------------------------------------
int32_t bx_omon::reset_condition() // update histograms with current event
{
  get_message (bx_message::debug) << "reset condition" << dispatch;
  // pick up rule
  // reset param accord to rule
  if (law == 1) {count = after_e;};
  if (law == 2) {last_condition_time = time(0);}
  return 1; // 0=false
}

//-----------------------------------------------------------------------------
void bx_omon::export_ascii() // update histograms with current event
/*!!!
  ofstream fpccd;
  fpccd.open(fname.c_str());
  fpccd.close();
*/
{
  get_message (bx_message::debug) << "export_ascii" << dispatch;
  std::string fname = "";
  // === list to export ===
  // == laben part ==

  // charge histograms and index // files ccd.txt and cci.txt
  fname = path + "/ccd.txt";
//  get_message (bx_message::debug) << fname << dispatch;
  FILE* fpccd = fopen(fname.c_str(),"wt+");
  fname = path + "/cci.txt";
  FILE* fpcci = fopen(fname.c_str(),"wt+");
  fprintf(fpccd,"# nchannels\n2240\n");
  for (int32_t ch=1; ch<MAXCHID; ch++) {
    fprintf(fpccd,"# Logical Channel %04d\n",ch);
    fprintf(fpcci,"%04ld\n",ftell(fpccd) );
    MyDump1DHisto(RawChargeLCh[ch], fpccd);
  }
  fclose(fpccd);
  fclose(fpcci);

  // charge histograms and index for high Voltage monitor // files vcd.txt and vci.txt
  fname = path + "/vcd.txt";
  FILE* fpvcd = fopen(fname.c_str(),"wt+");
  fname = path + "/vci.txt";
  FILE* fpvci = fopen(fname.c_str(),"wt+");
  fprintf(fpvcd,"# nchannels\n2240\n");
  for (int32_t ch=1; ch<MAXCHID; ch++) {
    fprintf(fpvcd,"# Logical Channel %04d\n",ch);
    fprintf(fpvci,"%04ld\n",ftell(fpvcd) );
    MyDump1DHisto(ChargeCutLCh[ch], fpvcd);
  }
  fclose(fpvcd);
  fclose(fpvci);

  // base sampling histos and index  // files bcd.txt and bci.txt
  fname = path + "/bcd.txt";
  FILE* fpbcd = fopen(fname.c_str(),"wt+");
  fname = path + "/bci.txt";
  FILE* fpbci = fopen(fname.c_str(),"wt+");
  fprintf(fpbcd,"# nchannels\n2240\n");
  for (int32_t ch=1; ch<MAXCHID; ch++) {
    fprintf(fpbcd,"# Logical Channel %04d\n",ch);
    fprintf(fpbci,"%04ld\n",ftell(fpbcd) );
    MyDump1DHisto(RawBaseLCh[ch], fpbcd);
  }
  fclose(fpbcd);
  fclose(fpbci);

  // time histos and index (time relative to trigger) // files tcd.txt and tci.txt
  fname = path + "/tcd.txt";
  FILE* fptcd = fopen(fname.c_str(),"wt+");
  fname = path + "/tci.txt";
  FILE* fptci = fopen(fname.c_str(),"wt+");
  fprintf(fptcd,"# nchannels\n2240\n");
  for (int32_t ch=1; ch<MAXCHID; ch++) {
    fprintf(fptcd,"# Logical Channel %04d\n",ch);
    fprintf(fptci,"%04ld\n",ftell(fptcd) );
    MyDump1DHisto(TimeLCh[ch], fptcd);
  }
  fclose(fptcd);
  fclose(fptci);

  // countrate history -- histos and index // files hcd.txt and hci.txt
  fname = path + "/hcd.txt";
  FILE* fphcd = fopen(fname.c_str(),"wt+");
  fname = path + "/hci.txt";
  FILE* fphci = fopen(fname.c_str(),"wt+");
  fprintf(fphcd,"# nchannels\n2240\n");
  for (int32_t ch=1; ch<MAXCHID; ch++) {
    fprintf(fphcd,"# Logical Channel %04d\n",ch);
    fprintf(fphci,"%04ld\n",ftell(fphcd) );
    for(int32_t i=0; i < HistoryLCh[ch]->GetNbinsX(); i++) {
      fprintf(fphcd, "%s %06.2f\n", HistoryLCh[ch]->GetXaxis()->GetBinLabel(i), (float)HistoryLCh[ch]->GetBinContent(i));
    }
  }
  fclose(fphcd);
  fclose(fphci);

  // occupancy
  fname = path + "/hpt.txt";
  FILE* fpoccupancy = fopen(fname.c_str(),"wt+");
  fprintf(fpoccupancy,"## pmt occupancy\n");
  fprintf(fpoccupancy,"2240\n");
  for(int32_t i=1; i<MAXCHID; i++) fprintf(fpoccupancy, "%d %6.2f\n", i, HisPmtOccupancy[i]);
  fclose(fpoccupancy);

  // main trigger history // file mth.txt
  fname = path + "/mth.txt";
  FILE* fptriggerhistory = fopen(fname.c_str(),"wt+");
  for(int32_t i=0; i < TriggerRateHistory->GetNbinsX(); i++) {
//    if (TriggerRateHistory->GetXaxis()->GetBinLabel(i) != "") { // this line does not work... ?!
    if (strlen(TriggerRateHistory->GetXaxis()->GetBinLabel(i)) != 0) { // this line does not work... ?!
      fprintf(fptriggerhistory, "%s %06.2f\n", TriggerRateHistory->GetXaxis()->GetBinLabel(i), (float)TriggerRateHistory->GetBinContent(i));
    }
  }
  fclose(fptriggerhistory);

  // triggers history // file mthm.txt
  fname = path + "/mthm.txt";
  FILE* mthm = fopen(fname.c_str(),"wt+");
  for(int32_t i=0; i < TriggerRateHistory->GetNbinsX(); i++) {
//    if ((TriggerRateHistory->GetXaxis()->GetBinLabel(i) != "") 
    if ((strlen(TriggerRateHistory->GetXaxis()->GetBinLabel(i)) != 0) 
		&& (TriggerRateHistory->GetBinContent(i) != 0)) { // this line does not work... ?!
      fprintf(mthm,
              "%s %06.2f %06.2f %06.2f %06.2f %06.2f %06.2f \n",
              TriggerRateHistory->GetXaxis()->GetBinLabel(i),
              (float)Trigger1RateHistory->GetBinContent(i),
              (float)Trigger2RateHistory->GetBinContent(i),
              (float)Trigger3RateHistory->GetBinContent(i),
              (float)Trigger4RateHistory->GetBinContent(i),
              (float)Trigger5RateHistory->GetBinContent(i),
              (float)Trigger6RateHistory->GetBinContent(i) //,
//              (float)Trigger7RateHistory->GetBinContent(i),
//              (float)Trigger8RateHistory->GetBinContent(i),
//              (float)Trigger9RateHistory->GetBinContent(i)
             );
    }
  }
  fclose(mthm);

  // main trigger distribution // file mtd.txt
//  TriggerRateHistogram->Clear();
  TriggerRateHistogram->Reset();
  for(int32_t i=0; i < TriggerRateHistory->GetNbinsX(); i++)
    if (TriggerRateHistory->GetBinContent(i) != 0)
    {
      TriggerRateHistogram->Fill((float)TriggerRateHistory->GetBinContent(i));
    }

  int32_t n = 22;
  fname = path + "/mtd.txt";
  FILE* fptriggerhistogram = fopen(fname.c_str(),"wt+");
  fprintf(fptriggerhistogram,"%d\n",n);
  for(int32_t i=0; i < n; i++) {
    fprintf(fptriggerhistogram, "%06.2f %06.2f\n", (float)TriggerRateHistogram->GetBinCenter(i), (float)TriggerRateHistogram->GetBinContent(i));
  }
  fprintf(fptriggerhistogram,"%d\n%5.1f\n%5.1f\n", (int32_t)TriggerRateHistogram->GetEntries(), (float)TriggerRateHistogram->GetMean(), (float)TriggerRateHistogram->GetRMS() );
  fclose(fptriggerhistogram);

  // event size history // file esh.txt
  fname = path + "/esh.txt";
  FILE* fpes = fopen(fname.c_str(),"wt+");
  for(int32_t i=0; i < EventSizeHistory->GetNbinsX(); i++)
//    if (EventSizeHistory->GetXaxis()->GetBinLabel(i) != "") // this line does not work... ?!
    if (strlen(EventSizeHistory->GetXaxis()->GetBinLabel(i)) != 0) // this line does not work... ?!
    {
      fprintf(fpes,"%s %06.2f\n", TriggerRateHistory->GetXaxis()->GetBinLabel(i), (float)EventSizeHistory->GetBinContent(i) / (float)TriggerRateHistory->GetBinContent(i) / 60.);
    }
  fclose(fpes);

  // event size distribution // file esd.txt
  fname = path + "/esd.txt";
  FILE* fpesd = fopen(fname.c_str(),"wt+");
  fprintf(fpesd,"%d\n",EventSizeDistribution->GetNbinsX());
  for(int32_t i=0; i < EventSizeDistribution->GetNbinsX(); i++)
    fprintf(fpesd,"%06.2f %06.2f\n",(float)EventSizeDistribution->GetBinCenter(i), (float)EventSizeDistribution->GetBinContent(i));
  fprintf(fpesd,"%d\n%5.1f\n%5.1f\n",
    (int32_t)EventSizeDistribution->GetEntries(), (float)EventSizeDistribution->GetMean(), (float)EventSizeDistribution->GetRMS() );
  fclose(fpesd);

  // total hits -- history  // file nhh.txt
  fname = path + "/nhh.txt";
  FILE* fpnh = fopen(fname.c_str(),"wt+");
  for(int32_t i=0; i < HitsTotalHistory->GetNbinsX(); i++)
//    if (HitsTotalHistory->GetXaxis()->GetBinLabel(i) != "") {
    if (strlen(HitsTotalHistory->GetXaxis()->GetBinLabel(i)) != 0) {
//      fprintf(fpnh,"h %s %06.2f\n", TriggerRateHistory->GetXaxis()->GetBinLabel(i), (float)HitsTotalHistory->GetBinContent(i));
//      fprintf(fpnh,"t %s %06.2f\n", TriggerRateHistory->GetXaxis()->GetBinLabel(i), (float)TriggerRateHistory->GetBinContent(i)*60);
      fprintf(fpnh,"%s %06.2f\n", TriggerRateHistory->GetXaxis()->GetBinLabel(i), (float)HitsTotalHistory->GetBinContent(i)/((float)TriggerRateHistory->GetBinContent(i)*60));
//      fprintf(fpnh,"\n");
//      fprintf(fpnh,"3 %s %06.2f\n", TriggerRateHistory->GetXaxis()->GetBinLabel(i), (float)HitsTotalHistory->GetBinContent(i) / (float)TriggerRateHistory->GetBinContent(i));
//      fprintf(fpnh,"4 %s %06.2f\n", TriggerRateHistory->GetXaxis()->GetBinLabel(i), ((float)HitsTotalHistory->GetBinContent(i) / (float)TriggerRateHistory->GetBinContent(i)) / 3600.);
    }
  fclose(fpnh);

  // total hits -- distribution // file nhd.txt
  fname = path + "/nhd.txt";
  FILE* fpnhd = fopen(fname.c_str(),"wt+");
  fprintf(fpnhd,"%d\n",HitsTotalDistribution->GetNbinsX());
  for(int32_t i=0; i < HitsTotalDistribution->GetNbinsX(); i++)
    fprintf(fpnhd,"%06.2f %06.2f\n", (float)HitsTotalDistribution->GetBinCenter(i), (float)HitsTotalDistribution->GetBinContent(i));
  fprintf(fpnhd,"%d\n%5.1f\n%5.1f\n",
    (int32_t)HitsTotalDistribution->GetEntries(), (float)HitsTotalDistribution->GetMean(), (float)HitsTotalDistribution->GetRMS() );
  fclose(fpnhd);

  // hits with cones -- history // file wch.txt
  fname = path + "/wch.txt";
  FILE* fwch = fopen(fname.c_str(),"wt+");
  for(int32_t i=0; i < HitsWConeHistory->GetNbinsX(); i++)
//    if (HitsWConeHistory->GetXaxis()->GetBinLabel(i) != "") {
    if (strlen(HitsWConeHistory->GetXaxis()->GetBinLabel(i)) != 0) {
      fprintf(fwch,"%s %06.2f\n", TriggerRateHistory->GetXaxis()->GetBinLabel(i), (float)HitsWConeHistory->GetBinContent(i) / (float)TriggerRateHistory->GetBinContent(i) / 60.);
    }
  fclose(fwch);

  // hits with cones -- distribution // file wcd.txt
  fname = path + "/wcd.txt";
  FILE* fpwcd = fopen(fname.c_str(),"wt+");
  fprintf(fpwcd,"%d\n",HitsWConeDistribution->GetNbinsX());
  for(int32_t i=0; i < HitsWConeDistribution->GetNbinsX(); i++)
    fprintf(fpwcd,"%06.2f %06.2f\n",(float)HitsWConeDistribution->GetBinCenter(i), (float)HitsWConeDistribution->GetBinContent(i));
  fprintf(fpwcd,"%d\n%5.1f\n%5.1f\n",
    (int32_t)HitsWConeDistribution->GetEntries(), (float)HitsWConeDistribution->GetMean(), (float)HitsWConeDistribution->GetRMS() );
  fclose(fpwcd);

  // hits no cones -- history // file nch.txt
  fname = path + "/nch.txt";
  FILE* fnch = fopen(fname.c_str(),"wt+");
  for(int32_t i=0; i < HitsNConeHistory->GetNbinsX(); i++)
//    if (HitsNConeHistory->GetXaxis()->GetBinLabel(i) != "") {
    if (strlen(HitsNConeHistory->GetXaxis()->GetBinLabel(i)) != 0) {
      fprintf(fnch,"%s %06.2f\n", TriggerRateHistory->GetXaxis()->GetBinLabel(i), (float)HitsNConeHistory->GetBinContent(i) / (float)TriggerRateHistory->GetBinContent(i) / 60.);
    }
  fclose(fnch);

  // hits no cones -- distribution // file ncd.txt
  fname = path + "/ncd.txt";
  FILE* fpncd = fopen(fname.c_str(),"wt+");
  fprintf(fpncd,"%d\n", HitsNConeDistribution->GetNbinsX());
  for(int32_t i=0; i < HitsNConeDistribution->GetNbinsX(); i++)
    fprintf(fpncd,"%06.2f %06.2f\n", (float)HitsNConeDistribution->GetBinCenter(i), (float)HitsNConeDistribution->GetBinContent(i));
  fprintf(fpncd,"%d\n%5.1f\n%5.1f\n",
    (int32_t)HitsNConeDistribution->GetEntries(), (float)HitsNConeDistribution->GetMean(), (float)HitsNConeDistribution->GetRMS() );
  fclose(fpncd);

  // trigger type distribution // file ttd.txt
  fname = path + "/ttd.txt";
  FILE* fpttd = fopen(fname.c_str(),"wt+");
  fprintf(fpttd,"%d\n", TriggerTypeDistribution->GetNbinsX());
  for(int32_t i=0; i < TriggerTypeDistribution->GetNbinsX(); i++)
    fprintf(fpttd, "%s %06.2f\n", TriggerTypeDistribution->GetXaxis()->GetBinLabel(i), (float)TriggerTypeDistribution->GetBinContent(i));
  fprintf(fpttd,"%d\n%5.1f\n%5.1f\n",
    (int32_t)TriggerTypeDistribution->GetEntries(), (float)TriggerTypeDistribution->GetMean(), (float)TriggerTypeDistribution->GetRMS() );
  fclose(fpttd);

  // events time distribution // file tbd.txt
  fname = path + "/tbd.txt";
  FILE* fptbd = fopen(fname.c_str(),"wt+");
  fprintf(fptbd,"%d\n", EventsTimeDistribution->GetNbinsX());
  for(int32_t i=1; i <= EventsTimeDistribution->GetNbinsX(); i++)
    fprintf(fptbd,"%s %06.2f\n", EventsTimeDistribution->GetXaxis()->GetBinLabel(i), (float)EventsTimeDistribution->GetBinContent(i));
  fprintf(fptbd,"%d\n%5.1f\n%5.1f\n",
    (int32_t)EventsTimeDistribution->GetEntries(), (float)EventsTimeDistribution->GetMean(), (float)EventsTimeDistribution->GetRMS() );
  fclose(fptbd);

  // gps time distribution // file gps.txt
  fname = path + "/gps.txt";
  FILE* fpgps = fopen(fname.c_str(),"wt+");
  fprintf(fpgps,"%d\n", GPS_ppc0_Distribution->GetNbinsX());
  for(int32_t i=1; i <= GPS_ppc0_Distribution->GetNbinsX(); i++)
    fprintf(fpgps, "%06.2f %06.2f %06.2f\n", (float)GPS_ppc0_Distribution->GetBinCenter(i), (float)GPS_ppc0_Distribution->GetBinContent(i), (float)GPS_build_Distribution->GetBinContent(i));
  fprintf(fpgps,"%d\n%5.1f\n%5.1f\n",
    (int32_t)GPS_ppc0_Distribution->GetEntries(), (float)GPS_ppc0_Distribution->GetMean(), (float)GPS_ppc0_Distribution->GetRMS() );
  fclose(fpgps);

  // average time distribution // file atd.txt
  fname = path + "/atd.txt";
  FILE* fpatd = fopen(fname.c_str(),"wt+");
  fprintf(fpatd,"%d\n", DistribAverTime->GetNbinsX());
  for(int32_t i=1; i <= DistribAverTime->GetNbinsX(); i++)
    fprintf(fpatd, "%06.2f %06.2f\n", (float)DistribAverTime->GetBinCenter(i), (float)DistribAverTime->GetBinContent(i));
  fprintf(fpatd,"%d\n%5.1f\n%5.1f\n",
    (int32_t)DistribAverTime->GetEntries(), (float)DistribAverTime->GetMean(), (float)DistribAverTime->GetRMS() );
  fclose(fpatd);

  // reference time distribution // file rtd.txt
  fname = path + "/rtd.txt";
  FILE* fprtd = fopen(fname.c_str(), "wt+");
  fprintf(fprtd,"%d\n", DistribRefTime->GetNbinsX());
  for(int32_t i=1; i <= DistribRefTime->GetNbinsX(); i++)
    fprintf(fprtd, "%06.2f %06.2f\n", (float)DistribRefTime->GetBinCenter(i), (float)DistribRefTime->GetBinContent(i));
  fprintf(fprtd,"%d\n%5.1f\n%5.1f\n",
    (int32_t)DistribRefTime->GetEntries(), (float)DistribRefTime->GetMean(), (float)DistribRefTime->GetRMS() );
  fclose(fprtd);

  // file end.txt (energy --  sum of amplitudes -- distribution)
  fname = path + "/end.txt";
  FILE* fpmsad = fopen(fname.c_str(), "wt+");
  fprintf (fpmsad, "#sum of amplitudes -- distribution\n");
  fprintf(fpmsad,"%d\n", QuasyEnergy->GetNbinsX());
  for(int32_t i=1; i <= QuasyEnergy->GetNbinsX(); i++) {
    fprintf(fpmsad,"%06.2f %06.2f %06.2f %06.2f %06.2f %06.2f\n", (float)QuasyEnergy->GetBinCenter(i), 
                                             (float)QuasyEnergy->GetBinContent(i), 
                                             (float)QuasyEnergyT01->GetBinContent(i),
                                             (float)QuasyEnergyT08->GetBinContent(i),
                                             (float)QuasyEnergyT32->GetBinContent(i),
                                             (float)QuasyEnergyT64->GetBinContent(i));
  }
  fprintf(fpmsad,"%d\n%5.1f\n%5.1f\n",
    (int32_t)QuasyEnergy->GetEntries(), (float)QuasyEnergy->GetMean(), (float)QuasyEnergy->GetRMS() );
  fclose(fpmsad);

  // file endc.txt (sum of amplitudes (cutted) -- distribution)
  fname = path + "/endc.txt";
  FILE* fpmsac = fopen(fname.c_str(), "wt+");
  fprintf (fpmsac, "#sum of amplitudes (cutted) -- distribution\n");
  fprintf(fpmsac,"%d\n",QuasyEnergyCutted->GetNbinsX());
  for(int32_t i=1; i <= QuasyEnergyCutted->GetNbinsX(); i++) {
    fprintf(fpmsac,"%06.2f %06.2f\n", (float)QuasyEnergyCutted->GetBinCenter(i), (float)QuasyEnergyCutted->GetBinContent(i));
  }
  fprintf(fpmsac,"%d\n%5.1f\n%5.1f\n",
    (int32_t)QuasyEnergyCutted->GetEntries(), (float)QuasyEnergyCutted->GetMean(), (float)QuasyEnergyCutted->GetRMS() );
  fclose(fpmsac);

  // run number // file ofrn.txt
  fname = path + "/ofrn.txt";
  FILE* ofrn = fopen(fname.c_str(),"wt+");
  fprintf(ofrn,"%d\n",run_number);
  fclose(ofrn);
  // end of ofrn.txt

  // == muon part ==

  char comment[100];
  sprintf(comment, "# nchannels\n%d\n", MAXOD);

  // lead time histograms and index // file mltcd.txt and mltci.txt
  fname = path + "/mltcd.txt";
  FILE* fpmltcd = fopen(fname.c_str(), "wt+");
  fname = path + "/mltci.txt";
  FILE* fpmltci = fopen(fname.c_str(), "wt+");
  fprintf(fpmltcd,comment);
  for(int32_t iCh=0; iCh<MAXOD; iCh++) {
    fprintf(fpmltcd,"# OD Logical Channel %04d\n",iCh+3001);
    fprintf(fpmltci,"%04ld\n",ftell(fpmltcd) );
    MyDump1DHisto(ODLeadTimeLCh[iCh],fpmltcd);
  }
  fclose(fpmltcd);
  fclose(fpmltci);

  // trail time histograms and index // file mttcd.txt and mttci.txt
  fname = path + "/mttcd.txt";
  FILE* fpmttcd = fopen(fname.c_str(), "wt+");
  fname = path + "/mttci.txt";
  FILE* fpmttci = fopen(fname.c_str(), "wt+");
  fprintf(fpmttcd,comment);
  for(int32_t iCh=0; iCh<MAXOD; iCh++) {
    fprintf(fpmttcd,"# OD Logical Channel %04d\n",iCh+3001);
    fprintf(fpmttci,"%04ld\n",ftell(fpmttcd) );
    MyDump1DHisto(ODTrailTimeLCh[iCh],fpmttcd);
  }
  fclose(fpmttcd);
  fclose(fpmttci);

  // charge histograms and index // file mccd.txt and mcci.txt
  fname = path + "/mccd.txt";
  FILE* fpmccd = fopen(fname.c_str(), "wt+");
  fname = path + "/mcci.txt";
  FILE* fpmcci = fopen(fname.c_str(), "wt+");
  fprintf(fpmccd,comment);
  for(int32_t iCh=0; iCh<MAXOD; iCh++) {
    fprintf(fpmccd,"# OD Logical Channel %04d\n",iCh+3001);
    fprintf(fpmcci,"%04ld\n",ftell(fpmccd) );
    MyDump1DHisto(ODChargeLCh[iCh],fpmccd);
  }
  fclose(fpmccd);
  fclose(fpmcci);

  // countrate history histograms and index // file mhcd.txt and mhci.txt
  fname = path + "/mhcd.txt";
  FILE* fpmhcd = fopen(fname.c_str(), "wt+");
  fname = path + "/mhci.txt";
  FILE* fpmhci = fopen(fname.c_str(), "wt+");
  fprintf(fpmhcd,comment);
  for(int32_t iCh=0; iCh<MAXOD; iCh++) {
    fprintf(fpmhcd,"# OD Logical Channel %04d\n",iCh+3001);
    fprintf(fpmhci,"%04ld\n",ftell(fpmhcd) );
    for(int32_t i=0; i < ODHistoryLCh[iCh]->GetNbinsX(); i++)
      fprintf(fpmhcd,"%s %06.2f\n",ODHistoryLCh[iCh]->GetXaxis()->GetBinLabel(i),(float)ODHistoryLCh[iCh]->GetBinContent(i));
  }
  fclose(fpmhcd);
  fclose(fpmhci);

  // occupancy // file mhpt.txt
  fname = path + "/mhpt.txt";
  FILE* mfpoccupancy = fopen(fname.c_str(), "wt+");
  fprintf(mfpoccupancy,"## pmt occupancy\n");
  fprintf(mfpoccupancy,"%d\n",MAXOD-1);
  for(int32_t i=1; i<MAXOD; i++) {
    fprintf(mfpoccupancy,"%d %6.2f\n",i+3000,MuonPmtOccupancy[i]);
  }
  fclose(mfpoccupancy);

  // muon trigger
  // file mmth.txt
  fname = path + "/mmth.txt";
  FILE* mfptriggerhistory = fopen(fname.c_str(), "wt+");
//  for(int32_t i=0; i < muTriggerRateHistory->GetNbinsX(); i++)
//  if (muTriggerRateHistory->GetXaxis()->GetBinLabel(i) != "") {
//    fprintf(mfptriggerhistory,"%s %06.2f\n",muTriggerRateHistory->GetXaxis()->GetBinLabel(i), (float)muTriggerRateHistory->GetBinContent(i));
//  }
  for(int32_t i=0; i < TriggerRateHistory->GetNbinsX(); i++)
//  if (TriggerRateHistory->GetXaxis()->GetBinLabel(i) != "") {
  if (strlen(TriggerRateHistory->GetXaxis()->GetBinLabel(i)) != 0) {
    fprintf(mfptriggerhistory,"%s %06.2f\n",TriggerRateHistory->GetXaxis()->GetBinLabel(i), (float)Trigger2RateHistory->GetBinContent(i));
  }
  fclose(mfptriggerhistory);

  // muon trigger history
  n = 22;
  for(int32_t i=0; i < muTriggerRateHistory->GetNbinsX(); i++)
    if (muTriggerRateHistory->GetBinContent(i) != 0) {
      muTriggerRateHistogram->Fill((float)muTriggerRateHistory->GetBinContent(i));
    }
  fname = path + "/mmtd.txt";
  FILE* mfptriggerhistogram = fopen(fname.c_str(), "wt+");
  fprintf(mfptriggerhistogram,"%d\n",n);
  for(int32_t i=0; i < n; i++)
    fprintf(mfptriggerhistogram,"%06.2f %06.2f\n", (float)muTriggerRateHistogram->GetBinCenter(i), (float)muTriggerRateHistogram->GetBinContent(i));
  fprintf(mfptriggerhistogram,"%d\n%5.1f\n%5.1f\n",
    (int32_t)muTriggerRateHistogram->GetEntries(), (float)muTriggerRateHistogram->GetMean(), (float)muTriggerRateHistogram->GetRMS() );
  fclose(mfptriggerhistogram);
//  muTriggerRateHistogram->Delete(); // не уничтожать, а обнулять... или вообще отказаться "временно"

  // file mmth.txt
  fname = path + "/mnhh.txt";
  FILE* mfpnh = fopen(fname.c_str(), "wt+");
  for(int32_t i=0; i < TriggerRateHistory->GetNbinsX(); i++)
//    if (TriggerRateHistory->GetXaxis()->GetBinLabel(i) != "") {
    if (strlen(TriggerRateHistory->GetXaxis()->GetBinLabel(i)) != 0) {
      fprintf(mfpnh,"%s %06.2f\n",TriggerRateHistory->GetXaxis()->GetBinLabel(i), (float)muHitsTotalHistory->GetBinContent(i) / (float)TriggerRateHistory->GetBinContent(i) / 60.);
    }

  fclose(mfpnh);

  // muon trigger distribution // file mtd.txt
  fname = path + "/mnhd.txt";
  FILE* fpmnhd = fopen(fname.c_str(), "wt+");
  fprintf(fpmnhd,"%d\n",muTriggerRateHistogram->GetNbinsX());
  for(int32_t i=0; i < muTriggerRateHistogram->GetNbinsX(); i++)
    fprintf(fpmnhd,"%06.2f %06.2f\n",(float)muTriggerRateHistogram->GetBinCenter(i), (float)muTriggerRateHistogram->GetBinContent(i));
  fprintf(fpmnhd,"%d\n%5.1f\n%5.1f\n",
    (int32_t)muTriggerRateHistogram->GetEntries(), (float)muTriggerRateHistogram->GetMean(), (float)muTriggerRateHistogram->GetRMS() );
  fclose(fpmnhd);

  // mu-Hits-Total-History // file nhh.txt
  // mu-Hits-Total-Distribution // file mnhd.txt
  // file msad.txt (muon sum of amplitudes -- distribution)
  fname = path + "/msad.txt";
  FILE* mfpmsad = fopen(fname.c_str(), "wt+");
  fprintf (mfpmsad, "#muon sum of amplitudes -- distribution\n");
  fprintf(mfpmsad,"%d\n",ODQuasyEnergy->GetNbinsX());
  for(int32_t i=1; i <= ODQuasyEnergy->GetNbinsX(); i++) {
    fprintf(mfpmsad,"%06.2f %06.2f\n",(float)ODQuasyEnergy->GetBinCenter(i), (float)ODQuasyEnergy->GetBinContent(i));
  }
  fprintf(mfpmsad,"%d\n%5.1f\n%5.1f\n",
    (int32_t)ODQuasyEnergy->GetEntries(), (float)ODQuasyEnergy->GetMean(), (float)ODQuasyEnergy->GetRMS() );
  fclose(mfpmsad);

  // == fwfd part ==

  // fwfd pedestal channel // files fphd.txt fphi.txt fmpm.txt
  fname = path + "/fphd.txt";
  FILE* fphd = fopen(fname.c_str(),"wt+");
  fname = path + "/fphi.txt";
  FILE* fphi = fopen(fname.c_str(),"wt+");
  fname = path + "/fmpm.txt";
  FILE* fmpm = fopen(fname.c_str(),"wt+");
  for(int32_t i=0; i<MAXFWFD; i++) {
    fprintf(fphd,"# fwfd pedestal channel %03d\n",i);
    fprintf(fphi,"%04ld\n",ftell(fphd));
    MyDump1DHisto(PedCh[i], fphd);
    fprintf(fmpm,"%03d %5.1f\n",i+1,(float)PedCh[i]->GetMean());
  };
  fclose(fphd);
  fclose(fphi);
  fclose(fmpm);
  
  // fwfd peak amplitude channel // file fpad.txt fpai.txt
  fname = path + "/fpad.txt";
  FILE* fpad = fopen(fname.c_str(),"wt+");
  fname = path + "/fpai.txt";
  FILE* fpai = fopen(fname.c_str(),"wt+");
  for(int32_t i=1; i<=MAXFWFD; i++) {
    fprintf(fpad,"# fwfd peak amplitude channel %03d\n",i);
    fprintf(fpai,"%04ld\n",ftell(fpad));
    MyDump1DHisto(AmplCh[i],fpad);
  };
  fclose(fpad);
  fclose(fpai);

  // fwfd peak position channel // files fppd.txt fppi.txt
  fname = path + "/fppd.txt";
  FILE* fppd = fopen(fname.c_str(),"wt+");
  fname = path + "/fppi.txt";
  FILE* fppi = fopen(fname.c_str(),"wt+");
  for(int32_t i=1; i<=MAXFWFD; i++) {
    fprintf(fppd,"# fwfd peak position channel %03d\n",i);
    fprintf(fppi,"%04ld\n",ftell(fppd));
    MyDump1DHisto(PeakCh[i],fppd);
  };
  fclose(fppd);
  fclose(fppi);
  
  // fwfd pealk charge channel // files fpqd.txt fpqi.txt
  fname = path + "/fpqd.txt";
  FILE* fpqd = fopen(fname.c_str(),"wt+");
  fname = path + "/fpqi.txt";
  FILE* fpqi = fopen(fname.c_str(),"wt+");
  for(int32_t i=1; i<=MAXFWFD; i++) {
    fprintf(fpqd,"# fwfd peak charge channel %03d\n",i);
    fprintf(fpqi,"%04ld\n",ftell(fpqd));
    MyDump1DHisto(ChargeCh[i],fpqd);
  };
  fclose(fpqd);
  fclose(fpqi);
  
  // peak Num Sum // file fsph.txt
  fname = path + "/fsph.txt";
  FILE* fsph = fopen(fname.c_str(),"wt+");
  MyDump1DHisto(PeakNSum,fsph);
  fclose(fsph);
  
  // ampl num sum // file fsah.txt
  fname = path + "/fsah.txt";
  FILE* fsah = fopen(fname.c_str(),"wt+");
  MyDump1DHisto(AmplNSum,fsah);
  fclose(fsah);
  
  // charge num sum // file fsqh.txt
  fname = path + "/fsqh.txt";
  FILE* fsqh = fopen(fname.c_str(),"wt+");
  MyDump1DHisto(ChargeNSum,fsqh);
  fclose(fsqh);
  
  // number of windows // file fowh.txt
  fname = path + "/fowh.txt";
  FILE* fowh = fopen(fname.c_str(),"wt+");
  MyDump1DHisto(NWin,fowh);
  fclose(fowh);
  
  // event flag // file fofh.txt
  fname = path + "/fofh.txt";
  FILE* fofh = fopen(fname.c_str(),"wt+");
  MyDump1DHisto(EvFlag,fofh);
  fclose(fofh);
  
  // t_trigger ? // file foth.txt
  fname = path + "/foth.txt";
  FILE* foth = fopen(fname.c_str(),"wt+");
  MyDump1DHisto(TTrig,foth);
  fclose(foth);
  
  // n_pulse ? // file fonh.txt
  fname = path + "/fonh.txt";
  FILE* fonh = fopen(fname.c_str(),"wt+");
  MyDump1DHisto(NPulse,fonh);
  fclose(fonh);
  
  // total size // file fooh.txt
  fname = path + "/fooh.txt";
  FILE* fooh = fopen(fname.c_str(),"wt+");
  MyDump1DHisto(TotSize,fooh);
  fclose(fooh);

  return;
}

//-----------------------------------------------------------------------------
// dump a root histogram into an ascii file
void bx_omon::MyDump1DHisto(TH1F* th, FILE* fp) {
  if (!th || !fp ) return;
//  get_message (bx_message::debug) << "Dump1D" << dispatch;
  fprintf(fp,"%6.2f\n%6.2f\n",(float)th->GetXaxis()->GetXmin(), (float)th->GetXaxis()->GetXmax());
  for(int32_t i=0; i < th->GetNbinsX(); i++)
    fprintf(fp,"%06.2f\n",(float)th->GetBinContent(i));
  fprintf(fp,"%d\n%5.1f\n%5.1f\n",(int32_t)th->GetEntries(), (float)th->GetMean(), (float)th->GetRMS() );
  fprintf(fp,"\n");
  return;
}

//-----------------------------------------------------------------------------
// reset histograms
void bx_omon::ResetHistos() {
  get_message (bx_message::info) << "OMON resetting histograms" << dispatch;

  for(int32_t i=0; i<MAXCHID; i++) {
    RawChargeLCh[i]->Reset("ICE");
    RawBaseLCh[i]->Reset("ICE");
    ChargeCutLCh[i]->Reset("ICE");
    TimeLCh[i]->Reset("ICE");
    HistoryLCh[i]->Reset("ICE");
  };

  DistribAverTime->Reset("ICE");
  DistribRefTime->Reset("ICE");
  TriggerRateHistory->Reset("ICE");
  EventSizeHistory->Reset("ICE");
  EventSizeDistribution->Reset("ICE");
  HitsTotalHistory->Reset("ICE");
  HitsTotalDistribution->Reset("ICE");
  HitsWConeHistory->Reset("ICE");
  HitsWConeDistribution->Reset("ICE");
  HitsNConeHistory->Reset("ICE");
  HitsNConeDistribution->Reset("ICE");
  TriggerTypeDistribution->Reset("ICE");
  EventsTimeDistribution->Reset("ICE");
  GPS_ppc0_Distribution->Reset("ICE");
  GPS_build_Distribution->Reset("ICE");
  QuasyEnergy->Reset("ICE");
  QuasyEnergyT01->Reset("ICE");
  QuasyEnergyT08->Reset("ICE");
  QuasyEnergyT32->Reset("ICE");
  QuasyEnergyT64->Reset("ICE");
  QuasyEnergyCutted->Reset("ICE");
  TriggerRateHistogram->Reset("ICE");

  for(int32_t i=0; i<MAXOD; i++) {
    ODLeadTimeLCh[i]->Reset("ICE");
    ODTrailTimeLCh[i]->Reset("ICE");
    ODChargeLCh[i]->Reset("ICE");
    ODHistoryLCh[i]->Reset("ICE");
  }
  ODQuasyEnergy->Reset("ICE");
  muTriggerRateHistory->Reset("ICE");
  muTriggerRateHistogram->Reset("ICE");
  muHitsTotalHistory->Reset("ICE");
  muHitsTotalDistribution->Reset("ICE");
  EvFlag->Reset("ICE");
  TTrig->Reset("ICE");
  for(int32_t i=0; i<MAXFWFD; i++) {
    PedCh[i]->Reset("ICE");
    PeakCh[i]->Reset("ICE");
    AmplCh[i]->Reset("ICE");
    ChargeCh[i]->Reset("ICE");
  }
  PeakNSum->Reset("ICE");
  AmplNSum->Reset("ICE");
  ChargeNSum->Reset("ICE");
  NWin->Reset("ICE");
  NPulse->Reset("ICE");
  TotSize->Reset("ICE");
}

/*
 * $Log: bx_omon.cc,v $
 * Revision 1.27  2009/09/16 16:17:23  sabelnik
 * probably fix warnings
 *
 * Revision 1.26  2008-09-01 17:12:11  sabelnik
 * OD: trg rate history & hits/trg history
 *
 * Revision 1.25  2008-08-29 14:42:00  sabelnik
 * OD  Channel Histograms
 *
 * Revision 1.24  2008-08-23 11:48:04  sabelnik
 * OD mapping errors fixed
 *
 * Revision 1.23  2008-06-30 09:12:48  sabelnik
 * OD map in OMON
 *
 * Revision 1.22  2008-03-31 15:20:32  sabelnik
 * trigger types history
 *
 * Revision 1.21  2008-02-20 15:44:02  pallas
 * Fix a bug. The single channel histogram was displaying "peak" instead
 * of the correct charge ( peak - base )
 *
 * Revision 1.20  2008-01-04 09:26:52  pallas
 * A small change to allow trigger type selection in online monitor histograms
 * Also, update law based on time implemented
 *
 * Revision 1.19  2007-04-02 18:14:38  sabelnik
 * triger types
 *
 * Revision 1.17  2007-04-01 19:46:28  sabelnik
 * energy nhits
 *
 * Revision 1.14  2006-11-24 11:52:40  sabelnik
 * Trigger Types Distribution filled
 *
 * Revision 1.13  2006/11/23 20:11:58  sabelnik
 * some more histogramms filled
 *
 * Revision 1.12  2006/11/22 19:41:48  sabelnik
 * 2 more histos filled
 *
 * Revision 1.11  2006/11/22 16:12:39  sabelnik
 * history of trigger rates filled
 *
 * Revision 1.10  2006/09/04 16:04:40  sabelnik
 * export all histos (not filled!)
 *
 * Revision 1.9  2006/09/01 17:25:56  sabelnik
 * export all laben and fwfd histograms (most not filled!)
 *
 * Revision 1.8  2006/08/31 16:35:26  sabelnik
 * all histos allocated and most saved (not filled!)
 *
 * Revision 1.6  2006/08/26 15:11:29  sabelnik
 * some histos init added
 *
 * Revision 1.5  2006/08/21 11:27:00  razeto
 * Updated to new barn_interface
 *
 * Revision 1.4  2006/08/18 15:22:46  sabelnik
 * charge histograms added
 *
 * Revision 1.3  2006/08/17 15:45:08  sabelnik
 * complete list of output files
 *
 * Revision 1.2  2006/08/08 11:21:56  sabelnik
 * detector and electronics maps implemented.
 * OMON web pages are generated from now on
 * by means of Echidna-based histogram server.
 *
 * Revision 1.1  2006/07/13 16:59:53  sabelnik
 * first entry
 *
 * Revision 1.0  2006/05/25 16:21:03 sabelnik
 * Derivated from bx_test_module
 */
