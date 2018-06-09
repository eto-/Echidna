// Usage: ROOT>.x wf.C(RUN, EVENT, FADC_OR_LABEN_EVNUM, SUM)
// RUN - run number, EVENT - event number,
// FADC_OR_LABEN_EVNUM = 0 - FADC evnum
// FADC_OR_LABEN_EVNUM = 1 - LABEN evnum
// SUM = 0 - draw ASUM
// SUM = 1 - draw DSUM

#include <cstdio>
#include <cstring>
#include <iostream>
#include <fstream>

//#include "TFile.h"
//#include "TH1F.h"
//#include "TTree.h"

#define RUN 		22010
#define EVENT 		790

//#include "/storage/gpfs_data/borexino/users/litvinov/Echidna_c16/event/BxEvent.hh"
#include "fadc_aux.cc"
#include "muons_ps.cc"
#include "noise_ps.cc"

void wf(int run = RUN, int event = EVENT, int FADC_OR_LABEN_EVNUM = 0)  
{

  TH1F* Awf = new TH1F("Awf","Awf",512,0,512);
  TH1F* Dwf = new TH1F("Dwf","Dwf",512,0,512);
  TCanvas *c = new TCanvas("c","ADSum",1000,800);
  c->Divide(1,2);


    string run_name = get_fadc_run_path(run);
    char tmp[1000];
    //sprintf(tmp, "/home/lukyanch/tmp/Run%06d_fadc_c16.root", run);
    //string run_name = tmp;
    
    if (run_name.empty())
    {
        printf("Run is not validated!\n");
        return;
    }
    
    TFile* f = TFile::Open(run_name.c_str());
    
    if (!f || !f->IsOpen())
        return;
    
//    BxEvent* ev = new BxEvent();
    BxFwfd* ev = new BxFwfd();
    
    TTree* tree = (TTree*)f->Get("BxFwfd");
    tree->SetBranchAddress("bxfwfd",&ev);
    
    int evnum_bx, evnum, trgtype, nclu;
    double time_prev, bc;
    float charge;
    char name[10], title[80];
    
//    float tail12 = 0., tail16 = 0., tail20 = 0.;

    int nentries = tree->GetEntries();
    
    if (FADC_OR_LABEN_EVNUM == 1)  {// search FADC waveform through evnum_bx

    for (int i = 0; i < nentries; i++)
    {
    	tree->GetEntry(i);
    	evnum = ev->GetEvNum();
    	evnum_bx = ev->GetEvNumBx()+1;
    	trgtype = ev->GetTrgType();
    	nclu = ev->GetNClusters();
	
	if (evnum_bx == event)
	{
    	    sprintf(title,"run %d, FADC ev. %d, BX ev. %d",run,evnum,evnum_bx);

	    printf("FADC ev. %d, BX ev. %d, FADC trg=%d, nclu = %d\n",evnum,evnum_bx,trgtype,nclu);
	    
	    for (int j = 0; j < 512; j++){ 
	      Awf->Fill(j,ev->GetWFormAsum(j));
	      Dwf->Fill(j,ev->GetWFormDsum(j));
	    }
    
	float dcharge = ev->GetCluster(0).GetDSumCharge();
    
    	for (int k = 0; k < nclu; k++)  {
             printf("Clu %d: pos=%d raw_time=%u ns",k,ev->GetCluster(k).GetPeakPos(),ev->GetRawTime());
             printf(" EA=%2.1f MeV ED=%2.1f MeV\n",ev->GetCluster(k).GetASumCharge()/378.,ev->GetCluster(k).GetDSumCharge()/352.);
            }

        if (is_muon(ev,tree) > 0)  printf("\n***** MUON-like event *****\n");
        if (is_noise(ev,0) > 0)  printf("\n***** NOISE-like event *****\n");
	
	Awf->SetTitle(title);
	Dwf->SetTitle(title);	
	c->cd(1);
    	Awf->Draw();
	c->cd(2);
    	Dwf->Draw();
	break;
     	}
	else if (evnum_bx > event)
	{
	     printf("NO LABEN event %d in FADC data found\n",event);
	     break;
	}
    }
    }

    else  {// search FADC waveform through FADC evnum 

    tree->GetEntry(event-1);
    
    //run = ev->GetRun();
    
    evnum = ev->GetEvNum();
    evnum_bx = ev->GetEvNumBx()+1;
    trgtype = ev->GetTrgType();
    nclu = ev->GetNClusters();
    
    sprintf(title,"run %d, FADC ev. %d, BX ev. %d",run,evnum,evnum_bx);
    printf("FADC ev. %d, BX ev. %d, FADC trg=%d, nclu = %d\n",evnum,evnum_bx,trgtype,nclu);
    
    for (int j = 0; j < 512; j++){
      Awf->Fill(j,ev->GetWFormAsum(j));
      Dwf->Fill(j,ev->GetWFormDsum(j));
    }
    float dcharge = ev->GetCluster(0).GetDSumCharge();
    
    for (int k = 0; k < nclu; k++)  {
        printf("Clu %d: pos=%d raw_time=%u ns",k,ev->GetCluster(k).GetPeakPos(),ev->GetRawTime());
        printf(" EA=%2.1f MeV ED=%2.1f MeV\n",ev->GetCluster(k).GetASumCharge()/378.,ev->GetCluster(k).GetDSumCharge()/352.);
    }
        
    if (is_muon(ev,tree) > 0)  printf("\n***** MUON-like event *****\n");
    if (is_noise(ev,0) > 0)  printf("\n***** NOISE-like event *****\n");
 
    ev->Clear();
     Awf->SetTitle(title);
    Dwf->SetTitle(title);
    c->cd(1);
    Awf->Draw();
    c->cd(2);
    Dwf->Draw();
   }
}

