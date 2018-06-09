/*
FADC data validation script

Once compiled (.L fadc_validation.C+), this macro is launched:

ROOT>fadc_validation("12345"), where 12345 is FADC run number
ROOT>save() // if you decide to validate the run

mailto: litvinovich@lngs.infn.it
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include <iostream>
#include <fstream>
#include <iomanip>

#include "TROOT.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TPaveStats.h"

//#include "/home/evgeny/include/FWFDEvent.hh"

#include "muons_ps.cc" // by Egor
#include "Conversion_factors.cc" // by Maxim


using namespace std;

int is_muon (BxFwfd *ev, TTree* tree);
double GetConversionFactor(int start_time, int sum); // returns ASum if sum=1, DSum if sum=2

//const int DISC_CHANNELS = 16; // declared in muons_ps.cc
const int MAJ_LEVEL = 4;

// Global variables to be saved in FADC_run_validation_out.txt
int run, nentries, n_antinues=0, n_errors=0, n_dcode_errors=0;
int start_time, end_time, runtime;
float n_muons_daily, percent_muons_ps, percent_1280, percent_noise;
float n_1MeV_Hz;

// Function seeking for antineutrino candidates
int is_antinue(int i, BxFwfd *ev, TTree* tree) {

  double time_to_muon = 0;
  float first, last;

  tree->GetEntry(i-1);

  if (ev->GetNClusters() == 0)  return 0;
  
  float asum_factor = GetConversionFactor(start_time,1);
  float dsum_factor = GetConversionFactor(start_time,2);
  
  float eneA1 = ev->GetCluster(0).GetASumCharge() / asum_factor; // MeV
  float eneD1 = ev->GetCluster(0).GetDSumCharge() / dsum_factor; // MeV

  int nclu1 = ev->GetNClusters();
  int trgtype1 = ev->GetTrgType();
//  int evnum = ev->GetEvNum();
//  int evnum_bx = ev->GetEvNumBx()+1;

  ev->Clear();

  // Positron?
  unsigned int raw_time = ev->GetRawTime(); // ns
  if (trgtype1 == 1 && nclu1 == 1 && raw_time > 1280 && eneA1 > 1. && eneD1 > 1. && eneD1 < 50.)  {

// Check for fake event with 1 cluster
      bool fake_event_1clu = false;
      if (nclu1 == 1 && raw_time < 1280 && ev->GetCluster(0).GetPeakAmpl() > 210)
          fake_event_1clu = true;

      if (fake_event_1clu)  return 0;

      // Quality cut on waveform
      int peak_pos = ev->GetCluster(0).GetPeakPos();
      float peak_ampl = ev->GetCluster(0).GetPeakAmpl();
      if (peak_pos > 9)  first = ev->GetWFormAsum(peak_pos-10);
      else first = ev->GetWFormAsum(0);
      if (peak_pos < 492)  last = ev->GetWFormAsum(peak_pos+20);
      else last = ev->GetWFormAsum(511);

      if (first < peak_ampl/4. && last < peak_ampl/3.)  {

      time_to_muon += ev->GetRawTime()*1e-6; // ms
      if (time_to_muon > 3.)  return 1; // 3 ms
      
      for (int k = i-2; k > i-500; k--)  {
      
           if (k < 0)  return 0;

	   tree->GetEntry(k);
	   if (is_muon(ev,tree) > 0)  return 0;

           time_to_muon += ev->GetRawTime()*1e-6; // ms
/*           evnum = ev->GetEvNum();
           evnum_bx = ev->GetEvNumBx()+1;
           eneD1 = ev->GetCluster(0).GetDSumCharge() / dsum_factor; // MeV
*/
	   ev->Clear();

           if (time_to_muon > 3.)  return 1; // 3 ms
	  }
     }
      else return 0;

     }
      else return 0;

  return 1; // no muon found within 3 ms before, so it's probably a real antinue candidate
}

//____________________________________________________

int is_dcode_error(BxFwfd* ev)  {

  int dcode_ch = 0;
            
  for (int c = 0; c < DISC_CHANNELS; ++c)
      {
       if (/*c != 7&&*/ ev->GetDCode(c))
           ++dcode_ch;
      }
            
  if (dcode_ch < MAJ_LEVEL + 1)  return 1; // +1 because ch. 8 is always fired
  else return 0;

}

//____________________________________________________

void fadc_validation(string runnum) {

  gStyle->SetMarkerStyle(7);    // small squared bullet
  gStyle->SetMarkerColor(50);   // red

  string fnameFadc;
//  fnameFadc = "/data1/fadc/rootfiles/Run0";
  fnameFadc = "/home/production/Echidna/Run0";

  fnameFadc += runnum;
  fnameFadc += "_fadc_c16.root";
  TFile* current_file = TFile::Open(fnameFadc.c_str(),"READ");
  
// Check if run already validated
  string str = "grep Run0";
  str += runnum;
  str += " /home/production/FADC_run_validation_out.txt";
//  str += " /home/evgeny/fadc/FADC_run_validation_out.txt";

  char is_validated[200];
  strcpy(is_validated,str.c_str());

  char no[3];
  if ( !system(is_validated) )  {
      printf("ATTENTION: run already validated. Continue? (y or n)\n");
      gets(no);
      if (no[0] == 'n')  exit(1);
     }

  if (!current_file)  {
      cerr << "Couldn't open file " << fnameFadc << endl;
      exit(1);
     }

  if (current_file->IsZombie())   {
      if (current_file->IsOpen())  current_file->Close();
      cout << "File is zombi: " << fnameFadc << endl;
      exit(1);
     }

  BxFwfd* ev = new BxFwfd();

  TTree* tree = (TTree*)current_file->Get("BxFwfd");
  tree->SetBranchAddress("bxfwfd",&ev);

  int evnum, evnum_bx, trgtype, nclu, ntrg1=0, n_1280=0, n_noise = 0, n_muons=0, n_muons_ps=0;
  int n_1MeV=0;
  UInt_t current_time, time_to_muon;

  int evnum_list[100],evnum_bx_list[100]; // of positron
  float eneD1_list[100]; // of positron

  int dcode[16] = {0};

  nentries = tree->GetEntries();

// Get time of the last event
  tree->GetEntry(nentries - 1);
  end_time = ev->GetUnixTime();
  ev->Clear();

// Get time of the last event
  tree->GetEntry(0);
  start_time = ev->GetUnixTime();
  ev->Clear();

  float asum_factor = GetConversionFactor(start_time,1);
  float dsum_factor = GetConversionFactor(start_time,2);

  int nbins = int((end_time - start_time) / 30); // 30 s bin
  int nbins_mu = int((end_time - start_time) / 60); // 60 s bin

  TH1F* hTemp = new TH1F("hTemp","Temporaty histo",nbins,start_time,end_time);
  TH1F* hTempMu = new TH1F("hTempMu","Temporaty histo for muons",nbins_mu,start_time,end_time);
//  TH1F* hTrg_temp = new TH1F("hTrg_temp","temporary histo for trigger type",130,0,130);
  TH1F* hTrg = new TH1F("","",13,0,13);

  TH2F* hEv = new TH2F("","",nbins,start_time,end_time,nentries/10,0,nentries);
  TH1F* hEvBx = new TH1F("","",nentries+1,0,nentries+1);

  TH1F* hEneA = new TH1F("","",200,0,5);
  TH1F* hEneD = new TH1F("","",2000,0,50);

//  int max_trgtype = 0;
  char name[10], title[80];
  
  for (int i = 0; i < nentries; i++)  {

       tree->GetEntry(i);

// First event in run
       if (i == 0)  {

           run = ev->GetRun();
/*	   cout << endl;
	   cout << "Run number: " << run << endl;
           cout << "Events: " << nentries << endl;*/
           printf("\nRun Number: %d\n",run);
           printf("Events: %d\n",nentries);

	   TDatime da(start_time);
//	   cout << "Run started: " << da.GetDay() << "/" << da.GetMonth() << "/" << da.GetYear();
//	   cout << " " << da.GetHour() << ":" << da.GetMinute() << " GMT" << endl;
printf("Run started: %d",da.GetDay()); printf("/");
printf("%d",da.GetMonth()); printf("/");
printf("%d",da.GetYear());
printf(" %d",da.GetHour()); printf(":");
printf("%d",da.GetMinute()); printf(" GMT\n");

         }

	  evnum = ev->GetEvNum();
	  evnum_bx = ev->GetEvNumBx();//+1
	  nclu = ev->GetNClusters();

	  UInt_t raw_time = ev->GetRawTime();
	  if (raw_time < 1280)  n_1280++;

	  trgtype = ev->GetTrgType();
//	  if (trgtype > max_trgtype)  max_trgtype = trgtype;
//	  hTrg_temp->Fill(trgtype);
	  hTrg->Fill(trgtype);
	  
// Calculate the number of ADC memory overflows (error)
          if (ev->GetError() == 1)  n_errors++;

// Check if all discriminator channels are OK
          if (trgtype == 1)  {
	      ntrg1++;
	      for (int channel = 0; channel < 16; channel++) {
	           dcode[channel] += ev->GetDCode(channel);
		  }
             }

// Check for dcode_errors (if N fired discriminator channels < MAJ)
  if ((ev->GetTrgType() & 1) && (raw_time > 1280))  n_dcode_errors += is_dcode_error(ev);

// Fill counting rate
	  current_time = ev->GetUnixTime();
	  hTemp->Fill(current_time);
	  hEv->Fill(current_time,evnum);
	  hEvBx->Fill(evnum,evnum_bx);

// Check muons
          int muon_value = is_muon(ev,tree);
	  if (muon_value > 0 && muon_value != 1)  n_muons_ps++;
          if (muon_value > 0)  {
	      n_muons++;
	      time_to_muon = 0;

	      hTempMu->Fill(current_time);
//	      hEvMu->Fill(current_time,evnum);
	     }
	  if (muon_value == -3)  n_noise++; // !! noise finder adjusted to muons is used. Replace in future !!

// Draw energy spectrum
//	  else if ((trgtype == 1 || trgtype == 16 || trgtype == 17) && (nclu != 0) && (muon_value <= 0))  {
	  else if ((trgtype == 1 || trgtype == 17) && (nclu != 0) && (muon_value <= 0))  {

	      time_to_muon += ev->GetRawTime();

	      float first, last;
	      for (int clu = 0; clu < nclu; clu++)  {
	           float eneA = ev->GetCluster(clu).GetASumCharge() / asum_factor; // MeV
	           float eneD = ev->GetCluster(clu).GetDSumCharge() / dsum_factor; // MeV

                        // Quality cut on waveform
                   int peak_pos = ev->GetCluster(0).GetPeakPos();
                   float peak_ampl = ev->GetCluster(0).GetPeakAmpl();
                   if (peak_pos > 9) first = ev->GetWFormAsum(peak_pos-10);
                   else first = ev->GetWFormAsum(0);
                   if (peak_pos < 492) last = ev->GetWFormAsum(peak_pos+20);
                   else last = ev->GetWFormAsum(511);

//                   if (first < peak_ampl/4. && last < peak_ampl/3.)  {
                   if (first < peak_ampl/4. && last < peak_ampl/3. && raw_time > 1280.)  {

		       hEneA->Fill(eneA);
		       hEneD->Fill(eneD);
		       if (eneA > 1.0)  n_1MeV++;

// BiPo212 candidate:
//if (clu == 1 && nclu == 2 && eneD > 0.7 && eneD < 1.4)  cout << evnum << endl;
		      }
		  }

// Check if antineutrino candidate?
//          float eneD2 = ev->GetCluster(0).GetDSumCharge() / dsum_factor; // MeV
//          double time_prev = ev->GetRawTime()*1e-6; // ms

          // Neutron?
/*	  if (nclu == 1 && time_prev < 3. && eneD2 > 1.3 && eneD2 < 3.)  {

	      // Candidate found
	      if (is_antinue(i, ev, tree) == 1)  {

		  tree->GetEntry(i-1);

                  eneD1_list[n_antinues] = ev->GetCluster(0).GetDSumCharge() / dsum_factor; // MeV
                  evnum_list[n_antinues] = ev->GetEvNum();
                  evnum_bx_list[n_antinues] = ev->GetEvNumBx()+1;

	          n_antinues += 1;
                  ev->Clear();
		 }
	     }*/
	     }

       ev->Clear();
      }

// Unixtime of the last event in run
  TDatime da2(current_time);
//  cout << "Run ended: " << da2.GetDay() << "/" << da2.GetMonth() << "/" << da2.GetYear();
//  cout << " " << da2.GetHour() << ":" << da2.GetMinute() << " GMT" << endl;

printf("Run ended: %d",da2.GetDay()); printf("/");
printf("%d",da2.GetMonth()); printf("/");
printf("%d",da2.GetYear());
printf(" %d",da2.GetHour()); printf(":");
printf("%d",da2.GetMinute()); printf(" GMT\n");

  runtime = current_time - start_time;
//  cout << "Run duration: " << runtime/3600. << " hours";
  printf("Run duration: %f hours ",runtime/3600.);
  if (runtime < 900)  printf("ATTENTION: run is too short and should not be validated\n\n");
  else if (runtime < 43200.)  printf("OK\n\n");
  else if (runtime > 43200.)  printf("ATTENTION: run is more than 12h long!\n\n");
/*  if (runtime < 900)  {
      cout << endl;
      cout << "=================" << endl;
      cout << "ATTENTION: run is too short and should not be validated" << endl;
      cout << "=================" << endl;
     }      
  else if (runtime < 43200.)  cout << " OK" << endl; // 12 hours

  else if (runtime > 43200.)  {
      cout << endl;
      cout << "=================" << endl;
      cout << "ATTENTION: run is more than 12h long!" << endl;
      cout << "=================" << endl;
     }
*/
//  cout << endl;
//  cout << "===== MUON'S RATE CHECK =====" << endl;
  printf("===== MUON'S RATE CHECK =====\n");
  n_muons_daily = n_muons*86400./runtime;
//  cout << n_muons_daily << " mu/day (2800 to 5800 expected)";
  printf("%.1f mu/day (2800 to 5800 expected)",n_muons_daily);

  if (n_muons_daily > 2800. && n_muons_daily < 5800.)
//  cout << " OK" << endl;
  printf(" OK\n");
  else printf(" !! Unexpected muon rate !!\n");

//  cout << n_muons_ps << " muons without OD trigger (identified by its pulse shape)" << endl;
  percent_muons_ps = (float)100.*n_muons_ps/n_muons;
  printf("%d (%.1f%) muons without OD trigger (identified by its pulse shape)\n",n_muons_ps,percent_muons_ps);

//____________ Discriminator channels_____________

//  cout << endl;
//  cout << "DISCRIMINATOR CHANNELS (8th may contain 1.0, 14th and 16th - 0.0):" << endl;
//  cout << " 1    2    3    4    5    6    7    8    9    10   11   12   13   14   15   16  " << endl;

printf("\nDISCRIMINATOR CHANNELS (8th may contain 1.0, 14th and 16th - 0.0):\n");
printf(" 1    2    3    4    5    6    7    8    9    10   11   12   13   14   15   16  \n");

       printf("%.1f  %.1f  %.1f  %.1f  %.1f  %.1f  %.1f  %.1f  %.1f  %.1f  %.1f  %.1f  %.1f  %.1f  %.1f  %.1f  hits/event\n",
(float)dcode[0]/ntrg1,
(float)dcode[1]/ntrg1,
(float)dcode[2]/ntrg1,
(float)dcode[3]/ntrg1,
(float)dcode[4]/ntrg1,
(float)dcode[5]/ntrg1,
(float)dcode[6]/ntrg1,
(float)dcode[7]/ntrg1,
(float)dcode[8]/ntrg1,
(float)dcode[9]/ntrg1,
(float)dcode[10]/ntrg1,
(float)dcode[11]/ntrg1,
(float)dcode[12]/ntrg1,
(float)dcode[13]/ntrg1,
(float)dcode[14]/ntrg1,
(float)dcode[15]/ntrg1);

 
//___________ N dcode errors ____________

  printf("\nN dcode errors: %d ",n_dcode_errors);
  if (n_dcode_errors > 5)  printf(" !! ATTENTION !! Too many dcode errors\n");
  else printf(" OK\n");

//___________ N memory overflows ____________

  printf("\nN memory overflows: %d ",n_errors);
  if (n_errors > 0)  printf(" !! ATTENTION !! Some events could be lost\n");
  else printf(" OK\n");

//___________ Noise events (is_muon == -3) ___________

  percent_noise = (float)100.*n_noise/nentries;
  printf("\nNoise events: %d (%.1f%)",n_noise,percent_noise);
  if (percent_noise > 50.)  printf(" !! ATTENTION !! Too high noise level\n");
  else printf(" OK\n");

//___________ Overlapped events (raw_time < 1280 ns) ___________

/*  percent_1280 = (float)100.*n_1280/nentries;
  printf("\nOverlapped events: %d (%.1f%)",n_1280,percent_1280);
  if (percent_1280 > 50.)  printf(" !! ATTENTION !! Too many overlapped windows\n");
  else printf(" OK\n");
*/

//___________ Rate above 1 MeV ___________

  n_1MeV_Hz = (float)n_1MeV/runtime;
  printf("\nRate above 1 MeV: %f Hz",n_1MeV_Hz);
  if (n_1MeV_Hz < 0.01)  printf(" !! ATTENTION !! Too low rate\n"); // down to 0.005 from 0.01 on July 9, 2014 (EL), back to 0.01 on Sep 19, 2014 (EL)
  else if (n_1MeV_Hz > 0.04)  printf(" Too high rate?\n");
  else printf(" OK\n");

//___________________Antineutrino candidates___________________

//  cout << endl;
//  cout << "ANTINEUTRINO CANDIDATES:" << endl;
//  printf("\nANTINEUTRINO CANDIDATES:\n");

//  for (int d = 0; d < n_antinues; d++)  printf("Positron ev. %d, BX ev. %d, ED=%.1f MeV\n",evnum_list[d],evnum_bx_list[d],eneD1_list[d]);
//  cout << n_antinues << " candidates found" << endl;
//  printf("%d candidates found\n",n_antinues);

//_________________________ Histograms ________________________

  sprintf(name,"hRate");
  sprintf(title,"Run %d, FADC counting rate, Hz",run);
  TH1F* hRate = new TH1F("","",nbins,0,nbins);
  hRate->SetNameTitle(name,title);

  sprintf(name,"hRateMu");
  sprintf(title,"Run %d, Muons counting rate, Hz",run);
  TH1F* hRateMu = new TH1F("","",nbins_mu,0,nbins_mu);
  hRateMu->SetNameTitle(name,title);

  for (int r = 0; r < nbins; r++)  hRate->Fill(r,hTemp->GetBinContent(r)/30.); // 30 s bins
  for (int rr = 0; rr < nbins_mu; rr++)  hRateMu->Fill(rr,hTempMu->GetBinContent(rr)/60.); // 60 s bins

  hRate->Sumw2();
  hRateMu->Sumw2();

  TCanvas* c1 = new TCanvas("c1","Dsum spectrum");
  c1->SetLogy();
  c1->SetGridx();
  sprintf(name,"hEneD");
  sprintf(title,"Run %d, dsum energy, MeV",run);
  hEneD->SetNameTitle(name,title);
  hEneD->GetXaxis()->SetTitle("Energy, MeV");
  hEneD->Draw();

  TCanvas* c2 = new TCanvas("c2","Asum spectrum");
  c2->SetLogy();
  c2->SetGridx();
  sprintf(name,"hEneA");
  sprintf(title,"Run %d, asum energy, MeV",run);
  hEneA->SetNameTitle(name,title);
  hEneA->GetXaxis()->SetTitle("Energy, MeV");
  hEneA->Draw();

  TCanvas* c3 = new TCanvas("c3","Trigger type");
//  c3->SetLogy();
  c3->SetGridx();

  sprintf(name,"hTrg");
  sprintf(title,"Run %d, FADC trigger type",run);
  hTrg->SetNameTitle(name,title);
  hTrg->GetXaxis()->SetTitle("type of the FADC trigger");

/*  TH1F* hTrg = new TH1F("hTrg","FADC trigger type",max_trgtype+2,0,max_trgtype+2);
  hTrg->SetNameTitle(name,title);

  hTrg->GetXaxis()->SetTitle("type of the FADC trigger");
  for (int t = 1; t < max_trgtype+2; t++)  {
       double bin_content = hTrg_temp->GetBinContent(t);
       hTrg->Fill(t-1,bin_content);
      }*/
  hTrg->Draw();

  TCanvas* c4 = new TCanvas("c4","BX evnum vs. FADC evnum");
//  c4->SetGridy();

  TPad* pad1 = new TPad("pad1","",0,0,1,1);
  pad1->SetGridy();
  pad1->Draw();
  pad1->cd();

  sprintf(name,"hEvBx");
  sprintf(title,"Run %d, BX evnum vs. FADC evnum",run);
  hEvBx->SetNameTitle(name,title);
  hEvBx->GetXaxis()->SetTitle("FADC event number");
  hEvBx->GetYaxis()->SetTitle("BX event number");
  hEvBx->Draw();

  pad1->Update();
  TPaveStats* ps1 = (TPaveStats*)hEvBx->GetListOfFunctions()->FindObject("stats");
  ps1->SetX1NDC(0.2); ps1->SetX2NDC(0.4);
//  ps1->SetY1NDC(0.65); ps1->SetY2NDC(0.85);
  pad1->Modified();

  TCanvas* c5 = new TCanvas("c5","evnum vs. unix_time");

  TPad* pad2 = new TPad("pad2","",0,0,1,1);
  pad2->SetGridy();
  pad2->Draw();
  pad2->cd();

  sprintf(name,"hEv");
  sprintf(title,"Run %d, evnum vs. unix_time",run);
  hEv->SetNameTitle(name,title);
  hEv->GetXaxis()->SetTitle("unix time");
  hEv->GetYaxis()->SetTitle("FADC event number");
  hEv->Draw();

  pad2->Update();
  TPaveStats* ps2 = (TPaveStats*)hEv->GetListOfFunctions()->FindObject("stats");
  ps2->SetX1NDC(0.2); ps2->SetX2NDC(0.4);
//  ps2->SetY1NDC(0.65); ps2->SetY2NDC(0.85);
  pad2->Modified();

  TCanvas* c6 = new TCanvas("c6","Counting rate");
  c6->SetGridy();
  hRate->GetXaxis()->SetTitle("Time counts, 1 bin = 30 s");
  hRate->Draw();

  TCanvas* c7 = new TCanvas("c7","Muons counting rate");
  c7->SetGridy();
  hRateMu->GetXaxis()->SetTitle("Time counts, 1 bin = 60 s");
  hRateMu->Draw();

}

void save()  {

  char zero[8];
  printf("Type 0 if you want to validate all events, or N events being validated\n");
  gets(zero);
  int validated_events = atoi(zero);
  if (validated_events == 0)  validated_events = nentries;

  stringstream ss;
  ss.fill('0');
  ss << setw(6) << run;

  string filename = "Run";
  filename += ss.str();
  filename += "_fadc_c16.root";

//  ofstream file;
//  file.open("test.txt",ios::app);

  FILE* fp;
//  char name[80] = "/home/evgeny/fadc/test.txt";
  char name[80] = "/home/production/FADC_run_validation_out.txt";

  if ((fp=fopen(name,"a")) == NULL) {
      printf("Failed opening txt file\n");
      exit(1);
     }

//  file << filename << setw(7) << nentries << setw(4) << n_antinues << endl;
/*
  file << filename << setw(12) << start_time << setw(12) << end_time << setw(7) << runtime << setw(7) << nentries << setw(4) << n_antinues;
  file << endl;
*/

 
fprintf(fp,"%s %12d %12d %7d %7d %7d %9.1f %7.1f %5d %6.1f %7.3f %4d %4d\n",filename.c_str(),start_time,end_time,runtime,nentries,validated_events,n_muons_daily,percent_muons_ps,n_antinues,percent_noise,n_1MeV_Hz,n_errors,n_dcode_errors);

  fclose(fp);

  if (validated_events == nentries)
//      cout << "Run " << run << " validated" << endl;
      printf("Run %d validated.\n",run);
  else printf("Run %d, validated events: %d\n",run,validated_events);

}


