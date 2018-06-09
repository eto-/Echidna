// Run Validation Macro
// MP 12/2/07
// BC, LL 20/2/07
// GT  18/3/07
// GT  01/08/07
// PR  07/09/07
#include <iostream>
#include <vector>
#include "TPaveLabel.h"
#include "TNetFile.h"
#include "TFile.h"
#include "TH1F.h"
#include "TF1.h"
#include "TH2F.h"
#include "TTree.h"
#include "event/BxEvent.hh"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TMath.h"
#include "TROOT.h"
#include <stdlib.h>
#include <TSQLServer.h>
#include <TSQLResult.h>
#include <TTreeResult.h>
#include <TSQLRow.h>
#include <string>
#include <sstream>
#include <iostream>
#include <time.h>

#define muon
#define laben

class bx_Validation_SAVE_Class {
public:
  
  Int_t run, cycle, RunEvents, evnum_start, evnum_end, t_diff, btb_thresh;
  Float_t VerticalMean1, VerticalMean2, VerticalMean3, VerticalMean4, VerticalMean5;  //laben
  Float_t MuonRateID, MuonRateOD, MuonRateOD25; //muon
  Int_t   LastMcrEvent, LastAlignedEvent;  //muon
  Int_t good_laben_channels, good_muon_channels;
  char *loc_validation_time, *loc_valid, *loc_group, *loc_root_files;
  const char *loc_start_time;
  
  TH1F *histo1[4], *histo2[4], *histo3[4], *histo4[4], *histo5[4], *histo6[2], *histo7;
  TH1F *histo8[6], *histo9[6], *histo10[6], *histo11[4], *histo12[3], *histo12a[2], *histo13[6], *histo13a[2], *histo14[6], *histo15[6];
  TCanvas *canvas1,*canvas2,*canvas3,*canvas4,*canvas5,*canvas6,*canvas7;   //laben
  TCanvas *canvas8,*canvas9,*canvas10,*canvas11,*canvas12,*canvas13,*canvas14,*canvas15; //muon
  
private:  
};

bx_Validation_SAVE_Class *save_validation_parameters = new bx_Validation_SAVE_Class;

void SaveValidation();

long RunEvents(int RunNumber) {
  TSQLServer* db = TSQLServer::Connect("pgsql://bxdb.lngs.infn.it/daq_config","borex_guest","xyz");
  TSQLResult* result = db->Query(Form("SELECT \"Events\" FROM \"Run\" WHERE \"RunNumber\" = %d", RunNumber) );
  if( result->GetRowCount() == 1) {
    TSQLRow* row = result->Next(); 
    return atol(row->GetField(0));
  } else return -1;
}

const char * StartDate(int RunNumber) {
  TSQLServer* db = TSQLServer::Connect("pgsql://bxdb.lngs.infn.it/daq_config","borex_guest","xyz");
  TSQLResult* result = db->Query(Form("SELECT \"StartTime\" FROM \"Run\" WHERE \"RunNumber\" = %d", RunNumber) );
  if( result->GetRowCount() == 1) {
    TSQLRow* row = result->Next(); 
    return row->GetField(0);
  }  else {
    const char * bla = 0;
    return bla;
  }  
}
 
char * TimeNow( ) {
        time_t t = time(0);
	static char s[200];
        struct tm *ts = localtime(&t);
        strftime(s,200,"%Y-%m-%d %H:%M:%S",ts);
	//        printf("Tempo: %s\n",s);
        return s;
}


long int StartDate_int(int RunNumber) {
  TSQLServer* db = TSQLServer::Connect("pgsql://bxdb.lngs.infn.it/daq_config","borex_guest","xyz");
  TSQLResult* result = db->Query(Form("SELECT \"StartTime\" FROM \"Run\" WHERE \"RunNumber\" = %d", RunNumber) );
  if( result->GetRowCount() == 1) {
    TSQLRow* row = result->Next(); 
    result = db->Query(Form("SELECT EXTRACT(EPOCH FROM TIMESTAMP \'%s\')", row->GetField(0)));
    if (result->GetRowCount() == 1) {
      row = result->Next();
      return ::strtol (row->GetField(0), 0, 10);
    }
  }
  return 0;
}

std::string GetPath (long int start_time, int cycle) {
  struct tm *start_date = ::localtime (&start_time);
  time_t week_second = start_time - ((start_date->tm_wday * 24 + start_date->tm_hour) * 60 + start_date->tm_min) * 60;
  struct tm *week_date = ::localtime (&week_second);
  char tmp_str[100];
  ::strftime (tmp_str, 99, "/%Y/%b_%d/", week_date);
  std::ostringstream str;
  str << "http://bxmaster-data.lngs.infn.it//bxstorage/rootfiles/cycle_" << cycle << tmp_str;
  return str.str ();
}

std::string GetWeek (long int start_time) {
  struct tm *start_date = ::localtime (&start_time);
  time_t week_second = start_time - ((start_date->tm_wday * 24 + start_date->tm_hour) * 60 + start_date->tm_min) * 60;
  struct tm *week_date = ::localtime (&week_second);
  char tmp_str[100];
  ::strftime (tmp_str, 99, "%b_%d", week_date);
  std::ostringstream str;
  str << tmp_str;
  return str.str ();
}


void SetGeneralStyle() {

  int color1 = 9;               // main color: dark blue
	gStyle->SetAxisColor(color1);    
	gStyle->SetAxisColor(color1,"Y");    
	gStyle->SetLabelColor(color1);  
	gStyle->SetLabelColor(color1,"Y");  
	gStyle->SetTitleColor(color1,"X"); 
	gStyle->SetTitleColor(color1,"Y"); 
	gStyle->SetTitleTextColor(color1); 
	gStyle->SetStatTextColor(color1); 
	gStyle->SetOptStat(111111); 

  // bullet
	gStyle->SetMarkerStyle(7);    // small squared bullet
	gStyle->SetMarkerColor(50);   // red

}


//
// compute mean of vertical axis in one-dim histos 
// through an option, you can include zero bins or not
// points that are out mean more than 5 sigmas are skipped in second iteration
//
float VerticalMean(TH1F* h, float rms_limit, bool IgnoreZeros = true ) {
  int size = h->GetSize() - 2; //number of bins, without under/over flow
  float mean = 0.,rms=0;
  float S = 0.;
  int cnt=0;
  for(int i=1; i<=size; i++) {
    double value= h->GetBinContent (i); 
    if (value >0.0000001 || !IgnoreZeros) {
      cnt ++;
      double delta = value - mean;   
      if(cnt) mean = mean + delta / cnt;
      S = S + delta * (value - mean);
    }
  }
  if(cnt) rms = TMath::Sqrt(double(S / cnt));

  cout << "mean " << mean << " rms " << rms << endl;
  
  
  double mean_old = mean;
  mean = 0.;
  S = 0.;
  cnt=0;
  
  for(int i=1; i<=size; i++) {
    double value= h->GetBinContent (i); 
    //cout << "value " << value << " mean_old" << mean_old << " rms_limit " << rms_limit << " rms " << rms << endl; 
    if ( TMath::Abs(value - mean_old) < rms_limit* rms && (value >0.0000001 || !IgnoreZeros)) {
      cnt ++;
      double delta = value - mean;   
      if(cnt) mean = mean + delta / cnt;
      S = S + delta * (value - mean);
    }
  }	
  if(cnt) rms = ::sqrt(S / cnt);
  std::cout << "mean " << mean << " rms " << rms << " (truncated) " << std::endl;
  return mean;
}



//
// ----------- first canvas : event size vs event number
//
void FirstCanvas(TCanvas * c1, TH1F** h) {
  gStyle->SetOptLogy(0);
  float global[4][4];
  
	// fit these histograms with a line
  double low = 1.01*h[0]->GetBinLowEdge(1);
  double high = 0.99*h[0]->GetBinLowEdge(h[0]->GetSize());
  TF1 *line = new TF1("line","[0] + [1]*x",low,high);
	
  c1->Divide(2,2);
  //h[0]->SetAxisRange(50.,250.,"Y");
  //h[1]->SetAxisRange(0.,200.,"Y");
  //h[2]->SetAxisRange(0.,50.,"Y");
  //h[3]->SetAxisRange(1800.,2300.,"Y");
  for(int i=0; i<4; i++) {
    // reasonable starting point
    float mean = 2000;
    //float mean = VerticalMean(h[i], 3., true);
    line->SetParameter(0,mean);
    line->SetParameter(1,0.);
    c1->cd(i+1);
    h[i]->Fit("line","QR");
    global[i][0] = line->GetParameter(0);
    global[i][1] = line->GetParameter(1);
    global[i][2] = line->GetChisquare() / line->GetNDF();
    
    cout << "Laben: Fit: " << i <<  " const=" << line->GetParameter(0) << " slope=" << line->GetParameter(1) << " chi2/ndf=" << global[i][2] << "\n";
   
    if(h[i]->GetMaximum () > 9999) h[i]->SetLabelOffset(-0.01,"Y");
    h[i]->SetTitleSize(0.045,"X");
    h[i]->SetTitleSize(0.045,"Y");
    h[i]->SetTitleOffset(0.8,"X");
    h[i]->SetTitleOffset(1.1,"Y");

    h[i]->SetXTitle("Event number");
    h[i]->SetYTitle("Event size (N decoded hits)");
    h[i]->Draw("simple");
    //c1->Update ();
  }
  c1->Update();
}

//
// ----------- second canvas : event size vs event number with rough binning
//
void SecondCanvas(TCanvas * c2, TH1F** h) {
  gStyle->SetOptLogy(0);
  float global[4][4];
  
	// fit these histograms with a line
  TF1 *line = new TF1("line","[0] + [1]*x",0.,h[0]->GetSize());

  c2->Divide(2,2);

  for(int i=0; i<4; i++) {
    c2->cd(i+1);
  
    // reasonable starting point
    float mean = TMath::Median(h[i]->GetSize(),h[i]->GetArray()); // VerticalMean(h[i], true);
    line->SetParameter(0,mean);
    line->SetParameter(1,0.);
    
    h[i]->Fit("line","Q");
    //h[i]->Fit("line","QME");
    global[i][0] = line->GetParameter(0);
    global[i][1] = line->GetParameter(1);
    global[i][2] = line->GetChisquare() / line->GetNDF();
    
    cout << "Laben: Fit: " << i <<  " const=" << line->GetParameter(0) << " slope=" << line->GetParameter(1) 
	 << " chi2/ndf=" << global[i][2] << "\n";  

    if(h[i]->GetMaximum () > 9999) h[i]->SetLabelOffset(-0.01,"Y");
    h[i]->SetTitleSize(0.045,"X");
    h[i]->SetTitleSize(0.045,"Y");
    h[i]->SetTitleOffset(0.8,"X");
    h[i]->SetTitleOffset(1.1,"Y");
	
    h[i]->SetXTitle("Run duration in 200 bins");
    h[i]->SetYTitle("Event size (decoded hits)");	
    h[i]->Draw("simple");
  }
  c2->Update();
}

//
// ----------- third canvas : number of neutrino trigger vs time 
//
void ThirdCanvas(TCanvas * c3, TH1F** h, double last_time) {
  
    // fit this histogram with a line
  TF1 *line = new TF1("line","[0] + [1]*x",0.,last_time);
  TF1 *line_100 = new TF1("line_100","[0] + [1]*x",0.,last_time);

  float global[3];
  
  h[0]->SetAxisRange (0,1.1*last_time,"X");
  h[1]->SetAxisRange (0,1.1*last_time,"X");
  c3->Divide(1,2);
  c3->cd (1);
 
  // reasonable starting point 1
  float mean = TMath::Median(h[0]->GetSize(),h[0]->GetArray()); // VerticalMean(h[i], true);
  line->SetParameter(0,mean);
  line->SetParameter(1,0.);
  h[0]->Fit("line","QR");  


   //h->Fit("line","QME");
  global[0] = line->GetParameter(0);
  global[1] = line->GetParameter(1);
  global[2] = line->GetChisquare() / line->GetNDF();
  cout << "Neutrino rate: Fit: const = " << line->GetParameter(0) << " slope = " << line->GetParameter(1) 
       << " chi2/ndf = " << global[2] << "\n";

  h[0]->SetTitleSize(0.045,"X");
  h[0]->SetTitleSize(0.055,"Y");
  h[0]->SetTitleOffset(0.8,"Y");
  h[0]->SetXTitle("Time [seconds], 1 min/bin");
  h[0]->SetYTitle("Trigger rate (counts per second)");
  h[0]->Draw("simple");
 
  c3->cd(2);
  // reasonable starting point 2
  float mean_100 = TMath::Median(h[1]->GetSize(),h[1]->GetArray()); // VerticalMean(h[i], true);
  line_100->SetParameter(0,mean_100);
  line_100->SetParameter(1,0.);
  h[1]->Fit("line_100","QR");  

  h[1]->SetTitleSize(0.045,"X");
  h[1]->SetTitleSize(0.055,"Y");
  h[1]->SetTitleOffset(0.8,"Y");
  h[1]->SetXTitle("Time [seconds], 1 min/bin");
  h[1]->SetYTitle("Trigger rate (counts per second)");

  h[1]->Draw ("simple");
  c3->Update ();
}


void FourthCanvas(TCanvas * c4, TH1F** h,TH1F* h1) {
  int max = h[0]->GetMaximumBin();
  h[0]->SetAxisRange(0,1000.);
  h[1]->SetAxisRange(0,1.4*max);
  h[2]->SetAxisRange(0,1.4*max);
  h[3]->SetAxisRange(0,1.4*max);
  
  for (int i = 0; i < 4; i++) {
    h[i]->SetTitleSize(0.045,"X");
    h[i]->SetTitleSize(0.045,"Y");
    h[i]->SetTitleOffset(0.8,"X");
    h[i]->SetTitleOffset(1.,"Y");
  }

  char lbl[13];
  sprintf (lbl,"Decoded Hits");
  TPaveLabel * label = new TPaveLabel(0.6,0.65,0.7,0.75,lbl,"NDC");
  label->SetBorderSize(0);
  label->SetTextColor(2);
  label->SetTextSize(0.5);
  char lbl2[13];
  sprintf (lbl2,"Clustered Hits");
  TPaveLabel * label2 = new TPaveLabel(0.6,0.55,0.7,0.65,lbl2,"NDC");
  label2->SetBorderSize(0);
  label2->SetTextSize(0.5);

  c4->Divide(2,2);
  c4->cd (1);
  h[0]->SetXTitle("N raw hits");
  h[0]->Draw ();
  c4->cd (2);
  h[1]->Draw ();
  h[1]->SetXTitle("N decoded hits");
  c4->cd (3);
  h1->SetXTitle("N clusters");
  h1->Draw ();
  c4->cd (4);
  h[2]->SetXTitle("N hits");
  h[2]->Draw ();
  h[3]->SetLineColor (2);
  h[3]->Draw("same");
  label->Draw ();
  label2->Draw ();
  c4->Update ();

}

void FifthCanvas(TCanvas* c5, TH1F** h) {
  gStyle->SetOptLogy(1);
  for (int i = 0; i < 4; i++) {
    h[i]->SetTitleSize(0.045,"X");
    h[i]->SetTitleSize(0.045,"Y");
    h[i]->SetTitleOffset(0.8,"X");
    h[i]->SetTitleOffset(1.,"Y");
    h[i]->SetXTitle("Time [ns]");
  }
  c5->Divide(2,2);
  c5->cd (1);
  h[0]->Draw ();
  c5->cd (2);
  h[1]->Draw ();
  c5->cd (3);
  h[2]->Draw ();
  c5->cd (4);
  h[3]->Draw ();
  gStyle->SetOptLogy(0);
  c5->Update ();
}

void SixthCanvas(TCanvas * c6, TH1F** h) {
  gStyle->SetOptLogy(0);
  for (int i = 0; i < 4; i++) {
    h[i]->SetTitleSize(0.045,"X");
    h[i]->SetTitleSize(0.045,"Y");
    h[i]->SetTitleOffset(0.8,"X");
    h[i]->SetTitleOffset(1.,"Y");
    h[i]->SetXTitle("Logical channel");
  }
  c6->Divide(2,2);
  c6->cd (1);
  h[0]->Draw ();
  c6->cd (2);
  h[1]->Draw ();
  c6->cd (3);
  h[2]->Draw ();
  c6->cd (4);
  h[3]->Draw ();
  c6->Update ();
}

void SeventhCanvas(TCanvas * c7, TH1F* raw, TH1F* dec, TH1F* clus, TH1F* charge, double duration, int run) {
  gStyle->SetOptLogy(1);
  c7->Divide(2,2);
  c7->cd (1);
  raw->SetAxisRange(0,1000);
  raw->Draw ();
  c7->cd (2);
  dec->SetAxisRange(0,400);
  dec->Draw ();
  c7->cd (3);
  clus->SetTitleSize(0.045,"X");
  clus->SetTitleOffset(0.8,"X");
  clus->SetXTitle("N clustered hits");
  clus->SetAxisRange(0,400);
  clus->Draw ();
  c7->cd (4);
  charge->SetTitleSize(0.045,"X");
  charge->SetTitleOffset(0.8,"X");
  charge->SetXTitle("Charge");
  charge->SetAxisRange(0,400);
  charge->Draw ();

  //to plot run 5031
  TFile *vzor = new TFile("/home/production/run_plots/6384.root");
  //  TFile *vzor = new TFile("/home/production/run_plots/6472.root");
  //  TFile *vzor = new TFile("/home/production/run_plots/5031.his");
  double duration_norm =  19956; //(6384), BTB=20, 7mus
  //  double duration_norm =  13112; //(6472), BTB=25, 7mus
  //  double duration_norm = 21594; (run 5031)
  TH1F *nraw_5031 = (TH1F*)vzor->Get("n_raw_hits");
  nraw_5031->SetLineStyle(2);
  nraw_5031->SetLineColor(2);
  nraw_5031->SetAxisRange(0,400);
  nraw_5031->Scale(duration/duration_norm);
  TH1F *ndec_5031 = (TH1F*)vzor->Get("n_dec_hits");
  ndec_5031->SetLineStyle(2);
  ndec_5031->SetLineColor(2);
  ndec_5031->SetAxisRange(0,400);
  ndec_5031->Scale(duration/duration_norm);
  TH1F *nclus_5031 = (TH1F*)vzor->Get("n_clus_hits");
  nclus_5031->SetLineStyle(2);
  nclus_5031->SetLineColor(2);
  nclus_5031->SetAxisRange(0,400);
  nclus_5031->Scale(duration/duration_norm);
  TH1F *charge_5031 = (TH1F*)vzor->Get("charge");
  charge_5031->SetLineStyle(2);
  charge_5031->SetLineColor(2);
  charge_5031->SetAxisRange(0,400);
  charge_5031->Scale(duration/duration_norm);
  c7->cd (1);
  nraw_5031->Draw("same");
  char lbl[13];
  sprintf (lbl,"Run 6384, BTB=20");
  //  sprintf (lbl,"Run 6472, BTB=25");
  //  sprintf (lbl,"Run 5031");
  TPaveLabel * label = new TPaveLabel(0.6,0.65,0.7,0.75,lbl,"NDC");
  label->SetBorderSize(0);
  label->SetTextColor(2);
  label->SetTextSize(0.5);
  label->Draw ();
  char lbl2[13];
  sprintf (lbl2,"Run %d",run);
  TPaveLabel * label2 = new TPaveLabel(0.6,0.55,0.7,0.65,lbl2,"NDC");
  label2->SetBorderSize(0);
  label2->SetTextSize(0.5);
  label2->Draw ();
  c7->cd (2);
  ndec_5031->Draw("same");
  label->Draw ();
  label2->Draw ();
  c7->cd (3);
  nclus_5031->Draw("same");
  label->Draw ();
  label2->Draw ();
  c7->cd (4);
  charge_5031->Draw("same");
  //vzor->Close ();
  label->Draw ();
  label2->Draw ();
  c7->Update ();
  gStyle->SetOptLogy(0);
}

//
// -------- Eighth Canvas

void EighthCanvas(TCanvas * c8, TH1F** m, int n_muons_id, int n_muons_od) {
  gStyle->SetOptLogy(0);
  double maxscale_1 = 0;
  for (int i=0; i<3; i++) {
    if (m[i]->GetBinContent(m[i]->GetMaximumBin())>maxscale_1) maxscale_1 = m[i]->GetBinContent(m[i]->GetMaximumBin());
  }
  maxscale_1 *= 1.1;
  double maxscale_2 = 0;
  for (int i=3; i<5; i++) {
    if (m[i]->GetBinContent(m[i]->GetMaximumBin())>maxscale_2) maxscale_2 = m[i]->GetBinContent(m[i]->GetMaximumBin());
  }
  maxscale_2 *= 1.1;

  TH1F *h_regrate = new TH1F("hrr","hrr",10,0,10);
  TH1F *h_1s_low  = new TH1F("h1l","h1l",10,0,10);
  TH1F *h_1s_high = new TH1F("h1h","h1h",10,0,10);
  TH1F *h_2s_low  = new TH1F("h2l","h2l",10,0,10);
  TH1F *h_2s_high = new TH1F("h2h","h2h",10,0,10);
  for (int i=0; i<10; i++) {
    h_regrate->Fill(i,4250);
    h_1s_low ->Fill(i,4250.*(1.-1./sqrt(n_muons_id/10.)));
    h_1s_high->Fill(i,4250.*(1.+1./sqrt(n_muons_id/10.)));
    h_2s_low ->Fill(i,4250.*(1.-2./sqrt(n_muons_id/10.)));
    h_2s_high->Fill(i,4250.*(1.+2./sqrt(n_muons_id/10.)));
  }

  TH1F *h_regrate_od = new TH1F("hrrod","hrr",10,0,10);
  TH1F *h_1s_low_od  = new TH1F("h1lod","h1l",10,0,10);
  TH1F *h_1s_high_od = new TH1F("h1hod","h1h",10,0,10);
  TH1F *h_2s_low_od  = new TH1F("h2lod","h2l",10,0,10);
  TH1F *h_2s_high_od = new TH1F("h2hod","h2h",10,0,10);
  for (int i=0; i<10; i++) {
    h_regrate_od->Fill(i,4000);
    h_1s_low_od ->Fill(i,4000.*(1.-1./sqrt(n_muons_od/10.)));
    h_1s_high_od->Fill(i,4000.*(1.+1./sqrt(n_muons_od/10.)));
    h_2s_low_od ->Fill(i,4000.*(1.-2./sqrt(n_muons_od/10.)));
    h_2s_high_od->Fill(i,4000.*(1.+2./sqrt(n_muons_od/10.)));
  } 

  c8->Divide(2,1);
  c8->cd (1);
  m[0]->SetLineColor(2);
  m[0]->SetAxisRange(0, maxscale_1, "Y");
  m[0]->SetStats(kFALSE);
  m[0]->SetMarkerColor(2);
  m[1]->SetLineColor(3);
  m[1]->SetLineStyle(2);
  m[1]->SetMarkerColor(3);
  m[2]->SetLineStyle(5);
  m[2]->SetMarkerColor(1);
  m[0]->GetXaxis()->SetTitle("time bins");
  m[0]->GetYaxis()->SetTitle("events/day");
  m[0]->Draw("P0,L");
  m[1]->Draw("P0,L,same");
  m[2]->Draw("P0,L,same");
  h_regrate->SetLineColor(4);
  h_regrate->Draw("L,same");
  h_1s_low ->SetLineStyle(2);
  h_1s_low ->SetLineColor(4);
  h_1s_low ->Draw("L,same");
  h_1s_high->SetLineStyle(2);
  h_1s_high->SetLineColor(4);
  h_1s_high->Draw("L,same");
  h_2s_low ->SetLineStyle(2);
  h_2s_low ->SetLineColor(33);
  h_2s_low ->Draw("L,same");
  h_2s_high->SetLineStyle(2);
  h_2s_high->SetLineColor(33);
  h_2s_high->Draw("L,same");
  TLegend *legend1 = new TLegend(.725,.79,.975,.99);
  legend1->SetTextColor(9);
  legend1->SetLineColor(9);
  legend1->AddEntry(m[0],"MTB");
  legend1->AddEntry(m[1],"MCR");
  legend1->AddEntry(m[2],"ID flag");
  legend1->AddEntry(h_regrate,"exp.rate","L");
  legend1->AddEntry(h_1s_high,"1-sig.dev.","L");
  legend1->AddEntry(h_2s_high,"2-sig.dev.","L");
  legend1->SetTextSize(0.04);
  legend1->Draw();
  c8->cd (2);
  m[3]->SetLineColor(2);
  m[3]->SetAxisRange(0, maxscale_2, "Y");
  m[3]->SetStats(kFALSE);
  m[3]->SetMarkerColor(2);
  m[4]->SetLineColor(3);
  m[4]->SetLineStyle(2);
  m[4]->SetMarkerColor(3);
  m[5]->SetLineStyle(5);
  m[5]->SetMarkerColor(1);
  m[3]->GetXaxis()->SetTitle("time bins");
  m[3]->GetYaxis()->SetTitle("events/day"); 
  m[3]->Draw("P0,L");
  m[4]->Draw("P0,L,same");
  m[5]->Draw("P0,L,same");
  h_regrate_od->SetLineColor(4);
  h_regrate_od->Draw("L,same");
  h_1s_low_od ->SetLineStyle(2);
  h_1s_low_od ->SetLineColor(33);
  h_1s_low_od ->Draw("L,same");
  h_1s_high_od->SetLineStyle(2);
  h_1s_high_od->SetLineColor(33);
  h_1s_high_od->Draw("L,same");
  h_2s_low_od ->SetLineStyle(2);
  h_2s_low_od ->SetLineColor(33);
  h_2s_low_od ->Draw("L,same");
  h_2s_high_od->SetLineStyle(2);
  h_2s_high_od->SetLineColor(33);
  h_2s_high_od->Draw("L,same");
  TLegend *legend2 = new TLegend(.7,.79,.975,.99);
  legend2->SetTextColor(9);
  legend2->SetLineColor(9);
  legend2->AddEntry(m[0],"MTB");
  legend2->AddEntry(m[2],"&nhits>25");
  legend2->AddEntry(m[1],"MCR");
  legend2->AddEntry(h_regrate_od,"exp.rate","L");
  legend2->AddEntry(h_1s_high_od,"1-sig.dev.","L");
  legend2->AddEntry(h_2s_high_od,"2-sig.dev.","L");
  legend2->SetTextSize(0.04);
  legend2->Draw();
  c8->Update ();
}

//
// --------- ninth Canvas
//

void NinthCanvas(TCanvas * c9, TH1F** m) {

  float global[6][4];
  gStyle->SetOptLogy(0);

  // fit these histograms with a line
  double low = 1.01*m[0]->GetBinLowEdge(1);
  double high = 0.99*m[0]->GetBinLowEdge(m[0]->GetSize());
  TF1 *line = new TF1("line","[0] + [1]*x",low,high);

  c9->Divide(3,2);

  for(int i=0; i<6; i++) {
    if (m[i]->GetEntries()==0) {
      cout << "Muon:  event size histogram #" << i << " is empty!" << endl;
      continue;
    }
  // reasonable starting point
    float mean = 200;
    line->SetParameter(0,mean);
    line->SetParameter(1,0.);
    c9->cd(i+1);
    if (i==5) line->SetRange(1250.,high);
    m[i]->Fit("line","QR");
    global[i][0] = line->GetParameter(0);
    global[i][1] = line->GetParameter(1);
    global[i][2] = line->GetChisquare() / line->GetNDF();
    cout << "Muon:  Fit: " << i <<  " const=" << line->GetParameter(0) << " slope=" << line->GetParameter(1) << " chi2/ndf=" << global[i][2] << "\n";
    m[i]->GetYaxis()->SetTitle("event size");
    m[i]->GetXaxis()->SetTitle("event number");
    m[i]->Draw("simple");
  }
  c9->Update();
}

//
// --------- tenth canvas
//

void TenthCanvas(TCanvas * c10, TH1F** m, double evnum_end) {

  float global[6][4];
  gStyle->SetOptLogy(0);

  // fit these histograms with a line
  TF1 *line = new TF1("line","[0] + [1]*x",0,m[0]->GetSize());

  c10->Divide(3,2);

  for(int i=0; i<6; i++) {
    if (m[i]->GetEntries()==0) {
      cout << "Muon:  event size histogram #" << i << " is empty!" << endl;
      continue;
    }
     c10->cd(i+1);
    // reasonable starting point
    float mean = TMath::Median(m[i]->GetSize(),m[i]->GetArray()); // VerticalMean(h[i], true);
    line->SetParameter(0,mean);
    line->SetParameter(1,0.);
    if (i==5) line->SetRange( (double) (m[5]->GetSize()*1250.)/evnum_end+2, m[5]->GetSize() );
    m[i]->Fit("line","QR");
    //h[i]->Fit("line","QME");
    global[i][0] = line->GetParameter(0);
    global[i][1] = line->GetParameter(1);
    global[i][2] = line->GetChisquare() / line->GetNDF();
    cout << "Muon:  Fit: " << i <<  " const=" << line->GetParameter(0) << " slope=" << line->GetParameter(1)
         << " chi2/ndf=" << global[i][2] << "\n";
    m[i]->GetYaxis()->SetTitle("event size");
    m[i]->GetXaxis()->SetTitle("event number (200 bins)");
    m[i]->Draw("simple");
  }
  c10->Update();
}

//
// --------- eleventh canvas
//

void EleventhCanvas(TCanvas * c11, TH1F** m, double last_time) {

  gStyle->SetOptLogy(0);

  // fit this histogram with a line
  TF1 *line = new TF1("line","[0] + [1]*x",0.,last_time);
  TF1 *line_25 = new TF1("line_25","[0] + [1]*x",0.,last_time);
  float global[6];

  for (int i=0; i<4; i++) m[i]->SetAxisRange (0,1.1*last_time,"X");
  c11->Divide(2,2);

  for (int i=0; i<4; i++) {
    m[i]->GetXaxis()->SetTitle("time [10 minute bins]");
    m[i]->GetYaxis()->SetTitle("event rate [/day]");
  }  

  c11->cd (1);
  // reasonable starting point 1
  float mean = TMath::Median(m[0]->GetSize(),m[0]->GetArray()); // VerticalMean(h[i], true);
  line->SetParameter(0,mean);
  line->SetParameter(1,0.);
  m[0]->Fit("line","QR");
  global[0] = line->GetParameter(0);
  global[1] = line->GetParameter(1);
  global[2] = line->GetChisquare() / line->GetNDF();
  cout << "Muon tt1 rate: Fit: const = " << line->GetParameter(0) << " slope = " << line->GetParameter(1)
       << " chi2/ndf = " << global[2] << "\n";
  m[0]->Draw("simple");
 
  c11->cd(2);
  mean = TMath::Median(m[1]->GetSize(),m[1]->GetArray()); // VerticalMean(h[i], true);
  line->SetParameter(0,mean);
  line->SetParameter(1,0.);
  m[1]->Fit("line","QR");
  m[1]->Draw("simple");
  if (m[2]->GetEntries()>0) {
    c11->cd(3);
  // reasonable starting point 2
    float mean_25 = TMath::Median(m[2]->GetSize(),m[2]->GetArray()); // VerticalMean(h[i], true);
    line_25->SetParameter(0,mean_25);
    line_25->SetParameter(1,0.);
    m[2]->Fit("line_25","QR");
    m[2]->Draw ("simple");
  
    c11->cd(4);
    mean_25 =  TMath::Median(m[3]->GetSize(),m[3]->GetArray()); // VerticalMean(h[i], true);
    line_25->SetParameter(0,mean_25);
    line_25->SetParameter(1,0.);
    m[3]->Fit("line_25","QR");
    global[3] = line->GetParameter(0);
    global[4] = line->GetParameter(1);
    global[5] = line->GetChisquare() / line->GetNDF();
    cout << "Muon tt2 rate (>25): Fit: const = " << line->GetParameter(0) << " slope = " << line->GetParameter(1)
         << " chi2/ndf = " << global[5] << "\n";
    m[3]->Draw ("simple");
  }
  c11->Update ();
}

//
// ---------- Twelfth Canvas
//

void TwelfthCanvas(TCanvas * c12, TH1F** m, TH1F** mc) {
  
  char lbl[13];
  sprintf (lbl,"tt2");
  TPaveLabel * label = new TPaveLabel(0.6,0.65,0.7,0.75,lbl,"NDC");
  label->SetBorderSize(0);
  label->SetTextColor(2);
  label->SetTextSize(0.5);
  char lbl2[13];
  sprintf (lbl2,"tt1");
  TPaveLabel * label2 = new TPaveLabel(0.6,0.55,0.7,0.65,lbl2,"NDC");
  label2->SetBorderSize(0);
  label2->SetTextSize(0.5);
//  m[0]->SetAxisRange(0,200.);
//  m[1]->SetAxisRange(0,200.);
//  m[2]->SetAxisRange(0,200.);
//  mc[1]->SetAxisRange(0,15.);
    
  gStyle->SetOptLogy(1);  
  c12->Divide(2,2);
  c12->cd (1);
  m[0]->GetXaxis()->SetTitle("n_raw_hits");
  m[0]->Draw ();
  c12->cd (2);
  m[1]->GetXaxis()->SetTitle("n_dec_hits");
  m[1]->Draw ();
  c12->cd (3);
  m[2]->GetXaxis()->SetTitle("n_clu_hits");
  m[2]->Draw ();
  c12->cd (4);
  mc[1]->SetLineColor(2);
  mc[1]->GetXaxis()->SetTitle("n_clusters");
  mc[0]->Draw ();
  if (mc[1]->GetEntries()>0) { mc[1]->Draw(); mc[0]->Draw("same"); label->Draw (); }
  else mc[0]->Draw();
  label2->Draw ();
  c12->Update ();
}

//
// ---------- Thirteenth Canvas
//

void ThirteenthCanvas(TCanvas * c13, TH1F** mo, TH1F** mc) {
  
   gStyle->SetOptLogy(1);

  char lbl[13];
  sprintf (lbl,"tt2");
  TPaveLabel * label = new TPaveLabel(0.6,0.65,0.7,0.75,lbl,"NDC");
  label->SetBorderSize(0);
  label->SetTextColor(2);
  label->SetTextSize(0.5);
  char lbl2[13];
  sprintf (lbl2,"tt1");
  TPaveLabel * label2 = new TPaveLabel(0.6,0.55,0.7,0.65,lbl2,"NDC");
  label2->SetBorderSize(0);
  label2->SetTextSize(0.5);
  mo[1]->SetAxisRange(0,200.);
  mo[3]->SetAxisRange(0,200.);
  mo[5]->SetAxisRange(0,200.);
  mc[1]->SetAxisRange(0,15.);
 
  gStyle->SetOptLogy(1);
  c13->Divide(2,2);
  c13->cd (1);
  if (mo[1]->GetEntries()>0) {
    mo[1]->SetLineColor(2);
    mo[1]->GetXaxis()->SetTitle("n_raw_hits");
    mo[1]->Draw();
    mo[0]->Draw("same");
    label->Draw();
    label2->Draw();
    c13->cd (2);
    mo[3]->SetLineColor(2);
    mo[3]->GetXaxis()->SetTitle("n_dec_hits");
    mo[3]->Draw();
    mo[2]->Draw("same");
    label->Draw();
    label2->Draw();
    c13->cd (3);
    mo[5]->SetLineColor(2);
    mo[5]->GetXaxis()->SetTitle("n_clu_hits");
    mo[5]->Draw();
    mo[4]->Draw("same");
    label->Draw();
    label2->Draw();
    c13->cd (4);
    mc[1]->SetLineColor(2);
    mc[1]->Draw ();
    mc[0]->Draw("same");
    label->Draw();
    label2->Draw();
  }
  else {
    mo[0]->GetXaxis()->SetTitle("n_raw_hits");
    mo[0]->Draw();
    label2->Draw();
    c13->cd (2);
    mo[2]->GetXaxis()->SetTitle("n_dec_hits");
    mo[2]->Draw();
    label2->Draw();
    c13->cd (3);
    mo[4]->GetXaxis()->SetTitle("n_clu_hits");
    mo[4]->Draw();
    label2->Draw();
    c13->cd (4);
    mc[0]->Draw ();
    label2->Draw();
  }
  c13->Update();
}

//
// --------- Fourteenth Canvas
//

void FourteenthCanvas(TCanvas* c14, TH1F** m) {
  gStyle->SetOptLogy(1);
  c14->Divide(2,3);
  for (int i=0; i<6; i++) {
    if (m[i]->GetEntries()==0) continue;
    c14->cd (i+1);
    m[i]->GetXaxis()->SetTitle("time [ns]");
    m[i]->Draw ();
  }
  c14->Update ();
}

//
// --------- Fifteenth Canvas
//

void FifteenthCanvas(TCanvas * c15, TH1F** m) {
  gStyle->SetOptLogy();
  c15->Divide(2,3);
  for (int i=0; i<6; i++) {
    if (m[i]->GetEntries()==0) continue;
    c15->cd(i+1);
    m[i]->GetXaxis()->SetTitle("mch");
    m[i]->GetYaxis()->SetTitle("# of decoded_hits");
    m[i]->Draw();
  }
  c15->Update();
  gStyle->SetOptLogy(0);
}


//
// --------- main macro
// 

void RunValidation(const char* fName,  int start=1, int calib = 1, int production = 1) {
  //gROOT->Reset();
  
  SetGeneralStyle();

  //open file
  TFile *f = TFile::Open(fName);
  if (f->IsZombie()) {
    cout << "File " << fName << " not found\n";
    return;
  }
  
  //find the cycle number from the root file name
  int cycle;
  int run;
  int length = strlen(fName);
  fName += (length - 18);
  printf("File name %s \n ",fName);
  if (sscanf (fName, "Run0%5d_c%d.root", &run, &cycle) != 2) { printf ("wrong file name\n"); exit (1); }
  else{
    printf("Cycle is %2d \n", cycle);
  }


  long int start_date = StartDate_int (run);
  
  char * group = new char[20];
  //  char group[20];
  std::string week = GetWeek (start_date);
  sprintf (group, "%s", week.c_str() );
  printf("Group %s \n", group);    

  
    //check electronics calibration
  TH1F *trigref = new TH1F;
  TH1F *lasref = new TH1F;
  TH2F *pulser_charge_vs_lg = new TH2F;
  TH2F *laser_charge_vs_lg = new TH2F;
  TH1F *raw_lg_pulser = new TH1F;
  TH1F *raw_lg_neutrino = new TH1F;

  if(calib){
    char el_calib_file[1000];
    if (run < 9402) sprintf(el_calib_file,"http://bxmaster-data.lngs.infn.it//bxstorage/rootfiles/cycle_11/2008/%s/ancillary/Run0%5d_electronics_calibrations_c11.root",week.c_str(),run);
    if (run > 9401) sprintf(el_calib_file,"http://bxmaster-data.lngs.infn.it//bxstorage/rootfiles/cycle_11/2009/%s/ancillary/Run0%5d_electronics_calibrations_c11.root",week.c_str(),run);
    cout << "File " << el_calib_file;

    TFile *fcalib = TFile::Open(el_calib_file);
    if (fcalib->IsZombie()) {
      cout << " not found\n";
      return;
    } else {
      cout << " found\n";
    }
    
    trigref = (TH1F*)fcalib->Get("barn/bx_calib_laben_electronics/trigref");
    lasref = (TH1F*)fcalib->Get("barn/bx_calib_laben_electronics/lasref");
    pulser_charge_vs_lg = (TH2F*)fcalib->Get("barn/bx_calib_laben_electronics/pulser_charge_vs_lg");
    laser_charge_vs_lg = (TH2F*)fcalib->Get("barn/bx_calib_laben_electronics/laser_charge_vs_lg");
    raw_lg_pulser = (TH1F*)fcalib->Get("barn/bx_calib_laben_electronics/raw_Nhits_pulser_vs_lg");
    raw_lg_neutrino = (TH1F*)fcalib->Get("barn/bx_calib_laben_electronics/raw_Nhits_neutrino_vs_lg");

    raw_lg_pulser->SetTitleSize(0.045,"X");
    raw_lg_neutrino->SetTitleSize(0.045,"X");
    raw_lg_pulser->SetXTitle("Logical channel");
    raw_lg_neutrino->SetXTitle("Logical channel");
    raw_lg_pulser->SetXTitle("Logical channel");
    raw_lg_neutrino->SetXTitle("Logical channel");

    //raw hits non normalizzati vanno visualizzati qua solo nel caso production = 0
    //altrimenti histo si normalizza e plotta later (dopo contare Nevents)
    if (!production) {
      char raw_lg_title[50];
      sprintf(raw_lg_title,"Run0%5d lg distribution of raw hits",run);
      TCanvas *raw_lg = new TCanvas("raw_lg",raw_lg_title);
      raw_lg->Divide(1,2);
      raw_lg->cd(1);
      raw_lg_pulser->Draw ();
      raw_lg->cd(2);
      raw_lg_neutrino->Draw (); 
      raw_lg->Modified();
      raw_lg->Update ();
    }

    char trigger_title[50];
    sprintf(trigger_title,"Run0%5d Charge trg reference, pulser triggers",run);
    TCanvas *trigger = new TCanvas("trigger",trigger_title);
    trigger->Divide (4,4);    
    TH1D *trg[16];
    TPaveLabel *label[16];
    for (int i = 1; i < 17; i ++){
      int trgref_lg = (int) trigref->GetBinContent(i);
      if (trgref_lg) {
	char name[13];
	sprintf(name,"trgref %d",trgref_lg);
	trg[i] = pulser_charge_vs_lg->ProjectionY(name,trgref_lg,trgref_lg);
	trigger->cd(i);
	trg[i]->SetAxisRange(0,100);
	trg[i]->Draw ();
	label[i] = new TPaveLabel(0.1,0.6,0.4,0.95,name,"NDC");
	label[i]->SetBorderSize(0);
	label[i]->Draw (); trigger->Modified();
	//	trigger->Update ();
      }
    }
    trigger->Update();

    char laser_title[50];
    sprintf(laser_title,"Run0%5d Charge laser ref, laser trigger",run);
    TCanvas * laser = new TCanvas("laser",laser_title);
    laser->Divide (3,3);    
    TH1D *las[16];
    TPaveLabel *label_las[16];
    for (int i = 1; i < 17; i ++){
      int lasref_lg = (int) lasref->GetBinContent(i);
      if (lasref_lg) {
	char name[13];
	sprintf(name,"lasref %d",lasref_lg);
	las[i] = laser_charge_vs_lg->ProjectionY(name,lasref_lg,lasref_lg);
	laser->cd(i);
	las[i]->SetAxisRange(0,100);
	las[i]->Draw ();
	label_las[i] = new TPaveLabel(0.1,0.6,0.4,0.95,name,"NDC");
	label_las[i]->SetBorderSize(0);
	label_las[i]->Draw ();laser->Modified();
//	laser->Update ();
      }
    }
  laser->Update();
  }   

  if (production) {

    f->cd ();
    //back to the production root file

    #ifdef laben
      int good_laben_channels = 0;
    #endif
    int good_muon_channels = 0;
    
    
    TTree *tree = (TTree*)f->Get("bxtree");
    if (!tree) {
      cout << "Invalid Borexino root file" << fName << endl;
      return;
    }
    
    
    
    
    TBranch *b = tree->GetBranch("events");
    BxEvent *ev = new BxEvent();
    b->SetAddress(&ev);
    
    tree->SetBranchStatus("*",0);
    tree->SetBranchStatus("run",1);
    tree->SetBranchStatus("evnum",1);
    tree->SetBranchStatus("trigger.trgtype",1);
    tree->SetBranchStatus("trigger.gpstimes[2]",1);
    tree->SetBranchStatus("trigger.btb_threshold",1);
    tree->SetBranchStatus("trigger.btb_inputs",1);
    tree->SetBranchStatus("laben.trigger_time",1);
    tree->SetBranchStatus("laben.laser_time",1);
    tree->SetBranchStatus("laben.n_decoded_hits",1);
    tree->SetBranchStatus("laben.n_raw_hits",1);
    tree->SetBranchStatus("laben.n_clusters",1);
    tree->SetBranchStatus("laben.n_clustered_hits",1);
    tree->SetBranchStatus("laben.decoded_hits.raw_time",1);
    tree->SetBranchStatus("laben.decoded_hits.lg",1);
    tree->SetBranchStatus("laben.clusters.charge",1);
    tree->SetBranchStatus("laben.clusters.mean_time",1);
    tree->SetBranchStatus("laben.clusters.peak_times",1);
    tree->SetBranchStatus("laben.rec_clusters.gatti",1);
    tree->SetBranchStatus("laben.n_live_pmts",1);
    tree->SetBranchStatus("muon.npmts",1);
    // added for muon validation
    #ifdef muon
      tree->SetBranchStatus("muon.n_decoded_hits",1);
      tree->SetBranchStatus("muon.n_raw_hits",1);
      tree->SetBranchStatus("muon.charge_sss",1);
      tree->SetBranchStatus("muon.charge_floor",1);
      tree->SetBranchStatus("muon.decoded_hits.charge",1);
      tree->SetBranchStatus("muon.decoded_hits.time",1);
      tree->SetBranchStatus("muon.decoded_hits.mch",1);
      tree->SetBranchStatus("muon.clusters.charge",1);
      tree->SetBranchStatus("muon.n_clusters",1);
      tree->SetBranchStatus("muon.n_clustered_hits_floor",1);
      tree->SetBranchStatus("muon.n_clustered_hits_sss",1);
      tree->SetBranchStatus("enabled_crates",1);
    #endif
    int entries = tree->GetEntries();

    //info about the run and the  first and last event
    tree->GetEntry(1);
    int evnum_first_entry = ev->GetEvNum ();
    double first_event_time = ev->GetTrigger().GetGpsTimeSec();
    
    tree->GetEntry(entries - 1);
    int evnum_last_entry = ev->GetEvNum ();
    //  double last_event_time = ev->GetTrigger().GetGpsTimeSec();
    
    
    if (run != ev->GetRun ())  std::cout << "WARNING: run number mismatch between file name and root event" << std::endl;
    else {  std::cout << "Opened run has the run number " << run << std::endl;}
    
    //get the correct number of events 
    int end;
    
    cout << "IF YOU KNOW FOR SURE: insert the root file entry  number of the last event you want to process" << endl;
    cout << " OR " << endl;
    cout << "if you expect me to find evnum of the last event written in raw data, hit 0" << endl;
  
    cin >> end ;
    if (end != 0){
      cout << "You want to check " << end << " entries in the root file" << endl;
    }
    if (end == 0){
      end = (int) RunEvents(run);
      cout << "The raw file contains " << end  << " events and you want to check them all" << endl;
    }
    
    
    if (entries < end ){ // in root files less events than the number user wants to analyse
      if (abs(entries - end) > 100){ // a big difference
	cout << "*************************************************************************************************" << endl;
	cout << "ERROR: big difference in the number of events in root file and in raw data" << endl;
	cout << "all the events from the root file will be processed, N entries is : " << entries << endl;
	cout << "last entry event number is " << evnum_last_entry << endl;
	cout << "*************************************************************************************************" << endl;
      }
      end = entries; //silently, just a small discrepancy
    }
    else if( (entries - end) > 100){ // entries > end: in root files is more  events than the number user wants to analyse
      cout << "*************************************************************************************************" << endl;
      cout << "WARNING: the number of events in file is " << entries << " and you are going to analyze only " << end << endl;
      cout << "*************************************************************************************************" << endl;
    }
    
    
    cout << "Info: the number of entries in the root file " << fName << " is " << entries << endl;
    
    
    tree->GetEntry(start);
    int evnum_start = ev->GetEvNum ();
    tree->GetEntry(end - 1);
    int evnum_end = ev->GetEvNum ();
    
    const char * start_time = StartDate(run);
    //char * start_time = (char) StartDate(run);
    char * validation_time = new char;
    validation_time = TimeNow (); 
    printf("Start time of this run: %s\n",start_time);
    printf("Validation time of this run: %s\n",validation_time);
    
    #ifdef laben

    // first canvas : event size stability (skipping precalib events, 1 bin per event)
    int histo1_size = end-start;
    char name1[100];
    sprintf (name1, "RUN_%d_Size vs Event Number - Trigger Laser",run); 
    TH1F *h1[4];
    h1[0]=new TH1F("h1_0",name1,histo1_size,(double)evnum_start,(double)evnum_end);
    h1[1]=new TH1F("h1_1","Size vs Event Number - Trigger Neutrino",histo1_size,(double)evnum_start,(double)evnum_end);
    h1[2]=new TH1F("h1_2","Size vs Event Number - Trigger Random",histo1_size,(double)evnum_start,(double)evnum_end);
    h1[3]=new TH1F("h1_3","Size vs Event Number - Trigger Pulser",histo1_size,(double)evnum_start,(double)evnum_end);
     	
    h1[0]->SetDirectory(0); h1[1]->SetDirectory(0); h1[2]->SetDirectory(0); h1[3]->SetDirectory(0);


    // second canvas: same as before averaging over periods
    int histo2_size = 200;
    TH1F *h2[4];
    char name2[100];
    sprintf (name2, "RUN_%d_Size vs Event Number - Trigger Laser",run); 
    h2[0]=new TH1F("h2_0",name2,histo2_size,0.,(double)histo2_size);
    h2[1]=new TH1F("h2_1","Size vs Event Number - Trigger Neutrino",histo2_size,0.,(double)histo2_size);
    h2[2]=new TH1F("h2_2","Size vs Event Number - Trigger Random",histo2_size,0.,(double)histo2_size);
    h2[3]=new TH1F("h2_3","Size vs Event Number - Trigger Pulser",histo2_size,0.,(double)histo2_size);

    h2[0]->SetDirectory(0); h2[1]->SetDirectory(0); h2[2]->SetDirectory(0); h2[3]->SetDirectory(0);   

  
    //3rd canvas
    double max_time = 3600. * 10;
    int    time_bin_rate = 60; //in seconds
    char name3[100];
    sprintf (name3, "RUN_%d_trigger rate (Hz) vs time in seconds (1 min / bin)",run); 
    char name3b[100];
    sprintf (name3b, "RUN_%d_trigger rate (Hz) for ev nhit > 100 vs time in seconds (1 min / bin)",run); 
    TH1F* neutrino_rate[2];
    neutrino_rate[0] = new TH1F("neutrino_rate all  ",name3, (int) max_time/ time_bin_rate,0,max_time); 
    neutrino_rate[1] = new TH1F("neutrino_rate nhit > 100",name3, (int) max_time/ time_bin_rate,0,max_time); 
    neutrino_rate[0]->Sumw2 ();
    neutrino_rate[1]->Sumw2 ();
    neutrino_rate[0]->SetDirectory(0); neutrino_rate[1]->SetDirectory(0);
  
    //lg
    TH1F* hlg[4];
    char namelg[100];
    sprintf (namelg, "RUN_%d_LASER: decoded hits logical channel",run); 
    hlg[0] = new TH1F("laser_lg",namelg,2240,1,2241);
    hlg[1] = new TH1F("neutrino_lg","NEUTRINO: decoded hits logical channel",2240,1,2241);
    hlg[2] = new TH1F("random_lg","RANDOM: decoded hits logical channel",2240,1,2241);
    hlg[3] = new TH1F("pulser_lg","PULSER: decoded hits logical channel",2240,1,2241);
    TH1F* rawlg;
    rawlg = new TH1F("rawlg","NEUTRINO: raw hits lg",2240,1,2241);
    
    hlg[0]->SetDirectory(0); hlg[1]->SetDirectory(0); hlg[2]->SetDirectory(0); hlg[3]->SetDirectory(0);

    // 4th canvas
    int hits_bins = 6000;
    char name4[100];
    TH1F* nhits[4];
    sprintf (name4, "RUN_%d_N raw hits in neutrino triggers",run); 
    nhits[0] = new TH1F("n_raw_hits", name4, hits_bins, 0, hits_bins);
    nhits[1] = new TH1F("n_dec_hits","N decoded hits in neutrino triggers",hits_bins,0,hits_bins);
    nhits[2] = new TH1F("n_clus_hits","N clustered hits, n_clusters = 1",hits_bins,0,hits_bins);
    nhits[3] = new TH1F("n_dec_hits_cl1","N decoded hits, n_clusters = 1",hits_bins,0,hits_bins);
    TH1F* n_clusters = new TH1F("n_clusters","clusters in neutrino triggers",4,0,4);	
    TH1F* charge = new TH1F("charge","charge for n_clusters == 1",hits_bins,0,hits_bins);

    charge->SetDirectory(0);

    nhits[0]->SetDirectory(0); nhits[1]->SetDirectory(0); nhits[2]->SetDirectory(0); nhits[3]->SetDirectory(0);

    // 5th canvas
    TH1F *htimes[4];
    int gate_length = 17000;
    int time_bin = 5; //ns
    char name5[100];
    sprintf (name5, "RUN_%d_LASER:decoded_hit_time - laser_time",run); 
    htimes[0] = new TH1F("times_laser",name5, 9000/time_bin,-1000,8000); 
    htimes[1] = new TH1F("times_neutrino","NEUTRINO:decoded_hit_time - trigger_time",(gate_length+1000)/time_bin,-1 * gate_length,1000); 
    htimes[2] = new TH1F("times_random","RANDOM:decoded_hit_time - trigger_time",(gate_length+1000)/time_bin,-1 * gate_length,1000); 
    htimes[3] = new TH1F("times_pulser","PULSER:decoded_hit_time - trigger_time",(gate_length+1000)/time_bin,-1 * gate_length,1000); 
  
    htimes[0]->SetDirectory(0); htimes[1]->SetDirectory(0); htimes[2]->SetDirectory(0); htimes[3]->SetDirectory(0);

    #endif

    int neutrino_cnt = 0;
    int ID_muon_cnt = 0;
    int OD_muon_cnt = 0;
    int OD25_muon_cnt = 0;
    int laser_cnt=0;
    int random_cnt=0;
    int pulser_cnt=0;
    
    int t_start=0;
    int t_end=0;
   
    //info about the last analyzed event
    tree->GetEntry(end - 1);
    double tot_time_minutes =  (ev->GetTrigger().GetGpsTimeSec() - first_event_time) / 60.;
    #ifdef laben
    int nbins_neutrino_rate = (int) tot_time_minutes;
    #endif
    #ifdef muon
    int nbins_muon_rate = (int) tot_time_minutes/10;
    #endif
 
    #ifdef muon
    // eighth canvas
    char name_m8[100];
    sprintf (name_m8, "RUN_%d_muon_rate in tt1 : muon flag stability",run); 
    char name_m8a[100];
    sprintf (name_m8a, "RUN_%d_muon_rate in tt2 : muon flag stability",run); 
    int nbins_muon_flags = 10;
    TH1F* m_flag_rate[6];
    m_flag_rate[0] = new TH1F("tt1 muon rate : MTB", name_m8, nbins_muon_flags, 0, nbins_muon_flags); 
    m_flag_rate[1] = new TH1F("tt1 muon rate : MCR", "Muon rate in tt1 : MCR", nbins_muon_flags, 0, nbins_muon_flags); 
    m_flag_rate[2] = new TH1F("tt1 muon rate : IDF", "Muon rate in tt1 : ID flags", nbins_muon_flags, 0, nbins_muon_flags); 
    m_flag_rate[3] = new TH1F("OD muon rate : MTB", name_m8a, nbins_muon_flags, 0, nbins_muon_flags);
    m_flag_rate[4] = new TH1F("OD muon rate : MCR", "Muon rate in tt2 : MCR", nbins_muon_flags, 0, nbins_muon_flags);
    m_flag_rate[5] = new TH1F("OD muon rate : MTB25", "Muon rate in tt2 : MTB > 25 nhits", nbins_muon_flags, 0, nbins_muon_flags); 
    for (int i=0; i<6; i++) m_flag_rate[i]->SetDirectory(0);

    // ninth canvas: muon event size stability
    char name_m1[100];
    sprintf (name_m1, "RUN_%d_LASER",run);
    TH1F *m1[6];
    int m1_size = end-start;
    m1[0]=new TH1F("m1_0", name_m1, m1_size,(double)evnum_start,(double)evnum_end);
    m1[1]=new TH1F("m1_1","MUON tt1", m1_size,(double)evnum_start,(double)evnum_end);
    m1[2]=new TH1F("m1_2","MUON tt2", m1_size,(double)evnum_start,(double)evnum_end);
    m1[3]=new TH1F("m1_3","NEUTRINO", m1_size,(double)evnum_start,(double)evnum_end);
    m1[4]=new TH1F("m1_4","RANDOM", m1_size,(double)evnum_start,(double)evnum_end);
    m1[5]=new TH1F("m1_5","PULSER", m1_size,(double)evnum_start,(double)evnum_end);
   
    m1[0]->SetDirectory(0); m1[1]->SetDirectory(0); m1[2]->SetDirectory(0); m1[3]->SetDirectory(0);
    m1[4]->SetDirectory(0); m1[5]->SetDirectory(0);
   
    // tenth canvas: muon event size stability, 200 bins
    TH1F *m2[6];
    char name_m2[100];
    int m2_size = 200;
    sprintf (name_m2, "RUN_%d_LASER",run);
    m2[0]=new TH1F("m2_0",name_m2, m2_size,0.,(double) m2_size);
    m2[1]=new TH1F("m2_1","MUON tt1", m2_size,0.,(double) m2_size);
    m2[2]=new TH1F("m2_2","MUON tt2", m2_size,0.,(double) m2_size);
    m2[3]=new TH1F("m2_3","NEUTRINO", m2_size,0.,(double) m2_size);
    m2[4]=new TH1F("m2_4","RANDOM", m2_size,0.,(double) m2_size);
    m2[5]=new TH1F("m2_5","PULSER", m2_size,0., (double) m2_size);

    for (int i=0; i<6; i++) m2[i]->SetDirectory(0);

    // 11th canvas: muon trigger rate
    char name_11[100];
    TH1F *muon_rate[4];
    sprintf (name_11, "RiUN_%d_Muon tt1 rate (/day) vs time, total",run);
    double max_muon_time = 3600.*10.;
    int muon_time_bin_rate = 600;
    muon_rate[0] = new TH1F("muon_rate tt1 all",name_11, (int) max_muon_time/ muon_time_bin_rate,0,max_muon_time);
    muon_rate[1] = new TH1F("muon_rate tt1 nhits>25", "Muon tt1 rate (/day) vs time, OD NHits>25", (int) max_muon_time/ muon_time_bin_rate, 0, max_muon_time);
    muon_rate[2] = new TH1F("muon_rate tt2 all", "Muon tt2 rate (/day) vs time, total", (int) max_muon_time/ muon_time_bin_rate,0,max_muon_time);
    muon_rate[3] = new TH1F("muon_rate tt2 nhits>25", "Muon tt2 rate (/day) vs time, OD NHits>25", (int) max_muon_time/ muon_time_bin_rate,0,max_muon_time);
    for (int i=0; i<4; i++) {
      muon_rate[i]->Sumw2();
      muon_rate[i]->SetDirectory(0);
    }
    
    // 15th canvas: mch distributions
    TH1F* hmch[6];
    char namemch[100];
    sprintf (namemch, "RUN_%d_LASER",run);
    hmch[0] = new TH1F("laser_mch",namemch,255,1,256);
    hmch[1] = new TH1F("tt1_muon_mch","MUON tt1",255,1,256);
    hmch[2] = new TH1F("tt2_muon_mch","MUON tt2",255,1,256);
    hmch[3] = new TH1F("neutrino_mch","NEUTRINO",255,1,256);
    hmch[4] = new TH1F("random_mch","RANDOM",255,1,256);
    hmch[5] = new TH1F("pulser_mch","PULSER",255,1,256);
    hmch[0]->SetDirectory(0); hmch[1]->SetDirectory(0); hmch[2]->SetDirectory(0);
    hmch[3]->SetDirectory(0); hmch[4]->SetDirectory(0); hmch[5]->SetDirectory(0);
    
    //12th canvas ID nhits distributions
    TH1F *m_idhits[3];
    int hits_bins_idm = 200;
    char name_m12[100];
    sprintf (name_m12, "RUN_%d_N laben raw hits for tt1 muons",run);
    m_idhits[0] = new TH1F("m_id_n_raw_hits", name_m12, hits_bins_idm, 0, hits_bins_idm*250);
    m_idhits[1] = new TH1F("m_id_n_dec_hits","N laben decoded hits for tt1 muons",hits_bins_idm,0,hits_bins_idm*250);
    m_idhits[2] = new TH1F("m_id_n_clus_hits","N laben clustered hits for tt1 muons",hits_bins_idm,0,hits_bins_idm*250);
    TH1F* m_id_nclusters[2];
    m_id_nclusters[0] = new TH1F("idm_id_n_clusters","laben clusters in tt1 muons",4,0,4); 
    m_id_nclusters[1] = new TH1F("odm_id_n_clusters","laben clusters in tt2 muons",4,0,4);
    m_idhits[0]->SetDirectory(0); m_idhits[1]->SetDirectory(0); m_idhits[2]->SetDirectory(0);
    m_id_nclusters[0]->SetDirectory(0); m_id_nclusters[1]->SetDirectory(0);

    //13th canvas OD nhits distributions
    TH1F *m_odhits[6];
    int hits_bins_odm = 500;
    char name_13[100];
    sprintf (name_13, "RUN_%d_N muon raw hits in tt1 muons",run);
    m_odhits[0] = new TH1F("idm_od_n_raw_hits", name_13, hits_bins_odm, 0, hits_bins_odm);
    m_odhits[1] = new TH1F("odm_od_n_raw_hits", "N muon raw hits in tt2 muons",hits_bins_odm,0,hits_bins_odm);
    m_odhits[2] = new TH1F("idm_od_n_dec_hits", "N muon decoded hits in tt1 muons",hits_bins_odm,0,hits_bins_odm);
    m_odhits[3] = new TH1F("odm_od_n_dec_hits", "N muon decoded hits in tt2 muons",hits_bins_odm,0,hits_bins_odm);
    m_odhits[4] = new TH1F("idm_od_n_clus_hits","N muon clustered hits in tt1 muons",hits_bins_odm,0,hits_bins_odm);
    m_odhits[5] = new TH1F("odm_od_n_clus_hits","N muon clustered hits in tt2 muons",hits_bins_odm,0,hits_bins_odm);
    TH1F* m_od_nclusters[2];
    m_od_nclusters[0] = new TH1F("idm_od_n_clusters","muon clusters in tt1 muons",14,0,14);
    m_od_nclusters[1] = new TH1F("odm_od_n_clusters","muon clusters in tt2 muons",14,0,14);
    for (int i=0; i<6; i++) m_odhits[i]->SetDirectory(0); 
    m_od_nclusters[0]->SetDirectory(0); m_od_nclusters[1]->SetDirectory(0);


    // 14th muon pulse shapes
    TH1F *mtimes[6];
    int muon_gate_length = 18000;
    int muon_time_bin = 5;
    char name_14[100];
    sprintf (name_14, "RUN_%d_LASER",run);
    mtimes[0] = new TH1F("mtimes_laser", name_14, (int) muon_gate_length/2/muon_time_bin, 0, muon_gate_length/2);
    mtimes[1] = new TH1F("mtimes_muonID", "MUON tt1", muon_gate_length/muon_time_bin, -1 * muon_gate_length, 0);
    mtimes[2] = new TH1F("mtimes_muonOD", "MUON tt2", muon_gate_length/muon_time_bin, -1 * muon_gate_length, 0);
    mtimes[3] = new TH1F("mtimes_neutrino", "NEUTRINO", muon_gate_length/muon_time_bin, -1 * muon_gate_length, 0);
    mtimes[4] = new TH1F("mtimes_random","RANDOM", muon_gate_length/muon_time_bin, -1 * muon_gate_length, 0);
    mtimes[5] = new TH1F("mtimes_pulser","PULSER", muon_gate_length/muon_time_bin, -1 * muon_gate_length, 0);
    mtimes[0]->SetDirectory(0); mtimes[1]->SetDirectory(0); mtimes[2]->SetDirectory(0);
    mtimes[3]->SetDirectory(0); mtimes[4]->SetDirectory(0); mtimes[5]->SetDirectory(0);

    #endif


    int neutrino_found = 0;
    #ifdef muon
      int muon_found = 0;
      int tt2_found = 0;
    #endif
		  
    // loop on all events 
    for(int iev=start; iev <= end; iev++) {    

      //check for broaked entries
      if (tree->GetEntry(iev - evnum_first_entry) <= 0) {
        std::cout << "unable to read entry " << iev - evnum_first_entry << std::endl;
        continue;
      }

      // load iev-event in memory )
      tree->GetEntry(iev - evnum_first_entry);

      // count total number of seconds from first to last event
      t_end = ev->GetTrigger().GetGpsTimeSec();    
      int trgtype = ev->GetTrigger().GetTrgType();
      #ifdef muon
        bool mtb = false;
        if(ev->GetTrigger().GetBtbInputs()==4) mtb = true;
      #endif
      
      //N found neutrino events
      if(trgtype == 1) neutrino_found ++;
      #ifdef muon
        if(mtb) muon_found++;
	if(trgtype == 2) tt2_found++;
      #endif
      
      //starting time is the first neutrino event
      if (neutrino_found == 1) {
        t_start = t_end;
        run = ev->GetRun();
      }


      int evnum_for_this_entry = ev->GetEvNum ();
      #ifdef laben
      double trgtime = ev->GetLaben().GetTriggerTime();
      double lasertime = ev->GetLaben().GetLaserTime();
      #endif
      double muontime = 0;
      double dsize = ev->GetLaben().GetNDecodedHits();
      double rsize = ev->GetLaben().GetNRawHits();
      #ifdef muon
        double mdsize = ev->GetMuon().GetNDecodedHits();
        double mrsize = ev->GetMuon().GetNRawHits();
	double csize  = ev->GetLaben().GetNClusteredHits();
	double mcsize = ev->GetMuon().GetNClusteredHits();
      #endif
 
      
      #ifdef laben
     //times and lg
      for(int ihit = 0; ihit < dsize; ihit ++){
        if(trgtype == 8) {
          htimes[0]->Fill(ev->GetLaben ().GetDecodedHits ()[ihit].GetRawTime () - lasertime);
    	  hlg[0]->Fill(ev->GetLaben ().GetDecodedHits ()[ihit].GetLg ());
        }       
        if(trgtype == 1) {
          htimes[1]->Fill(ev->GetLaben ().GetDecodedHits ()[ihit].GetRawTime () - trgtime);
          hlg[1]->Fill(ev->GetLaben ().GetDecodedHits ()[ihit].GetLg ());
        }
        if(trgtype == 64) {
          htimes[2]->Fill(ev->GetLaben ().GetDecodedHits ()[ihit].GetRawTime () - trgtime);
          hlg[2]->Fill(ev->GetLaben ().GetDecodedHits ()[ihit].GetLg ());
        }
        if(trgtype == 32) {
           htimes[3]->Fill(ev->GetLaben ().GetDecodedHits ()[ihit].GetRawTime () - trgtime);
           hlg[3]->Fill(ev->GetLaben ().GetDecodedHits ()[ihit].GetLg ());
        }
      }
      #endif

      #ifdef muon
      //OD times and mch
      for(int ihit = 0; ihit < mdsize; ihit ++){
        if(trgtype == 8) {
          mtimes[0]->Fill(ev->GetMuon().GetDecodedHits()[ihit].GetTime() - muontime);
          hmch[0]->Fill(ev->GetMuon().GetDecodedHits()[ihit].GetMch());
        }
        if (mtb && trgtype == 1) {
          mtimes[1]->Fill(ev->GetMuon().GetDecodedHits()[ihit].GetTime() - muontime);
          hmch[1]->Fill(ev->GetMuon().GetDecodedHits()[ihit].GetMch());
        }
	if (trgtype == 2) {
          mtimes[2]->Fill(ev->GetMuon().GetDecodedHits()[ihit].GetTime() - muontime);
    	  hmch[2]->Fill(ev->GetMuon().GetDecodedHits()[ihit].GetMch());
	}
        if (trgtype == 1 && !mtb) {
          mtimes[3]->Fill(ev->GetMuon().GetDecodedHits()[ihit].GetTime() - muontime);
          hmch[3]->Fill(ev->GetMuon().GetDecodedHits()[ihit].GetMch());
        }
        if(trgtype == 64) {
           mtimes[4]->Fill(ev->GetMuon().GetDecodedHits()[ihit].GetTime() - muontime);
    	   hmch[4]->Fill(ev->GetMuon().GetDecodedHits()[ihit].GetMch());
        }
        if(trgtype == 32) {
           mtimes[5]->Fill(ev->GetMuon().GetDecodedHits()[ihit].GetTime() - muontime);
    	   hmch[5]->Fill(ev->GetMuon().GetDecodedHits()[ihit].GetMch());
        }
      }
      #endif

      #ifdef laben
      //for neutrino triggers
      if(trgtype == 1) {
      
      //for(int ihit = 1; ihit <= rsize; ihit ++)
      //rawlg->Fill(GetLaben ().GetRawHits ()[ihit].GetLg ());
      
        double n_cl = (double) ev->GetLaben ().GetNClusters ();
        n_clusters->Fill (n_cl);
	
        //not to fill for the events from the last bin (if I do not have the whole minute of data)
        if ( (t_end - t_start)  < nbins_neutrino_rate * 60.){
          neutrino_rate[0]->Fill((t_end-t_start), 1./(double) time_bin_rate);
          if(dsize > 100)   neutrino_rate[1]->Fill((t_end-t_start), 1./(double) time_bin_rate);
        }      
        nhits[0]->Fill(rsize);
        nhits[1]->Fill(dsize);
        if(n_cl == 1.){
          nhits[2]->Fill(ev->GetLaben ().GetNClusteredHits ());
          nhits[3]->Fill(dsize);
          charge->Fill(ev->GetLaben ().GetClusters ()[0].GetCharge ());
        }
      }
      #endif
   
      #ifdef muon
      //for muons
      if(mtb) {
        double od_n_cl = (double) ev->GetMuon().GetNClusters();
	double id_n_cl = (double) ev->GetLaben().GetNClusters();
	if (trgtype==1) {
	  m_id_nclusters[0]->Fill(id_n_cl);
	  m_od_nclusters[0]->Fill(od_n_cl);
	}
	if (trgtype==2) {
          m_id_nclusters[1]->Fill(id_n_cl);
	  m_od_nclusters[1]->Fill(od_n_cl);
	}
	//not to fill for the events from the last bin (if I do not have the whole minute of data)
        if ( (t_end - t_start)  < nbins_muon_rate *muon_time_bin_rate){
	  if (trgtype==1) {
            muon_rate[0]->Fill((t_end-t_start), 24.*3600./(double) muon_time_bin_rate);
	    if (mdsize > 25) muon_rate[1]->Fill((t_end-t_start), 24.*3600./(double) muon_time_bin_rate);
	  }
	  if (trgtype==2) {
            muon_rate[2]->Fill((t_end-t_start), 24.*3600./(double) muon_time_bin_rate);
            if (mdsize > 25) { 
	      muon_rate[3]->Fill((t_end-t_start), 24.*3600./(double) muon_time_bin_rate);
	      OD25_muon_cnt++;
	    }
	  }
        }
	if (trgtype==1) {
          m_idhits[0]->Fill(rsize);
	  m_odhits[0]->Fill(mrsize);
          m_idhits[1]->Fill(dsize);
	  m_odhits[2]->Fill(mdsize);
	}
	if (trgtype==2) {
	  m_odhits[1]->Fill(mrsize);
	  m_odhits[3]->Fill(mdsize);
	}
	if (trgtype==1) {
	  if (id_n_cl) m_idhits[2]->Fill(csize);
          if (od_n_cl) m_odhits[4]->Fill(mcsize);
        }
	if (trgtype==2 && od_n_cl) m_odhits[5]->Fill(mcsize);
      }

      // muon flag rates
      double bin_rate_factor = (double) nbins_muon_flags/tot_time_minutes*1440.;
      bool mcr = (ev->GetMuon().GetNClusters()>0);
      double timebin_flag = (t_end-t_start)/(tot_time_minutes*60.)*nbins_muon_flags;
      if (trgtype==1) {
        bool idf = false;
        if (mtb) m_flag_rate[0]->Fill(timebin_flag, bin_rate_factor);
        if (mcr) m_flag_rate[1]->Fill(timebin_flag, bin_rate_factor);
        if (csize>2100 && ev->GetLaben().GetCluster(0).GetMeanTime()>100. && ev->GetLaben().GetRecClusters()[0].GetGatti()<0.55) idf = true;
	if (ev->GetLaben().GetCluster(0).GetPeakTimes()[0]>30. && csize>900 && csize<2100) idf = true;
	if (ev->GetLaben().GetCluster(0).GetPeakTimes()[0]>40. && csize>100 && csize<900) idf = true;
        if (idf) m_flag_rate[2]->Fill(timebin_flag, bin_rate_factor);
      }
      if (trgtype==2) {
	if (mtb) m_flag_rate[3]->Fill(timebin_flag, bin_rate_factor);
	if (mtb && ev->GetMuon().GetNDecodedHits()>25) m_flag_rate[4]->Fill(timebin_flag, bin_rate_factor);
	if (mcr) m_flag_rate[5]->Fill(timebin_flag, bin_rate_factor);
      }
      #endif

      #ifdef laben
      // fill first canvas histos
      switch(trgtype) {
      case 8: {
        h1[0]->Fill(evnum_for_this_entry,dsize);
        laser_cnt++;
        break;
      }
      case 1: {
        h1[1]->Fill(evnum_for_this_entry,dsize);
        neutrino_cnt++;
        break;
      }
      case 64: {
        h1[2]->Fill(evnum_for_this_entry,dsize);
        random_cnt++;
        break;
      }
      case 32: {
        h1[3]->Fill(evnum_for_this_entry,dsize);
        pulser_cnt++;
        break;
      }
      default: break;
      }
      #endif

      #ifdef muon
      switch(trgtype) {
      case 8: {
        m1[0]->Fill(evnum_for_this_entry,mdsize);
        break;
      }
      case 1: {
        if (mtb) {
	  m1[1]->Fill(evnum_for_this_entry,mdsize);
	  ID_muon_cnt++;
	}
	else m1[3]->Fill(evnum_for_this_entry,mdsize);
        break;
      }
      case 2: {
        m1[2]->Fill(evnum_for_this_entry,mdsize);
	OD_muon_cnt++;
	break;
      }
      case 64: {
        m1[4]->Fill(evnum_for_this_entry,mdsize);
        break;
      }
      case 32: {
        m1[5]->Fill(evnum_for_this_entry,mdsize);
        break;
      }
      default: break;
      }
      #endif

      
      //print out each 10000 events
      if (iev%10000==0 && iev) cout << "Event= " << iev << endl; 
     
      if (iev == end) {
        #ifdef laben
        good_laben_channels = ev->GetLaben ().GetNLivePmts (); 
	#endif
        //good_muon_channels =  0; 
        //good_muon_channels =  ev->GetMuon ().GetNPmts (); 
      }
      
    }
    int steplength;
    #ifdef laben
      steplength = histo2_size;
    #endif
    #ifdef muon
      steplength = m2_size;
    #endif

    int step = (end-start) / steplength;
    for(int hindex=0; hindex<6; hindex++) {
      for(int i=0; i<steplength; i++) {
        int start_his = i*step+1;
        int end_his = start_his + step;
      #ifdef laben
        int cnt=0;
        float tot = 0.;
	if (hindex<4) {
          for(int j=start_his; j < end_his; j++) {
            if (j<h1[hindex]->GetSize()) {
              float val = h1[hindex]->GetBinContent(j);
              if (val>0.) {
                tot += val;
                cnt++;
              }
            }
  	  }  
          if (cnt>0) tot /= cnt;
          else tot=0.;
          h2[hindex]->SetBinContent(i+1,tot);
          if(cnt) h2[hindex]->SetBinError(i+1,TMath::Sqrt(tot)/TMath::Sqrt(cnt));     // poisson errors
          if (hindex==3 && cnt ) h2[hindex]->SetBinError(i+1,TMath::Sqrt(tot)/(5.*TMath::Sqrt(cnt)));   // pulser is NOT poisson
	}
      #endif
      #ifdef muon
        int cnt_m=0;
        float tot_m = 0.;
        for (int j=start_his; j < end_his; j++) {
          if (j<m1[hindex]->GetSize()) {
            float val_m = m1[hindex]->GetBinContent(j);
	    if(val_m>0.) {
	      tot_m += val_m;
  	      cnt_m++;
  	    }
	  }																		            }
	if (cnt_m>0) tot_m /= cnt_m;
        else tot_m=0.;
        m2[hindex]->SetBinContent(i+1,tot_m);
	if(cnt_m) m2[hindex]->SetBinError(i+1,TMath::Sqrt(tot_m)/TMath::Sqrt(cnt_m));     // poisson errors
        if (hindex==3 && cnt_m ) m2[hindex]->SetBinError(i+1,TMath::Sqrt(tot_m)/(5.*TMath::Sqrt(cnt_m)));   // pulser is NOT poisson
      #endif
      }
    }

  #ifdef laben
    char c1_title[50];
    sprintf(c1_title,"Run0%5d Event size stability plots",run);
    TCanvas *c1 = new TCanvas("c1",c1_title,800,600);
    FirstCanvas(c1, h1);
    char c2_title[50];
    sprintf(c2_title,"Run0%5d Event size stability plots (2)",run);
    TCanvas *c2 = new TCanvas("c2",c2_title,800,600);
    SecondCanvas(c2, h2);
    char c3_title[50];
    sprintf(c3_title,"Run0%5d Neutrino trigger vs time",run);
    TCanvas *c3 = new TCanvas("c3",c3_title,800,600);
    ThirdCanvas(c3, neutrino_rate, (t_end-t_start));
    char c4_title[50];
    sprintf(c4_title,"Run0%5d Nhits and clusters",run);
    TCanvas *c4 = new TCanvas("c4",c4_title,800,600);
    FourthCanvas(c4, nhits,n_clusters);	
    char c5_title[50];
    sprintf(c5_title,"Run0%5d decoded hits times",run);
    TCanvas *c5 = new TCanvas("c5",c5_title,800,600);
    FifthCanvas(c5, htimes);	
    char c6_title[50];
    sprintf(c6_title,"Run0%5d Lg distribution",run);
    TCanvas *c6 = new TCanvas("c6",c6_title,800,600);
    SixthCanvas(c6, hlg);
    char c7_title[50];
    sprintf(c7_title,"Run0%5d Nhits and charge",run);
    TCanvas *c7 = new TCanvas("c7",c7_title,800,600);
    SeventhCanvas(c7, nhits[0], nhits[1], nhits[2], charge, (t_end-t_start), run);
  #endif
  #ifdef muon
    char c8_title[50];
    sprintf(c8_title,"Run0%5d Muon Flag Stability",run);
    TCanvas *c8 = new TCanvas ("c8", c8_title,800,600);
    EighthCanvas(c8, m_flag_rate, ID_muon_cnt, OD25_muon_cnt);
    char c9_title[50];
    sprintf(c9_title,"Run0%5d Outer Detector Event Size Stability",run);
    TCanvas *c9 = new TCanvas ("c9", c9_title,800,600);
    NinthCanvas(c9, m1);
    char c10_title[50];
    sprintf(c10_title,"Run0%5d Outer Detector Event Size Stability (2)",run);
    TCanvas *c10 = new TCanvas ("c10", c10_title,800,600);
    TenthCanvas(c10, m2, evnum_end);
    char c11_title[50];
    sprintf(c11_title,"Run0%5d Muon Trigger Rates",run);
    TCanvas *c11 = new TCanvas ("c11", c11_title,800,600);
    EleventhCanvas(c11, muon_rate, (t_end-t_start) );
    char c12_title[50];
    sprintf(c12_title,"Run0%5d Muon NHits Distribution in Laben",run);
    TCanvas *c12 = new TCanvas ("c12", c12_title,800,600);
    for (int i=0; i<3; i++) m_idhits[i]->Rebin(Int_t(1000/ID_muon_cnt));
    TwelfthCanvas(c12, m_idhits, m_id_nclusters);
    char c13_title[50];
    sprintf(c13_title,"Run0%5d Muon NHits Distribution in Muon",run);
    TCanvas *c13 = new TCanvas ("c13", c13_title,800,600);
    for (int i=0; i<6; i++) m_odhits[i]->Rebin(Int_t(1000/ID_muon_cnt));
    ThirteenthCanvas(c13, m_odhits, m_od_nclusters);
    char c14_title[50];
    sprintf(c14_title,"Run0%5d Muon Pulse Shape",run);
    TCanvas *c14 = new TCanvas ("c14", c14_title,800,600);
    FourteenthCanvas(c14, mtimes);
    char c15_title[50];
    sprintf(c15_title,"Run0%5d MuonChannel Distributions",run);
    TCanvas *c15 = new TCanvas ("c15", c15_title,800,600);
    FifteenthCanvas(c15, hmch);
  #endif 

    cout << "Run " << run << " statistics:\n";
    cout << "Total events="<<(end-start)<<" start="<<start<<"  end="<<end<<"\n";
    cout << "Total number of secs="<<(t_end-t_start)<<"\n";
    cout << "Neutrino=" << neutrino_cnt;
    #ifdef muon
      cout << "  tt1_Muons=" << ID_muon_cnt;
      cout << "  tt2_Muons=" << OD25_muon_cnt;
    #endif
    cout << "  Laser394="<<laser_cnt<<"  Random="<<random_cnt;
    cout << "  Pulser="<<pulser_cnt<<"\n";
    #ifdef laben
    cout << "Laben: Average Laser Size = " <<    endl;
    VerticalMean(h1[0],3.,true);	
    cout << "Laben: Average Neutrino Size = " << endl;
    VerticalMean(h1[1],3.,true);
    cout << "Laben: Average Random Size = " << endl;
    VerticalMean(h1[2],3.,true);
    cout << "Laben: Average Pulser Size = " << endl;
    VerticalMean(h1[3],3.,true);
    #endif
    #ifdef muon
    cout << "Muon:  Average Laser Size = " <<    endl;
    VerticalMean(m1[0],3.,true);	
    cout << "Muon:  Average Muon Size = " << endl;
    VerticalMean(m1[1],3.,true);
    cout << "Muon:  Average Random Size = " << endl;
    VerticalMean(m1[2],3.,true);
    cout << "Muon:  Average Pulser Size = " << endl;
    VerticalMean(m1[3],3.,true);
    #endif
      
    int btb_thresh = ev->GetTrigger().GetBtbThreshold ();
    
    //NEW PART THAT CHECKS FOR AMAZING OD BEHAVIOUR
    
    #ifdef muon
      int last_event = ev->GetEvNum();
      int mcr_down_event = last_event;
      TTreeResult* res1 = (TTreeResult*)tree->Query("evnum", "!((enabled_crates>>14)&1)");
//      cout << endl << "number of events with disabled mcr: " << res1->GetRowCount() << endl;
      // check if mcr died at some point
      if (res1->GetRowCount()) {
        TSQLRow* row = res1->Next();
        mcr_down_event = atoi(row->GetField(0)) - 1;
//	cout << "check of mcr: mcr_down_event = " << mcr_down_event << endl;
      }
      // check if mcr was disabled from profile and still looks enabled
      else {
	TTreeResult* res2 = (TTreeResult*)tree->Query("evnum", "muon.n_decoded_hits>0 && trigger.trgtype==32");
	if (!res2->GetRowCount()) mcr_down_event = -2;
//	cout << "mcr disabled from profile? number of tt32 events with non-zero muon hits: " << res2->GetRowCount() << endl;
	res2->Delete();
      }
      res1->Delete();

      // case where electronics was working fine but PMTs were off or had problems
      if(mcr_down_event>0) {
        TTreeResult* res3 = (TTreeResult*)tree->Query("evnum", "trigger.trgtype!=32", "", mcr_down_event);
        TTreeResult* res4 = (TTreeResult*)tree->Query("evnum", "muon.n_decoded_hits==0 && trigger.trgtype!=32", "", mcr_down_event);
        if ( res4->GetRowCount() > .4*res3->GetRowCount() ) mcr_down_event = -3;
        res3->Delete();
	res4->Delete();
      }

      int last_aligned_event = mcr_down_event;
      if (last_aligned_event>0) {
        TTreeResult* res3 = (TTreeResult*)tree->Query("evnum", "trigger.trgtype==32 && muon.n_decoded_hits<180", "", mcr_down_event);
        if (res3->GetRowCount()>=5) {
          TSQLRow* row = (TSQLRow*)res3->GetRows()->First();
          int first_disaligned_calib = atoi(row->GetField(0));
          TTreeResult* res4 = (TTreeResult*)tree->Query("evnum", "trigger.trgtype==32 && muon.n_decoded_hits>180", "", first_disaligned_calib);
          if (res4->GetRowCount()) {
            TSQLRow* row = (TSQLRow*)res4->GetRows()->Last();
            last_aligned_event = atoi(row->GetField(0));
          }
          else last_aligned_event = 0;
          res4->Delete();
        }
        res3->Delete();
      }
      cout << endl << "CHECK OF OUTER DETECTOR:" << endl;  
      if (mcr_down_event==last_event) cout << "MCR was working fine for the whole run." << endl;
      else {
        if (mcr_down_event>0) {
          cout << "MCR was down after event " << mcr_down_event << " (of " << last_event << " events, ";
          cout << int(mcr_down_event*100./last_event) << "%)." << endl;
        }
        if (mcr_down_event==-2) cout << "MCR was disabled from profile!" << endl;
        if (mcr_down_event==-3) cout << "Muon PMTs had problems!" << endl;
      }
      if (last_aligned_event > 0.99* mcr_down_event) cout << "No problems with ID/OD-disalignments." << endl;
      else if (mcr_down_event>=0)cout << "Disalignment of ID/OD events, beginning from event " << last_aligned_event+1 << endl;
      cout << endl;

      if (tt2_found==0) cout << "No trigger type 2! Check if MTB is enabled!" << endl;
    #endif
    
    //NEW PART THAT RESCALE THE ELECTRONICS CALIBRATION HISTOGRAMS IF THE CALIBRATION IS PRESENT
    #ifdef laben
    if(calib){
      char raw_lg_title2[50];
      sprintf(raw_lg_title2,"Run0%5d: lg of raw hits,  NORMALIZED to Ntriggers",run);
      TCanvas *raw_lg_renorm = new TCanvas("raw_lg_renorm",raw_lg_title2);
      raw_lg_renorm->Divide(1,2);
      raw_lg_renorm->cd(1);
      raw_lg_pulser->Scale( 1.0 / (float)pulser_cnt);
      raw_lg_pulser->Draw ();
      raw_lg_renorm->cd(2);
      raw_lg_neutrino->Scale( 1.0 / (float)neutrino_cnt);
      raw_lg_neutrino->Draw (); 
      raw_lg_renorm->Modified();
      raw_lg_renorm->Update ();
    }
    #endif

    //NEW PART THAT COLLECTS THE DATA TO BE SAVED

    char *valid = new char[5];
    sprintf(valid,"true");
    
    char *root_files = new char[100];

    std::string path = GetPath (start_date, cycle);
    int modE5 = run/100000;
    int modE4 = run/10000;
    int modE3 = run/1000;
    if( modE5 && modE5 < 10 ) sprintf (root_files, "%sRun%6d_c%d.root", path.c_str(),run,cycle);
    if( modE4 && modE4 < 10 ) sprintf (root_files, "%sRun0%5d_c%d.root", path.c_str(),run,cycle);
    if( modE3 && modE3 < 10 ) sprintf (root_files, "%sRun0%5d_c%d.root", path.c_str(),run,cycle);
    
    char comment[50];
    sprintf(comment,"comment");

    save_validation_parameters->run=run;
    save_validation_parameters->cycle=cycle;
    save_validation_parameters->RunEvents=(int) RunEvents(run);
    save_validation_parameters->evnum_start=evnum_start;
    save_validation_parameters->evnum_end=evnum_end;
    save_validation_parameters->t_diff=(t_end-t_start);
    save_validation_parameters->btb_thresh=btb_thresh;
    #ifdef laben
    save_validation_parameters->VerticalMean1=VerticalMean(neutrino_rate[0],3.,true );
    save_validation_parameters->good_laben_channels=good_laben_channels;
    #else
    save_validation_parameters->VerticalMean1=0;
    save_validation_parameters->good_laben_channels=0;
    #endif
    save_validation_parameters->good_muon_channels=good_muon_channels;
    #ifdef laben
    save_validation_parameters->VerticalMean2=VerticalMean(h1[0],3.,true); 
    save_validation_parameters->VerticalMean3=VerticalMean(h1[2],3.,true);
    save_validation_parameters->VerticalMean4=VerticalMean(h1[3],3.,true);
    save_validation_parameters->VerticalMean5=VerticalMean(h1[1],3.,true);
    #else
    save_validation_parameters->VerticalMean2=0;
    save_validation_parameters->VerticalMean3=0;
    save_validation_parameters->VerticalMean4=0;
    save_validation_parameters->VerticalMean5=0;
    #endif
    #ifdef muon
    save_validation_parameters->MuonRateID = VerticalMean(muon_rate[0],3.,true);
//    save_validation_parameters->MuonRateOD = VerticalMean(muon_rate[2],3.,true);
    save_validation_parameters->MuonRateOD25 = VerticalMean(muon_rate[3],3,true);
    save_validation_parameters->LastMcrEvent = mcr_down_event;
    save_validation_parameters->LastAlignedEvent = last_aligned_event;
    #else
    save_validation_parameters->MuonRateID = 0;
    save_validation_parameters->MuonRateOD = 0;
    save_validation_parameters->MuonRateOD25 = 0;
    save_validation_parameters->LastMcrEvent = 0;
    save_validation_parameters->LastAlignedEvent = 0;
    #endif
    save_validation_parameters->loc_start_time = start_time;
    save_validation_parameters->loc_validation_time = validation_time;
    save_validation_parameters->loc_valid = valid;
    save_validation_parameters->loc_group = group;
    save_validation_parameters->loc_root_files = root_files;    
   
    #ifdef laben
    for(int u=0;u<4;u++){
      save_validation_parameters->histo1[u] = h1[u];
      save_validation_parameters->histo2[u] = h2[u];
      save_validation_parameters->histo3[u] = hlg[u];
      save_validation_parameters->histo4[u] = nhits[u];
      save_validation_parameters->histo5[u] = htimes[u];
    }
    for(int u=0;u<2;u++) save_validation_parameters->histo6[u]=neutrino_rate[u];
    save_validation_parameters->histo7=charge;
    #endif
    #ifdef muon
    for (int u=0; u<6; u++) {
      save_validation_parameters->histo8[u]  = m_flag_rate[u];
      save_validation_parameters->histo9[u]  = m1[u];
      save_validation_parameters->histo10[u] = m2[u];
      save_validation_parameters->histo13[u] = m_odhits[u];
      save_validation_parameters->histo14[u] = mtimes[u];  
      save_validation_parameters->histo15[u] = hmch[u];
    }
    for (int u=0; u<4; u++) save_validation_parameters->histo11[u] = muon_rate[u];
    for (int u=0; u<3; u++) save_validation_parameters->histo12[u] = m_idhits[u];
    for (int u=0; u<2; u++) {
      save_validation_parameters->histo12a[u] = m_id_nclusters[u];
      save_validation_parameters->histo13a[u] = m_od_nclusters[u];
    }
    #endif
    #ifdef laben
    save_validation_parameters->canvas1=c1;
    save_validation_parameters->canvas2=c2;
    save_validation_parameters->canvas3=c3;
    save_validation_parameters->canvas4=c4;
    save_validation_parameters->canvas5=c5;
    save_validation_parameters->canvas6=c6;
    save_validation_parameters->canvas7=c7;
    #endif
    #ifdef muon
    save_validation_parameters->canvas8=c8;
    save_validation_parameters->canvas9=c9;
    save_validation_parameters->canvas10=c10;
    save_validation_parameters->canvas11=c11;
    save_validation_parameters->canvas12=c12;
    save_validation_parameters->canvas13=c13;
    save_validation_parameters->canvas14=c14;
    save_validation_parameters->canvas15=c15;
    #endif
  }

  cout << "Check if everything is ok following the normal procedure." << endl;
  cout << "If the Run is acceptable, save it by typing: SaveValidation()" << endl; /**/
  return; 
}


void SaveValidation(){

  cout << "Run Validation saving procedure starting." << endl;

  cout << "Are you sure to save this run: " << save_validation_parameters->run << endl;
  cout << "To do it please enter again the RunNumber: ";
  int shifter_choose;
  bool go_on_and_save=false;
  cin >> shifter_choose;
  if (shifter_choose==save_validation_parameters->run){
    go_on_and_save=true;
    cout << endl << "Start saving Run Validation of: " << save_validation_parameters->run << endl;
  }else{
    cout << "Error, RunNumber mismatch !!! Check what Run you want to save !!!" << endl;
    return;
  }
/*  cout << "Please enter your run comment. If everything is fine, just enter \"ok\": ";
  char comment[256];
  cin.getline(comment,256);
  cout << "You entered as comment: " << comment << endl;*/
  char comment[10] = "ok";

  if (go_on_and_save){
    char str_charge[100];
    #ifdef laben
    sprintf (str_charge, "/home/production/run_plots/C1_%d.ps", save_validation_parameters->run); 
    save_validation_parameters->canvas1->Print(str_charge);
    sprintf (str_charge, "/home/production/run_plots/C2_%d.ps", save_validation_parameters->run); 
    save_validation_parameters->canvas2->Print(str_charge);  
    sprintf (str_charge, "/home/production/run_plots/C3_%d.ps", save_validation_parameters->run); 
    save_validation_parameters->canvas3->Print(str_charge);
    sprintf (str_charge, "/home/production/run_plots/C4_%d.ps", save_validation_parameters->run); 
    save_validation_parameters->canvas4->Print(str_charge);
    sprintf (str_charge, "/home/production/run_plots/C5_%d.ps", save_validation_parameters->run); 
    save_validation_parameters->canvas5->Print(str_charge);
    sprintf (str_charge, "/home/production/run_plots/C6_%d.ps", save_validation_parameters->run); 
    save_validation_parameters->canvas6->Print(str_charge);
    sprintf (str_charge, "/home/production/run_plots/C7_%d.ps", save_validation_parameters->run); 
    save_validation_parameters->canvas7->Print(str_charge);
    #endif
    #ifdef muon
    /*sprintf (str_charge, "/home/wurm/offline/Echidna/rv/C8_%d.ps", save_validation_parameters->run); 
      save_validation_parameters->canvas8->Print(str_charge);*/
    sprintf (str_charge, "/home/production/run_plots/C8_%d.ps", save_validation_parameters->run); 
    save_validation_parameters->canvas8->Print(str_charge);
    sprintf (str_charge, "/home/production/run_plots/C9_%d.ps", save_validation_parameters->run); 
    save_validation_parameters->canvas9->Print(str_charge);  
    sprintf (str_charge, "/home/production/run_plots/C10_%d.ps", save_validation_parameters->run); 
    save_validation_parameters->canvas10->Print(str_charge);
    sprintf (str_charge, "/home/production/run_plots/C11_%d.ps", save_validation_parameters->run); 
    save_validation_parameters->canvas11->Print(str_charge);
    sprintf (str_charge, "/home/production/run_plots/C12_%d.ps", save_validation_parameters->run); 
    save_validation_parameters->canvas12->Print(str_charge);
    sprintf (str_charge, "/home/production/run_plots/C13_%d.ps", save_validation_parameters->run); 
    save_validation_parameters->canvas13->Print(str_charge);
    sprintf (str_charge, "/home/production/run_plots/C14_%d.ps", save_validation_parameters->run); 
    save_validation_parameters->canvas14->Print(str_charge);
    sprintf (str_charge, "/home/production/run_plots/C15_%d.ps", save_validation_parameters->run); 
    save_validation_parameters->canvas15->Print(str_charge);
    #endif
    
    cout << "Saving on run_validation_out file" << endl;
    FILE *fp=fopen("/home/production/run_validation_out.txt","a+");
    fprintf(fp,"%5d %3d %7d %7d %7d %6d %6d %8.1f %6.0f %6.0f %7d %6d %6d %6.0f %6.0f %7.0f %6.0f %7d %7d %20s %20s %5s %20s %s     %s\n",save_validation_parameters->run, save_validation_parameters->cycle, save_validation_parameters->RunEvents,save_validation_parameters->evnum_start,save_validation_parameters->evnum_end,save_validation_parameters->t_diff,save_validation_parameters->btb_thresh, save_validation_parameters->VerticalMean1, save_validation_parameters->MuonRateID, save_validation_parameters->MuonRateOD25, save_validation_parameters->good_laben_channels, save_validation_parameters->good_muon_channels, 0,save_validation_parameters->VerticalMean2, save_validation_parameters->VerticalMean3, save_validation_parameters->VerticalMean4,save_validation_parameters->VerticalMean5, save_validation_parameters->LastMcrEvent,save_validation_parameters->LastAlignedEvent, save_validation_parameters->loc_start_time, save_validation_parameters->loc_validation_time, save_validation_parameters->loc_valid, save_validation_parameters->loc_group, save_validation_parameters->loc_root_files,comment);   
    fclose(fp);
    
    cout << "    Saved" << endl;
    
//    TH1F **h1=save_validation_parameters->histo1;
    
    char histofile[100];
    sprintf (histofile,"/home/production/run_plots/%d.root", save_validation_parameters->run);
    //sprintf (histofile,"/home/wurm/offline/Echidna/rv/%d.root", save_validation_parameters->run);
    TFile *fh =new TFile (histofile, "RECREATE");
    fh->cd();
    for (int i=0;i<6;i++) { 
      //h1[1]->Write();
      #ifdef laben
      if (i<4) save_validation_parameters->histo1[i]->Write();
      if (i<4) save_validation_parameters->histo2[i]->Write();
      if (i<4) save_validation_parameters->histo3[i]->Write();
      if (i<4) save_validation_parameters->histo4[i]->Write();
      if (i<4) save_validation_parameters->histo5[i]->Write();	
      if (i<2) save_validation_parameters->histo6[i]->Write();
      if (i<1) save_validation_parameters->histo7->Write();
      #endif
      #ifdef muon
      save_validation_parameters->histo8[i]->Write();
      save_validation_parameters->histo9[i]->Write();
      save_validation_parameters->histo10[i]->Write();
      if (i<4)  save_validation_parameters->histo11[i]->Write();
      if (i<3)  save_validation_parameters->histo12[i]->Write();
      if (i<2)  save_validation_parameters->histo12a[i]->Write();
      save_validation_parameters->histo13[i]->Write();
      if (i<2)  save_validation_parameters->histo13a[i]->Write();
      save_validation_parameters->histo14[i]->Write();
      save_validation_parameters->histo15[i]->Write();
      #endif
      cout << "saving: " << i << endl; 
    }
    
    fh->Close();

    
    cout << "Run Validation saving procedure ended." << endl;
    
  }
  
}
