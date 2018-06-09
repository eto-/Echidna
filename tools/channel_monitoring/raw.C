{
 gROOT->Reset();


 #include <TMath.h>
 #include <string.h>
 #include <stdio.h>
 #include <fstream>
 gStyle->SetOptStat(0);
 gStyle->SetOptFit(0);
 gStyle->SetOptTitle(0);
 //gStyle->SetOptLogy(1);
 int dim = 2241;
 int  channel[dim], empty[dim], full[dim], neutrino[dim] , laser[dim], pulser[dim];
 int dim ;
 ifstream fdat("general.dat");
 int i = 1;
 TH1I *_empty = new TH1I("_empty","_toomanyhits",2240,1,2241); 
 TH1I *_full    = new TH1I("_full","_problems",2240,1,2241); 
 TH1I *_neutrino = new TH1I("_neutrino","_toomanyhits",2240,1,2241); 
 TH1I *_laser    = new TH1I("_laser","_problems",2240,1,2241); 
 TH1I *_pulser     = new TH1I("_pulser","_nohits",2240,1,2241); 
 while(!fdat.eof()) {
   
   fdat >> channel[i] >> full[i] >> empty[i] >> neutrino[i] >> laser[i] >> pulser[i];
   if(fdat.eof()) break ;
   _empty->Fill(channel[i],empty[i]);
   _empty->SetBinError(channel[i],0.05);   
  
   _full->Fill(channel[i],full[i]);
   _full->SetBinError(channel[i],0.05);   
   
   _neutrino->Fill(channel[i],neutrino[i]);
   _neutrino->SetBinError(channel[i],0.05);   
   
   _laser->Fill(channel[i],laser[i]);
   _laser->SetBinError(channel[i],0.05);   

   _pulser->Fill(channel[i],pulser[i]);
   _pulser->SetBinError(channel[i],0.05);   

   i++;
 } 
 dim -= 1;
 TCanvas *can = new TCanvas ("can","Channel monitoring" ,150,10,990,760);
 
 _empty->SetMarkerStyle(2);
 _empty->Draw("e");

 _full->SetMarkerColor(2);
 _full->SetMarkerStyle(3);
 _full->Draw("e,same");
 
 _neutrino->SetMarkerStyle(5);
 _neutrino->SetMarkerColor(3);
 _neutrino->Draw("e,same");
 
 _laser->SetMarkerStyle(5);
 _laser->SetMarkerColor(4);
 _laser->Draw("e,same");
 
 _pulser->SetMarkerStyle(5);
 _pulser->SetMarkerColor(5);
 _pulser->Draw("e,same");
 
 TLegend  *leg = new TLegend(0.88,0.7,1,1);
 leg->AddEntry(_empty,"fifo empty","p");
 leg->AddEntry(_full, "fifo full", "p");
 leg->AddEntry(_neutrino,"neutrino: loosing hits","p");
 leg->AddEntry(_laser,"laser: loosing hits","p");
 leg->AddEntry(_pulser,"pulser: loosing hits","p");
 leg->SetFillColor(0);
 leg->SetBorderSize(0);
 leg->Draw();
 
 fdat.close();
 fdat.clear();

}
