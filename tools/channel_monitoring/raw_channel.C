{
 gROOT->Reset();

 #include <vector>
 #include <TMath.h>
 #include <string.h>
 #include <stdio.h>
 #include <fstream>
 gStyle->SetOptStat(0);
 gStyle->SetOptFit(0);
 gStyle->SetOptTitle(1);
 //gStyle->SetOptLogy(1);
 int dim = 1000;
 string channel ;
 float  run[dim], empty[dim], full[dim];
 float neutrino[dim] , laser[dim], pulser[dim];
 int dim ;
 ifstream fdat("channel.dat");
 int i = 0;
 fdat >> channel ;
 while(!fdat.eof()) {   
   fdat >> run[i] >> full[i] >> empty[i] >> neutrino[i] >> laser[i] >> pulser[i];
   if(fdat.eof()) break ;

   i++;
 } 
 i--;
 TGraph *_empty    = new TGraph(i,run,empty);
 TGraph *_full     = new TGraph(i,run,full);
 TGraph *_neutrino = new TGraph(i,run,neutrino);
 TGraph *_laser    = new TGraph(i,run,laser);
 TGraph *_pulser   = new TGraph(i,run,pulser);

 dim -= 1;
 TCanvas *can = new TCanvas ("can","Channel monitoring" ,150,10,990,760);
  can.Divide(1,2);
  can.cd(1);
 string title1 = "FIFO - channel: " + channel; 
 _empty->SetTitle(title1.c_str());
 _empty->GetYaxis()->SetRangeUser(0.,2.);
 _empty->SetMarkerStyle(2);
 _empty->Draw("A*");

 _full->SetMarkerColor(2);
 _full->SetMarkerStyle(3);
 _full->Draw("*,same");
 
 can.cd(2);
 float maxy  = max(_neutrino->GetYaxis()->GetXmax(),_laser->GetYaxis()->GetXmax());
 maxy = max(maxy,_pulser->GetYaxis()->GetXmax());


 _neutrino->GetYaxis()->SetRangeUser(0.,maxy + 5.);
 
 string title2 = "Loosing Hits - channel:  " + channel; 
 _neutrino->SetTitle(title2.c_str());
 _neutrino->SetMarkerStyle(5);
 _neutrino->SetMarkerColor(3);
 _neutrino->Draw("A*");
 
 _laser->SetMarkerStyle(5);
 _laser->SetMarkerColor(4);
 _laser->Draw("*,same");
 
 _pulser->SetMarkerStyle(5);
 _pulser->SetMarkerColor(5);
 _pulser->Draw("*,same");
 
 TLegend  *leg = new TLegend(0.8,0.7,1,1);
 leg->AddEntry(_empty,"fifo empty","p");
 leg->AddEntry(_full, "fifo full", "p");
 leg->AddEntry(_neutrino,"neutrino: loosing hits","p");
 leg->AddEntry(_laser,"laser: loosing hits","p");
 leg->AddEntry(_pulser,"pulser: loosing hits","p");
 leg->SetFillColor(0);
 leg->SetBorderSize(0);
 can.cd(1);
 leg->Draw();
 string filename = "fifo" + channel + ".gif";
 cout <<filename<< endl ;
 can->SaveAs(filename.c_str());

 fdat.close();
 fdat.clear();

}
