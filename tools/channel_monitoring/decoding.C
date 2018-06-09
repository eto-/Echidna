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
 int  channel[dim], nohits[dim], toomanyhits[dim], problems[dim] ;
 int dim ;
 ifstream fdat("general.dat");
 int i = 1;
 TH1I *_nohits      = new TH1I("_nohits","_nohits",2240,1,2241); 
 TH1I *_toomanyhits = new TH1I("_toomanyhits","_toomanyhits",2240,1,2241); 
 TH1I *_problems    = new TH1I("_problems","_problems",2240,1,2241); 
 while(!fdat.eof()) {
   
   fdat >> channel[i] >> nohits[i] >> toomanyhits[i] >> problems[i] ;
   if(fdat.eof()) break ;
   _nohits->Fill(channel[i],nohits[i]);
   _nohits->SetBinError(channel[i],0.05);   
  
   _toomanyhits->Fill(channel[i],toomanyhits[i]);
   _toomanyhits->SetBinError(channel[i],0.05);   
   
   _problems->Fill(channel[i],problems[i]);
   _problems->SetBinError(channel[i],0.05);   
   
   i++;
 } 
 dim -= 1;
 TCanvas *can = new TCanvas ("can","Channel monitoring" ,150,10,990,760);
 
 _nohits->SetMarkerStyle(2);
 _nohits->Draw("e");

 _toomanyhits->SetMarkerColor(2);
 _toomanyhits->SetMarkerStyle(3);
 _toomanyhits->Draw("e,same");
 
 _problems->SetMarkerStyle(5);
 _problems->SetMarkerColor(3);
 _problems->Draw("e,same");
 TLegend  *leg = new TLegend(0.9,0.9,1,1);
 leg->AddEntry(_nohits,"no hit","p");
 leg->AddEntry(_toomanyhits, "hot channel", "p");
 leg->AddEntry(_problems,"problems","p");
 leg->SetFillColor(0);
 leg->SetBorderSize(0);
 leg->Draw();
 
 fdat.close();
 fdat.clear();

}
