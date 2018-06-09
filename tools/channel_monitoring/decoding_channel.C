{
 gROOT->Reset();
 #include <TMath.h>
 #include <string.h>
 #include <stdio.h>
 #include <fstream>
 gStyle->SetOptStat(0);
 gStyle->SetOptFit(0);
 gStyle->SetOptTitle(1);
 //gStyle->SetOptLogy(1);
 int dim = 1000;
 int  run[dim], nohits[dim], toohits[dim];
 float  peak[dim], rms[dim], runf[dim] ;
 int dim ;
 ifstream fdat;
 fdat.open("channel.dat",ios::in);
 int i = 0;
 string channel_no;
 int files;
 fdat >> channel_no >> files ;
 
 while(!fdat.eof()) {
   fdat >> run[i] >> nohits[i] >> toohits[i] >>peak[i] >> rms[i] ;
   runf[i] = float (run[i]);
   if(fdat.eof()) break ;
   i++;
 } 
 
 TCanvas *can = new TCanvas ("can",channel_no.c_str() ,150,10,990,760);
 TGraph *_nohits = new TGraph(i,run,nohits);
 TGraph *_toohits = new TGraph(i,run,toohits);
 TGraph *_peak = new TGraph(i,runf,peak);
 TGraph *_rms = new TGraph(i,runf,rms);
 TPaveText *txt = new  TPaveText(0.9,0.9,1,1,"ndc");
 txt->SetTextSize(0.04);
 txt->AddText(channel_no.c_str());
 txt->SetBorderSize(0);
 txt->SetFillColor(0);

 can->Divide(2,2);
 can->cd(1);
 _nohits->SetTitle("no hit");
 _nohits->Draw("A*");
 txt->Draw();

 can->cd(2);
 _toohits->SetTitle("hot channel");
 _toohits->Draw("A*");
 txt->Draw();

 can->cd(3);
 _peak->SetTitle("peak");
 _peak->SetMaximum(200.);
 _peak->Draw("A*");
  txt->Draw();

 can->cd(4);
 _rms->SetTitle("rms"); 
 _rms->SetMaximum(60.);
 _rms->Draw("A*");
 txt->Draw();
 
 string filename ;
 filename  = "channel_";
 filename  += channel_no;
  filename += ".gif";
 cout <<filename<< endl ;
 fdat.close();
 fdat.clear();
 can->SaveAs(filename.c_str());
}
