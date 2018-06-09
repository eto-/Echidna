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
 int  run[dim], 
  n_dead[dim], n_hot[dim], n_low[dim], n_retr[dim],
  p_dead[dim], p_hot[dim], p_low[dim], p_retr[dim],
  l_dead[dim], l_hot[dim], l_low[dim], l_retr[dim];
 int dim ;
 ifstream fdat("channel.dat");
 int i = 0;
 string channel;
 int nruns;
 fdat >> channel >> nruns;
 while(!fdat.eof()) {
   
   fdat >> run[i] >>   n_dead[i] >>p_dead[i] >> l_dead[i] 
        	      >> n_hot[i] >> p_hot[i] >> l_hot[i] 
		      >> n_low[i] >> p_low[i] >> l_low[i] 
        	      >> n_retr[i] >> p_retr[i] >> l_retr[i] ;
   if(fdat.eof()) break ;
   
   i++;
 } 

 i--;
 TGraph *_ndead = new TGraph(i,run,n_dead);
 TGraph *_nhot = new TGraph(i,run,n_hot);
 TGraph *_nlow = new TGraph(i,run,n_low);
 TGraph *_nretr = new TGraph(i,run,n_retr);
 
 TGraph *_pdead = new TGraph(i,run,p_dead);
 TGraph *_phot = new TGraph(i,run,p_hot);
 TGraph *_plow = new TGraph(i,run,p_low);
 TGraph *_pretr = new TGraph(i,run,p_retr);
 
 TGraph *_ldead = new TGraph(i,run,l_dead);
 TGraph *_lhot = new TGraph(i,run,l_hot);
 TGraph *_llow = new TGraph(i,run,l_low);
 TGraph *_lretr = new TGraph(i,run,l_retr);


  TCanvas can("can","can",150,10,990,760);

  TPaveText *txt = new  TPaveText(0.9,0.9,1,1,"ndc");
  txt->SetTextSize(0.04);
  txt->AddText(channel.c_str());
  txt->SetBorderSize(0);
  txt->SetFillColor(0);
  
  can.Divide(2,2);
 int MAX = 2 ;
 can.cd(1);
 _ndead->SetTitle("Dead");
 _ndead->GetYaxis()->SetRangeUser(0,MAX);
 _ndead->Draw("A*");
 _pdead->SetMarkerColor(2);
 _pdead->Draw("*,same");
 _ldead->SetMarkerColor(4);
 _ldead->Draw("*,same");
 txt->Draw();

 can.cd(2); 
 _nhot->SetTitle("Hot");
 _nhot->GetYaxis()->SetRangeUser(0,MAX);
 _nhot->Draw("A*");
 _phot->SetMarkerColor(2);
 _phot->SetMarkerStyle(2);
 _phot->Draw("*,same");
 _lhot->SetMarkerColor(4);
 _lhot->SetMarkerStyle(3);
 _lhot->Draw("*,same");
 txt->Draw();

 can.cd(3);
 _nlow->SetTitle("Low");
 _nlow->GetYaxis()->SetRangeUser(0,MAX);
 _nlow->Draw("A*");
 _plow->SetMarkerColor(2);
 _plow->Draw("*,same");
 _llow->SetMarkerColor(4);
 _llow->Draw("*,same");
 txt->Draw();
 
 can.cd(4);
 _nretr->SetTitle("Retriggers");
 _nretr->GetYaxis()->SetRangeUser(0,MAX);
 _nretr->Draw("A*");
 _pretr->SetMarkerColor(2);
 _pretr->Draw("*,same");
 _lretr->SetMarkerColor(4);
 _lretr->Draw("*,same");
 txt->Draw();
 
 can.cd(2);

 TLegend  *leg = new TLegend(0.9,0.9,1,1);
 leg->AddEntry(_ndead,"neutrino","p");
 leg->AddEntry(_pdead, "pulser", "p");
 leg->AddEntry(_ldead,"laser","p");
 leg->SetFillColor(0);
 leg->SetBorderSize(0);
 leg->Draw();
 
 string filename ;
 filename  = "channel_";
 filename  += channel;
 filename += ".gif";
 cout <<filename<< endl ;
 can->SaveAs(filename.c_str());
 
 fdat.close();
 fdat.clear();

}
