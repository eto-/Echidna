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
 int dim = 2242;
 int  channel[dim], 
  n_dead[dim], n_hot[dim], n_low[dim], n_retr[dim],
  p_dead[dim], p_hot[dim], p_low[dim], p_retr[dim],
  l_dead[dim], l_hot[dim], l_low[dim], l_retr[dim];
 int dim ;
 ifstream fdat("general.dat");
 int i = 1;

 while(!fdat.eof()) {
   
   fdat >> channel[i] >>   n_dead[i] >>p_dead[i] >> l_dead[i] 
        	      >> n_hot[i] >> p_hot[i] >> l_hot[i] 
		      >> n_low[i] >> p_low[i] >> l_low[i] 
        	      >> n_retr[i] >> p_retr[i] >> l_retr[i] ;
		  
   if(fdat.eof()) break ;
   
   i++;
 } 


 TGraph *_ndead = new TGraph(i,channel,n_dead);
 TGraph *_nhot  = new TGraph(i,channel,n_hot);
 TGraph *_nlow  = new TGraph(i,channel,n_low);
 TGraph *_nretr = new TGraph(i,channel,n_retr);
 
 TGraph *_pdead = new TGraph(i,channel,p_dead);
 TGraph *_phot  = new TGraph(i,channel,p_hot);
 TGraph *_plow  = new TGraph(i,channel,p_low);
 TGraph *_pretr = new TGraph(i,channel,p_retr);
 
 TGraph *_ldead = new TGraph(i,channel,l_dead);
 TGraph *_lhot  = new TGraph(i,channel,l_hot);
 TGraph *_llow  = new TGraph(i,channel,l_low);
 TGraph *_lretr = new TGraph(i,channel,l_retr);


  TCanvas can("can","can",150,10,990,760);

  
  can.Divide(2,2);
 int MAX = _ndead->GetYaxis()->GetXmax() +10 ;
 can.cd(1);
 _ndead->SetTitle("Dead");
 _ndead->GetYaxis()->SetRangeUser(0,MAX);
 _ndead->Draw("A*");
 _pdead->SetMarkerColor(2);
 _pdead->Draw("*,same");
 _ldead->SetMarkerColor(4);
 _ldead->Draw("*,same");

 can.cd(2); 
 _nhot->SetTitle("Hot");
 _nhot->GetYaxis()->SetRangeUser(0,MAX);
 _nhot->Draw("A*");
 _phot->SetMarkerColor(2);
 _phot->Draw("*,same");
 _lhot->SetMarkerColor(4);
 _lhot->Draw("*,same");

 can.cd(3);
 _nlow->SetTitle("Low");
 _nlow->GetYaxis()->SetRangeUser(0,MAX);
 _nlow->Draw("A*");
 _plow->SetMarkerColor(2);
 _plow->Draw("*,same");
 _llow->SetMarkerColor(4);
 _llow->Draw("*,same");
 
 can.cd(4);
 _nretr->SetTitle("Retriggers");
 _nretr->GetYaxis()->SetRangeUser(0,MAX);
 _nretr->Draw("A*");
 _pretr->SetMarkerColor(2);
 _pretr->Draw("*,same");
 _lretr->SetMarkerColor(4);
 _lretr->Draw("*,same");
 
 can.cd(2);

 TLegend  *leg = new TLegend(0.9,0.9,1,1);
 leg->AddEntry(_ndead,"neutrino","p");
 leg->AddEntry(_pdead, "pulser", "p");
 leg->AddEntry(_ldead,"laser","p");
 leg->SetFillColor(0);
 leg->SetBorderSize(0);
 leg->Draw();
 
 fdat.close();
 fdat.clear();

}
