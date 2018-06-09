#include <TTree.h>
#include <TBranch.h>
#include <TROOT.h>
#include <TFile.h>
#include <iostream>
#include <fstream>
#include <TH1F.h>
#include <event/BxEvent.hh>
#include <TCanvas.h>
#include <TStyle.h>

void create_histograms(char *listname="files_c11.lst", char* filename="histos_c11.root") {
	
	TFile *f = new TFile(filename,"recreate");
	
	TH1F *h_rawhits = new TH1F("hrh","n raw hits",2000,0,2000);
	TH1F *h_dechits = new TH1F("hdh","n decoded hits",2000,0,2000);
	TH1F *h_cluhits = new TH1F("hch","n clustered hits",2000,0,2000);
	
	TH1F *m_rawhits = new TH1F("mrh","n raw hits - muons",2500,0,50000);
	TH1F *m_dechits = new TH1F("mdh","n decoded hits - muons",2500,0,50000);
	TH1F *m_cluhits = new TH1F("mch","n clustered hits - muons",2500,0,50000);
	TH1F *m_costheta = new TH1F("mct","zenith angle - muons",60,-1,1);
	TH1F *m_phi = new TH1F("mph","azimuth angle - muons",60,0,360);
	
	
	ifstream listfile(listname);
	
	char dstname[256];
	
	listfile >> dstname;
	
	while (!listfile.eof()) {
		cout << dstname << endl;
		TFile *rf = TFile::Open(dstname,"read");
		TTree* tree = (TTree*) rf->Get("bxtree");
		Int_t events = tree->GetEntries();
		TBranch *b = (TBranch*) tree->GetBranch("events");
		BxEvent *ev = new BxEvent();
		b->SetAddress(&ev);	
		
		for (int iev=0; iev<events; iev++) {
			tree->GetEntry(iev);
			if (iev%10000==0) cout << iev << endl;
			if (ev->GetTrigger().GetTrgType()!=1) continue;
			if (ev->GetTrigger().GetBtbInputs()!=4) {
				h_rawhits->Fill(ev->GetLaben().GetNRawHits());
				h_dechits->Fill(ev->GetLaben().GetNDecodedHits());
				h_cluhits->Fill(ev->GetLaben().GetNClusteredHits());
			}
			/*if (ev->GetTrigger().GetBtbInputs()==4)*/ else {
				m_rawhits->Fill(ev->GetLaben().GetNRawHits());
				m_dechits->Fill(ev->GetLaben().GetNDecodedHits());
				m_cluhits->Fill(ev->GetLaben().GetNClusteredHits());
				if (ev->GetTrack().GetPhi()!=0) {
					m_phi->Fill(ev->GetTrack().GetPhi()/3.1416*180.);
					m_costheta->Fill(cos(ev->GetTrack().GetTheta()));
				}
			}
		}
		
		rf->Close();		
		listfile >> dstname;
		
	}
	
	f->Write();
	f->Close();
}



void compare_histograms(char *file_c11="histos_c11.root", char *file_c12="histos_c12.root", bool write=false) {

	gROOT->Reset();
	gROOT->SetStyle("Plain");

	// open histograms
	
	TFile *f11 = TFile::Open(file_c11,"read");
	
	TH1F *hrh11 = (TH1F*) f11->Get("hrh");
	TH1F *hdh11 = (TH1F*) f11->Get("hdh");
	TH1F *hch11 = (TH1F*) f11->Get("hch");
	TH1F *mrh11 = (TH1F*) f11->Get("mrh");
	TH1F *mdh11 = (TH1F*) f11->Get("mdh");
	TH1F *mch11 = (TH1F*) f11->Get("mch");
	TH1F *mph11 = (TH1F*) f11->Get("mph");
	TH1F *mct11 = (TH1F*) f11->Get("mct");
	
	TFile *f12 = TFile::Open(file_c12,"read");
	
	TH1F *hrh12 = (TH1F*) f12->Get("hrh");
	TH1F *hdh12 = (TH1F*) f12->Get("hdh");
	TH1F *hch12 = (TH1F*) f12->Get("hch");
	TH1F *mrh12 = (TH1F*) f12->Get("mrh");
	TH1F *mdh12 = (TH1F*) f12->Get("mdh");
	TH1F *mch12 = (TH1F*) f12->Get("mch");
	TH1F *mph12 = (TH1F*) f12->Get("mph");
	TH1F *mct12 = (TH1F*) f12->Get("mct");

	// residuals histograms
	
	TH1F *dhrh = new TH1F("dhrh","residuals in raw hits",2000,0,2000);
	TH1F *dhdh = new TH1F("dhdh","residuals in decoded hits",2000,0,2000);
	TH1F *dhch = new TH1F("dhch","residuals in clustered hits",2000,0,2000);
	TH1F *dmrh = new TH1F("dmrh","residuals in raw hits - muons",2500,0,50000);
	TH1F *dmdh = new TH1F("dmdh","residuals in decoded hits - muons",2500,0,50000);
	TH1F *dmch = new TH1F("dmch","residuals in clustered hits - muons",2500,0,50000);
	TH1F *dmph = new TH1F("dmph","residuals in azimuth angle",60,0,360);
	TH1F *dmct = new TH1F("dmct","residuals in zenith angle",60,-1,1);	
	
	dhrh->Add(hrh11);
	dhdh->Add(hdh11);
	dhch->Add(hch11);
	dmrh->Add(mrh11);
	dmdh->Add(mdh11);
	dmch->Add(mch11);
	dmph->Add(mph11);
	dmct->Add(mct11);
	dhrh->Add(hrh12,-1);
	dhdh->Add(hdh12,-1);
	dhch->Add(hch12,-1);
	dmrh->Add(mrh12,-1);
	dmdh->Add(mdh12,-1);
	dmch->Add(mch12,-1);
	dmph->Add(mph12,-1);
	dmct->Add(mct12,-1);
	
	TFile *fout;
	if (write) fout = TFile::Open("nhits_and_muon_histos.root","recreate");
	
	TCanvas *ch = new TCanvas("ch","nhits histograms laben",1);
	ch->Divide(3,2);
	
	ch->cd(1);
	hrh12->GetYaxis()->SetRangeUser(1,1e6);
	hrh11->GetXaxis()->SetRangeUser(0,1000);
	hrh12->SetLineColor(2);
	hrh11->SetXTitle("nhits");
	hrh11->Draw();
	hrh12->Draw("same");
	
	ch->cd(2);
	hdh12->GetYaxis()->SetRangeUser(1,1e6);
	hdh11->GetXaxis()->SetRangeUser(0,1000);
	hdh12->SetLineColor(2);
	hdh11->SetXTitle("nhits");
	hdh11->Draw();
	hdh12->Draw("same");
	
	ch->cd(3);
	hch12->GetYaxis()->SetRangeUser(1,1e6);
	hch11->GetXaxis()->SetRangeUser(0,1000);
	hch12->SetLineColor(2);
	hch11->SetXTitle("nhits");
	hch11->Draw();
	hch12->Draw("same");
	
	ch->cd(4);
	//gStyle->SetOptLogy(0);
	dhrh->GetYaxis()->SetRangeUser(-200,20);
	dhrh->GetXaxis()->SetRangeUser(0,1000);		
	dhrh->Draw();
	dhrh->SetXTitle("nhits");
	
	ch->cd(5);
	dhdh->GetYaxis()->SetRangeUser(-200,20);
	dhdh->GetXaxis()->SetRangeUser(0,1000);		
	dhdh->Draw();
	dhdh->SetXTitle("nhits");
	
	ch->cd(6);
	dhch->GetYaxis()->SetRangeUser(-2000,2000);
	dhch->GetXaxis()->SetRangeUser(0,1000);		
	dhch->Draw();
	dhch->SetXTitle("nhits");
	
	if (write) ch->Write();
	
	TCanvas *cm = new TCanvas("cm","nhits histograms - muons",1);
	cm->Divide(3,2);

	cm->cd(1);
	mrh12->SetLineColor(2);
	mrh11->Draw();
	mrh11->SetXTitle("nhits");
	mrh12->Draw("same");
	
	cm->cd(2);
	mdh12->SetLineColor(2);
	mdh11->SetXTitle("nhits");
	mdh11->Draw();
	mdh12->Draw("same");
	
	cm->cd(3);
	mch12->SetLineColor(2);
	mch11->SetXTitle("nhits");
	mch11->Draw();
	mch12->Draw("same");
	
	cm->cd(4);
	dmrh->SetXTitle("nhits");
	dmrh->Draw();
	
	cm->cd(5);
	dmdh->SetXTitle("nhits");
	dmdh->Draw();
	
	cm->cd(6);
	dmch->SetXTitle("nhits");
	dmch->Draw();
	
	if (write) cm->Write();
	
	TCanvas *ct = new TCanvas("ct","track histograms",1);
	ct->Divide(2,2);
	
	ct->cd(1);
	mph12->SetLineColor(2);
	mph11->SetXTitle("#phi");
	mph11->Draw();
	mph12->Draw("same");
	
	ct->cd(2);
	mct12->SetLineColor(2);
	mct11->SetXTitle("cos(#theta)");
	mct11->Draw();
	mct12->Draw("same");
	
	ct->cd(3);
	dmph->SetXTitle("#phi");
	dmph->Draw();
	
	ct->cd(4);
	dmct->SetXTitle("cos(#theta)");
	dmct->Draw();
	
	if (write) { ct->Write(); fout->Close(); }	
	
	
	
}

