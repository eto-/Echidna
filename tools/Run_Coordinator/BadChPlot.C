#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <sstream>
#include <string>
#include "TRandom.h"
#include "TH1F.h"
#include "TFile.h"
#include "TTree.h"
#include "TDirectory.h"
#include "TMath.h" 
#include <TH2F.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TDirectoryFile.h>
#include <TStyle.h>
#include <TSystem.h>
#include "TROOT.h"
#include <unistd.h>

/*
	Author for Muon Plots:	Quirin Meindl
	Date: 19/06/2012
*/


int find_last_calibration_run(const std::string& tempfile_path, int &calibration_run, char *calib_path, const int run){

	//Find latest laser calibrations
	ifstream calib_file(tempfile_path.c_str());
	if (!calib_file){
		cerr<<"Error opening "<<tempfile_path<<"\n";
		return 1;
	}

	char path[1024];
	char line[1024];
	char runnumber_char[6];
	int runnumber;

	bool find_last_calib_run = false;
	if (calibration_run == 0) {
		find_last_calib_run = true;
		std::cout << " --- No Calibration Run specified. Searching for last Calibration Run before Run " << run << ". --- " << std::endl;
	}

	bool calib_run_found = false;
	while(!calib_file.eof()){
		calib_file >> path;
		runnumber_char[0] = path[55];runnumber_char[1] = path[56];runnumber_char[2] = path[57];
		runnumber_char[3] = path[58];runnumber_char[4] = path[59];runnumber_char[5] = path[60];
		runnumber = atoi(runnumber_char);

		if (find_last_calib_run){	//if routine active, find last calibration run
			if (runnumber > calibration_run && runnumber < run) {	//calibration run is required to be before the Normal Run studied in the whole Macro
				calibration_run = runnumber;
				strcpy(calib_path,path);
				calib_run_found = true;
			}
		}
		else {				//if run is fixed, stop when detected
			if (runnumber == calibration_run) {
				strcpy(calib_path,path);
				calib_run_found = true;
				break;
			}
		}
		calib_file.getline(line,1024);
	}

	//delete temporary file
	calib_file.close();

	return !calib_run_found;
}


int BadChPlot(int run, int year, const std::string& week, const std::string& tempfile_path, int calibration_run){ 
	gStyle->SetCanvasColor(kWhite);

	string istoname1, istoname2, istoname3, istoname4, istoname5, istoname6;
	string istomuonname1, istomuonname2, istomuonname3, istomuonname4;
	string figname1, figname2, figname3, figname4, figname5, figname6;
	string figmuonname1, figmuonname2, figmuonname3, figmuonname4, figmuonname3_link, figmuonname4_link;
	stringstream isto1, isto2, isto3, isto4, isto5, isto6;
	stringstream istomuon1, istomuon2, istomuon3, istomuon4;
	stringstream fig1, fig2, fig3, fig4, fig5, fig6;
	stringstream figmuon1, figmuon2, figmuon3, figmuon4, figmuon3_link, figmuon4_link;

	ostringstream tab, rfpre, rfec;
	ostringstream tabmuon, prod;
	tab << "/home/production/dbMon/laben_broken_channels_0" << run << ".txt" ;
	tabmuon << "/home/production/dbMon/muon_broken_channels_0" << run << ".txt" ;
	rfpre << "http://bxmaster-data//bxstorage/rootfiles/cycle_15/" << year << "/" << week << "/ancillary/Run0" << run << "_precalibrations_c15.root";
	rfec << "http://bxmaster-data//bxstorage/rootfiles/cycle_14/" << year << "/" << week << "/ancillary/Run0" << run << "_electronics_calibrations_c14.root";
	prod << "http://bxmaster-data//bxstorage/rootfiles/cycle_14/" << year << "/" << week << "/Run0" << run << "_c14.root";

	std::cout << " --- Processing run " << run << " --- " << std::endl;

	//Find (lastest) Laser Calibration Run
	char calib_path[1024];
	int status = find_last_calibration_run(tempfile_path, calibration_run,calib_path,run);
	if (status == 1) {
		std::cerr << " --- Error : Laser Calibration run " << calibration_run << " not found --- " << std::endl;
		return 1;
	}
	else std::cout << " --- Laser Calibration run " << calibration_run << " found at path : " << calib_path << " --- "<< std::endl;

	const string tabname=tab.str();
	const string tabmuonname=tabmuon.str();
	const string rfprename=rfpre.str();
	const string rfecname=rfec.str();
	const string prodname=prod.str();

	int reas=0, lg=0;
	ifstream ff(tabname.c_str());
	if (!ff){
		cerr<<"Error opening "<<tabname<<"\n";
		return 1;
	}
 
	while(!ff.eof()){
		ff >> lg >> reas;
		//    ff >> reas >> lg;
		if(ff.eof()) break;
		isto1 << "PrecalibDecodingLg" << lg;
		istoname1=isto1.str();
    
		TFile *fil = TFile::Open(rfprename.c_str());
		if (fil->IsZombie()){
			cerr<<"Error opening "<<rfprename<<"\n";
			return 1;
		}
		TH2F* h = (TH2F*) fil->Get("/barn/bx_calib_laben_decoding/time_vs_channel");
		if (!h){
			cerr<<"Error getting TH2F \n";
			return 1;
		}
		TH1D* PrecalibDecoding = (TH1D*) h->ProjectionY(istoname1.c_str(),lg,lg);
		PrecalibDecoding->SetTitle(istoname1.c_str());
		fig1 << "/home/production/BadChannelsPlots/Run0" << run << "_Lg" << lg << "_1.jpg";
		figname1=fig1.str();
		TCanvas *c1=new TCanvas("c1","c1");
		PrecalibDecoding->Draw();
		c1->SaveAs(figname1.c_str());
		isto1.str( "" ); fig1.str( "" ); 


		if(reas==2){
			isto2 << "ChargePulserLg" << lg;
			isto3 << "BasePulserLg" << lg;
			isto4 << "PeakPulserLg" << lg;
			isto5 << "ChargeNeutrinoLg" << lg;
			isto6 << "TimingLaserLg" << lg;
			istoname2=isto2.str();
			istoname3=isto3.str();
			istoname4=isto4.str();
			istoname5=isto5.str();
			istoname6=isto6.str();

			fig2 << "/home/production/BadChannelsPlots/Run0" << run << "_Lg" << lg << "_2.jpg";
			fig3 << "/home/production/BadChannelsPlots/Run0" << run << "_Lg" << lg << "_3.jpg";
			fig4 << "/home/production/BadChannelsPlots/Run0" << run << "_Lg" << lg << "_4.jpg";
			fig5 << "/home/production/BadChannelsPlots/Run0" << run << "_Lg" << lg << "_5.jpg";
			fig6 << "/home/production/BadChannelsPlots/Run0" << run << "_Lg" << lg << "_6.jpg";
			figname2=fig2.str();
			figname3=fig3.str();
			figname4=fig4.str();
			figname5=fig5.str();
			figname6=fig6.str();

			TFile *fg = TFile::Open(rfecname.c_str());
			TH2F* ha = (TH2F*) fg->Get("/barn/bx_calib_laben_electronics/charge_vs_lg");
			TH2F* hb = (TH2F*) fg->Get("/barn/bx_calib_laben_electronics/base_vs_lg");
			TH2F* hc = (TH2F*) fg->Get("/barn/bx_calib_laben_electronics/peak_vs_lg");
			TH2F* hd = (TH2F*) fg->Get("/barn/bx_calib_laben_charge_tt1/bx_adc_charge_tt1");
     
			TH1D* ChargePulser = (TH1D*) ha->ProjectionY(istoname2.c_str(),lg,lg);
			ChargePulser->SetTitle(istoname2.c_str());
			TCanvas *c2=new TCanvas("c2","c2");
			ChargePulser->Draw();
			c2->SaveAs(figname2.c_str());
     
			TH1D* BasePulser = (TH1D*) hb->ProjectionY(istoname3.c_str(),lg,lg);
			BasePulser->SetTitle(istoname3.c_str());
			TCanvas *c3=new TCanvas("c3","c3");
			BasePulser->Draw();
			c3->SaveAs(figname3.c_str());
      
			TH1D* PeakPulser = (TH1D*) hc->ProjectionY(istoname4.c_str(),lg,lg);
			PeakPulser->SetTitle(istoname4.c_str());
			TCanvas *c4=new TCanvas("c4","c4");
			PeakPulser->Draw();
			c4->SaveAs(figname4.c_str());

			TH1D* ChargeNeutrino = (TH1D*) hd->ProjectionY(istoname5.c_str(),lg,lg);
			ChargeNeutrino->SetTitle(istoname5.c_str());
			TCanvas *c5=new TCanvas("c5","c5");
			ChargeNeutrino->Draw();
			c5->SaveAs(figname5.c_str());
      
			isto2.str( "" ); 
			isto3.str( "" );
			isto4.str( "" );
			isto5.str( "" );
			fig2.str( "" ); 
			fig3.str( "" );
			fig4.str( "" );
			fig5.str( "" );
			ha->Delete(); c2->Clear();
			hb->Delete(); c3->Clear();
			hc->Delete(); c4->Clear();
			hd->Delete(); c5->Clear();
			ChargePulser->Delete();
			BasePulser->Delete();
			PeakPulser->Delete();
			ChargeNeutrino->Delete();
    
		}
		h->Delete(); c1->Clear();
		PrecalibDecoding->Delete();
		fil->Close();
	}

 	ff.close();


	//--------------------
	//--------------------
	//NOW CHECK MUON DATA
	//--------------------
	//--------------------

	//(1) get broken muon channels
	ifstream ff_muon(tabmuonname.c_str());
	if (!ff_muon){
		cerr<<"Error opening "<<tabmuonname<<"\n";
		return 1;
	}

	int mch=0;
	char line[1024];
	std::vector<int> broken_muon_channels_vec;
	while(!ff_muon.eof()){
		mch = -1;
		ff_muon >> mch;
		if (mch == -1) break;	//stop at end of the file
		broken_muon_channels_vec.push_back(mch);
		ff_muon.getline(line,1024);
	}
	int n_broken_mch = broken_muon_channels_vec.size();
	if (n_broken_mch == 0) {	//if all muon channels are fine, stop
		std::cout<<" --- No problematic Muon Channels --- "<<std::endl;
		std::cout<<"NORMAL END"<<std::endl;
		return 0;
	}


	//------------
	//PULSER DATA
	//------------

	//(2) get Pulser Information from the normal production Run
	TFile *prodfile = TFile::Open(prodname.c_str());
	if (prodfile->IsZombie()){
		cerr<<"Error opening "<<prodname<<"\n";
		return 1;
	}
	TTree* prod_bxtree = (TTree*)prodfile->Get("bxtree");
	if (!prod_bxtree){
		cerr<<"Error getting BxTree \n";
		return 1;
	}

	int n_mchs = 256;	//first mch = 0 (=3001 in db); last mch = 255 (=3256 in db)
	double max_pulsercharge = 5;
	TCanvas *c0_muon=new TCanvas("c0_muon","c0_muon");
	TH2D *Hist_Mch_vs_PulserCharge = new TH2D("Hist_Mch_vs_PulserCharge","MCh vs PulserCharge",n_mchs,3001,n_mchs+3001,100,0,max_pulsercharge);	//Shift by 3001 to convert Echidna to database MCh notation
	prod_bxtree->Draw("muon.decoded_hits.charge : muon.decoded_hits.mch + 3001 >> Hist_Mch_vs_PulserCharge","trigger.trgtype == 32");


	//(3) save the Pulser Plots to jpg-files
 	for (int i=0; i<n_broken_mch;i++){

 		mch = broken_muon_channels_vec[i];

		istomuon1 << "ChargePulserMCh" << mch;
		istomuonname1=istomuon1.str();

		figmuon1 << "/home/production/BadChannelsPlots/Run0" << run << "_MCh" << mch << "_1.jpg";
		figmuonname1=figmuon1.str();

 		TH1D* ChargePulser_MCh = (TH1D*) Hist_Mch_vs_PulserCharge->ProjectionY(istomuonname1.c_str(),mch-3001+1,mch-3001+1);	//Shift of 3001 because of Echdina<->database; +1 because of binning (bin 0 is underflow bin)
		ChargePulser_MCh->SetTitle(istomuonname1.c_str());
		TCanvas *c1_muon=new TCanvas("c1_muon","c1_muon");
		ChargePulser_MCh->Draw();
		c1_muon->SaveAs(figmuonname1.c_str());

		istomuon1.str( "" );
		figmuon1.str( "" );
		ChargePulser_MCh->Delete(); c1_muon->Clear();
	}

	Hist_Mch_vs_PulserCharge->Delete();
	c0_muon->Clear();




	//------------
	//RANDOM DATA
	//------------

	//(4) get Random Information from the normal production Run
	double max_randomcharge = 5;
	TCanvas *c2_muon=new TCanvas("c2_muon","c2_muon");
	TH2D *Hist_Mch_vs_RandomCharge = new TH2D("Hist_Mch_vs_RandomCharge","MCh vs RandomCharge",n_mchs,3001,n_mchs+3001,100,0,max_randomcharge);	//Shift by 3001 to convert Echidna to database MCh notation
	prod_bxtree->Draw("muon.decoded_hits.charge : muon.decoded_hits.mch + 3001 >> Hist_Mch_vs_RandomCharge","trigger.trgtype == 1 && laben.n_clusters != 0 && trigger.btb_inputs != 4");


	//(5) save the Random Plots to jpg-files
 	for (int i=0; i<n_broken_mch;i++){

 		mch = broken_muon_channels_vec[i];

		istomuon2 << "ChargeRandomMCh" << mch;
		istomuonname2=istomuon2.str();

		figmuon2 << "/home/production/BadChannelsPlots/Run0" << run << "_MCh" << mch << "_2.jpg";
		figmuonname2=figmuon2.str();

 		TH1D* ChargeRandom_MCh = (TH1D*) Hist_Mch_vs_RandomCharge->ProjectionY(istomuonname2.c_str(),mch-3001+1,mch-3001+1);	//Shift of 3001 because of Echdina<->database; +1 because of binning (bin 0 is underflow bin)
		ChargeRandom_MCh->SetTitle(istomuonname2.c_str());
		TCanvas *c3_muon=new TCanvas("c3_muon","c3_muon");
		ChargeRandom_MCh->Draw();
		c3_muon->SaveAs(figmuonname2.c_str());

		istomuon2.str( "" );
		figmuon2.str( "" );
		ChargeRandom_MCh->Delete(); c3_muon->Clear();
	}

	Hist_Mch_vs_RandomCharge->Delete();
	c2_muon->Clear();

	prod_bxtree->Delete();
	prodfile->Close();






	//---------
	//LED DATA
	//---------

	//(6) Get the Timing Information from the last Laser-Calibration Run

	//if plots are not present, generate them
	TFile *calibfile = TFile::Open(calib_path);
	if (calibfile->IsZombie()){
		cerr<<"Error opening "<<calib_path<<"\n";
		return 1;
	}
	TTree* calib_bxtree = (TTree*)calibfile->Get("bxtree");
	if (!calib_bxtree){
		cerr<<"Error getting BxTree \n";
		return 1;
	}	

	int max_ledcharge = 150;
	TCanvas *c4_muon=new TCanvas("c4_muon","c4_muon");
	TH2D *Hist_Mch_vs_LEDCharge = new TH2D("Hist_Mch_vs_LEDCharge","MCh vs LEDCharge",n_mchs,3001,n_mchs+3001,max_ledcharge,0,max_ledcharge);//Shift by 3001 to convert Echidna to database MCh notation
	calib_bxtree->Draw("muon.decoded_hits.charge : muon.decoded_hits.mch + 3001 >> Hist_Mch_vs_LEDCharge","trigger.trgtype == 8 && abs(muon.decoded_hits.time - 717) < 40");

	double min_ledtiming = -8000;
	double max_ledtiming = 13000;
	TH2D *Hist_Mch_vs_LEDTime = new TH2D("Hist_Mch_vs_LEDTime","MCh vs LEDTime",n_mchs,3001,n_mchs+3001,1000,min_ledtiming,max_ledtiming);//Shift by 3001 to convert Echidna to database MCh notation
	calib_bxtree->Draw("muon.decoded_hits.time : muon.decoded_hits.mch + 3001 >> Hist_Mch_vs_LEDTime","trigger.trgtype == 8");


	//(7) save the LED Plots to jpg-files
 	for (int i=0; i<n_broken_mch;i++){

 		mch = broken_muon_channels_vec[i];

		istomuon3 << "ChargeLEDMCh" << mch << " , Calibration Run " << calibration_run;
		istomuon4 << "TimingLEDMCh" << mch << " , Calibration Run " << calibration_run;
		istomuonname3=istomuon3.str();
		istomuonname4=istomuon4.str();

		figmuon3 << "/home/production/BadChannelsPlots/CalibRun0" << calibration_run << "_MCh" << mch << "_1.jpg";
		figmuon4 << "/home/production/BadChannelsPlots/CalibRun0" << calibration_run << "_MCh" << mch << "_2.jpg";
		figmuonname3=figmuon3.str();
		figmuonname4=figmuon4.str();

		//Check if file already exists
		ifstream plot_file3(figmuonname3.c_str());
		ifstream plot_file4(figmuonname3.c_str());
		bool files_existent=false;
		if (plot_file3 && plot_file4) files_existent=true;
		plot_file3.close();plot_file4.close();

		if (!files_existent){	//if not existent => create

			TH1D* ChargeLED_MCh = (TH1D*) Hist_Mch_vs_LEDCharge->ProjectionY(istomuonname3.c_str(),mch-3001+1,mch-3001+1);//Shift of 3001 because of Echdina<->database; +1 because of binning (bin 0 is underflow bin)
			ChargeLED_MCh->SetTitle(istomuonname3.c_str());
			TCanvas *c5_muon=new TCanvas("c5_muon","c5_muon");
			ChargeLED_MCh->Draw();
			c5_muon->SaveAs(figmuonname3.c_str());
	
			TH1D* TimeLED_MCh = (TH1D*) Hist_Mch_vs_LEDTime->ProjectionY(istomuonname4.c_str(),mch-3001+1,mch-3001+1);//Shift of 3001 because of Echdina<->database; +1 because of binning (bin 0 is underflow bin)
			TimeLED_MCh->SetTitle(istomuonname4.c_str());
			TCanvas *c6_muon=new TCanvas("c6_muon","c6_muon");
			if (TimeLED_MCh->GetEntries() != 0) c6_muon->SetLogy();
			TimeLED_MCh->Draw();
			c6_muon->SaveAs(figmuonname4.c_str());

			ChargeLED_MCh->Delete(); c5_muon->Clear();
			TimeLED_MCh->Delete(); c6_muon->Clear();
		}
		else {
			std::cout << " --- Calibration plots for MCh " << mch << " via Calibration Run " << calibration_run << " already existent. Performing only linking. --- "<< std::endl;
		}

		//now make links
		figmuon3_link << "/home/production/BadChannelsPlots/Run0" << run << "_MCh" << mch << "_3.jpg";
		figmuon4_link << "/home/production/BadChannelsPlots/Run0" << run << "_MCh" << mch << "_4.jpg";
		figmuonname3_link=figmuon3_link.str();
		figmuonname4_link=figmuon4_link.str();
		symlink (figmuonname3.c_str(), figmuonname3_link.c_str());
		symlink (figmuonname4.c_str(), figmuonname4_link.c_str());

		istomuon3.str( "" );
		istomuon4.str( "" );
		figmuon3.str( "" );
		figmuon4.str( "" );
		figmuon3_link.str( "" );
		figmuon4_link.str( "" );
 	}

	Hist_Mch_vs_LEDCharge->Delete(); Hist_Mch_vs_LEDTime->Delete();
	c4_muon->Clear();
	calib_bxtree->Delete();
	calibfile->Close();

	tab.str( "" );
	tabmuon.str( "" );
	rfpre.str( "" );
	rfec.str( "" );
	prod.str( "" );

	std::cout<<"NORMAL END"<<std::endl;
	return 0;
}
