#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cmath>
#include <vector>

#include <TFile.h>
#include <TDirectory.h>
#include <TCanvas.h>
#include <TTree.h>
#include <TH1D.h>
#include <TF1.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TVector3.h>
#include <TString.h>
#include <TObjString.h>
#include <TObjArray.h>
#include <TKey.h>
#include <TIterator.h>
#include "event/BxEvent.hh"

using namespace std;

void set_branch_status(TTree*, bool);
std::string make_string(const std::string&, Float_t);
std::string rootfilename(const std::string&, Int_t, bool);
void draw_test(TDirectory* dir, TH1D&, TH1D*, const std::string&, const std::string&, TLegend&);
void hstyle(TH1D&, Int_t, bool);
void set_bxstyle();

struct channel_{
	UShort_t lg;
	bool conc;     // light concentrator, (populary said "cone")
	TVector3 pos;  // PMT position
};

int main(int argc, char** argv){
	set_bxstyle();
	bool data = false;

	for (int i=1; i<argc; i++){
		const std::string stmp = argv[i];
		const size_t sep = stmp.find('=');
		const std::string key   = stmp.substr(0,sep);
		const std::string value = stmp.substr(sep+1,stmp.length());
		if (key=="data")
			data = atoi(value.c_str());
		else {
			std::cerr<<"WARNING unknown option: "<<stmp<<'\n';
			exit(255);
		}
	}

	// read input geometry file
	const std::string ifilegeo_name = "tools/mc_validation/geometry_profile17.dat";
	ifstream ifilegeo(ifilegeo_name.c_str());
	if (!ifilegeo){
		std::cerr<<"ERROR: I can't open "<<ifilegeo_name<<std::endl;
		exit(-1);
	}
	const int nchannels = 2241;
	channel_ channel[nchannels];
	while (!ifilegeo.eof()){
		UShort_t lg;
		bool conc;
		Float_t x, y, z;
		ifilegeo>>lg>>conc>>x>>y>>z;
		if (ifilegeo.eof()) break;
		if (lg>= nchannels){
			std::cerr<<"Warning: logical channel > 2240\n";
			continue;
		}
		TVector3 pos(x, y, z);
		channel_ tmpchannel;
		tmpchannel.lg = lg;
		tmpchannel.conc = conc;
		tmpchannel.pos = pos;
		channel[lg] = tmpchannel;
	}
	ifilegeo.close();

	// read input calibration source file
	const std::string ifilesource_name = "tools/mc_validation/calibration_sources.dat";
	ifstream ifilesource(ifilesource_name.c_str());
	if (!ifilesource){
		std::cerr<<"ERROR: I can't open "<<ifilesource_name<<std::endl;
		exit(-1);
	}
	std::vector<std::string> source_id;
	std::vector<std::string> source_name;
	std::vector<int> source_run;
	std::vector<Float_t> source_x;
	std::vector<Float_t> source_y;
	std::vector<Float_t> source_z;
	std::vector<Float_t> source_ene;
	std::vector<Int_t> source_pmt_cut_low;
	std::vector<Int_t> source_pmt_cut_high;
	while (!ifilesource.eof()){
		string name;
		int run, cut_low, cut_high;
		Float_t x, y, z, ene;
		ifilesource >> name >> run >> ene >> x >> y >> z >> cut_low >> cut_high;
		if (ifilesource.eof()) break;
		source_name.push_back(name);
		source_run.push_back(run);
		source_x.push_back(x);
		source_y.push_back(y);
		source_z.push_back(z);
		source_ene.push_back(ene);
		source_pmt_cut_low.push_back(cut_low);
		source_pmt_cut_high.push_back(cut_high);
		ostringstream oss;
		oss<<name<<"Run0"<<run;
		source_id.push_back(oss.str());
	}
	ifilesource.close();
	const unsigned nfiles = source_id.size();

	// prepare output file
	const std::string outfile_name = (data) ? "calibration_sources.root" : "mc_validation.root";
	TFile outfile(outfile_name.c_str(), "RECREATE");
	if (outfile.IsZombie()){
		std::cerr<<"ERROR: I can't open "<<outfile_name<<std::endl;
		exit(-1);
	}

	TFile* rootfile_data;
	const std::string rootfile_data_name = "tools/mc_validation/calibration_sources.root";
	if (!data){
		// open data calibration source rootfile
		rootfile_data = new TFile(rootfile_data_name.c_str(), "READ");
		if (rootfile_data->IsZombie()){
			std::cerr<<"ERROR: I can't open "<<rootfile_data_name<<std::endl;
			exit(-1);
		}
	}


	// main loop on simulation root files
	for (unsigned i=0; i<nfiles; i++){

		// prepare output directory and histograms
		TDirectory* this_dir = outfile.mkdir(source_id[i].c_str());
		this_dir->cd();

		const string se = make_string("Calibration Source Energy =", source_ene[i]);
		TObjString SE(se.c_str());
		SE.Write();
		const string sx = make_string("Calibration Source Position X =", source_x[i]);
		TObjString SX(sx.c_str());
		SX.Write();
		const string sy = make_string("Calibration Source Position Y =", source_y[i]);
		TObjString SY(sy.c_str());
		SY.Write();
		const string sz = make_string("Calibration Source Position Z =", source_z[i]);
		TObjString SZ(sz.c_str());
		SZ.Write();
		TVector3 source_pos(source_x[i], source_y[i], source_z[i]);

		TDirectory* this_dir_data;
		if (!data){
			// copy all histograms from calibration_source.root (data rootfile)
			// in the folder source_id[i]/data in mc_validation.root (mc rootfile)
			TDirectory* data_dir_data = rootfile_data->GetDirectory(source_id[i].c_str());
			this_dir_data = this_dir->mkdir("data");
			this_dir_data->cd();
			TKey* key;
			TIter nextkey(data_dir_data->GetListOfKeys());
			while ((key = (TKey*) nextkey())){
				TObject* obj = key->ReadObj();
				if ((obj->IsA())->InheritsFrom("TH1")){
					((TH1D*) obj)->Write();
				}
			}
		}

		this_dir->cd();
		const Int_t tclus = 1600;
		TH1D h_nnpmts("nnpmts", "NNPMTS; nnpmts", 2000, 0., 2000.);
		TH1D h_nnhits("nnhits", "NNHITS; nnhits", 2000, 0., 2000.);
		TH1D h_ncharge("ncharge", "NCHARGE; ncharge", 2000, 0., 2000.);
		TH1D h_hit_times("hit_times", "Clustered Hit Time all Lg; time (ns)", tclus, 0., Double_t(tclus));
		TH1D h_rec_times("rec_times", "Rec Time all Lg; time (ns)", tclus, 0., Double_t(tclus));
		TH1D h_hit_times_conc("hit_times_conc", "Clustered Hit Time, Channels with Concentrators; time (ns)", tclus, 0., Double_t(tclus));
		TH1D h_hit_times_noconc("hit_times_noconc", "Clustered Hit Time, Channels with NO Concentrators; time (ns)", tclus, 0., Double_t(tclus));
		TH1D h_hit_times_1stchhit("hit_times_1stchhit", "Clustered Hit Time, First Channel Hit only; time (ns)", tclus, 0., Double_t(tclus));
		TH1D h_hit_times_no1stchhit("hit_times_no1stchhit", "Clustered Hit Time, all Hits except First Channel Hit; time (ns)", tclus, 0., Double_t(tclus));
		TH1D h_hit_times_near("hit_times_near", "Clustered Hit Time, PMT-Source distance < 4 m; time (ns)", tclus, 0., Double_t(tclus));
		TH1D h_hit_times_far("hit_times_far", "Clustered Hit Time, PMT-Source distance > 8 m; time (ns)", tclus, 0., Double_t(tclus));
		TH1D h_hit_times_near_1stchhit("hit_times_near_1stchhit", "Clustered Hit Time, PMT-Source distance < 4 m, First Channel Hit Only; time (ns)", tclus, 0., Double_t(tclus));
		TH1D h_hit_times_far_1stchhit("hit_times_far_1stchhit", "Clustered Hit Time, PMT-Source distance > 8 m, First Channel Hit Only; time (ns)", tclus, 0., Double_t(tclus));
		TH1D h_nnpmts_near("nnpmts_near", "NNPMTS, PMT-Source distance < 4m; nnpmts", 2000, 0., 2000.);
		TH1D h_nnpmts_far("nnpmts_far", "NNPMTS, PMT-Source distance > 8m; nnpmts", 2000, 0., 2000.);
		TH1D h_nnpmts_conc("nnpmts_conc", "NNPMTS, Channels with Concentrators; nnpmts", 2000, 0., 2000.);
		TH1D h_nnpmts_noconc("nnpmts_noconc", "NNPMTS, Channels with NO Concentrators; nnpmts", 2000, 0., 2000.);
		TH1D h_ncharge_pmts("ncharge_pmts", "NCHARGEPMTS (Charge of NPMTS); ncharge_pmts", 2000, 0., 2000.);
		TH1D h_raw_charge("raw_charge", "Peak - Base (without pileup corr.); raw_charge;", 255, 0, 255);
		TH1D h_raw_charge_1stchhit("raw_charge_1stchhit", "Peak - Base, First Channel Hit only", 255, 0 , 255);
		TH1D h_raw_charge_no1stchhit("raw_charge_no1stchhit", "Peak - Base, all Hits except First Channel Hit", 255, 0 , 255);
		TH1D h_charge_1stchhit("charge_1stchhit", "Clustered Hit Charge, First Channel Hit only", 100, 0 , 20);
		TH1D h_charge_no1stchhit("charge_no1stchhit", "Clustered Hit Charge, all Hits except First Channel Hit", 100, 0 , 20);
		TH1D h_x("x", "PositionLNGS X", 1200, -6., 6.);
		TH1D h_y("y", "PositionLNGS Y", 1200, -6., 6.);
		TH1D h_z("z", "PositionLNGS Z", 1200, -6., 6.);
		TH1D h_r0("r0", "PositionLNGS Source-Event distance", 50, 0, 0.5);
		TH1D h_gatti("gatti", "Gatti distribution", 200, -0.2, 0.2);
		TH1D h_peak_time("peak_time", "Peak Time", 300, 0, 300);
		TH1D h_mean_time_short("mean_time_short", "Mean Time (short cluster)", 500, 0, 500);
		TH1D h_rms_time_short("rms_time_short", "RMS Time (short cluster)", 300, 0, 300);
		TH1D h_hits_on_pmts("hits_on_pmts", "NHits/NPmts", 300, 1., 4.);
		TH1D h_charge_on_pmts("charge_on_pmts", "NCharge/NNPmts", 80, 0., 4.);
		TH1D h_charge_on_hits("charge_on_hits", "NCharge/NNHits", 80, 0., 4.);
		TH1D h_conc_on_noconc("conc_on_noconc", "NPmts Concentrators / NPmts NoConcentrators", 500, 0., 50);
		TH1D h_decoded_times("decoded_times", "Decoded Hit Times; time (ns) respect 1st clustered hit;", tclus + 500 , -500., Double_t(tclus));
		TH1D h_decoded_times_order1("decoded_times_order1", "Decoded Hit Times, order 1; time (ns) respect 1st clustered hit;", tclus+500, -500., Double_t(tclus));
		TH1D h_decoded_times_order2("decoded_times_order2", "Decoded Hit Times, order 2; time (ns) respect 1st clustered hit;", tclus+500, -500., Double_t(tclus));
		TH1D h_decoded_dt_12("decoded_dt_12", "Time difference between order 2 and 1 decoded hit; time (ns);", tclus, 0., Double_t(tclus));

		const Int_t color = (data) ? kBlack : 4;
		hstyle(h_nnpmts, color, data);
		hstyle(h_nnhits, color, data);
		hstyle(h_ncharge, color, data);
		hstyle(h_hit_times, color, data);
		hstyle(h_rec_times, color, data);
		hstyle(h_hit_times_conc, color, data);
		hstyle(h_hit_times_noconc, color, data);
		hstyle(h_hit_times_1stchhit, color, data);
		hstyle(h_hit_times_no1stchhit, color, data);
		hstyle(h_hit_times_near, color, data);
		hstyle(h_hit_times_far, color, data);
		hstyle(h_hit_times_near_1stchhit, color, data);
		hstyle(h_hit_times_far_1stchhit, color, data);
		hstyle(h_nnpmts_near, color, data);
		hstyle(h_nnpmts_far, color, data);
		hstyle(h_nnpmts_conc, color, data);
		hstyle(h_nnpmts_noconc, color, data);
		hstyle(h_ncharge_pmts, color, data);
		hstyle(h_raw_charge, color, data);
		hstyle(h_raw_charge_1stchhit, color, data);
		hstyle(h_raw_charge_no1stchhit, color, data);
		hstyle(h_charge_1stchhit, color, data);
		hstyle(h_charge_no1stchhit, color, data);
		hstyle(h_x, color, data);
		hstyle(h_y, color, data);
		hstyle(h_z, color, data);
		hstyle(h_r0, color, data);
		hstyle(h_gatti, color, data);
		hstyle(h_peak_time, color, data);
		hstyle(h_mean_time_short, color, data);
		hstyle(h_rms_time_short, color, data);
		hstyle(h_hits_on_pmts, color, data);
		hstyle(h_charge_on_pmts, color, data);
		hstyle(h_charge_on_hits, color, data);
		hstyle(h_conc_on_noconc, color, data);
		hstyle(h_decoded_times, color, data);
		hstyle(h_decoded_times_order1, color, data);
		hstyle(h_decoded_times_order2, color, data);
		hstyle(h_decoded_dt_12, color, data);

		const std::string this_source_rootfile_name = rootfilename(source_id[i], source_run[i], data);
		TFile rootfile(this_source_rootfile_name.c_str(), "READ");
		if (rootfile.IsZombie()){
			std::cerr<<"ERROR: I can't open "<<this_source_rootfile_name<<std::endl;
			exit(i+1);
		}
		cout<<this_source_rootfile_name<<" opened\n";
		
		// get bxtree, loop on it

		TTree* bxtree = (TTree*) rootfile.Get("bxtree");
		if (!bxtree){
			cerr<<"Error: No bxtree in "<<this_source_rootfile_name<<endl;
			continue;
		}

		// address and activate branch 
		BxEvent* ev   = new BxEvent();
		bxtree->SetBranchAddress("events", &ev);
		set_branch_status(bxtree, data);
		
		const Int_t cut_low  = source_pmt_cut_low.at(i);
		const Int_t cut_high = source_pmt_cut_high.at(i);
		Long64_t t_mu = 0;
		// event loop
		const Int_t nentries = bxtree->GetEntries();

		//Run Duration in Second
		bxtree->GetEntry(0);	
		const Long64_t gps_time_start = ev->GetTrigger().GetGpsTimeSec();
		bxtree->GetEntry(nentries-2);	
		const Long64_t gps_time_stop = ev->GetTrigger().GetGpsTimeSec();
		Long64_t duration = gps_time_stop - gps_time_start;

		for (Int_t ie=0; ie<nentries; ie++){
			bxtree->GetEntry(ie);

			const Int_t NClusters  = ev->GetLaben().GetNClusters();
			if (data){
				const Int_t nhits_short = (NClusters>0) ? ev->GetLaben().GetCluster(0).GetNHitsShort() : 0;
				const Int_t trg_type   = ev->GetTrigger().GetTrgType();
				const Int_t btb_inputs = ev->GetTrigger().GetBtbInputs();
				const bool is_muon_aligned = ev->GetMuon().IsAligned();
				const bool muon_has_cluster = ev->GetMuon().HasCluster();
				bool muon_mtb = (btb_inputs==4);
				bool muon_mcr = (is_muon_aligned) ? (muon_has_cluster) : false;
				bool muon_idf = false; // to be implemented ...
				bool muon_internal_strict =  (NClusters>0) ? ( ((trg_type==1) && (muon_mtb || muon_mcr || muon_idf)) || ((trg_type==2) && (nhits_short>80)) ) : false;
				bool muon_internal_special_d5 = ((trg_type==1) && (NClusters==0) && (muon_mtb || muon_mcr));

				bool muon_internal = muon_internal_strict || muon_internal_special_d5; 
				Long64_t gps_time_sec = ev->GetTrigger().GetGpsTimeSec();
				Long64_t gps_time_ns  = ev->GetTrigger().GetGpsTimeNs();
				Long64_t t_ev = gps_time_sec + 1000000000* gps_time_ns;

				if (muon_internal){
					t_mu = t_ev;
					continue;
				}
				Long_t dt_mu = t_ev - t_mu;
				if (dt_mu < 10000000){ // 10 ms
					//	cout<<"Tagged Muon Daughter!\n";
					continue;
				}
				if ((trg_type!=1) && (btb_inputs!=0)) continue;

			}

			if (NClusters!=1) continue;

			const Int_t NPmts     = ev->GetLaben().GetCluster(0).GetNPmts();
			const Int_t NLivePmts   = ev->GetLaben().GetNLivePmts();
			const Float_t NNPmts  = 2000./NLivePmts*NPmts;
			if ((NNPmts<cut_low) || (NNPmts>cut_high)) continue;

			const Float_t x = ev->GetLaben().GetCluster(0).GetPositionLNGS().GetX();
			const Float_t y = ev->GetLaben().GetCluster(0).GetPositionLNGS().GetY();
			const Float_t z = ev->GetLaben().GetCluster(0).GetPositionLNGS().GetZ();
			TVector3 ev_pos(x, y, z);
			if ((ev_pos - source_pos).Mag()>0.5) continue;

			h_x.Fill(x);
			h_y.Fill(y);
			h_z.Fill(z);
			h_r0.Fill((ev_pos-source_pos).Mag());

			const Int_t NLiveCharge = ev->GetLaben().GetNLiveCharge();
			const Int_t NHits       = ev->GetLaben().GetCluster(0).GetNHits();
			const Float_t Charge      = ev->GetLaben().GetCluster(0).GetCharge();
			const Float_t ChargePmts  = ev->GetLaben().GetCluster(0).GetChargeNpmts();
			const Float_t NNHits  = 2000./NLivePmts*NHits;
			const Float_t NCharge = 2000./NLiveCharge*Charge;
			const Float_t NChargePmts = 2000./NLiveCharge*ChargePmts;

			h_nnpmts.Fill(NNPmts);
			h_nnhits.Fill(NNHits);
			h_ncharge.Fill(NCharge);
			h_ncharge_pmts.Fill(NChargePmts);
			h_hits_on_pmts.Fill(NNHits/NNPmts);
			h_charge_on_pmts.Fill(NCharge/NNPmts);
			h_charge_on_hits.Fill(NCharge/NNHits);

			const Int_t npeaks = ev->GetLaben().GetCluster(0).GetNPeaks();
			const Float_t gatti = ev->GetLaben().GetRecCluster(0).GetGatti();
			const Float_t peak_time = (npeaks==1) ? (ev->GetLaben().GetCluster(0).GetPeakTimes()[0]) : 0;
			const Float_t mean_time_short = ev->GetLaben().GetCluster(0).GetMeanTimeShort();
			const Float_t rms_time_short = ev->GetLaben().GetCluster(0).GetRMSTimeShort();

			h_gatti.Fill(gatti);
			h_peak_time.Fill(peak_time);
			h_mean_time_short.Fill(mean_time_short);
			h_rms_time_short.Fill(rms_time_short);

			// time distribution
			const Int_t NDeco = ev->GetLaben().GetNDecodedHits();
			static std::vector<bool> ch_nohit(nchannels);  	// true = this channel has recieved no hit before
			ch_nohit.assign(nchannels, true);		// no channel has recieved a hit at this level
			Double_t ch_t_1[nchannels];			// 1st decoded hit time of each channel 

			// find first decoded hit time time 
			Double_t first_hit_time = 0.;
			for (int d=0; d<NDeco; d++){
				if (ev->GetLaben().GetDecodedHit(d).GetNumCluster() == 1){
					first_hit_time = ev->GetLaben().GetDecodedHit(d).GetRawTime();
					break;
				}
			}

			Int_t NPmts_near = 0, NPmts_far = 0, NPmts_conc = 0, NPmts_noconc = 0;
			// fill time distrib of clustered hits histos
			for (int d = 0; d < NDeco; d++){
				const UShort_t lg = ev->GetLaben().GetDecodedHit(d).GetLg();
				if (lg>=nchannels){
					std::cout<<"Warning: logical channel > 2240 \n";
					continue;
				}
				const Double_t hit_time = ev->GetLaben().GetDecodedHit(d).GetRawTime() - first_hit_time;
				const Int_t order = ev->GetLaben().GetDecodedHit(d).GetOrder(); 
				h_decoded_times.Fill(hit_time);
				if (order==1) {
					ch_t_1[lg] = hit_time;
					h_decoded_times_order1.Fill(hit_time);
					}
				if (order==2) {
					h_decoded_times_order2.Fill(hit_time);
					h_decoded_dt_12.Fill(hit_time - ch_t_1[lg]);
					}

				if(ev->GetLaben().GetDecodedHit(d).GetNumCluster() == 1){
				//clustered hits
					const bool has_conc = (channel[lg].conc);
					const bool has_nohit   = (ch_nohit[lg]);
					const Float_t pmt_source_distance = (source_pos - channel[lg].pos).Mag();
					const bool is_near = ((pmt_source_distance<4.));
					const bool is_far  = ((pmt_source_distance>8.));
					const Double_t rec_time = ev->GetLaben().GetDecodedHit(d).GetRecTime();
					const Float_t raw_charge = ev->GetLaben().GetDecodedHit(d).GetRawCharge();
					const Float_t charge = ev->GetLaben().GetDecodedHit(d).GetCharge();
					h_hit_times.Fill(hit_time);
					h_rec_times.Fill(rec_time);
					h_raw_charge.Fill(raw_charge);
					if (has_nohit){
						h_hit_times_1stchhit.Fill(hit_time);
						h_raw_charge_1stchhit.Fill(raw_charge);
						h_charge_1stchhit.Fill(charge);
						if (is_near) {
							h_hit_times_near_1stchhit.Fill(hit_time);
							NPmts_near++;
							}
						if (is_far) {
							h_hit_times_far_1stchhit.Fill(hit_time);
							NPmts_far++;
							}
						if (has_conc){
							NPmts_conc++;
						}
						if (!has_conc){
							NPmts_noconc++;
						}
						ch_nohit[lg] = false;
					}
					else if (!has_nohit){
						h_hit_times_no1stchhit.Fill(hit_time);
						h_raw_charge_no1stchhit.Fill(raw_charge);
						h_charge_no1stchhit.Fill(charge);
					} //
					if (is_near){
						h_hit_times_near.Fill(hit_time);
					}
					if (is_far){
						h_hit_times_far.Fill(hit_time);
					} //
					if (has_conc){
						h_hit_times_conc.Fill(hit_time);
					}
					else if (!has_conc){
						h_hit_times_noconc.Fill(hit_time);
					} //
				} 
			} // end decoded hits loop
			
			const Float_t NNPmts_near = NPmts_near*2000./NLivePmts;
			const Float_t NNPmts_far  = NPmts_far*2000./NLivePmts;
			if (NNPmts_near>0) h_nnpmts_near.Fill(NNPmts_near);
			if (NNPmts_far>0)  h_nnpmts_far.Fill(NNPmts_far);
			const Float_t NNPmts_conc   = NPmts_conc*2000./NLivePmts;
			const Float_t NNPmts_noconc = NPmts_noconc*2000./NLivePmts;
			h_nnpmts_conc.Fill(NNPmts_conc);
			h_nnpmts_noconc.Fill(NNPmts_noconc);
			if (NNPmts_noconc>0) h_conc_on_noconc.Fill(NNPmts_conc/NNPmts_noconc); 

		} // end bxtree event loop 
		rootfile.Close();
		outfile.cd(source_id[i].c_str());
		const string sd = make_string("Duration=", duration);
		TObjString SD(sd.c_str());
		SD.Write();
		h_nnpmts.Write();
		h_nnhits.Write();
		h_ncharge.Write();
		h_hit_times.Write();
		h_rec_times.Write();
		h_hit_times_1stchhit.Write();
		h_hit_times_no1stchhit.Write();
		h_hit_times_conc.Write();
		h_hit_times_noconc.Write();
		h_hit_times_near.Write();
		h_hit_times_far.Write();
		h_hit_times_near_1stchhit.Write();
		h_hit_times_far_1stchhit.Write();
		h_nnpmts_near.Write();
		h_nnpmts_far.Write();
		h_nnpmts_conc.Write();
		h_nnpmts_noconc.Write();
		h_ncharge_pmts.Write();
		h_raw_charge.Write();
		h_raw_charge_1stchhit.Write();
		h_raw_charge_no1stchhit.Write();
		h_charge_1stchhit.Write();
		h_charge_no1stchhit.Write();
		h_x.Write();
		h_y.Write();
		h_z.Write();
		h_r0.Write();
		h_gatti.Write();
		h_peak_time.Write();
		h_mean_time_short.Write();
		h_rms_time_short.Write();
		h_hits_on_pmts.Write();
		h_charge_on_pmts.Write();
		h_charge_on_hits.Write();
		h_conc_on_noconc.Write();
		h_decoded_times.Write();
		h_decoded_times_order1.Write();
		h_decoded_times_order2.Write();
		h_decoded_dt_12.Write();


		// canvas composition here:
		// c1 -> general (nnpmts, nnhits, hit_times, rec_times)
		// c2 -> near_far (nnpmts_near, nnpmts_far, hit_times_near, hit_times_far)
		// c3 -> conc_noconc (nnpmts_conc, nnpmts_noconc, hit_times_conc, hit_times_noconc)
		// c4 -> first_hit (hit_times_1stch, hit_times_no1stch, hit_times_near_1stch, hit_times_far_1stch)
		// c5 -> charge (ncharge, ncharge_pmts, raw_charge, raw_charge_1stchhit);
		// c6 -> ratios (hits_on_pmts, charge_on_hits, charge_on_pmts, conc_on_noconc);
		// c7 -> times (gatti, peak_time, mean_time_short, rms_time_short)
		// c8 -> position (r0, z, x, y);
		// c9 -> decoded (decoded_times, decoded_dt_12, decoded_times_order1, decoded_times_order2);
		//
		if (!data){
			TDirectory* this_dir_canvas = this_dir->mkdir("canvas");
			// c1
			this_dir_canvas->cd();
			ostringstream oss_title;
			oss_title<<"General comparison histograms for "<<source_id[i];
			TCanvas c1("c1_general", oss_title.str().c_str());
			c1.Divide(2,2);
			c1.cd(1);
			TLegend leg_1_1(0.65, 0.90, 0.98, 0.55);
			TH1D* data_nnpmts;
			draw_test(this_dir_data, h_nnpmts, data_nnpmts, "nnpmts", source_id[i], leg_1_1);
			c1.cd(2);
			TLegend leg_1_2(0.65, 0.90, 0.98, 0.55);
			TH1D* data_nnhits;
			draw_test(this_dir_data, h_nnhits, data_nnhits, "nnhits", source_id[i], leg_1_2);
			TVirtualPad* c1_3 = c1.cd(3);
			c1_3->SetLogy();
			TLegend leg_1_3(0.65, 0.90, 0.98, 0.55);
			TH1D* data_hit_times;
			draw_test(this_dir_data, h_hit_times, data_hit_times, "hit_times", source_id[i], leg_1_3);
			c1.cd(4);
			TLegend leg_1_4(0.65, 0.90, 0.98, 0.55);
			TH1D* data_rec_times;
			draw_test(this_dir_data, h_rec_times, data_rec_times, "rec_times", source_id[i], leg_1_4);
			c1.cd(0);
			this_dir_canvas->cd();
			c1.Write();

			// c2
			this_dir_canvas->cd();
			ostringstream oss_title2;
			oss_title2<<"Near/far PMTS comparison histograms for "<<source_id[i];
			TCanvas c2("c2_near_far", oss_title2.str().c_str());
			c2.Divide(2,2);
			c2.cd(1);
			TLegend leg_2_1(0.65, 0.90, 0.98, 0.55);
			TH1D* data_nnpmts_near;
			draw_test(this_dir_data, h_nnpmts_near, data_nnpmts_near, "nnpmts_near", source_id[i], leg_2_1);
			c2.cd(2);
			TLegend leg_2_2(0.65, 0.90, 0.98, 0.55);
			TH1D* data_nnpmts_far;
			draw_test(this_dir_data, h_nnpmts_far, data_nnpmts_far, "nnpmts_far", source_id[i], leg_2_2);
			TVirtualPad* c2_3 = c2.cd(3);
			c2_3->SetLogy();
			TLegend leg_2_3(0.65, 0.90, 0.98, 0.55);
			TH1D* data_hit_times_near;
			draw_test(this_dir_data, h_hit_times_near, data_hit_times_near, "hit_times_near", source_id[i], leg_2_3);
			TVirtualPad* c2_4 = c2.cd(4);
			c2_4->SetLogy();
			TLegend leg_2_4(0.65, 0.90, 0.98, 0.55);
			TH1D* data_hit_times_far;
			draw_test(this_dir_data, h_hit_times_far, data_hit_times_far, "hit_times_far", source_id[i], leg_2_4);
			c2.cd();
			this_dir_canvas->cd();
			c2.Write();

			// c3
			this_dir_canvas->cd();
			ostringstream oss_title3;
			oss_title3<<"Channels with/without concentrators comparison histograms for "<<source_id[i];
			TCanvas c3("c3_conc_noconc", oss_title3.str().c_str());
			c3.Divide(2,2);
			c3.cd(1);
			TLegend leg_3_1(0.65, 0.90, 0.98, 0.55);
			TH1D* data_nnpmts_conc;
			draw_test(this_dir_data, h_nnpmts_conc, data_nnpmts_conc, "nnpmts_conc", source_id[i], leg_3_1);
			c3.cd(2);
			TLegend leg_3_2(0.65, 0.90, 0.98, 0.55);
			TH1D* data_nnpmts_noconc;
			draw_test(this_dir_data, h_nnpmts_noconc, data_nnpmts_noconc, "nnpmts_noconc", source_id[i], leg_3_2);
			TVirtualPad* c3_3 = c3.cd(3);
			c3_3->SetLogy();
			TLegend leg_3_3(0.65, 0.90, 0.98, 0.55);
			TH1D* data_hit_times_conc;
			draw_test(this_dir_data, h_hit_times_conc, data_hit_times_conc, "hit_times_conc", source_id[i], leg_3_3);
			TVirtualPad* c3_4 = c3.cd(4);
			c3_4->SetLogy();
			TLegend leg_3_4(0.65, 0.90, 0.98, 0.55);
			TH1D* data_hit_times_noconc;
			draw_test(this_dir_data, h_hit_times_noconc, data_hit_times_noconc, "hit_times_noconc", source_id[i], leg_3_4);
			c3.cd(0);
			this_dir_canvas->cd();
			c3.Write();

			// c4
			this_dir_canvas->cd();
			ostringstream oss_title4;
			oss_title4<<"First/NoFirst PMT hit comparison histograms for "<<source_id[i];
			TCanvas c4("c4_1st_no1st", oss_title4.str().c_str());
			c4.Divide(2,2);
			TVirtualPad* c4_1 = c4.cd(1);
			c4_1->SetLogy();
			TLegend leg_4_1(0.65, 0.90, 0.98, 0.55);
			TH1D* data_hit_times_1stchhit;
			draw_test(this_dir_data, h_hit_times_1stchhit, data_hit_times_1stchhit, "hit_times_1stchhit", source_id[i], leg_4_1);
			TVirtualPad* c4_2 = c4.cd(2);
			c4_2->SetLogy();
			TLegend leg_4_2(0.65, 0.90, 0.98, 0.55);
			TH1D* data_hit_times_no1stchhit;
			draw_test(this_dir_data, h_hit_times_no1stchhit, data_hit_times_no1stchhit, "hit_times_no1stchhit", source_id[i], leg_4_2);
			TVirtualPad* c4_3 = c4.cd(3);
			c4_3->SetLogy();
			TLegend leg_4_3(0.65, 0.90, 0.98, 0.55);
			TH1D* data_hit_times_near_1stchhit;
			draw_test(this_dir_data, h_hit_times_near_1stchhit, data_hit_times_near_1stchhit, "hit_times_near_1stchhit", source_id[i], leg_4_3);
			TVirtualPad* c4_4 = c4.cd(4);
			c4_4->SetLogy();
			TLegend leg_4_4(0.65, 0.90, 0.98, 0.55);
			TH1D* data_hit_times_far_1stchhit;
			draw_test(this_dir_data, h_hit_times_far_1stchhit, data_hit_times_far_1stchhit, "hit_times_far_1stchhit", source_id[i], leg_4_4);
			c4.cd(0);
			this_dir_canvas->cd();
			c4.Write();

			// c5
			this_dir_canvas->cd();
			ostringstream oss_title5;
			oss_title5<<"Charge comparison histograms for "<<source_id[i];
			TCanvas c5("c5_charge", oss_title5.str().c_str());
			c5.Divide(2,2);
			c5.cd(1);
			TH1D* data_ncharge;
			TLegend leg_5_1(0.65, 0.90, 0.98, 0.55);
			draw_test(this_dir_data, h_ncharge, data_ncharge, "ncharge", source_id[i], leg_5_1);
			c5.cd(2);
			TH1D* data_ncharge_pmts;
			TLegend leg_5_2(0.65, 0.90, 0.98, 0.55);
			draw_test(this_dir_data, h_ncharge_pmts, data_ncharge_pmts, "ncharge_pmts", source_id[i], leg_5_2);
			c5.cd(3);
			TH1D* data_raw_charge;
			TLegend leg_5_3(0.65, 0.90, 0.98, 0.55);
			draw_test(this_dir_data, h_raw_charge, data_raw_charge, "raw_charge", source_id[i], leg_5_3);
			//c5.cd(5);
			//TH1D* data_raw_charge_1stchhit;
			//TLegend leg_5_3(0.65, 0.90, 0.98, 0.55);
			//draw_test(this_dir_data, h_raw_charge_1stchhit, data_raw_charge_1stchhit, "raw_charge_1stchhit", source_id[i], leg_5_3);
			//c5.cd(4);
			//TH1D* data_raw_charge_no1stchhit;
			//TLegend leg_5_4(0.65, 0.90, 0.98, 0.55);
			//draw_test(this_dir_data, h_raw_charge_no1stchhit, data_raw_charge_no1stchhit, "raw_charge_no1stchhit", source_id[i], leg_5_4);
			c5.cd(4);
			TH1D* data_charge_1stchhit;
			TLegend leg_5_4(0.65, 0.90, 0.98, 0.55);
			draw_test(this_dir_data, h_charge_1stchhit, data_charge_1stchhit, "charge_1stchhit", source_id[i], leg_5_4);
			//c5.cd(6);
			//TH1D* data_charge_no1stchhit;
			//TLegend leg_5_6(0.65, 0.90, 0.98, 0.55);
			//draw_test(this_dir_data, h_charge_no1stchhit, data_charge_no1stchhit, "charge_no1stchhit", source_id[i], leg_5_6);

			c5.cd(0);
			this_dir_canvas->cd();
			c5.Write();
		
			//c6
			this_dir_canvas->cd();
			ostringstream oss_title6;
			oss_title6<<"Ratio comparison histograms for "<<source_id[i];
			TCanvas c6("c6_ratio", oss_title6.str().c_str());
			c6.Divide(2,2);
			c6.cd(1);
			TH1D* data_hits_on_pmts;
			TLegend leg_6_1(0.65, 0.90, 0.98, 0.55);
			draw_test(this_dir_data, h_hits_on_pmts, data_hits_on_pmts, "hits_on_pmts", source_id[i], leg_6_1);
			c6.cd(2);
			TH1D* data_charge_on_pmts;
			TLegend leg_6_2(0.65, 0.90, 0.98, 0.55);
			draw_test(this_dir_data, h_charge_on_pmts, data_charge_on_pmts, "charge_on_pmts", source_id[i], leg_6_2);
			c6.cd(3);
			TH1D* data_charge_on_hits;
			TLegend leg_6_3(0.65, 0.90, 0.98, 0.55);
			draw_test(this_dir_data, h_charge_on_hits, data_charge_on_hits, "charge_on_hits", source_id[i], leg_6_3);
			c6.cd(4);
			TH1D* data_conc_on_noconc;
			TLegend leg_6_4(0.65, 0.90, 0.98, 0.55);
			draw_test(this_dir_data, h_conc_on_noconc, data_conc_on_noconc, "conc_on_noconc", source_id[i], leg_6_4);
			c6.cd(0);
			this_dir_canvas->cd();
			c6.Write();


			//c7
			this_dir_canvas->cd();
			ostringstream oss_title7;
			oss_title7<<"Charge comparison histograms for "<<source_id[i];
			TCanvas c7("c7_times", oss_title7.str().c_str());
			c7.Divide(2,2);
			c7.cd(1);
			TH1D* data_gatti;
			TLegend leg_7_1(0.65, 0.90, 0.98, 0.55);
			draw_test(this_dir_data, h_gatti, data_gatti, "gatti", source_id[i], leg_7_1);
			c7.cd(2);
			TH1D* data_peak_time;
			TLegend leg_7_2(0.65, 0.90, 0.98, 0.55);
			draw_test(this_dir_data, h_peak_time, data_peak_time, "peak_time", source_id[i], leg_7_2);
			c7.cd(3);
			TH1D* data_mean_time_short;
			TLegend leg_7_3(0.65, 0.90, 0.98, 0.55);
			draw_test(this_dir_data, h_mean_time_short, data_mean_time_short, "mean_time_short", source_id[i], leg_7_3);
			c7.cd(4);
			TH1D* data_rms_time_short;
			TLegend leg_7_4(0.65, 0.90, 0.98, 0.55);
			draw_test(this_dir_data, h_rms_time_short, data_rms_time_short, "rms_time_short", source_id[i], leg_7_4);
			c7.cd(0);
			this_dir_canvas->cd();
			c7.Write();

			//c8
			this_dir_canvas->cd();
			ostringstream oss_title8;
			oss_title8<<"Position comparison histograms for "<<source_id[i];
			TCanvas c8("c8_position", oss_title8.str().c_str());
			c8.Divide(2,2);
			c8.cd(1);
			TH1D* data_r0;
			TLegend leg_8_1(0.65, 0.90, 0.98, 0.55);
			draw_test(this_dir_data, h_r0, data_r0, "r0", source_id[i], leg_8_1);
			c8.cd(2);
			TH1D* data_z;
			TLegend leg_8_2(0.65, 0.90, 0.98, 0.55);
			draw_test(this_dir_data, h_z, data_z, "z", source_id[i], leg_8_2);
			c8.cd(3);
			TH1D* data_x;
			TLegend leg_8_3(0.65, 0.90, 0.98, 0.55);
			draw_test(this_dir_data, h_x, data_x, "x", source_id[i], leg_8_3);
			c8.cd(4);
			TH1D* data_y;
			TLegend leg_8_4(0.65, 0.90, 0.98, 0.55);
			draw_test(this_dir_data, h_y, data_y, "y", source_id[i], leg_8_4);
			c8.cd(0);
			this_dir_canvas->cd();
			c8.Write();

			// c9
			this_dir_canvas->cd();
			ostringstream oss_title9;
			oss_title9<<"Decoded Hit Times for "<<source_id[i];
			TCanvas c9("c0_decoded", oss_title9.str().c_str());
			c9.Divide(2,2);
			TVirtualPad* c9_1 = c9.cd(1);
			c9_1->SetLogy();
			TLegend leg_9_1(0.65, 0.90, 0.98, 0.55);
			TH1D* data_decoded_times;
			draw_test(this_dir_data, h_decoded_times, data_decoded_times, "decoded_times", source_id[i], leg_9_1);
			TVirtualPad* c9_2 = c9.cd(2);
			c9_2->SetLogy();
			TLegend leg_9_2(0.65, 0.90, 0.98, 0.55);
			TH1D* data_decoded_dt_12;
			draw_test(this_dir_data, h_decoded_dt_12, data_decoded_dt_12, "decoded_dt_12", source_id[i], leg_9_2);
			TVirtualPad* c9_3 = c9.cd(3);
			c9_3->SetLogy();
			TLegend leg_9_3(0.65, 0.90, 0.98, 0.55);
			TH1D* data_decoded_times_order1;
			draw_test(this_dir_data, h_decoded_times_order1, data_decoded_times_order1, "decoded_times_order1", source_id[i], leg_9_3);
			TVirtualPad* c9_4 = c9.cd(4);
			c9_4->SetLogy();
			TLegend leg_9_4(0.65, 0.90, 0.98, 0.55);
			TH1D* data_decoded_times_order2;
			draw_test(this_dir_data, h_decoded_times_order2, data_decoded_times_order2, "decoded_times_order2", source_id[i], leg_9_4);
			c9.cd(0);
			this_dir_canvas->cd();
			c9.Write();



		}

	} // end loop on simulation root files

	std::cout<<"Closing output rootfile "<<outfile_name<<'\n';
	outfile.Close();
	if (!data) {
		rootfile_data->Close();
	}
	std::cout<<"NORMAL END\n";
	return 0;
}

void set_branch_status(TTree* bxtree, bool data=0){
        bxtree->SetBranchStatus("*", 0);
	bxtree->SetBranchStatus("laben.n_decoded_hits", 1);
	bxtree->SetBranchStatus("laben.decoded_hits", 1);
	bxtree->SetBranchStatus("laben.decoded_hits.lg", 1);
	bxtree->SetBranchStatus("laben.decoded_hits.order", 1);
	bxtree->SetBranchStatus("laben.decoded_hits.raw_time", 1);
	bxtree->SetBranchStatus("laben.decoded_hits.rec_time", 1);
	bxtree->SetBranchStatus("laben.decoded_hits.num_cluster", 1);
	bxtree->SetBranchStatus("laben.decoded_hits.raw_charge", 1);
	bxtree->SetBranchStatus("laben.decoded_hits.charge", 1);
	bxtree->SetBranchStatus("laben.clusters.npmts", 1);
	bxtree->SetBranchStatus("laben.clusters.nhits", 1);
	bxtree->SetBranchStatus("laben.clusters.charge", 1);
	bxtree->SetBranchStatus("laben.clusters.charge_npmts", 1);
	bxtree->SetBranchStatus("laben.n_clusters", 1);
	bxtree->SetBranchStatus("laben.n_live_pmts", 1);
	bxtree->SetBranchStatus("laben.n_live_charge", 1);
	bxtree->SetBranchStatus("laben.rec_clusters.gatti");
	bxtree->SetBranchStatus("laben.clusters.peak_times");
	bxtree->SetBranchStatus("laben.clusters.mean_time_short");
	bxtree->SetBranchStatus("laben.clusters.rms_time_short");
	bxtree->SetBranchStatus("laben.clusters.position_lngs.x",1);
	bxtree->SetBranchStatus("laben.clusters.position_lngs.y",1);
	bxtree->SetBranchStatus("laben.clusters.position_lngs.z",1);
	if (data){
		// data rootfiles
		bxtree->SetBranchStatus("trigger.trgtype",1);
		bxtree->SetBranchStatus("trigger.gpstimes[2]",1);
		bxtree->SetBranchStatus("trigger.btb_inputs",1);
		bxtree->SetBranchStatus("muon.is_aligned", 1);
		bxtree->SetBranchStatus("muon.has_clustered", 1);
		bxtree->SetBranchStatus("muon.has_cluster_sss", 1);
		bxtree->SetBranchStatus("muon.has_cluster_floor", 1);
		bxtree->SetBranchStatus("laben.clusters.nhits_short", 1);
	}
}

std::string make_string(const std::string& key, Float_t val){
	std::ostringstream oss;
	oss<<key<<val;
	return oss.str();
}

std::string rootfilename(const std::string& source_id, Int_t run, bool is_data){
	if (is_data){
		const std::string runtype = (source_id.substr(0,6));
		const std::string where   = (runtype != "normal") ? "ancillary/"  : "";
		const std::string flag    = (runtype != "normal") ? "_source"     : "";
		stringstream basepath_ss;
		basepath_ss.str("");
		basepath_ss<<"/bxstorage/rootfiles/cycle_"<<CYCLE_NUMBER<<"/2009/";
		const std::string basepath = basepath_ss.str();
		std::string path;
		if ((run>=10301) && (run<=10413)) path = "Jun_14/" + where;
		else if ((run>=10414) && (run<=10481)) path = "Jun_21/" + where;
		else if ((run>=10537) && (run<=10585)) path = "Jul_12/" + where;
		else if ((run>=10586) && (run<=10663)) path = "Jul_19/" + where;
		ostringstream oss;
		oss<<basepath<<path<<"Run0"<<run<<flag<<"_c"<<CYCLE_NUMBER<<".root";
		return (oss.str());
	}
	else{
		return("../mc_validation/" + source_id + ".root");
	}
}

void hstyle(TH1D& h, Int_t color, bool data) {
	if (!data){ 
		h.SetFillStyle(3001);
		h.SetFillColor(color);
	}
	h.SetLineColor(color);
	h.SetLineWidth(2);
}

void draw_test(TDirectory* dir, TH1D& h_mc, TH1D* h_data, const std::string& hname, const std::string& source_id, TLegend& leg){
	if (h_mc.GetEntries()>0){
		h_data = (TH1D*) dir->Get(hname.c_str());
		if (!h_data){
			std::cerr<<"Warning: no "<<hname<<" in "<<source_id<<"/data\n";
			exit(40);
		}
		h_data->SetName(("data"+hname).c_str());
		// ranges
		if ((hname=="nnpmts") || (hname=="nnhits") || (hname=="ncharge") || (hname=="ncharge_pmts") || (hname=="nnpmts_conc") || (hname=="nnpmts_noconc") || (hname=="nnpmts_near") || (hname=="nnpmts_far") || (hname=="raw_charge") || (hname=="raw_charge_1stchhit") || (hname=="raw_charge_no1stchhit") || (hname=="charge_1stchhit") || (hname=="charge_no1stchhit") || (hname=="ncharge_far") || (hname=="gatti") || (hname=="peak_time") || (hname=="mean_time_short") || (hname=="rms_time_short") || (hname=="x") || (hname=="y") || (hname=="z") || (hname=="hits_on_pmts") || (hname=="charge_on_pmts") || (hname=="charge_on_hits") || (hname=="conc_on_noconc")){
			Double_t mean = h_data->GetMean();
			Double_t rms = h_data->GetRMS();
			h_data->GetXaxis()->SetRangeUser(mean - 3.5*rms, mean + 4.0*rms);
		}
		else if (hname=="rec_times") h_data->GetXaxis()->SetRangeUser(0, 60);
		
		// drawings
		h_data->SetStats(0);
		h_mc.SetStats(0);
		h_data->DrawNormalized();
		h_mc.DrawNormalized("same");
		h_data->DrawNormalized("same"); // black line over blue fill
		
		// tests
		Double_t chi2 = 0.;
		Int_t ndf = 0;
		const Int_t nbins = h_data->GetNbinsX()-1;
		const Double_t Nd = h_data->Integral(1, nbins);
		const Double_t Nm = h_mc.Integral(1, nbins);

		for (Int_t i=1; i<nbins; i++){
			Double_t nd = h_data->GetBinContent(i);
			Double_t nm = h_mc.GetBinContent(i);
			if ((nm+nd)>10){
				chi2 += TMath::Power((Nd*nm - Nm*nd), 2)/(nm + nd);
				ndf++;
				}
		}
		chi2 /= (Nd*Nm);
		Double_t ks_prob = h_mc.KolmogorovTest(h_data);
		
		// legend
		leg.SetX1NDC(0.7);
		leg.SetX2NDC(1.);
		leg.SetY1NDC(0.85);
		leg.SetY2NDC(0.55);
		leg.SetHeader(source_id.c_str());
		leg.AddEntry(h_data, "Data");
		leg.AddEntry(&h_mc, "MC");
		ostringstream oss_ks, oss_chi2;
		oss_ks.precision(3);
		oss_ks<<"KS Prob. "<<ks_prob;
		oss_chi2.precision(3);
		oss_chi2<<"#Chi^{2}/NDF = "<<chi2<<"/"<<ndf;
		leg.AddEntry(&h_mc, oss_chi2.str().c_str());
		leg.AddEntry(&h_mc, oss_ks.str().c_str());
		leg.Draw();
	}
}

void set_bxstyle(){
	gStyle->SetCanvasColor(kWhite);
	gStyle->SetTitleFillColor(kWhite);
	gStyle->SetOptTitle(0);
	gStyle->SetFuncWidth(2);
	gStyle->SetHistLineWidth(1);
	gStyle->SetLegendBorderSize(0);
	//gStyle->SetOptFit(1111);
	gStyle->SetStatBorderSize(0);
	gStyle->SetTitleBorderSize(0);
	gStyle->SetDrawBorder(0);
	gStyle->SetLabelSize(.04,"xyz");
	gStyle->SetTitleSize(.04,"xyz");
	gStyle->SetLabelFont(102,"xyz");
	gStyle->SetOptStat("ne");
	gStyle->SetStatFont(102);
	gStyle->SetTitleFont(102,"xyz");
	gStyle->SetTitleFont(102,"pad");
	gStyle->SetStatStyle(0);
	gStyle->SetStatX(1);
	gStyle->SetStatY(1);
	gStyle->SetStatW(.2);
	gStyle->SetStatH(.15);
	gStyle->SetTitleStyle(0);
	gStyle->SetTitleX(.2);
	gStyle->SetTitleW(.65);
	gStyle->SetTitleY(.98);
	gStyle->SetTitleH(.07);
	gStyle->SetStatColor(0);
	gStyle->SetFillColor(10);
	gStyle->SetFillStyle(0);
	gStyle->SetTextFont(102);
	gStyle->SetCanvasBorderMode(0);
	gStyle->SetPadBorderMode(0);
	gStyle->SetFrameBorderMode(0);
	gStyle->SetDrawBorder(0);
	gStyle->SetPalette(1,0);
	gStyle->SetTitleYOffset(1.20);
	gStyle->SetPadRightMargin(0.05);
	gStyle->SetPadTopMargin(0.07);
}
