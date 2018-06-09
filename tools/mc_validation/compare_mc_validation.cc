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

std::string make_string(const std::string&, Float_t);
std::string get_name(const std::string&); // remove path and extension
void draw_compare(TDirectory* data, std::vector<TDirectory*>& mc, const std::vector<std::string>&, const std::string&, const std::string&, TLegend&);
void hstyle(TH1D&, Int_t, bool);
void set_bxstyle();

int main(int argc, char** argv){
	set_bxstyle();
	std::vector<std::string> mc_validation_filename_vector; // full path + extension

	for (int i=1; i<argc; i++){
		const std::string stmp = argv[i];
		const size_t sep = stmp.find(':');
		const std::string key   = stmp.substr(0,sep);
		const std::string value = stmp.substr(sep+1,stmp.length());
		if (key=="add"){
			mc_validation_filename_vector.push_back(value);
			}
		else {
			std::cerr<<"WARNING unknown option: "<<stmp<<'\n';
			exit(255);
		}
	}

	// check: at least 1 file
	if (mc_validation_filename_vector.size()<1){
		cerr<<"WARNING: not enough files to compare. Use add:\n";
		exit(1);
	}

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
	const std::string outfile_name = "mc_validation_comparison.root";
	TFile outfile(outfile_name.c_str(), "RECREATE");
	if (outfile.IsZombie()){
		std::cerr<<"ERROR: I can't open "<<outfile_name<<std::endl;
		exit(-1);
	}

	TFile* rootfile_data;
	const std::string rootfile_data_name = "tools/mc_validation/calibration_sources.root";
	// open data calibration source rootfile
	rootfile_data = new TFile(rootfile_data_name.c_str(), "READ");
	if (rootfile_data->IsZombie()){
		std::cerr<<"ERROR: I can't open "<<rootfile_data_name<<std::endl;
		exit(-1);
	}

	TFile** rootfile_mc = new TFile* [mc_validation_filename_vector.size()];
	for (unsigned i=0; i<mc_validation_filename_vector.size(); i++){
		rootfile_mc[i] = new TFile(mc_validation_filename_vector.at(i).c_str(), "READ");
		if (rootfile_mc[i]->IsZombie()){
			std::cerr<<"ERROR: I can't open "<<mc_validation_filename_vector[i]<<std::endl;
			exit(-1);
		}
		cout<<mc_validation_filename_vector[i]<<" loaded\n";
	}

	// main loop on simulations
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

		// copy all histograms from calibration_source.root (data rootfile)
		// in the folder source_id[i]/data in mc_validation.root (mc rootfile)
		TDirectory* data_dir_data = rootfile_data->GetDirectory(source_id[i].c_str());
		TDirectory* this_dir_data = this_dir->mkdir("data");
		this_dir_data->cd();
		TKey* key;
		TIter nextkey(data_dir_data->GetListOfKeys());
		while ((key = (TKey*) nextkey())){
			TObject* obj = key->ReadObj();
			if ((obj->IsA())->InheritsFrom("TH1")){
				((TH1D*) obj)->Write();
			}
		}

		// do the same here for montecarlo validation histograms
		// TDirectory* are saved for future use
		std::vector<TDirectory*> mc_directory_vector;
		std::vector<std::string> mc_validation_vector; 		// shortened name
		for (unsigned ul=0; ul < mc_validation_filename_vector.size(); ul++){
			TDirectory* mc_dir_mc = rootfile_mc[ul]->GetDirectory(source_id[i].c_str());
			// remove path and .root from mc_validation_filename_vector
			const std::string mc_validation = get_name(mc_validation_filename_vector.at(ul));
			TDirectory* this_dir_mc = this_dir->mkdir(mc_validation.c_str());
			this_dir_mc->cd();
			TKey* key;
			TIter nextkey(mc_dir_mc->GetListOfKeys());
			while ((key = (TKey*) nextkey())){
				TObject* obj = key->ReadObj();
				if ((obj->IsA())->InheritsFrom("TH1")){
					((TH1D*) obj)->Write();
				}
			}
			mc_directory_vector.push_back(this_dir_mc);
			mc_validation_vector.push_back(mc_validation);
		}

		this_dir->cd();

		// canvas composition here:
		// c1 -> general (nnpmts, nnhits, hit_times, rec_times)
		// c2 -> near_far (nnpmts_near, nnpmts_far, hit_times_near, hit_times_far)
		// c3 -> conc_noconc (nnpmts_conc, nnpmts_noconc, hit_times_conc, hit_times_noconc)
		// c4 -> first_hit (hit_times_1stch, hit_times_no1stch, hit_times_near_1stch, hit_times_far_1stch)
		// c5 -> charge (ncharge, charge_pmts, raw_charge, raw_charge_1stchhit);
		// c6 -> ratios (hits_on_pmts, charge_on_hits, charge_on_pmts, conc_on_noconc);
		// c7 -> times (gatti, peak_time, mean_time_short, rms_time_short)
		// c8 -> position (r0, z, x, y);
		// c9 -> decoded (decoded_times, decoded_dt_12, decoded_times_order1, decoded_times_order2);
		//
		TDirectory* this_dir_canvas = this_dir->mkdir("canvas");
	
		// c1
		this_dir_canvas->cd();
		ostringstream oss_title;
		oss_title<<"General comparison histograms for "<<source_id[i];
		TCanvas c1("c1_general", oss_title.str().c_str());
		c1.Divide(2,2);
		c1.cd(1);
		TLegend leg_1_1(0.65, 0.90, 0.98, 0.55);
		draw_compare(this_dir_data, mc_directory_vector, mc_validation_vector, "nnpmts", source_id[i], leg_1_1);
		c1.cd(2);
		TLegend leg_1_2(0.65, 0.90, 0.98, 0.55);
		draw_compare(this_dir_data, mc_directory_vector, mc_validation_vector, "nnhits", source_id[i], leg_1_2);
		TVirtualPad* c1_3 = c1.cd(3);
		c1_3->SetLogy();
		TLegend leg_1_3(0.65, 0.90, 0.98, 0.55);
		draw_compare(this_dir_data, mc_directory_vector, mc_validation_vector, "hit_times", source_id[i], leg_1_3);
		c1.cd(4);
		TLegend leg_1_4(0.65, 0.90, 0.98, 0.55);
		draw_compare(this_dir_data, mc_directory_vector, mc_validation_vector, "rec_times", source_id[i], leg_1_4);
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
		draw_compare(this_dir_data, mc_directory_vector, mc_validation_vector, "nnpmts_near", source_id[i], leg_2_1);
		c2.cd(2);
		TLegend leg_2_2(0.65, 0.90, 0.98, 0.55);
		draw_compare(this_dir_data, mc_directory_vector, mc_validation_vector, "nnpmts_far", source_id[i], leg_2_2);
		TVirtualPad* c2_3 = c2.cd(3);
		c2_3->SetLogy();
		TLegend leg_2_3(0.65, 0.90, 0.98, 0.55);
		draw_compare(this_dir_data, mc_directory_vector, mc_validation_vector, "hit_times_near", source_id[i], leg_2_3);
		TVirtualPad* c2_4 = c2.cd(4);
		c2_4->SetLogy();
		TLegend leg_2_4(0.65, 0.90, 0.98, 0.55);
		draw_compare(this_dir_data, mc_directory_vector, mc_validation_vector, "hit_times_far", source_id[i], leg_2_4);
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
		draw_compare(this_dir_data, mc_directory_vector, mc_validation_vector, "nnpmts_conc", source_id[i], leg_3_1);
		c3.cd(2);
		TLegend leg_3_2(0.65, 0.90, 0.98, 0.55);
		draw_compare(this_dir_data, mc_directory_vector, mc_validation_vector, "nnpmts_noconc", source_id[i], leg_3_2);
		TVirtualPad* c3_3 = c3.cd(3);
		c3_3->SetLogy();
		TLegend leg_3_3(0.65, 0.90, 0.98, 0.55);
		draw_compare(this_dir_data, mc_directory_vector, mc_validation_vector, "hit_times_conc", source_id[i], leg_3_3);
		TVirtualPad* c3_4 = c3.cd(4);
		c3_4->SetLogy();
		TLegend leg_3_4(0.65, 0.90, 0.98, 0.55);
		draw_compare(this_dir_data, mc_directory_vector, mc_validation_vector, "hit_times_noconc", source_id[i], leg_3_4);
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
		draw_compare(this_dir_data, mc_directory_vector, mc_validation_vector, "hit_times_1stchhit", source_id[i], leg_4_1);
		TVirtualPad* c4_2 = c4.cd(2);
		c4_2->SetLogy();
		TLegend leg_4_2(0.65, 0.90, 0.98, 0.55);
		draw_compare(this_dir_data, mc_directory_vector, mc_validation_vector, "hit_times_no1stchhit", source_id[i], leg_4_2);
		TVirtualPad* c4_3 = c4.cd(3);
		c4_3->SetLogy();
		TLegend leg_4_3(0.65, 0.90, 0.98, 0.55);
		draw_compare(this_dir_data, mc_directory_vector, mc_validation_vector, "hit_times_near_1stchhit", source_id[i], leg_4_3);
		TVirtualPad* c4_4 = c4.cd(4);
		c4_4->SetLogy();
		TLegend leg_4_4(0.65, 0.90, 0.98, 0.55);
		draw_compare(this_dir_data, mc_directory_vector, mc_validation_vector, "hit_times_far_1stchhit", source_id[i], leg_4_4);
		c4.cd(0);
		this_dir_canvas->cd();
		c4.Write();

		// c5
		this_dir_canvas->cd();
		ostringstream oss_title5;
		oss_title5<<"Charge comparison histograms for "<<source_id[i];
		TCanvas c5("c5_charge", oss_title5.str().c_str());
		c5.Divide(2,3);
		c5.cd(1);
		TLegend leg_5_1(0.65, 0.90, 0.98, 0.55);
		draw_compare(this_dir_data, mc_directory_vector, mc_validation_vector, "ncharge", source_id[i], leg_5_1);
		c5.cd(2);
		TLegend leg_5_2(0.65, 0.90, 0.98, 0.55);
		draw_compare(this_dir_data, mc_directory_vector, mc_validation_vector, "ncharge_pmts", source_id[i], leg_5_2);
		c5.cd(3);
		TLegend leg_5_3(0.65, 0.90, 0.98, 0.55);
		draw_compare(this_dir_data, mc_directory_vector, mc_validation_vector, "raw_charge", source_id[i], leg_5_3);
		c5.cd(4);
		TLegend leg_5_4(0.65, 0.90, 0.98, 0.55);
		draw_compare(this_dir_data, mc_directory_vector, mc_validation_vector, "charge_1stchhit", source_id[i], leg_5_4);
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
		TLegend leg_6_1(0.65, 0.90, 0.98, 0.55);
		draw_compare(this_dir_data, mc_directory_vector, mc_validation_vector, "hits_on_pmts", source_id[i], leg_6_1);
		c6.cd(2);
		TLegend leg_6_2(0.65, 0.90, 0.98, 0.55);
		draw_compare(this_dir_data, mc_directory_vector, mc_validation_vector, "charge_on_pmts", source_id[i], leg_6_2);
		c6.cd(3);
		TLegend leg_6_3(0.65, 0.90, 0.98, 0.55);
		draw_compare(this_dir_data, mc_directory_vector, mc_validation_vector, "charge_on_hits", source_id[i], leg_6_3);
		c6.cd(4);
		TLegend leg_6_4(0.65, 0.90, 0.98, 0.55);
		draw_compare(this_dir_data, mc_directory_vector, mc_validation_vector, "conc_on_noconc", source_id[i], leg_6_4);
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
		TLegend leg_7_1(0.65, 0.90, 0.98, 0.55);
		draw_compare(this_dir_data, mc_directory_vector, mc_validation_vector, "gatti", source_id[i], leg_7_1);
		c7.cd(2);
		TLegend leg_7_2(0.65, 0.90, 0.98, 0.55);
		draw_compare(this_dir_data, mc_directory_vector, mc_validation_vector, "peak_time", source_id[i], leg_7_2);
		c7.cd(3);
		TLegend leg_7_3(0.65, 0.90, 0.98, 0.55);
		draw_compare(this_dir_data, mc_directory_vector, mc_validation_vector, "mean_time_short", source_id[i], leg_7_3);
		c7.cd(4);
		TLegend leg_7_4(0.65, 0.90, 0.98, 0.55);
		draw_compare(this_dir_data, mc_directory_vector, mc_validation_vector, "rms_time_short", source_id[i], leg_7_4);
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
		TLegend leg_8_1(0.65, 0.90, 0.98, 0.55);
		draw_compare(this_dir_data, mc_directory_vector, mc_validation_vector, "r0", source_id[i], leg_8_1);
		c8.cd(2);
		TLegend leg_8_2(0.65, 0.90, 0.98, 0.55);
		draw_compare(this_dir_data, mc_directory_vector, mc_validation_vector, "z", source_id[i], leg_8_2);
		c8.cd(3);
		TLegend leg_8_3(0.65, 0.90, 0.98, 0.55);
		draw_compare(this_dir_data, mc_directory_vector, mc_validation_vector, "x", source_id[i], leg_8_3);
		c8.cd(4);
		TLegend leg_8_4(0.65, 0.90, 0.98, 0.55);
		draw_compare(this_dir_data, mc_directory_vector, mc_validation_vector, "y", source_id[i], leg_8_4);
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
		draw_compare(this_dir_data, mc_directory_vector, mc_validation_vector, "decoded_times", source_id[i], leg_9_1);
		TVirtualPad* c9_2 = c9.cd(2);
		c9_2->SetLogy();
		TLegend leg_9_2(0.65, 0.90, 0.98, 0.55);
		draw_compare(this_dir_data, mc_directory_vector, mc_validation_vector, "decoded_dt_12", source_id[i], leg_9_2);
		TVirtualPad* c9_3 = c9.cd(3);
		c9_3->SetLogy();
		TLegend leg_9_3(0.65, 0.90, 0.98, 0.55);
		draw_compare(this_dir_data, mc_directory_vector, mc_validation_vector, "decoded_times_order1", source_id[i], leg_9_3);
		TVirtualPad* c9_4 = c9.cd(4);
		c9_4->SetLogy();
		TLegend leg_9_4(0.65, 0.90, 0.98, 0.55);
		draw_compare(this_dir_data, mc_directory_vector, mc_validation_vector, "decoded_times_order2", source_id[i], leg_9_4);
		c9.cd(0);
		this_dir_canvas->cd();
		c9.Write();


	} // end loop on simulation root files

	std::cout<<"Closing output rootfile "<<outfile_name<<'\n';
	outfile.Close();
		rootfile_data->Close();
	for (unsigned i=0; i<mc_validation_filename_vector.size(); i++){
		rootfile_mc[i]->Close();
	}
	std::cout<<"NORMAL END\n";
	return 0;
}


std::string make_string(const std::string& key, Float_t val){
	std::ostringstream oss;
	oss<<key<<val;
	return oss.str();
}

std::string get_name(const std::string& fullname){

	// remove path
	size_t begin = fullname.rfind("/");
	const string extname = (begin!=string::npos) ? fullname.substr(begin + 1): fullname;

	// remove extension
	size_t end   = extname.find(".root");
	const string name = extname.substr(0, end);

	return name;
}

void hstyle(TH1D& h, Int_t color, bool data) {
	if (!data){ 
		h.SetFillStyle(3001);
		h.SetFillColor(color);
	}
	h.SetLineColor(color);
	h.SetLineWidth(2);
}

void draw_compare(TDirectory* data_dir, std::vector<TDirectory*>& mc_dir, const std::vector<std::string>& mc_validation, const std::string& hname, const std::string& source_id, TLegend& leg){
	// get mc_histograms
	TH1D** h_mc = new TH1D*[mc_dir.size()];
	for (unsigned i=0; i<mc_dir.size(); i++){
		h_mc[i] = (TH1D*) mc_dir.at(i)->Get(hname.c_str());
		if (!h_mc[i]){
			std::cerr<<"Warning: no "<<hname<<" in "<<source_id<<" of "<<mc_validation[i]<<endl;
			exit(40);
		}
		h_mc[i]->SetName((mc_validation[i] + hname).c_str());
	}

	TH1D* h_data = (TH1D*) data_dir->Get(hname.c_str());
	if (!h_data){
		std::cerr<<"Warning: no "<<hname<<" in "<<source_id<<"/data\n";
		exit(40);
	}
	h_data->SetName(("data"+hname).c_str());

	if (h_data->GetEntries()>0){
		// ranges
		if ((hname=="nnpmts") || (hname=="nnhits") || (hname=="ncharge") || (hname=="ncharge_pmts") || (hname=="nnpmts_conc") || (hname=="nnpmts_noconc") || (hname=="nnpmts_near") || (hname=="nnpmts_far") || (hname=="raw_charge") || (hname=="raw_charge_1stchhit") || (hname=="raw_charge_no1stchhit") || (hname=="charge_1stchhit") || (hname=="charge_no1stchhit") || (hname=="ncharge_far") || (hname=="gatti") || (hname=="mean_time_short") || (hname=="rms_time_short") || (hname=="x") || (hname=="y") || (hname=="z") || (hname=="hits_on_pmts") || (hname=="charge_on_pmts") || (hname=="charge_on_hits") || (hname=="conc_on_noconc")){
			Double_t mean = h_data->GetMean();
			Double_t rms = h_data->GetRMS();
			h_data->GetXaxis()->SetRangeUser(mean - 3.5*rms, mean + 4.0*rms);
		}
		else if (hname=="rec_times") h_data->GetXaxis()->SetRangeUser(0, 60);
		
		// drawings
		h_data->SetStats(0);
		h_data->DrawNormalized();

		// legend
		leg.SetX1NDC(0.7);
		leg.SetX2NDC(1.);
		leg.SetY1NDC(0.85);
		leg.SetY2NDC(0.55);
		leg.SetHeader(source_id.c_str());
		leg.AddEntry(h_data, "Data");

		for (unsigned u=0; u<mc_validation.size(); u++){
			h_mc[u]->SetStats(0);
			//h_mc[u]->SetFillStyle(3305 + 10*u);
			h_mc[u]->SetFillStyle(0);
			h_mc[u]->SetLineColor(u+2);
			h_mc[u]->DrawNormalized("same");
			// tests
			Double_t chi2 = 0.;
			Int_t ndf = 0;
			const Int_t nbins = h_data->GetNbinsX()-1;
			const Double_t Nd = h_data->Integral(1, nbins);
			const Double_t Nm = h_mc[u]->Integral(1, nbins);
			for (Int_t i=1; i<nbins; i++){
				Double_t nd = h_data->GetBinContent(i);
				Double_t nm = h_mc[u]->GetBinContent(i);
				if ((nm+nd)>10){
					chi2 += TMath::Power((Nd*nm - Nm*nd), 2)/(nm + nd);
					ndf++;
				}
			}
			chi2 /= (Nd*Nm);
			leg.AddEntry(h_mc[u], mc_validation.at(u).c_str());
			ostringstream oss_chi2;
			oss_chi2.precision(3);
			oss_chi2<<"#Chi^{2}/NDF = "<<chi2<<"/"<<ndf;
			leg.AddEntry(h_mc[u], oss_chi2.str().c_str());
		}
		h_data->DrawNormalized("same"); // black line over blue fill
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
