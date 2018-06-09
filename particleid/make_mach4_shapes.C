//make_echidna_shapes.C
//Macro to produce shape files from the Echidna root file. Need to have saved reconstructed clusters in Echidna.
//Alvaro Chavarria, Princeton
//3 August 2009

#include <iostream>
#include <sstream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <list>
#include <vector>
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TVector3.h"

using namespace std;

//Method to find peak

Float_t get_peak(vector<Float_t>::iterator begin, vector<Float_t>::iterator end)
{
	if (begin == end) return 0;
	
	unsigned n = 0;
	for (vector<Float_t>::iterator i = begin; i != end; i++) n++;
	
	unsigned peak_sample = 2 * (int)std::sqrt(1.0 * n) + 1;
	vector<Float_t> window(peak_sample, 0);
	
	unsigned j = 0;
	vector<Float_t>::iterator i = begin;
	for (; j < peak_sample && j < n; i++, j++)
		window[j] = *i;
	
	// If 5 or fewer hits, just return the median
	if (n <= 5) switch (n) {
		case 1: case 3: case 5: return window[n/2]; break;
		case 2: return (window[0] + window[1]) / 2; break;
		case 4: return (window[1] + window[2]) / 2; break;
	}
	
	Float_t t_peak = window[(peak_sample - 1)/2], t_diff = window[peak_sample - 1] - window[0];
	
	for (; i != end; i++, j++) {
		Float_t current = *i;
		window[j % peak_sample] = current;
		Float_t t_diff_test = current - window[(j + 1) % peak_sample];
		if (t_diff_test < t_diff) {
			t_peak = window[(j + 1 + peak_sample/2) % peak_sample];
			t_diff = t_diff_test;
		}
	}
	return t_peak;
}


//Now the main method

void make_shapes(TString echidna_filename, TString output_filename, Int_t num_events){
	

	//read root files and set output file and tree.
	
	TFile* echidna_f = new TFile( echidna_filename );
	TTree* echidna_t = (TTree*) echidna_f->Get("bxtree");
	
	TFile* output_f = new TFile( output_filename, "RECREATE" );
	
	//to use same settings as the histograms in Echidna/particleid/data/mach4_shapes.root
	
	Float_t t_min = -50;
	Float_t t_max = 950;
	
	TH1F* beta_hist = new TH1F( "beta", "beta", 500, t_min, t_max );
	TH1F* alpha_hist = new TH1F( "alpha", "alpha", 500, t_min, t_max);
	
	//Make the BxEvent
	BxEvent *ev = 0;
	echidna_t->SetBranchAddress("events", &ev);
	
	//BiPo criteria
	Float_t bi_nhits_max = 800; //these nhits values obtained for run 10321
	Float_t bi_nhits_min = 400;
	Float_t po_nhits_max = 290;
	Float_t po_nhits_min =  200;
	Float_t ds_max = 0.8; //in m
	Double_t dt_max = 500000; //in ns
	Float_t r_max = 4;
	
	//get number of entries
	Int_t entries = echidna_t->GetEntries();
	Int_t counter = 0;
	
	//because of the way Echidna tree is structure it is easier to look only for sequence of events in different triggers. Make sure they are both single cluster triggers. Effectively making a min dt cut.
	
	for( Int_t i = 0; i < entries - 1; i++ ){
		
		if ((i % 10000) == 0)
			cout<<"Processing Cluster No: "<< i << endl;
			
		echidna_t->GetEntry(i);
		//store first variables.
		
		//check to make sure that the cluster has rec hits stored
		if( !( ev->GetLaben().HasRecHits() ) ){
			cout << "Event has no reconstructed hits" << endl;
			continue;
		}
		
		Int_t nclusters1 = ev->GetLaben().GetClusters().size();
		if( nclusters1 != 1 ) continue;
		
		TVector3 pos1(ev->GetLaben().GetClusters()[0].GetPositionMach4().GetX(), ev->GetLaben().GetClusters()[0].GetPositionMach4().GetY(), ev->GetLaben().GetClusters()[0].GetPositionMach4().GetZ() );

		Double_t time1 = ev->GetLaben().GetTriggerTime();
		Float_t nhits1 = ev->GetLaben().GetClusters()[0].GetNpmts();
		Int_t btb1 = ev->GetTrigger().GetBtbInputs();
		Int_t ttype1 = ev->GetTrigger().GetTrgType();
		
		echidna_t->GetEntry(i+1);
		//store second variables.
		
		Int_t nclusters2 = ev->GetLaben().GetClusters().size();
		if( nclusters2 != 1 ) continue;		
		
		TVector3 pos2( ev->GetLaben().GetClusters()[0].GetPositionMach4().GetX(), ev->GetLaben().GetClusters()[0].GetPositionMach4().GetY(), ev->GetLaben().GetClusters()[0].GetPositionMach4().GetZ() );

		Double_t time2 = ev->GetLaben().GetTriggerTime();
		Float_t nhits2 =  ev->GetLaben().GetClusters()[0].GetNpmts();
		Int_t btb2 = ev->GetTrigger().GetBtbInputs();
		Int_t ttype2 = ev->GetTrigger().GetTrgType();
		
		Double_t dt = time2 - time1;
		Float_t ds = (pos1 - pos2).Mag();
		
		//cout << btb1 << " " << btb2 << " " << ttype1 << " " << ttype2 << " " << dt << " " << ds << " " << pos1.Mag() << endl;
		
		Bool_t bipo = false;
		
		if(btb2 == 0 && btb1 == 0 && ttype1 == 1 && ttype2 == 1 && dt < dt_max && ds < ds_max && pos1.Mag() < r_max && pos2.Mag() < r_max && nhits1 > bi_nhits_min && nhits1 < bi_nhits_max && nhits2 < po_nhits_max && nhits2 > po_nhits_min ) bipo = true;
		
		//check if bipo.
		
		if( bipo ){
			
			///COPYING THE CODE HERE TWICE IS NASTY BUT CINT WONT STOP FUCKING WITH ME
			
			//ONE
			vector<Float_t> times;
			Int_t size =  ev->GetLaben().GetRecHits().size();
			
			for( Int_t j = 0; j <  size; j++){
				
				Float_t hit_time = ev->GetLaben().GetRecHits()[j].GetTime();
				times.push_back(hit_time);
			}
			
			sort(times.begin(), times.end());
			
			Float_t t_peak = get_peak(times.begin(), times.end());
			
			for(Int_t j = 0; j < size; j++) alpha_hist->Fill(times[j] - t_peak);
			
			///////////////////////////////////////////////////////////////////
			
			echidna_t->GetEntry(i);
			
			//TWO
			
			times.clear();
			size = ev->GetLaben().GetRecHits().size();
			
			for( Int_t j = 0; j <  size; j++){
				
				Float_t hit_time = ev->GetLaben().GetRecHits()[j].GetTime();
				times.push_back(hit_time);
			}
			
			sort(times.begin(), times.end());
			
			Float_t t_peak = get_peak(times.begin(), times.end());
			
			for(Int_t j = 0; j < size; j++) beta_hist->Fill(times[j] - t_peak);
			
			//////////////////////////////////////////////////////////////////
			
			counter++;
			i++;
		}
		
		if( counter >= num_events ) break;
		
	}
	
	//Write tree and close output file.
	beta_hist->Write();
	alpha_hist->Write();
	output_f->Close();

}
