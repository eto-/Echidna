// PlotExample2.C
// $Id: PlotExample2.C,v 1.2 2005/02/14 16:42:44 zavatare Exp $
// This is a very simple example on how to plot a variable using root with Echidna events
// The example plots the difference between the montecarlo true energy and the detected
// energy (assuming so far a very simple correnspondence of 1 MeV --> 300 pmt hits)
//
// The root file name is a parameter that must be given to this macro
// In case the parameter is not given, the default name "echidna.root" is used
// This macro can be obviously used on geneb mc events only!
//
//

#include <vector>

void PlotExample2(const char* name=0) {

	// reset root to have a neat starting point
	gROOT->Reset();
	
	// choose file name
	char fname[500];
	if (name) 
		strcpy(fname,name);            // use given name
	else
		strcpy(fname,"echidna.root");  // use default name
    
	// open root file by creating object "f" of type "TFile" and passing
	// file name as parameter to the contructor; root opens file automatically; 
	TFile* f=new TFile(fname);

	// get the event tree from the file
	TTree* tree = f->Get("bxtree");

	// ask how many events are stored in this tree
	int events = tree->GetEntries();

	// create a branch 
	TBranch *b = tree->GetBranch("events");

	// create an event (empty!)
	BxEvent *ev = new BxEvent();

	// tell the branch where this event is
	b->SetAddress(&ev);

	// create an empty histogram
	// the name of the histogram is "eee"
	// not to be confused with the histogram title, the second parameters
	// numbers are the number of bins, low edge and high edge
	TH1F *h1 = new TH1F("eee","Geneb Energy - Reconstructed Energy",400,-2.,2.);

	// loop on all events skipping the first 1000 (pre-calibration events)
	// (start from 0 if you wish to include pre-calibrations events in your analysis)
	for(int iev=1001; iev<events; iev++) {

		// load iev-event in memory
		tree->GetEntry( iev );

		// get the montecarlo truth object of the event
		BxMcTruth& truth = ev->GetMcTruth();

		// get the whole frame vector
		std::vector<BxMcTruthFrame>& frames = truth.GetFrames();

		// get first energy deposit frame of geneb
		if ( frames.size() > 0) {
			BxMcTruthFrame& frame = frames[0];

			// compute "energy"
			if ( ev->GetLaben().GetClusters().size() > 0 ) {
				Float_t energy = (ev->GetLaben().GetClusters()[0].GetNpmt())/300.;
		
				// fill histogram with events of trigger type 1 (normal triggers)
				if (ev->GetTrigger().GetTrgType() == 1 )
					h1->Fill( frame.GetGenebEee() - energy );
			}
		}
	}
	// draw the histogram
	h1->Draw();

	// do NOT delete or close the file
	// the histogram will be destroyed!

    return;
}


