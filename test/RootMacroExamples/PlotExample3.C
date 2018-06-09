// PlotExample3.C
// $Id: PlotExample3.C,v 1.1 2005/02/14 16:42:44 zavatare Exp $
// same as example 2 but with a gaussian fit added

#include <vector>

void PlotExample3(const char* name=0) {

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
	
	// create a function for fitting
	TF1 *g1 = new TF1("fit","gaus",-2.,2.);
 	
	// loop on all events skipping the first 1000 (pre-calibration events)
	// (start from 0 if you wish to include pre-calibrations events in your analysis)
	for(int iev=1001; iev<events; iev++) {

		// load iev-event in memory
		tree->GetEntry( iev );

		// get the montecarlo truth object of the event
		BxMcTruth& truth = ev->GetMcTruth();

		// get the whole frame vector
		std::vector<BxMcTruthFrame>& frames = truth.GetFrames();


		if ( frames.size() > 0 ) {
			// get first energy deposit frame of geneb
			BxMcTruthFrame& frame = frames[0];

			// compute "energy"
			BxLaben& laben=ev->GetLaben();
			std::vector<BxLabenCluster>& clusters = laben.GetClusters();
			Int_t Nsize = clusters.size();
			if (Nsize > 0)   {
		
		  		Float_t energy = (ev->GetLaben().GetClusters()[0].GetNpmt())/300.;
		
		  		// fill histogram with events of trigger type 1 (normal triggers)
		  
				if (ev->GetTrigger().GetTrgType() == 1 )
					h1->Fill( frame.GetGenebEee() - energy );
			}
		}
	}
	
	// draw the histogram

	h1->Draw();

	// set starting point parameters to mean value, rms and max height for histogram
    g1->SetParameter(0,h1->GetMaximum());
	g1->SetParameter(1,h1->GetMean());
	g1->SetParameter(2,h1->GetRMS());

	cout << "Fit starting point: (" << h1->GetMaximum() << "," << h1->GetMean() << "," << h1->GetRMS() << ")"<< endl; 

	h1->Fit(g1,"R");
	// do NOT delete or close the file
	// the histogram will be destroyed!

    return;
}
