// PlotExample1.C
//
// This is a very simple example on how to plot a variable using root with Echidna events
// The example plots the trigger type
// The root file name is a parameter that must be given to this macro
// In case the parameter is not given, the default name "echidna.root" is used
//

#include <vector>

void PlotExample1(const char* name=0) {

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
	// the name of the histogram is "trgtype"
	// not to be confused with the histogram title, the second parameter!
	// numbers are the number of bins, low edge and high edge
	TH1F *h1 = new TH1F("trgtype","Trigger Type ",65,0.,65.);

	// loop on all events 
	for(int iev=0; iev<events; iev++) {

		// load iev-event in memory
		tree->GetEntry( iev );

		// fill histogram retrieving variable from the event using event "Getters"
		h1->Fill(ev->GetTrigger().GetTrgType());
	}
	
	// draw the histogram
	h1->Draw();

	// do NOT delete or close the file
	// the histogram will be destroyed!

    return;
}


