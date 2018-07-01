// 
// class for alpha-beta discrimination
// MP 22-02-07
// Implemented methods:
//

#include <fstream>
#include "TMath.h"
#include "TH1F.h"
#include "bx_alphabeta.hh"
#include "bx_echidna_event.hh"
#include "BxEvent.hh"

#ifdef REPLACE_GETLINE
void getline(char** s, size_t* n, FILE *fp) { fgets(*s,200,fp);}
#endif

//ClassImp(bx_alphabeta)

// ctor
bx_alphabeta::bx_alphabeta(const std::string& name) : message_client("bx_alphabeta"), 
                               hHisto(0), hHistoCumu(0), 
			       fReady(false), fPositionCorrection(true), 
			       fErrorCounter(0), fSmooth(0), fName(name) {
	fAlphaSamples.clear();
	fBetaSamples.clear();

	// helper histograms
	
	std::string s = "h_"+GetName()+"MeanAlpha";
	hMeanAlpha = new TH1F(s.c_str(),"Truncated mean for Alpha candidates",
			100,0.,50.);

	s = "h_"+GetName()+"MeanBeta";
	hMeanBeta = new TH1F(s.c_str(),"Truncated mean for Beta candidates",
			100,0.,50.);

	s = "h_"+GetName()+"ModaAlpha";
	hModaAlpha = new TH1F(s.c_str(),"Moda for Alpha candidates",
			100,0.,50.);

	s = "h_"+GetName()+"ModaBeta";
	hModaBeta = new TH1F(s.c_str(),"Moda for Beta candidates",
			100,0.,50.);

	// define default binning
	/*DefaultBinning[0] = -42.;
	int i,j;

	for(i=1; i<30; i++)
		DefaultBinning[i] = DefaultBinning[i-1] + 3.;

	for(j=i; j<60; j++)
		DefaultBinning[j] = DefaultBinning[j-1] + 3. + 20.*(TMath::Exp((j - 30.)/35.)-0.999);
	*/
		
	// define default binning
	DefaultBinning[0] = 0.;
	int i,j;

	for(i=1; i<30; i++)
		DefaultBinning[i] = DefaultBinning[i-1] + 3.;

	for(j=i; j<60; j++)
		DefaultBinning[j] = DefaultBinning[j-1] + 3. + 20.*(TMath::Exp((j - 30.)/35.)-0.999);
	

	NDefaultBins = j-1;

	current_sample.clear();

	svm_energy_file = 0;
	svm_sample_file = 0;
	svm_alpha_file = 0;
	svm_beta_file = 0;

}


// dtor
bx_alphabeta::~bx_alphabeta() {
}


// allow definition of time histogram binning
void bx_alphabeta::set_time_histogram(int nbins, float* binning) {
	if ( binning ) {

		std::string s = "h_" + GetName() + "Histo";
		hHisto = new TH1F(s.c_str(),"",nbins,binning);

		s = "h_" + GetName() + "HistoCumu";
		hHistoCumu = new TH1F(s.c_str(),"",nbins,binning);

		s = "h_" + GetName() + "AlphaTimeDist";
		hAlphaTimeDist = new TH1F(s.c_str(),"",nbins,binning);

		s = "h_" + GetName() + "BetaTimeDist";
		hBetaTimeDist = new TH1F(s.c_str(),"",nbins,binning);

		s = "h_"+GetName()+"AlphaTimeDistCumu";
		hAlphaTimeCumulativeDist = new TH1F(s.c_str(),"",nbins,binning);

		s = "h_" + GetName() + "BetaTimeDistCumu";
		hBetaTimeCumulativeDist = new TH1F(s.c_str(),"",nbins,binning);

		NDefaultBins = nbins;
		for(int i=0; i<=nbins; i++)
			DefaultBinning[i] = binning[i];
	}
	else {
		std::string s = "h_" + GetName() + "Histo";
		hHisto = new TH1F(s.c_str(),"",NDefaultBins,DefaultBinning);

		s = "h_" + GetName() + "HistoCumu";
		hHistoCumu = new TH1F(s.c_str(),"",NDefaultBins,DefaultBinning);

		s = "h_" + GetName() + "AlphaTimeDist";
		hAlphaTimeDist = new TH1F(s.c_str(),"",NDefaultBins,DefaultBinning);

		s = "h_" + GetName() + "BetaTimeDist";
		hBetaTimeDist = new TH1F(s.c_str(),"",NDefaultBins,DefaultBinning);

		s = "h_"+GetName()+"AlphaTimeDistCumu";
		hAlphaTimeCumulativeDist = new TH1F(s.c_str(),"",NDefaultBins,DefaultBinning);

		s = "h_" + GetName() + "BetaTimeDistCumu";
		hBetaTimeCumulativeDist = new TH1F(s.c_str(),"",NDefaultBins,DefaultBinning);

		get_message(bx_message::info) << GetName() << " Setting histograms with default binning" << dispatch;
	}
	//for(int i=0; i < NDefaultBins; i++) {
	//	std::cout << DefaultBinning[i] << " ";
	//}
	get_message(bx_message::info) << GetName() << " Binning: Start="<<DefaultBinning[0]<<"   End="<<DefaultBinning[NDefaultBins] << dispatch;
	get_message(bx_message::info) << GetName() << " Num Bins = " << NDefaultBins << dispatch;
}

// get time hits from bx_echidna_event and build bx_sample
bx_sample* bx_alphabeta::get_sample(const bx_echidna_event& ev, int ncluster, bx_sample* sampl ) {
	
	// check cluster number range (0..size-1)
	if (ncluster < 0 || ncluster >= ev.get_laben().get_nclusters()) {
		get_message(bx_message::error) << GetName() << " Wrong cluster number in bx_alphabeta::get_sample(bx_echidna_event). Ignoring." << dispatch;
		return 0;
	}

	// if provided, use existing sample, otherwise create a new one
	bx_sample *sample = 0;
	if (sampl)
		sample = sampl;
	else
		sample = new bx_sample();
	sample->clear();

        // if position correction is not required, used cluster hits
	if ( !fPositionCorrection ) {
		const bx_laben_cluster &cluster = ev.get_laben().get_cluster(ncluster);
		int size = cluster.get_clustered_nhits();
		for (int  i=0; i < size; i++) {
			double t = cluster.get_clustered_hit(i).get_time();
			sample->push_back( t );
		}
	}
	// otherwise use reco hits with correction of time of flight
	else {
		const bx_laben_ab_cluster &cluster = ev.get_laben().get_ab_cluster(ncluster);
		int size = cluster.get_rec_nhits();
		for (int  i=0; i < size; i++) {
			double t = cluster.get_rec_hit(i).get_time();
			sample->push_back( t );
		}
	}

        // ignore tiny events
	if ( sample->size() < 20 ) {
		get_message(bx_message::warn) << GetName() << " Cluster has less than 20 photons. Ignored." << dispatch;
		if (!sampl) delete sample;
		return 0;
	}

	// check that hit times start from 0; if not, correct it.
	if ( (*sample)[0] < -0.1 || (*sample)[0] > 0.1 ) {
		if (fErrorCounter == 0)
			get_message(bx_message::error) << GetName() << " First hit of event does not have time 0!" << dispatch;
		fErrorCounter++;
		int size = sample->size();
		for(int i=0; i<size; i++) {
			(*sample)[i] -= (*sample)[0];
		}
	} 

	return sample;
}

// get time hits from BxEvent and build bx_sample
bx_sample* bx_alphabeta::get_sample(const BxEvent& ev, int ncluster, bx_sample* sampl) {

	if (ncluster <= 0 || ncluster > ev.GetLaben().GetNClusters()) {
		get_message(bx_message::error) << GetName() << " Wrong cluster number in bx_alphabeta::get_sample(BxEvent). Ignoring." << dispatch;
		return 0;
	}

	int size = ev.GetLaben ().GetNClusteredHits();
    
	// if provided, use existing sample, otherwise create a new one
	bx_sample *sample = 0;
	if (sampl)
		sample = sampl;
	else
		sample = new bx_sample();
	sample->clear();

	// get hits time
	// can use uncorrected or corrected time for comparison
	for (int  i=0; i < size; i++) {
	        if (!fPositionCorrection) {
			const BxLabenClusteredHit &hit = ev.GetLaben().GetClusteredHits()[i];
			if ( hit.GetNumCluster() == ncluster )
				sample->push_back( hit.GetTime() );
	        }
		else {
			const BxLabenRecHit &hit = ev.GetLaben().GetRecHits()[i];
			if ( hit.GetNumCluster() == ncluster )
				sample->push_back( hit.GetTime() );
		}
	}

        // ignore tiny events
	if ( sample->size() < 20 ) {
		get_message(bx_message::warn) << GetName() << " Cluster has less than 20 photons. Ignored. " << dispatch;
		if (!sampl) delete sample;
		return 0;
	}
	
	// check that hit times start from 0; if not, correct it.
	if ( (*sample)[0] < -0.1 || (*sample)[0] > 0.1 ) {
		get_message(bx_message::error) << GetName() << " First hit of event does not have time 0! Corrected!" << dispatch;
		int size = sample->size();
		for(int i=0; i<size; i++) {
			(*sample)[i] -= (*sample)[0];
		}
	}

	return sample;
	
}



// add methods: they must be called to prepare learning
// add an echidna_event
void bx_alphabeta::add(const std::string& cname,const bx_echidna_event& ev, int ncluster) {
	bx_sample *s = get_sample(ev,ncluster);
	if (s)
		add(cname,s);
}

// add a root event
void bx_alphabeta::add(const std::string& cname,const BxEvent& ev, int ncluster) {
	bx_sample *s = get_sample(ev,ncluster);
	if (s)
		add(cname,s);
}

// add a sample
void bx_alphabeta::add(const std::string& cname, bx_sample* sample) {
	
	if ( fReady) {
		get_message(bx_message::error) << GetName() << " You cannot add samples to an already READY classifier" << dispatch;
		return;
	}

	bool alpha=false;
	if (cname=="Alpha") {
		alpha = true;
	} else if ( cname!="Beta" ) {
		get_message(bx_message::error) << GetName() << " Invalid class name. Only Alpha and Beta accepted" << dispatch;
		return;
	}

	if (!sample)
		return;

	float moda = 0.;
	float mean = 0.; // sample->truncated_mean(moda);

	// fill helper histos
	if ( alpha ) {
		hMeanAlpha->Fill(mean);
		hModaAlpha->Fill(moda);
	} else {
		hMeanBeta->Fill(mean);
		hModaBeta->Fill(moda);
	}

	// refer all times to the truncated mean
	int size = sample->size();
	for(int i=0; i<size; i++)
		(*sample)[i] -= mean;
	
	// add to samples
	if (alpha) {
		fAlphaSamples.push_back(sample);
	}
	else {
		fBetaSamples.push_back(sample);
	}
	
	// std::cout << "Mean = " << mean << "   First sample time: " << (*sample)[0] << std::endl;
}

//
// dump gatti weights on screen or on stream
//
void bx_alphabeta::dump_gatti(std::ostream& os) {
	os << std::endl;
	os << "dump of gatti cumulative weights\n";
	for(int i=1; i<=NDefaultBins; i++) {
		os << GattiCumuWeight[i] << " ";
	}
	os << std::endl;
	os << "dump of gatti non cumulative weights\n";
	for(int i=0; i<NDefaultBins; i++) {
		os << GattiWeight[i] << " ";
	}
	os << std::endl;
}

//
// learning for all classifiers
//
bool bx_alphabeta::learn() {
	
	if ( fReady )
		return true;
	
	if ( !hHisto )
		set_time_histogram();

	int n_samples_alpha = fAlphaSamples.size();
	int n_samples_beta = fBetaSamples.size();

	get_message(bx_message::info) << GetName() << " Start learning with " << n_samples_alpha << " alpha samples and " <<
		n_samples_beta << " beta samples" << dispatch;
	
	// open svm file
	FILE *svmfile = fopen("svmfile_samples.ab","wt+");
	
	// build alpha reference histogram
	for(int i=0; i<n_samples_alpha; i++) {
		hHisto->Reset("ICE");
		bx_sample &s = *(fAlphaSamples[i]);
		int n = s.size();
		for(int j=0; j<n; j++) {
			hAlphaTimeDist->Fill(s[j]);
			hHisto->Fill(s[j]);            // svm
		}
		float area = hHisto->Integral();
		if (area > 0. ) {
			fprintf(svmfile,"-1 ");                     // alpha type -1
			for(int j=1; j<=NDefaultBins; j++) {
				float val = hHisto->GetBinContent(j)/area;
				fprintf(svmfile,"%d:%06.4f ",j,val);
			}
			fprintf(svmfile,"\n");
		}
	}

	// build beta reference histogram
	for(int i=0; i<n_samples_beta; i++) {
		bx_sample &s = *(fBetaSamples[i]);
		int n = s.size();
		for(int j=0; j<n; j++) {
			hBetaTimeDist->Fill(s[j]);
			hHisto->Fill(s[j]);            // svm
		}
		// write svm file
		float area = hHisto->Integral();
		if (area > 0. ) {
			fprintf(svmfile,"+1 ");                     // beta type +1
			for(int j=1; j<=NDefaultBins; j++) {
				float val = hHisto->GetBinContent(j)/area;
				fprintf(svmfile,"%d:%06.4f ",j,val);
			}
			fprintf(svmfile,"\n");
		}
	}
	
	fclose(svmfile); // close svm file

	float a = hAlphaTimeDist->Integral();
	float b = hBetaTimeDist->Integral();

	for(int i=1; i<=NDefaultBins; i++) {
		hAlphaTimeDist->SetBinContent(i,hAlphaTimeDist->GetBinContent(i)/a);
		hBetaTimeDist->SetBinContent(i,hBetaTimeDist->GetBinContent(i)/b);
	}
	
	for(int i=1; i<=NDefaultBins; i++) {
		hAlphaTimeCumulativeDist->SetBinContent(i,hAlphaTimeDist->Integral(1,i));
		hBetaTimeCumulativeDist->SetBinContent(i,hBetaTimeDist->Integral(1,i));
	}

	if (fSmooth) {
		hAlphaTimeCumulativeDist->Smooth(500);
		hBetaTimeCumulativeDist->Smooth(500);
		hAlphaTimeDist->Smooth(500);
		hBetaTimeDist->Smooth(500);
	}

	// compute gatti weights
	for(int i=1; i<=NDefaultBins; i++) {
		// non cumulative
		float a = hAlphaTimeDist->GetBinContent(i);
		float b = hBetaTimeDist->GetBinContent(i);
		float c = a+b;
		if ( c > 0.000001)
			GattiWeight[i-1] = (a-b)/c;
		else
			GattiWeight[i-1] = 0.;

		if ( (b > 1.E-10) && (a > 1.E-10) ) 
			LklMeanRatio[i-1] = TMath::Log(a/b);
		else {
			get_message(bx_message::error) << GetName() <<  " Mean value in histogram too small. Not enough samples or too many bins ! (1)" 
			                               << dispatch;
		}
		
		// cumulative
		a = hAlphaTimeCumulativeDist->GetBinContent(i);
		b = hBetaTimeCumulativeDist->GetBinContent(i);
		c = a+b;
		if ( c > 0.000001)
			GattiCumuWeight[i-1] = (a-b)/c;
		else
			GattiCumuWeight[i-1] = 0.;
		
		if ( (b > 1.E-10) && (a > 1.E-10) ) 
			LklMeanRatioCumu[i-1] = TMath::Log(a/b);
		else {
			get_message(bx_message::error) << GetName() << " Mean value in histogram too small. Not enough samples or too many bins ! (2)"
				<< dispatch;
		}
	}

	dump_gatti();

	fReady = true;
	return fReady;
}

// load classifiers from a file (default name available)
bool bx_alphabeta::load(const std::string& fname) {
	if (fReady) {
		get_message(bx_message::error) << GetName() << " You cannot load on a ready classifier" << dispatch;
		return false;
	}
	get_message(bx_message::info) << GetName() << " Loading alpha beta learning file " << fname << dispatch;
	FILE *fp=fopen(fname.c_str(),"rt");
	if (!fp)
		return false;
	char *s=0;
	size_t len=0;
	getline(&s,&len,fp);
	getline(&s,&len,fp);
	getline(&s,&len,fp);
	NDefaultBins = -1;
	if ( sscanf(s,"Histo size: %d",&NDefaultBins) != 1 ) {
		get_message(bx_message::error) << GetName() << " Invalid file " << fname << ". Missing histo size" << dispatch;
		fclose(fp);
		return false;
	}
	if (NDefaultBins < 0) {
		get_message(bx_message::error) << GetName() << " Wrong size in file " << fname << ". Size="<<NDefaultBins << dispatch;
		fclose(fp);
		return false;		
	}
	getline(&s,&len,fp);
	for(int i=0; i<=NDefaultBins; i++) {
		getline(&s,&len,fp);
		sscanf(s,"%f",&DefaultBinning[i]);
	}
	set_time_histogram(NDefaultBins,DefaultBinning);
	getline(&s,&len,fp);
	for(int i=0; i<NDefaultBins; i++) {
		getline(&s,&len,fp);
		sscanf(s,"%f",&GattiWeight[i]);
	}
	getline(&s,&len,fp);
	for(int i=0; i<NDefaultBins; i++) {
		getline(&s,&len,fp);
		sscanf(s,"%f",&GattiCumuWeight[i]);
	}
	getline(&s,&len,fp);
	for(int i=0; i<NDefaultBins; i++) {
		getline(&s,&len,fp);
		sscanf(s,"%f",&LklMeanRatio[i]);
	}
	getline(&s,&len,fp);
	for(int i=0; i<NDefaultBins; i++) {
		getline(&s,&len,fp);
		sscanf(s,"%f",&LklMeanRatioCumu[i]);
	}
	getline(&s,&len,fp);
	for(int i=0; i<NDefaultBins; i++) {
		getline(&s,&len,fp);
		float aa;
		sscanf(s,"%f",&aa);
		hAlphaTimeDist->SetBinContent(i+1,aa);
	}
	getline(&s,&len,fp);
	for(int i=0; i<NDefaultBins; i++) {
		getline(&s,&len,fp);
		float aa;
		sscanf(s,"%f",&aa);
		hBetaTimeDist->SetBinContent(i+1,aa);
	}
		
	fclose(fp);

	fReady = true;

	dump_gatti();

	return true;
}
	
// set a root event event for filtering
void bx_alphabeta::set(const BxEvent& ev, int ncluster) {
	if ( get_sample(ev,ncluster,&current_sample) )
		set_histos();
}

// set an echidna event event for filtering
void bx_alphabeta::set(const bx_echidna_event& ev, int ncluster) {
	if ( get_sample(ev,ncluster,&current_sample) ) 
		set_histos();
}

// convert the sample into histograms to be used by filters
void bx_alphabeta::set_histos() {
	// refer all times to the truncated mean
	// float moda;
	float mean = 0.; // current_sample.truncated_mean(moda);
	int size = current_sample.size();
	hHisto->Reset("ICE");
	hHistoCumu->Reset("ICE");
	for(int i=0; i<size; i++) {
		current_sample[i] -= mean;
		hHisto->Fill(current_sample[i]);
	}
	float integ = hHisto->Integral();
	for(int i=1; i<=size; i++) {
		hHisto->SetBinContent(i,hHisto->GetBinContent(i)/integ);
	}
	for(int i=1; i<=size; i++) {
		hHistoCumu->SetBinContent(i,hHisto->Integral(1,i));
		}
}


//
// ------- alpha-beta discriminating variables
//

// tail tot
float bx_alphabeta::tailtot(float timecut) {
	if (current_sample.size() == 0) {
		get_message(bx_message::error) << GetName() << " You must set an event before filtering!" << dispatch;
		return 0.;
	}
	int tail=0;
	int size = current_sample.size();
	for(int i=0; i<size; i++) {
		float t = current_sample[i];
		if (t>timecut)
			tail++;
	}
	return ((float)tail)/((float)size);
}

// compute gatti filter
// returns cumulative distribution weight (cumu)
// as reference the non cumulative one is also given
float bx_alphabeta::gatti(float& no_cumu) {
	if (current_sample.size()==0) {
		get_message(bx_message::error) << GetName() << " You must set an event before filtering!" << dispatch;
		return 0.;
	}
	if (!fReady) {
		get_message(bx_message::error) << GetName() << " You must learn first, or load a file!" << dispatch;
		return 0.;
	}
	float cumu = 0.;
	no_cumu = 0.;
	int size = hHisto->GetNbinsX();
	if ( size != NDefaultBins ) {
		get_message(bx_message::error) << GetName() << " Error in binning!! Size != NDefaultBins!!!!" << dispatch;
	}
	for(int i=0; i<size; i++) {
		float a = hHisto->GetBinContent(i+1);    // root bins start from 1
		no_cumu += GattiWeight[i]*a;
		a = hHistoCumu->GetBinContent(i+1);      // root bins start from 1
		cumu += GattiCumuWeight[i]*a;
	}

	return cumu;
}

float bx_alphabeta::lkl(float& no_cumu) {
	if (current_sample.size()==0) {
		get_message(bx_message::error) << GetName() << " You must set an event before filtering!" << dispatch;
		return 0.;
	}
	if (!fReady) {
		get_message(bx_message::error) << GetName() << " You must learn first, or load a file!" << dispatch;
		return 0.;
	}
	float cumu = 0.;
	no_cumu = 0.;
	int size = hHisto->GetNbinsX();
	for(int i=0; i<size; i++) {
		float a = hHisto->GetBinContent(i+1);    // root bins start from 1
		no_cumu += LklMeanRatio[i]*a;
		a = hHistoCumu->GetBinContent(i+1);      // root bins start from 1
		cumu += LklMeanRatioCumu[i]*a;
	}

	return cumu;
}

//
// kolmogorov test between sample and reference distribution
// the ratio between the alpha probability and the beta probability is returned
// this is done for both the time distributions and the cumulative time distributions
//
float bx_alphabeta::kolmo(float& no_cumu) {
	if (current_sample.size()==0) {
		get_message(bx_message::error) << GetName() << " You must set an event before filtering!" << dispatch;
		return 0.;
	}
	if (!fReady) {
		get_message(bx_message::error) << GetName() << " You must learn first, or load a file!" << dispatch;
		return 0.;
	} 
	
	int n = hAlphaTimeCumulativeDist->GetNbinsX();
	float *sh1=hHisto->GetArray();
	float *sh2=hHistoCumu->GetArray();
	float *sa1=hAlphaTimeCumulativeDist->GetArray();
	float *sb1=hBetaTimeCumulativeDist->GetArray();
	float *sa2=hAlphaTimeDist->GetArray();
	float *sb2=hBetaTimeDist->GetArray();

	// convert to double
	double h1[2000];
	double h2[2000];
	double a1[2000];
	double b1[2000];
	double a2[2000];
	double b2[2000];
    	for(int i=0; i<n; i++) {
		h1[i] = sh1[i];
		h2[i] = sh2[i];
		a1[i] = sa1[i];
		b1[i] = sb1[i];
		a2[i] = sa2[i];
		b2[i] = sb2[i];
	}

	// perform Kolmogorov test with the two reference distributions, cumulative and not
	float cumu_a = TMath::KolmogorovTest(n,a1,n,h1,"");
	float cumu_b = TMath::KolmogorovTest(n,b1,n,h1,"");
	float no_cumu_a = TMath::KolmogorovTest(n,a2,n,h2,"");
	float no_cumu_b = TMath::KolmogorovTest(n,b2,n,h2,"");

	// returns ratio as discriminating variable
	float cumu = cumu_a / cumu_b;
	no_cumu = no_cumu_a / no_cumu_b;
	
	return cumu;
}

float bx_alphabeta::knnr(int k) {
	if (current_sample.size()==0) {
		get_message(bx_message::error)<< GetName() << " You must set an event before filtering!" << dispatch;
		return 0.;
	}
	return 0.;
}

bool bx_alphabeta::save(const std::string& fname) {
	if (!fReady) {
		get_message(bx_message::error) << GetName() << " You must learn first! Nothing to save" << dispatch;
		return false;
	}
	FILE *fp=fopen(fname.c_str(),"wt+");
	
	fprintf(fp,"AlphaBeta classifier save file\n");
	time_t t = time(0);
	struct tm *ts = localtime(&t);
	char s[300];
	strftime(s,30,"%d-%m-%Y  %H:%M:%S",ts);
	fprintf(fp,"File saved: %s\n",s);
        fprintf(fp,"Histo size: %d\n",NDefaultBins);
	fprintf(fp,"Binning:\n");
	for(int i=0; i<=NDefaultBins; i++) 
		fprintf(fp,"%f\n",DefaultBinning[i]);
	fprintf(fp,"Gatti non-cumulative weights\n");
	for(int i=0; i<NDefaultBins; i++) 
		fprintf(fp,"%f\n",GattiWeight[i]);
	fprintf(fp,"Gatti cumulative weights\n");
	for(int i=0; i<NDefaultBins; i++) 
		fprintf(fp,"%f\n",GattiCumuWeight[i]);
	fprintf(fp,"Lkl non-cumulative mean values ratio\n");
	for(int i=0; i<NDefaultBins; i++) 
		fprintf(fp,"%f\n",LklMeanRatio[i]);
	fprintf(fp,"Lkl cumulative mean values ratio\n");
	for(int i=0; i<NDefaultBins; i++) 
		fprintf(fp,"%f\n",LklMeanRatioCumu[i]);
	fprintf(fp,"Alpha Time Distribution\n");
	for(int i=1; i<=NDefaultBins; i++)
		fprintf(fp,"%f\n",hAlphaTimeDist->GetBinContent(i));
	fprintf(fp,"Beta Time Distribution\n");
	for(int i=1; i<=NDefaultBins; i++)
		fprintf(fp,"%f\n",hBetaTimeDist->GetBinContent(i));
	fclose(fp);
	return true;
}

void bx_alphabeta::exportsvm(const BxEvent& ev, int ncluster) {
	if ( get_sample(ev,ncluster,&current_sample) ) {
		set_histos();
		
		float area = hHisto->Integral();
		float charge = ev.GetLaben().GetCluster(ncluster-1).GetCharge();
		if (area > 0.001 && charge>60. && charge < 2000. ) {

			if (svm_energy_file==0)
				svm_energy_file=fopen("svm_energies.dat","a+");

			if (svm_sample_file==0)
				svm_sample_file=fopen("svm_sample.dat","a+");

			fprintf(svm_energy_file,"%f\n",charge);
			fprintf(svm_sample_file,"+1 ");                     // beta type +1
			for(int j=1; j<=NDefaultBins; j++) {
				float val = hHisto->GetBinContent(j)/area;
				fprintf(svm_sample_file,"%d:%06.4f ",j,val);
			}
			fprintf(svm_sample_file,"\n");
		}
	}
}

void bx_alphabeta::exportsvm_mc(const BxEvent& ev, int ncluster) {
	if ( get_sample(ev,ncluster,&current_sample) ) {
		set_histos();

		float area = hHisto->Integral();
		float charge = ev.GetLaben().GetCluster(ncluster-1).GetCharge();
		int geneb_num = ev.GetMcTruth().GetFrame(0).GetEventId();

		if (area > 0.001 && charge>60. && charge < 2000. ) {

			if (svm_energy_file==0)
				svm_energy_file=fopen("svm_energies.dat","a+");

			if (svm_sample_file==0)
				svm_sample_file=fopen("svm_sample.dat","a+");

			if (svm_alpha_file==0)
				svm_alpha_file=fopen("svm_alpha_mc.dat","a+");

			if (svm_beta_file==0)
				svm_beta_file=fopen("svm_beta_mc.dat","a+");

                        // charge file
			fprintf(svm_energy_file,"%f\n",charge);
			
			// sample file
			fprintf(svm_sample_file,"0 ");
			for(int j=1; j<=NDefaultBins; j++) {
				float val = hHisto->GetBinContent(j)/area;
				fprintf(svm_sample_file,"%d:%06.4f ",j,val);
			}
			fprintf(svm_sample_file,"\n");

			// FIXME: trucco temporaneo per separare i BiPo
			// alpha file
			if ( geneb_num%2 == 0 ) {
				fprintf(svm_alpha_file,"-1 ");                     
				for(int j=1; j<=NDefaultBins; j++) {
					float val = hHisto->GetBinContent(j)/area;                                           
					fprintf(svm_alpha_file,"%d:%06.4f ",j,val);                                       
				} 
				fprintf(svm_alpha_file,"\n");
			}
			else {

				// beta file
				fprintf(svm_beta_file,"+1 ");                     
				for(int j=1; j<=NDefaultBins; j++) {
					float val = hHisto->GetBinContent(j)/area;                                           
					fprintf(svm_beta_file,"%d:%06.4f ",j,val);                                       
				} 
				fprintf(svm_beta_file,"\n");
			}
		}
	}
}

