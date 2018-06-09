// 
// class for alpha-beta discrimination
// MP 22-02-07
// Implemented methods:
//

#ifndef _BX_ALPHABETA_HH_
#define _BX_ALPHABETA_HH_
#ifndef CYCLE_NUMBER
#define CYCLE_NUMBER 18
#endif

#include <vector>


class bx_echidna_event;
class TH1F;

#include "event/BxEvent.hh"
#include "interface/messenger.hh"
#include "bx_sample.hh"

class bx_alphabeta : public message_client {
public:
	bx_alphabeta(const std::string& name = "default");

	virtual ~bx_alphabeta();

	// set methods
	// default is used
	void set_time_histogram(int nbins=0, float* binning=0);
	
	// add methods: they must be called to prepare learning
	void add(const std::string& cname, const bx_echidna_event&, int ncluster=1);
	void add(const std::string& cname, const BxEvent&, int ncluster=1);

	// learning for all classifiers
	bool learn();

	// load classifiers from a file (default name available)
	bool load(const std::string& = "");
	
	// save classifiers to a file
	bool save(const std::string&);

	// set an event for filtering
	void set(const BxEvent&, int ncluster=1);                   // cluster number 1 based in root
	void set(const bx_echidna_event&, int ncluster=0);          // cluster number 1 based in echidna
	
	// alpha-beta discriminating variables

	// Tail over total: tail definition is parameter
	float tailtot(float timecut);

	// Gatti filter. Return value using cumulative distr. Reference is not cumulative one. 
	float gatti(float&);

	// Likelihood ratio assumung poisson statistics. Return value using cumulative distr. Reference is not cumulative one. 
	float lkl(float&);
	
	// Kolmogorov test between histograms. Return value using cumulative distr. Reference is not cumulative one. 
	float kolmo(float&);

	// Nearest neighbor tecnique. Return value using cumulative distr. Reference is not cumulative one. 
	float knnr(int k);
	
	// write current event for performing svm analysis externally
	// the sample and the energy are written in separate files
	// useful for svm development and studies
	void exportsvm(const BxEvent& ev, int ncluster = 1);
	
	// same for mc bipo events
	void exportsvm_mc(const BxEvent& ev, int ncluster = 1);
	
	// Some public histograms for monitoring
	TH1F *hMeanAlpha;
	TH1F *hMeanBeta;
	TH1F *hModaAlpha;
	TH1F *hModaBeta;

	TH1F *hAlphaTimeDist;
	TH1F *hBetaTimeDist;
	TH1F *hAlphaTimeCumulativeDist;
	TH1F *hBetaTimeCumulativeDist;

	bool isReady() const {return fReady;}

	const std::string& GetName() const {return fName;}

	void SetSmooth(bool val=true) {fSmooth=val;}
	
private:
	//ClassDef(bx_alphabeta,CYCLE_NUMBER)

	// scratch histograms used for filtering
	TH1F *hHisto;
	TH1F *hHistoCumu;

	// fill histograms and normalize them correctly for filtering
	void set_histos(); 

	// alpha and beta good samples
	std::vector<bx_sample*> fAlphaSamples;
	std::vector<bx_sample*> fBetaSamples;

	// function that convert the events into bx_sample
	bx_sample* get_sample(const bx_echidna_event& ev, int nclus, bx_sample* sampl=0);
	bx_sample* get_sample(const BxEvent& ev, int nclus, bx_sample* sampl=0);

	// add a sample to the learning list
	void add(const std::string& cname, bx_sample*);

	// histogram binning structure
	int NDefaultBins;
	float DefaultBinning[5000];
	
	// true if classifier ready for filtering
	bool fReady;
	
	// parameter to correct for event position or not
	bool fPositionCorrection;

	// scratch sample to avoid new and delete
	bx_sample current_sample;
	

	// gatti weights
	float GattiWeight[1000];          // non cumulative distribution weights
	float GattiCumuWeight[1000];      // cumulative distribution weights
	float LklMeanRatio[1000];         // log ratio of mean values for LKL
	float LklMeanRatioCumu[1000];     // same for cumulative distribution
	
	void dump_gatti(std::ostream& os = std::cout);

	int fErrorCounter;
	bool fSmooth;

	FILE *svm_energy_file;   //! do not store on file
	FILE *svm_sample_file;   //! do not store on file
	FILE *svm_alpha_file;    //! do not store on file
	FILE *svm_beta_file;     //! do not store on file

	std::string fName;       // this object name
};


#endif
