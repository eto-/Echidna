/* Port of Mach4 ProcessAlphaBeta Madule
 * Alvaro Chavarria, Princeton
 */

#ifndef _BX_PID_AB_MACH4_H
#define _BX_PID_AB_MACH4_H

#include "bx_base_module.hh"
#include "bx_rec_general.hh"

#include "bx_alphabeta.hh"

#include <TFile.h>

class bx_pid_ab_mach4: public bx_base_module {
public:
    bx_pid_ab_mach4 ();
    virtual ~bx_pid_ab_mach4 () {}
	
    virtual void begin ();
	
    virtual bx_echidna_event* doit (bx_echidna_event *ev);
	
    virtual void end ();
    
private:
	
	std::vector<std::pair<std::string, std::vector<double> > > gatti_weights;
	///< Vector of histograms of Gatti weights (for several possible filters.)
	///< The first element in the pair is the name of the filter.
	
	bool skip_gatti_filter;
	///< Boolean that sets if Gatti parameter should be calculated.
	
	double gatti_hist_min;
	///< Start time for Gatti histograms.
	double gatti_hist_max;
	///< End time for Gatti histograms.
	double gatti_bin_size;
	///< Gatti histogram bin size.
	bool gatti_weights_loaded;
	///< Flag showing whether we have Gatti weights loaded.
	
	double xpmt[2241];
	double ypmt[2241];
	double zpmt[2241];
	///< Variables to keep the position of the PMTs for the TOF subtraction.
	
	double GetGattiWeight(size_t i, double t, bool * = NULL);
	///< Returns the Gatti weight associated with the given time
	///< using the @a i th histogram of Gatti weights.
	void DoGattiFilter(bx_echidna_event *ev, int cluster, double ref_time, std::vector<double>*);
	///< Use the Gatti filters on the given cluster, given the peak time offset.
	void LoadGattiShapes(TFile *f, const std::string &ahist, const std::string &bhist, const std::string &name);
	///< Load the Gatti shapes corresponding to the given histograms from
	///< the given ROOT file.
		
};

#endif
