/* Port of Mach4 ProcessAlphaBeta Madule
 * Alvaro Chavarria, Princeton
 */

#include "bx_pid_ab_mach4.hh"
#include "bx_echidna_event.hh"
#include "barn_interface.hh"
#include "bx_alphabeta.hh"
#include "db_channel.hh"
#include "bx_dbi.hh"
#include "db_run.hh"
#include <algorithm>
#include "TVector3.h"
#include "mach4/BxMath.h"

#define SPEED_OF_LIGHT 0.299792458 /* m/ns */

// ctor
bx_pid_ab_mach4::bx_pid_ab_mach4 (): bx_base_module("bx_pid_ab_mach4", bx_base_module::main_loop) {
	require_event_stage (	bx_detector::laben, bx_base_event::reconstructed); 
	require_trigger_type (bx_trigger_event::neutrino);
	require_trigger_type (bx_trigger_event::neutron);
	gatti_weights_loaded = false;
}

void bx_pid_ab_mach4::LoadGattiShapes(TFile *f, const std::string & ahist, const std::string & bhist, const std::string & name)
{
	TH1F *alpha_shape, *beta_shape;
	
	alpha_shape = (TH1F *)f -> Get(ahist.c_str()) -> Clone();
	beta_shape = (TH1F *)f -> Get(bhist.c_str()) -> Clone();
	
	alpha_shape -> Scale(1 / alpha_shape -> Integral());
	beta_shape -> Scale(1 / beta_shape -> Integral());
	
	int min_bin, max_bin;
	TAxis *tmp = alpha_shape -> GetXaxis();
	
	min_bin = tmp -> GetFirst();
	max_bin = tmp -> GetLast();
	
	// we overwrite these every time, so they'd better be the same for all histos
	gatti_bin_size = tmp -> GetBinWidth(min_bin);
	gatti_hist_min = tmp -> GetBinLowEdge(min_bin);
	gatti_hist_max = tmp -> GetBinUpEdge(max_bin);
	
	gatti_weights.push_back(std::pair<std::string, std::vector<double> >());
	std::pair<std::string, std::vector<double> > &ref = gatti_weights.back();
	std::vector<double> &weights = ref.second;
	ref.first = name;
	
	weights.resize(int((gatti_hist_max - gatti_hist_min) / gatti_bin_size));
	for (int i = min_bin; i <= max_bin; i++) {
		double a = alpha_shape -> GetBinContent(i);
		double b = beta_shape -> GetBinContent(i);
		weights[i - min_bin] = (a - b) / (a + b);
	}
	
	alpha_shape->Delete();
	beta_shape->Delete();
}

void bx_pid_ab_mach4::begin () {
	
	TFile *f;
	
	gatti_weights_loaded = false;
	gatti_weights.clear();
	
	std::string shapes_fname = get_parameter("shapes_file").get_string ();
	f = new TFile( shapes_fname.c_str(), "READ");
	
	//Read from config file
	skip_gatti_filter = get_parameter("skip_gatti_filter").get_bool();
	
	if (!f || !f -> IsOpen()) {
		get_message(bx_message::error) << "Couldn't find alpha/beta shapes file for use with the Gatti filter." << dispatch;
		return;// we still don't want to fail, just don't use the Gatti filter.
	}
	
	gROOT -> cd();
	
	LoadGattiShapes(f, "alpha_bipo", "beta_bipo", "bipo");
	LoadGattiShapes(f, "alpha_pos", "beta_pos", "pos");
	LoadGattiShapes(f, "alpha_ttr", "beta_ttr", "ttr");
	LoadGattiShapes(f, "alpha_skew", "beta_skew", "skew");
	gatti_weights_loaded = (gatti_weights.size() != 0);
	f->Close("R");
	delete f;
	
	//Now that the time of flight subtraction will be done in this module we need to get the PMT positions
  	for(size_t channel = 1; channel <= constants::laben::channels; channel++)
	{
		const db_channel_laben& channel_info = dynamic_cast<const db_channel_laben&>(bx_dbi::get()->get_channel(channel));
		
		xpmt[channel] = channel_info.pmt_x();
		ypmt[channel] = channel_info.pmt_y();
		zpmt[channel] = channel_info.pmt_z();
	}
	
	
}

void bx_pid_ab_mach4::end () {
}

double bx_pid_ab_mach4::GetGattiWeight(size_t i, double t, bool *flag)
{ // careful: no checks to see if gatti_weights is valid here!
	if (i >= gatti_weights.size()) return 0;
	std::vector<double> &weights = gatti_weights[i].second;
	int j = int((t - gatti_hist_min) / gatti_bin_size);
	if (j < 0 || j >= int(weights.size())) {
		if (flag) *flag = false;
		return 0;
	}
	if (flag) *flag = true;
	return weights[j];
}

void bx_pid_ab_mach4::DoGattiFilter(bx_echidna_event *ev, int iclus, double ref_time, std::vector<double>* times)
{
	if (!gatti_weights_loaded) return; // nothing to do without the weights...
	
	//bx_laben_ab_mach4_cluster& ab = ev->get_laben().get_ab_mach4_cluster(iclus);
	
	for (size_t i = 0; i < gatti_weights.size(); i++) {
		
		bool is_valid_bin;
		double gatti_index = 0;
		unsigned nhits = 0;
		
		size_t n = times->size();	
		
		for (unsigned j = 0; j < n; j++) {
		
			double t = times->at(j) - ref_time;
			gatti_index += GetGattiWeight(i, t, &is_valid_bin);
			if (is_valid_bin) nhits++;
		}
		
		if (nhits != 0) gatti_index /= nhits;
		else get_message(bx_message::debug) << "No hits inside pulse window, event " << ev -> get_event_number() << ", cluster " << iclus << "." << dispatch;
		
//		ab.v_gatti_mach4[i] = gatti_index;
		
	}
}

bx_echidna_event* bx_pid_ab_mach4::doit (bx_echidna_event *ev) {
	
	// get number of clusters
	int nc = ev->get_laben().get_nclusters();
	
	if (nc <= 0)
		return ev;
	
	// loop on clusters and compute alpha-beta variables
	for(int iclus=0; iclus<nc; iclus++) {
		
		// get echidna ab_cluster
//		bx_laben_ab_mach4_cluster& ab = ev->get_laben().get_ab_mach4_cluster(iclus);
		
		// set variables to invalid values
//		for(int index=0; index<17; index++) ab.v_tailtot_mach4[index] = -100.;
//		for(int index=0; index<4; index++) ab.v_gatti_mach4[index] = -100.;		
		
/*		ab.f4_peak_mach4 = -100.;
		ab.f4_mean_mach4 = -100.;
		ab.f4_rms_mach4 = -100.;
		ab.f4_skew_mach4 = -100.;
		ab.f4_kurt_mach4 = -100.;
*/		
		// Re-calculate hit time distribution parameters now that we
		// in principle know where the hit occurred.
		
		size_t n = ev->get_laben().get_cluster(iclus).get_clustered_nhits();
		std::vector<double> times;
		// <-- an array to store emission times of photons relative to event time
		
		for (size_t i = 0; i < n; i++) {

			//get time of clustered hit and channel number.
			double t = ev->get_laben().get_cluster(iclus).get_clustered_hit(i).get_time();
			size_t channel = ev->get_laben().get_cluster(iclus).get_clustered_hit(i).get_decoded_hit().get_raw_hit().get_logical_channel();

			//perform TOF subtraction to obtain emission time.
			double n_c = ev->get_laben().get_cluster(iclus).get_position_mach4().get_ref_index() / SPEED_OF_LIGHT;
						
			TVector3 ev_pos (ev->get_laben().get_cluster(iclus).get_position_mach4().get_x(), ev->get_laben().get_cluster(iclus).get_position_mach4().get_y(), ev->get_laben().get_cluster(iclus).get_position_mach4().get_z());
			TVector3 pmt_pos (xpmt[channel], ypmt[channel], zpmt[channel]);
			
			t -=  n_c * (ev_pos - pmt_pos).Mag();
			
			//sanity check
			if (std::fabs(t) > 1e4 || std::isnan(t) || std::isinf(t))
				// something is wrong; most probably position reconstruction didn't
				// converge and hence the position module set an invalid emission time
				continue;
			
			times.push_back(t);
		}
		
		n = times.size();
		if (!n) continue;
		std::sort(times.begin(), times.end());
		
		// Calculate peak time of cluster relative to first hit.
		double t_peak = BxMath::get_peak(times.begin(), times.end());
//		ab.f4_peak_mach4 = t_peak;
		
		// Calculate other parameters of the emission time distribution.
		std::vector<double> moments =
		BxMath::get_distribution_moments(times.begin(), times.end(), 4);
		
/*		ab.f4_mean_mach4 = moments[1];
		ab.f4_rms_mach4 = std::sqrt(moments[2]);
		ab.f4_skew_mach4 = moments[3] / (std::sqrt(moments[2]) * moments[2]);
		ab.f4_kurt_mach4 = moments[4] / (moments[2] * moments[2]) - 3.;
*/		
		// Calculate tail-to-total ratios
		for (size_t bin = 0; bin < 17; bin++) {
			// Fill in tail/total's for tail offsets that are multiples of 5 ns in
			// the range [30 ns, 110 ns], as these cover either side of the optimum
			// figure of merit, which is achieved by an offset around 45-55 ns.
			size_t offset = 5 * (bin + 6);
			double ttr = 0;
			
			for (int i = (int)times.size() - 1; i >= 0; i--) {
				if (times[i] - t_peak >= (double)offset)
					ttr++;
				else
					break;
			}
			ttr /= n;
			
//			ab.v_tailtot_mach4[bin] = ttr;
		}
		
		if ( !skip_gatti_filter ) DoGattiFilter(ev, iclus, t_peak, &times);
		
	}
	
	return ev;
	
}
