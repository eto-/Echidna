/* BOREXINO Reconstruction program
 *
 * Author: Stefano Davini <Stefano.Davini@ge.infn.it>
 * Maintainer: Stefano Davini <Stefano.davini@ge.infn.it>
 *
 * Calculation of several positron(c11, oPs)/electron discrimination variables
 *
 */

#include "bx_pid_positron.hh"
#include "bx_echidna_event.hh"
#include "barn_interface.hh"
#include "bx_dbi.hh"
#include "db_run.hh"
#include <algorithm>
#include <TFile.h>
#include <TH2D.h>
#include <TH1D.h>
#include <iostream>

// ctor
bx_pid_positron::bx_pid_positron (): bx_base_module("bx_pid_positron", bx_base_module::main_loop) {
  require_event_stage (	bx_detector::laben, bx_base_event::reconstructed); 
  require_trigger_type (bx_trigger_event::neutrino);
  require_trigger_type (bx_trigger_event::neutron);
}


// module interface
void bx_pid_positron::begin () {

	gatti_weights_loaded = false;
	bool noerr = true; // no error occurred

	nnhits_min  = get_parameter("nnhits_min").get_int();
	nnhits_max  = get_parameter("nnhits_max").get_int();
	ene_bins    = nnhits_max-nnhits_min;
	/*
	 * ene bins is the number of energy bins for the gatti shapes: 1 bin = 1 hit
	 */
	time_max    = get_parameter("time_max").get_int();
	smear_nops  = get_parameter("smear_nops").get_int();
	smear_ops   = get_parameter("smear_ops").get_int();
	smear_beta  = get_parameter("smear_beta").get_int();
	smear_c11   = get_parameter("smear_c11").get_int();
	sum_min     = get_parameter("sum_min").get_float();

	const std::string shapes_fname = get_parameter("shapes_file").get_string ();
	TFile* f_shapes = new TFile(shapes_fname.c_str(), "READ");
	
	if ( (!f_shapes) || (f_shapes->IsZombie()) || (!f_shapes->IsOpen()) ) {
		get_message(bx_message::error) << "Couldn't find shapes file for use with the Gatti filter." << dispatch;
		return; // don't use gatti filter
	}
	else get_message(bx_message::info) << "Shapes file " << shapes_fname << " successfully opened."<< dispatch;

	noerr *= load_shape(f_shapes, "nops", shape_nops, smear_nops);
	noerr *= load_shape(f_shapes,  "ops", shape_ops,  smear_ops);
	noerr *= load_shape(f_shapes, "beta", shape_beta, smear_beta);
	noerr *= load_shape(f_shapes,  "c11", shape_c11,  smear_c11);
	f_shapes->Close("R");
	delete f_shapes;
	if (noerr){
		get_message(bx_message::info) << "Shapes file " << shapes_fname << " successfully closed."<< dispatch;
	}
	else{
		get_message(bx_message::error) << "Some error occurred when loading shapes " << dispatch;
		return;
	}

	noerr*=compute_gatti_weights(gatti_weight_ops_beta, shape_ops, shape_beta);
	noerr*=compute_gatti_weights(gatti_weight_c11_beta, shape_c11, shape_beta);
	noerr*=compute_gatti_weights(gatti_weight_ops_nops, shape_ops, shape_nops);

	mdelete(shape_nops);
	mdelete(shape_ops);
	mdelete(shape_beta);
	mdelete(shape_c11);

	if (noerr){
		gatti_weights_loaded = true;
		get_message(bx_message::info) << "Gatti weights computed."<< dispatch;
	}
	else{
		get_message(bx_message::error) << "Some error occurred when computing Gatti weights " << dispatch;
		return;
	}
	
	hSample = new TH1F("pid_positron_time_sample", "", time_max, 0., time_max);
}

// doit
bx_echidna_event* bx_pid_positron::doit (bx_echidna_event *ev) {

	
	const int nc = ev->get_laben().get_nclusters();
	if (nc <= 0) return ev;


	// loop on clusters and compute gatti variables
	for(int iclus=0; iclus<nc; iclus++) {
		
		// get echidna positron_cluster
		bx_laben_positron_cluster& ps = ev->get_laben().get_positron_cluster(iclus);

		// set variables to invalid values
		ps.f4_gatti_ops_beta     = -100.;
		ps.f4_gatti_c11_beta     = -100.;
		ps.f4_gatti_ops_nops     = -100.;
		
		if (!gatti_weights_loaded) continue;

		const int nhits     = ev->get_laben().get_cluster(iclus).get_clustered_nhits();
		const int nlivepmts = ev->get_laben().get_n_live_pmts();
		// TODO: better normalization (invalid pmts??)
		const float nnhits = nhits*2000./nlivepmts;

		if (nnhits < 300) continue;

		/*
		 * nnhits_min, nnhits_max should be larger than the c11/pep region, in order to catch event with bad energy resolution
		 */

		const int ebin =  ( (nnhits < nnhits_min) || (nnhits >= nnhits_max) ) ? ( (nnhits < nnhits_min) ? 0 : ene_bins - 1  )  : int(nnhits - nnhits_min);

		//if ( (nnhits < nnhits_min) || (nnhits >= nnhits_max) )
		//	continue; 
		//const int ebin = int(nnhits - nnhits_min);

		hSample->Reset();

		// Fill histo of times (skip the first 2 and shift to avoid noise hit effects)
		const float shift_time = ps.get_rec_hit(2).get_time();
		for (int i=2; i<nhits; i++) hSample->Fill( ps.get_rec_hit(i).get_time()-shift_time );

		// Normalize histo
		const Double_t integral = hSample->Integral(1, time_max);
		if (integral>0)
			hSample->Scale(1./integral);
		else
			get_message(bx_message::warn) << "Sample histogram area too small " << dispatch;


		double gatti_ops_beta = 0.;
		double gatti_c11_beta = 0.;
		double gatti_ops_nops = 0.;
		for (int it=0; it<time_max; it++) {
			const double sample = hSample->GetBinContent(hSample->FindBin(it));
			gatti_ops_beta  +=  sample * gatti_weight_ops_beta[ebin][it];
			gatti_c11_beta  +=  sample * gatti_weight_c11_beta[ebin][it];
			gatti_ops_nops  +=  sample * gatti_weight_ops_nops[ebin][it];
		}


		ps.f4_gatti_ops_beta     = gatti_ops_beta;
		ps.f4_gatti_c11_beta     = gatti_c11_beta;
		ps.f4_gatti_ops_nops     = gatti_ops_nops;

	} // end of loop on clusters

	return ev;
}

void bx_pid_positron::end () {

	if (gatti_weights_loaded){
		mdelete(gatti_weight_ops_beta);
		mdelete(gatti_weight_c11_beta);
		mdelete(gatti_weight_ops_nops);
		delete hSample;
	}
}


bool bx_pid_positron::load_shape(TFile *f, const std::string & histo_name, double** & shape, int smear){

	// allocation of rec-times shapes matrix -> 1st index energy, 2nd index time
	bool noerror = allocate(shape);
	if (!noerror){
		get_message(bx_message::error) << "No memory for allocating " << histo_name << " matrix"<< dispatch;
		return false;
	}

	// 2D histograms: x axis: energy (1 bin = 1 nhits), y axis: rec-time at time t (1 bin = 1 ns)
	TH2D* h2shape = (TH2D*)f -> Get(histo_name.c_str()) -> Clone();
	if (!h2shape){
		get_message(bx_message::error) << "Rec-times Shapes " << histo_name << " not found"<< dispatch;
		return false;
	}

	for (int ie = 0; ie < ene_bins; ie++){
		const int proj_start = nnhits_min + ie - smear;
		const int proj_stop  = nnhits_min + ie + smear;

		TH1D* h1shape = (TH1D*) h2shape->ProjectionY("h1_ene_proj", proj_start+1, proj_stop); 
		if (!h1shape) return false;

		// normalization; 1 because Integral() arguments are bin
		const Double_t integral = h1shape->Integral(1, time_max);
		if (integral==0){
			get_message(bx_message::warn) << "Rec-times Shapes " << histo_name << " empty beetween bins "<< proj_start <<" and "<< proj_stop << dispatch;
			continue;
		}
		h1shape->Scale(1./integral);
	
		for (int it = 0; it< time_max; it++){
			shape[ie][it] = h1shape->GetBinContent(h1shape->FindBin(it));
		}
		h1shape->Delete();
	}

	h2shape->Delete();
	get_message(bx_message::info) << "Shapes " << histo_name << " successfully loaded."<< dispatch;
	return true;
}

bool bx_pid_positron::compute_gatti_weights(double** &w, double** p, double** e){
	bool noerror= allocate(w);
	if (!noerror) return false;
	for (int ie=0; ie < ene_bins; ie++){
		for (int it=0; it < time_max; it++){
			const double pp = p[ie][it];
			const double ee = e[ie][it];
			const double sum  = pp + ee; 
			const double diff = pp - ee;
			w[ie][it] = (sum > sum_min) ? diff/sum : 0 ;
		}
	}
	return true;
}

bool bx_pid_positron::allocate(double** &shape){
	shape = new double* [ene_bins];
	if (!shape) return false;   
	for (int i=0; i<ene_bins; i++){
		shape[i] = new double [time_max];
		if (!shape[i]) return false;
	}
	return true;
}

void bx_pid_positron::mdelete(double** &shape){
	for (int i=0; i<ene_bins; i++){
		delete (shape[i]);	
	}
	delete [] (shape);
}

/*
 *  $Log: bx_pid_positron.cc,v $
 *  Revision 1.6  2011/04/13 14:40:57  davini
 *  enlarged range: nnhits >= 300; using ebin=0 for low nnhits events, ene_bin=ene_bins-1 for high nnhits events
 *
 *  Revision 1.5  2011-02-23 10:30:30  davini
 *  smearing implemented;
 *
 *  Revision 1.4  2011-02-21 10:04:59  davini
 *  small changes; logs;
 *
 *  Revision 1.3  2011-02-20 15:47:48  davini
 *  Stable version. Under Alessandro's suggestion, the real lord of C++ and master of pointers and references, double** -> double** & in allocation, load, compute delete function, in order to avoid allocation/passing problems. Debug done.
 *
 *  Revision 1.2  2011-02-19 13:32:33  davini
 *  Basic features added; the module writes in the event; extensive testingstill needed;
 *
 *  Revision 1.1  2011-02-18 15:30:59  davini
 *  added module
 *
 */
