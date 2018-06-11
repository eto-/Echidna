/* BOREXINO Reconstruction program
 *
 * Author: Stefano Davini <Stefano.Davini@ge.infn.it>
 * Maintainer: Stefano Davini <Stefano.Davini@ge.infn.it>
 *
 * Calculation of several positron(c11, oPs)/electron discrimination variables
 *
 */

#ifndef _BX_PID_POSITRON_H
#define _BX_PID_POSITRON_H

#include "bx_base_module.hh"
#include "bx_rec_general.hh"
#include "bx_alphabeta.hh"
#include <TFile.h>

class bx_pid_positron: public bx_base_module {
  public:
	bx_pid_positron ();
	virtual ~bx_pid_positron () {}

	virtual void begin ();

	virtual bx_echidna_event* doit (bx_echidna_event *ev);

	virtual void end ();
    
  private:

	bool load_shape(TFile *f, const std::string &hist_name,  double** &rec_time_shapes, int32_t smear);
	/*
	 * Fill the rec_time_shapes (energy-time) matrix with the content of the hisogram named hist_name in TFile f.
	 * The rec-times are normalized within the function.
	 * The matrix rec_time_shapes is allocated in the function;
	 * smear is used to take the content of near energy bins (entry correspondend to energy n will be filled with [n-smear, n+smear]);
	 */
	bool compute_gatti_weights(double** &weights, double** positron_shape, double** electron_shape);
	/*
	 * Fill the (energy, time) matrix (gatti) weigths using the  matrix positron_shape and electron_shape.
	 * The matrix weights is allocated in the function.
	 */
	bool allocate(double** &shape);
	/*
	 * Allocate the energy-time matrix;
	 * returns false if error occurred
	 */
	void mdelete(double** &shape);
	/*
	 * Free the memory of the energy-time matrix
	 */

	bool gatti_weights_loaded;
    
	int32_t nnhits_min;
	int32_t nnhits_max;
	int32_t ene_bins;
	int32_t smear_ops;
	int32_t smear_nops;
	int32_t smear_c11;
	int32_t smear_beta;
	int32_t time_max;
	float sum_min;

	double** shape_nops;
	double** shape_ops;
	double** shape_c11;
	double** shape_beta;

	double** gatti_weight_ops_beta;
	double** gatti_weight_c11_beta;
	double** gatti_weight_ops_nops;

	TH1F* hSample;
    
};

#endif
/*
 * $Log: bx_pid_positron.hh,v $
 * Revision 1.5  2011/02/23 10:30:30  davini
 * smearing implemented;
 *
 * Revision 1.4  2011-02-20 15:47:47  davini
 * Stable version. Under Alessandro's suggestion, the real lord of C++ and master of pointers and references, double** -> double** & in allocation, load, compute delete function, in order to avoid allocation/passing problems. Debug done.
 *
 * Revision 1.3  2011-02-19 13:32:33  davini
 * Basic features added; the module writes in the event; extensive testingstill needed;
 *
 * Revision 1.2  2011-02-18 15:36:04  ddangelo
 * debugging
 *
 * Revision 1.1  2011-02-18 15:30:59  davini
 * added module
 *
 */


