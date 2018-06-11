/* BOREXINO Reconstruction program
 *
 * Author: Alessandro Razeto <Alessandro.Razeto@ge.infn.it>
 * Maintainer: Alessandro Razeto <Alessandro.Razeto@ge.infn.it>
 *
 * $Id: bx_splitting_filter.hh,v 1.9 2008/02/21 11:22:00 razeto Exp $
 *
 * Splitting module: filter the cluster signal with a sample 
 * using FFT for evaluating the cross-correlation
 * 
 */
#ifndef _BX_SPLITTING_FILTER_H
#define _BX_SPLITTING_FILTER_H

#include "bx_base_module.hh"
#include "bx_rec_general.hh"
#include "cmap.hh"
#include <fftw3.h>

class TH1F;
class TH2F;
class TSpectrum;

class bx_splitting_filter: public bx_base_module {
  public:
    bx_splitting_filter ();
    virtual ~bx_splitting_filter () {}

    virtual void begin ();
    virtual bx_echidna_event* doit (bx_echidna_event *ev);
    virtual void end ();
  private:
      // FFT buffers, plans and normalizations
    double *shape_in, *shape_out;
    double *sample_F, *shape_in_F, *shape_out_F;
    fftw_plan plan_in, plan_out;
    double f8_sample_area, f8_sample_square_area;

      // Parameters
    static const int32_t max_npeaks = 50;
    int32_t i4_window_width, i4_zero_bin;
    int32_t i4_sample_ramp_lenght, i4_sample_width;
    float f4_sample_exp_tao;
    int32_t i4_charge_threshold, i4_time_threshold;
    enum peak_mode_t {
      tspectrum_simple,
      tspectrum_highres,
    } peak_mode;
    std::cmap<std::string, peak_mode_t> peak_mode_names;

      // Histograms
    TH2F *peak_time_vs_energy;

      // TSpectrum stuff
    TSpectrum *S;
    TH1F *h_shape_out;
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,0,0)
    double *tmp_in, *spectrum;
#else
    float *tmp_in, *spectrum;
#endif

    template<typename FLOAT> int32_t search_peak_base (FLOAT *v, int32_t& peak, int32_t low_limit);
    void ignore_peak (const std::string& msg, int32_t current_peak, float t, float charge, float duration);
    void reset_histograms (int32_t evnum);
};

#endif
/*
 * $Log: bx_splitting_filter.hh,v $
 * Revision 1.9  2008/02/21 11:22:00  razeto
 * Check on NaN added on fft output
 *
 * Revision 1.8  2007-04-13 13:56:14  razeto
 * BE quiter
 *
 * Revision 1.7  2007-03-29 14:27:10  razeto
 * Better splitting algorithm
 *
 * Revision 1.6  2005/09/19 15:19:23  razeto
 * Removed a useless histogram
 *
 * Revision 1.5  2004/12/19 17:03:05  razeto
 * Added generation of another debugging histogram
 *
 * Revision 1.4  2004/12/17 15:40:34  razeto
 * Indroduced peak finding algorithm.
 * The identification of peaks work but need to be tested; some magic numbers
 * still hardcoded.
 * Waiting for the presence in the event of tagging fields.
 *
 * Revision 1.3  2004/12/13 12:39:17  razeto
 * Changed to use fftw version 2 instead of version 3 (not available on cluster)
 *
 * Revision 1.2  2004/12/10 13:53:47  razeto
 * A first working filter
 *
 * Revision 1.1  2004/12/09 17:10:06  razeto
 * First code (not working)
 *
 */
