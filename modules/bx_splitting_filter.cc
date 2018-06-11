/* BOREXINO Reconstruction program
 *
 * Author: Alessandro Razeto <Alessandro.Razeto@ge.infn.it>
 * Maintainer: Alessandro Razeto <Alessandro.Razeto@ge.infn.it>
 *
 * $Id: bx_splitting_filter.cc,v 1.26 2009/10/23 12:40:46 razeto Exp $
 *
 * Implementation of bx_splitting_filter; some note on the optimal filter algorithm.
 * 1) naming convenction: lets define SAMPLE the reference shaple, SHAPE_IN the hit 
 * time distribution of the cluster (with 1ns binning) and SHAPE_OUT the results from
 * the filtering procedure (with the same binning).
 * 2) SAMPLE will be time reversed from the real shape; this is done to keep
 * normalization (point 3).
 * 3) normalization: consider f and g where f is the SAMPLE and g is the SHAPE_IN;
 * suppose Max(f) = 1 and define N1 = I(f) and N2 = I(f * f) as the integral of f 
 * and f^2; then the following normalizations can be demostrated:
 * 	- considering g as one single hit, I (g * f) = 1
 * 	- considering g as a platou g = k, I(k * f) = k * N1
 * 	- consider G as the responce function of the detector, than g1 ϵ G and g2 ϵ G
 * 	means g2 = c * g1. I (g * f) depend only on the actual charge of g. If
 * 	we assume that g = a * f for f != 0, than I (g * f) = N2.
 * 	This is the optimum filter hypothesys.
 * Currently no normalization is done in the output, which is close to consider
 * short peaks (which is false). The optimum filter hypotesys is not respected
 * since the shape is slightly different from G. A normalization for extended
 * peaks depends on the shape of the data and should be calibrated.
 */
#include "bx_splitting_filter.hh"
#include "bx_echidna_event.hh"
#include "barn_interface.hh"
#include <math.h>
#include <algorithm>
#include <TH1F.h>
#include <TH2F.h>
#include <TSpectrum.h>

// ctor
bx_splitting_filter::bx_splitting_filter (): bx_base_module("bx_splitting_filter", bx_base_module::main_loop) {
  require_event_stage (bx_detector::laben, bx_base_event::clustered);
  peak_mode_names["tspectrum_simple"] = tspectrum_simple;
  peak_mode_names["tspectrum_highres"] = tspectrum_highres;
}

// module interface
void bx_splitting_filter::begin () {
  i4_window_width = get_parameter ("window_width").get_int ();
  i4_zero_bin = get_parameter ("zero_bin").get_int ();
  i4_sample_width = get_parameter ("sample_width").get_int ();
  i4_sample_ramp_lenght = get_parameter ("sample_ramp_lenght").get_int ();
  f4_sample_exp_tao = get_parameter ("sample_exp_tao").get_float ();
  i4_charge_threshold = get_parameter ("peak_charge_threshold").get_int ();
  i4_time_threshold = get_parameter ("peak_time_threshold").get_int ();
  const std::string &pm = get_parameter ("peak_mode").get_string ();
  if (!peak_mode_names.check (pm)) get_message (bx_message::critic) << "unknown peak mode " << pm << dispatch;
  else peak_mode = peak_mode_names[pm];
  
    // Histograms and TSpectrum stuff
  h_shape_out = new TH1F ("h_shape_out", "Histogram of shape_out", i4_window_width, -i4_zero_bin, i4_window_width - i4_zero_bin);
  peak_time_vs_energy = new TH2F ("peak_time_vs_energy", "Peak time vs peak energy", 300, -60, 540, 600, 0, 600);
  barn_interface::get ()->store (barn_interface::file, h_shape_out, this);
  barn_interface::get ()->store (barn_interface::file, peak_time_vs_energy, this);
  S = new TSpectrum (max_npeaks);
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,0,0)
  tmp_in = new double [i4_window_width];
  spectrum = new double [i4_window_width];
#else
  tmp_in = new float [i4_window_width];
  spectrum = new float [i4_window_width];
#endif

    // Internal buffer for FFT operations
  double *sample = (double *)fftw_malloc (sizeof(double) * i4_window_width);
  sample_F = (double *)fftw_malloc (sizeof(double) * i4_window_width);
  shape_in = (double *)fftw_malloc (sizeof(double) * i4_window_width);
  shape_in_F = (double *)fftw_malloc (sizeof(double) * i4_window_width);
  shape_out_F = (double *)fftw_malloc (sizeof(double) * i4_window_width);
  shape_out = (double *)fftw_malloc (sizeof(double) * i4_window_width * 2); // larger to allow offset   tmp_in = new float [i4_window_width];

     // FFTW plans
  fftw_plan sample_plan = fftw_plan_r2r_1d (i4_window_width, sample, sample_F, FFTW_R2HC, FFTW_ESTIMATE);
  plan_in = fftw_plan_r2r_1d (i4_window_width, shape_in, shape_in_F, FFTW_R2HC, FFTW_PATIENT);
  plan_out = fftw_plan_r2r_1d (i4_window_width, shape_out_F, shape_out + i4_sample_ramp_lenght, FFTW_HC2R, FFTW_PATIENT); // shift the output on left to have correct peak allignement

    // Fill the sample; time invert and translate the sample to the end of the histogram to use time aliasing
    // Calculate the integrals
  f8_sample_area = f8_sample_square_area = 0;
  for (int32_t i = 0; i < i4_window_width; i++) {
    if (i < i4_sample_ramp_lenght) sample[i4_window_width - 1 - i] = double(i + 1) / i4_sample_ramp_lenght; // Have maximun at 1
    else if (i < i4_sample_width) sample[i4_window_width - 1 - i] = ::expf ((i4_sample_ramp_lenght - i - 1) / f4_sample_exp_tao);
    else sample[i4_window_width - 1 - i] = 0;
    f8_sample_area += sample[i4_window_width - 1 - i];
    f8_sample_square_area += sample[i4_window_width - 1 - i] * sample[i4_window_width - 1 - i];
  }
    // Precalculate sample FFT
  fftw_execute(sample_plan);

    // Sample in time space and plan no more in use
  fftw_free (sample);
  fftw_destroy_plan(sample_plan);
}

bx_echidna_event* bx_splitting_filter::doit (bx_echidna_event *ev) {
    // Loop on every cluster
  for (int32_t i = 0; i < ev->get_laben ().get_nclusters (); i++) {
      // Get cluster reference
    bx_laben_cluster& cluster = ev->get_laben ().get_cluster (i);
    
      // 1) Filter the hit time distribution with the best filter obtaining SHAPE_OUT 

      // Fill the SHAPE_IN histogram
    std::fill_n (shape_in, i4_window_width, 0);
    for (int32_t j = 0; j < cluster.get_clustered_nhits (); j++) {
      const bx_laben_clustered_hit& hit = cluster.get_clustered_hit (j);
      int32_t bin = int32_t (::roundf (hit.get_time ()) + i4_zero_bin);
      if (bin >= i4_window_width) continue;
      shape_in[bin] ++;
    }

      // Do the FFT into shape_F
    fftw_execute (plan_in);
    
      // Multiply the 2 complex vector (see fftw reference for explanation of the math)
    shape_out_F[0] = sample_F[0] * shape_in_F[0];
    shape_out_F[i4_window_width / 2] *= sample_F[i4_window_width / 2];
    for (int32_t j = 1; j < i4_window_width / 2; j++) {
      double re = shape_in_F[j] * sample_F[j] - shape_in_F[i4_window_width - j] * sample_F[i4_window_width - j];
      double im = shape_in_F[j] * sample_F[i4_window_width - j] + shape_in_F[i4_window_width - j] * sample_F[j];
      shape_out_F[j] = re;
      shape_out_F[i4_window_width - j] = im;
    }

      // Do FFT back to the time space into the SHAPE_OUT histogram
    std::fill_n (shape_out, i4_sample_ramp_lenght, 0);
    fftw_execute (plan_out);

      // Normalize (since fftw does not) and strip ripple + undershoots
    for (int32_t i = 0; i < i4_window_width; i++) {
      shape_out[i] /= i4_window_width;
      if (shape_out[i] < 0.1) shape_out[i] = 0;
    }

      // 2) Search for peaks
    reset_histograms (ev->get_event_number ());
    int32_t npeaks = 0;
    float peaks[max_npeaks];
    if (peak_mode == tspectrum_simple) {
      h_shape_out->SetAxisRange (-i4_zero_bin, 500);
      S->Search (h_shape_out);
      npeaks = S->GetNPeaks ();
      std::copy (S->GetPositionX (), S->GetPositionX () + npeaks, peaks);
    } else if (peak_mode == tspectrum_highres) {
      std::copy (shape_out, shape_out + i4_window_width, tmp_in);
      npeaks = S->SearchHighRes (tmp_in, spectrum, i4_window_width, 8., 2, false, 3, false, 3);
      for (int32_t i = 0; i < npeaks; i++) {
        int32_t t = int32_t(S->GetPositionX ()[i]);
        if (t >= i4_window_width) t = i4_window_width - 1;
        peaks[i] = t - i4_zero_bin;
      }
    }
    std::sort (peaks, peaks + npeaks);

      // 3) Validate peaks
    int32_t end_prev_peak = 0;
    for (int32_t i = 0; i < npeaks; i++) {
      if (peaks[i] < i4_sample_width / -2) continue; // ignore fluctuation at the beginning of the shape
      int32_t peak_bin = int32_t(peaks[i]) + i4_zero_bin;
      int32_t start_peak = 0;
      switch (peak_mode) {
        case tspectrum_simple:
          start_peak = search_peak_base (shape_out, peak_bin, end_prev_peak);
          break;
        case tspectrum_highres:
          start_peak = search_peak_base (spectrum, peak_bin, end_prev_peak);
          break;
      };
      float charge = shape_out[peak_bin] - shape_out[start_peak];
      int32_t duration = peak_bin - start_peak;
      int32_t t = peak_bin - i4_zero_bin;
      if (duration < i4_time_threshold) ignore_peak ("short peak", i, t, charge, duration);
      else if (charge < i4_charge_threshold) ignore_peak ("under-threshold", i, t, charge, duration);
      else if (charge < ::sqrtf (shape_out[peak_bin])) ignore_peak ("poissoninan", i, t, charge, duration);
      else {
        //if (start_peak != end_prev_peak) get_message (bx_message::debug) << "found peak " << i << " at time " << t << "ns with charge " << charge << "hits and duration " << duration << "ns" << dispatch;
        //else get_message (bx_message::debug) << "found malformed peak " << i << " at time " << t << "ns with charge " << charge << "hits and duration " << duration << "ns" << dispatch;
        peak_time_vs_energy->Fill (t, charge);
	cluster.get_split_peaks ().push_back (bx_split_peak (t, charge));
        end_prev_peak = peak_bin;
      }
    }
  }
  return ev;
}

void bx_splitting_filter::end () {
  fftw_free (sample_F);
  fftw_free (shape_in);
  fftw_free (shape_in_F);
  fftw_free (shape_out);
  fftw_free (shape_out_F);
  fftw_destroy_plan (plan_in);
  fftw_destroy_plan (plan_out);
}

template<typename FLOAT> int32_t bx_splitting_filter::search_peak_base (FLOAT *v, int32_t& peak, int32_t low_limit) {
  if (peak - low_limit < 4) return low_limit;
  int32_t new_peak = peak;
  for (int32_t i = peak; i > low_limit + 4; i--) {
    float d = ::sqrtf (v[i]) * -0.03;
    if ((v[i - 1] - v[i]) < d || v[i] < 0.2) break;
    new_peak = i - 1;
  }
  //if (new_peak != peak) std::cout << "peak moved from " << peak << " to " << new_peak << std::endl;
  peak = new_peak;

  for (int32_t i = peak; i > low_limit + 4; i--) {
    bool negative_derivate = true;
    bool zero = true;
    for (int32_t j = 0; j < 4; j++) {
      float d = v[i - j] - v[i - j - 1];
      if (d > -0.2 || v[i] < 0.2) negative_derivate = false;
      if (v[i] > 0.2) zero = false;
    }
    //cout << "t " << i - i4_zero_bin << "ns " << negative_derivate; for (int32_t j = 0; j < 4; j++) cout << " " << v[i - j]; cout << endl;
    if (zero) return i;

    if (negative_derivate) return i;
  }
  return low_limit;
}

void bx_splitting_filter::ignore_peak (const std::string& msg, int32_t current_peak, float t, float charge, float duration) {
  bx_message::message_level level = bx_message::debug;
  if (current_peak == 0) level = bx_message::log;
  //get_message (level) << "ignoring " << msg << " peak " << current_peak << " at time " << t << "ns with charge " << charge << "hits and duration " << duration << "ns"<< dispatch;
}

void bx_splitting_filter::reset_histograms (int32_t evnum) {
  h_shape_out->Reset ();
  for (int32_t i = 0; i < i4_window_width; i++) {
    if (isnan (shape_out[i])) get_message (bx_message::critic) << "NaN found in shape_out for event " << evnum << " for bin " << i << dispatch;
    if (isnan (shape_out_F[i])) get_message (bx_message::critic) << "NaN found in shape_out_F for event " << evnum << " for bin " << i << dispatch;
    h_shape_out->Fill (i - i4_zero_bin, shape_out[i]);
  }
}
/*
 * $Log: bx_splitting_filter.cc,v $
 * Revision 1.26  2009/10/23 12:40:46  razeto
 * Muted
 *
 * Revision 1.25  2008-02-21 14:03:36  ludhova
 * debug NAN check
 *
 * Revision 1.24  2008-02-21 11:50:25  ludhova
 * NAN check also for h_shape_out_F
 *
 * Revision 1.23  2008-02-21 11:22:00  razeto
 * Check on NaN added on fft output
 *
 * Revision 1.22  2007-07-08 09:53:31  razeto
 * Depend only on clustered event but neutrino trigger
 *
 * Revision 1.21  2007-06-03 14:02:23  razeto
 * quiter
 *
 * Revision 1.20  2007-04-13 13:56:14  razeto
 * BE quiter
 *
 * Revision 1.19  2007-04-11 15:39:54  razeto
 * Do not use cout for messages
 *
 * Revision 1.18  2007-03-29 15:05:38  razeto
 * Fixed a typo
 *
 * Revision 1.17  2007-03-29 15:04:43  razeto
 * New comment on normalization
 *
 * Revision 1.16  2007-03-29 14:27:10  razeto
 * Better splitting algorithm
 *
 * Revision 1.15  2007-03-20 14:58:17  razeto
 * Some upgrades: comment, timing and thresholds
 *
 * Revision 1.14  2006/09/12 13:58:39  razeto
 * Do not return 0 unless when skipping event
 *
 * Revision 1.13  2006/08/21 11:19:02  razeto
 * Updated to new barn_interface
 *
 * Revision 1.12  2006/01/02 21:23:46  razeto
 * Changed test target to file target for db_barn
 *
 * Revision 1.11  2005/09/19 15:19:23  razeto
 * Removed a useless histogram
 *
 * Revision 1.10  2005/03/18 16:31:39  razeto
 * Removed a warn
 *
 * Revision 1.9  2005/03/15 09:41:42  razeto
 * Some upgrades in comments
 *
 * Revision 1.8  2004/12/21 17:18:14  razeto
 * Added event writing support
 *
 * Revision 1.7  2004/12/19 17:03:05  razeto
 * Added generation of another debugging histogram
 *
 * Revision 1.6  2004/12/17 15:40:34  razeto
 * Indroduced peak finding algorithm.
 * The identification of peaks work but need to be tested; some magic numbers
 * still hardcoded.
 * Waiting for the presence in the event of tagging fields.
 *
 * Revision 1.5  2004/12/13 15:49:53  razeto
 * Fixed few bugs and added some comments
 *
 * Revision 1.4  2004/12/13 12:39:17  razeto
 * Changed to use fftw version 2 instead of version 3 (not available on cluster)
 *
 * Revision 1.3  2004/12/10 13:53:47  razeto
 * A first working filter
 *
 * Revision 1.2  2004/12/09 17:45:04  razeto
 * Some fixing
 *
 * Revision 1.1  2004/12/09 17:10:06  razeto
 * First code (not working)
 *
 *
 */
