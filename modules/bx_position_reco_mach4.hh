#ifndef _BX_POSITION_RECO_MACH4_HH
#define _BX_POSITION_RECO_MACH4_HH

#include "bx_base_module.hh"
#include "bx_scint_pdf_mach4.hh"
#include "bx_mach4_minuit_wrapper.hh"
#include "constants.hh"

class bx_echidna_event;

class bx_laben_cluster;

class bx_position_reco_mach4 : public bx_base_module {

  ScintPDF scint_pdf[3];
  MinuitWrapper *mn;

  double xpmt[2241];
  double ypmt[2241];
  double zpmt[2241];

  size_t iteration;
  bx_laben_cluster *cluster;

  mutable size_t first_hit, end_hit;
  
  float refidx; //refractive index for fixed-index reconstruction to be read from config file
  
  bool Reconstruct(double n_idx, bx_echidna_event * ev,
		   size_t cluster_idx,
		   double *pos_info);


public:
  bx_position_reco_mach4();
  virtual ~bx_position_reco_mach4();

  virtual void begin();
  virtual bx_echidna_event* doit(bx_echidna_event *ev);
  virtual void end();

  double Likelihood(double x, double y, double z, double t) const;

  ScintPDF *get_pdf(int32_t i = 0) {return &scint_pdf[i];}

};


#endif
