
#ifndef _BX_POSITION_RECO_LNGS_HH
#define _BX_POSITION_RECO_LNGS_HH
#include "bx_base_module.hh"
#include "interpolator.hh"
#include "light_guide.hh"
#include "TVector3.h"


class bx_laben_cluster;
class TMinuit;

class bx_position_reco_lngs : public bx_base_module {
  public:
    bx_position_reco_lngs ();
    virtual ~bx_position_reco_lngs () {}
  
    // Operations for the framework 
    virtual void begin ();
    virtual bx_echidna_event* doit (bx_echidna_event *ev);
    virtual void end ();
    double my_fcn (double *x);

  private:
    double f8_tmax;
    float f4_ref_index;
    float f4_t_0;
    int f4_minuit_ok;
    int f4_minuit_blank;
    int f4_minuit_unreadable;
    int f4_minuit_unknown;
    int f4_minuit_abnormal;

    int f4_iconv_0;    
    int f4_iconv_1;    
    int f4_iconv_2;    
    int f4_iconv_3;    

    interpolator* pdf1;
    interpolator* pdf2;
    interpolator* pdf3;
    interpolator* pdf4;
    interpolator* pdf5;
    interpolator* pdf6;
    interpolator* pdf7;
    interpolator* pdf8;
    interpolator* pdf9;
    interpolator* pdf10;
    
    std::vector<TVector3> *f4_pmt_positions;
    std::vector<TVector3> *f4_pmt_normal_vectors;
    std::vector<bool> *b_pmt_cone;

    // to be used in my_fcn (see minuit root workaround)
    TMinuit *p_minuit;
 //   bx_echidna_event *p_fit_ev;
    bx_laben_cluster *p_fit_clu;
//    int i4_fit_cluster;
    int i4_n_iterations;
    int i4_hesse;
    
    bool b_iconv;

};

#endif
