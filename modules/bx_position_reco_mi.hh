/* BOREXINO Reconstruction program
 *
 * Author: Daniela Manuzio <dmanuzio@ge.infn.it>
 * Maintainer: Barbara Caccianiga <barbara.caccianiga@mi.infn.it>
 *
 * $Id: bx_position_reco_mi.hh,v 1.16 2011/03/02 17:21:13 razeto Exp $
 *
 * Event's position minimization 
*/

#ifndef _BX_POSITION_RECO_MI_HH
#define _BX_POSITION_RECO_MI_HH
#include "bx_base_module.hh"
#include "interpolator.hh"
#include "light_guide.hh"
#include "TVector3.h"
#include "TH1F.h"

class bx_echidna_event;
class TMinuit;

class bx_position_reco_mi : public bx_base_module {
  public:
    bx_position_reco_mi ();
    virtual ~bx_position_reco_mi () {}
  
    // Operations for the framework 
    virtual void begin ();
    virtual bx_echidna_event* doit (bx_echidna_event *ev);
    virtual void end ();
    double my_fcn (double *x);

  private:
    int32_t i4_keep_ratio;

    double f8_tmax;
    float f4_ref_index;
    float f4_t_0;
    int32_t f4_minuit_ok;
    int32_t f4_minuit_blank;
    int32_t f4_minuit_unreadable;
    int32_t f4_minuit_unknown;
    int32_t f4_minuit_abnormal;

    int32_t f4_iconv_0;    
    int32_t f4_iconv_1;    
    int32_t f4_iconv_2;    
    int32_t f4_iconv_3;    

    interpolator* pdf;
    
    std::vector<TVector3> *f4_pmt_positions;
    std::vector<TVector3> *f4_pmt_normal_vectors;
    std::vector<bool> *b_pmt_cone;

    // to be used in my_fcn (see minuit root workaround)
    TMinuit *p_minuit;
    bx_echidna_event *p_fit_ev;
    int32_t i4_fit_cluster;
    int32_t i4_n_iterations;
    int32_t i4_hesse;
    
    bool b_iconv;
    
    TH1F *p_iconv_ok;
    TH1F *p_iconv_nook;
    TH1F *p_edm_ok;
    TH1F *p_edm_nook;
    TH1F *p_histo_x;
    TH1F *p_histo_y;
    TH1F *p_histo_z;
    TH1F *p_histo_t;
    TH1F *p_histo_x_no;
    TH1F *p_histo_y_no;
    TH1F *p_histo_z_no;
    TH1F *p_histo_t_no;
};

#endif
