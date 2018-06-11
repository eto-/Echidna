/* BOREXINO Reconstruction program
 *
 * Author: Alessandro Razeto <Alessandro.Razeto@ge.infn.it>
 * Maintainer: Alessandro Razeto <Alessandro.Razeto@ge.infn.it>
 *
 * $Id: bx_position_reco_dbn.hh,v 1.8 2005/10/30 11:51:23 razeto Exp $
 *
 */
#ifndef _BX_POSITION_RECO_DBN_H
#define _BX_POSITION_RECO_DBN_H

#include "bx_base_module.hh"
#include "bx_rec_general.hh"
#include "TH2F.h"
#include "flight_path.hh"
#include <string>


class TMinuit;
class bx_position_reco_dbn: public bx_base_module {
  public:
    bx_position_reco_dbn ();
    virtual ~bx_position_reco_dbn () {}

    virtual void begin ();
    virtual bx_echidna_event* doit (bx_echidna_event *ev);
    virtual void end ();
    double my_fcn (int32_t npar, double *x, double *grad, int32_t iflag);

    static const uint16_t n_shell = 5;
    static const float shell_external_radius[n_shell];
  private:
    double f8_first_hits_gate;
    float f4_gauss_sigma[n_shell];
    float f4_t0_shift[n_shell];
    float f4_current_gauss_sigma_square;
    bool b_use_grad;
    
    float f4_sphere_radius;
    TMinuit *p_minuit;
    TH2F *p_iconv;
    flight_path path;

    int32_t i4_minuit_error_count;

      // to be used in my_fcn (see minuit root workaround)
    bx_echidna_event *p_fit_ev;
    int32_t i4_fit_cluster;
    int32_t i4_n_iterations;
};

#endif
/*
 * $Log: bx_position_reco_dbn.hh,v $
 * Revision 1.8  2005/10/30 11:51:23  razeto
 * Upgrade to have a better starting t
 *
 * Revision 1.7  2005/10/13 13:42:07  razeto
 * Added shell dependent pdf
 *
 * Revision 1.6  2005/10/12 14:44:26  razeto
 * Added derivative calculation, still in testing
 *
 * Revision 1.5  2005/10/06 21:30:16  razeto
 * Working toward analitic derivative
 *
 * Revision 1.4  2005/10/04 20:24:13  razeto
 * Added iterations value during minimization
 *
 * Revision 1.3  2005/09/22 11:54:21  razeto
 * Now the pdf can be slightly modified with 2 parameters
 *
 * Revision 1.2  2005/07/07 15:14:09  razeto
 * Reduced the minimization steps and added a debugging histogram
 *
 * Revision 1.1  2005/06/21 12:06:56  razeto
 * Added to repository
 *
 *
 */
