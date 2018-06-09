/* BOREXINO Reconstruction program
 *
 * Author: Daniela Manuzio <dmanuzio@lngs.infn.it>
 * Maintainer: Barbara Caccianiga <barbara.caccianiga@mi.infn.it>
 *
 * $Id: interpolator.hh,v 1.1 2004/12/23 10:28:58 razeto Exp $
 *
 * Implementation of an interpolation algorithm 
 *
 */

#ifndef _INTERPOLATOR_H
#define _INTERPOLATOR_H

#include "cmap.hh"
#include "messenger.hh" 
#include "bx_base_module.hh"
#include "vdt.hh"
#include "bx_named.hh"
#include <vector>

class interpolator: public bx_named {
  public:
    interpolator (int n, const double* x, const double* y);
    virtual ~interpolator () { delete [] f8_deriv; delete [] f8_x_coord; delete [] f8_y_coord; }

    double get_value (double val);
  private:
    const int i4_order;
    double* f8_deriv;
    double* f8_x_coord;
    double* f8_y_coord;
};

#endif
