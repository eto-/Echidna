/* BOREXINO Reconstruction program
 *
 * Author: Daniela Manuzio <dmanuzio@lngs.infn.it>
 * Maintainer: Barbara Caccianiga <barbara.caccianiga@mi.infn.it>
 *
 * $Id: interpolator.cc,v 1.3 2005/11/08 15:53:59 misiaszek Exp $
 *
 * Implementation of bx_interpolator
 *
 */
#include "interpolator.hh"
#include "messenger.hh"
#include <iostream> 
#include <math.h>

// ctor
interpolator::interpolator (int n, const double* x, const double* y): bx_named("interpolator"),  i4_order(n) {
  if (i4_order < 2) {
    get_message(bx_message::error) << "no interpolation done (< 2 points)" << dispatch;
    return;
  }

  double* par = new double[i4_order];
  f8_deriv = new double[i4_order];

  f8_x_coord = new double[i4_order];
  f8_y_coord = new double[i4_order];
  for (int i = 0; i < i4_order; i++) {
    f8_x_coord[i] = x[i];
    f8_y_coord[i] = y[i];
  }


  f8_deriv[0] = 0.;
  par[0] = 0.;

  for (int i = 1; i < i4_order - 1; i++) {
    double sig = (f8_x_coord[i] - f8_x_coord[i - 1])/(f8_x_coord[i + 1] - f8_x_coord[i - 1]);
    double p = sig * f8_deriv[i - 1] + 2.;
    f8_deriv[i] = (sig - 1.) / p;
    par[i] = (f8_y_coord[i + 1] - f8_y_coord[i]) / (f8_x_coord[i + 1] - f8_x_coord[i]) - (f8_y_coord[i] - f8_y_coord[i - 1]) / (f8_x_coord[i] - f8_x_coord[i - 1]);
    par[i] = (6 * par[i] / ( f8_x_coord[i + 1] - f8_x_coord[i - 1]) - sig * par[i - 1]) / p;
  }
  
  f8_deriv[i4_order - 1] = 0.;
      
  for (int i = i4_order - 2; i >= 0; i--) {
    f8_deriv[i] = f8_deriv[i] * f8_deriv[i + 1] + par[i];
  }
  
  delete [] par;
}

// return the interpolated value
double interpolator::get_value (double val) {
  int low = 0;
  int up = i4_order - 1;

  if (val < f8_x_coord[low]) return f8_y_coord[low];
  if (val > f8_x_coord[up]) return f8_y_coord[up];
    
  while (up-low > 1) {
    int index = (low+up)/2;
    if ( f8_x_coord[index] > val ) up = index;
    else low = index;
  }

  double h = f8_x_coord[up] - f8_x_coord[low];
  if (h == 0.){
    get_message(bx_message::error) << "no interpolation done (bad input)" << dispatch;
    return 0.;
  }
  
  double a = (f8_x_coord[up] - val) / h;
  double b = (val - f8_x_coord[low]) / h;
//
// for now, use the linear interpolation since cubic spline is not fully understood yet
//
  double interpolation = a * f8_y_coord[low] + b * f8_y_coord[up];
//
//  double interpolation = a * f8_y_coord[low] + b * f8_y_coord[up]+ 
//                         ((a * a * a - a) * f8_deriv[low] + (b * b * b - b) * f8_deriv[up]) * h * h / 6.;
  
  return interpolation;
}
