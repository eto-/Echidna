/* BOREXINO Reconstruction program
 *
 * Author: Daniela Manuzio <dmanuzio@ge.infn.it>
 * Maintainer: Barbara Caccianiga <barbara.caccianiga@mi.infn.it>
 *
 * $Id: light_guide.cc,v 1.2 2007/07/24 20:12:19 bcaccian Exp $
 *
 * Implemenentation of the Borexino light guides. 
 * It computes the mean path lenght travelled by each photon
 * inside the LG ( by interpolation from the values tabled in 
 * "Light Collectors for Borexino" note, by I.Manno)
 * 
*/
#include "light_guide.hh"

light_guide *light_guide::me = 0;

light_guide* light_guide::get () {
  if (!me) me = new light_guide();                                                                                          
  return me;
}

light_guide::light_guide (): bx_named("light_guide") {

  // Mean path lenght of the photon in the LG (distance) as a function of the incident angle (theta).
  // Data points taken from the relative table in the note "Light Collectors for Borexino", by I.Manno
  int points = 18;
  double theta[18] = {0.,5.,10.,15.,20.,25.,30.,35.,36.,37.,38.,39.,40.,41.,42.,43.,44.,47.};
  double path[18] = {0.2648, 0.2652, 0.2685, 0.2723, 0.2796, 0.2901, 0.3069, 0.3307, 0.3375, 0.3434,
    		     0.3490, 0.3510, 0.3561, 0.3622, 0.3691, 0.3770, 0.3853, 0.4177};		
		     
  lg_int = new interpolator(points, theta, path); 

  get_message(bx_message::debug) << "Data interpolation done" << dispatch;
  return;
}

double light_guide::get_path (double angle) {
  return lg_int->get_value(angle);
}
