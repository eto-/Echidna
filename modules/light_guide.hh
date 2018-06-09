/* BOREXINO Reconstruction program
 *
 * Author: Daniela Manuzio <dmanuzio@ge.infn.it>
 * Maintainer: Barbara Caccianiga <barbara.caccianiga@mi.infn.it>
 *
 * $Id: light_guide.hh,v 1.1 2004/12/23 10:28:58 razeto Exp $
 *
 * Implemenentation of the Borexino light guides. 
 * It computes the mean path lenght travelled by each photon
 * inside the LG ( by interpolation from the values tabled in 
 * "Light Collectors for Borexino" note, by I.Manno)
 * 
*/

#ifndef _LIGHT_GUIDE_HH
#define _LIGHT_GUIDE_HH
#include "interpolator.hh"
#include "bx_base_module.hh"
#include "bx_named.hh"
#include <vector>

class light_guide: public bx_named {
  public:
    static light_guide* get ();
    
    double get_path (double angle);
  private:
    light_guide ();
    ~light_guide () { delete lg_int; }
    static light_guide *me;

    interpolator* lg_int;
};


#endif
