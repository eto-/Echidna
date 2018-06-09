/* BOREXINO Reconstruction program
 *
 * Author: Stefano Davini <stefano.davini@ge.infn.it>
 * Maintainer: Davide D'Angelo <davide.dangelo@lngs.infn.it>
 *
 * $Id: bx_neutron_event.cc,v 1.5 2008/08/26 18:05:38 ddangelo Exp $
 *
 * Implementation of bx_neutron_event
 *
 */

#include "bx_neutron_event.hh"

bx_neutron_pulse::bx_neutron_pulse () : f4_charge(0.), f4_amplitude(0.), f4_peak_time(0.), f4_rise_time(0.), f4_fall_time(0.), f4_x(0.), f4_y(0.), f4_z(0.), f4_dx(0.), f4_dy(0.), f4_dz(0.) {
}

bx_neutron_event::bx_neutron_event () : b_is_enabled(false), b_is_associated(false), i4_n_neutrons(0) {
  pulses.clear();
}

/*
 * $Log: bx_neutron_event.cc,v $
 * Revision 1.5  2008/08/26 18:05:38  ddangelo
 * missing inizialization
 *
 * Revision 1.4  2008-08-26 15:30:30  ddangelo
 * added neutron enabled and association flag, muon aligned flags
 *
 * Revision 1.3  2008-08-05 16:30:52  ddangelo
 * added fall time, some variables renamed
 *
 * Revision 1.2  2008-07-17 16:30:30  ddangelo
 * debugging
 *
 * Revision 1.1  2008-07-11 17:05:45  ddangelo
 * added
 *
 *
 */
