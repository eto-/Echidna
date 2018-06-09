/* BOREXINO Reconstruction program
 *
 * Author: Davide D'Angelo <davide.dangelo@lngs.infn.it>
 * Maintainer: Davide D'Angelo <davide.dangelo@lngs.infn.it>
 *
 * $Id: bx_track.cc,v 1.9 2014/12/11 21:27:12 wurm Exp $
 *
 * Implementation of bx_track.hh
 *
 */

#include "bx_track.hh"

bx_track_by_points::bx_track_by_points () : f4_x1  (0.), f4_y1  (0.), f4_z1  (0.), f4_x2  (0.), f4_y2  (0.), f4_z2  (0.),
	f4_dx1 (0.), f4_dy1 (0.), f4_dz1 (0.), f4_dx2 (0.), f4_dy2 (0.), f4_dz2 (0.), f4_theta(0.), f4_phi(0.), 
	f4_dtheta(0.), f4_dphi(0.), f4_impact(0.), f4_dimpact(0.), f4_labennormhits(-1), i4_error(0), b_downward(true) { 
}

bx_track_fitted::bx_track_fitted () : f8_alpha (0.), f8_beta (0.), f8_gamma (0.), f8_delta (0.),
                                      f8_alpha_error (0.), f8_beta_error (0.), f8_gamma_error (0.), f8_delta_error (0.),
				      f4_chi2 (0.),
				      f4_phi(0.), f4_theta(0.), f4_impact(0.), f4_dphi(0.), f4_dtheta(0.), f4_dimpact(0.),
				      f4_labennormhits(-1), b_downward(true), u1_points(0), b_is_valid (false) {
}
/*
 * $Log: bx_track.cc,v $
 * Revision 1.9  2014/12/11 21:27:12  wurm
 * added getters and variables to BxEvent that concern C-11 likelihood tagging
 * * for tracks: new distance-to- and projection-on-track getters, dE/dx calculation
 * * counting of neutron-like clusters in tt128, new flag of neutron clusters
 *
 * Revision 1.8  2010/08/06 23:00:15  ddangelo
 * fixed initialization order
 *
 * Revision 1.7  2010-08-03 15:58:05  wurm
 * introduced theta, phi, impact variables for fitted tracks
 *
 * Revision 1.6  2010-05-21 15:55:14  ddangelo
 * different things on muon tracks
 *
 * 1.a) old laben_track renamed as laben_track_energy
 * new laben_track_tof added
 *
 * 1.b) (global) track renemed as track_global at base event level
 * track_cmt added at base event level (track by points)
 *
 * 1) all getters updated/integrated
 * is_tracked variable updated/integrated accordingly. inizialization.
 * job ported to root event as well. copy done.
 * friendship with old/new module updated
 *
 * 2) bxtrack_by_points class:
 * - theta, phi and impact added as variables.
 * - errors added on all of the above.
 * - error code variable requested by cmt tracker added
 *
 * Revision 1.5  2009-07-31 15:39:50  ddangelo
 * debugging the work of the [previous commit
 *
 * Revision 1.4  2009-04-15 17:12:38  ddangelo
 * n_point changed into a bitfield, internal and root event (BxTrackFitted class)
 *
 * Revision 1.3  2008-07-17 15:55:00  ddangelo
 * added track direction, removed unused getter
 *
 * Revision 1.2  2008-05-09 15:46:23  ddangelo
 * muon track implemented as "by_points" and "fitted" classes inheriting from common interface (pure virtual class bx_track)
 * laben and muon host an istance "by_points"
 * echidna event hosts an istance "fitted"
 * structure duplicated at root level and filled
 * GetDistance duplicated in BxPosition class
 *
 * Revision 1.1  2008-02-21 15:58:31  ddangelo
 * added muon track as a separate class.
 * Currently used twice: laben and muon (later also globally fit).
 * both inner and roor event
 * added tracked stage in base event enum
 * bx_muon_rec_event renamed as bx_muon_tracked_event
 *
 *
 */
