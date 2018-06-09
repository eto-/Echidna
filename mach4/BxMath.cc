#include "mach4/BxMath.h"
#include "mach4/Parameters.h"
#include <cmath>

/**
 * @file BxMath.cc
 * @brief Common mathematical functions.
 */

/**
 * @brief Estimate the number of photoelectrons, given the number
 * @brief of observed PMT hits and the number of live PMTs in the set.
 */
double BxMath::npe_from_nhits(size_t nhits, size_t npmts)
{
  if (nhits == 0 || npmts == 0)
    return 0.0;

  if (nhits >= npmts)
    // complete saturation; no idea what should be returned, but it should
    // be large
    nhits = npmts - 1;
  
  return std::log(1. - (1. * nhits) / npmts) / std::log(1. - 1. / npmts);
}

void BxMath::DecodeFlag(int flag, bool field[8])
{
  
  int remaining=flag;
  for(int i=0;i<8;i++) {
    field[i]=remaining%2;
    remaining=remaining/2;
  }
 
}

#ifdef DEBUG_RECON
# include <cassert>
# include <iostream>
  using std::cerr; using std::endl;
# define Assert(args) do { if (! (args)) { \
  cerr << "*** x1 = " << x1 << ", y1 = " << y1 << ", z1 = " << z1 << endl; \
  cerr << "*** x2 = " << x2 << ", y2 = " << y2 << ", z2 = " << z2 << endl; \
  cerr << "*** dist_to_pmt   = " << dist_to_pmt << endl; \
  cerr << "*** dist_in_scint = " << dist_in_scintillator << endl; \
  cerr << "*** dist_in_buff  = " << dist_in_buffer << endl; \
  cerr << "*** dplus         = " << dplus << endl; \
  cerr << "*** diff          = " << diff << endl; \
  cerr.flush(); assert(args); } } while (0)
#else
# define Assert(args) /* nothing */
#endif

static const double r_IV2 = VESSEL_RADIUS * VESSEL_RADIUS;

/**
 * @brief Determine the distance in scintillator and in buffer between two points in space.
 * @param x1 X position of first point, in meters.
 * @param y1 Y position of first point.
 * @param z1 Z position of first point.
 * @param x2 X position of second point. (Must be outside the Inner Vessel!!!)
 * @param y2 Y position of second point. (Must be outside the Inner Vessel!!!)
 * @param z2 Z position of second point. (Must be outside the Inner Vessel!!!)
 * @param dist_to_pmt Distance in meters between two points. (return value)
 * @param dist_in_scintillator Distance traveled in scintillator. (return value)
 * @param dist_in_buffer Distance traveled in buffer. (return value)
 */
void BxMath::get_distances(double x1, double y1, double z1,
			   double x2, double y2, double z2,
			   double & dist_to_pmt, double & dist_in_scintillator,
			   double & dist_in_buffer)
{
  // x1, y1, z1 are test point coordinates
  // x2, y2, z2 are PMT coordinates
  // 
  // WARNING WARNING WARNING:
  // We assume for simplicity that (x2, y2, z2) is always outside the IV.
  double dplus = 0, diff = 0;
  Assert(distance(x2, y2, z2, 0, 0, VESSEL_Z_CENTER) > VESSEL_RADIUS);

  // Initialize return variables assuming that the path goes only through
  // the buffer (later we will test this assumption and fix if needed).
  dist_to_pmt = distance(x1, y1, z1, x2, y2, z2);
  dist_in_buffer = dist_to_pmt;
  dist_in_scintillator = 0;

  // Transform to test point coordinates in order to use the formula
  // at http://en.wikipedia.org/wiki/Ray-sphere_intersection
  const double lx = (x2 - x1) / dist_to_pmt,
               ly = (y2 - y1) / dist_to_pmt,
               lz = (z2 - z1) / dist_to_pmt;

  const double sx = -x1, sy = -y1, sz = VESSEL_Z_CENTER - z1;
  const double innerprod = lx*sx + ly*sy + lz*sz;
  const double discriminant = innerprod * innerprod
                              - (sx*sx + sy*sy + sz*sz - r_IV2);
  // N.B. by construction, lx*lx + ly*ly + lz*lz must equal 1.

  if (discriminant <= 0)
    // no intersection or only tangential intersection of path with IV;
    // path from test point to PMT must only see buffer
    return;

  const double sqrt_disc = sqrt(discriminant);
  dplus = innerprod + sqrt_disc; // larger solution to the eqn
  if (dplus <= 0)
    // line through test point and PMT intersects IV, but only in the
    // direction opposite the test point, not in the path between, which
    // must see only buffer
    return;

  if (dplus >= dist_to_pmt)
    // test point is on far side of PMT from IV!  Minuit, what are you
    // doing way out there??
    return;

  // Check whether the smaller solution of the quadratic
  // (dminus := dplus - diff) is negative (dplus < diff).  If yes,
  // the path from test point to PMT crosses the IV only once (and
  // the test point is inside the IV).  If no, the path crosses the
  // IV twice (and the test point is outside the IV and on the other side
  // of the IV from the PMT).
  diff = 2 * sqrt_disc;
  dist_in_scintillator = (dplus <= diff) ? dplus : diff;
  dist_in_buffer = dist_to_pmt - dist_in_scintillator;

  // Distance in scintillator must be less than total distance
  Assert(dist_in_scintillator < dist_to_pmt);
  // Distance in scintillator must be at most the vessel diameter
  // (plus epsilon to account for floating point imprecision)
  Assert(dist_in_scintillator <= 2 * VESSEL_RADIUS + 0.01);
  // If test point is inside IV, then distance in buffer must be at
  // least the minimum distance between IV surface and PMT sphere
  // (minus epsilon to account for floating point imprecision)
  Assert(distance(x1, y1, z1, 0, 0, VESSEL_Z_CENTER) >= VESSEL_RADIUS ||
         dist_in_buffer >= 6.55 - VESSEL_RADIUS - VESSEL_Z_CENTER);
}

