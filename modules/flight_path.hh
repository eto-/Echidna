/* BOREXINO Reconstruction program
 *
 * Author: Alessandro Razeto <Alessandro.Razeto@ge.infn.it>
 * Maintainer: Alessandro Razeto <Alessandro.Razeto@ge.infn.it>
 *
 * $Id: flight_path.hh,v 1.14 2007/06/13 11:08:56 razeto Exp $
 *
 * Given a source point and a clustered hit flight_path calculates
 * some parameters of the path the photon followed to reach the 
 * pmt.
 * 
 */
#ifndef _FLIGHT_PATH
#define _FLIGHT_PATH

#include "bx_named.hh"
#include "bx_rec_general.hh"

class db_channel_laben;
class bx_base_position;
class TH2F;
class flight_path: public bx_named {
  public:
    flight_path ();
    void init (double v[3], const db_channel_laben* ch_info); 	// v[3] = { x, y, z }
    void init (double x, double y, double z, const db_channel_laben* ch_info);
    void init (const bx_base_position& position, const db_channel_laben* ch_info);
    void set_refraction_index (float ref_index);
    float get_refraction_index () const { return f4_ref_index; }
    float get_c_medium () const { return f4_c_medium_m_ns; } 
    
    float get_time () { distance (); return f4_distance / f4_c_medium_m_ns; }
    float get_time_no_light_guide () { straight_distance (); return f4_straight_distance / f4_c_medium_m_ns; }
    float get_distance () { distance (); return f4_distance; }
    float get_distance_no_light_guide () { straight_distance (); return f4_straight_distance; }
    float get_theta () { theta (); return f4_theta; }
    void get_distance_derivative (float v[3]); 			// v[3] = { x, y, z }
    void get_time_derivative (float v[3]) { get_distance_derivative (v); v[0] /= f4_c_medium_m_ns; v[1] /= f4_c_medium_m_ns; v[2] /= f4_c_medium_m_ns; }
  private:
    float f4_v[3];
    float f4_pmt[3];
    bool b_has_cone;
    bool inited;
    float f4_ref_index;
    float f4_c_medium_m_ns;
    static int32_t create_histo;
    static TH2F *path;

      // some cached values
    float f4_distance, f4_straight_distance, f4_theta;
      // and their calculators
    void straight_distance ();
    void theta ();
    void distance ();

    void _init ();
};

#endif
/*
 * $Log: flight_path.hh,v $
 * Revision 1.14  2007/06/13 11:08:56  razeto
 * Do not rely on the sphere collector radius from db, instead calculate it pmt by pmt
 *
 * Revision 1.13  2007-05-20 12:48:11  razeto
 * Fixed initialization bug (fixes refidx problem)
 *
 * Revision 1.12  2006-05-08 17:31:31  razeto
 * Added db_channel patch for fadc (sent to the mailing list)
 *
 * Revision 1.11  2005/10/30 11:17:50  razeto
 * Now the refraction index can be overridden from user(module) but not directly from its own parameters
 *
 * Revision 1.10  2005/10/17 10:15:08  razeto
 * Added a missimg member
 *
 * Revision 1.9  2005/10/11 14:45:03  razeto
 * Added time derivative
 *
 * Revision 1.8  2005/10/11 14:21:22  razeto
 * Upgraded to cache step results and to calculate derivative
 *
 * Revision 1.7  2005/10/06 21:27:56  razeto
 * Added a small optimization (cached pmt position)
 *
 * Revision 1.6  2005/09/06 11:58:45  razeto
 * Added refraction index parameter
 *
 * Revision 1.5  2005/06/20 14:09:39  razeto
 * Updated
 *
 * Revision 1.4  2005/05/16 12:40:00  razeto
 * Added an initializer and a simple method for getting c
 *
 * Revision 1.3  2005/03/20 09:02:55  razeto
 * Upgraded to cache some variables and to fill on demand an histogram
 *
 * Revision 1.2  2005/03/18 16:30:09  razeto
 * Upgraded to working version
 *
 * Revision 1.1  2005/03/18 11:22:49  razeto
 * Added initial version
 *
 */
