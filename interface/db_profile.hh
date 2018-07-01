/* BOREXINO Reconstruction program
 *
 * Author: Alessandro Razeto <Alessandro.Razeto@ge.infn.it>
 * Maintainer: Alessandro Razeto <Alessandro.Razeto@ge.infn.it>
 *
 * $Id: db_profile.hh,v 1.30 2015/01/09 15:03:08 misiaszek Exp $
 *
 * The database interface for profile informations
 *
 */
#ifndef _BD_PROFILE_H
#define _BD_PROFILE_H
#ifndef CYCLE_NUMBER
#define CYCLE_NUMBER 18
#endif
#include "db_acl.hh"

#include <TObject.h>

#include <string>
#include <map>
#include <vector>
//#include "constants.hh"

class db_profile: public db_acl, public TObject {
  public:
    enum channel_description_type {
      // general description
      empty = 0,
      ordinary = 1,
      trigger = 2,
      laser = 3,
      // laben description
      pulser = 4,
      cngs_id = 5,
      cngs_od = 6,
      cngs_trg = 7,
      // nuon descriptions
      integrity = 20,
      clock = 21,
      omt = 22,
    };

    db_profile() : db_acl(), TObject() {};
    
    class pmt_coordinates: public TObject { 
      public:
        float x, y, z, theta, phi; 
	pmt_coordinates() : TObject () {};
	ClassDef(db_profile::pmt_coordinates,CYCLE_NUMBER);
    };
    
    channel_description_type logical_channel_description (int32_t lg) const { return map_get(lg,logical_channel_description_v,"logical_channel_description_v"); }

    float laben_crate_delay		(int32_t crate) const { return map_get(crate,laben_crate_delay_v,"laben_crate_delay_v"); }
    
    int32_t neutrino_trigger_tag 	() const { return i4_neutrino_trigger_tag; }
    int32_t muon_trigger_tag 		() const { return i4_muon_trigger_tag; }
    int32_t neutron_trigger_tag	() const { return i4_neutron_trigger_tag; }
    int32_t laser266_trigger_tag 	() const { return i4_laser266_trigger_tag; }
    int32_t laser355_trigger_tag 	() const { return i4_laser355_trigger_tag; }
    int32_t laser394_trigger_tag 	() const { return i4_laser394_trigger_tag; }
    int32_t pulser_trigger_tag 	() const { return i4_pulser_trigger_tag; }
    int32_t random_trigger_tag 	() const { return i4_random_trigger_tag; }
    int32_t muon_tdc_range     	() const { return i4_muon_tdc_range; }
//    float  muon_tdc_range_ns     	() const { return i4_muon_tdc_range*constants::muon::tdc::ns_per_clock; }
    float trigger_start_gate		() const { return f_trigger_start_gate; }
    float trigger_end_gate		() const { return f_trigger_end_gate; }
    float laser394_trigger_start_gate	() const { return f_laser394_trigger_start_gate; }
    float laser394_trigger_end_gate	() const { return f_laser394_trigger_end_gate; }
    
    int32_t number_of_internal_pmt 		() const { return i_number_of_internal_pmt; }
    int32_t number_of_external_pmt 		() const { return i_number_of_external_pmt; }

    bool pmt_has_cone 			(int32_t lg) const { return map_get(lg,pmt_has_cone_v,"pmt_has_cone_v"); }
    int32_t pmt_hole_id 			(int32_t lg) const { return map_get(lg,pmt_hole_id_v,"pmt_hole_id_v"); }
    const pmt_coordinates& get_pmt_coordinates (int32_t lg) const { return map_get(lg,pmt_coordinates_v,"pmt_coordinates_v"); }
    int32_t pmt_fiber_bundle		(int32_t lg) const { return map_get(lg,pmt_fiber_bundle_v,"pmt_fiber_bundle_v"); }
    
    float stainless_steel_sphere_radius	() const { return f_stainless_steel_sphere_radius; }
    float inner_vessel_radius 		() const { return f_inner_vessel_radius; }
    float light_collector_radius 	() const { return f_light_collector_radius; }
    float light_guide_exit_aperture 	() const { return f_light_guide_exit_aperture; }
    float light_guide_entry_aperture	() const { return f_light_guide_entry_aperture; }
    float light_guide_lenght		() const { return f_light_guide_lenght; }
    float pmt_chathode_radius		() const { return f_pmt_chathode_radius; }
    float pmt_bulb_radius		() const { return f_pmt_bulb_radius; }

    float pmt_time_jitter		() const { return f_pmt_time_jitter; }
    float pmt_reflectivity		() const { return f_pmt_reflectivity; }
    float light_guide_reflectivity	() const { return f_light_guide_reflectivity; }
    float sss_specular_reflectivity	() const { return f_sss_specular_reflectivity; }
    float sss_diffusion_reflectivity	() const { return f_sss_diffusion_reflectivity; }
    float sss_uni_reflectivity		() const { return f_sss_uni_reflectivity; }

  private:
    std::map<std::string, db_profile::channel_description_type> channel_description_type_map;
    std::map<int32_t, channel_description_type> logical_channel_description_v;

    
    std::map<int32_t, float> laben_crate_delay_v;
    
    int32_t i4_neutrino_trigger_tag;
    int32_t i4_muon_trigger_tag;
    int32_t i4_neutron_trigger_tag;
    int32_t i4_laser266_trigger_tag;
    int32_t i4_laser355_trigger_tag;
    int32_t i4_laser394_trigger_tag;
    int32_t i4_pulser_trigger_tag;
    int32_t i4_random_trigger_tag;
    int32_t i4_muon_tdc_range;
    float f_trigger_start_gate;
    float f_trigger_end_gate;
    float f_laser394_trigger_start_gate;
    float f_laser394_trigger_end_gate;

    int32_t i_number_of_internal_pmt;
    int32_t i_number_of_external_pmt;
    
    float f_stainless_steel_sphere_radius;
    float f_inner_vessel_radius;
    float f_light_collector_radius;
    float f_light_guide_exit_aperture;
    float f_light_guide_entry_aperture;
    float f_light_guide_lenght;
    float f_pmt_chathode_radius;
    float f_pmt_bulb_radius;

    float f_pmt_time_jitter;
    float f_pmt_reflectivity;
    float f_light_guide_reflectivity;
    float f_sss_specular_reflectivity;
    float f_sss_diffusion_reflectivity;
    float f_sss_uni_reflectivity;
    
    
    std::map<int32_t, bool> pmt_has_cone_v;
    std::map<int32_t, int32_t> pmt_hole_id_v;
    std::map<int32_t, pmt_coordinates> pmt_coordinates_v;
    std::map<int32_t, int32_t> pmt_fiber_bundle_v;
    
    
    void m_load_laben_channel_mapping (int32_t profile_id);
    void m_load_muon_channel_mapping (int32_t profile_id);
    void m_load_muon_general_setting (int32_t profile_id);
    void m_load_detector_geometry (int32_t profile_id);
    void m_load_profile_table (int32_t profile_id);
    void m_load_physics_constants (int32_t profile_id);

    db_profile (int32_t profile_id);
    friend class bx_dbi;
    ClassDef(db_profile,CYCLE_NUMBER)
};
#endif
/*  
 *  $Log: db_profile.hh,v $
 *  Revision 1.30  2015/01/09 15:03:08  misiaszek
 *  cycle_18 new unstable
 *
 *  Revision 1.29  2013/06/18 18:56:37  razeto
 *  cycle_17 new unstable
 *
 *  Revision 1.28  2013-02-02 09:01:50  razeto
 *  Incremented to cycle_16 (cycle 15 was lost)
 *
 *  Revision 1.27  2012-03-22 19:08:36  razeto
 *  Added cngs reference channels
 *
 *  Revision 1.26  2011-04-19 05:54:58  razeto
 *  Moved to cycle 15 unstable
 *
 *  Revision 1.25  2011-02-18 17:10:05  ddangelo
 *  major code cleanup: removed fadc throughout the program
 *
 *  Revision 1.24  2010-08-06 17:20:16  razeto
 *  Moving to cycle 14
 *
 *  Revision 1.23  2009-12-03 16:11:14  misiaszek
 *  Change for new streamer and rootcint in channel_description_type_map
 *
 *  Revision 1.22  2009-11-26 13:42:52  razeto
 *  Moved to cycle_13_unstable
 *
 *  Revision 1.21  2009-10-26 11:19:48  ddangelo
 *  trgtype 'muon_tct' (or 'muon_totc') renamed as 'neutron' after real use
 *  trgtype 'muon_mtb' renamed as 'muon' as duplicity with the above one disappeared and 'mtb' is associated to btb_inputs flag
 *
 *  Revision 1.20  2008-12-15 17:13:55  razeto
 *  New cycle (12)
 *
 *  Revision 1.19  2008-10-17 13:41:12  razeto
 *  new development cycle (11)
 *
 *  Revision 1.18  2008-02-27 20:46:13  razeto
 *  new development cycle (10)
 *
 *  Revision 1.17  2008-02-27 20:26:30  razeto
 *  New clasdef(9) version and new cycle version
 *
 *  Revision 1.16  2007-12-07 17:33:02  ddangelo
 *  added loading of muon_general_stting table
 *
 *  Revision 1.15  2007-10-11 10:49:54  razeto
 *  Cycle 8 deployed
 *
 *  Revision 1.14  2007-06-22 15:15:27  razeto
 *  Moved to cycle 7
 *
 *  Revision 1.13  2007-05-07 15:47:29  razeto
 *  Cycle number in root classdef
 *
 *  Revision 1.12  2006/11/05 10:27:30  razeto
 *  Remove some useless include
 *
 *  Revision 1.11  2006/07/20 14:53:34  ludhova
 *  Added FadcChannelMapping
 *
 *  Revision 1.10  2006/01/25 13:23:54  misiaszek
 *  Moved from cmap to simple map (to work with root)
 *
 *  Revision 1.9  2004/11/26 15:25:10  razeto
 *  Thanks to Davide who so kindly underlined a syntactical error, I fixed the "Maintainer" word
 *
 *  Revision 1.8  2004/11/26 14:01:20  razeto
 *  Added Mantainer field
 *
 *  Revision 1.7  2004/11/24 09:46:41  razeto
 *  Moved some prototyping in db_acl from db_*
 *
 *  Revision 1.6  2004/09/27 13:31:06  razeto
 *  Fibre bundle, theta and phi read from database
 *
 *  Revision 1.5  2004/06/10 10:07:42  razeto
 *  Added FADC handling
 *
 *  Revision 1.4  2004/05/20 10:49:08  razeto
 *  Applied Davide's patch for loading muon tables
 *
 *  Revision 1.3  2004/04/26 13:48:29  razeto
 *  Added db_acl to check set calls for calling module to have the right privileges.
 *  Fixed the names of set/get methods.
 *  Modifications decided at the software Paris meeting.
 *
 *  Revision 1.2  2004/04/12 16:08:19  razeto
 *  Added crate delay fixed map
 *
 *  Revision 1.1  2004/04/09 08:05:00  razeto
 *  Added
 *
 */
