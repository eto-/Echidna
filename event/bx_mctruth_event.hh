/* BOREXINO Reconstruction program
 *
 * Author: Davide D'Angelo <davide.dangelo@lngs.infn.it>
 * Maintainer: Davide D'Angelo <davide.dangelo@lngs.infn.it>
 *
 * $Id: bx_mctruth_event.hh,v 1.10 2009/07/17 17:05:41 ddangelo Exp $
 *
 * The MonteCarlo Truth object
 * This portion of the event is read by bx_reader and written to file by bx_writer
 * No other module should write to this area.
 * 
 */

#ifndef _BX_MCTRUTH_EVENT_HH
#define _BX_MCTRUTH_EVENT_HH

#include "bx_rec_general.hh"
#include "bx_base_event.hh"
#include "constants.hh"
#include "bx_event_disk_format.h"

class bx_mctruth_hit {
public:
  bx_mctruth_hit ( uint16_t lg, float time);

  uint16_t get_lg () const { return u2_lg; }
  float        get_time () const { return f4_time; } 

private:
  uint16_t u2_lg;
  float          f4_time;
};

class bx_mctruth_daughter {
public:
  bx_mctruth_daughter ( mctruth_daughter_disk_format *raw);

  int    get_id        ()      const { return i4_id;             }
  int    get_pdg       ()      const { return i4_pdg;            }
  double get_time      ()      const { return f8_time;           } 
  float  get_energy    ()      const { return f4_energy;         }
  float  get_position  (int i) const { return f4_position_v[i];  }
  float  get_direction (int i) const { return f4_direction_v[i]; }

private:
  int    i4_id;
  int    i4_pdg;
  double f8_time;
  float  f4_energy;
  float  f4_position_v[3];
  float  f4_direction_v[3];
};

class bx_mctruth_deposit {
public:
  bx_mctruth_deposit ( mctruth_deposit_disk_format *raw);

  int   get_pdg_parent      () const { return i4_pdg_parent;    }
  float get_energy          () const { return f4_energy;        } 
  float get_position ( int i ) const { return f4_position_v[i]; } 

private:
  int   i4_pdg_parent;
  float f4_energy;
  float f4_position_v[3];
};

class bx_mctruth_user {
public:
  bx_mctruth_user ( mctruth_user_disk_format *raw);

  int    get_int1   () const { return i4_int1;   }
  int    get_int2   () const { return i4_int2;   }
  float  get_float1 () const { return f4_float1; }
  float  get_float2 () const { return f4_float2; }
  double get_double () const { return f8_double; } 

private:
  int i4_int1, i4_int2;
  float f4_float1, f4_float2;
  double f8_double;
};

class bx_mctruth_frame {
public:
  enum mctruth_subframe_type { /* the values of the info_type field of every subframe */
    scattered_info = 1, // this is skipped at the moment
    pe_id_list     = 2,
    pe_od_list     = 3
//    deposits_list = 4,
//    daughters_list = 5,
//    users_list = 6
};

  bx_mctruth_frame ( const char* disk_event);

  uint16_t get_file_id () const { return u2_file_id;         }
  double get_elec_event_time () const { return f8_elec_event_time; }
  int    get_event_id        () const { return i4_event_id;        }
  int    get_n_sequence      () const { return i4_n_sequence;      }
  int    get_isotope_coinc   () const { return i4_isotope_coinc;   }
  int    get_pdg             () const { return i4_pdg;             }
  double get_time            () const { return f8_time;            }
  float  get_energy          () const { return f4_energy;          }
  float  get_visible_energy  () const { return f4_visible_energy;  }
  float  get_position   (int i) const { return f4_position_v  [i]; }
  float  get_baricenter (int i) const { return f4_baricenter_v[i]; }
  float  get_direction  (int i) const { return f4_direction_v [i]; }
  int    get_id_npe          () const { return i4_id_npe;          } 
  int    get_od_npe          () const { return i4_od_npe;          } 
  int    get_n_daughters     () const { return i4_n_daughters;     } 
  int    get_n_deposits      () const { return i4_n_deposits;      } 
  int    get_n_users         () const { return i4_n_users;         } 
  int    get_n_id_photons    () const { return i4_n_id_photons;    } 
  int    get_n_od_photons    () const { return i4_n_od_photons;    } 

  const bx_mctruth_hit&      get_hit_id   (int i) const { return hit_id_v[i];    }
  const bx_mctruth_hit&      get_hit_od   (int i) const { return hit_od_v[i];    }
  const bx_mctruth_daughter& get_daughter (int i) const { return daughters_v[i]; }
  const bx_mctruth_deposit&  get_deposit  (int i) const { return deposits_v[i];  }
  const bx_mctruth_user&     get_user     (int i) const { return users_v[i];     }

private:
  uint16_t u2_file_id;
  double         f8_elec_event_time;    // relative to trg time
  int            i4_event_id;           // Event ID
  int            i4_n_sequence;         // Sequence number of isotope in the chain
  int            i4_isotope_coinc;      // 1 if it's a chain
  int            i4_pdg;                // Particle Data Group Code
  double         f8_time;               // Time
  float          f4_energy;             // Kinetic energy
  float          f4_visible_energy;     // Quenched energy
  float          f4_position_v[3];      // Position
  float          f4_baricenter_v[3];    // Baricenter
  float          f4_direction_v[3];     // Direction
  int            i4_id_npe;             // Number of photoelectrons in the ID (size of the vector hit_id_v)
  int            i4_od_npe;             // Number of photoelectrons in the OD (size of the vector hit_od_v)
  int            i4_n_daughters;        // Number of daughters (size of the vector daughters_v)
  int            i4_n_deposits;         // Number of deposits  (size of the vector deposits_v)
  int            i4_n_users;            // Size of the vector users_v
  int            i4_n_id_photons;       // Number of generated photons in the ID
  int            i4_n_od_photons;       // Number of generated photons in the OD
  std::vector<bx_mctruth_hit>      hit_id_v;
  std::vector<bx_mctruth_hit>      hit_od_v;
  std::vector<bx_mctruth_daughter> daughters_v;
  std::vector<bx_mctruth_deposit>  deposits_v;
  std::vector<bx_mctruth_user>     users_v;
};

class bx_mctruth_event: 
  public bx_base_event {
  public:
    bx_mctruth_event (const char *disk_event);
    virtual ~bx_mctruth_event () {}

    bool is_data() const { return b_is_data; }
    int16_t get_trigger_jitter   () const { return i2_trigger_jitter; }
    uint16_t get_nframes () const { return frame_v.size(); }
    const bx_mctruth_frame& get_frame (int i) const { return frame_v[i]; }

  private:
    bool b_is_data;    
    int16_t i2_trigger_jitter;
    std::vector<bx_mctruth_frame> frame_v;
};

#endif
/*
 * $Log: bx_mctruth_event.hh,v $
 * Revision 1.10  2009/07/17 17:05:41  ddangelo
 * added a flag for real data
 *
 * Revision 1.9  2007-11-06 14:45:11  ddangelo
 * debugging mctruth
 *
 * Revision 1.8  2007-10-31 15:42:14  ddangelo
 * added quenched energy (mctruth)
 * disk format, internal and root event
 *
 * Revision 1.7  2007-10-29 16:28:49  ddangelo
 * reading new format.
 * subframe types for daughters, deposits and users abandoned.
 * Hierarchy assumed instead.
 * file_id added.
 * hits reintroduced, sorting added.
 *
 * Revision 1.6  2007-10-11 11:23:53  ddangelo
 * new mctruth format (internal AND root event)
 *
 * Revision 1.5  2005-05-11 14:00:29  ddangelo
 * fixed a bug
 *
 * Revision 1.4  2005/01/19 13:15:44  ddangelo
 * added elec_event_time
 *
 * Revision 1.3  2004/11/26 15:25:10  razeto
 * Thanks to Davide who so kindly underlined a syntactical error, I fixed the "Maintainer" word
 *
 * Revision 1.2  2004/11/26 14:06:22  razeto
 * Added Mantainer field
 *
 * Revision 1.1  2004/07/12 17:03:56  ddangelo
 * added mctruth.
 * Still undebugged, tmp commit to allow debug.
 *
 *
 */
