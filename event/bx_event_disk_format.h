#ifndef _BX_EVENT_DISK_FORMAT_H
#define _BX_EVENT_DISK_FORMAT_H

struct event_header_disk_format { /* 40 bytes */
  uint32_t event_size_bytes;
  uint32_t run_number;
  uint32_t event_number;
  uint32_t enabled_crates;
  uint32_t laben_offset;
  uint32_t fadc_offset;
  uint32_t muon_offset;
  uint32_t mctruth_offset;
  uint32_t builder_time_seconds;
  uint32_t builder_time_nanoseconds;
};

struct trigger_disk_format { /* 40 bytes */
  uint16_t length;
  uint16_t error;
  uint16_t blocks;
  uint16_t version;
  uint32_t  dsp_registers;
  uint16_t btb_firmware;
  uint16_t btb_threshold;
  uint32_t  trg_time;
  uint32_t  evid;
  uint8_t  trgtype;
  uint8_t  btb_inputs;
  uint16_t tab_sum;
  uint32_t  gps1;
  uint32_t  gps2;
  uint32_t  gps3;
};

struct laben_header_disk_format { /* 8 bytes */
  uint32_t  entries;
  uint32_t  errors;
};

struct laben_hit_disk_format { /* 12 bytes */
  uint16_t channel;
  uint8_t  time_s[2];
  uint16_t time_l;
  uint8_t  charge[2];
  uint16_t flags;
  uint16_t errors;
};

struct muon_header_disk_format { /* 12 bytes */
  uint32_t  nwords;
  uint32_t  trgid;
  uint32_t  error_flag;
};

struct muon_edge_disk_format { /* 4 bytes */
  uint32_t tdc_word;
};

struct mctruth_header_disk_format { /* 8 bytes */
  uint32_t  length;
  uint16_t frames;
  int16_t   trigger_jitter;
  /* >= 1 frames follow */
};

struct mctruth_daughter_disk_format {
  int    id;
  int    pdg;
  double time;
  float  energy;
  float  position[3];
  float  direction[3];
};

struct mctruth_deposit_disk_format {
  int    pdg_parent;
  float  energy;
  float  position[3];
};

struct mctruth_user_disk_format {
  int     i1;
  int     i2;
  float   f1;
  float   f2;
  double  d;
};

struct mctruth_frame_disk_format { /* 64 bytes */
  uint16_t length;
  uint16_t file_id;
  double         elec_event_time;
  int            event_id;           // Event ID
  int            n_sequence;         // Sequence number of isotope in the chain
  int            isotope_coinc;      // 1 if it's a chain
  int            pdg;                // Particle Data Group Code
  double         time;               // Time
  float          energy;             // Kinetic energy
  float          visible_energy;     // Kinetic energy
  float          position[3];      // Position
  float          baricenter[3];    // Baricenter
  float          direction[3];     // Direction
  int            id_npe;             // Number of photoelectrons in the ID (size of the vector hit_id_v)
  int            od_npe;             // Number of photoelectrons in the OD (size of the vector hit_od_v)
  int            n_daughters;        // Number of daughters (size of the vector daughters_v)
  int            n_deposits;         // Number of deposits  (size of the vector deposits_v)
  int            n_users;            // Size of the vector users_v
  int            n_id_photons;       // Number of generated photons in the ID
  int            n_od_photons;       // Number of generated photons in the OD
  /* subframes may follow */
};

struct mctruth_subframe_disk_format { /* > 4 bytes */
  uint16_t length;
  uint16_t type;		/* enum mctruth_subframe_type = 2 */
  uint16_t tracking_iphs[0];	/* in profile_id 0 (alligned at 4 bytes boundary) */
  float 	 tracking_thps[0];
};

/* enum mc_file_mode { 	/\* bit field with 8 bit to be used in the field mode *\/ */
/*   plain = 0, */
/*   overriding_times = 1, */
/*   copying_tracking_hits = 2, */
/* }; */

/* struct mctruth_file_scattered_info_subframe_disk_format { /\* 100 bytes *\/ */
/*   uint16_t length; */
/*   uint16_t info_type;		/\* enum mctruth_subframe_type = 1 *\/ */
/*   uint16_t file_id; */
/*   uint8_t  mode;                  /\* uses enum mctruth_file_mode *\/ */
/*   char 		 filename[77]; */
/*   uint32_t  geneb_nrun; */
/*   float		 geneb_tempostat; */
/*   uint32_t  geneb_nevents; */
/*   float		 overriding_rate; */
/* }; */

#endif

