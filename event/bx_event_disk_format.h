#ifndef _BX_EVENT_DISK_FORMAT_H
#define _BX_EVENT_DISK_FORMAT_H

struct event_header_disk_format { /* 40 bytes */
  unsigned long event_size_bytes;
  unsigned long run_number;
  unsigned long event_number;
  unsigned long enabled_crates;
  unsigned long laben_offset;
  unsigned long fadc_offset;
  unsigned long muon_offset;
  unsigned long mctruth_offset;
  unsigned long builder_time_seconds;
  unsigned long builder_time_nanoseconds;
};

struct trigger_disk_format { /* 40 bytes */
  unsigned short length;
  unsigned short error;
  unsigned short blocks;
  unsigned short version;
  unsigned long  dsp_registers;
  unsigned short btb_firmware;
  unsigned short btb_threshold;
  unsigned long  trg_time;
  unsigned long  evid;
  unsigned char  trgtype;
  unsigned char  btb_inputs;
  unsigned short tab_sum;
  unsigned long  gps1;
  unsigned long  gps2;
  unsigned long  gps3;
};

struct laben_header_disk_format { /* 8 bytes */
  unsigned long  entries;
  unsigned long  errors;
};

struct laben_hit_disk_format { /* 12 bytes */
  unsigned short channel;
  unsigned char  time_s[2];
  unsigned short time_l;
  unsigned char  charge[2];
  unsigned short flags;
  unsigned short errors;
};

struct muon_header_disk_format { /* 12 bytes */
  unsigned long  nwords;
  unsigned long  trgid;
  unsigned long  error_flag;
};

struct muon_edge_disk_format { /* 4 bytes */
  unsigned long tdc_word;
};

struct mctruth_header_disk_format { /* 8 bytes */
  unsigned long  length;
  unsigned short frames;
  signed short   trigger_jitter;
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
  unsigned short length;
  unsigned short file_id;
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
  unsigned short length;
  unsigned short type;		/* enum mctruth_subframe_type = 2 */
  unsigned short tracking_iphs[0];	/* in profile_id 0 (alligned at 4 bytes boundary) */
  float 	 tracking_thps[0];
};

/* enum mc_file_mode { 	/\* bit field with 8 bit to be used in the field mode *\/ */
/*   plain = 0, */
/*   overriding_times = 1, */
/*   copying_tracking_hits = 2, */
/* }; */

/* struct mctruth_file_scattered_info_subframe_disk_format { /\* 100 bytes *\/ */
/*   unsigned short length; */
/*   unsigned short info_type;		/\* enum mctruth_subframe_type = 1 *\/ */
/*   unsigned short file_id; */
/*   unsigned char  mode;                  /\* uses enum mctruth_file_mode *\/ */
/*   char 		 filename[77]; */
/*   unsigned long  geneb_nrun; */
/*   float		 geneb_tempostat; */
/*   unsigned long  geneb_nevents; */
/*   float		 overriding_rate; */
/* }; */

#endif

