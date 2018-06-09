/**
 * @file Parameters.h
 * @brief Contains constant parameters used by the various modules of the 
 * @brief reconstruction program. Can also contain variables that need to
 * @brief be globally accessible.
 * @author Kevin, Tibi, Ben, Richard, Josh
 */
#ifndef _PARAMETERS_H
#define _PARAMETERS_H

#include <stdint.h>

/**
 * @defgroup params Global parameters
 * @{
 */
/**
 *   @defgroup gen_params General parameters
 *   @{
 */
/** 
 * @brief Number of events saved in circular buffers of InputReader and 
 * DataProcessor
 */
#define CIRCULAR_BUFFER_SIZE 1000

/** 
 * @brief Maximum expected size of raw event on disk, in bytes
 * If the event is any larger we will do a (slow) dynamic memory allocation
 * to hold it.
 */
#define EXPECTED_MAX_EVENT_SIZE 350000

/** 
 * @brief Maximum size of raw event we are willing to believe, in bytes
 */
#define ABSOLUTE_MAX_EVENT_SIZE 1000000

/** 
 * @brief Minimum run number to expect
 */
#define MIN_RUN_NUMBER 1000

/**
 * @brief Maximum run number to expect
 */
#define MAX_RUN_NUMBER 999999

/** 
 * @brief Frequency with which to print that we are processing event @a N
 */
#define EVENTS_BETWEEN_OUTPUT_MARK 2500

/** 
 * @brief Minimum number of calibration events required per run
 */
#define MIN_CALIBRATION_EVENTS 400 

/**
 * @brief Desired number of calibration events; if less read from previous
 */
#define DESIRED_CALIBRATION_EVENTS 9950

#define OLDEST_CALIBRATION_RUN 100

#ifndef DEFAULT_CFGFILE
/** 
 * @brief Name of config file to use
 */
#define DEFAULT_CFGFILE "bx.cfg"
#endif

/** 
 * @brief Default verbosity
 */
#define DEFAULT_VERBOSITY DEBUG

/** 
 * @brief Minimum PMT channel number
 */
#define MIN_PMT 0
/**
 * @brief Maximum PMT channel number
 */
#define MAX_PMT 2300

/** 
 * @brief Minimum muon channel number
 */
#define MIN_Muon 3000
/**
 * @brief Maximum muon channel number
 */
#define MAX_Muon 3300

/**
 * @brief Minimum FADC channel number
 */
#define MIN_FADC 4000
/**
 * @brief Maximum FADC channel number
 */
#define MAX_FADC 4200

/** 
 * @brief Maximum value of a channel number (PMT, muon or FADC)
 */
#define MAX_Channel 4500

/**
 *   @defgroup recon Position Reconstruction parameters
 *   @{
 */
/// @name High time cut to use for raw hits wrt. first event time
//   @{
#define RAW_TIME_HIGH_CUT   54
//   @}

/// @name Low and high cuts for PDF to use during position recon wrt. peak time
//   @{
#define PDF_LOW_CUT   -5
#define PDF_HIGH_CUT   15
//   @}

/** @brief Time between hits to define start of cluster, ns */
#define		 FIRST_HIT_DT 3

/** @brief Speed of light in a vacuum, in meters/ns */
#define	   SPEED_OF_LIGHT 0.299792458

/** @brief Index of refraction of pseudocumene. */
#define    IDX_REF_PC 1.53  /* n at 390 nm */

/** @brief Z-offset of Inner Vessel center from center of PMTs, meters. */
#define    VESSEL_Z_CENTER 0.20

/** @brief Average radius of Inner Vessel, meters. */
#define    VESSEL_RADIUS 4.25

/** @brief Average radius of PMT (center of IV to front of photocathode?), meters */
#define    PMT_RADIUS 6.559

/** @brief Radius of light concentrator at entry apeture, meters */
#define    LIGHT_CONC_RADIUS 0.16

/** @brief Optical coverage of inner detector, fractional */
#define OPTICAL_COVERAGE 0.32

/** @brief PMT Detection efficiency */
#define DET_EFFICIENCY 0.16

/** @brief Attenuation length in scintillator at 440 nm, meters */
#define SCINT_ATTENUATION 12.0

/** @brief Attenuation length in buffer at 440 nm, meters */
#define BUFFER_ATTENUATION 2.9
// @}

/** 
 * @brief Command line arguments recognized by the program (getopt syntax)
 */
#define CMDLINE_ARGS ":vq-:hpsbcd:f:l:"

/**
 * @brief Maximum number of calibration / processing modules that can
 * be defined in the DataProcessor.
 */
#define MAX_MODULES 25

class CommandLine;
/**
 * @brief The command line parser object. Might change from a global variable.
 */
extern CommandLine *cmdlineparser;
//   @}

/** 
 *   @defgroup muon_hardware Muon Hardware Parameters
 *   @{
 */
const int muon_chips_per_board = 4; ///< Number of muon chips per board
const int muon_channels_per_chip = 32;///< Number of channels per chip
const int muon_channel_offset = 3000;///< Offset for muon channel numbers
const int muon_time_buffer = 1<<13;///< Muon time buffer @todo More explicitly?
const int muon_num_channels=256; ///< Number of muon channels
const double muon_ns_per_tick=1.0416;
///< Muon ns per tick
const double muon_pe_per_tick=0.07295;
///< Muon PE per tick
//   @}

/**
 *   @defgroup muon_cal Muon Calibration Parameters
 *   @{
 */
const bool muon_use_mean = true;
///< Use mean or mode in muon pedestal calibration
const int muon_pulse_tolerance = 128;///< Pulse tolerance @todo More explicitly?
const int muon_min_pedestal_width = 610;///< Minimum width of muon pedestal
const int muon_max_pedestal_width = 810;///< Maximum width of muon pedestal

// Muon Charge and Time Calibration Parameters
const double muon_min_efficiency = .005;
///< Determines minimum number of hits to consider a channel for calibration
const int muon_time_fit_range = 20;
///< Time fit range for muon calibration @todo Explain this?
const int muon_time_max = 5;
///< Muon time max
const int muon_time_rms_min = 0;
///< Muon time RMS min
const int muon_time_rms_max = 5;
///< Muon time RMS max
const int muon_charge_max = 30;
///< Muon chage max
const int muon_charge_rms_min = 1;
///< Muon charge RMS min
const int muon_charge_rms_max = 10;
///< Muon charge RMS max
//   @}

/**
 *   @defgroup laben Laben Parameters
 *   @{
 */
#define N_LABEN_CHANNELS 2240
///< The total number of logical laben channels
#define N_LABEN_CRATES   14
///< The number of laben crates
#define N_CHANNELS_PER_CRATE (N_LABEN_CHANNELS / N_LABEN_CRATES)
///< Number of channels per laben crate
#define N_CHANNELS_PER_BOARD 8
///< Number of channels on each laben board
#define N_LABEN_BOARDS (N_LABEN_CHANNELS / N_CHANNELS_PER_BOARD)
///< Total number of laben boards in the inner detector
#define MIN_HITS_PER_CRATE_NEEDED 4
///< Num of hits needed for a crate to count as triggering in the cluster module
#define MIN_CRATES_NEEDED 3
///< Num of triggering crates needed for the cluster module to accept event
#define N_REFERENCE_CHANNELS 14
///< The number of dedicated trigger reference channels
const int laben_reference_channel[N_REFERENCE_CHANNELS]=
  {126,297,464,604,719,939,1084,1221,1324,1453,1627,1802,1955,2105};
///< Array containing the reference channel numbers
#define N_LASER_CHANNELS 6
///< The number of dedicated laser reference channels
const int laser_reference_channel[N_LASER_CHANNELS]=
  {130,266,1086,1457,1634,2104};
///< Array containing the laser reference channels
#define LABEN_REFERENCE_WIDTH 50
///< Maximum width of reference times before convergence
#define INTEGRATOR_DECAY_TIME 500
///< Decay time of ADC in ns
#define NOMINAL_LASER_PEAK 650
///< Expected center of the laser peak relative to the triggers
#define GRAY_WINDOW 3276800
///< Width of possible PMT hit times in ns for matching hits to triggers
//   @}

/**
 *   @defgroup splitting  Splitting Parameters
 *   @{
 */
const int n_pmts_ref = 2000;
///< Number of PMT's to normalize hits
const int dark_bin_size = 1600;
///< Size of bins for fitting dark noise (in ns)
const int cluster_bin_size = 16;
///< Size of bins for finding clusters (in ns)
const int cluster_start_bins = 3;
///< Number of bins in start_gate
const double  cluster_start_threshold = 3;
///< Average number of hits per bin - used as start threshold
const int cluster_end_bins = 6;
///< Number of bins in end_gate
const double cluster_end_sigma = 3;
///< Number of sigma above dark noise - used for end threshold
const double neutron_start_sigma = 3;
///< Number of sigma above dark noise - used for start threshold in neutron clusters
const double neutron_end_sigma = 1;
///< Number of sigma above dark noise - used for end threshold in neutron clusters
const int cluster_tail_bins = 32;
///< Minimum number of hits to add to end of cluster. IMP: Do NOT reduce this number - needed for good Gatti parameter separation between alphas and betas
const int cluster_integral = 20;
///< Minimum number of hits in each cluster
const int cluster_log_en = 100;
///< Hit value above which a longer tail is added
const int neutron_actual_start_limit = 20;
///< Maximum number of hits (from first hit) to look for actual start of event
const int neutron_actual_start_width = 32;
///< Size of gate to look for actual start
const int neutron_actual_start_threshold = 3;
///< Number of hits that have to fall within "neutron_actual_start_width" to determine actual start
//   @}

/**
 *   @defgroup G4BX  Parameters used for incorporating montecarlo
 *   @{
 */
#define G4_N_CHANNELS 6000
///< The number of channles in G4BX - not sure where this number comes from
#define G4_CHANNEL_INT_TIME 80
///< Integration time of channels in ns.
#define G4_CHANNEL_DEAD_TIME_START 80
///< Start of dead time of channel, in ns, measured from first hit.
#define G4_CHANNEL_DEAD_TIME_END 140
///< End of dead time of channel, in ns, measured from first hit.
// @}

#endif
