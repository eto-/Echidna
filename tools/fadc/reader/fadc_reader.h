// Header file of BOREXINO FADC rawdata to ROOT production program

#include "TApplication.h"
#include "Rtypes.h"
#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TSpectrum.h"
#include "TStopwatch.h"

#include <sstream>
#include <cstdio>
#include <zlib.h>
#include <byteswap.h>

#include "Riostream.h"
#include <cmath>
#include <ctime>
#include <vector>
#include <algorithm>

#include <event/BxEvent.hh>

#include "FWFDHeaders.h"

// ****************************** DEFINES ******************************
#define READER_VERSION  "030314"

#define RAWDATA_PATH    "/bxstorage/fadc/rawdata"
//#define RAWDATA_PATH    "/data1/fadc/rawdata"
#define GPS_SANE_PATH   "/home/lukyanch/fadc"

#define ECHIDNA_STYLE_TIME  0   // 0 - use 1970-based UTC, 1 - use Echidna GPS time format
#define CORRECT_UTC_PPS     0   // 1 - substract PPS to UTC time difference
#define WRITE_GPS_SANE_LIST 1   // 1 - runs with sane GPS timestamps will be writen to fadc_gps_sane.txt


// *********************** FUNCTION DECLARATIONS ***********************

extern "C" void* gzopen64(const char*, const char*);

unsigned int convert(unsigned int);
unsigned short convert16(unsigned short t);
int read_raw_event(gzFile gzFin, TEvent *raw, int rawdata_ver);
int guess_raw(gzFile gzf);

int gps_decode(int verbose, int rawdata_ver, TRawHeader raw_header, TEvent raw, BxFwfd* ev);
int add_to_gps_sane_list(string RunName);


// ****************************** CONSTANTS ****************************

const double WND_LEN_NS = NTICK * NS_PER_TICK;

const int CHANNEL_CHECK_MIN_EVNUM       = 1000; // Minimum number of events in run to run a channel check
const int CHANNEL_CHECK_MIN_DIF         = 30;   // Lower difference between min and max value in channel data
                                                // in muon event will increase suspicious counter
const double CHANNEL_CHECK_MISS_RATE    = 0.5;  // Max rate of suspicious events in channel to consider it Ok

const int NOCHECK_CH_CNT                = 7;    // Number of channels ignored on channel check
// this channels will be ignored in channel check:
const int nocheck_channels[NOCHECK_CH_CNT] = {2, 65, 67, 68, 69, 70, 71};
