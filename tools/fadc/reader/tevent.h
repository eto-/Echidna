//=== Borexino FADC rawdata file format =========================

#include <stdint.h>

#define RAWDATA_MAGIC   0x46414443
#define EVENT_MAGIC     0xAAAA5555
#define ADC_MAGIC       0xBEEFBEEF
#define ERROR_MAGIC     0x33333333

#define FLAG_PPS_ENA    0x00000001

#define GPS_WORDS       5
#define MAX_ADC         34

typedef uint32_t THalfChannel[64];

// Event info produced by PPC
typedef struct {
    uint32_t Version;
    uint32_t UnixTime; // 32 bits!
    uint32_t Counter;
    uint32_t Error;
} TSysEvent;

// Event info produced by V1495 trigger unit
typedef struct {
    // Trigger unit FIFO data
    uint16_t fifo80, fifo82, fifo84, fifo86, fifo88, fifo8A, fifo8C, fifo8E;
    // Trigger unit astronomical/dead timers
    uint16_t AsLo, AsMi, AsHi, DeLo, DeMi, DeHi;
    // GPS timestamp (SLAV2000 or PPS)
    uint16_t gps_time[GPS_WORDS];
} TTrgEvent;

// Per ADC event
typedef struct {
    THalfChannel HalfChannel[6];
    uint32_t LoHiPage;
    uint32_t TimeLog;
    uint32_t Pattern;
    uint32_t test2;
} TADCEvent;

// FADC rawdata event
typedef struct {
    TSysEvent sys;
    TTrgEvent trg;
    TADCEvent adc[MAX_ADC];
    uint32_t test1;
} TEvent;

// FADC rawdata file header
typedef struct {
    uint32_t Magic;
    uint32_t Online_ver;
    uint32_t Run_num;
    uint32_t Start_time;
    uint32_t flags;
    uint32_t pps_diff;
    uint32_t Trg_FW_ver;
    uint32_t reserved0;
} TRawHeader;


/* ************************* OLD EVENT TYPES  ************************ */

// Event info produced by V1495 trigger unit (ver 2)
typedef struct {
    // Trigger unit FIFO data
    uint16_t fifo80, fifo82, fifo84, fifo86, fifo88, fifo8A, fifo8C, fifo8E;
    // Trigger unit astronomical/dead timers
    uint16_t AsLo, AsMi, AsHi, DeLo, DeMi, DeHi;
} TTrgEvent_ver2;

// FADC rawdata event ver 2
typedef struct {
    TSysEvent sys;
    TTrgEvent_ver2 trg;
    TADCEvent adc[MAX_ADC];
    uint32_t test1;
} TEvent_ver2;

// Event info produced by V1495 trigger unit (ver 1)
typedef struct {
    // Trigger unit FIFO data
    uint16_t fifo80, fifo82, fifo84, fifo86, fifo88, fifo8A, fifo8C, fifo8E;
} TTrgEvent_ver1;

// FADC rawdata event ver 1
typedef struct {
    TSysEvent sys;
    TTrgEvent_ver1 trg;
    TADCEvent adc[MAX_ADC];
    uint32_t test1;
} TEvent_ver1;
