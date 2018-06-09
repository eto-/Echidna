#define Unix0           0x495C0780  //Time in second from 01/01/1970 00:00 up to 31/12/2008 24:00
#define TAU             195.0       // Time constant for Deconvolution algo

/* FADC Raw event format */
#define KEY1            ADC_MAGIC     // keyword 1
#define KEY2            EVENT_MAGIC   // keyword 2
#define NTICK           512           // number of digitizations
#define NS_PER_TICK     2.5           // ns in one ADC count
#define NBUF            128           // probably not used
#define NCH             102           // number of FADC channels
#define SIZE            64            //
#define NMOD            34            // number of FADC modules
#define NBASE           8             // number of pedestals

#define RAW_DATA_FORMATS 4  // number of rawdata formats to date

// Significant run numbers
#define RUN_FIRST_FADC          12006   // first FADC run synced to LABEN
#define RUN_NEW_BTB_FIRMWARE    12422
#define RUN_CNGS_START          18016   // first with CNGS signals
#define RUN_FIRST_RAW2          12302   // first run with second rawdata format (TU timers introduced)
#define RUN_FIRST_RAW3          17128   // first run with third rawdata format (GPS stamp introduced)
#define RUN_FIRST_RAW4          21212   // first run with fourth rawdata format (header introduced)
#define RUN_TEST_RAW4           505     // first test run with fourth rawdata format

//Shared lib names
#define EV1_LIB "librawevent_1.so"
#define EV2_LIB "librawevent_2.so"
#define EV3_LIB "librawevent_3.so"
