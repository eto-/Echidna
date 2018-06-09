// ****************** Rawdata related functions ******************

#include "fadc_reader.h"

// Read raw event from gziped rawdata file and convert to latest format if needed
int read_raw_event(gzFile gzFin, TEvent *raw, int rawdata_ver)
{
    // if event is in latest format - just read it
    if (rawdata_ver >= 3)
        return gzread(gzFin, raw, sizeof(*raw));
    else if (rawdata_ver == 2)
    {
        TEvent_ver2 raw_ver2;
        int ret = gzread(gzFin, &raw_ver2, sizeof(raw_ver2));

        // copy data from old event type
	raw->test1 = raw_ver2.test1;
        raw->sys = raw_ver2.sys;

        memcpy(raw->adc, raw_ver2.adc, NMOD*sizeof(TADCEvent));

        raw->trg.fifo80 = raw_ver2.trg.fifo80;
        raw->trg.fifo82 = raw_ver2.trg.fifo82;
        raw->trg.fifo84 = raw_ver2.trg.fifo84;
        raw->trg.fifo86 = raw_ver2.trg.fifo86;
        raw->trg.fifo88 = raw_ver2.trg.fifo88;
        raw->trg.fifo8A = raw_ver2.trg.fifo8A;
        raw->trg.fifo8C = raw_ver2.trg.fifo8C;
        raw->trg.fifo8E = raw_ver2.trg.fifo8E;

        raw->trg.AsLo = raw_ver2.trg.AsLo;
        raw->trg.AsMi = raw_ver2.trg.AsMi;
        raw->trg.AsHi = raw_ver2.trg.AsHi;
        raw->trg.DeLo = raw_ver2.trg.DeLo;
        raw->trg.DeMi = raw_ver2.trg.DeMi;
        raw->trg.DeHi = raw_ver2.trg.DeHi;

        memset(raw->trg.gps_time, 0, GPS_WORDS*sizeof(uint16_t));

        return ret;
    }
    else if (rawdata_ver == 1)
    {
        TEvent_ver1 raw_ver1;
        int ret = gzread(gzFin, &raw_ver1, sizeof(raw_ver1));

        // copy data from old event type
        raw->sys = raw_ver1.sys;
	raw->test1 = raw_ver1.test1;


        memcpy(raw->adc, raw_ver1.adc, NMOD*sizeof(TADCEvent));

        raw->trg.fifo80 = raw_ver1.trg.fifo80;
        raw->trg.fifo82 = raw_ver1.trg.fifo82;
        raw->trg.fifo84 = raw_ver1.trg.fifo84;
        raw->trg.fifo86 = raw_ver1.trg.fifo86;
        raw->trg.fifo88 = raw_ver1.trg.fifo88;
        raw->trg.fifo8A = raw_ver1.trg.fifo8A;
        raw->trg.fifo8C = raw_ver1.trg.fifo8C;
        raw->trg.fifo8E = raw_ver1.trg.fifo8E;

        raw->trg.AsLo = 0;
        raw->trg.AsMi = 0;
        raw->trg.AsHi = 0;
        raw->trg.DeLo = 0;
        raw->trg.DeMi = 0;
        raw->trg.DeHi = 0;

        memset(raw->trg.gps_time, 0, GPS_WORDS*sizeof(uint16_t));

        return ret;
    }

    return -1;
}


// DAQ machine to production machine endianess convert
unsigned int convert(unsigned int t)  {

  return bswap_32(t);
}


unsigned short convert16(unsigned short t)  {

  return bswap_16(t);
}


// Try to guess raw file type based on keyword position
int guess_raw(gzFile gzf)
{
    const int bufsize = sizeof(TEvent); // the largest of all
    char buf[bufsize];
    TRawHeader *header;
    TEvent *ev3;
    TEvent_ver2 *ev2;
    TEvent_ver1 *ev1;
    
    gzrewind(gzf);
    gzread(gzf, buf, bufsize);
    gzrewind(gzf);
    
    header = (TRawHeader*) buf;
    
    if (convert(header->Magic) == RAWDATA_MAGIC)
        return 4;
    
    ev3 = (TEvent*) buf;
    if (convert(ev3->test1) == EVENT_MAGIC)
        return 3;
    
    ev2 = (TEvent_ver2*) buf;
    if (convert(ev2->test1) == EVENT_MAGIC)
        return 2;
    
    ev1 = (TEvent_ver1*) buf;
    if (convert(ev1->test1) == EVENT_MAGIC)
        return 1;

    // we were unable to determine, probably file is wrong
    return -1;
}
