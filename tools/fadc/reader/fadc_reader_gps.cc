// ****************** GPS-related functions for FADC reader ******************

#include "fadc_reader.h"

// Decode GPS time of event and check that it is sane
int gps_decode(int verbose, int rawdata_ver, TRawHeader raw_header, TEvent raw, BxFwfd* ev)
{
    static unsigned long secev, prev_secev;
    static unsigned long nsecev, prev_nsecev;
    int gps_errors = 0;

    if (rawdata_ver == 3)
    {
        int year, mnth, day, hour, min, sec, msec, usec, nsec;

        // Decode GPS time from SLAV2000 format
        nsec = (convert16(raw.trg.gps_time[0]) & 0x000F) * 100;
        usec = ((convert16(raw.trg.gps_time[0]) & 0x00F0) >> 4) + ((convert16(raw.trg.gps_time[0]) & 0x0F00) >> 8)*10 + ((convert16(raw.trg.gps_time[0]) & 0xF000) >> 12)*100;
        msec = (convert16(raw.trg.gps_time[1]) & 0x000F) + ((convert16(raw.trg.gps_time[1]) & 0x00F0) >> 4)*10 + ((convert16(raw.trg.gps_time[1]) & 0x0F00) >> 8)*100;
        sec  = (convert16(raw.trg.gps_time[2]) & 0x000F) + ((convert16(raw.trg.gps_time[2]) & 0x0070) >> 4)*10;
        min  = ((convert16(raw.trg.gps_time[2]) & 0x0F00) >> 8) + ((convert16(raw.trg.gps_time[2]) & 0x7000) >> 12)*10;
        hour = (convert16(raw.trg.gps_time[3]) & 0x000F) + ((convert16(raw.trg.gps_time[3]) & 0x0030) >> 4)*10;
        day  = ((convert16(raw.trg.gps_time[3]) & 0x0F00) >> 8) + ((convert16(raw.trg.gps_time[3]) & 0x3000) >> 12)*10;
        mnth = (convert16(raw.trg.gps_time[4]) & 0x000F) + ((convert16(raw.trg.gps_time[4]) & 0x0010) >> 4)*10;
        year = ((convert16(raw.trg.gps_time[4]) & 0x0F00) >> 8) + ((convert16(raw.trg.gps_time[4]) & 0xF000) >> 12)*10;

        // compute number of seconds since 01-01-2000
        struct tm working_date;
        working_date.tm_sec = sec;
        working_date.tm_min = min;
        working_date.tm_hour = hour - 1;     // ugly hack to fix DAQ timezone
        working_date.tm_mday = day;
        working_date.tm_mon = mnth-1;
        working_date.tm_year = (year > 90) ? year : (year+100);
        working_date.tm_isdst = 0;

        time_t timet = timegm(&working_date);

        // patch for gps overflow error occurred at run 16789 (2011-10-09)
        if (year > 90) {
            timet += 1024*7*86400;
            time_t tmp_timet = timet;
            if (hour == 0) tmp_timet += 3600;
            struct tm *fixed_date = gmtime(&tmp_timet);
            day  = fixed_date->tm_mday;
            mnth = fixed_date->tm_mon+1;
            year = fixed_date->tm_year;
        }

        secev = timet;
        nsecev = (unsigned long)(msec*1000000 + usec*1000 + nsec);
    }
    else if ((rawdata_ver >= 4) && (convert(raw_header.flags) & FLAG_PPS_ENA))
    {
        long sec;
        double subsec_clk;

        const double REF_CLK_FREQ = 50000000.;

        sec = convert(raw_header.Start_time) + convert16(raw.trg.gps_time[2]) + 1;
        double PPS_clk_cor;

        PPS_clk_cor = REF_CLK_FREQ / (convert16(raw.trg.gps_time[3]) +
                             (convert16(raw.trg.gps_time[4]) << 16));

        if ((PPS_clk_cor > 1.01) || (PPS_clk_cor < 0.99))
        {
            if (verbose)
            {
                cerr << "Insane PPS clock correction coefficient: " << PPS_clk_cor;
                cerr << "(" << hex << convert16(raw.trg.gps_time[4]) << convert16(raw.trg.gps_time[3]) << ")" << dec;
                cerr << " in event " << ev->GetEvNum() << endl;
            }
            gps_errors++;
        }

        subsec_clk = (convert16(raw.trg.gps_time[0])*1. + (convert16(raw.trg.gps_time[1]) << 16)*1.) *
                     PPS_clk_cor;

        if ((subsec_clk > REF_CLK_FREQ*1.001) || (subsec_clk < 0))
        {
            if (verbose)
            {
                cerr << "Insane subsecond PPS counter: " << (long long)subsec_clk;
                cerr << " in event " << ev->GetEvNum() << endl;
            }
            gps_errors++;
        }

        #if CORRECT_UTC_PPS
        double pps_diff = convert(raw_header.pps_diff);

        if ((pps_diff < 0) || (pps_diff > REF_CLK_FREQ))
        {
            if (verbose)
            {
                cerr << "Insane PPS<->UTC difference: " << (long long)pps_diff;
                cerr << " in event " << ev->GetEvNum() << endl;
            }
            gps_errors++;
        }

        // compensate PPS<->UTC difference
        if (subsec_clk < pps_diff)
        {
            sec--;
            subsec_clk = REF_CLK_FREQ - (pps_diff - subsec_clk);
        }
        else
            subsec_clk -= pps_diff;
        #endif

        secev = sec;
        nsecev = (unsigned long)(subsec_clk * (1000000000L / REF_CLK_FREQ)); // nanoseconds
    }
    else
    {
        // there is no GPS timestamps in this rawfile
        ev->SetGpsTimeSec(0);
        ev->SetGpsTimeNs(0);
        return 1;
    }

    // check that PPC read event later than it occured
    long unixtime = ev->GetUnixTime();

    if (((unixtime - (long)secev) > 100) ||
        (((long)secev - unixtime) > 1)  ||
        ((((long)secev - unixtime) > 0) && (nsecev > 300000000)))
    {
        if (verbose)
        {
            cerr << "GPS seconds (" << secev << ") inconsistent with DAQ unix time (" << unixtime << ")";
            cerr << " in event " << ev->GetEvNum() << endl;
        }
        gps_errors++;
        
        secev = 0;
        nsecev = 0;
    }

    // check that GPS timestamps are continuous
    if ((secev < prev_secev) || ((secev == prev_secev) && ((long long)nsecev <= (long long)prev_nsecev)))
    {
        if (verbose)
        {
            cerr << "GPS time leap! Previous: " << prev_secev << "." << prev_nsecev;
            cerr << " Current: " << secev << "." << nsecev;
            cerr << " in event " << ev->GetEvNum() << endl;
        }
        gps_errors++;
        
        secev = 0;
        nsecev = 0;
    }

    long long prevtime_gps = (((long long)(secev-prev_secev) * 1000000000L) + (nsecev-prev_nsecev));

    if (((prev_secev - secev) < 4) && (abs(double(prevtime_gps - (long long)ev->GetRawTime())) > 100000))
    {
        if (verbose)
        {
            cerr << "Big prevtime diff: " << (long long)prevtime_gps - (long long)ev->GetRawTime();
            cerr << " rawtime: " << (long long)ev->GetRawTime() << " GPS: " << (long long)prevtime_gps;
            cerr << " in event " << ev->GetEvNum() << endl;
        }
        gps_errors++;
        
        secev = 0;
        nsecev = 0;
    }

    if (secev > 1577836800L) // year 2020? :)
    {
        if (verbose)
        {
            cerr << "Insane GPS seconds: " << secev;
            cerr << " in event " << ev->GetEvNum() << endl;
        }
        gps_errors++;
        
        secev = 0;
        nsecev = 0;
    }

    prev_secev = secev;
    prev_nsecev = nsecev;

    #if ECHIDNA_STYLE_TIME
    // convert to strange Echidna-style GPS time format
    if (secev > 1230768000L) secev += 2;        //31 Dec 2008 2 leap seconds
    else if (secev > 1136073600L) secev += 1;   //1 Jan 2005 1 leap second

    if (secev > 0) secev -= 946684800L;         // going from 1970 to 2000 based number
    #endif

    // write time to event
    ev->SetGpsTimeSec(secev);
    ev->SetGpsTimeNs(nsecev);

    return gps_errors;
}


// Add run to FADC list of runs with sane GPS times
int add_to_gps_sane_list(string RunName)
{
    const string GPS_SANE_LIST = GPS_SANE_PATH+(string)"/fadc_gps_sane.txt";
    const string GPS_INSANE_LIST = GPS_SANE_PATH+(string)"/fadc_gps_insane.txt";

    // check if run is already present in list

    ifstream infile(GPS_SANE_LIST.c_str());

    if (!infile)
        return -1;

    string str;

    while (infile >> str)
        if (str == RunName)
            return 1;

    infile.close();
    infile.open(GPS_INSANE_LIST.c_str());

    if (!infile)
        return -1;

    while (infile >> str)
        if (str == RunName)
            return 2;

    // write to sane list
    ofstream outfile;

    outfile.open(GPS_SANE_LIST.c_str(), std::ios_base::app);
    outfile << RunName << endl;
    return 0;
}
