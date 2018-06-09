/***************** Reading FWFD event.
Starting from Echidna cycle_13, rootechidna.so is required.
         mailto: litvinovich@lngs.infn.it ****************/

/* Revisions:
 *
 * 17.10.13 by Georgy Lukyanchenko:
 * Added partial run production with -e option and channel check
 *
 * 23.10.13 by Georgy Lukyanchenko:
 * Added new PPS+NTP timing decode and GPS sanity checks
 *
 * 27.10.13 by Georgy Lukyanchenko:
 * Added ability to process older FADC rawdata formats
 *
 * 11.02.14 by Georgy Lukyanchenko:
 * Small fixes
 *
 * 24.02.14 by Georgy Lukyanchenko:
 * Small changes in clustering and waveform decoding (small rawtime fix)
 *
 * 03.03.14 by Georgy Lukyanchenko:
 * Small fixes for reprocessing
 *
 * mailto: lukyanch@lngs.infn.it
 */


using namespace std;


#include "fadc_reader.h"


// Print usage help
void usage(char *name) {
    cout << "Usage: " << name << " RUN_NUMBER [-v] [-r RAWDATA_FORMAT] [-e LAST_EVENT_NUMBER]" << endl;
    exit(0);
}


// Main function
int main(int argc, char **argv)  {

    TEvent raw;

    int rawdata_ver = 0;

    int Last_event = 0;
    int Run_number = 0;

    unsigned int qq, res;
    unsigned char nclu=0;
    int i, j, k, m;
    int evnum = 0, prev_evnum;
    int current_ChID;
    unsigned int peak_pos0=0, peak_pos[6], peak_posD[6];
    int cluster_start[6], cluster_end[6], cluster_startD = 0, cluster_endD = 0;
    float Q[6], atten_Q[6], QD[6], peak_amp, peak_ampl[6];
    double time_rest=0.; // ms
    int gps_errors = 0;
    int overfl_errors = 0;

    #if CYCLE >= 16
    double time_prev = 0.0;
    #else
    unsigned int time_prev = 0;
    #endif

    bool pile_up[6];

    unsigned char wf[NMOD][3][NTICK];
    int susp_channels[NCH] = {0};

    int verbose = 0;
    int susp_checks = 0;

    Float_t *source = new float[NTICK];
    Float_t *dest = new float[NTICK];

    vector<float> base(NBASE); // 8
    vector<unsigned int> wform(NTICK);
    vector<float> wform_based(NTICK);
    vector<float> wform_based_sum(NTICK);
    vector<float> dd(NTICK);

    // Parse commandline parameters
    if ((argc < 2) || (argc > 7))
        usage(argv[0]);

    for (i = 2; i < argc; i++)
    {
        string parameter  = argv[i];

        if (parameter == "-v")
            verbose = 1;
        else if (parameter == "-e")
        {
            if (i == argc-1)
                usage(argv[0]);

            Last_event = atoi(argv[++i]);

            if (Last_event < 1)
            {
                cout << "Invalid last event supplied: " << Last_event << endl;
                usage(argv[0]);
            }
        }
        else if (parameter == "-r")
        {
            if (i == argc-1)
                usage(argv[0]);

            rawdata_ver = atoi(argv[++i]);
            if ((rawdata_ver < 1) || (rawdata_ver > RAW_DATA_FORMATS))
            {
                cout << "Invalid rawdata format supplied: " << rawdata_ver << endl;
                usage(argv[0]);
            }
        }
        else
        {
                cout << "Invalid option: " << parameter << endl;
                usage(argv[0]);
        }
    }

    Run_number = atoi(argv[1]);

    if (!Run_number)
        usage(argv[0]);

    // Print welcome message
    cout << "BOREXINO FADC rawdata to ROOT production program version ";
    cout << READER_VERSION << " started." << endl;
    cout << "Using Echidna cycle " << CYCLE << " data format" << endl;
    cout << "Processing run " << Run_number << " up to ";
    if (Last_event) {
        cout << "event number " << Last_event << "." << endl;
    } else
        cout << "last event." << endl;

    char run_string[10];
    sprintf(run_string, "%05d", Run_number);

    string logname = string(ROOT_FILES_PATH)+"/Run0";
    logname += run_string;
    logname += "_fadc_c";
    logname += Form("%d",CYCLE);
    logname += ".log";

    ofstream log(logname.c_str(), ios::out);
    if (!log)  {
        cerr << "Failed opening LOG file" << endl;
        exit(1);
    }

    string file_out = string(ROOT_FILES_PATH)+"/Run0";
    file_out += run_string;
    file_out += "_fadc_c";
    file_out += Form("%d",CYCLE);
    file_out += ".root";

    TFile* outputfile = new TFile(file_out.c_str(),"RECREATE","BX FWFD rootfile");

    log << "Rootfile: " << file_out << endl;
    if (Last_event)
        log << "Processing up to event " << Last_event << endl;

    TTree* tree = new TTree("BxFwfd","BxFwfd");

    BxFwfd* ev = new BxFwfd();

    Int_t bsize = 32000;
    Int_t split = 99;

    tree->Branch("bxfwfd","BxFwfd",&ev,bsize,split);
    ev->SetRun(Run_number);

    int evnum_prev = -1;
    int jump = 0;
    unsigned int n_neutrino_trg = 0;
    unsigned int keyword1_error = 0;
    unsigned int keyword2_error = 0;
    unsigned int tspectrum_warn = 0;

    int64_t first_event_time = 0, last_event_time = 0;

    string file_in = string(RAWDATA_PATH)+"/bx0";
    file_in += run_string;
    file_in += ".iop.gz";

    gzFile gzFin;

    std::string::size_type pos = file_in.find ("://");
    std::string method = file_in.substr (0, pos);
    std::string tmp_file = "";

    if (method == "http") {
      tmp_file = "/tmp/bx0";
      tmp_file += run_string;
      tmp_file += ".iop.gz";
      std::string command = "wget -nv -O " + tmp_file + " "  + file_in;
      log << "executing " << command << endl;
      system(command.c_str());
      gzFin = gzopen64 ( tmp_file.c_str(), "rb");
    } else gzFin = gzopen64(file_in.c_str(),"rb");
    
    if (gzFin == NULL) {
        cerr << "Rawdata file does not exist: " << file_in << endl;
        outputfile->Close();
        exit(1);
    }
    else  log << "FWFD rawdata file: " << file_in << endl;
    
    gzrewind(gzFin);

    TStopwatch* time = new TStopwatch();
    time->Start();

    string online_ver;
    int trg_ver = 0;

    // suppress ROOT warnings in non-verbose mode
    if (!verbose)
        gErrorIgnoreLevel = kError;

    if (!rawdata_ver)
    {
        // if rawdata format was not supplied - determine it by run number
        if (Run_number < RUN_FIRST_FADC)
        {
            cout << "Trying to guess rawdata format: ";
            rawdata_ver = guess_raw(gzFin);
            if (rawdata_ver > 0)
                cout << rawdata_ver << endl;
            else
            {
                cout << "invalid" << endl;
                cerr << "Rawdata file is corrupted or in unknown format. Exiting" << endl;
                log  << "Rawdata file is corrupted or in unknown format. Exiting" << endl;
                exit(-1);
            }
        }
        else if (Run_number < RUN_FIRST_RAW2) rawdata_ver = 1;
        else if (Run_number < RUN_FIRST_RAW3) rawdata_ver = 2;
        else if (Run_number < RUN_FIRST_RAW4) rawdata_ver = 3;
        else rawdata_ver = 4;
        
    }


    TRawHeader raw_header;

    if (rawdata_ver >= 4)
    {
        // read header (starting with rawdata version 4)

        if (!gzread(gzFin, &raw_header, sizeof(raw_header)))
        {
            cerr << "Unable to read rawdata header" << endl;
            outputfile->Close();
            exit(1);
        }

        if (convert(raw_header.Magic) != RAWDATA_MAGIC)
        {
            cerr << "Rawdata header magic error: got " << hex << convert(raw_header.Magic);
            cerr << " expected " << RAWDATA_MAGIC << endl;
            outputfile->Close();
            exit(1);
        }

        // after ver. 4 rawdata version equals to major version of online
        rawdata_ver = (convert(raw_header.Online_ver) & 0xFFFF0000) >> 16;

        stringstream online_ver_s;
        online_ver_s << hex << ((convert(raw_header.Online_ver) & 0xFFFF0000) >> 16)
                     << "." <<  (convert(raw_header.Online_ver) & 0x0000FFFF);
        online_ver = online_ver_s.str();
        
        // after online 4.3 header contains trigger firmware revision
        trg_ver = convert(raw_header.Trg_FW_ver);
    }
    else
    {
        stringstream online_ver_s;
        online_ver_s << rawdata_ver;
        online_ver = online_ver_s.str();
    }

    cout << "Reading rawdata produced by FADC DAQ ver. " << online_ver
        << " with trigger FW rev. " << hex << trg_ver << dec << endl;
    log  << "Reading rawdata produced by FADC DAQ ver. " << online_ver
        << "; Trigger FW rev. " << hex << trg_ver << dec << endl;

    while (read_raw_event(gzFin, &raw, rawdata_ver) > 0)  {

        ev->SetStatus(true);

        last_event_time = convert(raw.sys.UnixTime);
        ev->SetUnixTime(last_event_time);

        prev_evnum = evnum;
        evnum = convert(raw.sys.Counter); // evnum is now the counter

        if (evnum == 1)
            first_event_time = last_event_time;

        // check evnum sanity
        if (evnum != prev_evnum + 1)
        {
            cerr << "Erroneous evnum: " << evnum << ". Rawdata is corrupted"\
                 " or in wrong format - exiting." << endl;
            log  << "Erroneous evnum: " << evnum << ". Rawdata is corrupted"\
                 " or in wrong format - exiting." << endl;
            outputfile->Close();
            exit(-1);
        }

        if (evnum && !(evnum % 1000))
        {
            cout << "Processed " << evnum << " events." << endl;
            log  << "Processed " << evnum << " events." << endl;
        }

        ev->SetEvNum(evnum);

        int overfl_err = convert(raw.sys.Error);
        ev->SetError(overfl_err);
        
        if (overfl_err)
        {
            if (verbose)
                cout << "ADC buffer overflow error in event " << evnum \
                    << ", after it some events are lost" << endl;
                    
            log << "ADC buffer overflow error in event " << evnum \
                << ", after it some events are lost" << endl;
                
            overfl_errors++;
        }   

        qq  = raw.test1;
        res = convert(qq);

        if (res != KEY2)  {
            log << "FWFD event " << ev->GetEvNum() << ": wrong keyword " << hex << res;
            log << " instead of " << KEY2 << dec << endl;
            if (verbose) {
                cerr << "FWFD event " << ev->GetEvNum() << ": wrong keyword " << hex << res;
                cerr << " instead of " << KEY2 << dec << endl;
            }
            keyword2_error++;
        }

        unsigned short dcode;
        res = convert(raw.trg.fifo82);
        dcode = (res >> 16) & 0xFFFF;

        for (i = 0; i < 16; i++)  {
            if (((dcode >> i) & 1) == 1)  ev->SetDCode(i, 1);
            else  ev->SetDCode(i, 0);
        }

        int trgtype = (convert(raw.trg.fifo84) >> 16) & 0xFFFF;

        if (trgtype == 1)  n_neutrino_trg++;

        if ((trgtype & 9) == 9) susp_checks++; // channel check will be performed on this event

        // remove trigger 16 (CNGS-ID)
        if (trgtype & 16)
            trgtype &= ~16;

        // fix trigger 0
        if (trgtype == 0)
            trgtype = 1;
            
        ev->SetTrgType(trgtype);

        // Waveform decode
        
        wform_based_sum.assign(NTICK,0.);
        int wnd_shift = 0;  // first ADC count belonging to this window

        for (m = 0; m < NMOD; m++)  { //34

            for (i = 0; i < SIZE; i++)  { //64
                qq = raw.adc[m].HalfChannel[0][i];
                wf[m][0][8*i + 2] = (qq & 0xFF000000) >> 24; // A1
                wf[m][0][8*i + 6] = (qq & 0x0000FF00) >>  8; // A2
                wf[m][0][8*i + 0] = (qq & 0x00FF0000) >> 16; // C1
                wf[m][0][8*i + 4] = (qq & 0x000000FF)      ; // C2
                qq = raw.adc[m].HalfChannel[1][i];
                wf[m][0][8*i + 3] = (qq & 0xFF000000) >> 24; // B1
                wf[m][0][8*i + 7] = (qq & 0x0000FF00) >>  8; // B2
                wf[m][0][8*i + 1] = (qq & 0x00FF0000) >> 16; // D1
                wf[m][0][8*i + 5] = (qq & 0x000000FF)      ; // D2

                qq = raw.adc[m].HalfChannel[2][i];
                wf[m][1][8*i + 2] = (qq & 0xFF000000) >> 24; // A1
                wf[m][1][8*i + 6] = (qq & 0x0000FF00) >>  8; // A2
                wf[m][1][8*i + 0] = (qq & 0x00FF0000) >> 16; // C1
                wf[m][1][8*i + 4] = (qq & 0x000000FF)      ; // C2
                qq = raw.adc[m].HalfChannel[3][i];
                wf[m][1][8*i + 3] = (qq & 0xFF000000) >> 24; // B1
                wf[m][1][8*i + 7] = (qq & 0x0000FF00) >>  8; // B2
                wf[m][1][8*i + 1] = (qq & 0x00FF0000) >> 16; // D1
                wf[m][1][8*i + 5] = (qq & 0x000000FF)      ; // D2

                qq = raw.adc[m].HalfChannel[4][i];
                wf[m][2][8*i + 2] = (qq & 0xFF000000) >> 24; // A1
                wf[m][2][8*i + 6] = (qq & 0x0000FF00) >>  8; // A2
                wf[m][2][8*i + 0] = (qq & 0x00FF0000) >> 16; // C1
                wf[m][2][8*i + 4] = (qq & 0x000000FF)      ; // C2
                qq = raw.adc[m].HalfChannel[5][i];
                wf[m][2][8*i + 3] = (qq & 0xFF000000) >> 24; // B1
                wf[m][2][8*i + 7] = (qq & 0x0000FF00) >>  8; // B2
                wf[m][2][8*i + 1] = (qq & 0x00FF0000) >> 16; // D1
                wf[m][2][8*i + 5] = (qq & 0x000000FF)      ; // D2
            }

            // Decode timing information
            res = convert(raw.adc[m].TimeLog);
            int stopaddr = ((res & 0x3F) << 3) + 8;
            
            if (m == 0)  {

                time_prev = (res & 0x7FFFFFFF) * 20; // ns. Will be also kept in clusters, corrected to cluster's peak_pos
                ev->SetRawTime(time_prev);
                
                if (time_prev < WND_LEN_NS)
                    wnd_shift = (int)((WND_LEN_NS - time_prev) / NS_PER_TICK) - 8;
                
                res = convert(raw.adc[m].Pattern);

                // get LABEN event number
                int evnum_bx;
                if (Run_number >= RUN_NEW_BTB_FIRMWARE)
                    evnum_bx = res & 0x7FFF;
                else
                    evnum_bx = res & 0x0000FFFF;    // OLD BTB firmware (before run 12422)

                // After new LABEN firmware has been installed (run 12422)
                // evnum_bx can jump forward on 512, but next event come back.
                // That's why I introduce additional check for, e.g., 2048 (btw even 1024 should be enough):
                if ( evnum_bx < evnum_prev && (evnum_prev - evnum_bx) > 2048 )  jump++;

                evnum_prev = evnum_bx;
                evnum_bx += jump*32768;
                ev->SetEvNumBx(evnum_bx);

            }
            
            // Check per ADC keyword
            qq  = raw.adc[m].test2;
            res = convert(qq);
            if (res != KEY1)  {
                log << "FWFD event " << ev->GetEvNum() << " module " << m << ": wrong keyword ";
                log << hex << res << " instead of " << KEY1 << dec << endl;
                
                if (verbose) {
                    cerr << "FWFD event " << ev->GetEvNum() << " module " << m << ": wrong keyword ";
                    cerr << hex << res << " instead of " << KEY1 << dec << endl;
                }
                
                keyword1_error++;
            }

            for (k = 0; k < 3; k++)  { // cycle on each channel within the m-th module

                current_ChID = m*3 + k; // 0 - 101
                base.assign(NBASE,0.);
                wform.assign(NTICK,0);
                wform_based.assign(NTICK,0.);

                for (i = wnd_shift; i < NTICK; i++)  
                    wform[i-wnd_shift] = wf[m][k][(i + stopaddr) % NTICK];

                //__________________________________
                // Check that all channels digitize smth relevant

                if ((trgtype & 9) == 9) {
                    // check only muon and pulser events
                    int wform_max = *max_element(wform.begin(), wform.end());
                    int wform_min = *min_element(wform.begin(), wform.end());

                    if ((wform_max - wform_min) < CHANNEL_CHECK_MIN_DIF)
                        susp_channels[current_ChID]++;
                }

                //__________________________________

                // Calculate pedestals (8 (4) per each channel within first 64 (32) ticks)

                if (m == 0 && k == 0)  {
                    vector<unsigned int>::iterator it0 = min_element(wform.begin(), wform.begin()+128);
                    peak_pos0 = it0 - wform.begin();
                }

                if (peak_pos0 > 64)  { // in this case the calculation is 8-based
                    for (i = 0; i < NBASE; i++)  {
                        for (j = 0; j < NBASE; j++)  {
                            base[j] += (float)wform[i*8 + j] / 8.;
                        }
                    }

                    // Pedestals subtraction
                    for (i = 0; i < SIZE; i++)  { //64
                        for (j = 0; j < NBASE; j++)  {
                            wform_based[i*8 + j] = base[j] - (float)wform[i*8 + j];
                        }
                    }
                }
                else  {  // in this case the calculation is 4-based
                    for (i = 0; i < 4; i++)  {
                        for (j = 0; j < 4; j++)  {
                            base[j] += (float)wform[i*4 + j] / 4.;
                        }
                    }

                    // Pedestals subtraction
                    for (i = 0; i < 128; i++)  { //64 * 2! Else the last 256 ticks in the window will be lost
                        for (j = 0; j < 4; j++)  {
                            wform_based[i*4 + j] = base[j] - (float)wform[i*4 + j];
                        }
                    }
                }
                
                // Fill "fake" part of window with zeroes
                for (i = NTICK - wnd_shift; i < NTICK; i++)
                    wform_based[i] = 0;

                //_______________________________________

                // Clustering! There should always be at least one cluster
                vector<float>::iterator it1 = max_element(wform_based.begin(), wform_based.end());
                unsigned int peak_posi = it1 - wform_based.begin();
                peak_amp = *it1;

                if (m == 0 && k == 0)  { // zero channel - ASUM

                    for (i = 0; i < NTICK; i++)  source[i] = wform_based[i];

                    peak_pos[0] = peak_posi;
                    peak_ampl[0] = peak_amp;

                    for (i = 0; i < 6; i++)  Q[i] = 0.; // 6 clusters at max.

                    nclu = 0;

                    // Save analog waveform
                    ev->ClearWForms(); 
                    ev->SetWFormAsum(wform_based);

                    if (peak_amp > 25)  { // Call TSpectrum only in case of "high-energy" events with peak_ampl > 25

                        TSpectrum s(24);
                        nclu = s.SearchHighRes(source, dest, 512, 2.5, 33., kFALSE, 1, kFALSE, 1);

                        if (nclu == 0)  nclu = s.SearchHighRes(source, dest, 512, 2.5, 33., kTRUE, 1, kFALSE, 1);

                        if (nclu > 6)  {
                            log << "FWFD event " << evnum << ": Warning in <TSpectrum::SearchHighRes>: Peak buffer full" << endl;
                            tspectrum_warn++;
                            nclu = 6;
                        }

                        Float_t* xpeaks = s.GetPositionX();

                        for (int p = 0; p < nclu; p++)  {

                            peak_pos[p] = int(xpeaks[p] + 0.5);
                            int pos = peak_pos[p];
                            peak_ampl[p] = wform_based[pos];

                            if (peak_pos[p] < 15)          cluster_start[p] = 0;
                            else                           cluster_start[p] = peak_pos[p] - 15;
                            if ((peak_pos[p] + 40) > 511)  cluster_end[p] = 511;
                            else                           cluster_end[p] = peak_pos[p] + 40;

                            Q[p] = 0.;
                            for (int k1 = cluster_start[p]; k1 < cluster_end[p]; k1++)  Q[p] += wform_based[k1];
                        }
                        
                        // if TSpectrum found zero clusters we have to create a fa

                        //____ Check if clusters are piled-up? ____

                        if (nclu > 1)  {

                            // Sorting clusters chronologically
                            bool t = true;
                            while (t) {
                                t = false;
                                for (i = 0; i < nclu-1; i++)  {
                                    if (peak_pos[i] > peak_pos[i+1])  {
                                        int temp_pos = peak_pos[i+1];
                                        peak_pos[i+1] = peak_pos[i];
                                        peak_pos[i] = temp_pos;

                                        float temp_ampl = peak_ampl[i+1];
                                        peak_ampl[i+1] = peak_ampl[i];
                                        peak_ampl[i] = temp_ampl;

                                        float temp_Q = Q[i+1];
                                        Q[i+1] = Q[i];
                                        Q[i] = temp_Q;
                                        t = true;
                                    }
                                }
                            }

                            for (i = 0; i < 6; i++)  pile_up[i] = false;

                            for (i = 1; i < nclu; i++)  {
                                if (peak_pos[i] - peak_pos[i-1] < 50)  {

                                    Q[i] = 0.; Q[i-1] = 0.;

                                    for (unsigned int bb = peak_pos[i-1] - 15; bb < peak_pos[i] + 40; bb++)  Q[i] += wform_based[bb];

                                    Q[i-1] = Q[i]; // Q[i-1] is currently the total charge in two clusters
                                    Q[i-1] = peak_ampl[i-1] / 0.06;
                                    Q[i] -= Q[i-1];

                                    pile_up[i-1] = true;
                                    pile_up[i] = true;
                                }
                            }
                        }
                    }  // end calling TSpectrum
                    
                    // if we could not use TSpectrum - search for cluster manually
                    if (!nclu)  { // only one cluster is possible

                        nclu = 1;
                        if (peak_pos[0] < 15)          cluster_start[0] = 0;
                        else                           cluster_start[0] = peak_pos[0] - 15;
                        if ((peak_pos[0] + 40) > 511)  cluster_end[0] = 511;
                        else                           cluster_end[0] = peak_pos[0] + 40;

                        for (int k2 = cluster_start[0]; k2 < cluster_end[0]; k2++)  Q[0] += wform_based[k2];
                    }
                } // end analyzing zero-channel - ASUM

                //______________

                if (m == 0 && k == 1)  { // channel 1 - attenuated (factor 3.2) ASUM

                    for (i = 0; i < 6; i++)  atten_Q[i] = 0.;  // 6 clusters at max.

                    for (i = 0; i < nclu; i++)  {
                        if (pile_up[i])  atten_Q[i] = Q[i]/3.2;
                        else  {

                        if (peak_pos[i] < 15)          cluster_start[i] = 0;
                        else                           cluster_start[i] = peak_pos[i] - 15;
                        if ((peak_pos[i] + 40) > 511)  cluster_end[i] = 511;
                        else                           cluster_end[i] = peak_pos[i] + 40;

                        for (int kkk = cluster_start[i]; kkk < cluster_end[i]; kkk++)  atten_Q[i] += wform_based[kkk];
                        }
                    }
                }

                //_______________

                if (m > 0)  {

                    // Starting from run 18016 exclude this channels from the DSum, because
                    // they are used now for BX-CNGS project:
                    if (Run_number > RUN_CNGS_START) {
                        if (current_ChID == 66)  continue;
                        if (current_ChID == 67)  continue;
                        if (current_ChID == 69)  continue;
                        if (current_ChID == 70)  continue;
                    }

                    for (j = 0; j < NTICK; j++)  wform_based_sum[j] += wform_based[j];

                }

            }   // end of cycle on channels (k)
        }       // end of cycle on modules (m)

        //_______________ End of event ________________

        // Deconvoluted digital sum

        dd.assign(NTICK,0.);

        float bmp = 0.;
        for (int d = 0; d < NTICK; d++)  {
            bmp += wform_based_sum[d];
            dd[d] = wform_based_sum[d] + bmp/TAU;
        }

        ev->SetWFormDsum(dd);

        for (i = 0; i < 6; i++)  QD[i] = 0.; // 6 clusters at max.

        ev->ClearClusters();

        for (unsigned int cl = 0; cl < nclu; cl++)  {

            if (peak_pos[cl] > 13)  peak_posD[cl] = peak_pos[cl] - 13; // 13 bins shifted averagely
            else peak_posD[cl] = 0;

            if (peak_posD[cl] < 15)          cluster_startD = 0;
            else                             cluster_startD = peak_posD[cl] - 15;
            if ((peak_posD[cl] + 40) > 511)  cluster_endD = 511;
            else                             cluster_endD = peak_posD[cl] + 40;

            for (int hh = cluster_startD; hh < cluster_endD; hh++)  QD[cl] += dd[hh];

            QD[cl] = QD[cl]/10.; // charge of the deconvoluted DSum

            // How can I calculate digital charge in case of piled-up events ??
            if (pile_up[cl] && peak_ampl[cl] < 200)  QD[cl] = Q[cl] * 0.92;

            #if CYCLE >= 16
            if (cl == 0)
                ev->SetCluster(nclu, cl, peak_pos[cl], time_prev*1.0E-6 - 0.0000025*(511-peak_pos[cl]) + time_rest, peak_ampl[cl], Q[cl], atten_Q[cl], QD[cl], 0, 0, 0);
            else
                ev->SetCluster(nclu, cl, peak_pos[cl], 0.0000025*(peak_pos[cl] - peak_pos[cl-1]), peak_ampl[cl], Q[cl], atten_Q[cl], QD[cl], 0, 0, 0); // time_prev in ms
            #else
            if (cl == 0)
                ev->SetCluster(nclu, cl, peak_pos[cl], time_prev*1.0E-6 - 0.0000025*(511-peak_pos[cl]) + time_rest, peak_ampl[cl], Q[cl], atten_Q[cl], QD[cl]);
            else
                ev->SetCluster(nclu, cl, peak_pos[cl], 0.0000025*(peak_pos[cl] - peak_pos[cl-1]), peak_ampl[cl], Q[cl], atten_Q[cl], QD[cl]); // time_prev in ms
            #endif

            if (cl+1 == nclu)  time_rest = (511 - peak_pos[cl]) * 0.0000025; // ms

        }

        // decode GPS time stamp only if there were no GPS errors in previous events
        if (gps_errors == 0)
            gps_errors += gps_decode(verbose, rawdata_ver, raw_header, raw, ev);
        else
            gps_errors += gps_decode(verbose, 1, raw_header, raw, ev);  // will not attempt to decode GPS

        ev->SetNClusters(nclu);
        tree->Fill();

        // stop on requested event number
        if (Last_event && evnum == Last_event)  break;

    } // End of while

    tree->Write();
    outputfile->Close();

    delete ev;

    cout << evnum << " events written" << endl;
    log  << evnum << " events written" << endl;

    if (evnum < Last_event) {
        // number of events in data was smaller than supplied last event
        cout << "Warning: you requested to process " << Last_event;
        cout << " events, but there are only " << evnum;
        cout << " events in raw data file for this run." << endl;
    }

    unsigned bad_channels = 0;

    if (evnum > CHANNEL_CHECK_MIN_EVNUM) {
        // Report channels failed check
        for (i = 0; i < NCH; i++)
            if (susp_channels[i] > (susp_checks * CHANNEL_CHECK_MISS_RATE)) {
    
                bool ignore_ch = false;
    
                for (j = 0; j < NOCHECK_CH_CNT; j++)
                    if (i == nocheck_channels[j]) {
                        // do not report on ignored channels
                        ignore_ch = true;
                        break;
                    }
    
                if (!ignore_ch) {
                    // this channel is bad - report it
                    cout << "Channel " << i << " data check failed! Cabling problem?" << endl;
                    log  << "Channel " << i << " data check failed! Cabling problem?" << endl;
                    bad_channels++;
                }
            }
    
        if (!bad_channels) {
            cout << "All channels were working fine." << endl;
            log  << "All channels were working fine." << endl;
        }
    }

    time->Stop();

    #if WRITE_GPS_SANE_LIST
    // if GPS timestamps are fine - add run to GPS sane list
    if (!gps_errors)
    {
        int addtolist = add_to_gps_sane_list((string)run_string);

        if (addtolist < 0)
            cout << "Unable to open GPS sanity list files." << endl;
        else if (!addtolist)
            cout << "GPS times are sane, run added to GPS sane list file." << endl;
        else if (addtolist == 1)
            cout << "GPS times are sane, run is already present in GPS sane list file." << endl;
        else if (addtolist == 2)
            cout << "GPS times seem sane, but run is present in GPS insane list file." << endl;
    }
    #endif

    // make readable time representations
    char start_time[100], end_time[100];
    strftime(start_time, 50, "%d.%m.%y %H:%M:%S", localtime((time_t*)&first_event_time));
    strftime(end_time, 50, "%d.%m.%y %H:%M:%S", localtime((time_t*)&last_event_time));

    // Output report to stdout
    cout << "RUN duration: " << last_event_time - first_event_time << " s" << endl;
    cout << "First event time: " << start_time << " Last event time: " << end_time << endl;

    cout << n_neutrino_trg << " NEUTRINO triggers (" << n_neutrino_trg*100/evnum << "%)" << endl;
    cout << tspectrum_warn << " TSpectrum warnings (" << tspectrum_warn*100/evnum << "%)" << endl;
    cout << keyword1_error << " keyword1 errors" << endl;
    cout << keyword2_error << " keyword2 errors" << endl;
    if (evnum > CHANNEL_CHECK_MIN_EVNUM)
        cout << bad_channels   << " channel errors"  << endl;
    cout << gps_errors << " GPS errors " << endl;
    cout << overfl_errors << " overflow errors" << endl;

    cout << "time cpu: " << time->CpuTime() << " s" << endl;
    cout << "time real: " << time->RealTime() << " s" << endl;

    // Output report to log
    log << "RUN duration: " << last_event_time - first_event_time << " s" << endl;

    log << n_neutrino_trg << " NEUTRINO triggers (" << n_neutrino_trg*100/evnum << "%)" << endl;
    log << tspectrum_warn << " TSpectrum warnings (" << tspectrum_warn*100/evnum << "%)" << endl;
    log << keyword1_error << " keyword1 errors" << endl;
    log << keyword2_error << " keyword2 errors" << endl;
    log << bad_channels   << " channel errors"  << endl;
    log << gps_errors << " GPS errors " << endl;
    log << overfl_errors << " overflow errors" << endl;

    log << "time cpu: " << time->CpuTime() << " s" << endl;
    log << "time real: " << time->RealTime() << " s" << endl;

    log.close();
    
    std::string command = "rm ";
    command += tmp_file;
    system(command.c_str()); 

    delete time;
    return 0;
}
