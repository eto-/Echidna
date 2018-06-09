/*
 ***************** DETERMINE MUONS BY FADC PULSE SHAPE *****************
 */

int is_muon(BxFwfd* ev, TTree* tree);

const int       START_CHANNEL           = 0;
const int       FADC_EVENT_SIZE         = 512;
const int       DISC_CHANNELS           = 16;

const double    ACHARGE_TO_MEV          = 378.;
const double    DCHARGE_TO_MEV          = 352.;

const double    NOISE_TIME_AFTER_MUON   = 0.02; // 20 us

const int       LOW_AMPLITUDE           = 25;
const int       MUON_CLOSE_PEAKS        = 30;

const int       HALF_AMPL_GAP           = 30;

const int       MAX_NOISE_AMPL          = 80;

const double    MIN_ASUM_MEV            = 1.0;
const double    MIN_DSUM_MEV            = 1.0;

const UInt_t    MIN_TIME_PREV           = 1300;

const int       DSUM_INTERVAL           = 5;
const int       MUON_RISE_TIME_DSUM     = 24;
const int       MUON_MAX_RISE_TIME_DSUM = 600;
const int       MUON_HIGH_AMPLITUDE_DSUM= 18000;
const int       NOISE_THRESHOLD_DSUM    = 30;
const int       MUON_CLOSE_PEAKS_DSUM   = 30;

const int       DSUM_2PEAK              = 30;
const int       DSUM_MIN2PEAK           = 8;
const double    DSUM_PEAKDIF            = 0.70; // in peak ampl
const int       DSUM_2PEAK_MINAMPL      = 550;


/* Reasons */
const int       MUON_REASONS            = 5;
const int       REJECT_REASONS          = 7;

const int       REASON_NONMUON_NOTHING  =  0;
const int       REASON_NONMUON_LOWAMPL  = -1;
const int       REASON_NONMUON_SMALLTTP = -2;
const int       REASON_NONMUON_NOISE    = -3;
const int       REASON_NONMUON_PULSER   = -4;
const int       REASON_NONMUON_AFTMUON  = -5;
const int       REASON_NONMUON_DCODE    = -6;

const int       REASON_MUON_MTB         = 1;
const int       REASON_MUON_HIAMPL      = 2;
const int       REASON_MUON_TWOPEAKS    = 3;
const int       REASON_MUON_SLOWRISE    = 4;
const int       REASON_MUON_CNGS        = 5;


/* Get time from last muon event in ms */
double time_from_prev_muon(TTree* tree, BxFwfd* ev)
{
    int evnum = ev->GetEvNum()-1;
    double result = 0.;
    
    for (int i = evnum; ; i--)
    {
        // stop on first event
        if (i == 0)
        {
            result += 10000.;
            break;
        }

        // use raw time    
        result += ev->GetRawTime() / 1000000.;
        
        tree->GetEntry(i-1);
        
        // previous muon found
        if (is_muon(ev, tree))
            break;
    }
    
    tree->GetEntry(evnum);
    return result;
}


/* Determine noise FADC events */
int is_noise_mu(BxFwfd* ev)
{   
    if (!ev->GetNClusters())
        return 1;
        
    if (ev->GetTrgType() == 65)
        return 4;
    
    Float_t first_ampl = ev->GetCluster(0).GetPeakAmpl();
    int first_peak_pos = ev->GetCluster(0).GetPeakPos();
    
    if (   (first_ampl < MAX_NOISE_AMPL) 
        && ((first_peak_pos < HALF_AMPL_GAP) || ev->GetWFormAsum(first_peak_pos - HALF_AMPL_GAP) >  (first_ampl / 2)))
        return 2;
        
    if ((ev->GetNClusters() > 3) && (first_ampl < MAX_NOISE_AMPL/2) && (ev->GetCluster(1).GetPeakAmpl() < MAX_NOISE_AMPL/2))
        return 3;

    if ((ev->GetNClusters() > 3) && (first_ampl < MAX_NOISE_AMPL/2) && (ev->GetCluster(1).GetPeakAmpl() < MAX_NOISE_AMPL/2))
      return 5; // 

    return 0;    
}

/* Determine "muonish" FADC events */
int is_muon(BxFwfd* ev, TTree* tree)
{   
    int i;
    int trgtype = ev->GetTrgType();
    int reason = 0;
    
    if ((trgtype & 2) || (trgtype & 4))
        return REASON_NONMUON_PULSER; /* this is pulser or laser event */
        
    // count active discriminators channels
    int d_chan = 0;
    for (i = 0; i < DISC_CHANNELS; ++i)
        if (ev->GetDCode(i))
            ++d_chan;
            
    if (d_chan < 5)
        return REASON_NONMUON_DCODE;
        
    if (((trgtype & 15) == 9)) 
    {
        return REASON_MUON_MTB;     // FADC muon trigger from MTB
    }
        
    if (((trgtype & 33) == 33)) 
    {
        return REASON_MUON_CNGS;     // FADC muon trigger from CNGS OD-discr
    }
        
    if (is_noise_mu(ev))
        return REASON_NONMUON_NOISE; /* this is junk event, not a muon */    
        
    Float_t first_ampl = ev->GetCluster(0).GetPeakAmpl();
    int first_peak_pos = ev->GetCluster(0).GetPeakPos();
    
    if (first_ampl < LOW_AMPLITUDE)
        return REASON_NONMUON_LOWAMPL; /* too low amplitude to determine... */
    
    if (ev->GetRawTime() < MIN_TIME_PREV)
        return REASON_NONMUON_SMALLTTP; /* small time to previous - probably junk */
               

    // Manually find the peak position in DSUM waveform
    Float_t first_dsum_ampl = 0;
    int first_dsum_peak_pos = 0;
    Float_t second_dsum_ampl = 0;
    int second_dsum_peak_pos = FADC_EVENT_SIZE*2;  
    
    for (i = (/*(first_peak_pos > START_CHANNEL+DSUM_INTERVAL) ? (first_peak_pos - DSUM_INTERVAL) :*/ START_CHANNEL); 
         i < ((first_peak_pos < FADC_EVENT_SIZE-DSUM_INTERVAL) ? (first_peak_pos + DSUM_INTERVAL) : FADC_EVENT_SIZE);
         i++)
        if (ev->GetWFormDsum(i) > first_dsum_ampl)
        {
            first_dsum_ampl = ev->GetWFormDsum(i);
            first_dsum_peak_pos = i;
        }
        
    if (first_dsum_ampl > MUON_HIGH_AMPLITUDE_DSUM)
    {
        return REASON_MUON_HIAMPL; /* looks like really large muon signal */
    }

    for (i = first_dsum_peak_pos; i > START_CHANNEL; i--)
        if (ev->GetWFormDsum(i) <= NOISE_THRESHOLD_DSUM)
        {
            int risetime = (first_dsum_peak_pos - i);
            
            if ((risetime > MUON_RISE_TIME_DSUM) && (risetime < MUON_MAX_RISE_TIME_DSUM))
            {
                reason = REASON_MUON_SLOWRISE; // looks like muon with slow rise front, have to check prev muon
            }
            
            break;
        } 

    if (!reason && (first_dsum_ampl > DSUM_2PEAK_MINAMPL))
    {    
        // Find second DSUM peak if it is present
        for (i = ((first_dsum_peak_pos > START_CHANNEL+DSUM_2PEAK) ? (first_dsum_peak_pos - DSUM_2PEAK) : START_CHANNEL); 
             i < ((first_dsum_peak_pos < FADC_EVENT_SIZE-DSUM_2PEAK) ? (first_dsum_peak_pos + DSUM_2PEAK) : FADC_EVENT_SIZE);
             i++)
            if (   (ev->GetWFormDsum(i) > (first_dsum_ampl * DSUM_PEAKDIF))
                && (abs(i - first_dsum_peak_pos) > DSUM_MIN2PEAK)
                && (ev->GetWFormDsum(i) > second_dsum_ampl)
                && (i > START_CHANNEL + 3)
                && (i < FADC_EVENT_SIZE - 3)
                && (ev->GetWFormDsum(i) > ev->GetWFormDsum(i+3))
                && (ev->GetWFormDsum(i) > ev->GetWFormDsum(i-3)))
            {
                second_dsum_ampl = ev->GetWFormDsum(i);
                second_dsum_peak_pos = i;
            }
            
        if ((abs(first_dsum_peak_pos - second_dsum_peak_pos) < MUON_CLOSE_PEAKS_DSUM) && !(ev->GetTrgType() & 64))
        {
            reason = REASON_MUON_TWOPEAKS; /* two close peaks in DSUM - looks like muon signal, have to check prev muon */
        }
    } 
    #if 0   // do not use 2-peak method for low-amplitude events...
    else if (!reason && (ev->GetNClusters() > 1))
    {   
        if (   ((ev->GetCluster(1).GetPeakPos() - first_peak_pos) < MUON_CLOSE_PEAKS)
            && ((ev->GetCluster(1).GetASumCharge() / ACHARGE_TO_MEV) > MIN_ASUM_MEV)
            && ((ev->GetCluster(1).GetDSumCharge() / DCHARGE_TO_MEV) > MIN_DSUM_MEV)
           )
        {
            reason = REASON_MUON_TWOPEAKS; /* two close peaks - looks like muon signal, have to check prev muon */
        }
    }
    #endif
    
    // if we've found muon by pulseshape - we have to check time from previous muon to reject noise
    if (reason)
    {
        if (time_from_prev_muon(tree, ev) > NOISE_TIME_AFTER_MUON)
            return reason;
        else
            return REASON_NONMUON_AFTMUON;
    }
                    
    return 0; /* overwise this is not a muon */

}    
