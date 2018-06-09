/*
 ***************** DETERMINE MUONS BY FADC PULSE SHAPE *****************
 */
// the amplitude, last, first is defined as the integral of 5 bins on ASum

int is_noise(BxFwfd* ev);

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

/* Reasons */
const int       NOISE_REASONS               = 1;
const int       REASON_NOISE_NONSCINTEDGES  = 1;


/* Determine noise FADC events */
int is_noise(BxFwfd* ev, int cluster)
{
  int peak_pos = ev->GetCluster(cluster).GetPeakPos();
  float peak_ampl = ev->GetCluster(cluster).GetPeakAmpl();
  float last = 0;
  float first = 0;
  
  std::cout << "peak position " << peak_pos << std::endl;

  //integral of 5
  for(int i = 1; i < 3; i++) peak_ampl += (ev->GetWFormAsum(peak_pos + i) + ev->GetWFormAsum(peak_pos - i));
  peak_ampl = peak_ampl/5;

  if (peak_pos > 9)  {
    double firstsum = 0;
    first = ev->GetWFormAsum(peak_pos-10);
    for(int i = 1; i < 3; i++) first += (ev->GetWFormAsum(peak_pos -10 + i) + ev->GetWFormAsum(peak_pos -10 - i));
    first = first /5;
  }  
  else first = ev->GetWFormAsum(0);
  if (peak_pos < 492)  {
    last = ev->GetWFormAsum(peak_pos+20);
    for(int i = 1; i < 3; i++) last += (ev->GetWFormAsum(peak_pos +20 + i) + ev->GetWFormAsum(peak_pos + 20 - i));
    last = last /5;
  }
  else last = ev->GetWFormAsum(511);
  if (first > peak_ampl/4. && last > peak_ampl/3.)  return 1;

  return 0;
}
