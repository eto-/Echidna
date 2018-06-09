/*
 ***************** DETERMINE NOISE BY FADC PULSE SHAPE *****************
 *
 * $Id:
 *
 */

//the peak, last and first are defined as integrals of 5 bins
//integral is done on DSum, the peak position is shifted by 13 with respect to cluster found on ASUM
//offset of 10 of first is increased to 13

//int is_noiseD(BxFwfd* ev, int cluster);

/* Reasons */
const int       NOISE_REASONS_D               = 1;
const int       REASON_NOISE_NONSCINTEDGES_D  = 1;


/* Determine noise FADC events */
int is_noiseD(BxFwfd* ev, int cluster)
{
  int peak_pos = ev->GetCluster(cluster).GetPeakPos() - 13;
  float peak_ampl = ev->GetWFormDsum(peak_pos);
  float last = 0;
  float first = 0;
  
//  std::cout << "peak position " << peak_pos << std::endl;
//  std::cout << "peak amplitude in D " << peak_ampl << std::endl;

  //integral of 5
  for(int i = 1; i < 3; i++) peak_ampl += (ev->GetWFormDsum(peak_pos + i) + ev->GetWFormDsum(peak_pos - i));
  peak_ampl = peak_ampl/5;

  if (peak_pos > 9)  {
    double firstsum = 0;
    first = ev->GetWFormDsum(peak_pos-13);
    for(int i = 1; i < 3; i++) first += (ev->GetWFormDsum(peak_pos -13 + i) + ev->GetWFormDsum(peak_pos -13 - i));
    first = first /5;
  }  
  else first = ev->GetWFormDsum(0);
  if (peak_pos < 492)  {
    last = ev->GetWFormDsum(peak_pos+18);//20
    for(int i = 1; i < 3; i++) last += (ev->GetWFormDsum(peak_pos +18 + i) + ev->GetWFormDsum(peak_pos +18 - i));
    last = last /5;
  }
  else last = ev->GetWFormDsum(511);
//  std::cout << "last " << last << std::endl;
//  std::cout << "first " << first << std::endl;
  if (first > peak_ampl/4. && last > peak_ampl/3.)  return 1;

  return 0;
}

