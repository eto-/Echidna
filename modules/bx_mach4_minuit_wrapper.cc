#include "bx_mach4_minuit_wrapper.hh"
#include <TMinuit.h>
#include <cmath>

// Interface to MINUIT

MinuitWrapper::MinuitWrapper(void (* fcn)(int &, double *, double &,
			     double *, int))
{
  mn = new TMinuit;
  mn->SetFCN(fcn);
  
  // initialize MINUIT
  mn->mninit(5, 6, 7);
  
  // tell MINUIT to be quiet
  int errflag;
#ifndef DEBUG_RECON
  mn->mncomd("SET PRINTOUT -1", errflag);
  mn->mncomd("SET NOWARNINGS", errflag);
#endif

  // Set "UP" to 0.5 for log likelihood fit:
  mn->mncomd("SET ERRORDEF 0.5", errflag);
  
  // Set strategy to the most reliable (albeit more time-consuming)
  mn->mncomd("SET STRATEGY 2", errflag);

  std::fill_n(names, MAX_MINUIT_PARS, (std::string *)0);
  std::fill_n(values, MAX_MINUIT_PARS, 0);
  std::fill_n(errors, MAX_MINUIT_PARS, 0);
  std::fill_n(limits[0], MAX_MINUIT_PARS, 0);
  std::fill_n(limits[1], MAX_MINUIT_PARS, 0);
  std::fill_n(is_set, MAX_MINUIT_PARS, false);
}

MinuitWrapper::~MinuitWrapper()
{
  delete mn;
  for (int i = 0; i < MAX_MINUIT_PARS; i++)
    delete names[i];
}

void MinuitWrapper::SetParameter(size_t index, double value)
{
  is_set[index] = true;
  values[index] = value;
}

void MinuitWrapper::SetParameterName(size_t index, std::string name)
{
  is_set[index] = true;
  delete names[index];
  names[index] = new std::string(name);
}

void MinuitWrapper::SetParLimits(size_t index, double min, double max)
{
  is_set[index] = true;
  limits[0][index] = min;
  limits[1][index] = max;
}

double MinuitWrapper::GetParameter(size_t index)
{
  return values[index];
}

double MinuitWrapper::GetParError(size_t index)
{
  return errors[index];
}

void MinuitWrapper::Execute(int & errflag, int & convergence)
{
  // copy saved settings into parameters
  mn->mncler();
  for (size_t i = 0; i < MAX_MINUIT_PARS; i++)
    if (is_set[i]) {
      if (! names[i]) { // name "p0", "p1", etc.
	names[i] = new std::string("p");
	if (i > 9)
	  *(names[i]) += (char)(i / 10 + '0');
	*(names[i]) += (char)(i % 10 + '0');
      }
      
      double stepsize;
      if (limits[0][i] || limits[1][i])
	stepsize = 0.2 * (limits[1][i] - limits[0][i]);
      else if (values[i])
	stepsize = 0.2 * values[i];
      else
	stepsize = 1;
      
      int dummy;
      mn->mnparm(i, *(names[i]), values[i], stepsize,
		 limits[0][i], limits[1][i], dummy);
    }
  
  mn->mncomd("MIGRAD", errflag);
#if 0 
  if (!errflag) {
    // Try to find a better nearby minimum
    mn->mncomd("IMPROVE", errflag);
    // Improve error estimates
    mn->mncomd("HESSE", errflag);
  }
#endif
  TString s;
  double dummyf[3];
  int dummyi, npari, nparx;

  mn->mnstat(dummyf[0], dummyf[1], dummyf[2],
	     npari, nparx, convergence);
  
  for (int i = 0; i < nparx; i++)
    if (is_set[i])
      mn->mnpout(i, s, values[i], errors[i], limits[0][i], limits[1][i],
		 dummyi);
}
