#include "bx_scint_pdf_mach4.hh"
#include <stdint.h>
#include <cmath>

void ScintPDF::Load(const ScintPDF & p)
{
  have_reflection = p.have_reflection;
  q[0] = p.q[0]; q[1] = p.q[1]; q[2] = p.q[2];
  sigma = p.sigma;
  tau[0] = p.tau[0]; tau[1] = p.tau[1]; tau[2] = p.tau[2];
  refl[0] = p.refl[0]; refl[1] = p.refl[1]; refl[2] = p.refl[2];

  a[0] = p.a[0]; a[1] = p.a[1]; a[2] = p.a[2];
  b = p.b;
  c[0] = p.c[0]; c[1] = p.c[1]; c[2] = p.c[2];
  lambda[0] = p.lambda[0]; lambda[1] = p.lambda[1]; lambda[2] = p.lambda[2];
  lambdadiff[0] = p.lambdadiff[0]; lambdadiff[1] = p.lambdadiff[1];
  gcoeff = p.gcoeff;
}

template <typename T>
static T sqr (T x) { return x*x; }

void ScintPDF::CalculateCoefficients()
{
  static const double sqrt2_inv = 1.0 / std::sqrt(2.0);

  b = sqrt2_inv / sigma;
  for (uint32_t i = 0; i < 2; i++) {
    a[i] = q[i] * std::exp(0.5 * sqr(sigma / tau[i]));
    c[i] = -sigma * sqrt2_inv / tau[i];
    lambda[i] = 1.0 / tau[i];
    if (i > 0) lambdadiff[i - 1] = lambda[0] - lambda[i];
  }
  if (have_reflection) gcoeff = -0.5 / sqr(refl[2]);
  // careful to set have_reflection before calling this function
}

double ScintPDF::operator() (double t) const
{
  double result = 0;
  for (unsigned i = 0; i < 2; i++) 
    {
      double argerf = b * t + c[i];
      result += 0.5 * a[i] * erfc(-argerf) * std::exp(-t * lambda[i]) * lambda[i];
    }
  if (have_reflection) {
    result += refl[0]*std::exp(gcoeff*sqr(t - refl[1]));
  }
  return result;
}

double ScintPDF::Log(double t) const
{
  // Kevin's idea about reusing exponential worked nice, but with the new shape
  // we get larger values for lambdadiff, and std::exp( ) blows up into inf
  // so we stick to the pedestrian method for now. (ps: and it doesn't work
  // if we add the gaussian for the reflected light.)

  return std::log (operator() (t));

  // We *could* just take the log of the result of operator(), but we can
  // instead be clever, removing one call to std::exp() in order to reduce
  // the number of libm function calls from 5 to 4:

  /*  double result, argerf = b * t + c[0];
  if (argerf >= 0)
    result = a[0] * (1 + erf(argerf));
  else
    result = a[0] * erfc(-argerf);

  argerf = b * t + c[1];
  if (argerf >= 0)
    result += a[1] * (1 + erf(argerf)) * std::exp(t * lambdadiff[0]);
  else
    result += a[1] * erfc(-argerf) * std::exp(t * lambdadiff[0]);

  argerf = b * t + c[2];
  if (argerf >= 0)
    result += a[2] * (1 + erf(argerf)) * std::exp(t * lambdadiff[1]);
  else
    result += a[2] * erfc(-argerf) * std::exp(t * lambdadiff[1]);

    return -lambda[0] * t + std::log(result);*/
}
