/**
 * @file ScintPDF.h
 * @brief Defines the scintillator PDF function.
 */
#ifndef _SCINT_PDF_H
#define _SCINT_PDF_H

class ScintPDF {
public:
  enum pdf_type { STANDARD, DIFFUSER_BALL };
private:
  double q[3];
  ///< Amplitudes of the fast and slow components of scintillator light decay.
  ///< Must sum to one.

  double sigma;
  ///< Gaussian spreading of the scintillator PDF due to scattering, PMt
  ///< jitter, etc.  In ns.

  double tau[3];
  ///< Time constants of the fast and slow components of scintillator light
  ///< decay.  In ns.

  double refl[3];
  ///< Coefficients for the gaussian representing the reflection peak.
  
  double a[3], b, c[3], lambda[3], lambdadiff[2];
  ///< Linear coefficients that are functions of q[i], sigma and tau[i].
  ///< Calculated once for efficiency.

  double gcoeff;
  ///< Coefficient for the gaussian representing the reflection peak.
  ///< Calculated once for efficiency.
  
  pdf_type type_to_use;
  ///< Allows to select within a run between different pdf types

public:
  bool have_reflection;

  ScintPDF()  : type_to_use(STANDARD), have_reflection(false) { }
  ScintPDF(const ScintPDF & pdf) { Load(pdf); }
  ScintPDF & operator= (const ScintPDF & pdf) { Load(pdf); return *this; }
  ~ScintPDF() { }

  void SetQ0(double Q0) { q[0] = Q0; q[1] = 1. - Q0; q[2] = 0;}
  void SetQ01(double Q0, double Q1) {
    q[0] = Q0; q[1] = Q1; q[2] = 1 - Q0 - Q1;
  }
  void SetSigma(double Sigma) { sigma = Sigma; }
  void SetTau(double Tau0, double Tau1, double Tau2 = 1.0)
  { tau[0] = Tau0; tau[1] = Tau1; tau[2] = Tau2; }
  // Tau2 defaults to 1.0 to avoid division by zero.
  void CalculateCoefficients();
  ///< MUST set Q0 [Q1], Sigma & Tau, then run this, before using operator()
  ///< or Log
  void Load(const ScintPDF & pdf); ///< Copy coefficients from another PDF

  void SetReflectionSigma(double Sigma) { refl[2] = Sigma; }
  void SetReflectionCenter(double Center) { refl[1] = Center; }
  void SetReflectionAmplitude(double Ampl) { refl[0] = Ampl; }

  double GetQ0() const { return q[0]; }
  double GetQ1() const { return q[1]; }
  double GetQ2() const { return q[2]; }
  double GetSigma() const { return sigma; }
  double GetTau0() const { return tau[0]; }
  double GetTau1() const { return tau[1]; }
  double GetTau2() const { return tau[2]; }

  double GetReflectionSigma() const { return refl[2]; }
  double GetReflectionCenter() const { return refl[1]; }
  double GetReflectionAmplitude() const { return refl[0]; }
  
  void SetPDFType(pdf_type newtype) { type_to_use = newtype; }
  pdf_type GetPDFType() const {return type_to_use; }

  double operator() (double t) const;
  ///< @brief Return the value of the normalized PDF at a given time,
  ///< for given path lengths in scintillator and buffer and given NPE

  double Value      (double t) const
  { return operator() (t); }
  ///< @brief Synonym for operator() for people who don't want that syntax

  double Log        (double t) const;
  ///< @brief Return the natural log of the normalized PDF at a given time,
  ///< for given path lengths in scintillator and buffer and given NPE
};

#endif
