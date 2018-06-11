/* BOREXINO Reconstruction program
 *
 * Authors: A.Derbin; V.Muratova; O.Smirnov 
 * Maintainers: A.Derbin; V.Muratova; O.Smirnov
 *
 * Radial Energy correction module :
 * Q_corr=Q/fQ(Q,r)
 * Npe_corr=Npe/fNpe(Npe,r)
 * Nhits_corr=Nhits/fNhits(Nhits,r)
 * r - reconstructed event radius
 * 1dim correction functions are obtained from the MC for 4 
 * energies (200,400,800 and 1200 keV)
 */
#ifndef __BX_ENERGY_RECO_DBN_HH__
#define __BX_ENERGY_RECO_DBN_HH__

class TH1F;
#include "bx_base_module.hh"

#define nbins 55
#define NearlyZero 1E-10

class bx_energy_reco_dbn: public bx_base_module {
  public:
    bx_energy_reco_dbn ();
    virtual ~bx_energy_reco_dbn ();

    virtual void begin ();
    virtual bx_echidna_event* doit (bx_echidna_event *ev);
    virtual void end ();

  private:
    float f4_pe_dbn;
    float f4_npe_dbn;
    float f4_nhits_dbn;
    const static double SSpline_Y[3][55]; 
    const static double SSpline_M[3][55]; 
    float Spline(int32_t N, float X , const double *data,const double *m,double *y);
    TH1F *pHisto_npe;
    TH1F *pHisto_pe;
    TH1F *pHisto_nhits;
};

#endif
