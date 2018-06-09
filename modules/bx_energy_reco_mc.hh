/* BOREXINO Reconstruction program
 *
 * Author: Marcin Misiaszek <marcin.misiaszek@lngs.infn.it>
 * Maintainer: Marcin Misiaszek <marcin.misiaszek@lngs.infn.it>
 *
 * Energy reconstruction module based on 2dim calibration 
 * tables for nhits, charge and npe. 
 *
 * 
 */
#ifndef _BX_ENERGY_RECO_MC_H
#define _BX_ENERGY_RECO_MC_H

#include "bx_base_module.hh"
#include "cmap.hh"
#include <string>

class interpolator2d;

class bx_energy_reco_mc: public bx_base_module {
  public:
    bx_energy_reco_mc ();
    virtual ~bx_energy_reco_mc ();

    virtual void begin ();
    virtual bx_echidna_event* doit (bx_echidna_event *ev);
    virtual void end ();

  private:
    const static double nhits[5][3], charge[5][3], npe[5][3];    
    interpolator2d *interpol; 
    float m_get_energy(float radius, float nhits);
    int method;
};

#endif
