/* BOREXINO Reconstruction program
 *
 * Author: Sergey Sukhotin <sukhotin@in2p3.fr>
 *         Evgeny Litvinovich <litvinov@lngs.infn.it>
 * Maintainer: Evgeny Litvinovich <litvinov@lngs.infn.it>
 *
 * $Id: bx_position_reco_msk.hh,v 1.12 2007/10/25 16:20:01 litvinov Exp $
 *
 * Moscow event's position & energy reconstruction
*/

#ifndef _BX_POSITION_RECO_MSK_HH
#define _BX_POSITION_RECO_MSK_HH
#include "bx_base_module.hh"
#include "TVector3.h"
#include "TH1F.h"

class bx_echidna_event;
class TMinuit;

class bx_position_reco_msk : public bx_base_module {
  public:
    bx_position_reco_msk();
    virtual ~bx_position_reco_msk() {}
  
    // Operations for the framework 
    virtual void begin();
    virtual bx_echidna_event* doit(bx_echidna_event *ev);
    virtual void end();
     double msk_fcn_t(double *x);

  private:
     double f8_time_cut;
     double f8_landau_mpv;
     double f8_landau_sigma;

      float f4_ref_index;

     std::vector<TVector3> *f4_pmt_positions;
     std::vector<bool> *b_pmt_cone;
     std::vector<int> i4_disabled_channels;

     TMinuit *p_minuit;
     bx_echidna_event *p_fit_ev;
     int i4_fit_cluster;
     int i4_n_minuit_warns;

     TH1F *msk_x, *msk_y, *msk_z, *msk_t;
     TH1F *msk_r;
};

#endif

/*
 * $Log: bx_position_reco_msk.hh,v $
 * Revision 1.12  2007/10/25 16:20:01  litvinov
 * decoupling "Moscow" energy and position reconstruction onto 2 independent modules
 *
 * Revision 1.11  2007-05-24 09:20:14  litvinov
 * getting status of Migrad convergency and putting it into log if problems
 *
 * Revision 1.10  2007-05-05 16:27:42  litvinov
 * Now the module returns reconstructed energy in terms of i4_npe and f4_charge.
 * i4_nhits still missed as IMHO useless
 *
 * Revision 1.9  2007-05-02 16:28:26  litvinov
 * change getter name from get_bad_channels to get_disabled_channels
 * according to last modification in bx_detector
 *
 * Revision 1.8  2007-05-02 11:41:44  litvinov
 * optimization for faster calculations
 *
 * Revision 1.7  2007-04-29 15:10:19  litvinov
 * moving parameters used by module into echidna.cfg
 *
 * Revision 1.6  2007-04-27 18:17:30  litvinov
 * two fcn-functions instead of one;
 * coming back po Poissonian statistics;
 * preparation to return npe/nhits/charge
 *
 * Revision 1.5  2007-04-23 11:16:28  litvinov
 * check for bad channels via detector_interface
 *
 * Revision 1.4  2006-11-17 17:51:29  litvinov
 * Module has been fully rewritten. No more groups of PMTs again.
 * Many new features. Work is in progress.
 *
 * Revision 1.3  2006-01-11 12:01:15  litvinov
 * Implementation of geometrical sharing of the PMTs to finite number of groups. This approx. 6-7 times reduces machine time.
 *
 * Revision 1.2  2005/10/16 16:25:34  litvinov
 * updated to call the Poisson predefined in ROOT
 *
 * Revision 1.1  2004/12/13 14:04:39  litvinov
 * Moscow event's position & energy reconstruction module is implemented
 *
 *
 */
