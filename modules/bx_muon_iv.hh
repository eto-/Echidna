/* BOREXINO Reconstruction program
 *
 * Author: Livia Ludhova<Livia.Ludhova@mi.infn.it>
 * 
 * Maintainer: Livia Ludhova<Livia.Ludhova@mi.infn.it>
 */

#ifndef _BX_MUON_IV_HH
#define _BX_MUON_IV_HH

#include "bx_base_module.hh"
#include "bx_rec_general.hh"

class bx_laben_event;
class TH1F;
class TH2F;

class bx_muon_iv: public bx_base_module {
  public:
  // this section if fixed; do not edit
    bx_muon_iv ();
    virtual ~bx_muon_iv () {}

    virtual void begin ();
    virtual bx_echidna_event* doit (bx_echidna_event *ev);
    virtual void end ();
  private:
  
  //std::vector<double> *my_vector;
  //bool m_check_this_and_that(const bx_laben_event& er);
  int32_t muon_event;
  int32_t npe_satur_thresh;
  double hits_satur_thresh; 
  double pmt_hitted_thresh;
  double pmt_satur_thresh;
  double mean_time_thresh;
  double notcone_cone_thresh;
  double N_ordinary_lg ; 
  TH1F* notcone_cone;
  TH1F* hit_lg_junk; 
  TH1F* satur_hit_lg_junk; 
};

#endif
