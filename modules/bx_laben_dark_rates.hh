/* BOREXINO Reconstruction program
 *
 * Author: Livia Ludhova<livia.ludhova@mi.infn.it>
 * 
*/

#ifndef _BX_LABEN_DARK_RATES_HH
#define _BX_LABEN_DARK_RATES_HH
#include "bx_rec_general.hh"
#include "bx_base_module.hh"
#include "TROOT.h"
#include "TH1F.h"
#include "TH2F.h"
#include <map>

class bx_echidna_event;
class bx_laben_cluster;
class bx_laben_event;

class bx_laben_dark_rates : public bx_base_module {
  public:
    bx_laben_dark_rates ();
    virtual ~bx_laben_dark_rates () {}
      
    // Operations for the framework 
    virtual void begin ();
    virtual bx_echidna_event* doit (bx_echidna_event *ev);
    virtual void end ();

  private:
  
  bool db_write;
  int N_tt64;
  double sum_dec_hits_per_lg;
  double error_sum_dec_hits_per_lg;


  //to be read from DB
  float gate_width;

  //histos
  TH1F* n_dechits; 
  TH1F* npmts_win1; 
  TH1F* npmts_win2; 

};

#endif
