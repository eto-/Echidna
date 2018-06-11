/* BOREXINO Reconstruction program
 *
 * Author: Koun Choi <kounchoi@hawaii.edu>
 *
 * $Id: bx_new_charge_weight.hh,v 1.3 2015/08/17 17:20:06 misiaszek Exp $
 *
 */

#ifndef _BX_NEW_CHARGE_WEIGHT_HH
#define _BX_NEW_CHARGE_WEIGHT_HH

#include "bx_echidna_event.hh"
#include "bx_base_module.hh"
#include "bx_rec_general.hh"
#include <list>
#include <vector>

class db_channel_laben;

class bx_new_charge_weight: public bx_base_module {
  public:
    bx_new_charge_weight ();
    virtual ~bx_new_charge_weight () {}

    virtual void begin ();
    virtual bx_echidna_event* doit (bx_echidna_event *ev);
    virtual void end ();
    
  private:

    uint8_t *p_disabled_lg;
    uint32_t i4_n_disabled_channels, i4_n_disabled_charge;
    const db_channel_laben **ch_info_v;
    bool *p_disabled_pmts_lg;
    bool *p_disabled_charge_lg;

    int32_t nlg;
    float xp,yp,zp;
    float QE;
    float R_pmt;
    float d_pmt;
	float angt_pmts;
	float angt_charge;
	float QEt_pmts;
	float QEt_charge;
	float allt_pmts;
	float allt_charge;

};
#endif

/*  
 *  $Log: bx_new_charge_weight.hh,v $
 *  Revision 1.3  2015/08/17 17:20:06  misiaszek
 *  Id & log keywords added
 *
 */
