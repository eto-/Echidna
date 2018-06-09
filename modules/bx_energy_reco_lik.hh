/* BOREXINO Reconstruction program
 *
 * Author: Marcin Misiaszek <marcin.misiaszek@lngs.infn.it>
 * Maintainer: Marcin Misiaszek <marcin.misiaszek@lngs.infn.it>
 *
 * Energy reconstruction module based on maximum likelihood  
 * method. Temporarily using Igor Machulin - Geant4 fit for the cone
 * (ported from bx_position_reco_msk).
 *
 * 
 */
#ifndef _BX_ENERGY_RECO_lik_H
#define _BX_ENERGY_RECO_lik_H

#include "bx_base_module.hh"
#include <string>

class TMinuit;
class bx_dbi;

class vector3 {
    public:
	vector3() {};
	vector3( double xx, double yy, double zz) { x[0] = xx; x[1] = yy; x[2] = zz;};
	virtual ~vector3() {};
	double Mag() { return sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]); }
	double& operator[] (const int i) {return x[i]; }
	vector3 operator- (const vector3 &rhs) {
	    vector3 temp(this->x[0] - rhs.x[0], this->x[1] - rhs.x[1], this->x[2] - rhs.x[2]);
	    return temp;
	    }
	vector3 operator* (const double c) { 
	    vector3 temp(this->x[0] * c, this->x[1] * c, this->x[2] * c);
	    return temp; 
	    }   
	double operator* (const vector3 &rhs) {
	    return ((this->x[0] * rhs.x[0]) + (this->x[1] * rhs.x[1]) + (this->x[2] * rhs.x[2]));
	    }
    private:
    	double x[3];
};


class bx_energy_reco_lik: public bx_base_module {
  public:
    bx_energy_reco_lik ();
    virtual ~bx_energy_reco_lik ();

    virtual void begin ();
    virtual bx_echidna_event* doit (bx_echidna_event *ev);
    virtual void end ();

    double lik_fcn (double *x);
       
  private:
    float f4_ref_index;
    float f4_lg_entry_radius;
    float f4_cathode_radius;
    int method;
    
    double f8_prob, expected_charge, f8_result;
    
    std::vector<vector3> *f4_pmt_positions;
    std::vector<vector3> *f4_pmt_normal_vectors;
    std::vector<unsigned char> *u1_collected_charge;
    std::vector<double> *v_omega_attenuation;
    std::vector<bool> *b_pmt_cone;

    // to be used in lik_fcn (see minuit root workaround)
    TMinuit *p_minuit;
    bx_echidna_event *p_fit_ev;
    int i4_fit_cluster;
};

#endif
