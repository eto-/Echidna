/* BOREXINO Reconstruction program
 *
 * Author: Marco Pallavicini <Marco.Pallavicini@ge.infn.it>
 * Maintainer: Marco Pallavicini <Marco.Pallavicini@ge.infn.it>
 *
 * $Id: bx_pid_alphabeta.hh,v 1.9 2015/07/13 15:08:04 misiaszek Exp $
 * Calculation of several alpha/beta discrimination variables
 * Classifiers learning should be done outside echidna
 *
 */

#ifndef _BX_PID_ALPHABETA_H
#define _BX_PID_ALPHABETA_H

#include "bx_base_module.hh"
#include "bx_rec_general.hh"

#include "bx_alphabeta.hh"

class bx_pid_alphabeta: public bx_base_module {
  public:
    bx_pid_alphabeta ();
    virtual ~bx_pid_alphabeta () {}

    virtual void begin ();

    virtual bx_echidna_event* doit (bx_echidna_event *ev);

    virtual void end ();
    
  private:

    // helper function to search for nearest run number file with alpha beta parameters
    // for each run, the nearest previous run is used
    bool search_for_right_file(const std::string& path, int run, std::string& fname);

    // tailtot calculation used for mlp algorithm
    float tailtot_ab_mlp(float timecut, std::vector<double> &current_sample);
    
    // calculate rms & kurtosis
    float rec_rms(std::vector<double> &current_sample, float& kurtosis);

    // alpha/beta object
    bx_alphabeta fAB;
    
    // minimum clusters size to consider alpha/beta 
    int fMinClusterSize;
    
    // coeff for gatti and lkl
    int fNBins;
    int fBinSizeNs;
    float fGattiWeights[5000];
    float fGattiWeightsC[5000];
    float fLkl[5000];
    float fLklC[5000];   
    
    // histos helper for gatti and lkl
    TH1F *hSample;
    TH1F *hSampleC;
};

#endif
