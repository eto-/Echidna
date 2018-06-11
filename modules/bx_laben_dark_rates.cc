/* BOREINO Reconstruction program
 *
 * Author:  Livia Ludhova<livia.ludhova@mi.infn.it>
 *
 */

#include <cmath>
#include <iostream>
#include <vector>
#include <algorithm>
#include "bx_laben_dark_rates.hh"
#include "bx_echidna_event.hh"
#include "bx_dbi.hh"
#include "db_run.hh"
#include "db_profile.hh"
#include "db_channel.hh"
#include "barn_interface.hh"
#include "constants.hh"
#include "TF1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"

bx_laben_dark_rates::bx_laben_dark_rates () : bx_base_module("bx_laben_dark_rates", bx_base_module::main_loop) {
  require_event_stage (bx_detector::laben, bx_base_event::clustered);
  require_trigger_type (bx_trigger_event::random);
}
      
  //BEGIN
void bx_laben_dark_rates::begin () {

  db_write = get_parameter ("db_write").get_bool ();

  N_tt64 = 0; 
  sum_dec_hits_per_lg  = 0;
  error_sum_dec_hits_per_lg  = 0;

  //reda gate width from DB
  db_run& run_info = bx_dbi::get ()->get_run ();
  gate_width   = run_info.get_laben_gate_width ();	     
  get_message(bx_message::log) << "Gate width " << gate_width << " ns" << dispatch; 
  
  //define histos and their storing in barn
  n_dechits = new TH1F ("n_dechits","N decoded hits spectrum of TT64",300,-0.5,299.5);
  barn_interface::get ()->store (barn_interface::file, n_dechits, this);

  npmts_win1 = new TH1F ("npmts_win1","npmts_win1 spectrum of TT64",150,-0.5,149.5);
  barn_interface::get ()->store (barn_interface::file, npmts_win1, this);
   
  npmts_win2 = new TH1F ("npmts_win2","npmts_win2 spectrum of TT64",150,-0.5,149.5);
  barn_interface::get ()->store (barn_interface::file, npmts_win2, this);
      
  
 
}  



  //DOIT
bx_echidna_event* bx_laben_dark_rates::doit (bx_echidna_event *ev) {
  
  //  //from random triggers
  N_tt64 ++;
    
  double Nlg = ev->get_laben ().get_n_live_pmts () - ev->get_laben ().get_invalid_pmts ();

  double Nhits = ev->get_laben ().get_decoded_nhits ();
  n_dechits->Fill(Nhits);
  if(Nlg) {
    sum_dec_hits_per_lg += (double) Nhits / (double) Nlg;
   
    double error_1_event = sqrt((double) Nhits) / (double) Nlg;
   
    error_sum_dec_hits_per_lg += pow(error_1_event,2);
  }

  const std::vector<int32_t>& NPmts_w1 = ev->get_laben ().get_npmts_win1 ();
  const std::vector<int32_t>& NPmts_w2 = ev->get_laben ().get_npmts_win2 ();

  for(int32_t i = 0; i < (int32_t) NPmts_w1.size (); i ++) npmts_win1->Fill (NPmts_w1[i]);
  for(int32_t i = 0; i < (int32_t) NPmts_w2.size (); i ++) npmts_win2->Fill (NPmts_w2[i]);


  return ev;     
  
}


double fitfunction(double *x, double *par) {
  return (par[0] * TMath::Poisson(x[0],par[1]) + par[2] * exp(- x[0]/par[3]));
}
 
double fitexp(double *x, double *par) {
  return (par[0] * exp(- x[0]/par[1]));
}
 

  //END
void bx_laben_dark_rates::end () {

  
  double mean_dark_rate_per_used_pmt = -10;
  double error_mean_dark_rate_per_used_pmt = -10;

  double mu_win1 = -10;
  double mu_win2 = -10;

  double pois_const_win1 = -10;
  double mu_fit_win1 = -10;
  double exp_const_win1 = -10;
  double tau_win1 = -10;
				  
  double error_pois_const_win1 = -10;
  double error_mu_fit_win1 = -10;
  double error_exp_const_win1 = -10;
  double error_tau_win1 = -10;

  double pois_const_win2 = -10;
  double mu_fit_win2 = -10;
  double exp_const_win2 = -10;
  double tau_win2 = -10;
				  
  double error_pois_const_win2 = -10;
  double error_mu_fit_win2 = -10;
  double error_exp_const_win2 = -10;
  double error_tau_win2 = -10;

  //we analyse only when TT64 ispresent and histos have some inputs
  if(N_tt64 && npmts_win1->Integral(1,3) && npmts_win2->Integral(1,3)) { 
    
    mean_dark_rate_per_used_pmt = sum_dec_hits_per_lg / (N_tt64 * gate_width * 1e-9);
    error_mean_dark_rate_per_used_pmt = sqrt(error_sum_dec_hits_per_lg) / (N_tt64 * gate_width * 1e-9);

    mu_win1 = -1 * log(npmts_win1->Integral(1,1)/npmts_win1->Integral(1,3)); ////log (1. - N(bin 0) /N(bin 0 - bin 4)
    mu_win2 = -1 * log(npmts_win2->Integral(1,1)/npmts_win2->Integral(1,3)); ////log (1. - N(bin 0) /N(bin 0 - bin 4)

    //fit functions
    TF1 *fitexpo = new TF1("fitexpo",fitexp,5.,15.,2);

    TF1 *fitfun = new TF1("fitfun",fitfunction,0.,15.,4);
    fitfun->SetParName(0,"Pois Const");
    fitfun->SetParName(1,"Mu");
    fitfun->SetParName(2,"Exp Const");
    fitfun->SetParName(3,"Tau");
    
    //FIT NPMTS_WIN1
    //fitting exp part to get input values for exponential
    fitexpo->SetParameter(0,npmts_win1->Integral(3,3));
    fitexpo->SetParameter(1,2);
    npmts_win1->Fit("fitexpo","LQN","",5,15);
    
    //fitting Poiss + exp
    fitfun->SetParameter(0,npmts_win1->Integral(1,3));
    fitfun->SetParameter(1,mu_win1);
    fitfun->SetParameter(2,fitexpo->GetParameter(0));
    fitfun->SetParameter(3,fitexpo->GetParameter(1));
    
    npmts_win1->Fit("fitfun","LQN","",0,15);
    
    pois_const_win1 = fitfun->GetParameter(0);
    mu_fit_win1 = fitfun->GetParameter(1);
    exp_const_win1 = fitfun->GetParameter(2);
    tau_win1 = fitfun->GetParameter(3);
    
    error_pois_const_win1 = fitfun->GetParError(0);
    error_mu_fit_win1 = fitfun->GetParError(1);
    error_exp_const_win1 = fitfun->GetParError(2);
    error_tau_win1 = fitfun->GetParError(3);

    //FIT NPMTS_WIN2
    //fitting exp part to get input values for exponential
    fitexpo->SetParameter(0,npmts_win2->Integral(3,3));
    fitexpo->SetParameter(1,2);
    npmts_win2->Fit("fitexpo","LQN","",5,15);
    
    //fitting Poiss + exp
    fitfun->SetParameter(0,npmts_win2->Integral(1,3));
    fitfun->SetParameter(1,mu_win2);
    fitfun->SetParameter(2,fitexpo->GetParameter(0));
    fitfun->SetParameter(3,fitexpo->GetParameter(1));
    
    npmts_win2->Fit("fitfun","LQN","",0,15);

    pois_const_win2 = fitfun->GetParameter(0);
    mu_fit_win2 = fitfun->GetParameter(1);
    exp_const_win2 = fitfun->GetParameter(2);
    tau_win2 = fitfun->GetParameter(3);
    
    error_pois_const_win2 = fitfun->GetParError(0);
    error_mu_fit_win2 = fitfun->GetParError(1);
    error_exp_const_win2 = fitfun->GetParError(2);
    error_tau_win2 = fitfun->GetParError(3);
  }


  
    get_message(bx_message::log) << " mean_dark_rate_per_used_pmt = " <<  mean_dark_rate_per_used_pmt  << " +- " << error_mean_dark_rate_per_used_pmt << dispatch; 

    get_message(bx_message::log) << " mu_win1 = " <<  mu_win1 << dispatch; 
    get_message(bx_message::log) << " mu_win2 = " <<  mu_win2 << dispatch; 

    get_message(bx_message::log) << " pois_const_win1 = " << pois_const_win1 << " +- " << error_pois_const_win1  << dispatch;
    get_message(bx_message::log) << " mu_fit_win1 = " << mu_fit_win1 << " +- " << error_mu_fit_win1 << dispatch;
    get_message(bx_message::log) << " exp_const_win1 = " << exp_const_win1 << " +- " <<  error_exp_const_win1 << dispatch;
    get_message(bx_message::log) << " tau_win1 = " << tau_win1 << " +- " << error_tau_win1 << dispatch;

    get_message(bx_message::log) << " pois_const_win2 = " << pois_const_win2 << " +- " << error_pois_const_win2  << dispatch;
    get_message(bx_message::log) << " mu_fit_win2 = " << mu_fit_win2 << " +- " << error_mu_fit_win2 << dispatch;
    get_message(bx_message::log) << " exp_const_win2 = " << exp_const_win2 << " +- " <<  error_exp_const_win2 << dispatch;
    get_message(bx_message::log) << " tau_win2 = " << tau_win2 << " +- " << error_tau_win2 << dispatch;

    db_run& run_info = bx_dbi::get()->get_run ();

  if(db_write && N_tt64 > 3000 && !run_info.is_laben_dark_rates_parametrisation_present   ()){ 
    //    db_run& run_info = bx_dbi::get()->get_run ();

    //set visitors
    run_info.set_mean_dark_rate_per_used_pmt (mean_dark_rate_per_used_pmt, this);
    run_info.set_error_mean_dark_rate_per_used_pmt (error_mean_dark_rate_per_used_pmt, this);

    run_info.set_mu_win1 (mu_win1, this);
    run_info.set_mu_win2 (mu_win2, this);

    run_info.set_pois_const_win1 (pois_const_win1, this);            
    run_info.set_error_pois_const_win1 (error_pois_const_win1, this);            
    run_info.set_mu_fit_win1(mu_fit_win1, this);                
    run_info.set_error_mu_fit_win1(error_mu_fit_win1, this);                
    run_info.set_exp_const_win1(exp_const_win1, this);
    run_info.set_error_exp_const_win1(error_exp_const_win1, this);
    run_info.set_tau_win1(tau_win1, this);
    run_info.set_error_tau_win1(error_tau_win1, this);

    run_info.set_pois_const_win2 (pois_const_win2, this);            
    run_info.set_error_pois_const_win2 (error_pois_const_win2, this);            
    run_info.set_mu_fit_win2(mu_fit_win2, this);                
    run_info.set_error_mu_fit_win2(error_mu_fit_win2, this);                
    run_info.set_exp_const_win2(exp_const_win2, this);
    run_info.set_error_exp_const_win2(error_exp_const_win2, this);
    run_info.set_tau_win2(tau_win2, this);
    run_info.set_error_tau_win2(error_tau_win2, this);

    
    //write to the DB
    run_info.write_laben_dark_rates_parametrisation (true, this);
    get_message(bx_message::log) << "Writing to DB" << dispatch; 
  }

}
 




 
 
 
 



