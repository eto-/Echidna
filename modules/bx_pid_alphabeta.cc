/* BOREXINO Reconstruction program
 *
 * Author: Marco Pallavicini <Marco.Pallavicini@ge.infn.it>
 * Maintainer: Marco Pallavicini <Marco.Pallavicini@ge.infn.it>
 *
 * $Id: bx_pid_alphabeta.cc,v 1.26 2015/08/26 11:21:18 misiaszek Exp $
 * Calculation of several alpha/beta discrimination variables
 *
 */

#include "bx_pid_alphabeta.hh"
#include "bx_echidna_event.hh"
#include "barn_interface.hh"
#include "bx_alphabeta.hh"
#include "bx_dbi.hh"
#include "db_run.hh"
#include "PSD_MLPv8.C"
#include <algorithm>

// ctor
bx_pid_alphabeta::bx_pid_alphabeta (): bx_base_module("bx_pid_alphabeta", bx_base_module::main_loop) {
  require_event_stage (	bx_detector::laben, bx_base_event::reconstructed); 
//  require_trigger_type (bx_trigger_event::neutrino);
//  require_trigger_type (bx_trigger_event::neutron);
}


// module interface
void bx_pid_alphabeta::begin () {
	
  const std::string path = get_parameter("alphabeta_file_path").get_string();
  int32_t run = bx_dbi::get()->get_run().get_number();
	
  std::string fName="";
  if ( search_for_right_file(path, run, fName) ) {
        
	  // fAB.load(fName);
	  fAB.set_time_histogram();
	
	
	  FILE *fp2 = fopen(fName.c_str(),"rt");
        
	  if (!fp2) {
		  get_message(bx_message::critic) << "Alpha beta file" << fName << " not found for run=" << run << dispatch;
	  } else {
		  get_message(bx_message::info) << "Alpha beta file" << fName << " successfully opened for run=" << run << dispatch;
	  }

	  fscanf(fp2,"%d",&fNBins);
	  fscanf(fp2,"%d",&fBinSizeNs);
		
	  int32_t dummy; 

	  for(int32_t i=1; i<=fNBins; i++) {
		  fscanf(fp2,"%d %f",&dummy, &(fGattiWeights[i]));
		  if ( dummy != i ) {
			  get_message(bx_message::warn) << "1 Error reading file i=" << dummy << " != " << i << dispatch;
			  return;
		  }
	  }

	  for(int32_t i=1; i<=fNBins; i++) {
		  fscanf(fp2,"%d %f",&dummy, &(fGattiWeightsC[i]));
		  if ( dummy != i ) {
			  get_message(bx_message::warn) << "2 Error reading file i=" << dummy << " != " << i << dispatch;
			  return;
	  	}
	  }

	  for(int32_t i=1; i<=fNBins; i++) {
		  fscanf(fp2,"%d %f",&dummy, &(fLkl[i]));
		  if ( dummy != i ) {
			  get_message(bx_message::warn) << "3 Error reading file i=" << dummy << " != " << i << dispatch;
			  return;
		  }
	  }

  	for(int32_t i=1; i<=fNBins; i++) {
	  	fscanf(fp2,"%d %f",&dummy, &(fLklC[i]));
		  if ( dummy != i ) {
			  get_message(bx_message::warn) << "4 Error reading file i=" << dummy << " != " << i << dispatch;
			  return;
		  }
	  }

	  fclose(fp2);
	  hSample = new TH1F("TimeSample","xx",fNBins,0.,(float)fNBins*fBinSizeNs);
	  hSampleC = new TH1F("TimeSampleTimeC","xx",fNBins,0.,(float)fNBins*fBinSizeNs);	
  }
	
  /*
  if (fAB.isReady() == false) {
	  get_message(bx_message::warn) << "Cannot load alpha beta file. File " << fName << " not found " << dispatch;
	  get_message(bx_message::warn) << "Tailtot only will be computed " << dispatch;
  }
  */

  fMinClusterSize = get_parameter("minimum_cluster_size").get_int();
}

// doit
bx_echidna_event* bx_pid_alphabeta::doit (bx_echidna_event *ev) {

        // get number of clusters
	int32_t nc = ev->get_laben().get_nclusters();
	if (nc <= 0) return ev;

	// loop on clusters and compute alpha-beta variables
	for(int32_t iclus=0; iclus<nc; iclus++) {

		// get echidna ab_cluster
		bx_laben_ab_cluster& ab = ev->get_laben().get_ab_cluster(iclus);

		// set variables to invalid values
		for(int32_t i=0; i<10; i++) ab.v_tailtot[i] = 2.;
		ab.f4_lkl       = -100.;
		ab.f4_gatti     = -100.;
		ab.f4_lklc      = -100.;
		ab.f4_gattic    = -100.;
		ab.f4_rise_time = -100.;
		hSample->Reset();
		hSampleC->Reset();

		// check that the cluster has some hits
		int32_t nhits = ev->get_laben().get_cluster(iclus).get_clustered_nhits();
		if ( nhits < fMinClusterSize )
			continue;

		// set current cluster
		fAB.set(*ev, iclus);

		// compute alpha/beta variables: tailtot with tail between 40 and 130, step 10
		for(int32_t index=0; index<10; index++) {
			float tail = 40. + index*10.;
			float aa = fAB.tailtot( tail );
			if (aa > 0. && aa <= 1. )
				ab.v_tailtot[index] = aa;
		}

		// Fill histo of times (skip the first 2 and shift to avoid noise hit effects)
		float shift_time = ab.get_rec_hit(2).get_time();
		for (int32_t  i=2; i < nhits; i++) hSample->Fill( ab.get_rec_hit(i).get_time()-shift_time );

		// Normalize histo
		float aa = hSample->Integral();
		if (aa>0.000000001)
			hSample->Scale(1./aa);
		else
			get_message(bx_message::warn) << "Sample histogram area too small " << dispatch;


		// Fill Cumulative histo
		for(int32_t i=1; i<=hSample->GetNbinsX(); i++) hSampleC->SetBinContent(i, hSample->Integral(1,i));

		// Normalize cumulative histo
		aa = hSampleC->Integral();
		//if (aa>0.000000001)
		//hSampleC->Scale(1./aa);
		//else
		//get_message(bx_message::warn) << "SampleC histogram area too small " << dispatch;

		float gatti     = 0.;
		float gattic    = 0.;
		float lkl       = 0.;
		float lklc      = 0.;
		for(int32_t i=1; i<=hSample->GetNbinsX(); i++) {
			gatti  +=  hSample ->GetBinContent(i) * fGattiWeights [i];
			gattic +=  hSampleC->GetBinContent(i) * fGattiWeightsC[i];
			lkl    +=  hSample ->GetBinContent(i) * fLkl          [i];
			lklc   +=  hSampleC->GetBinContent(i) * fLklC         [i];
		}

		// Compute rise time of the integrated curve
		int32_t bin = 0;
		while (hSampleC->GetBinContent(++bin) < 0.9) ;
		float rise_time = (bin-1);

		ab.f4_gatti     = gatti;
		ab.f4_lkl       = lkl;
		ab.f4_gattic    = gattic;
		ab.f4_lklc      = lklc;
		ab.f4_rise_time = rise_time;

		/////////////////////////////////////////////////////////////////////////////////////////////
	        //new tailtot for a/b discrimination with MLP algorithm
	        //
 	
 		//fill sample of rec_hit times
    		std::vector<double> sample;
    		std::vector<double> sample_c11;
		for (int32_t  i=2; i < nhits; i++) {
                  double time = ab.get_rec_hit(i).get_time()-shift_time;
                  if( time > 800 ) continue;
 		  if( time < 0.1 ) continue;
                  sample.push_back( time );
    		}

		for (int32_t  i=0; i < nhits; i++) {
                  double time = ab.get_rec_hit(i).get_time();
                  if( time > 110 ) continue;
                  sample_c11.push_back( time );
    		}

                //calculate tailtot for a/b;
	        ab.v_tailtot_ab_mlp[0] = tailtot_ab_mlp(5,   sample);
                ab.v_tailtot_ab_mlp[1] = tailtot_ab_mlp(10,  sample);
                ab.v_tailtot_ab_mlp[2] = tailtot_ab_mlp(20,  sample);
                ab.v_tailtot_ab_mlp[3] = tailtot_ab_mlp(40,  sample);
                ab.v_tailtot_ab_mlp[4] = tailtot_ab_mlp(80,  sample);
                ab.v_tailtot_ab_mlp[5] = tailtot_ab_mlp(120, sample);
                ab.v_tailtot_ab_mlp[6] = tailtot_ab_mlp(180, sample);
                ab.v_tailtot_ab_mlp[7] = tailtot_ab_mlp(240, sample);
                ab.v_tailtot_ab_mlp[8] = tailtot_ab_mlp(310, sample);
                ab.v_tailtot_ab_mlp[9] = tailtot_ab_mlp(380, sample);

                //calculate tailtot for c11/b;
	        ab.v_tailtot_c11_mva[0] = tailtot_ab_mlp(5,  sample_c11);
                ab.v_tailtot_c11_mva[1] = tailtot_ab_mlp(10, sample_c11);
                ab.v_tailtot_c11_mva[2] = tailtot_ab_mlp(15, sample_c11);
                ab.v_tailtot_c11_mva[3] = tailtot_ab_mlp(20, sample_c11);
                ab.v_tailtot_c11_mva[4] = tailtot_ab_mlp(25, sample_c11);
                ab.v_tailtot_c11_mva[5] = tailtot_ab_mlp(30, sample_c11);
                ab.v_tailtot_c11_mva[6] = tailtot_ab_mlp(35, sample_c11);
                ab.v_tailtot_c11_mva[7] = tailtot_ab_mlp(40, sample_c11);
                ab.v_tailtot_c11_mva[8] = tailtot_ab_mlp(45, sample_c11);
                ab.v_tailtot_c11_mva[9] = tailtot_ab_mlp(50, sample_c11);
                
                //calculate rms & kurtosis
		ab.rms = rec_rms(sample, ab.kurtosis);
		ab.rms_c11 = rec_rms(sample_c11, ab.kurtosis_c11);

		//calculate MLP for ab ver. 8
		ab.mlp_ab = PSD_MLPv8(	ab.v_tailtot_ab_mlp[0], 
					ab.v_tailtot_ab_mlp[1], 
					ab.v_tailtot_ab_mlp[2], 
					ab.v_tailtot_ab_mlp[3], 
					ab.v_tailtot_ab_mlp[4], 
					ab.v_tailtot_ab_mlp[5], 
					ab.v_tailtot_ab_mlp[6], 
					ab.v_tailtot_ab_mlp[7], 
					ab.v_tailtot_ab_mlp[8], 
					ab.v_tailtot_ab_mlp[9], 
					ab.rms,
					ab.kurtosis, 
					ev->get_laben().get_cluster(iclus).get_mean_time() );

	} // end of loop on clusters

	return ev;
}

void bx_pid_alphabeta::end () {
}

// search right learning file
// returns true if found
// store file name into variable fname
bool bx_pid_alphabeta::search_for_right_file(const std::string& path, int32_t run, std::string& fname) {
	int32_t tmp = run;
	while( tmp > 2000 ) {
		char fs[1000];
		sprintf(fs,"%s/Run%06d.w",path.c_str(),tmp);
		FILE *fp=fopen(fs,"rt");
		if (fp) {
			fclose(fp);
			fname = fs;
			get_message(bx_message::info) << "Found File =  " << fname << dispatch;
			return true;
		}
		tmp--;
	}
	return false;
}

float bx_pid_alphabeta::tailtot_ab_mlp(float timecut, std::vector<double> &current_sample) {
        if (current_sample.size() == 0) {
                return 0.;
        }
        int32_t tail=0;
        int32_t size = current_sample.size();
        for(int32_t i=0; i<size; i++) {
                float t = current_sample[i];
                if ((t)>=timecut)
                        tail++;
        }
        return ((float)tail)/((float)size);
}

float bx_pid_alphabeta::rec_rms(std::vector<double> &current_sample, float& kurtosis) {
        if (current_sample.size() == 0) {
                //_message(bx_message::error) << GetName() << " You must set an event before filtering!" << dispatch;
                return 0.;
        }
        float sum = 0.;
        int32_t size = current_sample.size();
        for(int32_t i=0; i<size; i++) {
                sum += current_sample[i];
        }
        float mean = ((float)sum)/((float)size);
        float m2 = 0.;
        for(int32_t i=0; i<size; i++) {
                m2 += pow(mean - current_sample[i], 2);
        }
        float m4 = 0.;
        for(int32_t i=0; i<size; i++) {
                m4 += pow(mean - current_sample[i], 4);
        }
        m2 =  m2 / ((float) size);
        m4 =  m4 / ((float) size);

        kurtosis = m4 / (m2 * m2) - 3.;

        return sqrt(m2);
}
