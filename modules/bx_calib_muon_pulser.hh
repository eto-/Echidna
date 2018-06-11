/*
Author: Timo Lewke <timo.lewke@ph.tum.de>, Michael Wurm <mwurm@ph.tum.de>
Maintainer: Timo lewke <timo.lewke@ph.tum.de>

Module for muon veto pulser, to distinguish good and bad channels
*/



#ifndef _BX_CALIB_MUON_PULSER_HH
#define _BX_CALIB_MUON_PULSER_HH
#include "bx_rec_general.hh"
#include "bx_base_module.hh"
#include "TROOT.h"
#include "TH2F.h"
#include "TH1F.h"

#include <map>

class bx_echidna_event;

class bx_calib_muon_pulser : public bx_base_module {
	public:
		bx_calib_muon_pulser ();
		virtual ~bx_calib_muon_pulser () {}

		virtual void begin ();
		virtual bx_echidna_event* doit (bx_echidna_event *ev);
		virtual void end ();

	private:
    // internal ev/hit ctrs
		int32_t i4_trg_ctr;
		int32_t a[256];
		std::vector<int32_t> v_nhits;
	
    // histograms for barn
		TH1F *pulser_hits_vs_lg;
		TH1F *working_channels;

};

#endif
