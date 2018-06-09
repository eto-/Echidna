/*
 * Author: Juergen Winter <juergen.winter@ph.tum.de>
 * Maintainer: Juergen Winter <juergen.winter@ph.tum.de>
 *
 * $Id: bx_calib_muon_alignment.hh,v 1.5 2008/09/26 14:53:49 winter Exp $
 *
 * Module for tagging misalignment between muon data and mcr
*/


#ifndef _BX_CALIB_MUON_ALIGNMENT_HH
#define _BX_CALIB_MUON_ALIGNEMNT_HH
#include "bx_rec_general.hh"
#include "bx_base_module.hh"
#include "TROOT.h"
#include "TH2F.h"
#include "TH1F.h"

#include <map>

class bx_echidna_event;

class bx_calib_muon_alignment : public bx_base_module {
	public:
		bx_calib_muon_alignment ();
		virtual ~bx_calib_muon_alignment () {}

		virtual void begin ();
		virtual bx_echidna_event* doit (bx_echidna_event *ev);
		virtual void end ();

	private:
    	// internal ev ctr
		int num;
		bool b_check_last;
		bool b_check_mcr;
		
		unsigned long u4_evnum;
		time_t time;	
		unsigned long u4_evnum_up;
		time_t time_up;
		unsigned long u4_evnum_aligned;
		time_t time_aligned;
	
	// to read from echidna.cfg
		int nmis_thresh;
		int nhits_thresh;
};

#endif

/*
 * $Log: bx_calib_muon_alignment.hh,v $
 * Revision 1.5  2008/09/26 14:53:49  winter
 * changed time format
 *
 * Revision 1.4  2008-09-11 16:18:23  winter
 * write to db, check mcr crash
 *
 * Revision 1.3  2008-08-13 14:50:42  winter
 * completed
 *
 * Revision 1.2  2008-05-13 12:42:55  ddangelo
 * compiles
 *
 * Revision 1.1  2008-05-13 12:33:13  ddangelo
 * added
 *
 */
