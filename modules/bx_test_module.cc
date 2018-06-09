/* BOREXINO Reconstruction program
 *
 * Author: Alessandro Razeto <Alessandro.Razeto@ge.infn.it>
 * 	   Marco Pallavicini <pallas@ge.infn.it>
 * Maintainer: Davide D'Angelo <Davide.Dangelo@lngs.infn.it>
 *
 * $Id: bx_test_module.cc,v 1.16 2013/01/29 06:15:37 mosteiro Exp $
 *
 * Implementation of bx_test_module
 * If you cut and paste from this file, 
 * remember to remove comments
 * 
 * A few lines of sample code are commented out
 * to allow clean compilation.
 * C-style comments are used in this cases, 
 * while C++style are used for explanations
 */

#include "bx_test_module.hh"
#include "bx_echidna_event.hh"
#include "barn_interface.hh"
#include "db_run.hh"
#include "db_profile.hh"
#include "bx_dbi.hh"
#include "TH1F.h"
#include "TH2F.h"

// ctor
bx_test_module::bx_test_module (): bx_base_module("bx_test_module", bx_base_module::main_loop) {
  // Option: you can require to run only on events which have reached a given stage (for a given detector segment)
  require_event_stage (bx_detector::laben, bx_base_event::clustered);
  // Option: you can require to run only on events with a given trigger type
  require_trigger_type (bx_trigger_event::neutrino);
}


// module interface
void bx_test_module::begin () {
  // Sample message
  get_message (bx_message::debug) << "test begin" << dispatch;

  // resource initialization
  i4_times = 0;
  my_vector = new std::vector<double>(1000);

  // Option: histogram to test how the module work, to be written to file
  my_histo_check = new TH1F ("my_histo_check", "My useful plot of this and that", 255, 0, 255);
  barn_interface::get ()->store (barn_interface::file, my_histo_check, this);

  // Option: histogram for internal calculations, not to be written to file
  my_histo_compute = new TH2F ("my_histo_compute", "My histo to compute my stuff", constants::laben::channels, 1, constants::laben::channels+1, 255, 0, 255);
  barn_interface::get ()->store (barn_interface::junk, my_histo_compute, this);
}


bx_echidna_event* bx_test_module::doit (bx_echidna_event *ev) {
  // Get a const ptr to read from (or more)
  const bx_laben_event& er = ev->get_laben ();

  // Get a non-const ptr to write to (calib and precalib modules won't need it)
  /*  bx_laben_decoded_event& ew = dynamic_cast<bx_laben_decoded_event&>(ev->get_laben ()); */

  // Check if data are present, otherwise return
  if (!er.get_decoded_nhits()) 
    return ev;

  // Precalib and calib modules only: say you have data to enable end() method computation
  b_has_data = true;

  // Get const references to database visitors with 1 or more keys (run, profile, calib_profile)
  const db_run& run_info = bx_dbi::get ()->get_run ();
  const db_profile& profile_info = bx_dbi::get ()->get_profile ();

  // Option: call a private method to delegate some computation
  bool my_result = m_check_this_and_that (er);

  // Sending a message with "critic" level ends the program. 
  // The msg handler does this for you, 
  // don't call exit(), abort() or similar functions by your own
  if (my_result != true) 
    get_message (bx_message::critic) << "Fatal error, check not passed!" << dispatch;
    
  // Do semthing ...
  if ( (i4_times++ % 100) == 0)
    if (i4_times!=1) 
      get_message (bx_message::debug) << "test doit called " << i4_times-1 << " times" << dispatch;

  // Example: loop on hits
  for (int i=0 ; i<er.get_decoded_nhits(); i++) {    
    // get a const reference to the hit you want to process
    const bx_laben_decoded_hit& h = er.get_decoded_hit(i);
    int lg = h.get_raw_hit().get_logical_channel();

    // Option: you can select only channels of a given type 
    if (profile_info.logical_channel_description(lg) != db_profile::ordinary)
      continue;

    // Get some calib or precalib info from database visitor(s).
    double pedestal = run_info.get_muon_precalib_pedestal(lg);

    // Get a user parameter
    unsigned long my_sample_parameter = (unsigned long)(get_parameter("my_sample_parameter").get_int());

    // Compute something ...
    double my_computed_value = (h.get_raw_time()-pedestal)*my_sample_parameter;
    
    // Fill histograms
    my_histo_compute->Fill(lg, my_computed_value);
  } 

  // Write to event level (calib and precalib modules won't need it)
  // Contact event maintainer to know wehre to write
  /* ew.my_variable = my_histo_compute->GetMean(); */

  // Event based modules: mark achieved event stage
  /*  er.mark_stage (bx_base_event::my_stage); */
  return ev;
}

void bx_test_module::end () {
  get_message (bx_message::debug) << "test end" << dispatch;

  // Calib and precalib only: work only if data were present
  if(b_has_data) {
    
    // Get a non-const reference to a db visitor to write to
    /*    db_run& run_info = bx_dbi::get ()->get_run (); */

    // Example: Loop on channels
    for(int i=0; i<constants::muon::channels; i++) { 
      // compute your stuff...

      // Save it to db; contact db interface maintainer for appropriate setters.
      /* run_info.set_my_variable(channel, value, this); */

      // Fill a test histo
      my_histo_check->Fill(5.);
    }

  }

  // If user wanted it, ask db interface to write to db:
  // Policies on db writing have to be agreed in working group
  if (get_parameter ("write_my_calib").get_bool ()) 
    /* run_info.write_my_calib (true, this) */;

  // Deallocate resourses (but not the ROOT objects!!!)
  delete my_vector;
}

// Private method(s)
bool bx_test_module::m_check_this_and_that(const bx_laben_event& er) {
  // ...
  return true;
}

/*
 * $Log: bx_test_module.cc,v $
 * Revision 1.16  2013/01/29 06:15:37  mosteiro
 * channel count starts at 1, not 0
 *
 * Revision 1.15  2006-08-21 11:19:02  razeto
 * Updated to new barn_interface
 *
 * Revision 1.14  2006/01/09 16:21:03  razeto
 * Updated to the new root_barn target
 *
 * Revision 1.13  2005/02/10 17:37:39  ddangelo
 * removed afew compilation warning
 *
 * Revision 1.12  2005/02/03 19:00:18  ddangelo
 * maintainer changed (Razeto->D'Angelo)
 * Module really implemented with many possible examples for developers.
 * It features all aspects described in "Module's programmers guide" in the docs.
 *
 * Revision 1.11  2004/11/26 15:25:11  razeto
 * Thanks to Davide who so kindly underlined a syntactical error, I fixed the "Maintainer" word
 *
 * Revision 1.10  2004/11/26 14:01:21  razeto
 * Added Mantainer field
 *
 * Revision 1.9  2004/09/22 13:28:37  razeto
 * Changed bx_reco_event to bx_echidna_event
 *
 * Revision 1.8  2004/04/06 12:42:17  razeto
 * Added module role to base_module, and the code to handle everywhere. Some other minor fixes
 *
 * Revision 1.7  2004/04/03 09:20:35  razeto
 * Added messenger
 *
 * Revision 1.6  2004/04/02 13:40:01  razeto
 * Added messenger
 *
 * Revision 1.5  2004/03/21 18:55:05  razeto
 * Some cosmetic changes
 *
 * Revision 1.4  2004/03/20 19:10:03  pallas
 * Debugging
 *
 * Revision 1.3  2004/03/20 19:05:17  pallas
 * Debugging
 *
 * Revision 1.2  2004/03/20 18:56:46  pallas
 * Change counter
 *
 * Revision 1.1  2004/03/20 17:39:45  pallas
 * First release
 *
 * Revision 1.1  2004/03/20 14:16:47  razeto
 * Added sources
 *
 */
