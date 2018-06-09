/* BOREXINO Reconstruction program
 *
 * Author: Alessandro Razeto <Alessandro.Razeto@ge.infn.it>
 * 	   Marco Pallavicini <pallas@ge.infn.it>
 * Maintainer: Alessandro Razeto <Alessandro.Razeto@ge.infn.it>
 *
 * $Id: bx_reco_framework.cc,v 1.38 2008/08/11 12:39:26 ddangelo Exp $
 *
 * Implementation of bx_reco_framework
 *
 */
#include "bx_reco_framework.hh"
#include "bx_module_factory.hh"
#include "messenger.hh"
#include "bx_dbi.hh"
#include "db_run.hh"

#include <iostream>
#include <algorithm>

// ctor
bx_reco_framework::bx_reco_framework(): bx_named ("bx_reco_framework"), modules(bx_module_factory::get ()->get_modules ()) {
  p_reader = dynamic_cast<bx_reader *>(bx_module_factory::get ()->get_module (bx_base_module::reader));

  i4_max_events = get_parameter ("max_events").get_int ();
  i4_skip_events = get_parameter ("skip_events").get_int ();
  i4_precalib_events = get_parameter ("precalib_events").get_int ();
}

// dtor
bx_reco_framework::~bx_reco_framework () {
  // modules has not to be deleted since are owned by bx_module_factory
  bx_module_factory::get ()->delete_module ();
}

// execute all program; basically, this is the Borexino Reco Main
void bx_reco_framework::run () {
  
  // init reader module 
  p_reader->begin ();

    // Do precalibration (4 cycles), only if not yet present on db
  if (get_parameter ("do_precalib").get_bool ()) {
    if (bx_dbi::get ()->get_run ().is_precalib_present () && !get_parameter ("redo_precalib").get_bool ())
      get_message (bx_message::warn) << "run already precalibrated, set bx_reco_framework.redo_precalib to 1 to rerun them" << dispatch;
    else {
      bx_dbi::get ()->get_run ().reset_precalibrations (this);
      m_do_precalib (bx_base_module::precalib_cycle1);
      m_do_precalib (bx_base_module::precalib_cycle2);
      m_do_precalib (bx_base_module::precalib_cycle3);
      m_do_precalib (bx_base_module::precalib_cycle4);
    }
  } else if (!bx_dbi::get ()->get_run ().is_precalib_present () && !get_parameter ("ignore_precalib_absence").get_bool ())
    get_message (bx_message::critic) << "run not yet precalibrated, run echidna with \"precalibrations\" configuration first" << dispatch;

    // Check for electronics calibration
  if (get_parameter ("require_electronics_calibrations").get_bool () && !bx_dbi::get ()->get_run ().is_laben_electronic_channel_present ())
    get_message (bx_message::critic) << "run not yet electronics calibrated, run echidna with \"electronics_calibrations\" configuration first" << dispatch;

    // Init main loop modules and writer
  std::for_each (modules.begin (), modules.end (), bx_base_module::begin_operation (bx_base_module::main_loop));
  std::for_each (modules.begin (), modules.end (), bx_base_module::begin_operation (bx_base_module::writer));

    // Skip first events
  bx_echidna_event *ev_skip;
  while (i4_skip_events-- > 0 && (ev_skip = p_reader->doit (0))) delete ev_skip;

    // Do main loop
  long int ev_count = 0;
  while (bx_echidna_event *ev = p_reader->doit (0)) {
    // check here for maximum number of events; to be done
    if (i4_max_events >= 0 && ev_count++ >= i4_max_events) break;

    try {
      // loop on modules and filter event
      std::for_each (modules.begin (), modules.end (), bx_base_module::doit_operation (ev, bx_base_module::main_loop));
  
      // write event (root file or other)
      std::for_each (modules.begin (), modules.end (), bx_base_module::doit_operation (ev, bx_base_module::writer));
    } catch (bx_base_module::skip_event &s) {
      get_message (bx_message::debug) << "skipped event by " << s.name << dispatch;
    }

    delete ev;
  }

  // end of job
  p_reader->end ();
  std::for_each (modules.begin (), modules.end (), bx_base_module::end_operation (bx_base_module::main_loop));
  std::for_each (modules.begin (), modules.end (), bx_base_module::end_operation (bx_base_module::writer));

  // write visitors data to database
  bx_dbi::get ()->flush_visitors (this);
}

void bx_reco_framework::m_do_precalib (bx_base_module::module_role cycle) {
  std::for_each (modules.begin (), modules.end (), bx_base_module::begin_operation (cycle));
  for (int i = 0; i < i4_precalib_events; i++) {
    bx_echidna_event *ev = p_reader->doit (0);
    if (!ev) break;
    try {
      std::for_each (modules.begin (), modules.end (), bx_base_module::doit_operation (ev, cycle));
    } catch (bx_base_module::skip_event &s) {
      get_message (bx_message::debug) << "skipped event by " << s.name << dispatch;
    }
    delete ev;
  }
  p_reader->rewind ();
  std::for_each (modules.begin (), modules.end (), bx_base_module::end_operation (cycle));
}

/*
 * $Log: bx_reco_framework.cc,v $
 * Revision 1.38  2008/08/11 12:39:26  ddangelo
 * complying to a modified variable name
 *
 * Revision 1.37  2007-11-26 11:08:00  razeto
 * Now echidna will not run unless electronics_calibrations data are present (required by livia)
 *
 * Revision 1.36  2007-01-30 11:03:07  razeto
 * When doing precalib alway start from scratch
 *
 * Revision 1.35  2006/10/13 15:30:04  razeto
 * Added a flag
 *
 * Revision 1.34  2006-08-21 11:01:07  razeto
 * Added event skipping support + new precalib policy
 *
 * Revision 1.33  2006/05/10 12:14:39  razeto
 * Removed some useless include
 *
 * Revision 1.32  2006/04/02 11:28:49  razeto
 * Fixed a small leak in skipping events
 *
 * Revision 1.31  2006/04/02 11:19:41  razeto
 * Added skip first event option (-j)
 *
 * Revision 1.30  2005/12/30 11:37:47  razeto
 * Added a log message
 *
 * Revision 1.29  2005/03/04 11:51:11  razeto
 * Transformed the framework in a named to pass more options.
 * Added some options.
 *
 * Revision 1.28  2005/03/02 15:45:43  razeto
 * Fixed a test (buggy for older runs) from framework to precalib modules
 *
 * Revision 1.27  2004/11/26 15:25:10  razeto
 * Thanks to Davide who so kindly underlined a syntactical error, I fixed the "Maintainer" word
 *
 * Revision 1.26  2004/11/26 14:01:19  razeto
 * Added Mantainer field
 *
 * Revision 1.25  2004/10/25 02:48:37  razeto
 * Added a better "end of echidna" message diplaying the number of errors and warns
 *
 * Revision 1.24  2004/10/19 18:41:58  razeto
 * Integrated visitors writing in the framework
 *
 * Revision 1.23  2004/09/30 14:36:12  razeto
 * Added a dctor which delete modules
 *
 * Revision 1.22  2004/09/23 10:09:24  razeto
 * Changed internal reader name (to follow writer)
 *
 * Revision 1.21  2004/09/22 14:04:13  razeto
 * Changed bx_reco_event to bx_echidna_event
 *
 * Revision 1.20  2004/06/07 10:59:20  razeto
 * Upgraded to handle multiple writer at the same time
 *
 * Revision 1.19  2004/05/26 09:24:53  razeto
 * Upgraded to check if precalib are needed
 *
 * Revision 1.18  2004/05/26 08:27:26  razeto
 * Added max number of events (negative = unlimited)
 *
 * Revision 1.17  2004/05/21 08:41:15  razeto
 * Fixed the trgtype code; still to be removed
 *
 * Revision 1.16  2004/05/19 08:43:37  razeto
 * Updated to follow new trigger raw event syntax
 *
 * Revision 1.15  2004/05/18 14:54:52  razeto
 * Added a simple test, to be removed
 *
 * Revision 1.14  2004/04/12 16:09:42  razeto
 * Added precalib loop handler. Added the loops for the precalibs
 *
 * Revision 1.13  2004/04/06 12:43:55  razeto
 * Added module role to base_module, and the code to handle everywhere. Some other minor fixes
 *
 * Revision 1.12  2004/04/03 09:24:42  razeto
 * Added messenger
 *
 * Revision 1.11  2004/04/01 12:08:53  razeto
 * Moved the base_module_operation operators in the bx_base_module class
 * scope; moved the bx_base_module_vector typedef in the bx_base_module
 * class.
 * Added a event delete to the frame.
 *
 * Revision 1.10  2004/03/24 16:22:08  razeto
 * Moved reader,writer,precalib creation from framework to factory
 *
 * Revision 1.9  2004/03/24 14:23:35  razeto
 * Moved module list ownership from framewrok to module_factory, this allows:
 *   - a better syntax for bx_options which can store parameters even if the
 * framework is not created
 *   - inter module comunications
 * Since the module factory was already a singleton this does not affect
 * the modules array lifecycle.
 *
 * Revision 1.8  2004/03/23 13:01:05  razeto
 * Moved base_module_operatos to bx_base_module.hh.
 * Added compare_operation.
 * Added get_module by name to framework.
 *
 * Revision 1.7  2004/03/22 14:50:08  razeto
 * Introduced algorithms for looping on modules with predefined operations
 *
 * Revision 1.6  2004/03/22 14:29:25  razeto
 * Some cosmetic changes
 *
 * Revision 1.5  2004/03/20 19:05:17  pallas
 * Debugging
 *
 * Revision 1.4  2004/03/20 18:55:02  pallas
 * Debugging framework
 *
 * Revision 1.3  2004/03/20 18:46:45  pallas
 * debugging
 *
 * Revision 1.2  2004/03/20 17:50:54  pallas
 * Debugging
 *
 * Revision 1.1  2004/03/20 17:39:45  pallas
 * First release
 *
 * Revision 1.1  2004/03/20 14:16:47  razeto
 * Added sources
 *
 *
 */
