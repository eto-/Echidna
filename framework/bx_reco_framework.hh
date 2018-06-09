/* BOREXINO Reconstruction program
 *
 * Author: Alessandro Razeto <Alessandro.Razeto@ge.infn.it>
 * 	   Marco Pallavicini <pallas@ge.infn.it>
 * Maintainer: Alessandro Razeto <Alessandro.Razeto@ge.infn.it>
 *
 * $Id: bx_reco_framework.hh,v 1.17 2006/10/13 15:30:04 razeto Exp $
 *
 * Main reconstruction class
 * It handles the whole framework by getting a list of modules through the
 * bx_module_factory class and executing it sequentially
 * 
 */
#ifndef _BX_RECO_FRAMEWORK_H
#define _BX_RECO_FRAMEWORK_H

#include "bx_named.hh"
#include "bx_base_module.hh"
#include "bx_reader.hh"
#include "bx_rec_general.hh"

class bx_reco_framework: public bx_named {
  public:
    bx_reco_framework ();
    virtual ~bx_reco_framework ();

    virtual void run ();
    
  private:
    bx_base_module::bx_base_module_vector &modules;
    bx_reader *p_reader;
    int i4_max_events, i4_skip_events;
    int i4_precalib_events;
    void m_do_precalib (bx_base_module::module_role cycle);
};

#endif
/*
 * $Log: bx_reco_framework.hh,v $
 * Revision 1.17  2006/10/13 15:30:04  razeto
 * Added a flag
 *
 * Revision 1.16  2006-08-21 11:01:07  razeto
 * Added event skipping support + new precalib policy
 *
 * Revision 1.15  2006/04/02 11:19:41  razeto
 * Added skip first event option (-j)
 *
 * Revision 1.14  2005/03/04 11:51:11  razeto
 * Transformed the framework in a named to pass more options.
 * Added some options.
 *
 * Revision 1.13  2004/11/26 15:25:10  razeto
 * Thanks to Davide who so kindly underlined a syntactical error, I fixed the "Maintainer" word
 *
 * Revision 1.12  2004/11/26 14:01:19  razeto
 * Added Mantainer field
 *
 * Revision 1.11  2004/09/23 10:09:24  razeto
 * Changed internal reader name (to follow writer)
 *
 * Revision 1.10  2004/06/07 10:59:20  razeto
 * Upgraded to handle multiple writer at the same time
 *
 * Revision 1.9  2004/04/27 16:58:00  ddangelo
 * removed old includes
 *
 * Revision 1.8  2004/04/12 16:09:42  razeto
 * Added precalib loop handler. Added the loops for the precalibs
 *
 * Revision 1.7  2004/04/06 12:43:55  razeto
 * Added module role to base_module, and the code to handle everywhere. Some other minor fixes
 *
 * Revision 1.6  2004/04/01 12:08:53  razeto
 * Moved the base_module_operation operators in the bx_base_module class
 * scope; moved the bx_base_module_vector typedef in the bx_base_module
 * class.
 * Added a event delete to the frame.
 *
 * Revision 1.5  2004/03/24 14:23:35  razeto
 * Moved module list ownership from framewrok to module_factory, this allows:
 *   - a better syntax for bx_options which can store parameters even if the
 * framework is not created
 *   - inter module comunications
 * Since the module factory was already a singleton this does not affect
 * the modules array lifecycle.
 *
 * Revision 1.4  2004/03/23 13:01:05  razeto
 * Moved base_module_operatos to bx_base_module.hh.
 * Added compare_operation.
 * Added get_module by name to framework.
 *
 * Revision 1.3  2004/03/22 14:29:25  razeto
 * Some cosmetic changes
 *
 * Revision 1.2  2004/03/20 18:55:02  pallas
 * Debugging framework
 *
 * Revision 1.1  2004/03/20 17:39:45  pallas
 * First release
 *
 * Revision 1.1  2004/03/20 14:16:47  razeto
 * Added sources
 *
 *
 */
