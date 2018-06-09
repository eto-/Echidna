/* BOREXINO Reconstruction program
 *
 * Author: Alessandro Razeto <Alessandro.Razeto@ge.infn.it>
 * 	   Marco Pallavicini <pallas@ge.infn.it>
 * Maintainer: Alessandro Razeto <Alessandro.Razeto@ge.infn.it>
 *
 * $Id: bx_module_factory.hh,v 1.13 2004/11/26 15:25:10 razeto Exp $
 *
 * Factory used to build all modules
 * Singleton
 * It is the only place in the program where modules are created
 */
#ifndef _BX_MODULE_FACTORY_H
#define _BX_MODULE_FACTORY_H

#include "bx_rec_general.hh"
#include "bx_base_module.hh"

#include <string>

class bx_module_factory {
  public:
    static bx_module_factory *get ();
    void create_modules (const std::vector<std::string>& module_name_list);

    // Return the module list (for the framework)
    bx_base_module::bx_base_module_vector &get_modules () { return modules; }
    
    bx_base_module *get_module (const std::string& module_name);
    bx_base_module *get_module (bx_base_module::module_role role);
    
    void delete_module (const std::string& module_name = "ALL");
  private:
    bx_module_factory ();
    ~bx_module_factory ();
    static bx_module_factory *me;
    
    bx_base_module::bx_base_module_vector modules;
    void m_check_single_presence (bx_base_module::module_role role);  	// check if there is just one active
    								      	// module with the indicated role.
};

#endif
/*
 * $Log: bx_module_factory.hh,v $
 * Revision 1.13  2004/11/26 15:25:10  razeto
 * Thanks to Davide who so kindly underlined a syntactical error, I fixed the "Maintainer" word
 *
 * Revision 1.12  2004/11/26 14:01:19  razeto
 * Added Mantainer field
 *
 * Revision 1.11  2004/09/30 14:35:36  razeto
 * Added delete_module to module_factory to destroy one or all modules
 *
 * Revision 1.10  2004/07/22 10:13:01  razeto
 * Jumbo patch for the configuration manager applied, read the documentation for usage
 *
 * Revision 1.9  2004/05/26 08:15:41  razeto
 * Removed unused headers (Davide)
 *
 * Revision 1.8  2004/04/27 16:58:00  ddangelo
 * removed old includes
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
 * Revision 1.5  2004/03/26 16:31:14  razeto
 * Fixed a bug: reader, writer and precalibrator were pointers, but the could
 * be even not initialized since the initialization were left to a procedure
 * which would allocate them if they where present on the modules.cfs file.
 * This could not be true.
 * Now these special modules are in a dedicated vector which however has
 * the same paradigm of the standard module vector.
 * This solves the problem; maybe a role for the module could be introduced
 * in future, allowing to have just one vector.
 *
 * Revision 1.4  2004/03/24 16:22:08  razeto
 * Moved reader,writer,precalib creation from framework to factory
 *
 * Revision 1.3  2004/03/24 14:23:35  razeto
 * Moved module list ownership from framewrok to module_factory, this allows:
 *   - a better syntax for bx_options which can store parameters even if the
 * framework is not created
 *   - inter module comunications
 * Since the module factory was already a singleton this does not affect
 * the modules array lifecycle.
 *
 * Revision 1.2  2004/03/22 14:29:25  razeto
 * Some cosmetic changes
 *
 * Revision 1.1  2004/03/20 17:39:45  pallas
 * First release
 *
 * Revision 1.1  2004/03/20 14:16:47  razeto
 * Added sources
 *
 *
 */
