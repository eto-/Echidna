/* BOREXINO Reconstruction program
 * 
 * Author: Alessandro Razeto <Alessandro.Razeto@ge.infn.it>
 * Maintainer: Alessandro Razeto <Alessandro.Razeto@ge.infn.it>
 *
 * $Id: barn_interface.hh,v 1.6 2008/10/01 16:55:32 razeto Exp $
 *
 * A repository for root objects, to be used from within modules to 
 * register histograms and other root objects (TNamed).
 * Registered root objects will be then written to the root file
 * According to a policy.
 *
 */

#ifndef _BARN_INTERFACE_H
#define _BARN_INTERFACE_H
#include "bx_named.hh"

#include <TROOT.h>
#include <TNamed.h>
#include <TFile.h>
#include <TDirectory.h>
#include <TObject.h>
#include <TObjString.h>
#include <list>
#include <string>

class root_item;
class bx_barn;
class TSocket;

class barn_interface: public bx_named {
  public:
    enum item_type {
      file,
      junk,
      net,
    };

    static barn_interface *get ();
    static void disable_root_signals ();

    void store (item_type type, TNamed *obj, const bx_named *owner);
    void write (TFile* file);

    void network_send (TObject *obj, const bx_named *owner);
  private:
    barn_interface ();
    static barn_interface* me;
    TDirectory* m_get_dir (TDirectory* cwd, const std::string& dir_name);
    TObject* m_get_obj (const root_item& item);
    void m_obj_write (TDirectory *dir, const root_item& item);
    void m_db_write  (TDirectory *dir, bx_barn *obj); 	 

    std::list<root_item> items; 
    TSocket *socket;  
    bool send;
    int retries;
};

class root_item {
  public:
    root_item (barn_interface::item_type type, TNamed *obj, const bx_named *owner);
    barn_interface::item_type get_type () const { return i_type; }
    const std::string& get_name () const { return name; }
    const std::string& get_owner_name () const { return owner_name; } 
    
  private:
    barn_interface::item_type i_type;
    std::string name;
    std::string owner_name;
};

#endif
/*
 * $Log: barn_interface.hh,v $
 * Revision 1.6  2008/10/01 16:55:32  razeto
 * Fixed a bug
 *
 * Revision 1.5  2008-10-01 16:38:52  razeto
 * Auto-reconnect to viewer
 *
 * Revision 1.4  2007-02-01 16:56:25  razeto
 * Added network dispatching to messenger
 *
 * Revision 1.3  2006/11/09 14:01:29  razeto
 * Fixed a warning
 *
 * Revision 1.2  2006/11/09 13:47:12  razeto
 * Added new messagge dispatcher for online echidna to gviewer
 *
 * Revision 1.1  2006/08/21 11:05:50  razeto
 * Moved bx_root_barn to barn_interface and db_barn to bx_barn
 *
 * Revision 1.12  2006/06/13 14:42:14  razeto
 * Removed duoble qualification (to compile with gcc 4.1)
 *
 * Revision 1.11  2006/05/24 15:11:26  razeto
 * Upgraded to send owner name to network server
 *
 * Revision 1.10  2006/05/17 14:59:40  razeto
 * Added experimental interface sending object to network
 *
 * Revision 1.9  2006/01/02 21:21:45  razeto
 * Removed test target, renamed barn folder, moved and renamed db_barn
 *
 * Revision 1.8  2005/11/04 00:31:25  misiaszek
 *
 * ----------------------------------------------------------------------
 * Modified Files:
 *  	bx_root_barn.hh
 *
 * db_barn object & m_db_write function added
 * ----------------------------------------------------------------------
 *
 * Revision 1.7  2005/02/25 16:46:40  razeto
 * Upgraded to store histos even from bx_named
 *
 * Revision 1.6  2004/11/26 15:25:10  razeto
 * Thanks to Davide who so kindly underlined a syntactical error, I fixed the "Maintainer" word
 *
 * Revision 1.5  2004/11/26 14:01:20  razeto
 * Added Mantainer field
 *
 * Revision 1.4  2004/08/06 10:30:28  razeto
 * cycle_1 branch merged in the main trunk.
 *
 * Revision 1.3.2.1  2004/07/29 08:53:25  razeto
 * Added junk object type
 *
 * Revision 1.3  2004/05/21 11:34:48  razeto
 * Disabled root signal handlers
 *
 * Revision 1.2  2004/05/21 08:34:52  razeto
 * Updated to working
 *
 * Revision 1.1  2004/05/18 15:07:01  razeto
 * Added
 *
 */
