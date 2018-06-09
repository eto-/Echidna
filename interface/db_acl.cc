/* BOREXINO Reconstruction program
 *
 * Author: Alessandro Razeto <Alessandro.Razeto@ge.infn.it>
 * Maintainer: Alessandro Razeto <Alessandro.Razeto@ge.infn.it>
 *
 * $Id: db_acl.cc,v 1.7 2006/02/02 10:38:55 razeto Exp $
 *
 * The database table acl implementation.
 *
 */
#include "db_acl.hh"
#include "bx_dbi.hh"
#include "bx_named.hh"
#include <algorithm>

#ifndef _ECHIDNA_ROOTLIB_
bool db_acl::check_acl (const std::string& method, const bx_named* obj) {
  if (map_check (method, acl)) {
    const bx_named_name_list& list = acl[method];
    if (std::find (list.begin (), list.end (), obj->get_name ()) != list.end ()) 
      return true;
  }
  bx_dbi::get ()->get_message (bx_message::error) << "bx_table: acl failed for bx_named " << obj->get_name () << " with " << method << dispatch;
  return false;
}
#else
bool db_acl::check_acl (const std::string& method, const bx_named* obj) { return false; }
#endif
/*  
 *  $Log: db_acl.cc,v $
 *  Revision 1.7  2006/02/02 10:38:55  razeto
 *  Fixed unresolved symbol of db_acl in bx_phys
 *
 *  Revision 1.6  2006/01/25 13:23:05  misiaszek
 *  Added cmap functionality to db_acl
 *
 *  Revision 1.5  2004/11/26 15:25:10  razeto
 *  Thanks to Davide who so kindly underlined a syntactical error, I fixed the "Maintainer" word
 *
 *  Revision 1.4  2004/11/26 14:01:20  razeto
 *  Added Mantainer field
 *
 *  Revision 1.3  2004/08/14 21:18:43  razeto
 *  Upgraded to use bx_named instead of only modules
 *
 *  Revision 1.2  2004/05/31 11:35:30  razeto
 *  Removed a dependancy
 *
 *  Revision 1.1  2004/04/26 13:48:29  razeto
 *  Added db_acl to check set calls for calling module to have the right privileges.
 *  Fixed the names of set/get methods.
 *  Modifications decided at the software Paris meeting.
 *
 */
