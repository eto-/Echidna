/* BOREXINO Reconstruction program
 *
 * Author: Alessandro Razeto <Alessandro.Razeto@ge.infn.it>
 * Maintainer: Alessandro Razeto <Alessandro.Razeto@ge.infn.it>
 *
 * $Id: db_acl.hh,v 1.13 2009/12/03 16:08:47 misiaszek Exp $
 *
 * The database table acls.
 * This class implement template meber function the same as cmap; this 
 * replication is seeded because of root which is not able to store
 * template map (like cmap). The map_get and map_check are shorthand to
 * preserve the db_* visitors functionality while being able to store them
 * in root files.
 *
 */
#ifndef _BD_TABLE_H
#define _BD_TABLE_H
#include <map>
#include <vector>
#include <sstream>
#include <string>
#include <stdexcept>


class bx_named;
class bx_dbi;

class db_acl {
  protected:
    db_acl () {}    // NO one will be allowed to create any db_* (except bx_dbi)
    ~db_acl () {}		// NO one will be allowed to destroy any db_* (except bx_dbi)
    void set_acl (const std::string& method, const std::string& bx_named_name) { acl[method].push_back(bx_named_name); }
    bool check_acl (const std::string& method, const bx_named* module);
  private:
//    db_acl (const db_acl& r)  {} // NO copy constructor for db_* classes
    typedef std::vector<std::string> bx_named_name_list;
    //std::map<std::string, bx_named_name_list> acl;
    std::map<std::string, std::vector<std::string> > acl;
  protected:
    template <typename key_t, typename value_t> const value_t& map_get (key_t k, const std::map<key_t, value_t>& map, const std::string& map_name) const { 
      typename std::map<key_t, value_t>::const_iterator item = map.find (k);
      if (item == map.end ()) {
	std::ostringstream msg;
	msg << "map: element " << k << " not found";
	if (map_name.size ()) msg << " in map \"" << map_name << "\"";
#ifdef CMAP_ABORT
	abort ();
#else
        throw std::runtime_error(msg.str ());
#endif	
      }
      return item->second;
    }

    template <typename key_t, typename value_t> bool map_check (key_t k, const std::map<key_t, value_t>& map) const { 
      return map.find(k) != map.end (); 
    }
    
    friend class bx_dbi;
};
#endif
/*  
 *  $Log: db_acl.hh,v $
 *  Revision 1.13  2009/12/03 16:08:47  misiaszek
 *  Changes for new streamer and rootcint in acl map
 *
 *  Revision 1.12  2007-09-14 04:00:59  razeto
 *  Allow copy ctor
 *
 *  Revision 1.11  2006-11-05 10:27:30  razeto
 *  Remove some useless include
 *
 *  Revision 1.10  2006/09/11 15:06:36  razeto
 *  Now the map name is required in map_get
 *
 *  Revision 1.9  2006/02/02 09:59:03  razeto
 *  Added some std:: qualifier
 *
 *  Revision 1.8  2006/01/25 13:23:05  misiaszek
 *  Added cmap functionality to db_acl
 *
 *  Revision 1.7  2004/11/26 15:25:10  razeto
 *  Thanks to Davide who so kindly underlined a syntactical error, I fixed the "Maintainer" word
 *
 *  Revision 1.6  2004/11/26 14:01:20  razeto
 *  Added Mantainer field
 *
 *  Revision 1.5  2004/11/24 13:05:03  razeto
 *  Upgraded to the new cmap ctor with name
 *
 *  Revision 1.4  2004/11/24 09:46:41  razeto
 *  Moved some prototyping in db_acl from db_*
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
