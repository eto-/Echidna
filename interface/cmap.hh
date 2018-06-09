/* BOREXINO Reconstruction program
 *
 * Author: Alessandro Razeto <Alessandro.Razeto@ge.infn.it>
 * Maintainer: Alessandro Razeto <Alessandro.Razeto@ge.infn.it>
 *
 * $Id: cmap.hh,v 1.12 2006/06/13 14:42:51 razeto Exp $
 *
 * In this file a new template is defined: cmap is a surrogate of map.
 * The advantages of cmap is that the [] operator is defined const
 * and can be used in a const environment (like a const reference to cmap).
 * On the contrary if [] operator does not find the required element an
 * exception (runtime_error) is generated.
 * !!! To use cmap the key has to be capable to be pushed on a ostream
 * (a ostream& operator<< (ostream&, const key &) has been defined);
 * most builtin types (numbers. string, and so on) have this capability.
 * This class is diccult to read and to understand, but the syntax is 
 * defined by STL, since cmap inherits from map :-(
 * cmap is in the std namespace to be alligned with every container 
 * (map, vector, ...).
 * No bx_message is used since cmap is intended to trow errors like
 * other standard containers.
 *
 */
#ifndef _CMAP_H
#define _CMAP_H

#include <map>
#include <stdexcept>
#include <utility>
#include <memory>
#include <sstream>
#include <string>

namespace std {
  template <typename key, typename value, 
    typename compare_operation = less<key>, 
    typename alloc = allocator<pair<const key, value> > >	// >>> is not handled in c++ for templates
  								// since >> is an operator
  class cmap: public map<key, value, compare_operation, alloc> {
    private:
      typedef map<key, value, compare_operation, alloc> _map;	// these definitions are shortcut for the
      typedef cmap<key, value, compare_operation, alloc> _cmap; // following code
      std::string map_name;
    public:
      typedef typename _map::const_iterator const_iterator;
      typedef typename _map::iterator iterator;
        // Some default ctors
      cmap (const std::string& name="unknown"): _map(), map_name(name) {}
      cmap (const _map& x, const std::string& name = "unknown"): _map(x), map_name(name) {}
      cmap (const _cmap& x): _map(x), map_name(x.map_name) {}
        // Check if the element with key is present. 
      bool check (const key& k) const { return _map::find (k) != _map::end (); }
        // Operators[] : 2 version are present non-const and const
	// Both are needed since then the compiler can choose the one which
	// fits the usage, The non-const operator[] has to be rewritten since 
	// the one inherited from map always has lower precedence than the 
	// const operator[] (and defacto is never used as is). A simple wrapper
	// restabilishes the precedence between the 2 [] operators.
      value& operator[] (const key& k) { return (dynamic_cast<_map &>(*this))[k]; }
      const value& operator[] (const key& k) const {
        const_iterator item = _map::find (k);

        if (item == _map::end ()) {
	  ostringstream msg;
	  msg << "cmap: element " << k << " not found";
	  if (map_name.size ()) msg << " in map \"" << map_name << "\"";
#ifdef CMAP_ABORT
	  abort ();
#else
	  throw runtime_error(msg.str ());
#endif
	}
        return item->second;
      }
        // rfind is an utility to find the first element whose value is equal to the
	// given test value: it's a find on the values instead of the keys.
	// Since rfind can find more than one element the iterator to the first
	// found element is returned; the search can be restarted passing that
	// iterator + 1 back to rfind as second element.
	// Two version, const and non-const are present.
	// To uset rfind the == operator has to be meaningfull for the value type
	// (all builtin type have support for it).
	// A version with 1 argument is needed since "start = begin ()" is not
	// allowed (since begin () is not a value).
      const_iterator rfind (const value& v, const_iterator start) const {
	const_iterator item = start;
	for (; item != _map::end (); item++) if (item->second == v) break;
	return item;
      }
      const_iterator rfind (const value& v) const { return rfind (v, _map::begin ()); }
      iterator rfind (const value& v, iterator start) {
	iterator item = start;
	for (; item != _map::end (); item++) if (item->second == v) break;
	return item;
      }
      iterator rfind (const value& v) { return rfind (v, _map::begin ()); }
   };
}
        




#endif
/*
 * $Log: cmap.hh,v $
 * Revision 1.12  2006/06/13 14:42:51  razeto
 * Added template specification (to compile with gcc 4.1)
 *
 * Revision 1.11  2004/11/30 13:20:21  razeto
 * Added an ifdef for debugging
 *
 * Revision 1.10  2004/11/26 15:25:10  razeto
 * Thanks to Davide who so kindly underlined a syntactical error, I fixed the "Maintainer" word
 *
 * Revision 1.9  2004/11/26 14:01:20  razeto
 * Added Mantainer field
 *
 * Revision 1.8  2004/11/24 13:03:03  razeto
 * Added name to cmap
 *
 * Revision 1.7  2004/08/06 10:30:28  razeto
 * cycle_1 branch merged in the main trunk.
 *
 * Revision 1.6.4.1  2004/07/23 16:07:22  razeto
 * Added some template qualification to compile even with g++ 3.4
 *
 * Revision 1.6  2004/04/13 12:07:49  razeto
 * Fixed a print
 *
 * Revision 1.5  2004/04/05 13:11:06  razeto
 * Updated a comment
 *
 * Revision 1.4  2004/04/02 13:35:03  razeto
 * Added some comments. Added typedef for iterators. Added rfind methods
 *
 * Revision 1.3  2004/04/01 12:01:47  razeto
 * Fixed 3 bugs:
 *   - nonconst [] operator added
 *   - equality test using = (instea of ==)
 *   - returnign the iterator value with the wrong syntax
 * Added a check method
 *
 * Revision 1.2  2004/03/31 09:40:58  razeto
 * Added comments
 *
 * Revision 1.1  2004/03/31 09:34:20  razeto
 * Added to the repository
 *
 */
