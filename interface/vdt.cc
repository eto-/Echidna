/* BOREXINO Reconstruction program
 *
 * Author: Alessandro Razeto <Alessandro.Razeto@ge.infn.it>
 * Maintainer: Alessandro Razeto <Alessandro.Razeto@ge.infn.it>
 *
 * $Id: vdt.cc,v 1.20 2008/09/24 12:13:08 razeto Exp $
 *
 * Implementation of vdt
 *
 */
#include "vdt.hh"
#include "messenger.hh"

#include <sstream>
#include <iterator>
#include <ctype.h>
#include <stdlib.h>

vdt::vdt_type vdt::m_search_type (const std::string& str) {
    // empty string
  if (!str.size ()) return string_vdt;

    // vector first since are simple
  if (str[0] == '{' && str[str.size () - 1] == '}') return vector_vdt;

    // true and false as litteral are integers
  if (str == "true" || str == "TRUE" || str == "false" || str == "FALSE") return int_vdt;

    // then check for number (int or float)
  char *error_ptr = 0;
  ::strtod (str.c_str (), &error_ptr);
  if (error_ptr && *error_ptr) return string_vdt;  // not a number then is a string

    // then check for float
  *error_ptr = 0;
  ::strtol (str.c_str (), &error_ptr, 0);
  if (error_ptr && *error_ptr) return float_vdt;  // not integer then is a float

    // last is integer
  return int_vdt;
}

void vdt::m_assign_value (const std::string& str) {
  bx_message msg(bx_message::critic, "vdt: ");
  if (name_s.size ()) msg << "\"" << name_s << "\" ";

  char *error_ptr;
  value_s = str; // ALWAYS save original string
  switch (type_internal) {
    case vdt::float_vdt:
      value_f = ::strtod (str.c_str (), &error_ptr);
      if (value_f == 0 && error_ptr && *error_ptr) msg << "float conversion failed for " << str << dispatch;
      break;
    case vdt::int_vdt:
      if (str == "true" || str == "TRUE") value_i = 1;
      else if (str == "false" || str == "FALSE") value_i = 0;
      else {
	value_i = ::strtol (str.c_str (), &error_ptr, 0);  // 0 means autodetect the base, which is fine for 0x starting hex
	if (value_i == 0 && error_ptr && *error_ptr) msg << "integer conversion failed for \"" << str << "\"" << dispatch;
      }
      break;
    case vdt::vector_vdt:
      do {
	  // Fill the vector
        std::string sub_str = str.substr (1, str.size () - 2);  	// copy to local (modificable) discarding '{' and '}'
	std::string::size_type colon = 0;
        while (sub_str.size ()) {
          colon = sub_str.find (',', 0);			  	// Find first ',' (since string is erased from previous fields)

          if (colon == std::string::npos) colon = sub_str.size (); 	// If no more ',' are present sub_str is the last field

          std::string::size_type start_field = 0;			// look for spaces around the element like " 12 ,"
	  for (; start_field < colon; start_field++) if (!isspace (sub_str[start_field])) break;
          std::string::size_type end_field = colon - 1;
	  for (; end_field > start_field; end_field--) if (!isspace (sub_str[end_field])) break;

	  std::string value = sub_str.substr (start_field, end_field - start_field + 1); // Copy the value

	  if (!value.size ()) break;					// Stop on last void element

	  value_v.push_back (vdt (value, name_s));			// Fill the vector

	  sub_str.erase (0, colon + 1);  				// Move forward
        }	
      } while (0); // This do while is present to allow the declaration of variables
      break;
    case vdt::string_vdt:  // nothing to do for strings
      break;
    default:
      msg << "internal vdt error: cannot call m_assign_value with a undefined type value" << dispatch;
  }
} 


vdt::vdt (const std::string& value, const std::string& name, vdt_type type): name_s(name) {
  vdt_type tmp_type = m_search_type (value);

  type_internal = tmp_type;
  if (type != unassigned_vdt) {
    m_check_type (type);
    type_internal = type;
  }

  m_assign_value (value);
}

const vdt& vdt::operator= (const std::string& value) {
  vdt_type tmp_type = m_search_type (value);

  if (type_internal != unassigned_vdt) m_check_type (tmp_type);

  type_internal = tmp_type;
  m_assign_value (value);

  return *this;
}

const vdt& vdt::operator= (long int value) {
  
  if (type_internal != unassigned_vdt) m_check_type (vdt::int_vdt);

  type_internal = int_vdt; 
  value_i = value; 
  
  return *this; 
}
  
const vdt& vdt::operator= (double value) {
  
  if (type_internal != unassigned_vdt) m_check_type (vdt::float_vdt);

  type_internal = float_vdt; 
  value_f = value; 
  
  return *this; 
}

long int vdt::get_int () const {
  m_check_type (vdt::int_vdt);

  return value_i;
}

double vdt::get_float () const {
  m_check_type (vdt::float_vdt);

  return (type_internal == float_vdt) ? value_f : value_i;	// int->float cast always allowed
}
  
const std::string& vdt::get_string () const {
  //m_check_type (vdt::string_vdt); RETURNING a string is always allowed

  return value_s;
}

const vdt::vdt_vector& vdt::get_vector () const {
  m_check_type (vdt::vector_vdt);

  return value_v;
}

void vdt::m_check_type (vdt_type new_type) const {

  if (type_internal == new_type) return;
  
  if (new_type == float_vdt && type_internal == int_vdt) return;	// int->float cast always allowed
  
  bx_message msg(bx_message::critic, "vdt: ");
  if (name_s.size ()) msg << "\"" << name_s << "\" ";
  msg << "requested " << char(new_type) << " while curren type is " << char(type_internal) << " (" << *this << ") cast not allowed" << dispatch;
}


std::ostream& operator<< (std::ostream& out, const vdt& v) {
  switch (v.get_type ()) {
    case vdt::int_vdt:
      out << v.get_int ();
      break;
    case vdt::float_vdt:
      out << v.get_float ();
      break;
    case vdt::string_vdt:
      out << v.get_string ();
      break;
    case vdt::vector_vdt:
      out << "{ ";
      if (v.get_vector ().size ()) {
        std::copy (v.get_vector ().begin (), v.get_vector ().end () - 1, std::ostream_iterator<vdt>(out, ", "));
        out << *(v.get_vector ().end () - 1) << " ";
      }
      out << "}";
      break;
    default:
      out << "unknown type"; 
  }
  return out;
}
  
/*
 * $Log: vdt.cc,v $
 * Revision 1.20  2008/09/24 12:13:08  razeto
 * Minor upgrades after discussion with Marco Vignati
 *
 * Revision 1.19  2007-11-20 13:11:06  razeto
 * even better parsing algorithm (fix a small bug)
 *
 * Revision 1.18  2007-11-20 11:31:20  razeto
 * Better parsing (for type identification
 *
 * Revision 1.17  2006-08-26 10:32:48  razeto
 * Removed a useless include
 *
 * Revision 1.16  2006/07/18 12:39:43  razeto
 * Now ostream<<vdt support vector conforming to its input requirements
 *
 * Revision 1.15  2005/07/13 12:33:12  razeto
 * Fixed a bug
 *
 * Revision 1.14  2004/11/26 15:25:10  razeto
 * Thanks to Davide who so kindly underlined a syntactical error, I fixed the "Maintainer" word
 *
 * Revision 1.13  2004/11/26 14:01:20  razeto
 * Added Mantainer field
 *
 * Revision 1.12  2004/09/29 09:59:02  razeto
 * Fixed handling of base 10 expo number
 *
 * Revision 1.11  2004/06/14 11:51:59  razeto
 * Fixed a stamp
 *
 * Revision 1.10  2004/06/01 15:36:14  razeto
 * Upgrade to full support vdt vectors
 *
 * Revision 1.9  2004/06/01 09:09:25  razeto
 * Added array support (with { 1, 2, 3, 4 } syntax)
 *
 * Revision 1.8  2004/05/26 09:20:47  razeto
 * Upgraded: now a vdt object typically knows its name (usefull for error reporting)
 *
 * Revision 1.7  2004/05/25 17:14:27  razeto
 * Upgraded
 *
 * Revision 1.6  2004/05/18 14:22:48  razeto
 * upgraded
 *
 * Revision 1.5  2004/04/09 07:48:19  razeto
 * Fixed a bug for negative numbers. Added int->float read casting
 *
 * Revision 1.4  2004/04/05 13:09:07  razeto
 * Added messanger for exception handling. Updated a comment
 *
 * Revision 1.3  2004/04/01 12:05:01  razeto
 * Fixed a bug in the assignement operator.
 * Added a ostreamer function.
 * Added stupid assignement operator since casting is not done.
 *
 * Revision 1.2  2004/03/29 13:23:57  razeto
 * Added the float/int int/float casting
 *
 * Revision 1.1  2004/03/29 13:14:12  razeto
 * Added interface directory (where options and dbi reside).
 * Added vdt object for storing parameters
 *
 */
