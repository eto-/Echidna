/* BOREXINO Reconstruction program
 *
 * Author: Alessandro Razeto <Alessandro.Razeto@ge.infn.it>
 * Maintainer: Alessandro Razeto <Alessandro.Razeto@ge.infn.it>
 *
 * $Id: vdt.hh,v 1.13 2006/07/18 14:35:12 razeto Exp $
 *
 * Variable Data Type is an object which is constructed from
 * a string and holts the type of the assigne value.
 * vdt("1000") creates a vdt obj holding an integer with value
 * 1000. Later the obj value can be fetced only if the query
 * matches the vdt type, else an exception is generated.
 * A assignement operator is present, which does not allow
 * runtime type change (except float/int).
 *
 */
#ifndef _VDT_H
#define _VDT_H
#include <string>
#include <iostream>
#include <vector>

class bx_vdt {
  public:
    enum vdt_type {
      int_vdt = 'I',
      float_vdt = 'F',
      string_vdt = 'S',
      vector_vdt = 'V',
      unassigned_vdt = 'U'
    };
    typedef std::vector<bx_vdt> vdt_vector;
    
    bx_vdt (): type_internal(unassigned_vdt), name_s("") {}
    bx_vdt (const std::string& value, const std::string& name = "", vdt_type type = unassigned_vdt);
    bx_vdt (long int value): name_s("") { value_i = value; type_internal = int_vdt; }
    bx_vdt (double value): name_s("") { value_f = value; type_internal = float_vdt; }
    template<typename V> bx_vdt (const std::vector<V, std::allocator<V> >& value): name_s("") { // This should stay in .cc file, but
      										 	     // gcc does not have export for templates
      type_internal = vector_vdt;
      std::back_insert_iterator<std::vector<bx_vdt> > inserter(value_v);
      std::copy (value.begin (), value.end (), inserter);
    }
      // Since cast is not called in some onstructor operator these methods are necessary
    bx_vdt (int value): name_s("") { value_i = value; type_internal = int_vdt; }
    bx_vdt (short value): name_s("") { value_i = value; type_internal = int_vdt; }
    bx_vdt (char value): name_s("") { value_i = value; type_internal = int_vdt; }
    bx_vdt (unsigned long value): name_s("") { value_i = value; type_internal = int_vdt; }
    bx_vdt (unsigned short value): name_s("") { value_i = value; type_internal = int_vdt; }
    bx_vdt (unsigned char value): name_s("") { value_i = value; type_internal = int_vdt; }
    bx_vdt (float value): name_s("") { value_f = value; type_internal = float_vdt; }

    void set_name (const std::string& name) { name_s = name; }

    const bx_vdt& operator= (const std::string& value); 
    const bx_vdt& operator= (long int value);
    const bx_vdt& operator= (double value);
    template<typename V> const bx_vdt& operator= (const std::vector<V, std::allocator<V> >& value) { // This should stay in .cc file, but
      												  // gcc does not have export for templates
      if (type_internal != unassigned_vdt) m_check_type (bx_vdt::vector_vdt);
      type_internal = vector_vdt;
      std::back_insert_iterator<std::vector<bx_vdt> > inserter(value_v);
      std::copy (value.begin (), value.end (), inserter);
      return *this;
    }

      // Since cast is not called in assignement operator these methods are necessary
    const bx_vdt& operator= (int value)      { return (*this) = (long int)(value); }
    const bx_vdt& operator= (short value)	  { return (*this) = (long int)(value); }
    const bx_vdt& operator= (char value)     { return (*this) = (long int)(value); }
    const bx_vdt& operator= (unsigned long value) { return (*this) = (long int)(value); }
    const bx_vdt& operator= (unsigned short value) { return (*this) = (long int)(value); }
    const bx_vdt& operator= (unsigned char value) { return (*this) = (long int)(value); }
    const bx_vdt& operator= (float value)    { return (*this) = double(value); }

    bool is_valid ()	 const  { return (type_internal == unassigned_vdt) ? 0 : 1; }
    
    vdt_type get_type () const	{ return type_internal; }
    long int get_int () const;
    double get_float () const;
    const std::string& get_string () const;
    const vdt_vector& get_vector () const;
    bool get_bool () const { return bool (get_int ()); }
    
  private:
    vdt_type type_internal;
    
    long int value_i;
    double value_f;
    std::string value_s;
    vdt_vector value_v;
    
    std::string name_s;
    void m_check_type (vdt_type new_type) const;  // Throws an exception if the internal type 
    						  // is different from the required type
    vdt_type m_search_type (const std::string& str);
    void m_assign_value (const std::string& str);
};

typedef bx_vdt vdt;
std::ostream& operator<< (std::ostream&, const bx_vdt&);
#endif
/*
 * $Log: vdt.hh,v $
 * Revision 1.13  2006/07/18 14:35:12  razeto
 * Upgraded to handle vector in construction/assignement
 *
 * Revision 1.12  2005/08/02 14:01:49  razeto
 * Removed some useless include
 *
 * Revision 1.11  2004/11/26 15:25:10  razeto
 * Thanks to Davide who so kindly underlined a syntactical error, I fixed the "Maintainer" word
 *
 * Revision 1.10  2004/11/26 14:01:20  razeto
 * Added Mantainer field
 *
 * Revision 1.9  2004/10/19 16:12:47  razeto
 * Added a get_bool method and improved ctors
 *
 * Revision 1.8  2004/06/01 09:09:25  razeto
 * Added array support (with { 1, 2, 3, 4 } syntax)
 *
 * Revision 1.7  2004/05/26 09:20:47  razeto
 * Upgraded: now a vdt object typically knows its name (usefull for error reporting)
 *
 * Revision 1.6  2004/05/25 17:14:27  razeto
 * Upgraded
 *
 * Revision 1.5  2004/04/09 07:48:19  razeto
 * Fixed a bug for negative numbers. Added int->float read casting
 *
 * Revision 1.4  2004/04/01 12:05:01  razeto
 * Fixed a bug in the assignement operator.
 * Added a ostreamer function.
 * Added stupid assignement operator since casting is not done.
 *
 * Revision 1.3  2004/03/29 13:26:07  razeto
 * Fixed a typo
 *
 * Revision 1.2  2004/03/29 13:23:27  razeto
 * Fixed a bug, object could be uninitialized
 *
 * Revision 1.1  2004/03/29 13:14:12  razeto
 * Added interface directory (where options and dbi reside).
 * Added vdt object for storing parameters
 *
 */
