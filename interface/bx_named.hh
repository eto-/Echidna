/* BOREXINO Reconstruction program
 *
 * Author: Alessandro Razeto <Alessandro.Razeto@ge.infn.it>
 * Maintainer: Alessandro Razeto <Alessandro.Razeto@ge.infn.it>
 *
 * $Id: bx_named.hh,v 1.10 2007/02/07 13:13:52 razeto Exp $
 *
 * A bx_named in the ancestor of each named object; the properties of named 
 * object are:
 * 1) having a name
 * 2) parameter interface: the parameter stored in the parameter_broker are available
 * to each bx_named (using name as key)
 * 3) bx_message interface
 * Module are bx_named but not only; other objects can inherit from bx_named and 
 * acquire the bx_named properties.
 *
 * The bx_named::set_parameter is virtual and can be overloaded, but in such case the 
 * overloading method MUST call internally bx_named::set_parameter(...) to really set the
 * parameter.
 *
 * The parameter "print_level" and "log_level" are cached internally.
 *
 */
#ifndef _BX_NAMED_H
#define _BX_NAMED_H

#include "vdt.hh"
#include "messenger.hh"
#include "parameter_broker.hh"
#include "cmap.hh"

#include <string>

class bx_named {
  public:
    bx_named (const std::string& myname): s_name(myname), b_level_cached(false) { if (!p_broker) m_error_no_pbroker (); } // This is faster since can be inlined
    virtual ~bx_named () {}

      // Get the name
    const std::string& get_name () const { return s_name; }

      // Parameter interface
    const vdt& get_parameter (const std::string& name) const { return p_broker->get_parameter (s_name, name); }
    const vdt& get_parameter_other (const std::string& other, const std::string& name) const { return p_broker->get_parameter (other, name); }
    bool check_parameter     (const std::string& name) const { return p_broker->check_parameter (s_name, name); }
    bool check_parameter_other (const std::string& other, const std::string& name) const { return p_broker->check_parameter (other, name); }
    virtual void set_parameter (const std::string& name, const vdt& value);
      // This method is not virtual since it is simply a shorthand for calling 
      // the previus (even if this method specialize the previus, since a vdt
      // can always be built from a string this will create an unnamed vdt, 
      // this method will create directly the vdt assigning it a name -more
      // usefull for error messages from the vdt object-)
      // This method will call the current overloaded version of the previus
      // method, so can be used in child of bx_named overriding set_parameter.
    void set_parameter (const std::string& name, const std::string& value) { set_parameter (name, vdt (value, name)); }

    
      // Get a module message. 
      // !!! get_message give a initialized message, so if the message is used 
      // in multiline command a reference has been used
      // get_message () << line1;
      // get_message () << line2 << dispatch; only dumps line2 !!!!!!!
      // instead using a reference like bx_message &msg = get_message ();
      // msg << line1; msg << line2 << dispatch; is OK.
    bx_message &get_message (bx_message::message_level level = bx_message::none);

    const bx_named& operator= (const bx_named& r) { s_name = r.get_name (); return *this; }

    static bool parameter_broker_inited () { return p_broker; }
    static const std::string get_configuration_name () { return p_broker->get_configuration_name (); }
  private:
    std::string s_name;
    bx_message message; 
    bool b_level_cached;
    vdt log_level, print_level;
    void m_error_no_pbroker (); // internal routine used just to print an error on ctor.  This is done 
    				// to reduce the ctor as much as possible.

    static parameter_broker* p_broker;
    static void set_parameter_broker (parameter_broker* broker) { p_broker = broker; }
    friend class bx_options;
};

#endif
/*
 * $Log: bx_named.hh,v $
 * Revision 1.10  2007/02/07 13:13:52  razeto
 * get_configuration_name is static (since only depende upon parameter broker)
 *
 * Revision 1.9  2007/02/02 11:53:00  razeto
 * Parameter broker must be initialized before creating a bx_named
 *
 * Revision 1.8  2006-12-05 13:37:20  razeto
 * Added configuration name
 *
 * Revision 1.7  2006/10/18 15:48:31  razeto
 * Now bx_named can get/check parameters of other bx_named
 *
 * Revision 1.6  2005/08/02 14:01:49  razeto
 * Removed some useless include
 *
 * Revision 1.5  2004/11/26 15:25:10  razeto
 * Thanks to Davide who so kindly underlined a syntactical error, I fixed the "Maintainer" word
 *
 * Revision 1.4  2004/11/26 14:01:20  razeto
 * Added Mantainer field
 *
 * Revision 1.3  2004/08/30 17:10:31  razeto
 * Some optimization (mainly inlined ctor and retarded parameter initialization)
 *
 * Revision 1.2  2004/07/22 10:13:01  razeto
 * Jumbo patch for the configuration manager applied, read the documentation for usage
 *
 * Revision 1.1  2004/06/10 12:38:19  razeto
 * Added bx_named and parameter_broker interface
 *
 */
