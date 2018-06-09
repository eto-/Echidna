/* BOREXINO Reconstruction program
 *
 * Author: Alessandro Razeto <Alessandro.Razeto@ge.infn.it>
 * Maintainer: Alessandro Razeto <Alessandro.Razeto@ge.infn.it>
 *
 * $Id: parameter_broker.cc,v 1.10 2005/07/13 12:33:54 razeto Exp $
 *
 * Implementation of parameter_broker
 *
 */
#include "parameter_broker.hh"
#include "messenger.hh"
#include "vdt.hh"
#include <iomanip>

void parameter_broker::init_parameter (const std::string& obj_name, const std::string& param_name, const vdt& value) {
    // A bit complicated syntax follow to avoid calling cmap::cmap() in the map internals
    // The following code is equivalent to "parameter_matrix[obj_name][param_name] = value;"
  if (!parameter_matrix.check (obj_name)) 
    parameter_matrix.insert(std::cmap<std::string, parameter_namespace>::value_type(obj_name, parameter_namespace("parameter_namespace")));
  parameter_matrix.find(obj_name)->second[param_name] = value;
}

void parameter_broker::init_parameter (const std::string& obj_name, const std::string& param_name, const std::string& value) {
  init_parameter (obj_name, param_name, vdt(value, param_name));
}

void parameter_broker::set_parameter (const std::string& obj_name, const std::string& param_name, const std::string& value) {
  set_parameter (obj_name, param_name, vdt(value, param_name));
}

void parameter_broker::set_parameter (const std::string& obj_name, const std::string& param_name, const vdt& value) {
  if (!check_parameter (obj_name, param_name) && param_name != "print_level" && param_name != "log_level") {
    bx_message msg(bx_message::warn, obj_name);
    msg << ": assigning value to unknown parameter \"" << param_name << "\" (" << value << ") will initialize it. Check for syntax error" << dispatch;
  }

  init_parameter (obj_name, param_name, value);
}

void parameter_broker::set_parameter (const std::string& obj_name_dot_param_name, const std::string& value) {
  std::string::size_type dot_position = obj_name_dot_param_name.find (".");
  if (dot_position == std::string::npos || dot_position == 0) {
    bx_message msg(bx_message::critic, "parameter_broker: ");
    msg << "set_module_dot_parameter_value require module.parameter syntax (" << obj_name_dot_param_name << ")" << dispatch;
  }

  std::string obj_name = obj_name_dot_param_name.substr (0, dot_position);
  std::string parameter_name = obj_name_dot_param_name.substr (dot_position + 1);

  set_parameter (obj_name, parameter_name, value);
}

void parameter_broker::m_error_param_not_found (const std::runtime_error &er, const std::string& obj_name, const std::string& param_name) const {
  std::string str = er.what ();

  if (str.substr (0, 13) == "cmap: element") {
    bx_message msg(bx_message::critic, "parameter_broker: ");
    msg << "parameter " << param_name << " not found in the " << obj_name << " namespace" << dispatch;
  } else throw er;
}

void parameter_broker::log_configuration () const {
  bx_message msg(bx_message::log, "parameter_broker: configuration dump:");
  msg << std::endl;
  for (std::cmap<std::string, parameter_namespace>::const_iterator ns = parameter_matrix.begin (); ns != parameter_matrix.end (); ns++) {
    msg << "\t" << ns->first << " namespace:" << std::endl;
    for (parameter_namespace::const_iterator it = ns->second.begin (); it != ns->second.end (); it++) 
#if __GNUC__ >= 3
      msg << "\t    " << std::setw(30) << std::left << it->first << " = " << it->second << std::endl;
#else
      msg << "\t    " << std::setw(30) << it->first << " = " << it->second << std::endl;
#endif

  }
  msg << dispatch;
}
/*
 * $Log: parameter_broker.cc,v $
 * Revision 1.10  2005/07/13 12:33:54  razeto
 * Updated a warning
 *
 * Revision 1.9  2004/11/26 15:25:10  razeto
 * Thanks to Davide who so kindly underlined a syntactical error, I fixed the "Maintainer" word
 *
 * Revision 1.8  2004/11/26 14:01:20  razeto
 * Added Mantainer field
 *
 * Revision 1.7  2004/11/24 13:04:33  razeto
 * Upgraded to the new cmap ctor with name
 *
 * Revision 1.6  2004/10/06 11:04:54  razeto
 * Changed a printout to compile on the cluster
 *
 * Revision 1.5  2004/09/30 14:37:22  razeto
 * Added parameter_broker::log_configuration to log the current configuration
 *
 * Revision 1.4  2004/08/30 17:05:15  razeto
 * Upgraded to check for parameter existance error
 *
 * Revision 1.3  2004/07/22 10:13:01  razeto
 * Jumbo patch for the configuration manager applied, read the documentation for usage
 *
 * Revision 1.2  2004/07/16 07:33:45  razeto
 * Fixed a malformed warning message
 *
 * Revision 1.1  2004/06/10 12:38:19  razeto
 * Added bx_named and parameter_broker interface
 *
 */
