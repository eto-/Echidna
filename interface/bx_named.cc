/* BOREXINO Reconstruction program
 *
 * Author: Alessandro Razeto <Alessandro.Razeto@ge.infn.it>
 * Maintainer: Alessandro Razeto <Alessandro.Razeto@ge.infn.it>
 *
 * $Id: bx_named.cc,v 1.6 2006/02/12 11:28:29 razeto Exp $
 *
 * Implementation of bx_named
 *
 */
#include "bx_named.hh"

parameter_broker *bx_named::p_broker = 0;

void bx_named::m_error_no_pbroker () {
  bx_message msg(bx_message::critic, s_name);
  msg << "trying to initialize a bx_named (\"" << s_name << "\") without the parameter broker" << dispatch;
}

void bx_named::set_parameter (const std::string& name, const vdt& value) {
  p_broker->set_parameter (s_name, name, value);

    // keep alligned the cached values
  if (name == "log_level") log_level = value;
  else if (name == "print_level") print_level = value;
}

bx_message &bx_named::get_message (bx_message::message_level level) {
  message.clear (get_name () + ": ");


  if (!b_level_cached) { // Assigne cached values unless it has already been done
      // assign the print and log level values
    if (check_parameter ("log_level")) log_level = get_parameter ("log_level");
    else {
      log_level = int (bx_message::get_default_log_level ());
      log_level.set_name ("log_level");
    }
    if (check_parameter ("print_level")) print_level = get_parameter ("print_level");
    else {
      print_level = int (bx_message::get_default_print_level ());
      print_level.set_name ("print_level");
    }
    b_level_cached = true;
  }
  message.set_level (level); 
  message.set_log_level (log_level);
  message.set_print_level (print_level);

  return message;
}
/*
 * $Log: bx_named.cc,v $
 * Revision 1.6  2006/02/12 11:28:29  razeto
 * Fixed an harmless bug
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
