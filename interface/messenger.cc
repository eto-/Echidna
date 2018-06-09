/* BOREXINO Reconstruction program
 *
 * Author: Alessandro Razeto <Alessandro.Razeto@ge.infn.it>
 * Maintainer: Alessandro Razeto <Alessandro.Razeto@ge.infn.it>
 *
 * $Id: messenger.cc,v 1.27 2008/10/10 09:19:49 razeto Exp $
 *
 * Implementation of messenger and bx_message
 *
 */
#include "messenger.hh"
#ifndef _ECHIDNA_ROOTLIB_
#include "barn_interface.hh"
#include "bx_named.hh"
namespace {
  bx_named *my_messenger = 0;
}
#endif

messenger* messenger::me = 0;

messenger* messenger::get () {
  if (!me) 
    me = new messenger;

  return me;
}

messenger::messenger (): i4_error_count(0), i4_warn_count(0), b_line_buffered(false) {}


void messenger::open_log_file (const std::string& str) {
  // Open log file
  if (log.is_open () || !str.size ()) return;

  log.open (str.c_str ());
  if (!log) throw std::runtime_error("Can not open file " + str + " for logging");
}


void messenger::submit (bx_message &msg) throw(std::runtime_error) {
  if (msg.str ().substr (msg.str ().size() - 1, 1) != std::string("\n")) msg << std::endl;

  if (msg.check_print_level ()) {
    std::cout << ">>> " << msg.level () << ": " << msg.str ();
    if (b_line_buffered) std::cout << std::flush;
  }
  if (log && msg.check_log_level ()) {
    log << ">>> " << msg.level () << ": " << msg.str ();
    if (b_line_buffered) log << std::flush;
  }
 
#ifndef _ECHIDNA_ROOTLIB_
  if (msg.get_level () >= bx_message::error) {
    if (bx_named::parameter_broker_inited ()) {
      if (!my_messenger) my_messenger = new bx_named ("messenger");
      TObjString rmsg(msg.str ().c_str ());
      barn_interface::get ()->network_send (&rmsg, my_messenger);
    }
  }
#endif
  switch (msg.get_level ()) {
    case bx_message::critic:
      if (log) log.close ();
      throw std::runtime_error (msg.str ());
    case bx_message::error:
      i4_error_count++;
      break;
    case bx_message::warn:
      i4_warn_count++;
      break;
    default: break;
  }
} 

bx_message::message_level bx_message::default_print_level = bx_message::warn;
bx_message::message_level bx_message::default_log_level = bx_message::log;
std::cmap<bx_message::message_level, std::string> bx_message::message_codes("message_codes");

void bx_message::m_init () {
  if (message_codes.size ()) return;
  message_codes[bx_message::none] = "none";
  message_codes[bx_message::debug] = "debug";
  message_codes[bx_message::log] = "log";
  message_codes[bx_message::info] = "info";
  message_codes[bx_message::warn] = "warn";
  message_codes[bx_message::error] = "error";
  message_codes[bx_message::critic] = "critic";
}

bx_message::message_level bx_message::m_decode (const vdt& v) const {

  bx_message::message_level ret_val = bx_message::none;
  
  bx_message error_msg(bx_message::warn, "bx_message: ");

  if (v.get_type () == vdt::string_vdt) {
    std::cmap<message_level, std::string>::const_iterator it = message_codes.rfind (v.get_string ());
    if (it == message_codes.end ()) {
      error_msg << "can not find message level for \"" << v <<"\", assuming \"none\"";
      error_msg.dispatch ();
    } else ret_val = it->first;
  } else if (v.get_type () == vdt::int_vdt) {
    if (v.get_int () < bx_message::none || v.get_int () > bx_message::critic) {
      error_msg << "can not find message level for \"" << v <<"\", assuming \"none\""; 
      error_msg.dispatch ();
    } else ret_val = bx_message::message_level(v.get_int ());
  } else {
    error_msg << "can not find message level for \"" << v <<"\", assuming \"none\""; 
    error_msg.dispatch ();
  }
    
  return ret_val;
}

bx_message &message_client::get_message (bx_message::message_level level) {
  message.clear (my_name + ": ");
  message.set_level (level);
  if (my_print_level != bx_message::none) message.set_print_level (my_print_level);
  return message;
}

/*
 * $Log: messenger.cc,v $
 * Revision 1.27  2008/10/10 09:19:49  razeto
 * send only error and critic messages
 *
 * Revision 1.26  2007-06-04 13:26:54  razeto
 * Use default log level if none level is specified
 *
 * Revision 1.25  2007-03-30 13:39:28  razeto
 * Added line buffering option (-V) for online echidna
 *
 * Revision 1.24  2007-02-07 13:12:53  razeto
 * Reference instead of value
 *
 * Revision 1.23  2007/02/02 11:53:00  razeto
 * Parameter broker must be initialized before creating a bx_named
 *
 * Revision 1.22  2007-02-01 16:56:25  razeto
 * Added network dispatching to messenger
 *
 * Revision 1.21  2006/05/08 14:24:45  razeto
 * Added throw declarations to messanger dispatch chain
 *
 * Revision 1.20  2006/02/12 11:31:38  razeto
 * Added message_client class as an abstraction for bx_bhys
 *
 * Revision 1.19  2006/02/01 16:10:07  razeto
 * Fixed a small bug (check if log is open before close)
 *
 * Revision 1.18  2005/08/02 14:01:49  razeto
 * Removed some useless include
 *
 * Revision 1.17  2005/07/06 12:25:34  razeto
 * Added initialization for some count variables
 *
 * Revision 1.16  2004/11/26 15:25:10  razeto
 * Thanks to Davide who so kindly underlined a syntactical error, I fixed the "Maintainer" word
 *
 * Revision 1.15  2004/11/26 14:01:20  razeto
 * Added Mantainer field
 *
 * Revision 1.14  2004/11/24 13:04:10  razeto
 * Upgraded to the new cmap ctor with name
 *
 * Revision 1.13  2004/10/25 02:48:37  razeto
 * Added a better "end of echidna" message diplaying the number of errors and warns
 *
 * Revision 1.12  2004/09/30 12:53:10  razeto
 * Added bx_message::log level
 *
 * Revision 1.11  2004/08/30 17:07:05  razeto
 * Some simple optimization (mainly inlined ctor)
 *
 * Revision 1.10  2004/08/06 10:30:28  razeto
 * cycle_1 branch merged in the main trunk.
 *
 * Revision 1.9.2.1  2004/07/26 14:39:16  razeto
 * Fixed open log file routine.
 *
 * Revision 1.9  2004/07/22 10:13:01  razeto
 * Jumbo patch for the configuration manager applied, read the documentation for usage
 *
 * Revision 1.8  2004/05/18 14:23:06  razeto
 * Fixed a bug
 *
 * Revision 1.7  2004/04/18 10:30:14  razeto
 * Fixed to compile wit g++ 2.95 present on the cluster
 *
 * Revision 1.6  2004/04/12 16:07:40  razeto
 * Added configurable default print level for messages
 *
 * Revision 1.5  2004/04/09 07:47:28  razeto
 * Updated: changed the log level to info, added a string field to the bx_message ctor
 *
 * Revision 1.4  2004/04/06 12:39:28  razeto
 * Added logfile handling
 *
 * Revision 1.3  2004/04/05 13:10:23  razeto
 * Added comments, added an error level, changed some methods
 *
 * Revision 1.2  2004/04/02 14:04:42  razeto
 * Added some comments. Fixed a typo
 *
 */
