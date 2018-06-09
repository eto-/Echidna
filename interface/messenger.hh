/* BOREXINO Reconstruction program
 *
 * Author: Alessandro Razeto <Alessandro.Razeto@ge.infn.it>
 * Maintainer: Alessandro Razeto <Alessandro.Razeto@ge.infn.it>
 *
 * $Id: messenger.hh,v 1.23 2011/01/15 13:04:59 razeto Exp $
 *
 * 2 classes are defined: 
 *   messenger: which is the singleton which receives messages and dump
 *   them on monitor and on logfile
 *   bx_message: which is the message created from modules to diagnose
 *   something. They are ostream, so the standard cout syntax (<<) can be
 *   used. A specialized dispatch manipulator has to be called at the
 *   end of the message to dispatch it to the messenger.
 *   Module usage is simplified by the bx_base_module::get_message(message_level)
 *   which fills all the level fields.
 *   For other uses an istance of bx_message has to be initialized using
 *   set_level and set_*_level.
 *   A message is printed if its level is >= than the print_level (idem for log).
 *   
 */
#ifndef _MESSANGER_H
#define _MESSANGER_H
#include <string>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <stdexcept>
#include <typeinfo>

#include "cmap.hh"
#include "vdt.hh"

class bx_message;

class messenger {
  public:
    static messenger* get ();

    void open_log_file (const std::string& str);
    void submit (bx_message &) throw(std::runtime_error);
    void set_line_buffered (bool on) { b_line_buffered = on; }
    void flush () { if (log) log.close (); }

    int get_error_count () const { return i4_error_count; }
    int get_warn_count () const { return i4_warn_count; }
  private:
    messenger ();
    static messenger *me;
    std::ofstream log;
    int i4_error_count, i4_warn_count;
    bool b_line_buffered;
};

  // a bx_message inherits from ostringstream, so it is a stream similar to cout,
  // then the message text can be written with << operators.
  // Alternativelly a () oerator is present so that the message text can be set
  // whith bx_message msg; msg("text").dispatch();
class bx_message: public std::ostringstream {
  public:
    enum message_level { 
      none = 0,
      debug = 1,
      log = 2,
      info = 3,
      warn = 4,
      error = 5,
      critic = 6,
    };

    bx_message (message_level level = none, const std::string& initial_message = ""): i4_print_level(default_print_level), i4_log_level(default_log_level), i4_message_level(level) { if (!message_codes.size ()) m_init (); exceptions (std::ios::badbit); if (initial_message.size ()) *this << initial_message; }
    bx_message (const bx_message& a): std::ostringstream(a.str ()), i4_print_level(a.i4_print_level), i4_log_level(a.i4_log_level), i4_message_level(a.i4_message_level) {}

      // Access to the default values
    static void set_default_print_level (message_level level) { default_print_level = level; }
    static message_level get_default_print_level () { return default_print_level; }
    static message_level get_default_log_level () { return default_log_level; }


      // Clear the message both text and levels
    void clear (const std::string& txt="") { str(""); if (txt.size ()) *this << txt; i4_print_level = default_print_level; i4_log_level = default_log_level; i4_message_level = none; } 

      // Set the message text (equivalent to message << str) 
    bx_message& operator() (const std::string& str) { return (bx_message &)(*this << str); } // no need to intelligent cast since
    											     // this is known to be a bx_message!!!

      // Set the message level (can be specified in the ctor even)
    void set_level (message_level level) { i4_message_level = level; }
      // Set the verbosity of this message: if it has to be printed and logged.
      // These methods are meant to be used by bx_base_module::get which knows
      // the module verbosity levels.
    void set_log_level (message_level level) { i4_log_level = level; }
    void set_print_level (message_level level) { i4_print_level = level; }
      // These methos are equivalent to the previus but accept a generic
      // vdt value (as got from the modules parameters map)
    void set_print_level (const vdt& level) { i4_print_level = m_decode(level); }
    void set_log_level (const vdt& level) { i4_log_level = m_decode(level); }

      // Get the string level representation
    const std::string& level () const { return message_codes[i4_message_level]; }
    
      // Some methods to chech the message verbosity for logs and printout
    bool check_log_level () const { return (i4_log_level == none) ? false : (i4_message_level >= i4_log_level); }
    bool check_print_level () const { return (i4_print_level == none) ? false : (i4_message_level >= i4_print_level); }

      // Get the message level
    message_level get_level () const { return i4_message_level; }
    
      // Dispatch the message to the messenger. Equivalent to message << dispatch
    void dispatch () throw(std::runtime_error) { messenger::get ()->submit (*this); }
  private:
    message_level i4_print_level;
    message_level i4_log_level;
    message_level i4_message_level;
    message_level m_decode (const vdt& v) const;
    static message_level default_print_level;
    static message_level default_log_level;
    static void m_init ();
    static std::cmap<message_level, std::string> message_codes;
};

  // function to use bx_message << dispatch sintax. Equivalent to message.dispatch () 
inline std::ostream& dispatch (std::ostream& stream) throw(std::runtime_error, std::bad_cast) {
  (dynamic_cast<bx_message&>(stream)).dispatch (); 
  return stream; 
}

inline std::ostream& newline (std::ostream& stream) {
  stream << std::endl;
  bx_message *msg = dynamic_cast<bx_message *>(&stream);
  if (msg) {
    std::ios::fmtflags old_opt = stream.flags ();
    stream.setf(std::ios::left, std::ios::adjustfield);
    stream << std::setw (msg->level ().size () + 3 + 3) << ">>>";
    stream.flags (old_opt);
  }
  return stream; 
}

  // To be used outside echidna (by bx_phys) without parameter broker 
class message_client {
  protected:
      // Using none will allow to use default message level
    message_client (const std::string& client_name, bx_message::message_level print_level = bx_message::none): my_name(client_name), my_print_level(print_level) {}
    bx_message &get_message (bx_message::message_level level);
    const std::string& get_my_name () const { return my_name; }
  private:
    std::string my_name;
    bx_message::message_level my_print_level;
    bx_message message;
};

#endif
/*
 * $Log: messenger.hh,v $
 * Revision 1.23  2011/01/15 13:04:59  razeto
 * Added my_name getter
 *
 * Revision 1.22  2007-06-04 13:28:35  razeto
 * Use default log level if none level is specified
 *
 * Revision 1.21  2007-03-30 13:39:28  razeto
 * Added line buffering option (-V) for online echidna
 *
 * Revision 1.20  2007-02-07 13:12:53  razeto
 * Reference instead of value
 *
 * Revision 1.19  2006/05/10 12:03:35  razeto
 * Added 2 include for exceptions
 *
 * Revision 1.18  2006/05/08 14:24:45  razeto
 * Added throw declarations to messanger dispatch chain
 *
 * Revision 1.17  2006/02/12 11:31:38  razeto
 * Added message_client class as an abstraction for bx_bhys
 *
 * Revision 1.16  2006/01/09 16:17:43  razeto
 * Fixed a typo
 *
 * Revision 1.15  2005/06/10 13:08:34  razeto
 * Added a new bx_message manipulator: newline. It should be used in the place
 * of std::endl to have multiple lines messages.
 *
 * Revision 1.14  2004/11/26 15:25:10  razeto
 * Thanks to Davide who so kindly underlined a syntactical error, I fixed the "Maintainer" word
 *
 * Revision 1.13  2004/11/26 14:01:20  razeto
 * Added Mantainer field
 *
 * Revision 1.12  2004/10/25 02:48:37  razeto
 * Added a better "end of echidna" message diplaying the number of errors and warns
 *
 * Revision 1.11  2004/09/30 13:36:40  razeto
 * Added bx_message::log level
 *
 * Revision 1.10  2004/09/16 16:00:35  razeto
 * fixed a typo
 *
 * Revision 1.9  2004/08/30 17:07:05  razeto
 * Some simple optimization (mainly inlined ctor)
 *
 * Revision 1.8  2004/07/22 10:13:01  razeto
 * Jumbo patch for the configuration manager applied, read the documentation for usage
 *
 * Revision 1.7  2004/04/12 16:07:40  razeto
 * Added configurable default print level for messages
 *
 * Revision 1.6  2004/04/09 07:47:28  razeto
 * Updated: changed the log level to info, added a string field to the bx_message ctor
 *
 * Revision 1.5  2004/04/06 12:39:28  razeto
 * Added logfile handling
 *
 * Revision 1.4  2004/04/05 13:10:23  razeto
 * Added comments, added an error level, changed some methods
 *
 * Revision 1.3  2004/04/03 09:21:53  razeto
 * fixed a bug
 *
 * Revision 1.2  2004/04/02 14:04:42  razeto
 * Added some comments. Fixed a typo
 *
 */
