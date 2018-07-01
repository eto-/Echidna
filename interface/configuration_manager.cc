/* BOREXINO Reconstruction program
 * 
 * Author: Alessandro Razeto <Alessandro.Razeto@ge.infn.it>
 * Maintainer: Alessandro Razeto <Alessandro.Razeto@ge.infn.it>
 *
 * $Id: configuration_manager.cc,v 1.8 2006/12/05 12:17:18 razeto Exp $
 *
 * Implementation of configuration_manager classes
 * Some handler classes are defined (configuration_{file|stanza|line} to help
 * parsing files.
 *
 */
#include "configuration_manager.hh"

#include <iterator>
#include <functional>
#include <algorithm>
#include <unistd.h>
#include <string.h>
#include <errno.h>


// ==================== HELPER CLASSES (definition) ======================= //
namespace {
    // The configuration_file is the main interface to the configuration files: dependig on the
    // user_file flags it expects a main config file or a user config file. The first is a list of
    // stanzas while the second is just a list of statements. All the stuff is done at the constructor.
  class configuration_file {
    public:
      configuration_file (const std::string& filename, bool user_file, configuration_manager *cm);
    private:
  };

  class configuration_line;

    // The configuration_stanza is the low level interface for the main configuration file, which is 
    // a list of stanzas. It parses the stanza and is send_to_cm is true sends the results (a configuration_item)
    // to the configuration_manager. All the stuff is done at the constructor.
  class configuration_stanza {
    public:
      configuration_stanza (const configuration_line& first_line, std::istream& in, configuration_manager *cm);
  };

    // The configuration_line is the interface to a line of configuration file; it 
    // reads the file untill a valid line is found (or EOF is reached).
    // When a non-empty line is found (after removing comments and spaces) it is parsed
    // to find which line type it is, the results are stored in the bool flags and in the
    // token[12] strings. All the stuff is done at the constructor.
  class configuration_line {
    public:
      configuration_line (std::istream& in);

      bool is_valid () const { return b_valid ;}
      bool is_stanza_begin () const { return b_stanza_begin; }
      bool is_stanza_end () const { return b_stanza_end ;}
      const std::string& get_token1 () const { return s_token1; }
      const std::string& get_token2 () const { return s_token2; }

      static void reset_filename (const std::string& filename) { s_filename = filename; i4_line_count = 0; }
    private:
      bool b_stanza_begin;
      bool b_stanza_end;
      bool b_valid;
      std::string s_token1, s_token2, s_line;

      static std::string s_filename;
      static int i4_line_count;
    friend std::ostream& operator<< (std::ostream& out, const configuration_line& line); // must be friend of both classes
  };
  std::ostream& operator<< (std::ostream& out, const configuration_line& line) { 
    return out << "\"" << line.s_line << "\" at " << line.s_filename << ":" << line.i4_line_count;
  }
};

// ==================== HELPER CLASSES (implementation) ======================= //
// --------------------- Configuration File ------------------------
configuration_file::configuration_file (const std::string& filename, bool user_file, configuration_manager *cm) {

  bx_message msg(bx_message::critic, "configuration_file: ");

    // Check for the file presence
  if (::access (filename.c_str (), R_OK) < 0) msg << "file \"" << filename << "\" access: " << ::strerror (errno) << dispatch;

    // Open the file
  std::ifstream in (filename.c_str ());
  if (!in) msg << "open file \"" << filename << "\" failed" << dispatch;

    // Reset the line counter and set the filename for the line helper class
  configuration_line::reset_filename (filename);

    // Create a configuration item for the user file
  configuration_manager::configuration_namespace user_item (configuration_manager::configuration_namespace::user_file);

    // Loop till end of file
  while (in) {
      // Get line
    configuration_line line(in);

      // End if not valid
    if (!line.is_valid ()) {
      break;
      // If begin of sanza create a valid stanza (<foo bar> not <end>)
    } else if (line.is_stanza_begin ()) {
      configuration_manager *_cm = cm;
      if (user_file) { 
	  // In a local user file sanzas are not allowed, set it as invalid
        bx_message msg2(bx_message::error, "configuration_file: ");
	msg2 << "trying to open a stanza in user file: " << line << dispatch;
	_cm = 0;
      }
        // Create the stanza (which will read lines up to <end> included)
      configuration_stanza(line, in, _cm);
      // Manage syntax error (<end> outside a stanza)
    } else if (line.is_stanza_end ()) {
      msg << "syntax error, end of stanza outside a stanza: " << line << dispatch;
      // If here line is a single assegnation statement outside a stanza
    } else {
      if (!user_file) {
    	  // Single statement not allowed on gloabal config file, skip line
  	bx_message msg2(bx_message::error, "configuration_file: ");
	msg2 << "statement otuside a stanza for a main file: " <<  line << dispatch;
      } else {
	  // Handle the statements
        user_item.add_parameter (configuration_manager::configuration_parameter (line.get_token1 (), line.get_token2 ()));
      }
    }
  }

    // If this is a user file, than user_item is filled with some values and so add it to the cm.
  if (user_file && cm) cm->add_item (user_item);
}

// --------------------- Configuration Stanza ------------------------
configuration_stanza::configuration_stanza (const configuration_line& first_line, std::istream& in, configuration_manager* cm) {

  bx_message msg(bx_message::critic, "configuration_file: ");

  if (!first_line.is_stanza_begin ()) 
    msg << "internal error, configuration_stanza created with invalid line: " << first_line << dispatch;

  configuration_manager::configuration_namespace::item_type type = configuration_manager::configuration_namespace::module;
  if (first_line.get_token1 () == "module") type = configuration_manager::configuration_namespace::module;
  else if (first_line.get_token1 () == "named") type = configuration_manager::configuration_namespace::named;
  else if (first_line.get_token1 () == "config") type = configuration_manager::configuration_namespace::config;
  else {
    bx_message msg2(bx_message::error, "configuration_stanza: ");
    msg2 << "unknown stanza type " << first_line << dispatch;
    cm = 0;
  }

  configuration_manager::configuration_namespace item (type, first_line.get_token2 ());

  //bool stanza_ended = false;
    // Loop till end of stanza
  while (in) {
      // Get line
    configuration_line line(in);
  
      // End if not valid
    if (!line.is_valid ()) {
      break;
      // If begin of sanza create a valid stanza (<foo bar> not <end>)
    } else if (line.is_stanza_begin ()) {
      if (line.get_token1 () == "include" && type == configuration_manager::configuration_namespace::config) {
	const configuration_manager::configuration_namespace& config = cm->get_configuration (line.get_token2 ());
	item += config;
      } else {
          // Nested stanzas not allowed
	bx_message msg2(bx_message::error, "configuration_stanza: ");
	msg2 << "found a nested stanza: " << line << dispatch;
	configuration_stanza(line, in, 0);
      }
      // Finish at the end of stanza
    } else if (line.is_stanza_end ()) {
      //stanza_ended = true;
      break;
    } else {
        // Handle the statements
      item.add_parameter (configuration_manager::configuration_parameter (line.get_token1 (), line.get_token2 ()));
    }
  }
    // Add to the configuration_manager
  if (cm) cm->add_item (item);
}

// --------------------- Configuration Line ------------------------
std::string configuration_line::s_filename;
int configuration_line::i4_line_count;

configuration_line::configuration_line (std::istream& in): b_stanza_begin(false), b_stanza_end(false), b_valid(false) {

  std::string::size_type pos;

  std::string line, line_saved;
    // Loop till end of file
  while (in) {
      // Get line
    std::getline (in, line);
    line_saved = line;

      // Increment internal line counter
    i4_line_count++;

      // Erase inline comments
    pos = line.find ('#');
    if (pos != std::string::npos) line.erase (pos);

      // Erase white space at the beginning of the line
    pos = line.find_first_not_of ("\t ");
    if (pos != 0) line.erase (0, pos);

      // Erase white space at the end of the line
    pos = line.find_last_not_of ("\t ");
    if (pos != std::string::npos) line.erase (pos + 1);

      // Copy line internally
    s_line = line;

      // Skip blank line
    if (!line.empty ()) break;
  }

    // Check validity (end of file reached)
  if (line.empty ()) return;  // b_valid is already false

    // Fine
  b_valid = true;

    // If <end> simply set flag
  if (line == "<end>") {
    b_stanza_end = true;
    return;
  }

    // If here the line can contain or a stanza begin "<foo bar>" or a assegnation line "foo.bar 1"
  if (line[0] == '<' && line[line.size () - 1] == '>') {
    b_stanza_begin = true;
  
      // Remove <> chearcters
    line.erase (0, 1);
    line.erase (line.size () - 1);
  
      // And the possible spaces after and before them
    if (line[0] == ' ' || line[0] == '\t') {
      pos = line.find_first_not_of ("\t ");
      line.erase (0, pos);
    }
    if (line[line.size () - 1] == ' ' || line[line.size () - 1] == '\t') {
      pos = line.find_last_not_of ("\t ");
      line.erase (pos + 1);
    }
  }

    // Now the line should contain the 2 tokens separated by spaces, fill the local variables
  bx_message msg (bx_message::critic, "configuration_line: ");
  msg << "on line \"" << line_saved << "\"";
  pos = line.find_first_of ("\t ");
  if (pos == std::string::npos) msg << " token1 not found" << dispatch;
  s_token1 = line.substr (0, pos);
  pos = line.find_first_not_of ("\t ", pos);
  if (pos == std::string::npos) msg << " token2 not found" << dispatch;
  s_token2 = line.substr (pos);
}
// ==================== HELPER CLASSES (end) ======================= //


// --------------------- Configuration Manager ------------------------
void configuration_manager::add_main_config_file (const std::string& filename) {
  configuration_file (filename, false, this);
}

void configuration_manager::add_user_config_file (const std::string& filename) {
  configuration_file (filename, true, this);
}

void configuration_manager::add_comman_line_parameter (const std::string& name_dot_parameter, const std::string& value) {
  configuration_namespace item (configuration_namespace::command_line_parameters);
  item.add_parameter (configuration_parameter (name_dot_parameter, value));

  add_item (item);
}

void configuration_manager::add_item (const configuration_namespace& item) {
  switch (item.get_type ()) {
    case configuration_namespace::module:
      module_list.push_back (item);
      break;
    case configuration_namespace::named:
      named_list.push_back (item);
      break;
    case configuration_namespace::config:
      config_list.push_back (item);
      break;
    case configuration_namespace::user_file:
      user_file_params.push_back (item);
      break;
    case configuration_namespace::command_line_parameters:
      command_line_params.push_back (item);
      break;
  }
}

const configuration_manager::configuration_namespace& configuration_manager::get_configuration (const std::string& configuration_name) {
  std::vector<configuration_namespace>::const_iterator config = std::find_if (config_list.begin (), config_list.end (), check_ns_name (configuration_name));

  if (config == config_list.end ()) {
    bx_message msg(bx_message::critic, "configuration_manager: ");
    msg << "configuration \"" << configuration_name << "\" not found" << dispatch;
  }

  return *config;
}

std::vector<std::string> configuration_manager::get_module_name_list () const {
  std::vector<std::string> module_name_list;

  for (unsigned long i = 0; i < module_list.size (); i++) module_name_list.push_back (module_list[i].get_namespace_name ());

  return module_name_list;
}
  
  
void configuration_manager::upload_configuration (const std::string& configuration_name) {
  const configuration_namespace& config = get_configuration (configuration_name);

  for (unsigned long i = 0; i < module_list.size (); i++) module_list[i].upload_parameter_broker (this);
  for (unsigned long i = 0; i < named_list.size (); i++) named_list[i].upload_parameter_broker (this);
  config.upload_parameter_broker (this);
  set_configuration_name (configuration_name);
  for (unsigned long i = 0; i < user_file_params.size (); i++) user_file_params[i].upload_parameter_broker (this);
  for (unsigned long i = 0; i < command_line_params.size (); i++) command_line_params[i].upload_parameter_broker (this);
 
}
    

// --------------------- Configuration Item ------------------------
void configuration_manager::configuration_namespace::upload_parameter_broker (parameter_broker *broker) const {
  bx_message msg(bx_message::debug, "configuration_namespace: ");
  if (i_type == module || i_type == named) msg << "initing " << s_ns_name << " with ";
  else msg << "uploading ";
  for (unsigned long i = 0; i < param_v.size (); i++) {
    if (i_type == module || i_type == named) {
      broker->init_parameter (s_ns_name, param_v[i].get_name (), param_v[i].get_value ());
    } else {
      broker->set_parameter (param_v[i].get_name (), param_v[i].get_value ());
    }
  }
  std::copy (param_v.begin (), param_v.end (), std::ostream_iterator<configuration_parameter>(msg, " "));
  msg << dispatch;
}



/*
 * $Log: configuration_manager.cc,v $
 * Revision 1.8  2006/12/05 12:17:18  razeto
 * Add <include configuration> inside an other configuration
 *
 * Revision 1.7  2004/12/10 12:19:45  razeto
 * Removed some compiler warnings
 *
 * Revision 1.6  2004/12/10 12:15:28  razeto
 * Added a error check
 *
 * Revision 1.5  2004/11/26 15:25:10  razeto
 * Thanks to Davide who so kindly underlined a syntactical error, I fixed the "Maintainer" word
 *
 * Revision 1.4  2004/11/26 14:01:20  razeto
 * Added Mantainer field
 *
 * Revision 1.3  2004/08/29 14:33:02  razeto
 * Fixed a compiler warning
 *
 * Revision 1.2  2004/08/06 10:30:28  razeto
 * cycle_1 branch merged in the main trunk.
 *
 * Revision 1.1.2.1  2004/07/22 11:09:59  razeto
 * Fixed and include
 *
 * Revision 1.1  2004/07/22 10:13:01  razeto
 * Jumbo patch for the configuration manager applied, read the documentation for usage
 *
 */

