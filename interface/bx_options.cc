/* BOREXINO Reconstruction program
 * 
 * Author: Alessandro Razeto <Alessandro.Razeto@ge.infn.it>
 * Maintainer: Alessandro Razeto <Alessandro.Razeto@ge.infn.it>
 *
 * $Id: bx_options.cc,v 1.28 2008/08/05 14:38:41 razeto Exp $
 *
 * Implementation of bx_options
 * Here bx_message messages are not used and direct throw are called on error,
 * this is done because every error is a line option/config file parsing error 
 * and have to be reported directly.
 *
 */
#include "bx_options.hh"
#include "bx_module_factory.hh"
#include "messenger.hh"
#include "configuration_manager.hh"

#include <iostream>
#include <fstream>
#include <sstream>
#include <unistd.h>
#include <time.h>
#include <string.h>
#include <bx_dbi.hh>

bx_options *bx_options::me = 0;



void bx_options::initialize_echidna (int argc, char * const argv[]) {
  if (me) return;
  me = new bx_options (argc, argv);
}

void bx_options::finish () {
  if (!me) return;
  if (!me->required_log_filename && me->log_filename.size () && bx_named::parameter_broker_inited ()) {
    std::ostringstream new_name;
    if (me->log_path.size ()) new_name << me->log_path << "/";
    new_name << "Run" << std::setw (6) << std::setfill ('0') << bx_dbi::get ()->get_current_run_number ();
    if (bx_named::get_configuration_name () != "default") new_name << "_" << bx_named::get_configuration_name ();
    new_name << "_c" << CYCLE_NUMBER << ".log";
    ::rename (me->log_filename.c_str (), new_name.str ().c_str ());
    std::cout << "Logfile renamed " << new_name.str () << std::endl;
  }
}

bx_options::bx_options (int argc, char * const argv[]): required_log_filename(false) {

  extern char *optarg;
  extern int optind, opterr, optopt;

  configuration_manager *cm = new configuration_manager;
  int verbosity = 0;
  std::string user_config_file = "user.cfg";
  std::string main_config_file = "echidna.cfg";
  std::string config = "default";
  std::string run;

  do {		// do {} while to avoid too many local variables
    time_t t = time (0);
    char time_string[31];
    ::strftime(time_string, 31, "echidna_%d:%m:%Y_%H%M%S.log", ::localtime(&t));
    log_filename = time_string;
  } while (0); 

  opterr = 0;
  while (1) {
    const char valid_options[] = ":vqc:C:f:e:j:p:lL:Vo:O:";
    int c = ::getopt (argc, argv, valid_options);
    if (c == -1) break;
    
    switch (c) {
      case 'c':
	user_config_file = optarg;
	break;
      case 'C':
	main_config_file = optarg;
	break;
      case 'f':
	cm->add_comman_line_parameter ("bx_reader.file_url", optarg);
	run = optarg;
	break;
      case 'e':
	cm->add_comman_line_parameter ("bx_reco_framework.max_events", optarg);
	break;
      case 'j':
	cm->add_comman_line_parameter ("bx_reco_framework.skip_events", optarg);
	break;
      case 'p':
	  // optind is the index to the next option (so argv[optind] is the second opt parameter)
	if (!argv[optind] || argv[optind][0] == '-') 
	  for (size_t k = 0; k < ::strlen (valid_options); k++) {
	    if (valid_options[k] != ':' && valid_options[k] == argv[optind][0]) 
	      throw std::runtime_error ("option -p requires 2 parameters with the \"module.param_name value\" syntax"); 
	  }
	cm->add_comman_line_parameter (optarg, argv[optind++]);
	break;
      case 'l':
        if (optind < argc && argv[optind] && argv[optind][0] != '-') log_filename = argv[optind++];
	else log_filename = "";
	required_log_filename = true;
	break;
      case 'L':
	log_path = optarg;
	break;
      case 'V':
	messenger::get ()->set_line_buffered (true);
	break;
      case 'v':
	verbosity ++;
	break;
      case 'q':
	verbosity --;
	break;
      case 'o':
	cm->add_comman_line_parameter ("bx_writer.file_name", optarg);	
	break;
      case 'O':
	cm->add_comman_line_parameter ("bx_writer.directory", optarg);	
	break;
      case ':':
        throw std::runtime_error ("missing value on command line option -" + char (optopt));
	break;
      default:
        throw std::runtime_error ("unknown command line option: -" + char (optopt));
    }
  }

    // Check for extra option comman line fields
  if (optind < argc) {
      // First argument = config
    if ((optind + 1) == argc) config = argv[optind++];
    if (optind < argc) {
      std::ostringstream msg;
      msg << "unknown fields on the command line:";
      while (optind < argc) msg << " \"" << argv[optind++] << "\"";
      throw std::runtime_error (msg.str ());
    }
  }

    // Set verbosity level and log file
  if (verbosity == -1) bx_message::set_default_print_level (bx_message::error);
  else if (verbosity == 1) bx_message::set_default_print_level (bx_message::info);
  else if (verbosity == 2) bx_message::set_default_print_level (bx_message::log);
  else if (verbosity < -1) bx_message::set_default_print_level (bx_message::critic);
  else if (verbosity > 2) bx_message::set_default_print_level (bx_message::debug);
  else bx_message::set_default_print_level (bx_message::warn);
  if (!required_log_filename && run.size () > 6 && run.substr (0, 6) == "run://" && vdt (run.substr (6, run.size ())).get_type () == vdt::int_vdt) {
    int run_number = vdt (run.substr (6, run.size ())).get_int ();
    std::ostringstream new_name;
    if (log_path.size ()) new_name << log_path << "/";
    new_name << "Run" << std::setw (6) << std::setfill ('0') << run_number;
    if (config != "default") new_name << "_" << config;
    new_name << "_c" << CYCLE_NUMBER << ".log";
    log_filename = new_name.str ();
    required_log_filename = true;
  } else if (log_path.size ()) log_filename = log_path + "/" + log_filename;
  messenger::get ()->open_log_file (log_filename);
  bx_message msg(bx_message::log);
  msg << "Running Echidna version " << ECHIDNA_VERSION << " on logfile " << log_filename << " with: ";
  for (int i = 0; i < argc; i++) msg << argv[i] << " ";
  if (verbosity < 1) std::cout << msg.str () << std::endl;
  msg << dispatch;

    // Send config files to the configuration manager
  cm->add_main_config_file (main_config_file);
  if (user_config_file.size ()) cm->add_user_config_file (user_config_file);
  
    // Upload the parameter broker with parameters from the current config
  cm->upload_configuration (config);

    // Set the parameter broker to bx_named (indeed the configuration manager is a parameter broker)
  bx_named::set_parameter_broker (cm);

    // Log the configuration
  cm->log_configuration ();

    // Create modules
  bx_module_factory::get ()->create_modules (cm->get_module_name_list());
}

/*
 * $Log: bx_options.cc,v $
 * Revision 1.28  2008/08/05 14:38:41  razeto
 * Writing logfile on output
 *
 * Revision 1.27  2007-06-23 11:14:46  razeto
 * -O added to specify output directory
 *
 * Revision 1.26  2007-06-07 10:29:46  razeto
 * Add a printout on cout with echidna arguments (to track pbs logs)
 *
 * Revision 1.25  2007-06-07 10:10:55  razeto
 * Use final log file name when using run:// input file specification
 *
 * Revision 1.24  2007-05-09 20:27:34  razeto
 * Add cycle index to log (conforming to root files)
 *
 * Revision 1.23  2007-03-30 13:39:27  razeto
 * Added line buffering option (-V) for online echidna
 *
 * Revision 1.22  2007-02-07 13:15:02  razeto
 * Added -L to specify a path to logfiles. Added finish to move log file name to canonical name.
 *
 * Revision 1.21  2006/05/10 12:14:39  razeto
 * Removed some useless include
 *
 * Revision 1.20  2006/04/02 11:19:42  razeto
 * Added skip first event option (-j)
 *
 * Revision 1.19  2005/06/28 15:54:12  razeto
 * Fixed a check to allow negative paramenters on line command
 *
 * Revision 1.18  2005/03/04 11:50:45  razeto
 * Transformed the framework in a named to pass more options
 *
 * Revision 1.17  2004/11/26 15:25:10  razeto
 * Thanks to Davide who so kindly underlined a syntactical error, I fixed the "Maintainer" word
 *
 * Revision 1.16  2004/11/26 14:01:20  razeto
 * Added Mantainer field
 *
 * Revision 1.15  2004/10/06 09:43:36  razeto
 * Added -q to decrement verbosity
 *
 * Revision 1.14  2004/09/30 14:40:17  razeto
 * Added some log printout: echidna version, command line and current configuration
 *
 * Revision 1.13  2004/09/30 13:04:17  razeto
 * Added bx_message::log level handling
 *
 * Revision 1.12  2004/09/23 10:10:25  razeto
 * Updated reader and writer name
 *
 * Revision 1.11  2004/08/06 10:30:28  razeto
 * cycle_1 branch merged in the main trunk.
 *
 * Revision 1.10.2.1  2004/07/26 14:16:29  razeto
 * Fixed a output file name parameter name.
 *
 * Revision 1.10  2004/07/22 10:13:01  razeto
 * Jumbo patch for the configuration manager applied, read the documentation for usage
 *
 * Revision 1.9  2004/05/26 09:22:37  razeto
 * Upgraded to last version of vdt
 *
 * Revision 1.8  2004/05/26 08:27:15  razeto
 * Added max number of events (negative = unlimited)
 *
 * Revision 1.7  2004/04/24 17:40:22  razeto
 * Updated to the new bx_reader file access
 *
 * Revision 1.6  2004/04/18 09:16:52  razeto
 * Added -o option, to select the output root file
 *
 * Revision 1.5  2004/04/12 16:07:39  razeto
 * Added configurable default print level for messages
 *
 * Revision 1.4  2004/04/06 12:39:28  razeto
 * Added logfile handling
 *
 * Revision 1.3  2004/04/05 13:12:27  razeto
 * Updated a comment, changed the behaviour if no user option file is present
 *
 * Revision 1.2  2004/03/29 13:17:30  razeto
 * Updated bx_options to the vdt usage
 *
 * Revision 1.1  2004/03/26 16:36:11  razeto
 * Introduced bx_options; a lot of code modified to read options and parameters.
 * Bx_event_reader interface changed: now there is a standard constructor and
 * the file opening is done at begin using the parameters for the file name.
 *
 *
 */
