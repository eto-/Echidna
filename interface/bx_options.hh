/* BOREXINO Reconstruction program
 *
 * Author: Alessandro Razeto <Alessandro.Razeto@ge.infn.it>
 * Maintainer: Alessandro Razeto <Alessandro.Razeto@ge.infn.it>
 *
 * $Id: bx_options.hh,v 1.12 2007/02/07 13:15:02 razeto Exp $
 *
 * Options interface: bx_options::initialize_echidna is the first line of
 * code that should be called in echidna, since initializes the parameter
 * broker (aka parameter and configuration management), the messenger (log
 * file and print level) and the module factory.
 * 
 *
 */
#ifndef _BX_OPTIONS_H
#define _BX_OPTIONS_H
#include <string>
#include <vector>

class bx_options {
  public:
    static void initialize_echidna (int argc, char *const argv[]);
    static void finish ();

  private:
    static bx_options *me;

    bx_options (int argc, char *const argv[]);
    bool required_log_filename;
    std::string log_filename, log_path;
};
 
#endif
/*
 * $Log: bx_options.hh,v $
 * Revision 1.12  2007/02/07 13:15:02  razeto
 * Added -L to specify a path to logfiles. Added finish to move log file name to canonical name.
 *
 * Revision 1.11  2005/08/02 14:01:49  razeto
 * Removed some useless include
 *
 * Revision 1.10  2005/03/04 11:50:45  razeto
 * Transformed the framework in a named to pass more options
 *
 * Revision 1.9  2004/11/26 15:25:10  razeto
 * Thanks to Davide who so kindly underlined a syntactical error, I fixed the "Maintainer" word
 *
 * Revision 1.8  2004/11/26 14:01:20  razeto
 * Added Mantainer field
 *
 * Revision 1.7  2004/07/22 10:13:01  razeto
 * Jumbo patch for the configuration manager applied, read the documentation for usage
 *
 * Revision 1.6  2004/06/10 12:40:01  razeto
 * Removed an unuseful include
 *
 * Revision 1.5  2004/05/26 08:27:15  razeto
 * Added max number of events (negative = unlimited)
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
