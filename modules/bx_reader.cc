/* BOREXINO Reconstruction program
 *
 * Author: Alessandro Razeto <Alessandro.Razeto@ge.infn.it>
 * Maintainer: Alessandro Razeto <Alessandro.Razeto@ge.infn.it>
 *
 * $Id: bx_reader.cc,v 1.47 2012/04/02 11:40:37 razeto Exp $
 *
 * Implementation of bx_reader.
 *
 */
#include "bx_reader.hh"
#include "bx_event_disk_format.h"
#include "bx_dbi.hh"
#include "db_run.hh"
#include "vdt.hh"
#include "bx_detector.hh"

#include <iostream>
#include <string>
#include <sstream>
#include <stdexcept>
#include <unistd.h>
#include <errno.h>
#include <stdio.h>
#include <time.h>
#include <string.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <string.h>

void sighandler (int32_t sig) {
  bx_message msg (bx_message::critic, "bx_reader: ");
  int32_t pid, status;
  if ((pid = waitpid (-1, &status, WUNTRACED)) <0) {
    bx_message msg (bx_message::critic, "bx_reader: ");
    msg << "waitpid returned " << strerror (errno) << dispatch;
  }
  
  if (WIFSIGNALED (status)) msg << "Process " << pid << " signaled with " << WTERMSIG(status) << dispatch;
  else {
    if (!WEXITSTATUS (status)) msg.set_level (bx_message::log);
    msg << "Process " << pid << " exit with " << WEXITSTATUS(status) << dispatch;
  }
}  

bx_reader::bx_reader (): bx_base_module ("bx_reader", bx_base_module::reader), u4_run_number(0), p_popen_connection(0), first_event(true) {
}

void bx_reader::m_open_gzfile (const std::string &file_url) {
  
  std::string::size_type pos = file_url.find ("://");
  
  std::string name, method;
  if (pos == std::string::npos) {
    name = file_url;
    method = "file";
  } else {
    name = file_url.substr (pos + 3);
    method = file_url.substr (0, pos);
  }

  if (method == "file") {  
    if (::access (name.c_str (), R_OK) < 0) 
      get_message (bx_message::critic) << "file \"" << name << "\" access: " << ::strerror (errno) << dispatch;
    gzfile = ::gzopen (name.c_str (), "r");
    if (!gzfile) get_message (bx_message::critic) << "gzopen: Z_MEM_ERROR" << dispatch; 
  } else if (method == "run") {
    if (vdt(name).get_type () != vdt::int_vdt) 
      get_message (bx_message::critic) << "run method requires an integer and \"" << name << "\" is not" << dispatch; 
    
      // Fill the url using the repository_url and a fixed rule for repository hierarchy 
        // Get repository url
    std::string repository_url = get_parameter ("repository_url").get_string ();
    
        // Calculate path in the repository
    int32_t run_number = vdt(name).get_int ();
    bx_dbi::get ()->set_current_run_number (u4_run_number = run_number);
    const db_run& run_info = bx_dbi::get ()->get_run ();
    time_t start_time = run_info.get_start_time ();
    struct tm *start_date = ::localtime (&start_time);
    time_t week_second = start_time - ((start_date->tm_wday * 24 + start_date->tm_hour) * 60 + start_date->tm_min) * 60 + 4800;
    struct tm *week_date = ::localtime (&week_second);
    std::string full_url;
    if (!run_info.get_number_of_files ()) get_message (bx_message::critic) << "No files for run " << run_number << dispatch;
    for (int32_t slice = 1; slice <= run_info.get_number_of_files (); slice++) {
      char tmp_str[100];
      ::strftime (tmp_str, 99, "/rawdata/%Y/%b_%d/", week_date);
      std::string file_path = tmp_str;

        // Get the file name
      std::ostringstream f_name;
      f_name << "Run" << std::setfill('0') << std::setw (6) << run_number << "_" << std::setfill('0') << std::setw (2) << slice << ".out.gz";

      // Get the url
      std::string url = repository_url + file_path + f_name.str ();
      if (slice > 1) full_url += " ";
      full_url += url;
      ::strftime (tmp_str, 99, "%a %d %b %Y", start_date);
      get_message (bx_message::info) << "run " << run_number << " taken on " << tmp_str << " is at \"" << url << "\"" << dispatch;
    }
    m_open_gzfile (full_url);
  } else if (method == "http") {
    std::string command = get_parameter ("wget_command").get_string () + " " + file_url;
    get_message (bx_message::log) << "executing " << command << dispatch;
    if (signal (SIGCHLD, sighandler) < 0) get_message (bx_message::critic) << "signal returned " << strerror (errno) << dispatch;
    p_popen_connection = ::popen(command.c_str (), "r");
    if (!p_popen_connection) get_message (bx_message::critic) << "error opening wget pipe" << dispatch;
    gzfile = ::gzdopen (::fileno (p_popen_connection), "r");
    if (!gzfile) get_message (bx_message::critic) << "gzopen: Z_MEM_ERROR" << dispatch; 
  } else {
    get_message (bx_message::critic) << "unknun file access method \"" << method <<"\"" << dispatch;
  }
}

void bx_reader::begin () {
  if (!check_parameter ("file_url")) 
    get_message (bx_message::critic) << "file_url parameter not specified" << dispatch;
  else m_open_gzfile (get_parameter ("file_url").get_string ());

  // Some default constant, later to be read from options
  u4_max_pool_size = get_parameter ("max_pool_size").get_int ();
  u4_pool_max_event_count = get_parameter ("pool_max_event_count").get_int ();
  u2_pool_readhaed_step = get_parameter ("pool_readhaed_step").get_int ();

  u4_pool_size = 0;
  current_event = disk_pool.end ();
  int32_t i = 0;
  for (;i < u2_pool_readhaed_step; i++) if (!m_feed_event (-1)) break; 
  if (!i) get_message (bx_message::critic) << "empty file" << dispatch; 

  current_event = disk_pool.begin ();

  detector_interface::init (*this);
  bx_echidna_event e(current_event->second);

  if (!u4_run_number) bx_dbi::get ()->set_current_run_number (u4_run_number = e.get_run_number ()); 
  if (e.get_run_number () != u4_run_number) get_message (bx_message::critic) << "error initialization failed: set " << u4_run_number << " file is " << e.get_run_number () << dispatch;

  bx_message &msg =  get_message (bx_message::info);
  msg << "echidna initialized with run number " << u4_run_number << newline;
  detector_interface::post_init (e, *this); 
}

bx_echidna_event *bx_reader::get_event (int32_t event_number, uint32_t trg_type) {
  if (event_number < -1) get_message (bx_message::critic) << "required negative event " << event_number << dispatch;
  
    // If the required event is no more in the pool throw an exception
  if (event_number != -1 && event_number < disk_pool.begin ()->first) 
    get_message (bx_message::critic) << "required event " << event_number << " while spooled events start at " << disk_pool.begin ()->first << dispatch;

    // Look for the next event
  disk_event_pool::const_iterator past_event = current_event;
  if (event_number == -1) {
    if (first_event) first_event = false; 
    else current_event ++; 
  } else current_event = disk_pool.find (event_number);

    // If the event is not present load it from file
  if (current_event == disk_pool.end ()) {
     if (!m_feed_event (event_number)) return 0;
       // Some readhaed
     for (uint16_t i = 0; i < u2_pool_readhaed_step; i++) if (!m_feed_event (-1)) break;

      // Now the event should be in memory pool, look for it again
    if (event_number == -1) current_event = ++past_event;
    else current_event = disk_pool.find (event_number);
    if (current_event == disk_pool.end ()) 
      get_message(bx_message::critic) << "internal error: event " << event_number << " not found in the pool" << dispatch;
  }

    // Strip the pool if it excedes size (both event counts and total ram size)
  while (disk_pool.size () > u4_pool_max_event_count && u4_pool_size > u4_max_pool_size) {
    disk_event_pool::iterator to_delete = disk_pool.begin ();
      // Stop if the present event is encountered  
    if (to_delete->first == current_event->first) break;
    //get_message (bx_message::debug) << "Deleting " << to_delete->first << dispatch;
    u4_pool_size -= ((event_header_disk_format *)(to_delete->second))->event_size_bytes;
    delete [] to_delete->second;
    disk_pool.erase (to_delete);
  }
  
  if (u4_run_number != ((event_header_disk_format *)(current_event->second))->run_number)
    get_message (bx_message::critic) << "echidna only support fixed run number" << dispatch;
  
  return new bx_echidna_event (current_event->second);
}
     
char *bx_reader::m_gzfile_exception () {
  gzgetc (gzfile); // This way a read is forced to check for eof
  if (::gzeof (gzfile)) return 0;
  
  bx_message &message = get_message (bx_message::critic);

  int32_t error_code;
  const char *msg = ::gzerror (gzfile, &error_code);
  
  if (error_code != Z_ERRNO) message << "gzlib: " << msg;
  else message << "gzlib: " << ::strerror (errno);

  message << dispatch;

  return 0;
}


char *bx_reader::m_feed_event (int32_t event_number) {
    // A dedicated filed i4_last_read_evnum is necessary since in case events are skiped 
    // (spool.end() - 1)->first is not always the last read event.
  //if (event_number != -1 && event_number < i4_last_read_evnum) get_message (bx_message::debug) <<  "back seek not supported " << event_number << " minor than last read event " << i4_last_read_evnum << dispatch;
  
    // Read the event header
  event_header_disk_format head;
  if (::gzread (gzfile, &head, sizeof (head)) != sizeof (head)) return m_gzfile_exception ();

    // Seek forward looking for requested event if event_number != -1
  if (event_number != -1) {
    if (head.event_number < (uint32_t)event_number) {
      if (::gzseek (gzfile, head.event_size_bytes - sizeof (head), SEEK_CUR) < 0) return m_gzfile_exception ();
      return m_feed_event (event_number);
    } else if (head.event_number > (uint32_t)event_number) {
        // Seeked too much, requested event not present
      if (::gzseek (gzfile, head.event_size_bytes - sizeof (head), SEEK_CUR) < 0) return m_gzfile_exception ();
        // Update the last event info, this means the next event could be read is > i4_last_read_evnum
      i4_last_read_evnum = head.event_number;
      return 0;
    }
  }
  
    // This point is reached is the event on disk is the event to read
  char * ptr = new char[head.event_size_bytes];
  if (::gzread (gzfile, ptr + sizeof (head), head.event_size_bytes - sizeof (head)) != int32_t(head.event_size_bytes - sizeof (head))) {
    delete [] ptr;
    return m_gzfile_exception ();
  }
    // Copy header to the buffer too
  *(event_header_disk_format *)ptr = head;

    // Update some internal info
  i4_last_read_evnum = head.event_number;
  u4_pool_size += head.event_size_bytes;

    // Update pool hash
  disk_pool[head.event_number] = ptr;
  
  //if (! (head.event_number % 100)) get_message (bx_message::debug) << "pool event " << disk_pool.size () << " size " << u4_pool_size << dispatch;
  return ptr;
}

void bx_reader::end () {
  if (p_popen_connection) {
    get_message(bx_message::info) << "closing network data stream, a transport layer warning in the next line is normal" << dispatch;
    if (signal (SIGCHLD, SIG_DFL) < 0) get_message (bx_message::critic) << "signal returned " << strerror (errno) << dispatch;
    ::pclose (p_popen_connection);
  }
  get_message(bx_message::log) << "processed up to event " << current_event->first << dispatch;
}
/*
 * $Log: bx_reader.cc,v $
 * Revision 1.47  2012/04/02 11:40:37  razeto
 * gzgetc is a macro
 *
 * Revision 1.46  2011-03-21 13:25:02  razeto
 * Set bx_dbi current_run_number (fix mc bug)
 *
 * Revision 1.45  2011-03-08 19:04:22  razeto
 * wget error ar citical
 *
 * Revision 1.44  2010-08-18 09:10:20  razeto
 * Added a check on wget exit
 *
 * Revision 1.43  2010-01-08 16:05:03  razeto
 * Fixed DST directory naming in storage
 *
 * Revision 1.42  2009-10-23 12:40:46  razeto
 * Muted
 *
 * Revision 1.41  2008-06-20 16:21:10  razeto
 * Added an include to compile with gcc 4.3
 *
 * Revision 1.40  2007-09-30 21:26:12  razeto
 * Fixed not sequential events bug
 *
 * Revision 1.39  2007-03-10 15:18:26  ddangelo
 * Now reader init the bx_detector (even with the event detector status).
 *
 * Revision 1.38  2007-02-06 18:02:25  razeto
 * New mapping of weeks, to follow bx_repository
 *
 * Revision 1.37  2006/11/06 18:13:27  razeto
 * Zero event is no more lost in war
 *
 * Revision 1.36  2006/10/22 23:58:49  razeto
 * Fixed end of file (for broken data) handling
 *
 * Revision 1.35  2006/03/21 14:22:52  razeto
 * Added a check on the number of slices
 *
 * Revision 1.34  2006/01/06 13:01:44  razeto
 * Now the reader download all run slices (with run:// method)
 *
 * Revision 1.33  2005/06/29 12:07:59  razeto
 * Upgraded the reader to do some printing
 *
 * Revision 1.32  2005/06/27 16:13:52  razeto
 * Added few parameter in echidna.cfg
 *
 * Revision 1.31  2004/12/22 15:49:50  razeto
 * Fixed a bug in message error level
 *
 * Revision 1.30  2004/11/26 16:20:46  razeto
 * Changed the level of a message
 *
 * Revision 1.29  2004/11/26 15:25:11  razeto
 * Thanks to Davide who so kindly underlined a syntactical error, I fixed the "Maintainer" word
 *
 * Revision 1.28  2004/11/26 14:01:21  razeto
 * Added Mantainer field
 *
 * Revision 1.27  2004/09/30 14:38:33  razeto
 * Added some cleanup and a warning in bx_reader::end
 *
 * Revision 1.26  2004/09/30 12:53:45  razeto
 * Changed a printout levele from info to log
 *
 * Revision 1.25  2004/09/27 14:20:57  razeto
 * Fixed a bug, added a debug printout
 *
 * Revision 1.24  2004/09/23 10:09:01  razeto
 * Changed internal reader name (to follow writer)
 *
 * Revision 1.23  2004/09/22 13:28:10  razeto
 * Changed bx_reco_event to bx_echidna_event
 *
 * Revision 1.22  2004/09/15 15:01:07  razeto
 * Fixed the english syntax to match the Davide's suggestions :-D
 *
 * Revision 1.21  2004/08/11 13:47:07  razeto
 * Updated to follow the new repository convention (month_day)
 *
 * Revision 1.20  2004/08/05 09:07:29  razeto
 * Upgraded reader to read from repository (even over network) using run:// url syntax
 *
 * Revision 1.19  2004/05/26 09:24:32  razeto
 * Upgraded to set the run number to db during begin (usefull for any db operation in modules::begin)
 *
 * Revision 1.18  2004/05/18 14:28:37  razeto
 * Removed unused stuff
 *
 * Revision 1.17  2004/04/27 12:14:48  ddangelo
 * updated to match bx_event_disk_format.h variable name changes
 *
 * Revision 1.16  2004/04/26 13:49:54  razeto
 * Added bx_run include
 *
 * Revision 1.15  2004/04/24 17:39:08  razeto
 * Updated the file access method with url like sintax. Still to be developed
 *
 * Revision 1.14  2004/04/18 10:31:38  razeto
 * Updated
 *
 * Revision 1.13  2004/04/13 14:52:15  razeto
 * Added run number setting for bx_dbi
 *
 * Revision 1.12  2004/04/06 12:42:17  razeto
 * Added module role to base_module, and the code to handle everywhere. Some other minor fixes
 *
 * Revision 1.11  2004/04/05 13:07:45  razeto
 * Added messanger for exception handling. Updated a comment
 *
 * Revision 1.10  2004/04/03 09:20:35  razeto
 * Added messenger
 *
 * Revision 1.9  2004/03/29 13:15:40  razeto
 * Updated bx_base_module to the vtd usage
 *
 * Revision 1.8  2004/03/26 16:36:11  razeto
 * Introduced bx_options; a lot of code modified to read options and parameters.
 * Bx_event_reader interface changed: now there is a standard constructor and
 * the file opening is done at begin using the parameters for the file name.
 *
 * Revision 1.7  2004/03/22 14:57:11  razeto
 * Updated code to avoid g++ -Wall warnings
 *
 * Revision 1.6  2004/03/22 13:44:18  razeto
 * Upgraded the handling of pool and introduced some other fixes
 *
 * Revision 1.5  2004/03/20 19:07:27  razeto
 * Fixed a bug
 *
 * Revision 1.4  2004/03/20 18:48:23  razeto
 * Fixed a bug
 *
 * Revision 1.3  2004/03/20 17:47:28  razeto
 * Using corrected bx_base_module constructor
 *
 * Revision 1.2  2004/03/20 17:38:46  razeto
 * Added a first working object
 *
 * Revision 1.1  2004/03/20 14:16:47  razeto
 * Added sources
 *
 */
