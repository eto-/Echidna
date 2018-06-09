/* BOREXINO Reconstruction program
 *
 * Author: Stefano Davini (stefano.davini@ge.infn.it)
 * Maintainer: Stefano Davini (stefano.davini@ge.infn.it)
 *
 * $Id: bx_v1731sys.hh,v 1.6 2008/12/15 13:37:16 davini Exp $
 *
 * Module for analysis of V1731 data (28 channels system for neutron detection and other)
 * 
*/

#ifndef _BX_V1731SYS_HH
#define _BX_V1731SYS_HH

#include "bx_rec_general.hh"
#include "bx_base_module.hh"
#include "bx_neutron_event.hh"
#include "v1731_event_correlator.hh"
#include "v1731_event.hh"
#include "v1731_neutron_analyzer.hh"
#include "zlib.h"

class bx_echidna_event;

class bx_v1731sys : public bx_base_module {
public:
  bx_v1731sys ();
  virtual ~bx_v1731sys () {}
  
  // Operations for the framework 
  virtual void begin ();
  virtual bx_echidna_event* doit (bx_echidna_event *ev);
  virtual void end ();
  
private:
  v1731_event_correlator::transfer transfer_opt;  // http or file
  v1731_event_correlator::mapper_time_file* map_time_file;
  v1731_event_correlator::v1731event_indexer v1731event_index;
  bool b_map_created;
  gzFile      currNfile;            // used in v1731_event_getter
  FILE*       p_connect;            // used in v1731 event_getter
  std::string s_curr_file_name;     
  long        i4_curr_file_position;
  unsigned long last_neutron_gpsstime;  // used for run w0 trgID aq
  long          last_neutron_file_position;
  int           last_neutron_laben_nclusters;
  void m_load_file_map();
  void m_clear_ev(bx_neutron_event&);
  void m_write_ev(bx_neutron_event&, const v1731_neutron_analyzer::neutron_analyzer&);
};

#endif
/*
 * $Log: bx_v1731sys.hh,v $
 * Revision 1.6  2008/12/15 13:37:16  davini
 * event index
 *
 * Revision 1.5  2008-11-27 14:18:42  davini
 * addes vars and lines for sequential readout; gzFile and pipe are private vars;
 *
 * Revision 1.4  2008-08-06 14:58:56  davini
 * added TRGID alignement
 *
 * Revision 1.3  2008-07-29 17:04:34  davini
 * improved alignement
 *
 * Revision 1.2  2008-07-14 09:28:31  davini
 * module ready
 *
 * Revision 1.1  2008-04-02 13:39:42  pallas
 * New module for the analysis of the new neutron system data
 * Developer and maintainer: S. Davini
 * Empty module now. Do NOT use until further notice.
 *
 *
 */
