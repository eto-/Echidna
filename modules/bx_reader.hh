/* BOREXINO Reconstruction program
 *
 * Author: Alessandro Razeto <Alessandro.Razeto@ge.infn.it>
 * Maintainer: Alessandro Razeto <Alessandro.Razeto@ge.infn.it>
 *
 * $Id: bx_reader.hh,v 1.18 2006/11/06 18:13:27 razeto Exp $
 *
 * This is the pure base class for the event hierarchy
 *
 */
#ifndef _BX_READER_H
#define _BX_READER_H

#include "bx_echidna_event.hh"
#include "bx_base_module.hh"
#include "bx_rec_general.hh"

#include <map>
#include <string>
#include <stdio.h>
#include <zlib.h>

class bx_reader: public bx_base_module {
  public:
    bx_reader ();
    virtual ~bx_reader () {}

    virtual void begin ();
    virtual bx_echidna_event* doit (bx_echidna_event *ev) { return get_event (); }
    virtual void end ();
    
    bx_echidna_event *get_event (int32_t event_number = -1, uint32_t trg_type = 0);
    void rewind () { current_event = disk_pool.begin (); }

  private:
    uint32_t u4_run_number;
    gzFile gzfile;
    void m_open_gzfile (const std::string &filename);
    FILE *p_popen_connection;
    bool first_event;

    uint32_t u4_pool_size, u4_max_pool_size, u4_pool_max_event_count;
    uint16_t u2_pool_readhaed_step;
    int32_t i4_last_read_evnum;
    
    typedef std::map<int32_t, char*> disk_event_pool;
    disk_event_pool disk_pool;

    disk_event_pool::const_iterator current_event;
    char *m_feed_event (int32_t event_number);
    char *m_gzfile_exception ();
};

#endif
/*
 * $Log: bx_reader.hh,v $
 * Revision 1.18  2006/11/06 18:13:27  razeto
 * Zero event is no more lost in war
 *
 * Revision 1.17  2005/06/29 12:07:59  razeto
 * Upgraded the reader to do some printing
 *
 * Revision 1.16  2004/11/26 15:25:11  razeto
 * Thanks to Davide who so kindly underlined a syntactical error, I fixed the "Maintainer" word
 *
 * Revision 1.15  2004/11/26 14:01:21  razeto
 * Added Mantainer field
 *
 * Revision 1.14  2004/09/30 14:38:33  razeto
 * Added some cleanup and a warning in bx_reader::end
 *
 * Revision 1.13  2004/09/23 10:09:01  razeto
 * Changed internal reader name (to follow writer)
 *
 * Revision 1.12  2004/09/22 13:28:10  razeto
 * Changed bx_reco_event to bx_echidna_event
 *
 * Revision 1.11  2004/05/18 14:28:37  razeto
 * Removed unused stuff
 *
 * Revision 1.10  2004/04/18 10:31:38  razeto
 * Updated
 *
 * Revision 1.9  2004/04/12 15:57:20  razeto
 * Added rewind method, to be used after each precalib loop
 *
 * Revision 1.8  2004/04/01 12:09:55  razeto
 * Moved a typedef (used only internally to the class) in bx_event_reader
 *
 * Revision 1.7  2004/03/26 16:36:11  razeto
 * Introduced bx_options; a lot of code modified to read options and parameters.
 * Bx_event_reader interface changed: now there is a standard constructor and
 * the file opening is done at begin using the parameters for the file name.
 *
 * Revision 1.6  2004/03/22 14:57:12  razeto
 * Updated code to avoid g++ -Wall warnings
 *
 * Revision 1.5  2004/03/22 13:44:18  razeto
 * Upgraded the handling of pool and introduced some other fixes
 *
 * Revision 1.4  2004/03/21 18:55:46  razeto
 * Some cosmetic changes
 *
 * Revision 1.3  2004/03/20 18:00:38  razeto
 * Added some virtual methods: TO BE completed
 *
 * Revision 1.2  2004/03/20 17:38:46  razeto
 * Added a first working object
 *
 * Revision 1.1  2004/03/20 14:16:47  razeto
 * Added sources
 *
 */
