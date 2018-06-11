/* BOREXINO Reconstruction program
 *
 * Author: Davide D'Angelo <davide.dangelo@lngs.infn.it>
 * based on previous work by Razeto&Pallavicini
 * Maintainer: Davide D'Angelo <davide.dangelo@lngs.infn.it>
 *
 * $Id: bx_writer.hh,v 1.13 2006/08/21 11:25:33 razeto Exp $
 *
 * Special module to write ROOT file
 * 
 */

#ifndef _BX_WRITER_HH
#define _BX_WRITER_HH

#include "bx_base_module.hh"
#include "bx_rec_general.hh"

class TFile;
class TTree;
class bx_echidna_event;
class bx_write_opts;
class BxEvent;

// module
class bx_writer: public bx_base_module {
  public:
    bx_writer ();
    virtual ~bx_writer () {}

    virtual void begin ();
    virtual bx_echidna_event* doit (bx_echidna_event *ev);
    virtual void end ();
	
  private:
    void m_parse_options(bx_write_opts&);   // helper method to retrieve user parameters and fill an opt bitfield 
    TFile *p_root_file;       // ptr to TFile
    TTree *p_root_tree;       // ptr to TTree 
    BxEvent *p_root_event;    // ptr to BxEvent, used for TTree::Branch() 
    uint32_t u4_written_events;
};

#endif
/*
 * $Log: bx_writer.hh,v $
 * Revision 1.13  2006/08/21 11:25:33  razeto
 * Updated to new barn_interface + count written events
 *
 * Revision 1.12  2004/11/26 15:25:11  razeto
 * Thanks to Davide who so kindly underlined a syntactical error, I fixed the "Maintainer" word
 *
 * Revision 1.11  2004/11/26 14:06:23  razeto
 * Added Mantainer field
 *
 * Revision 1.10  2004/09/22 15:02:22  ddangelo
 * class name moved from bx_event_writer to bx_writer to match file name
 *
 * Revision 1.9  2004/09/22 13:58:59  ddangelo
 * updated bx_reco event into bx_echidna_event
 *
 * Revision 1.8  2004/06/07 12:57:23  ddangelo
 * fixed some dependencies
 *
 * Revision 1.7  2004/06/03 14:59:36  ddangelo
 * new/delete on BxEvent replaced with operator=().
 * This should solve the problem of TTree::Fill() caching the addresses of TClonesArray.
 * BxEvent still in tmp state.
 *
 * Revision 1.6  2004/05/30 11:58:14  ddangelo
 * Redisign to use new Root file classes.
 * Filler visitor removed.
 * User control over TTree::Branch() parameters introduced.
 * User parameter control over hit lists writing introduced.
 *
 * Revision 1.5  2004/04/18 09:15:57  razeto
 * Fixed indentation
 *
 * Revision 1.4  2004/04/16 15:02:22  pallas
 * Start developing of root event structure
 * Not defined yet, but preprared infrastructure
 *
 * Revision 1.3  2004/03/21 18:53:01  razeto
 * Some cosmetic changes
 *
 * Revision 1.2  2004/03/20 18:55:02  pallas
 * Debugging framework
 *
 * Revision 1.1  2004/03/20 17:39:45  pallas
 * First release
 *
 * Revision 1.1  2004/03/20 14:16:47  razeto
 * Added sources
 *
 */
