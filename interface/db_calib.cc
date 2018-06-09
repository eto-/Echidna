/* BOREXINO Reconstruction program
 *
 * Author: Alessandro Razeto <Alessandro.Razeto@ge.infn.it>
 * Maintainer: Alessandro Razeto <Alessandro.Razeto@ge.infn.it>
 *
 * $Id: db_calib.cc,v 1.7 2007/11/09 22:52:35 razeto Exp $
 *
 * The database interface for run informations
 *
 */
#include "db_calib.hh"
#include "bx_dbi.hh"
#include "messenger.hh"

#include <sstream>

ClassImp(db_calib)

#ifndef _ECHIDNA_ROOTLIB_

db_calib::db_calib (int calib_profile): db_acl(), TObject() {

  bx_dbi *dbi = bx_dbi::get ();
  bx_message &msg = bx_dbi::get ()->get_message (bx_message::critic);
  msg << "db_calib: ";

  std::ostringstream where_str;
  where_str << "\"CalibProfile\"=" << calib_profile;
  const bx_dbi::table &table = dbi->query (bx_dbi::bx_calib, "\"ScintillatorParameters\"", where_str.str (), "*", "", -1);

  f4_refraction_index_data = table["SourceDataRefidx"][0].get_float ();
  f4_refraction_index_mc = table["MonteCarloRefidx"][0].get_float ();

  const bx_dbi::table &table2 = dbi->query (bx_dbi::bx_calib, "\"PDFHistogram\"", where_str.str (), "*", "", 0);
  const bx_dbi::column& pdf_time_c = table2["Time"];
  const bx_dbi::column& pdf_mc_c = table2["MonteCarloPDF"];
  const bx_dbi::column& pdf_data_c = table2["SourceDataPDF"];
  for (unsigned long i = 0; i < pdf_time_c.size (); i++) {
    pdf_time_v.push_back (pdf_time_c[i].get_float ());
    pdf_mc_v.push_back (pdf_mc_c[i].get_float ());
    pdf_data_v.push_back (pdf_data_c[i].get_float ());
  }

  bx_dbi::get ()->close_db_connections ();
  bx_dbi::get ()->get_message (bx_message::info) << "db_calib: db_calib (" << calib_profile << ") initialized" << dispatch;
}

#endif

/*  
 *  $Log: db_calib.cc,v $
 *  Revision 1.7  2007/11/09 22:52:35  razeto
 *  Fixed a printout
 *
 *  Revision 1.6  2007-06-03 16:13:30  razeto
 *  Do not keep db connection opened (else bxdb has lots of problems)
 *
 *  Revision 1.5  2006-01-25 13:24:05  misiaszek
 *  Moved from cmap to simple map (to work with root)
 *
 *  Revision 1.4  2004/11/26 15:25:10  razeto
 *  Thanks to Davide who so kindly underlined a syntactical error, I fixed the "Maintainer" word
 *
 *  Revision 1.3  2004/11/26 14:01:20  razeto
 *  Added Mantainer field
 *
 *  Revision 1.2  2004/09/29 09:59:25  razeto
 *  Fixed some bugs in readout
 *
 *  Revision 1.1  2004/09/16 11:44:48  razeto
 *  Added
 *
 */
