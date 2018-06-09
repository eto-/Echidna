/* BOREXINO Reconstruction program
 *
 * Author: Alessandro Razeto <Alessandro.Razeto@ge.infn.it>
 * Maintainer: Alessandro Razeto <Alessandro.Razeto@ge.infn.it>
 *
 * $Id: db_calib.hh,v 1.20 2015/01/09 15:03:08 misiaszek Exp $
 *
 * The database interface for calib profile dependent informations
 *
 */
#ifndef _BD_CALIB_H
#define _BD_CALIB_H
#ifndef CYCLE_NUMBER
#define CYCLE_NUMBER 18
#endif
#include "db_acl.hh"

#include <TObject.h>

#include <string>
#include <map>
#include <vector>
#include <time.h>

class db_calib: public db_acl, public TObject {
  public:
    db_calib() : db_acl(), TObject() {};
    float get_refraction_index_data 	() const { return f4_refraction_index_data; }
    float get_refraction_index_mc 	() const { return f4_refraction_index_mc; }

    typedef std::vector<float> pdf_v;
    const pdf_v& get_pdf_time_vector	() const { return pdf_time_v; }
    const pdf_v& get_pdf_mc_vector 	() const { return pdf_mc_v; }
    const pdf_v& get_pdf_data_vector 	() const { return pdf_data_v; }
  
  private:
    float f4_refraction_index_data, f4_refraction_index_mc;
    pdf_v pdf_time_v, pdf_mc_v, pdf_data_v;

      // Private ctors and dctor
    db_calib (int calib_profile);

      // Only bx_dbi can istantiate and destroy db_calib objects
    friend class bx_dbi;
    ClassDef(db_calib,CYCLE_NUMBER)
};
#endif
/*  
 *  $Log: db_calib.hh,v $
 *  Revision 1.20  2015/01/09 15:03:08  misiaszek
 *  cycle_18 new unstable
 *
 *  Revision 1.19  2013/06/18 18:56:37  razeto
 *  cycle_17 new unstable
 *
 *  Revision 1.18  2013-02-02 09:01:49  razeto
 *  Incremented to cycle_16 (cycle 15 was lost)
 *
 *  Revision 1.17  2011-04-19 05:54:58  razeto
 *  Moved to cycle 15 unstable
 *
 *  Revision 1.16  2010-08-06 17:20:16  razeto
 *  Moving to cycle 14
 *
 *  Revision 1.15  2009-11-26 13:42:51  razeto
 *  Moved to cycle_13_unstable
 *
 *  Revision 1.14  2008-12-15 17:13:55  razeto
 *  New cycle (12)
 *
 *  Revision 1.13  2008-10-17 13:41:12  razeto
 *  new development cycle (11)
 *
 *  Revision 1.12  2008-02-27 20:46:13  razeto
 *  new development cycle (10)
 *
 *  Revision 1.11  2008-02-27 20:26:30  razeto
 *  New clasdef(9) version and new cycle version
 *
 *  Revision 1.10  2007-10-11 10:49:54  razeto
 *  Cycle 8 deployed
 *
 *  Revision 1.9  2007-06-22 15:15:26  razeto
 *  Moved to cycle 7
 *
 *  Revision 1.8  2007-05-07 15:47:31  razeto
 *  Cycle number in root classdef
 *
 *  Revision 1.7  2006/11/05 10:27:30  razeto
 *  Remove some useless include
 *
 *  Revision 1.6  2006/01/25 13:24:05  misiaszek
 *  Moved from cmap to simple map (to work with root)
 *
 *  Revision 1.5  2004/11/26 15:25:10  razeto
 *  Thanks to Davide who so kindly underlined a syntactical error, I fixed the "Maintainer" word
 *
 *  Revision 1.4  2004/11/26 14:01:20  razeto
 *  Added Mantainer field
 *
 *  Revision 1.3  2004/11/24 09:46:41  razeto
 *  Moved some prototyping in db_acl from db_*
 *
 *  Revision 1.2  2004/09/29 10:04:01  razeto
 *  Changed the interface to return vector for pdf
 *
 *  Revision 1.1  2004/09/16 11:44:48  razeto
 *  Added
 *
 */
