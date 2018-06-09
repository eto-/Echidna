/* BOREXINO Reconstruction program
 * 
 * Author: Alessandro Razeto <Alessandro.Razeto@ge.infn.it>
 *	   Marco Pallavicini <pallas@ge.infn.it>
 * Maintainer: Alessandro Razeto <Alessandro.Razeto@ge.infn.it>
 *
 * $Id: bx_rec_general.hh,v 1.28.2.1 2015/08/26 11:35:01 misiaszek Exp $
 *
 * A general file included by everything for macro definition
 *
 */

#ifndef _BX_REC_GENERAL_H
#define _BX_REC_GENERAL_H
// ECHIDNA_VERSION is a string indicating the version of the
// Project: versioning is based on cycles. Versions can be one
// of the following strings (for a cycle N):
// 	cycle_N_unstable  (the free developement before the cycle testing branch)
// 	cycle_N_testing	  (the version during of the testing branch)
// 	cycle_N_release   (the stable version, end of the branch)
#define CYCLE_NUMBER 18
#define ECHIDNA_VERSION "cycle_18_testing"


#include <math.h>
extern "C" {
  double round(double x);
  float roundf(float x);
};

#endif
/*
 * $Log: bx_rec_general.hh,v $
 * Revision 1.28.2.1  2015/08/26 11:35:01  misiaszek
 * cycle_18 testing
 *
 * Revision 1.28  2015/01/09 15:03:07  misiaszek
 * cycle_18 new unstable
 *
 * Revision 1.27  2013/06/18 18:41:46  razeto
 * cycle_17 new unstable
 *
 * Revision 1.26  2013-02-02 09:01:49  razeto
 * Incremented to cycle_16 (cycle 15 was lost)
 *
 * Revision 1.25  2011-04-19 05:54:58  razeto
 * Moved to cycle 15 unstable
 *
 * Revision 1.24  2011-04-19 05:51:53  razeto
 * Cycle 14 testing released
 *
 * Revision 1.23  2010-08-06 17:20:16  razeto
 * Moving to cycle 14
 *
 * Revision 1.22  2009-11-26 13:42:51  razeto
 * Moved to cycle_13_unstable
 *
 * Revision 1.21  2008-12-15 17:12:10  razeto
 * New cycle (12)
 *
 * Revision 1.20  2008-10-17 13:41:12  razeto
 * new development cycle (11)
 *
 * Revision 1.19  2008-02-27 20:46:13  razeto
 * new development cycle (10)
 *
 * Revision 1.18  2008-02-27 20:26:29  razeto
 * New clasdef(9) version and new cycle version
 *
 * Revision 1.17  2007-10-11 10:49:53  razeto
 * Cycle 8 deployed
 *
 * Revision 1.16  2007-06-22 15:10:47  razeto
 * Moved to cycle 7
 *
 * Revision 1.15  2007-03-29 13:51:30  razeto
 * Added numerical cycle
 *
 * Revision 1.14  2007-03-27 13:08:09  razeto
 * Set branch name
 *
 * Revision 1.13  2006-06-23 12:03:52  razeto
 * Updated version to new development cycle (5)
 *
 * Revision 1.12  2005/09/02 12:48:55  razeto
 * Set the cycle_4_unstable version tag
 *
 * Revision 1.11  2005/03/01 15:05:46  razeto
 * Merged with cycle_2
 *
 * Revision 1.10  2004/12/13 12:38:14  razeto
 * Added a roundf definitions since gcc 2.95 is missing it
 *
 * Revision 1.9  2004/12/03 14:53:20  razeto
 * New tag for cycle_3_unstable added
 *
 * Revision 1.8.2.3  2005/03/01 15:01:15  razeto
 * Version cycle 2 released
 *
 * Revision 1.8.2.2  2004/12/03 14:54:25  razeto
 * Fixed a typo
 *
 * Revision 1.8.2.1  2004/12/03 14:47:54  razeto
 * Updated tag
 *
 * Revision 1.8  2004/12/03 09:48:39  razeto
 * Changed status names to a more symmetric pattern
 *
 * Revision 1.7  2004/11/26 15:25:09  razeto
 * Thanks to Davide who so kindly underlined a syntactical error, I fixed the "Maintainer" word
 *
 * Revision 1.6  2004/11/26 14:01:18  razeto
 * Added Mantainer field
 *
 * Revision 1.5  2004/09/30 14:39:29  razeto
 * Added a string representing the version (updated at the current version) and a comment
 *
 * Revision 1.4  2004/09/28 13:43:51  razeto
 * Removed the ifdef for root barn histos
 *
 * Revision 1.3  2004/05/21 11:22:14  razeto
 * Added a macro used by laben precalib modules. This fixes a bug
 *
 * Revision 1.2  2004/03/20 12:25:16  razeto
 * Added comment header and footer. Added ifdef
 *
 * Revision 1.1.1.1  2004/03/19 18:23:49  razeto
 * Imported sources
 *
 */
