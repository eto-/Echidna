/* BOREXINO Reconstruction program
 * 
 * Author: Davide D'Angelo <davide.dangelo@lngs.infn.it>
 * Maintainer: Davide D'Angelo <davide.dangelo@lngs.infn.it>
 *
 * $Id: Mach4EventLinkDef.hh,v 1.2 2010/03/09 13:07:35 razeto Exp $
 *
 * CINT LinkDef for event 
 *
 */

#ifndef _MACH4EVENTLINKDEF_HH
#define _MACH4EVENTLINKDEF_HH
#if defined(__ROOTCLING__) || defined(__CINT__)
#define _ROOT_CINT_
#endif

#ifdef _ROOT_CINT_

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ class Mach4Event;
#pragma link C++ class std::vector<Mach4Event>;

#endif
#endif
/*
 * $Log: Mach4EventLinkDef.hh,v $
 * Revision 1.2  2010/03/09 13:07:35  razeto
 * Added dictionary for vector<Mach4Event>
 *
 * Revision 1.1  2010-03-08 20:38:23  razeto
 * Added merging tools
 *
 *
 */
