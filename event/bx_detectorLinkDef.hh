/* BOREXINO Reconstruction program
 * 
 * Author:  Alessandro Razeto <Alessandro.Razeto@ge.infn.it>
 * Maintainer:  Alessandro Razeto <Alessandro.Razeto@ge.infn.it>
 *
 * $Id: bx_detectorLinkDef.hh,v 1.3 2008/09/26 09:56:05 razeto Exp $
 *
 * CINT LinkDef for bx_detector
 *
 */

#ifndef _BX_DETECTOR_LINKDEF_HH
#define _BX_DETECTOR_LINKDEF_HH
#if defined(__ROOTCLING__) || defined(__CINT__)
#define _ROOT_CINT_
#endif

#ifdef _ROOT_CINT_

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ class bx_detector;
#if ROOT_VERSION_CODE < 327680
#pragma link C++ class std::map<std::string, int>;
#pragma link C++ class std::vector<int>;
#endif

#endif
#endif
/*
 * $Log: bx_detectorLinkDef.hh,v $
 * Revision 1.3  2008/09/26 09:56:05  razeto
 * Do not create standard containre linkage for new root
 *
 * Revision 1.2  2006-08-21 14:42:06  razeto
 * Added few types for map and vector
 *
 * Revision 1.1  2006/08/21 11:11:59  razeto
 * Improved bx_detector and created detector_interface
 *
 */
