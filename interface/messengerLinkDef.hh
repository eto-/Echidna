/* BOREXINO physics tools
 *
 * linkdef definition for messages
 */
#ifndef _MESSENGERLINKDEF_HH
#define _MESSENGERLINKDEF_HH
#if defined(__ROOTCLING__) || defined(__CINT__)
#define _ROOT_CINT_
#endif

#ifdef _ROOT_CINT_

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#include "RVersion.h"
#if ROOT_VERSION_CODE >= 327680 // version newer than 5
#pragma link C++ class bx_message;
#endif

#pragma link C++ class messenger;
#endif
#endif
