/* BOREXINO Reconstruction program
 *
 * vdt linkdef
 */
#ifndef _VDT_LINKDEF_HH
#define _VDT_LINKDEF_HH
#if defined(__ROOTCLING__) || defined(__CINT__)
#define _ROOT_CINT_
#endif

#ifdef _ROOT_CINT_

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ class bx_vdt;
#pragma link C++ class std::vector<bx_vdt>;
#pragma link C++ function operator<< (std::ostream&, const bx_vdt&);
#endif
#endif
