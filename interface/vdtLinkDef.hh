/* BOREXINO Reconstruction program
 *
 * vdt linkdef
 */
#ifndef _VDT_LINKDEF_HH
#define _VDT_LINKDEF_HH

#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ class vdt;
#pragma link C++ class std::vector<vdt>;
#pragma link C++ function operator<< (std::ostream&, const vdt&);
#endif
#endif
