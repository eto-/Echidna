/* BOREXINO Reconstruction program
 *
 * Author: Marcin Misiaszek <marcin.misiaszek@lngs.infn.it>
 * Maintainer: Marcin Misiaszek <marcin.misiaszek@lngs.infn.it>
 *
 *
 * Implementation of an 2d interpolation algorithm 
 * Ported from "Numerical Recipes in C++"
 * by William H. Press (Editor), Saul A. Teukolsky (Editor), 
 * William T. Vetterling, Brian P. Flannery 
 */


#ifndef _INTERPOLATOR2D_H
#define _INTERPOLATOR2D_H

#include "bx_named.hh"

template <class T> class NRMat; 
template <class T> class NRVec;

class interpolator2d: public bx_named {
  public:
	typedef double DP;
	typedef const NRVec<DP> Vec_I_DP;
	typedef const NRMat<DP> Mat_I_DP;
	typedef NRVec<DP> Vec_DP, Vec_O_DP;
	typedef NRMat<DP> Mat_DP, Mat_O_DP;
	
	interpolator2d(Vec_I_DP &x1a, Vec_I_DP &x2a, Mat_I_DP &ya);
	virtual ~interpolator2d ();
	
	double get_value (double xx1, double xx2);

  private:
        Vec_DP*  x1c;
	Vec_DP*  x2c;
        Mat_DP*  yc;
	Mat_DP*  y2c;
  
	void splie2(Vec_I_DP &x1a, Vec_I_DP &x2a, Mat_I_DP &ya, Mat_O_DP &y2a);
	void splin2(Vec_I_DP &x1a, Vec_I_DP &x2a, Mat_I_DP &ya, Mat_I_DP &y2a,const DP x1, const DP x2, DP &y);
	void splint(Vec_I_DP &xa, Vec_I_DP &ya, Vec_I_DP &y2a, const DP x, DP &y);
	void spline(Vec_I_DP &x, Vec_I_DP &y, const DP yp1, const DP ypn, Vec_O_DP &y2);
};


template <class T>
class NRVec {
private:
        int nn; // size of array. upper index is nn-1
        T *v;
public:
        NRVec();
        NRVec(int n) : nn(n), v(new T[n]) {};   // Zero-based array
        NRVec(const T &a, int n);       //initialize to constant value
        NRVec(const T *a, int n);       // Initialize to array
        NRVec(const NRVec &rhs);        // Copy constructor
        NRVec & operator=(const NRVec &rhs);    //assignment
        NRVec & operator=(const T &a);  //assign a to every element
        T & operator[](const int i) {return v[i]; }     //i'th element
        const T & operator[](const int i) const;
        int size() const;
        ~NRVec() { if (v != 0) delete[] (v); }
};


template <class T>
class NRMat {
private:
        int nn;
        int mm;
        T **v;
public:
        NRMat();
        NRMat(int n, int m) : nn(n), mm(m), v(new T*[n]) {
							v[0] = new T[m*n];
							for (int i=1; i< n; i++)
							v[i] = v[i-1] + m;
							};                    // Zero-based array
        
	NRMat(const T &a, int n, int m);        //Initialize to constant
        NRMat(const T *a, int n, int m);        // Initialize to array
        NRMat(const NRMat &rhs);                // Copy constructor
        NRMat & operator=(const NRMat &rhs);    //assignment
        NRMat & operator=(const T &a);          //assign a to every element
        T* operator[](const int i) {return v[i];}      //subscripting: pointer to row i
        const T* operator[](const int i) const;
        int nrows() const;
        int ncols() const;
        ~NRMat() { if (v != 0) { delete[] (v[0]); delete[] (v); }; }
};




#endif



