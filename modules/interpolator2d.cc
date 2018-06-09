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


#include "interpolator2d.hh"

template class NRMat<double>;
template class NRVec<double>;

interpolator2d::interpolator2d(Vec_I_DP &x1a, Vec_I_DP &x2a, Mat_I_DP &ya): bx_named("interpolator2d")
	    {
	    y2c = new Mat_DP(x1a.size(), x2a.size());
	    yc  = new Mat_DP(x1a.size(), x2a.size());
	    x1c = new Vec_DP(x1a.size());
	    x2c = new Vec_DP(x2a.size());
	    
	    *x1c = x1a;
	    *x2c = x2a;
	    *yc  = ya;
	    splie2( *x1c, *x2c, *yc, *y2c);
	    };
	    
double interpolator2d::get_value (double xx1, double xx2)
	    {
	    double f;
	    splin2(*x1c,*x2c,*yc,*y2c,xx1,xx2,f);
	    return f;
	    }

interpolator2d::~interpolator2d () { delete x1c; delete x2c; delete yc; delete y2c; }

	    

template <class T>
NRVec<T>::NRVec() : nn(0), v(0) {}

/*
template <class T>
NRVec<T>::NRVec(int n) : nn(n), v(new T[n]) {}
*/

template <class T>
NRVec<T>::NRVec(const T& a, int n) : nn(n), v(new T[n])
{
	for(int i=0; i<n; i++)
		v[i] = a;
}


template <class T>
NRVec<T>::NRVec(const T *a, int n) : nn(n), v(new T[n])
{
        for(int i=0; i<n; i++)
                v[i] = *a++;
}

template <class T>
NRVec<T>::NRVec(const NRVec<T> &rhs) : nn(rhs.nn), v(new T[nn])
{
        for(int i=0; i<nn; i++)
                v[i] = rhs[i];
}

template <class T>
NRVec<T> & NRVec<T>::operator=(const NRVec<T> &rhs)
// postcondition: normal assignment via copying has been performed;
//              if vector and rhs were different sizes, vector
//              has been resized to match the size of rhs
{
        if (this != &rhs)
        {
                if (nn != rhs.nn) {
                        if (v != 0) delete [] (v);
                        nn=rhs.nn;
                        v= new T[nn];
                }
                for (int i=0; i<nn; i++)
                        v[i]=rhs[i];
        }
        return *this;
}

template <class T>
NRVec<T> & NRVec<T>::operator=(const T &a)      //assign a to every element
{
        for (int i=0; i<nn; i++)
                v[i]=a;
        return *this;
}

/*
template <class T>
inline T & NRVec<T>::operator[](const int i)    //subscripting
{
        return v[i];
}
*/

template <class T>
inline const T & NRVec<T>::operator[](const int i) const        //subscripting
{
        return v[i];
}

template <class T>
inline int NRVec<T>::size() const
{
        return nn;
}

/*
template <class T>
NRVec<T>::~NRVec()
{
        if (v != 0)
                delete[] (v);
}
*/

template <class T>
NRMat<T>::NRMat() : nn(0), mm(0), v(0) {}


/*
template <class T>
NRMat<T>::NRMat(int n, int m) : nn(n), mm(m), v(new T*[n])
{
	v[0] = new T[m*n];
	for (int i=1; i< n; i++)
		v[i] = v[i-1] + m;
}
*/


template <class T>
NRMat<T>::NRMat(const T &a, int n, int m) : nn(n), mm(m), v(new T*[n])
{
        int i,j;
        v[0] = new T[m*n];
        for (i=1; i< n; i++)
                v[i] = v[i-1] + m;
        for (i=0; i< n; i++)
                for (j=0; j<m; j++)
                        v[i][j] = a;
}

template <class T>
NRMat<T>::NRMat(const T *a, int n, int m) : nn(n), mm(m), v(new T*[n])
{
        int i,j;
        v[0] = new T[m*n];
        for (i=1; i< n; i++)
                v[i] = v[i-1] + m;
        for (i=0; i< n; i++)
                for (j=0; j<m; j++)
                        v[i][j] = *a++;
}

template <class T>
NRMat<T>::NRMat(const NRMat &rhs) : nn(rhs.nn), mm(rhs.mm), v(new T*[nn])
{
        int i,j;
        v[0] = new T[mm*nn];
        for (i=1; i< nn; i++)
                v[i] = v[i-1] + mm;
        for (i=0; i< nn; i++)
                for (j=0; j<mm; j++)
                        v[i][j] = rhs[i][j];
}

template <class T>
NRMat<T> & NRMat<T>::operator=(const NRMat<T> &rhs)
// postcondition: normal assignment via copying has been performed;
//              if matrix and rhs were different sizes, matrix
//              has been resized to match the size of rhs
{
        if (this != &rhs) {
                int i,j;
                if (nn != rhs.nn || mm != rhs.mm) {
                        if (v != 0) {
                                delete[] (v[0]);
                                delete[] (v);
                        }
                        nn=rhs.nn;
                        mm=rhs.mm;
                        v = new T*[nn];
                        v[0] = new T[mm*nn];
                }
                for (i=1; i< nn; i++)
                        v[i] = v[i-1] + mm;
                for (i=0; i< nn; i++)
                        for (j=0; j<mm; j++)
                                v[i][j] = rhs[i][j];
        }
        return *this;
}

template <class T>
NRMat<T> & NRMat<T>::operator=(const T &a)      //assign a to every element
{
        for (int i=0; i< nn; i++)
                for (int j=0; j<mm; j++)
                        v[i][j] = a;
        return *this;
}

/*
template <class T>
inline T* NRMat<T>::operator[](const int i)     //subscripting: pointer to row i
{
        return v[i];
}
*/

template <class T>
inline const T* NRMat<T>::operator[](const int i) const
{
        return v[i];
}

template <class T>
inline int NRMat<T>::nrows() const
{
        return nn;
}

template <class T>
inline int NRMat<T>::ncols() const
{
        return mm;
}

/*
template <class T>
NRMat<T>::~NRMat()
{
        if (v != 0) {
                delete[] (v[0]);
                delete[] (v);
        }
}
*/

void interpolator2d::splie2(Vec_I_DP &x1a, Vec_I_DP &x2a, Mat_I_DP &ya, Mat_O_DP &y2a)
{
	int m,n,j,k;

	m=x1a.size();
	n=x2a.size();
	Vec_DP ya_t(n),y2a_t(n);
	for (j=0;j<m;j++) {
		for (k=0;k<n;k++) ya_t[k]=ya[j][k];
		spline(x2a,ya_t,1.0e30,1.0e30,y2a_t);
		for (k=0;k<n;k++) y2a[j][k]=y2a_t[k];
	}
}

void interpolator2d::spline(Vec_I_DP &x, Vec_I_DP &y, const DP yp1, const DP ypn,
	Vec_O_DP &y2)
{
	int i,k;
	DP p,qn,sig,un;

	int n=y2.size();
	Vec_DP u(n-1);
	if (yp1 > 0.99e30)
		y2[0]=u[0]=0.0;
	else {
		y2[0] = -0.5;
		u[0]=(3.0/(x[1]-x[0]))*((y[1]-y[0])/(x[1]-x[0])-yp1);
	}
	for (i=1;i<n-1;i++) {
		sig=(x[i]-x[i-1])/(x[i+1]-x[i-1]);
		p=sig*y2[i-1]+2.0;
		y2[i]=(sig-1.0)/p;
		u[i]=(y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]);
		u[i]=(6.0*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p;
	}
	if (ypn > 0.99e30)
		qn=un=0.0;
	else {
		qn=0.5;
		un=(3.0/(x[n-1]-x[n-2]))*(ypn-(y[n-1]-y[n-2])/(x[n-1]-x[n-2]));
	}
	y2[n-1]=(un-qn*u[n-2])/(qn*y2[n-2]+1.0);
	for (k=n-2;k>=0;k--)
		y2[k]=y2[k]*y2[k+1]+u[k];
}


void interpolator2d::splin2(Vec_I_DP &x1a, Vec_I_DP &x2a, Mat_I_DP &ya, Mat_I_DP &y2a,
	const DP x1, const DP x2, DP &y)
{
	int j,k;

	int m=x1a.size();
	int n=x2a.size();
	Vec_DP ya_t(n),y2a_t(n),yytmp(m),ytmp(m);
	for (j=0;j<m;j++) {
		for (k=0;k<n;k++) {
			ya_t[k]=ya[j][k];
			y2a_t[k]=y2a[j][k];
		}
		splint(x2a,ya_t,y2a_t,x2,yytmp[j]);
	}
	spline(x1a,yytmp,1.0e30,1.0e30,ytmp);
	splint(x1a,yytmp,ytmp,x1,y);
}


void interpolator2d::splint(Vec_I_DP &xa, Vec_I_DP &ya, Vec_I_DP &y2a, const DP x, DP &y)
{
	int k;
	DP h,b,a;

	int n=xa.size();
	int klo=0;
	int khi=n-1;
	while (khi-klo > 1) {
		k=(khi+klo) >> 1;
		if (xa[k] > x) khi=k;
		else klo=k;
	}
	h=xa[khi]-xa[klo];
//	if (h == 0.0) nrerror("Bad xa input to routine splint"); -->  message
	a=(xa[khi]-x)/h;
	b=(x-xa[klo])/h;
	y=a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2a[klo]
		+(b*b*b-b)*y2a[khi])*(h*h)/6.0;
}
