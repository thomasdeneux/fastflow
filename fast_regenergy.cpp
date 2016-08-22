#include <math.h>
#include <mex.h>

class vecteur {
  int n;
  double* t;
public:
  vecteur(int N) {n=N; t=new double[n];}
  vecteur(const vecteur& v) {n=v.n; t=new double[n]; for (int i=0;i<n;i++) t[i]=v.t[i];}
  ~vecteur() {delete [] t;}
  const vecteur& operator=(const vecteur& v) {
    delete [] t; n=v.n; t=new double[n]; 
    for (int i=0;i<n;i++) t[i]=v.t[i];
    return (*this);
  }

  inline double operator()(int i) const {
    if (i<0 || i>=n) mexErrMsgTxt("vector access outside range"); 
    return t[i];
  } 
  inline double& operator()(int i) {
    if (i<0 || i>=n) mexErrMsgTxt("vector writing outside range"); 
    return t[i];
  } 

  inline int size() const{return n;}

  vecteur operator+ (const vecteur& v) const {
    if (n!=v.n) mexErrMsgTxt("inconsistent sizes in vector + vector");
    vecteur w(n);
    for (int i=0; i<n; i++) w(i)=t[i]+v(i);
    return w;
  }
  vecteur operator- (const vecteur& v) const {
    if (n!=v.n) mexErrMsgTxt("inconsistent sizes in vector - vector");
    vecteur w(n);
    for (int i=0; i<n; i++) w(i)=t[i]-v(i);
    return w;
  }
  
  friend class matrice;
};
  
class matrice {
  int m,n;
  double* t;
  bool mxfree;
public:
  matrice(int M, int N) {m=M; n=N; mxfree=true; t=new double[m*n];}
  matrice(const matrice& v) {
    m=v.m; n=v.n; mxfree=true; t=new double[m*n]; for (int i=0;i<m*n;i++) t[i]=v.t[i];}
  matrice(const mxArray* A) {
    if (!mxIsDouble(A) || mxIsComplex(A))
      mexErrMsgTxt("mxArray2matrice used on non-real array");
    m=mxGetM(A); n=mxGetN(A);
    t=mxGetPr(A);
    mxfree=false;
  }
  ~matrice() {if (mxfree) delete [] t;}
  const matrice& operator=(const matrice& v) {
    delete [] t; m=v.m; n=v.n; t=new double[m*n]; 
    for (int i=0;i<m*n;i++) t[i]=v.t[i]; 
    return (*this);
  }

  inline double operator()(int i, int j) const {
    if (i<0 || i>=m || j<0 || j>=n) mexErrMsgTxt("matrix access outside range"); 
    return t[i+j*m];
  } 
  inline double& operator()(int i, int j) {
    if (i<0 || i>=m || j<0 || j>=n) mexErrMsgTxt("matrix writing outside range"); 
    return t[i+j*m];
  } 
  inline double operator[](int i) const {
    if (i<0 || i>=m*n) mexErrMsgTxt("matrix vector-access outside range"); 
    return t[i];
  } 
  inline double& operator[](int i) {
    if (i<0 || i>=m*n) mexErrMsgTxt("matrix vector-writing outside range"); 
    return t[i];
  } 

  inline int nlin() const {return m;}
  inline int ncol() const {return n;}
  inline double* pointer() const {return t;}
  inline double*& pointer() {return t;}

  vecteur operator*(const vecteur& v) const {
    if (n!=v.n) mexErrMsgTxt("inconsistent sizes in matrix x vector");
    vecteur w(m);
    for (int i=0; i<m; i++) {
      w(i)=0;
      for (int j=0; j<n; j++) w(i)+=(*this)(i,j)*v(j);
    }
    return w; 
  }
  matrice operator*(const matrice& a) const {
    if (n!=a.m) mexErrMsgTxt("inconsistent sizes in matrix x matrix");
    matrice b(m,a.n);
    for (int i=0; i<m; i++) {
      for (int j=0; j<a.n; j++) {
	b(i,j)=0;
	for (int k=0; k<n; k++) b(i,j)+=(*this)(i,k)*a(k,j);
      }
    }
    return b; 
  }
};
  
void checkargs(int nargout, int nargin, int out, int in){
	if (nargin!=in)
		mexErrMsgTxt("Wrong number of input arguments");
	if (nargout!=out) 
		mexErrMsgTxt("Wrong number of output arguments");
}
void checkargs(int nargout, int nargin, int minout, int maxout, int minin, int maxin){
	if (nargin<minin || nargin>maxin)
		mexErrMsgTxt("Wrong number of input arguments");
	if (nargout<minout || nargout>maxout)
		mexErrMsgTxt("Wrong number of output arguments");
}

vecteur mxArray2vecteur(const mxArray* A){
    if (!mxIsDouble(A) || mxIsComplex(A))
		mexErrMsgTxt("mxArray2vecteur used on non-real array");
	if (mxGetM(A)==0 || mxGetN(A)==0) return vecteur(0);
	if (mxGetM(A)!=1 && mxGetN(A)!=1)
		mexErrMsgTxt("mxArray2vecteur used on non-vecteur array");
	int n=mxGetM(A)*mxGetN(A);
	double* p=mxGetPr(A);
	vecteur a(n); 
	for (int i=0; i<n; i++) a(i)=p[i];
	return a;
}

matrice mxArray2matrice(const mxArray* A){
    if (!mxIsDouble(A) || mxIsComplex(A))
		mexErrMsgTxt("mxArray2matrice used on non-real array");
	if (mxGetM(A)==0 || mxGetN(A)==0) return matrice(0,0);
	int m=mxGetM(A), n=mxGetN(A);
	double* p=mxGetPr(A);
	matrice a(m,n);
	for (int i=0; i<m*n; i++) a[i]=p[i];
	return a;
}

mxArray* vecteur2mxArray(const vecteur& x){
	mxArray* B = mxCreateDoubleMatrix(x.size(),(x.size())?1:0,mxREAL);
	double* pr = mxGetPr(B);
	for (unsigned int i=0; i<x.size(); i++) pr[i] = x(i);
	return B;
}

mxArray* matrice2mxArray(const matrice& x){
	mxArray* B = mxCreateDoubleMatrix(x.nlin(),x.ncol(),mxREAL);
	double* pr = mxGetPr(B);
	for (unsigned int i=0; i<x.nlin(); i++) 
	  for (unsigned int j=0; j<x.ncol(); j++) 
	    pr[i+j*x.nlin()] = x(i,j);
	return B;
}

/////////////////////////////////////////////////////////////////////
// END OF TOOLS
///////////////////////////////////////////////////////////////////

matrice makematrice(vecteur iso, int nj, int ni) {
	double theta=iso(0);
	vecteur shift(2); shift(0)=iso(1); shift(1)=iso(2);
	matrice R(2,2); R(0,0)=cos(theta); R(0,1)=-sin(theta); R(1,0)=sin(theta); R(1,1)=cos(theta);
	vecteur C(2); C(0)=ni/2; C(1)=nj/2;
	vecteur T = shift + C - R*C;
	matrice M(2,3);
	M(0,0)=R(0,0); M(0,1)=R(0,1); M(1,0)=R(1,0); M(1,1)=R(1,1); 
	M(0,2)=T(0); M(1,2)=T(1);
	return M;
}

void floormodulo(double u, int& p, double& x) {
	p = (int) u;
	x = u-p;
}

void mexFunction( int nargout, mxArray* pargout[], 
		  int nargin, const mxArray* pargin[] )    
{

	// Cas fast_regenergy(iso,CS2) -> calcul transform�e inverse 
	// (flemme d'�crire une autre fonction qui devrait plut�t s'appeler 'fast_resample')
	if (nargin==2) {
		vecteur iso = mxArray2vecteur(pargin[0]);
		matrice CS2(pargin[1]); int nj=CS2.nlin(), ni=CS2.ncol();
		matrice CS2inv(nj,ni); 
		for (int i=0; i<ni*nj; i++) CS2inv[i]=0;
		matrice M = makematrice(iso,nj,ni); 
		for (int j=0; j<nj; j++) {
			for (int i=0; i<ni; i++) {
				vecteur w(3); w(0)=i+1; w(1)=j+1; w(2)=1;
				w = M*w;
				double u=w(0)-1, v=w(1)-1;
				int p,q; double x,y;
				floormodulo(u,p,x); floormodulo(v,q,y);
				// linear interpolation
				if (p<0 || p>=ni-1 || q<0 || q>=nj-1) CS2inv(j,i)=0;
				else CS2inv(j,i) = (1-x)*((1-y)*CS2(q,p)+y*CS2(q+1,p))+x*((1-y)*CS2(q,p+1)+y*CS2(q+1,p+1));
			}
		}
		pargout[0] = matrice2mxArray(CS2inv);
		return;
	}
 
	// Input
	checkargs(nargout,nargin,0,1,4,4);
	vecteur iso = mxArray2vecteur(pargin[0]);
	if (iso.size()==2) { // translation only -> set rotation = 0
		vecteur tmp = iso;
		iso = vecteur(3);
		iso(0)=0; iso(1)=tmp(0); iso(2)=tmp(1);
	}
	if (iso.size()!=3) mexErrMsgTxt("first argument (iso) should be length 2 or 3 vector");
	matrice CS1(pargin[1]); int nj=CS1.nlin(), ni=CS1.ncol();
	matrice CS2(pargin[2]); 
	if (CS2.nlin()!=nj || CS2.ncol()!=ni) mexErrMsgTxt("second and third arguments should be same size images");
	vecteur indices = mxArray2vecteur(pargin[3]);

	matrice M = makematrice(iso,nj,ni); 

	int nind=indices.size();
	matrice coords(3,nind);
	for (int k=0; k<nind; k++) {
		int ind=((int) indices(k))-1;
		int j=ind%nj, i=ind/nj;
		coords(0,k) = i+1;
		coords(1,k) = j+1;
		coords(2,k) = 1;
	} 
	matrice coordsM = M * coords;

	//double E = 0; // energy
	vecteur E(nind);
	for (int k=0; k<nind; k++) {
	  int j=((int) coords(1,k))-1, i=((int) coords(0,k))-1;
	  double u=coordsM(0,k)-1, v=coordsM(1,k)-1;
	  int p,q; double x,y;
	  floormodulo(u,p,x); floormodulo(v,q,y);
	  // linear interpolation
	  double CS2invji;
	  if (p<0 || p>=ni-1 || q<0 || q>=nj-1) CS2invji=0;
	  else CS2invji = (1-x)*((1-y)*CS2(q,p)+y*CS2(q+1,p))+x*((1-y)*CS2(q,p+1)+y*CS2(q+1,p+1));
	  double diff = (CS1(j,i)-CS2invji);
	  //E += diff*diff;
	  E(k) = diff/100;
	}
	//E = sqrt(E/nind);
 
	// Output
	//pargout[0] = double2mxArray(E);
	pargout[0] = vecteur2mxArray(E);

	// return
	return;
}   
