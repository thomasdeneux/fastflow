#include <math.h>

#include <Balloon\debuggage.h>
#include <MEX_files\mex_tools.h> 

void mexFunction( int nargout, mxArray* pargout[], 
		  int nargin, const mxArray* pargin[] )    
{

	//double* p = mxGetPr(pargin[0]);
	//matlab << p << " -> " << p[0] << "\n";
	//p[0]++;
	//matlab << p << " -> " << p[0] << "\n";


	
	// Input
	checkargs(nargout,nargin,0,1,5,5);
	const mxArray* Y = pargin[0];
	const int* ydims = mxGetDimensions(Y);
	int nj=ydims[0], ni=ydims[1], nt=ydims[2];
	double* Yp = mxGetPr(Y);
	const mxArray* H = pargin[1];
	double* Hp = mxGetPr(H);
	int j = mxArray2int(pargin[2]), i = mxArray2int(pargin[3]), filterhw = mxArray2int(pargin[4]);
	int filterwidth = 2*filterhw+1; // c'est aussi la taille de la matrice carrée H

	// Code -> equivalent to Matlab:
	// cross = repmat(H(:,:,theta),[1 1 nt]) .* Y(j-filterhw:j+filterhw,i-filterhw:i+filterhw,:);
	// cross = reshape(cross,filterwidth^2,nt);
	// Y(j,i-filterhw-1,:) = shiftdim(sum(cross),-1);
	int yidx = (j-1)+nj*(i-1);
	int yidxdec = (j-1)+nj*(i-filterhw-2);
	int nij = nj*ni;
	for (int t=0; t<nt; t++, yidx+=nij, yidxdec+=nij) {
		double sum = 0;
		int hidx = 0;
		int yidx2 = yidx - filterhw*nj;
		for (int i=-filterhw; i<=filterhw; i++, yidx2+=nj) {
			int yidx3 = yidx2 - filterhw;
			for (int j=-filterhw; j<=filterhw; j++, hidx++, yidx3++)
				sum += Hp[hidx] * Yp[yidx3];
		}
		Yp[yidxdec] = sum;
	}

	// Voilà c'est fini (Y changé par ce passage par référence normalement interdit sous Matlab !)

}