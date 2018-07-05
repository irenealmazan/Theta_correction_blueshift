/* orginal from Larry Lurio (July 2005). */
/* substantially modified  Mark Sutton (June 2006) 
 * rewrote peakfind to handle arbitrarily complicated droplets
 * changed img_anal to handle new ids
 * Change names to dropletfind and dropletanal
 *
 * (Nov 2010) split into separate pieces for matlab and modified for mex 
 * see expandimg, dropletfind, dropletanal and photonize
 */

/* Matlab usage: img_out = expandimg(img_in);
 * expand bitmap image so all pixels neighbouring a on bit are turned on
 */
#include <stdio.h>
#include <stdlib.h>
#include "mex.h"

extern void expandimg(int *img_out, int *img_in, int ncol, int nrow);

/* Input Arguments */
#define IMG_IN  prhs[0]

/* Output Arguments */
#define IMG_OUT  plhs[0]

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray*prhs[] )
{ 
    int *img_out,*img_in; 
    mwSize m,n; 
    int dims[2];
    
    /* Check for proper number of arguments */
    if (nrhs != 1) { 
	mexErrMsgTxt("Two input arguments required."); 
    } else if (nlhs > 2) {
	mexErrMsgTxt("Too many output arguments."); 
    	} 
    
    /* Get dimensions of Y. (watch out for implied transpose) */ 
    dims[0] = mxGetM(IMG_IN); 
    dims[1] = mxGetN(IMG_IN);
    if(mxINT32_CLASS!=mxGetClassID(IMG_IN)) mexErrMsgTxt("img_in wrong type should be int32()");
    /* printf("dims=(%d,%d)\n",dims[0],dims[1]); */
    
    /* Create a matrix for the return argument */ 
    IMG_OUT= mxCreateNumericMatrix((mwSize)dims[0], (mwSize)dims[1], mxINT32_CLASS,mxREAL);

    /* Assign pointers to the various parameters */ 
    img_out = mxGetData(IMG_OUT);
    img_in = mxGetData(IMG_IN);

    /* Do the actual computations in a subroutine */
    expandimg(img_out,img_in, dims[0],dims[1]);
    return;
}

void expandimg(int *img_out,int *img_in, int ncol, int nrow)
{
	/*This is not always correct, seems to be problems at edges only.*/
	int i,m,n;

	for (i=0;i<nrow*ncol;i++) {
		if (img_in[i]) {
			m=i%ncol; n=i/ncol;/* index of matrix */
			/* three below */
			if (n>0) {
				if(m>0) img_out[i-ncol-1] = 1;
				img_out[i-ncol] = 1;
				if (m<(ncol-1)) img_out[i-ncol+1] = 1;
			        }
			/* three at level */
			if (m>0) img_out[i-1] = 1;
			img_out[i] = 1;
			if (m<(ncol-1)) img_out[i+1] = 1;
			/* three above */
			if (n < (ncol-1)) {
				if (m>0) img_out[i+ncol-1] = 1;
				img_out[i+ncol] = 1;
				if (m<(ncol-1)) img_out[i+ncol+1] = 1;
			        }
			}
		}
}
