/* orginal from Larry Lurio (July 2005). */
/* modified  Mark Sutton (June 2006) 
 * rewrote peakfind to handle arbitrarily complicated droplets
 * changed img_anal to handle new ids
 * Change names to dropletfind and dropletanal
 *
 * (Nov 2010) split into separate pieces for matlab and modified for mex 
 * see expandimg, dropletfind, dropletanal and photonize
 */
/* Matlab usage: pimg=photonize(a,dimg,npix,photsnum,idlist,np);
 * Use results from dropletanal to place photons into image pimg.
 * photsnum is number of photons in a droplet (done in matlab from adus).
 * NOTE: a is modified by this code, so keep a copy.
 */
#include <stdio.h>
#include <stdlib.h>
#include "mex.h"

void photonize(int *img_out, int *img_in, int *bimg, int *npix, int *photsnum, int *idlist, int npeak,int ncol, int nrow, int adupphot);

/* Input Arguments */
#define IMG_IN     prhs[0]
#define DROPIMG_IN prhs[1]
#define NPIX       prhs[2]
#define PHOTSNUM   prhs[3]
#define IDLIST     prhs[4]
#define NOPEAK     prhs[5]
#define ADUPPHOT   prhs[6]

/* Output Arguments */
#define IMG_OUT   plhs[0]

void mexFunction( int nlhs, mxArray *plhs[], 
		  int nrhs, const mxArray*prhs[] )
     
{ 
    int *img_in,*dropimg_in,npeak,*nopeak; 
    int *img_out;
    int *npix, *photsnum, *idlist;
    int adupphot,*app;
    int dims[2];
    
    /* Check for proper number of arguments */
    if (nrhs != 7) { 
	mexErrMsgTxt("input arguments are img_in,drop_img,npix,photsnum,idlist,npeak"); 
    } else if (nlhs != 1) {
	mexErrMsgTxt("Not the correct number of output arguments."); 
        } 
    if(mxINT32_CLASS!=mxGetClassID(IMG_IN)) mexErrMsgTxt("img_in wrong type should be int32()");
    if(mxINT32_CLASS!=mxGetClassID(DROPIMG_IN)) mexErrMsgTxt("dropimg_in wrong type should be int32()");
    if(mxINT32_CLASS!=mxGetClassID(NPIX)) mexErrMsgTxt("npix wrong type should be int32()");
    if(mxINT32_CLASS!=mxGetClassID(PHOTSNUM)) mexErrMsgTxt("photsnum wrong type should be int32()");
    if(mxINT32_CLASS!=mxGetClassID(NOPEAK)) mexErrMsgTxt("nopeak wrong type should be int32()");
    if(mxINT32_CLASS!=mxGetClassID(ADUPPHOT)) mexErrMsgTxt("adupphot wrong type should be int32()");
    
    dims[0] = mxGetM(IMG_IN); 
    dims[1] = mxGetN(IMG_IN);
    /*printf("dims=(%d,%d)\n",dims[0],dims[1]);*/
    
    /* Create a matrix for the return argument */ 
    IMG_OUT= mxCreateNumericMatrix((mwSize)dims[0], (mwSize)dims[1], mxINT32_CLASS,mxREAL);
    img_out = mxGetData(IMG_OUT);

    /* Assign pointers to the various parameters */ 
    img_in = mxGetData(IMG_IN);
    dropimg_in = mxGetData(DROPIMG_IN);
    npix = mxGetData(NPIX);
    photsnum = mxGetData(PHOTSNUM);
    idlist = mxGetData(IDLIST);

    nopeak = ((int *) mxGetPr(NOPEAK));
    npeak=nopeak[0];
    app = ((int *) mxGetPr(ADUPPHOT));
    adupphot=app[0];

    /* Do the actual computations in a subroutine */
    photonize(img_out,img_in,dropimg_in, npix, photsnum, idlist, npeak, dims[0], dims[1], adupphot);
    return;
    
}

void photonize(int *img_out, int *img_in, int *bimg, int *npix, int *photsnum,int *idlist, int npeak,int ncol, int nrow, int adupphot) {
/* Convert droplets to photons. Uses simple rule to place first photon at
 * largest pixel in droplet then reducing this pixel. It then repeats this
 * for the number of photons in droplet.
 */

	int i,j,id,nophots,pos,pmax,max;
	/*printf("photonize %d %d\n",minadu,adupphot);*/
	/* for each droplet */
	for(i=0;i<npeak;i++) {
	    id=idlist[i]-1;
	    /*printf("id=%d\n",id+1);*/
	    nophots = photsnum[i];
	    /* if(nophots==0) printf("zero %d\n",i);*/
	    for(j=0;j<nophots;j++) {
	       /*find max in droplet and delete for each photon*/
	       pos=id;
	       max=-100000000; /* negative infinity*/
	       pmax=pos;
	       do {
		  /*printf("pos=%d\n",pos+1);*/
	  /*Are we building in correlations by using first max?*/
	  /*if(img_in[pos]==max) printf("double max %d,np=%d\n",pos+1,nophots);*/
	          if(img_in[pos]>max) {
	             pmax=pos;
	             max=img_in[pos];
	             }
		  pos=bimg[pos]-1;
	          }
	       while(pos!=id);
	       /*printf("pmax=%d\n",pmax+1);*/
	       img_out[pmax]++;
	       img_in[pmax] -= adupphot;
	       }
	}
}
