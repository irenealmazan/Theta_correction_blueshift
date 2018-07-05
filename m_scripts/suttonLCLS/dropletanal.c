/* orginal from Larry Lurio (July 2005). */
/* modified  Mark Sutton (June 2006) 
 * rewrote peakfind to handle arbitrarily complicated droplets
 * changed img_anal to handle new ids
 * Change names to dropletfind and dropletanal
 *
 * (Nov 2010) split into separate pieces for matlab and modified for mex 
 * see expandimg, dropletfind, dropletanal and photonize
 */
/* Matlab usage: [npix,adus,xcen,ycen,idlist]=dropletanal(a,dimg,np);
* a is image, dimg is returned from dropletfind (it is modified).
* np is number of droplets.
*/
#include <stdio.h>
#include <stdlib.h>
#include "mex.h"

extern void dropletanal(int *img, int *bimg, int *npix, double *xcen,double *ycen, int *adus, int *idlist, int npeak,int ncol, int nrow);

/* Input Arguments */
#define IMG_IN     prhs[0]
#define DROPIMG_IN prhs[1]
#define NOPEAK     prhs[2]

/* Output Arguments */
#define NPIX   plhs[0]
#define ADUS   plhs[1]
#define XCEN   plhs[2]
#define YCEN   plhs[3]
#define IDLIST plhs[4]

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray*prhs[] )
{ 
    int *img_in,*dropimg_in,npeak,*nopeak; 
    int dims[2];
    int *npix, *adus, *idlist;
    double *xcen, *ycen;
    
    /* Check arguments */
    if (nrhs != 3) { 
	mexErrMsgTxt("input arguments are img, drop_img, nopeaks."); 
    } else if (nlhs != 5) {
	mexErrMsgTxt("Not the correct number of output arguments."); 
        } 
    if(mxINT32_CLASS!=mxGetClassID(IMG_IN)) mexErrMsgTxt("img_in wrong type should be int32()");
    if(mxINT32_CLASS!=mxGetClassID(NOPEAK)) mexErrMsgTxt("no_peak wrong type should be int32()");
    
    nopeak = ((int *) mxGetPr(NOPEAK));
    npeak=nopeak[0];
    dims[0] = mxGetM(IMG_IN); 
    dims[1] = mxGetN(IMG_IN);
    
    /* Create a matrix for the return argument */ 
    NPIX = mxCreateNumericMatrix(1,npeak, mxINT32_CLASS,mxREAL);
    ADUS = mxCreateNumericMatrix(1,npeak, mxINT32_CLASS,mxREAL);
    XCEN = mxCreateDoubleMatrix(1, npeak, mxREAL);
    YCEN = mxCreateDoubleMatrix(1, npeak, mxREAL);
    IDLIST=mxCreateNumericMatrix(1,npeak, mxINT32_CLASS,mxREAL);

    /* Assign pointers to the various parameters */ 
    img_in = mxGetData(IMG_IN);
    dropimg_in = mxGetData(DROPIMG_IN);
    npix = mxGetData(NPIX);
    adus = mxGetData(ADUS);
    xcen = mxGetPr(XCEN);
    ycen = mxGetPr(YCEN);
    idlist = mxGetData(IDLIST);

    /* Do the actual computations in a subroutine */
    dropletanal(img_in, dropimg_in, npix, xcen, ycen,  adus, idlist, npeak, dims[0], dims[1]);
    return;
}
/*for Numerical recipes, modified for ints
updated from my version by online version. June/06 */
void hunt(int xx[], unsigned int n, int x, unsigned int *jlo)
{
	unsigned int jm,jhi,inc;
	int ascnd;

	ascnd=(xx[n] >= xx[1]);
	if (*jlo <= 0 || *jlo > n) {
		*jlo=0;
		jhi=n+1;
	} else {
		inc=1;
		if (x >= xx[*jlo] == ascnd) {
			if (*jlo == n) return;
			jhi=(*jlo)+1;
			while (x >= xx[jhi] == ascnd) {
				*jlo=jhi;
				inc += inc;
				jhi=(*jlo)+inc;
				if (jhi > n) {
					jhi=n+1;
					break;
				}
			}
		} else {
			if (*jlo == 1) {
				*jlo=0;
				return;
			}
			jhi=(*jlo)--;
			while (x < xx[*jlo] == ascnd) {
				jhi=(*jlo);
				inc <<= 1;
				if (inc >= jhi) {
					*jlo=0;
					break;
				}
				else *jlo=jhi-inc;
			}
		}
	}
	while (jhi-(*jlo) != 1) {
		jm=(jhi+(*jlo)) >> 1;
		if (x >= xx[jm] == ascnd)
			*jlo=jm;
		else
			jhi=jm;
	}
  if (x == xx[n]) *jlo=n;/* special for droplet*/
  if (x == xx[1]) *jlo=1;
}
/* (C) Copr. 1986-92 Numerical Recipes Software ^. */
/*Should adus vector be float or double?*/
void dropletanal(int *img, int *bimg, int *npix, double *xcen,double *ycen, int *adus,int *idlist, int npeak,int ncol, int nrow) {
/* Loop through image and analyze droplets, for number of pixels, center of mass
 * and total intensity. Along the way relabel pixel ids so each pixel in a
 * droplet points to the next pixel, with the last pixel pointing to the first.
 * This makes looping through droplets easy. NOTE this modifies bimg.
 */
	int i,pnum,id,pos;
	int ct;
	double x,y;
	/*printf("dropletana\n");*/
	pnum=0; /*size of idlist*/
	ct=0;
	for (i=0;i<nrow*ncol;i++) {
		if (bimg[i]) { 
			id=bimg[i];
			if(id!=i-1){ /* this makes cycles so photonize can traverse droplets*/
				bimg[i]=bimg[id-1]; /*swap ptrs*/
				bimg[id-1]=i+1;
				}
			if(pnum==0 || id>idlist[pnum-1]) { /*must be if a new droplet.*/
			    idlist[pnum]=id;
			    pnum++;
			    pos=pnum;
			    /*printf("add %d,%d,%d,%d\n",id,pos,idlist[pos-1],idlist[pos]);*/
			    }	
			if( id!= idlist[pos-1]) {
				hunt(idlist-1,pnum,id,&pos);
			        /*printf("hunt %d,%d,%d,%d\n",id,pos,idlist[pos-1],idlist[pos]);*/
				if(idlist[pos-1]!=id) {
				    /*I don't think this should happen and I*/
				    /*will remove this at some point if it*/
				    /*doesn't happen in practice.*/
				    /*If it does happen, figure out why and fix.*/
				    mexPrintf("error should not be here, not in list\n");
			    mexPrintf("err %d,%d,%d,%d,%d\n",id,pos,idlist[pos-1],idlist[pos],idlist[pnum-1]);
				    idlist[pnum]=id; /* this is proabably wrong*/
				    pnum++;
				    ct++;
				    if(ct>10) return;
				    pos=pnum;
				    }	
				}
			x=(double)(i%ncol);
			y=(double)(i/ncol);
			npix[pos-1]+=1;
			xcen[pos-1]+=x*(double)(img[i]);
			ycen[pos-1]+=y*(double)(img[i]);
			adus[pos-1]+=img[i];
			}
		}
	for (i=0;i<npeak;i++) {
		if (adus[i]) { /* offset 1 for matlab indices */
			xcen[i] = xcen[i]/(double)(adus[i])+1.0;
			ycen[i] = ycen[i]/(double)(adus[i])+1.0;
			}
		else {
			mexPrintf("bad droplet for i=%d\n",i);
			mexPrintf("xcen=%d ycen=%d adus=%d npix=%d\n",xcen[i],ycen[i],adus[i],npix[i]);
		        ct++;
		        if(ct>10) return;
			}
	    }
}
