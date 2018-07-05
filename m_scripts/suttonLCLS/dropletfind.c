/* orginal from Larry Lurio (July 2005). */
/* modified  Mark Sutton (June 2006) 
 * rewrote peakfind to handle arbitrarily complicated droplets
 * changed img_anal to handle new ids
 * Change names to dropletfind and dropletanal
 *
 * (Nov 2010) split into separate pieces for matlab and modified for mex 
 * see expandimg, dropletfind, dropletanal and photonize
 */
/* Matlab usage: [dimg,np]=dropletfind(int32(b));
 * where b is bitmap for pixels to be in a droplet (ie b=image>100).
 * use b=expandimg(b) to include boundaries of droplets. dimg is
 * a droplet id map and np is number of droplets found.
*/
#include <stdio.h>
#include <stdlib.h>
#include "mex.h"

extern void dropletfind(int *img_out, int *img_in, int ncol, int nrow, int *npeak);

/* Input Arguments */
#define IMG_IN  prhs[0]

/* Output Arguments */
#define IMG_OUT  plhs[0]
#define NOPEAK   plhs[1]

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray*prhs[] )
{ 
    int *img_out,*img_in,*nopeak; 
    int t,dims[2];
    const mwSize scl1[]={1};
    mxClassID category;
    
    /* Check arguments */
    if (nrhs != 1) { 
	mexErrMsgTxt("Two input arguments required."); 
    } else if (nlhs > 2) {
	mexErrMsgTxt("Too many output arguments."); 
        } 
    if(mxINT32_CLASS!=mxGetClassID(IMG_IN)) mexErrMsgTxt("img_in wrong type should be int32()");
    dims[0] = mxGetM(IMG_IN); 
    dims[1] = mxGetN(IMG_IN);
    /*printf("dims=(%d,%d)\n",dims[0],dims[1]);*/
    
    /* Create a matrix for the return argument */ 
    IMG_OUT= mxCreateNumericMatrix((mwSize)dims[0], (mwSize)dims[1], mxINT32_CLASS,mxREAL);
    NOPEAK=mxCreateNumericArray(1,scl1,mxINT32_CLASS,mxREAL);

    /* Assign pointers to the various parameters */ 
    img_out = mxGetData(IMG_OUT);
    img_in = mxGetData(IMG_IN);
    nopeak = mxGetData(NOPEAK);

    /* Do the actual computations in a subroutine */
    dropletfind(img_out,img_in, dims[0],dims[1], nopeak);
    return;
}



void dropletfind(int *img_out,int *img_in, int ncol, int nrow, int *npeak) {
/* As we sweep through the bitmap, check left and above pixels to see
 * if in a droplet. If so add current point to droplet. If both
 * left and above are droplets may need to merge them.
 * Trick is to use the index of the "first" occurence (ie upper left-most) of
 * a pixel in a droplet as the droplet id (plus 1 because background will be
 * droplet 0).  Then each added pixel to a droplet will get this id. This
 * allows for a straightforward merging of complicated (hairy-like) droplets
 * with not much looping over droplet ids.
 */
	int i,j,curr,above,left,peakcount,complicated;
	int typ,oldtyp,pos,pos1,pos2;
        complicated=0;
	/*Have special cases of 1)first pixel      2)rest of first row
	 *		        3)start of columns 4) interior pixels
	 *First Pixel
	 */
	if(img_in[0]) img_out[0]=1;
	peakcount=img_out[0];
	/*Rest of first row*/
	for(i=1;i<ncol;i++) { /* id first row */
	    if(img_in[i]) {
		if(img_in[i-1]) img_out[i]=img_out[i-1];
		else {
		    img_out[i]=i+1;
		    peakcount++;
		    }
		}
	    }
	/* matrix layout(in C) [[0,1,2,...,ncol-1]
	 *                      [ncol,...,2ncol-1]
	 *                      ...
	 *                      [(nrow-1)ncol,...,nrow*ncol-1]]
	 * so above of i is i-ncol and left of i is i-1
	 */
	curr=ncol; /*current pixel to study*/
	above=0;
	left=ncol-1;
	for (i=1;i<nrow;i++) {
	    oldtyp=0;
	    /*Check first col*/
	    if(img_in[curr]) {
		    if(img_in[above]) {
			img_out[curr]=img_out[above];
			oldtyp=3;
			}
		    else {
			img_out[curr]=curr+1;
			peakcount++;
			oldtyp=1;
			}
		    }
	    above++; curr++; left++;
	    /*Interior pixels*/
	    for (j=1;j<ncol;j++) {
		typ=0;
		    /*printf("%d,%d,%d,%d,%d,%d,%d\n",i,j,typ,curr,above,left,img_in[curr]);*/
		if (img_in[curr]) {
		    typ=1+img_in[above]*2+img_in[left];
		    /*printf("---%d,%d,%d,%d,%d\n",curr,typ,img_in[curr],img_in[above],img_in[left]);*/
		    /* Only four possibilities for left and above */
		    switch(typ){
			case 1: img_out[curr]=curr+1;
			   /*printf("new id %d,%d,%d\n",curr,typ,oldtyp);*/
			   peakcount++; /*new droplet */
			   break;
			case 2: img_out[curr]=img_out[left];
			   break;
			case 3: img_out[curr]=img_out[above];
			   break;
			case 4: img_out[curr]=img_out[above];
			   /*printf("merge %d,%d,%d\n",curr,typ,oldtyp);*/
			   if(oldtyp==1) {
				img_out[left]=img_out[curr];
			        peakcount--; /* simple merge of two droplets */
			        }
			   if(oldtyp==2) { /*oops complicated */
			       /*printf("complicated\n");*/
			       complicated=1;
			       pos=img_out[left]-1; /*find left's minimum id*/
			       while(img_out[pos]!=pos+1) pos=img_out[pos]-1;
			       pos1=pos;
			       pos=img_out[above]-1; /*find above's minimum id*/
			       while(img_out[pos]!=pos+1) pos=img_out[pos]-1;
			       pos2=pos;
			       if(pos1<pos) pos=pos1; /*most minimum, helps img_anal*/
			       /*printf("pos1,pos2 %d,%d,%d\n",pos,pos1,pos2);*/
			       /*id them all*/
			       img_out[curr]=img_out[pos]; /*add new pixel*/
			       img_out[pos1]=img_out[pos]; /*merge droplets*/
			       img_out[pos2]=img_out[pos];
			       img_out[left]=img_out[pos]; /*change these ids too*/
			       img_out[above]=img_out[pos];
			       if(pos1!=pos2) peakcount--; /* may have merged two droplets */
			       /*might be able to handle some cases specially
			         and thus help avoid "complicated" loop at end.
			         For example track back along a list of 2's, if
			         it ends with a 1 could merge now.*/
			       }
			   break;
			   }
		    }
	       oldtyp=typ;
	       above++; curr++; left++;
	       }
	    }
	npeak[0]=peakcount;

	if(complicated) { /* may still need to merge middles of some droplets. */
	    for (curr=0;curr<ncol*nrow;curr++) {
		if(img_out[curr]){
		    /*find ultimate id*/
		    pos=img_out[curr]-1;
		    while(img_out[pos]!=pos+1) pos=img_out[pos]-1;
		    img_out[curr]=img_out[pos];
		    }
		}
	    }
}
