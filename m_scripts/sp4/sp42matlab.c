// Written by Stephan Hruszkewycz 6-25-09 (shrus@anl.gov)
//
// This function is meant to be compiled as a Matlab executable so that  
// .sp4 files can be imported into a Matlab work environment as a matrix.  
// This function takes input arguments from the Matlab command line and 
// returns its output variables directly into Matlab.  This function can 
// handle sp4 files that contain 2- or 3-dimensional data in either real
// or complex sp4 form.  The header information contained in the sp4 
// array is read but not stored numerically in matlab.  The header values
// are written to a character stream that is passed to Matlab and can be
// displayed in the Matlab work environment.  The sp4 data values are 
// passed to a complex Matlab environment variable regardless of whether
// or not the sp4 array was real or complex.  In the case of real sp4 data,
// the Matlab variable will have zero assigned to the complex part.  
//
// Input variables:
//	filename - Matlab string containing the name of the sp4 file you wish
//	to import.
//
// Output variables:
//	data - This is a Matlab matrix containing the 2D or 3D data stored in
//	the sp4 data file.
//
// 	printout - Rather than printing directly to the screen, sp42matlab 
// 	writes all its output to a buffer which is returned as a 
// 	Matlab string.  In this way, the user can choose to view the
// 	contents by entering 'fprintf(printout)' in Matlab or suppress
// 	output in cases where many sequential data sets are read.  
//			
//
// Requirements:
//  
//	the FFTW library must be installed on your system:
// 		http://www.fftw.org/
//	the PythonPhasing C source code and header files must reside somewhere
//		on your system, but do not need to be compiled (we do our own 
//		compilation here).  This software suite is maintained by Ross
//		Harder (rharder@aps.anl.gov)
//	Matlab installed with mex support (all releases I have used have it). Be
//		sure you know where the mex compiler and the 'mex.h' header live.
//
// Compilation:
// 	
//	Other compilatin schemes may work for this code.  For instance, Matlab
//	allows compilation to be done from the Matlab command line.  I have only
//	ever done it as it is described here.  I have absolutely no experience 
//	trying to get this to work on a Windows machine with C compilers 
//	available on Windows.
//
// 	This C file can be compiled (but not linked) using gcc (the free GNU
// 	compiler collection).  It is to be linked using the 'mex' compiler 
// 	that resides in '$MATLABROOT/bin/'.  This code calls 'mex.h', 
//	'sp4array.h', and 'sp4util.h' and uses functions from the mex  
// 	and sp4 libraries.  The 'mex.h' library resides in 
// 	'$MATLABROOT/extern/include/' and the sp4 libraries are currently 
//  being maintained by Ross Harder (rharder@aps.anl.gov)
// 	
// 	Compilation succeded using gcc 4.0.1 on OS 10.5.5 
//	runs on Dual and Quad Core Intel chips (probably ok for all x86 arch.)
// 	running Matlab version 7.6.0.324 (R2008a)
// 	
//  gcc sp42matlab.c $SP4UTILFOLDER/Sp4*.c -I$MATLABROOT/extern/include -I$SP4UTILFOLDER 
//  	-I$PPHASEFOLDER -I$FFTWHEADER -c
//
//  $MATLABROOT/bin/mex CC=gcc LD=gcc -L $FFTWLIB/libfftw3.a -lm 
// 		-output sp42matlab sp42matlab.o Sp4*.o
// 
//  where $VARIABLE are unix-like environment variable that point to the following
//		MATLABROOT - the path leading to matlab folder
//					e.g.: /Applications/Matlab_R2008a
//					(in r2009a, these files live in a folder like this:
//					/Applications/Matlab_R2009a.app/extern/)
//		SP4UTILFOLDER - path to folder containing C source code for sp4 functions
//					like Sp4ArrayInit.c, etc. 
//		PPHASEFOLDER - path to folder containing 'sp4array.h' and 'sp4util.h'
//		FFTWHEADER - path to folder containing 'fftw3.h'
//		FFTWLIB - path to folder containing 'libfftw3.a'
//		
// a unix shell script (written for bash) that is included in this folder
// 'compilemex.sh' runs the above compilation commands.  Edit it to reflect
// the appropriate paths on your system and type:
//		bash compilemex sp42matlab 
// in your prompt, and cross your fingers. This script takes one command line
// argument which is the name of a C file sans extension to be compiled.
//
// WARNING: Matlab executables compiled from external source code theoretically
// support exception handling that should return the function to the Matlab 
// prompt, however in practice these catches rarely work.  Especially in the 
// case of passing incorrect numbers of input or output variables, Matlab shuts
// down without warning and all workspace variables are lost.  This terribly 
// frustrating situation can be partially avoided by using the Matlab file 
// 'loadsp4.m'.  This function serves as a simple calling function for 
// sp42matlab that tries to catch errors and return gracefully to the Matlab
// command line rather than calling a doomed sp42matlab command.  
// 
// IT IS HIGHLY RECOMMENDED THAT USERS CALL 'LOADSP4.M' WITHIN MATLAB RATHER 
// THAN DIRECTLY CALLING SP42MATLAB TO AVOID CATASTROPHIC CRASHES!


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <float.h>

#include "sp4array.h"
#include "sp4util.h"
#include "mex.h"

void mexFunction(int nlhs, mxArray * plhs[], 
		int nrhs, const mxArray * prhs[])
{
//these values are used to determine if the user assigned too many input or 
//output variables from the command line in Matlab.  
const int numinputargs  = 1;  //correct number of command line inputs
const int numoutputargs = 2;  //correct number of command line outputs

//declare the local pointer that will point to the spot in memory where 
//Matlab will expect a returned matrix.  
double * output_real_ptr;
double * output_imag_ptr;

int i,j,k;
int abort=0;
char *filename;
int x,y,z;
int xmax, ymax, zmax;
int datatype;
int ndim;
unsigned int nelements;
unsigned int nn[3];
int dims3d[3];
const int *dims3dptr = &dims3d[0];
mwSize buflen;
mwSize ndims;
const mwSize *dims;
int temp;
ERRORCODE stat;
ERRORSTRUCT err;
int counter =0;
char buffer[5000];
char tempbuf[200];


if (nrhs != numinputargs)
	{
	mexWarnMsgTxt("Incorrect number of input arguments"); 
	abort = 1;
	}
if (nlhs != numoutputargs)
	{
	mexWarnMsgTxt("Incorrect number of output arguments");
	abort =	1;
	}

if(!abort)
	{
	
	if( mxIsChar(prhs[0]))
		{
		//the first and only command line argument is the filename 
		//from which we read the sp4 data.
		
		buflen = mxGetNumberOfElements(prhs[0]) + 1;
		filename = mxCalloc(buflen, sizeof(char));
		temp = mxGetString(prhs[0], filename, buflen);
		sprintf(buffer, "\nloading %s into matlab\n", filename); 
		}
			
	//initialize an sp4 structure
	Sp4Array sp4local;
	Sp4Array* sp4ptr = &sp4local;
	
	err=Sp4ArrayLoad(sp4ptr, filename);
	sprintf(tempbuf,"load errors? %i %s\n", err.code, err.message);strcat(buffer, tempbuf);
	if(err.code==SUCCESS)
		{
		sprintf(tempbuf,"load successful!\n");
		}

	sprintf(tempbuf,"sp4.ndim = %i \n" , sp4ptr->ndim);strcat(buffer, tempbuf);
	for(i=0; i<sp4ptr->ndim; i++) 
		{
		sprintf(tempbuf,"dim %i size is %i\n", i, sp4ptr->nn[i]);
		strcat(buffer, tempbuf);
		}
	sprintf(tempbuf,"sp4.datatype = %i \n" , sp4ptr->datatype);strcat(buffer, tempbuf);
	sprintf(tempbuf,"sp4.comments = %s \n" , sp4ptr->comments);strcat(buffer, tempbuf);


	if( sp4ptr->ndim ==2 || (sp4ptr->ndim==3 && sp4ptr->nn[2] ==1))

	//check to see if this data is two dimensional.  2d data can either
	//have ndim =2 or ndim=3, where the third dimension is unity

		{
		sprintf(tempbuf,"loading 2d data\n");strcat(buffer, tempbuf);

		xmax = sp4ptr->nn[1];
		ymax = sp4ptr->nn[0];

		plhs[0] = mxCreateDoubleMatrix(ymax,xmax,mxCOMPLEX);
		
		output_real_ptr = mxGetPr(plhs[0]);
		output_imag_ptr = mxGetPi(plhs[0]);
					
		counter=0;		
		
		if(sp4ptr->datatype==1)
		
		//fill in the matlab pointer with an sp4 REAL array.  this array
		//takes up half the space of an sp4 COMPLEX array, as real 
		//numerical values are stacked on top of one another.  here, the 
		//matlab variable is complex.  we fill in the real mx pointer with
		//the sp4 data while filling in the mx imaginary pointer with 0.
		//we cannot simply copy the sp4 data memory block byte-for-byte
		//to the matlab block.  in two dimensions, sp4 data is stored 
		//sequentially by rows, whereas matlab interprets arrays by
		//assigning sequential byte values to rows.  
		//
		//eg array: [1 2 3]
		//          [4 5 6]
		//
		//this is stored as [1 2 3 4 5 6] in sp4, as opposed to 
		//being stored as   [1 4 2 5 3 6] in matlab.
		//hence the re-ordering below.

			{
			for(x=0; x< xmax; x++)
				{
				for(y=0; y<ymax; y++)
					{
					*(output_real_ptr+counter) 
							= sp4ptr->data[ y*xmax + x ];
	
					*(output_imag_ptr+counter) 
							= 0;
	
 					counter=counter+1;
 					}

				}
			} //end if datatype ==1
		
		if(sp4ptr->datatype ==2)
		
		//then fill in a complex matlab pointer array with a COMPLEX 
		//sp4 array.  a COMPLEX sp4 array has alternating places in 
		//memory for the real and imaginary parts of each voxel.  so,
		//the Nth pixel can be accessed in the sp4 array by:
		//    real part: data[ 2*N ]
		//    imag part: dtat[ 2*N +1 ]
		//besides the reordering described above, we also need to 
		//account for this type of sp4 byte spacing and put the 
		//appropriate components into the real and imaginary 
		//matlab pointers. (matlab must use two distinct memory 
		//blocks to store the real and imaginary parts of an array.)
		
			{
			for(x=0; x< xmax; x++)
				{
				for(y=0; y<ymax; y++)
					{
					*(output_real_ptr+counter) 
							= sp4ptr->data[ 2* (y*xmax + x) ];
	
					*(output_imag_ptr+counter) 
							= sp4ptr->data[ 2* (y*xmax + x) +1 ];
							//=0;
							
					counter=counter+1;
					}
				}
			}//end if datatype ==2
		
		}// end if data is 2 dimensional	
	

	if( sp4ptr->ndim ==3 && sp4ptr->nn[2] > 1)

	//check to see if data is three dimensional.  agian, ndim =3 is 
	//not the minimum condition here.  the third dimension must be 
	//greater that unity in order for this to truly be a three 
	//dimensional data set.

		{
 		sprintf(tempbuf,"loading 3d data\n");strcat(buffer, tempbuf);
 		dims3d[0] = sp4ptr->nn[0];
 		dims3d[1] = sp4ptr->nn[1];
 		dims3d[2] = sp4ptr->nn[2];
 		xmax = dims3d[1];
 		ymax = dims3d[0];
 		zmax = dims3d[2];
 		
 		plhs[0] = mxCreateNumericArray(3, dims3dptr, mxDOUBLE_CLASS, mxCOMPLEX);
 		
 		output_real_ptr = mxGetPr(plhs[0]);
 		output_imag_ptr = mxGetPi(plhs[0]);
 		
 		counter=0;
 		
 		if(sp4ptr->datatype==1)

		//fill in the matlab pointer with an sp4 REAL array.  this array
		//takes up half the space of an sp4 COMPLEX array, as real 
		//numerical values are stacked on top of one another.  here, the 
		//matlab variable is complex.  we fill in the real mx pointer with
		//the sp4 data while filling in the mx imaginary pointer with 0.
		//
		//we cannot simply copy the sp4 data memory block byte-for-byte
		//to the matlab block.  in three dimensions, sp4 data is stored 
		//completely differently from matlab as demonstrated below.  
		//
		//eg 3d array:                [1 2 3] [7  8  9 ]
		//(2 sequential data slices)  [4 5 6] [10 11 12]
		//
		//sp4 order:    [1 7 2 8 3 9 4 10 5 11 6 12]
		//matlab order: [1 4 2 5 3 6 7 10 8 11 9 12] 
		//hence the re-ordering below.

 			{
 			for(z=0; z<zmax; z++)
 				{
 				for(x=0; x<xmax; x++)
 					{
 					for(y=0; y<ymax; y++)
 						{
 						
 						*(output_real_ptr + counter) =
 								sp4ptr->data[ y*xmax*zmax + x*zmax + z];
 	
 						*(output_imag_ptr + counter) = 0;
 						
 						counter = counter+1;
 						}
 					}
 				}
 			}//end if datatype =1
 
 		if(sp4ptr->datatype==2)

		//then fill in a complex matlab pointer array with a COMPLEX 
		//sp4 array.  a COMPLEX sp4 array has alternating places in 
		//memory for the real and imaginary parts of each voxel.  so,
		//the Nth pixel can be accessed in the sp4 array by:
		//    real part: data[ 2*N ]
		//    imag part: dtat[ 2*N +1 ]
		//besides the reordering described above, we also need to 
		//account for this type of sp4 byte spacing and put the 
		//appropriate components into the real and imaginary 
		//matlab pointers. 
		
 			{
 			for(z=0; z<zmax; z++)
 				{
 				for(x=0; x<xmax; x++)
 					{
 					for(y=0; y<ymax; y++)
 						{
 						
 						*(output_real_ptr + counter) =
 								sp4ptr->data[ 2*(y*xmax*zmax + x*zmax + z)];
 	
 						*(output_imag_ptr + counter) =
 								sp4ptr->data[ 2*(y*xmax*zmax + x*zmax + z) +1 ];
 						
 						counter = counter+1;
 						}
 					}
 				}
 			}//end if datatype =2

		}// end if data is 3 dimensional	
	
	//pipe all the output to a Matlab string so that displaying is optional.
	plhs[1] = mxCreateString(buffer);

	} // end if not abort
	

if(abort) //the user put too many variable on one side or the other in Matlab
	{
	mexPrintf("abort calculation\n");
	
	//one of the reasons for catostrophic crashes is that the left hand side arguments
	//in the Matlab environment are left dangling.  Here we fill in all the l.h.s. 
	//variables in the Matlab environment with zeros to prevent an unannounced crash.
	for(i=0; i<nlhs; i++)
		{
		plhs[i] = mxCreateDoubleMatrix(1,1,mxREAL);
		output_real_ptr = mxGetPr(plhs[i]);
		*output_real_ptr = 0;
		}
	}

}