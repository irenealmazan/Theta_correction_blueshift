// Written by Stephan Hruszkewycz 6-25-09 (shrus@anl.gov)
//
// This function is meant to be compiled as a Matlab executable so that  
// .sp4 files can be saved from a Matlab work environment matrix.  
// This function takes input arguments from the Matlab command line and 
// returns after saving a .sp4 file in the current dir.  This function can 
// create sp4 files that contain 2- or 3-dimensional data in complex sp4 
// form.  Even if the matrix input from Matlab only contains real values,
// a complex sp4 file will be created with imaginary values set to zero.
//
// Input variables:
//	data - a 2D or 3D Matlab matrix that you wish to save. The form of 
//	the data MUST be that of a DOUBLE.  This program writes an sp4 memory
//	block that expects doubles (aka long float) at each entry.  A catch
//	for this is written in the Matlab calling function 'savesp4.m' which
//	should be used when calling this function.  matlab2sp4 should never
//	be directly called from the Matlab command line (see below).
//
//	filename - Matlab string containing the name of the sp4 file you wish
//	to create.
//
//	comment - a Matlab string that contains the comment to be written into
//	the sp4 file.  Up to 2048 characters are allowed.
//
// Output variables:
//	stat - 0 in the case of a successful sp4 save, 1 in case of fail.
//
// 	printout - Rather than printing directly to the screen, matlab2sp4 
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
// 	Compilation succeded using gcc 4.0.1 on Mac OS 10.5.5 
//	runs on Dual and Quad Core Intel chips (probably ok for all x86 arch.)
// 	running Matlab version 7.6.0.324 (R2008a)
// 	
//  gcc matlab2sp4.c $SP4UTILFOLDER/Sp4*.c -I$MATLABROOT/extern/include -I$SP4UTILFOLDER 
//  	-I$PPHASEFOLDER -I$FFTWHEADER -c
//
//  $MATLABROOT/bin/mex CC=gcc LD=gcc -L $FFTWLIB/libfftw3.a -lm 
// 		-output matlab2sp4 matlab2sp4.o Sp4*.o
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
//		bash compilemex matlab2sp4 
// in your prompt, and cross your fingers. This script takes one command line
// argument which is the name of a C file sans extension to be compiled.
//
// WARNING: Matlab executables compiled from external source code theoretically
// support exception handling that should return the function to the Matlab 
// prompt, however in practice these catches rarely work.  Especially in the 
// case of passing incorrect numbers of input or output variables, Matlab shuts
// down without warning and all workspace variables are lost.  This terribly 
// frustrating situation can be partially avoided by using the Matlab file 
// 'savesp4.m'.  This function serves as a simple calling function for 
// matlab2sp4 that tries to catch errors and return gracefully to the Matlab
// command line rather than calling a doomed matlab2sp4 command.  
// 
// IT IS HIGHLY RECOMMENDED THAT USERS CALL 'SAVESP4.M' WITHIN MATLAB RATHER 
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
//take care of some mex stuff first

//these values are used to determine if the user assigned too many input or 
//output variables from the command line in Matlab.  
const int numinputargs  = 3;  //correct number of command line inputs
const int numoutputargs = 2;  //correct number of command line outputs

//declare the local pointer that will point to the spot in memory where 
//Matlab will expect a returned matrix.  
double * output_real_ptr;

//declare the local pointers that point to the spots in memory where Matlab
//puts the command line arguments
double * datar;
double * datai;

int i,j,k;
int abort=0;
char *sp4outname;
char *comment;
int x,y,z;
int xmax, ymax, zmax;
int datatype;
int ndim;
unsigned int nelements;
unsigned int nn[3];
mwSize buflen;
mwSize ndims;
int dimsprod = 1;
const mwSize *dims;
int temp;
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

	//check to see if the value passed from Matlab is not a string.
	if( mxIsNumeric(prhs[0]) )
		{
		//dereference the pointer for the first cmd line argument from 
		//the Matlab prompt.  prhs represents pointers to the array of
		//matrices plugged into the right hand side of the command on the 
		//Matlab command.  here, we want the first one.
		
		ndims = mxGetNumberOfDimensions(prhs[0]);
		dims = mxGetDimensions(prhs[0]);
		sprintf(buffer, "number of dimensions: %i \n", ndims);
		for(i=0;i<ndims; i++) 
			{
			dimsprod=dimsprod * *(dims+i);
			sprintf(tempbuf,"dim %i size is %i\n", i+1, *(dims+i));
			strcat(buffer, tempbuf);
			}
		sprintf(tempbuf,"dimensions product is %i\n", dimsprod);strcat(buffer, tempbuf);
		datar = mxGetPr(prhs[0]);
		datai = mxGetPi(prhs[0]);
		}
	
	if( mxIsChar(prhs[1]))
		{
		//the second argument is the filename to be used to store the 
		//matrix in sp4 format.
		
		buflen = mxGetNumberOfElements(prhs[1]) + 1;
		sp4outname = mxCalloc(buflen, sizeof(char));
		temp = mxGetString(prhs[1], sp4outname, buflen);
		sprintf(tempbuf,"saving matrix to %s\n", sp4outname); strcat(buffer, tempbuf);
		}

	if( mxIsChar(prhs[2]))
		{
		//the third argument is the comment.
		
		buflen = mxGetNumberOfElements(prhs[2]) + 1;
		comment = mxCalloc(buflen, sizeof(char));
		temp = mxGetString(prhs[2], comment, buflen);
		sprintf(tempbuf,"comment is %s\n", comment); strcat(buffer, tempbuf);
		}
			
	//initialize an sp4 structure
	Sp4Array sp4local;
	Sp4Array* sp4ptr = &sp4local;
	
	if(ndims==2)
		{
		nn[0]= *(dims+0);
		nn[1]= *(dims+1);
		xmax = nn[1];
		ymax = nn[0];
		zmax = 1;
		}
		
	if(ndims==3)
		{
		nn[0]= *(dims+0);
		nn[1]= *(dims+1);
		nn[2]= *(dims+2);
		xmax = nn[1];
		ymax = nn[0];
		zmax = nn[2];
		}
	
	Sp4ArrayInit( sp4ptr, nn, ndims, DATATYPE_COMPLEX, comment);

	//we cannot simply copy the matlab data memory block byte-for-byte
	//to the sp4 block.  in three dimensions, sp4 data is stored 
	//completely differently from matlab as demonstrated below.  
	//
	//eg 3d array:                [1 2 3] [7  8  9 ]
	//(2 sequential data slices)  [4 5 6] [10 11 12]
	//
	//sp4 order:    [1 7 2 8 3 9 4 10 5 11 6 12]
	//matlab order: [1 4 2 5 3 6 7 10 8 11 9 12] 
	//
	//and in 2d: [1 2 3]
	//           [4 5 6]
	//
	//this is stored as [1 2 3 4 5 6] in sp4, as opposed to 
	//being stored as   [1 4 2 5 3 6] in matlab.
	//	
	//hence the re-ordering below.
	//
	//another significant difference between matlab data storage and
	//sp4 byte order is the way that they handle real and imag 
	//components.  a COMPLEX sp4 array has alternating places in 
	//memory for the real and imaginary parts of each voxel.  so,
	//the Nth pixel can be accessed in the sp4 array by:
	//    real part: data[ 2*N ]
	//    imag part: dtat[ 2*N +1 ]
	//besides the reordering described above, we also need to 
	//account for this type of sp4 byte spacing and write the 
	//appropriate components from the real and imaginary 
	//matlab pointers. (matlab must use two distinct memory 
	//blocks to store the real and imaginary parts of an array.)

	
	counter=0;
	for(y=0; y< ymax; y++)
		{
		for(x=0; x<xmax; x++)
			{
			for(z=0; z< zmax; z++)
				{
				
		 		sp4ptr->data[ 2*counter] = 
		 				*(datar + nn[0]*x + y + nn[0]*nn[1]*z); 				
 				if(datai)
 					sp4ptr->data[ 2*counter +1] = 
		 					*(datai + nn[0]*x + y + nn[0]*nn[1]*z);
				else sp4ptr->data[ 2*counter +1] = 0;
				
				counter=counter+1;
				
				}
			}
		}
		
	sprintf(tempbuf,"sp4test.datatype = %i \n" , sp4ptr->datatype);strcat(buffer, tempbuf);

	//start allocating memory to the left hand side matrix pointers.
	plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL);
	
	//get the address stored in plhs[0] and give it to the local output pointer
	output_real_ptr = mxGetPr(plhs[0]);
	
	//from now on, anything we do to the values stored at the address output_pointer
	//will end up showing up in the Matlab variable on the output side (left hand side).
	
	if(	Sp4ArraySave( sp4ptr, sp4outname) )
		*output_real_ptr=1;
	else
		*output_real_ptr=0;	
	
	sprintf(tempbuf,"\n");strcat(buffer, tempbuf);
	
	//pipe all the output to a Matlab string so that displaying is optional.
	plhs[1] = mxCreateString(buffer);

	} //end if not abort
	
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