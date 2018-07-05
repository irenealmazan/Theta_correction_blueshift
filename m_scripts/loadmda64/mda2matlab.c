// Written by Stephan Hruszkewycz 1-26-09
// 
// This function is meant to be compiled as a Matlab executable so that APS 
// .mda files can be imported into a Matlab work environment as a matrix.  
// This function takes input arguments from the Matlab command line and 
// returns its output variables directly into Matlab.  This function can
// handle line and area scans where fluorescence spectra are recorded at 
// every point.  These spectra are identified as such because they contain
// 2000 points, much more than a normal line or area scan would contain. 
// These fluorescence scans are not loaded or returned to the Matlab 
// environment at this point in development.  
// 
// Input variables:
// 
// 	filename - Matlab string containing the name of the file.  The name 
// 				must contain the entire path of the file if working outside 
// 				the data folder.
// 	detector - mda2matlab extracts ONE detector channel every time it is 
// 				invoked.  The 'detector' input variable can be either the 
// 				process varible channel number or a string containing the
// 				process varible detector name, e.g.: '26idcXMAP:mca8.R0'.
// 				
// Output variables:
// 
// 	data - This is a Matlab matrix containing the positioner and detector
// 			data from the .mda file.  Line scans will return a two-
// 			columned Matlab matrix where data(:,1) is the positioner and
// 			data(:,2) is the detector readout.  Area scans are returned as
// 			a 3D Matlab matrix where data(:,:,1) is the detector 
// 			readout, data(:,:,2) is the corresponding values of the outer
// 			loop positioner, and data(:,:,3) contains the positions of the
// 			inner loop positioner.
// 	printout - Rather than printing directly to the screen, mda2matlab 
// 			writes all its output to a buffer which is returned as a 
// 			Matlab string.  In this way, the user can choose to view the
// 			contents by entering 'fprintf(printout)' in Matlab or suppress
// 			output in cases where many sequential data sets are read.  
// 			
// Compilation:
// 	
// 	This C file is to be compiled (but not linked) using gcc (the free GNU
// 	compiler collection).  It is to be linked using the 'mex' compiler 
// 	that resides in '$MATLABDIRECTORY/bin/'.  This code calls both
// 	'mex.h' and 'mda-load.h' and utilizes functions from the mex and 
// 	mda-load libraries.  The 'mex.h' library resides in 
// 	'$MATLABDIRECTORY/extern/include/' and the mda-load library developed by
// 	Dohn Arms is avaiable at http://sector7.xor.aps.anl.gov/mda/ . 
// 	
// 	Compilation succeded using gcc 4.0.1 on Mac Dual Intel OS 10.5.5 
// 	running Matlab version 7.6.0.324 (R2008a)
// 	
// 	gcc mda2matlab.c -I/Applications/MATLAB_R2008a/extern/include 
// 		-I MDAUTILSDIRECTORY/ -c
// 
// 	/Applications/MATLAB_R2008a/bin/mex CC=gcc LD=gcc 
// 		-L MDAUTILSDIRECTORY/libmda-load.a -lm -output mda2matlab
// 		mda2matlab.o
// 		
// WARNING: Matlab executables compiled from external source code theoretically
// support exception handling that should return the function to the Matlab 
// prompt, however in practice these catches rarely work.  Especially in the 
// case of passing incorrect numbers of input or output variables, Matlab shuts
// down without warning and all workspace variables are lost.  This terribly 
// frustrating situation can be partially avoided by using the Matlab file 
// 'loadmda.m'.  This function serves as a simple calling function for 
// mda2matlab that tries to catch errors and return gracefully to the Matlab
// command line rather than calling a doomed mda2matlab command.  
// 
// IT IS HIGHLY RECOMMENDED THAT USERS CALL 'LOADMDA.M' WITHIN MATLAB RATHER 
// THAN DIRECTLY CALLING MDA2MATLAB TO AVOID CATASTROPHIC CRASHES!
		

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <string.h>
#include "mda-load.h"
#include "mex.h"

//function decalrations

void detectorwarning(struct mda_scan * pscan);

int finddetectorchan(int detectorchan, int lenargv, int lendetchanstr,
			struct mda_scan * pscan, char * detectorname);

//main function

void mexFunction(int nlhs, mxArray * plhs[], 
		int nrhs, const mxArray * prhs[])
{
const int numinputargs  = 2;
const int numoutputargs = 2;
int abort=0;

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
	
int16_t datarank;	
int pts, detectorchan, fluofilesize;
int finddetname=0;
int finddetnum=0; 
int outerpts, innerpts, linescan, areascan;
int i,j,k;
int dimsize[100];  //stores the dimension sizes.  large enough for a 100 nested loop scan
FILE *fptr;
struct mda_file *mda;
char buffer[5000];
char tempbuf[200];
int bufcnt = 0;
int founddetector = 0;
double * outputmat;
int scannedptsnum;

char * filename;
char * detectorname;
mwSize buflen;
mwSize dims[3];
int temp;

if(!abort)
	{

	
	if( mxIsChar(prhs[0]) )
		{
		buflen = mxGetNumberOfElements(prhs[0]) + 1;
		filename = mxCalloc(buflen, sizeof(char));
		temp = mxGetString(prhs[0], filename, buflen);
		sprintf(buffer, "\nfilename is %s\n", filename); 
		}
		
	if( mxIsChar(prhs[1]) )
		{
		buflen = mxGetNumberOfElements(prhs[1]) + 1;
		detectorname = mxCalloc(buflen, sizeof(char));
		temp = mxGetString(prhs[1], detectorname, buflen);
		sprintf(tempbuf,"reading detector %s\n", detectorname);strcat(buffer, tempbuf);
		finddetname = 1;
		}
		
	if( mxIsNumeric(prhs[1]) )
		{
		detectorchan = * (mxGetPr(prhs[1])) -1;
		sprintf(tempbuf, "reading detector number %i\n", detectorchan+1);strcat(buffer, tempbuf);		
		finddetnum =1;
		}
	
	
	fluofilesize = 2000;
	
	if( (fptr = fopen(filename, "r")) == NULL) 
		{
		mexWarnMsgTxt("error loading file pointer \n.");
		}
		
	if( (mda = mda_load(fptr))==NULL) 
		{
		mexWarnMsgTxt("error loading mda struct \n.");
		}
	
	fclose(fptr);
		
	datarank = mda->header->data_rank;
	
	sprintf(tempbuf, "\nscan number %i\n", 
			mda->header->scan_number);strcat(buffer, tempbuf);
	sprintf(tempbuf, "taken on %s\n", 
			mda->scan->time);strcat(buffer, tempbuf);
	sprintf(tempbuf, "\ndata is %i-dimensional\n", datarank);strcat(buffer, tempbuf);
	for(i=0; i< mda->header->data_rank; i++)
		{
		sprintf(tempbuf, "  intended scan size in dimention %i: %i \n", 
				i+1, mda->header->dimensions[i]);strcat(buffer, tempbuf);
		dimsize[i] = mda->header->dimensions[i];
		if(dimsize[i] == fluofilesize) 
			{
			sprintf(tempbuf, "     The size of dimension %i matches the number of\n", i+1); 
			strcat(buffer, tempbuf);
			sprintf(tempbuf, "     fluorescence detector chanels and is assumed to be\n"); strcat(buffer, tempbuf);
			sprintf(tempbuf, "     a fluorescence spectrum collected at every point.\n"); strcat(buffer, tempbuf);
			}
		}
	
	sprintf(tempbuf, "\nDIMENSION 1 INFO\n"); strcat(buffer, tempbuf);
	sprintf(tempbuf, "positioners: \n"); strcat(buffer, tempbuf);
	for(i=0;i<mda->scan->number_positioners; i++)
		{
		sprintf(tempbuf, "  (%i) %s  %s,\n", i+1, mda->scan->positioners[i]->name,
				mda->scan->positioners[i]->description);strcat(buffer, tempbuf);		
		}
	
	sprintf(tempbuf, "triggers: \n"); strcat(buffer, tempbuf);
	for(i=0;i<mda->scan->number_triggers; i++)
		{
		sprintf(tempbuf, "  (%i) %s,\n", i+1, mda->scan->triggers[i]->name);strcat(buffer, tempbuf);		
		}

	sprintf(tempbuf, "number of detectors: %i\n", 
			mda->scan->number_detectors);strcat(buffer, tempbuf);
	
	//determine what kind of data were taken in this scan
	linescan =0;
	areascan =0;
	if(datarank==1 || (datarank == 2 && dimsize[1]==fluofilesize)) linescan =1;
	if(!linescan && ( datarank==2 || (datarank==3 && dimsize[2]==fluofilesize)) ) areascan =1;
	
	//sprintf(tempbuf, "linescan %i, areascan %i\n", 
	//		linescan, areascan);strcat(buffer, tempbuf);
	//sprintf(tempbuf, "find num %i, find name %i\n", 
	//		finddetnum, finddetname);strcat(buffer, tempbuf);
	
// 1D scan ---------------------------------------------------------
	
	if(linescan) //then this is a 1D scan
		{
		
		pts = mda->scan->last_point;
		sprintf(tempbuf, "\n1D scan completed %d of %d points\n", pts, 
				mda->scan->requested_points); strcat(buffer, tempbuf);
		
		detectorchan = finddetectorchan(detectorchan, finddetname, finddetnum,
				mda->scan, detectorname);
		if(detectorchan >= 0) 
			{
			founddetector =1; 
			sprintf(tempbuf, "found detector: %i %s %s,\n", mda->scan->detectors[detectorchan]->number+1,
					mda->scan->detectors[detectorchan]->name,
					mda->scan->detectors[detectorchan]->description);
			strcat(buffer, tempbuf);
			}
			
		//sprintf(tempbuf, "detector chan is %i\n", detectorchan);strcat(buffer, tempbuf);
		sprintf(tempbuf,"\n");strcat(buffer, tempbuf);
	
		if(!founddetector) 
			{
			detectorwarning( mda->scan);
			plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL);
			outputmat = mxGetPr(plhs[0]);
			outputmat[0] = 0.0;
			}
		else
			{
			sprintf(tempbuf, "extracting line scan data\n\n");strcat(buffer, tempbuf);
			
			plhs[0] = mxCreateDoubleMatrix(pts, 2, mxREAL);		
			outputmat = mxGetPr(plhs[0]);
			
			for(i=0; i<pts; i++)
				{
				outputmat[i] = (mda->scan->positioners_data[0])[i];
				outputmat[i+pts] = (mda->scan->detectors_data[detectorchan])[i];
				}
			}
		}
	
	
// 2D scan ----------------------------------------------------------
	
 	if(areascan)
 		{
		sprintf(tempbuf, "\n  DIMENSION 2 INFO\n"); strcat(buffer, tempbuf);
		sprintf(tempbuf, "  positioners: \n"); strcat(buffer, tempbuf);
		for(i=0;i<mda->scan->sub_scans[0]->number_positioners; i++)
			{
			sprintf(tempbuf, "    (%i) %s  %s,\n", i+1, mda->scan->sub_scans[0]->positioners[i]->name,
					mda->scan->sub_scans[0]-> positioners[i]->description);strcat(buffer, tempbuf);		
			}

		sprintf(tempbuf, "  triggers: \n"); strcat(buffer, tempbuf);
		for(i=0;i<mda->scan->sub_scans[0]->number_triggers; i++)
			{
			sprintf(tempbuf, "    (%i) %s,\n", i+1, mda->scan->sub_scans[0]->triggers[i]->name);strcat(buffer, tempbuf);		
			}

		sprintf(tempbuf, "  number of detectors: %i\n", 
				mda->scan->sub_scans[0]->number_detectors);strcat(buffer, tempbuf);
		sprintf(tempbuf,"\n");strcat(buffer, tempbuf);
 		
		outerpts = mda->scan->last_point;
		innerpts = mda->scan->sub_scans[0]->requested_points;
		
 		scannedptsnum = innerpts * mda->scan->requested_points;
  		sprintf(tempbuf, "%i of %i subscans completed in this area scan\n", 
  				mda->scan->last_point, mda->scan->requested_points);strcat(buffer, tempbuf);
 				
		
		sprintf(tempbuf,"%i points in each line scan\n", innerpts);strcat(buffer, tempbuf);
		
		detectorchan = finddetectorchan(detectorchan, finddetname, finddetnum,
				mda->scan->sub_scans[0], detectorname);
		if(detectorchan >=0) 
			{
			founddetector =1; 
			sprintf(tempbuf, "found detector: %i %s %s,\n", 
					mda->scan->sub_scans[0]->detectors[detectorchan]->number+1,
					mda->scan->sub_scans[0]->detectors[detectorchan]->name,
					mda->scan->sub_scans[0]->detectors[detectorchan]->description);
			strcat(buffer, tempbuf);
			}
			
		if(!founddetector) 
			{
			detectorwarning( mda->scan->sub_scans[0] );
			plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL);
			outputmat = mxGetPr(plhs[0]);
			outputmat[0] = 0.0;
			}
		else
			{
			sprintf(tempbuf,"\nextracting area scan data\n\n");strcat(buffer, tempbuf);

			if(	mda->scan->last_point > 1)
				{
				dims[0]=outerpts;
				dims[1]=innerpts;
				dims[2]=3;
				plhs[0] = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL);
				outputmat = mxGetPr(plhs[0]);
		
				for(i=0; i<outerpts; i++)
					{
					for(j=0; j<innerpts; j++)
						{
						//fill in data here
						outputmat[i + j*outerpts] = 
								(mda->scan->sub_scans[outerpts-1-i]->detectors_data[detectorchan])[j];
						
						//fill in outer loop positioner here
						outputmat[i + j*outerpts + outerpts*innerpts] = 
								(mda->scan->positioners_data[0])[outerpts-1-i];
						
						//fill in inner loop positioner
						outputmat[i + j*outerpts + 2*outerpts*innerpts] = 
								(mda->scan->sub_scans[0]->positioners_data[0])[j];
						
						}
					}
				}
			else
				{
				//the scan was aborted before a single full line scan could be collected
				//the partial scan data will be reported as a 1D scan
				sprintf(tempbuf,"!!!!!PLEASE NOTE!!!!!\n");strcat(buffer, tempbuf);
				
				if((mda->scan->positioners_data[0])[0])
					{
					sprintf(tempbuf,"Only one line scan was completed in this area scan.\n");strcat(buffer, tempbuf);
					sprintf(tempbuf,"The data will be extracted and saved as a 1D scan.\n");strcat(buffer, tempbuf);
					sprintf(tempbuf,"The outer loop positioner for this scan is %f.\n", 
							(mda->scan->positioners_data[0])[0]);strcat(buffer, tempbuf);
					}
				else
					{
					sprintf(tempbuf,"Only a partial line scan was recorded in this area scan.\n");strcat(buffer, tempbuf);
					sprintf(tempbuf,"The data will be extracted and saved as a 1D scan.\n");strcat(buffer, tempbuf);
					sprintf(tempbuf,"The outer loop positioner for this scan is unknown.\n");
					strcat(buffer, tempbuf);
					}
					
				pts = mda->scan->sub_scans[0]->last_point;
				sprintf(tempbuf,"%i points in the line scan\n", pts);strcat(buffer, tempbuf);
						
				plhs[0] = mxCreateDoubleMatrix(pts, 2, mxREAL);		
				outputmat = mxGetPr(plhs[0]);
				
				for(i=0; i<pts; i++)
					{
					outputmat[i] = (mda->scan->sub_scans[0]->positioners_data[0])[i];
					outputmat[i+pts] = (mda->scan->sub_scans[0]->detectors_data[detectorchan])[i];
					}

				}//end else

			}// end else	
 		} // end if area scan
		

	mda_unload(mda);
	
	//pipe all the output to a Matlab string so that displaying is optional.
	plhs[1] = mxCreateString(buffer);
	
	} //end if !abort

if(abort) 
	{
	mexPrintf("abort load\n");
	for(i=0; i<nlhs; i++)
		{
		plhs[i] = mxCreateDoubleMatrix(1,1,mxREAL);
		outputmat = mxGetPr(plhs[i]);
		*outputmat = 0;
		}
	}
}



//function declarations

//DETECTORWARNING ----------------------------------------------

void detectorwarning(struct mda_scan * pscan)
{
int i;
mexPrintf("DETECTOR CHANNEL NOT FOUND\n");
mexPrintf("please choose another detector name or PV number from the list:\n\n");

for(i=0; i<pscan->number_detectors; i++)
	{
	mexPrintf("  name: %s,\tPV channel number: %i,\t description: %s\n", 
		 pscan->detectors[i]->name, pscan->detectors[i]->number+1,
		 pscan->detectors[i]->description);
	}

mexPrintf("\n");
}

//FINDDETECTORCHAN----------------------------------------------

int finddetectorchan(int detectorchan, int finddetname, int finddetnum,
			struct mda_scan * pscan, char * detectorname)
{

int i, counter=0;
int founddetector=0;

if(finddetname && !finddetnum) //then the user input a detector name rather than a chanel
 	{
	for(i=0; i<pscan->number_detectors; i++)
		{
		if(strcmp(pscan->detectors[i]->name, detectorname) ==0)
			{
			detectorchan = counter;
			founddetector = 1;
			}
		counter = counter+1;
		}
 	}
else //a number was input
	{
	for(i=0; i<pscan->number_detectors; i++)
		{
		if(pscan->detectors[i]->number == detectorchan)
			{
			detectorchan = counter;
			founddetector = 1;
			}
		counter = counter+1;
		}
	}

if(!founddetector) detectorchan = -1; 

return detectorchan;
}

