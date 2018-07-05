function savesp4(data, filename, comment, print2screen)

% savesp4(data, filename, {comment, print2screen})
%   
% Stephan Hruszkewycz 6-25-09
%
% This function serves as a gateway to the Matlab executable file compiled
% from 'matlab2sp4.c'.  Without the proper numbers input and output
% variables, running matlab2sp4 will crash Matlab causing it to exit with no
% warning.  To minimize such occurences, this function tries to catch errors
% and return gracefully to the Matlab prompt rather than calling a
% catastrophic matlab2sp4 command.  
%
%  Input variables:
%     data - a 2D or 3D Matlab matrix that you wish to save. The form of 
%     the data MUST be that of a DOUBLE.  This program writes an sp4 memory
% 	block that expects doubles (aka long float) at each entry.  A catch
% 	for this is written in the Matlab calling function 'savesp4.m' which
% 	should be used when calling this function.  matlab2sp4 should never
% 	be directly called from the Matlab command line (see below).
% 
% 	filename - Matlab string containing the name of the sp4 file you wish
% 	to create.
% 
% 	comment - a Matlab string that contains the comment to be written into
% 	the sp4 file.  Up to 2048 characters are allowed.
% 
%  Output variables:
% 	stat - 0 in the case of a successful sp4 save, 1 in case of fail.
% 
%  	printout - Rather than printing directly to the screen, matlab2sp4 
%  	writes all its output to a buffer which is returned as a 
%  	Matlab string.  In this way, the user can choose to view the
%  	contents by entering 'fprintf(printout)' in Matlab or suppress
%  	output in cases where many sequential data sets are read.  
% 			
%
% Details regarding compilation of matlab2sp4.c source code:
%
%  matlab2sp4.c is meant to be compiled as a Matlab executable so that  
%  .sp4 files can be saved from a Matlab work environment matrix.  
%  The code takes input arguments from the Matlab command line and 
%  returns after saving a .sp4 file in the current dir.  This function can 
%  create sp4 files that contain 2- or 3-dimensional data in complex sp4 
%  form.  Even if the matrix input from Matlab only contains real values,
%  a complex sp4 file will be created with imaginary values set to zero.
% 
%  Requirements:
%   
% 	the FFTW library must be installed on your system:
%  		http://www.fftw.org/
% 	the PythonPhasing C source code and header files must reside somewhere
% 		on your system, but do not need to be compiled (we do our own 
% 		compilation here).  This software suite is maintained by Ross
% 		Harder (rharder@aps.anl.gov)
% 	Matlab installed with mex support (all releases I have used have it). Be
% 		sure you know where the mex compiler and the 'mex.h' header live.
% 
%  Compilation:
%  	
% 	Other compilatin schemes may work for this code.  For instance, Matlab
% 	allows compilation to be done from the Matlab command line.  I have only
% 	ever done it as it is described here.  I have absolutely no experience 
% 	trying to get this to work on a Windows machine with C compilers 
% 	available on Windows.
% 
%  	This C file can be compiled (but not linked) using gcc (the free GNU
%  	compiler collection).  It is to be linked using the 'mex' compiler 
% 	that resides in '$MATLABROOT/bin/'.  This code calls 'mex.h', 
% 	'sp4array.h', and 'sp4util.h' and uses functions from the mex  
%  	and sp4 libraries.  The 'mex.h' library resides in 
%  	'$MATLABROOT/extern/include/' and the sp4 libraries are currently 
%     being maintained by Ross Harder (rharder@aps.anl.gov)
%  	
%  	Compilation succeded using gcc 4.0.1 on Mac OS 10.5.5 
% 	runs on Dual and Quad Core Intel chips (probably ok for all x86 arch.)
%  	running Matlab version 7.6.0.324 (R2008a)
%  	
%   gcc matlab2sp4.c $SP4UTILFOLDER/Sp4*.c -I$MATLABROOT/extern/include -I$SP4UTILFOLDER 
%   	-I$PPHASEFOLDER -I$FFTWHEADER -c
% 
%   $MATLABROOT/bin/mex CC=gcc LD=gcc -L $FFTWLIB/libfftw3.a -lm 
% 	-output matlab2sp4 matlab2sp4.o Sp4*.o
%  
%   where $VARIABLE are unix-like environment variable that point to the following
% 		MATLABROOT - the path leading to matlab folder
% 					e.g.: /Applications/Matlab_R2008a
% 					(in r2009a, these files live in a folder like this:
% 					/Applications/Matlab_R2009a.app/extern/)
% 		SP4UTILFOLDER - path to folder containing C source code for sp4 functions
% 					like Sp4ArrayInit.c, etc. 
% 		PPHASEFOLDER - path to folder containing 'sp4array.h' and 'sp4util.h'
% 		FFTWHEADER - path to folder containing 'fftw3.h'
% 		FFTWLIB - path to folder containing 'libfftw3.a'
% 		
%  a unix shell script (written for bash) that is included in this folder
%  'compilemex.sh' runs the above compilation commands.  Edit it to reflect
%  the appropriate paths on your system and type:
% 		bash compilemex matlab2sp4 
%  in your prompt, and cross your fingers. This script takes one command line
%  argument which is the name of a C file sans extension to be compiled.
% 
%  WARNING: Matlab executables compiled from external source code theoretically
%  support exception handling that should return the function to the Matlab 
%  prompt, however in practice these catches rarely work.  Especially in the 
%  case of passing incorrect numbers of input or output variables, Matlab shuts
%  down without warning and all workspace variables are lost.  This terribly 
%  frustrating situation can be partially avoided by using the Matlab file 
%  'savesp4.m'.  This function serves as a simple calling function for 
%  matlab2sp4 that tries to catch errors and return gracefully to the Matlab
%  command line rather than calling a doomed matlab2sp4 command.  
%  
%  IT IS HIGHLY RECOMMENDED THAT USERS CALL 'SAVESP4.M' WITHIN MATLAB RATHER 
%  THAN DIRECTLY CALLING SP42MATLAB TO AVOID CATASTROPHIC CRASHES!
%
% For more info on how matlab2sp4 works, see comments in sp42matlab.c file.

if(nargin<3) 
    comment ='';
    print2screen=1;
end

if(nargin<4) 
    print2screen=1;
end

%check that the data matrix is acceptable
if(~isa(data, 'numeric'))
    error('Input data must be numeric'); end

if(isempty(data))
    error('Cannot save an empty data matrix'); end

if(ndims(data)>3)
    error('Data matrix has more than three dimensions'); end

%check to see if the comment is ok
if(~isa(comment, 'char'))
    error('Comment must be a character string'); end

if(length(comment)+1 > 2048)
    error('Comment is too long. Keep it under 2048 characters.'); end

%check that the filename is ok
if(~isa(filename, 'char'))
    error('Filename must be a character string'); end

fid = fopen(filename, 'w');

if(fid==-1) 
    error('The sp4 file could not be openned.\nCheck your permissions');
else
    fclose(fid);
    
    [stat, printout] = matlab2sp4(double(data), filename, comment);
    
    if(print2screen) fprintf(printout); end
    if(stat) error('An error occured during saving.'); end
end
