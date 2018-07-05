function data = loadsp4(filename, print2screen)

% data = loadsp4(filename, {print2screen})
%   
% Stephan Hruszkewycz 6-25-09
%
% This function serves as a gateway to the Matlab executable file compiled
% from 'sp42matlab.c'.  Without the proper numbers input and output
% variables, running sp42matlab will crash Matlab causing it to exit with no
% warning.  To minimize such occurences, this function tries to catch errors
% and return gracefully to the Matlab prompt rather than calling a
% catastrophic sp42matlab command.  
% 
% Input variables:
%   filename - Matlab string containing the name of the sp4 file you wish
%   	to import.
%
%  	print2screen (optional)- Rather than printing directly to the screen, 
%       sp42matlab writes all its output to a buffer which is returned as a 
%       Matlab string.  This flag controls whether this output is
%       displayed.  0=no display, 1=print output.
%
% Output variables:
%   data - This is a Matlab matrix containing the 2D or 3D data stored in
%       the sp4 data file.
%
%
% Details regarding compilation of sp42matlab source code:
%
% sp42matlab.c is meant to be compiled as a Matlab executable so that  
% .sp4 files can be imported into a Matlab work environment as a matrix.  
% The code takes input arguments from the Matlab command line and 
% returns its output variables directly into Matlab.  This function can 
% handle sp4 files that contain 2- or 3-dimensional data in either real
% or complex sp4 form.  The header information contained in the sp4 
% array is read but not stored numerically in matlab.  The header values
% are written to a character stream that is passed to Matlab and can be
% displayed in the Matlab work environment.  The sp4 data values are 
% passed to a complex Matlab environment variable regardless of whether
% or not the sp4 array was real or complex.  In the case of real sp4 data,
% the Matlab variable will have zero assigned to the complex part.  
% 
% Requirements:
%  
%	the FFTW library must be installed on your system:
% 		http://www.fftw.org/
%	the PythonPhasing C source code and header files must reside somewhere
%		on your system, but do not need to be compiled (we do our own 
%		compilation here).  This software suite is maintained by Ross
%		Harder (rharder@aps.anl.gov)
%	Matlab installed with mex support (all releases I have used have it). Be
%		sure you know where the mex compiler and the 'mex.h' header live.
%
% Compilation:
% 	
%	Other compilatin schemes may work for this code.  For instance, Matlab
%	allows compilation to be done from the Matlab command line.  I have only
%	ever done it as it is described here.  I have absolutely no experience 
%	trying to get this to work on a Windows machine with C compilers 
%	available on Windows.
%
% 	This C file can be compiled (but not linked) using gcc (the free GNU
% 	compiler collection).  It is to be linked using the 'mex' compiler 
% 	that resides in '$MATLABROOT/bin/'.  This code calls 'mex.h', 
%	'sp4array.h', and 'sp4util.h' and uses functions from the mex  
% 	and sp4 libraries.  The 'mex.h' library resides in 
% 	'$MATLABROOT/extern/include/' and the sp4 libraries are currently 
%  being maintained by Ross Harder (rharder@aps.anl.gov)
% 	
% 	Compilation succeded using gcc 4.0.1 on Mac Dual Intel OS 10.5.5 
% 	running Matlab version 7.6.0.324 (R2008a)
% 	
%  gcc sp42matlab.c $SP4UTILFOLDER/Sp4*.c -I$MATLABROOT/extern/include -I$SP4UTILFOLDER 
%     	-I$PPHASEFOLDER -I$FFTWHEADER -c
%
%  $MATLABROOT/bin/mex CC=gcc LD=gcc -L $FFTWLIB/libfftw3.a -lm 
% 		-output sp42matlab sp42matlab.o Sp4*.o
%  
%  where $VARIABLE are unix-like environment variable that point to the following
%   	MATLABROOT - the path leading to matlab folder
%					e.g.: /Applications/Matlab_R2008a
%					(in r2009a, these files live in a folder like this:
%					/Applications/Matlab_R2009a.app/extern/)
%		SP4UTILFOLDER - path to folder containing C source code for sp4 functions
%					like Sp4ArrayInit.c, etc. 
%		PPHASEFOLDER - path to folder containing 'sp4array.h' and 'sp4util.h'
%		FFTWHEADER - path to folder containing 'fftw3.h'
%		FFTWLIB - path to folder containing 'libfftw3.a'
%		
%
% WARNING: Matlab executables compiled from external source code theoretically
% support exception handling that should return the function to the Matlab 
% prompt, however in practice these catches rarely work.  Especially in the 
% case of passing incorrect numbers of input or output variables, Matlab shuts
% down without warning and all workspace variables are lost.  This terribly 
% frustrating situation can be partially avoided by using the Matlab file 
% 'loadsp4.m'.  This function serves as a simple calling function for 
% sp42matlab that tries to catch errors and return gracefully to the Matlab
% command line rather than calling a doomed sp42matlab command.  
% 
% IT IS HIGHLY RECOMMENDED THAT USERS CALL 'LOADSP4.M' WITHIN MATLAB RATHER 
% THAN DIRECTLY CALLING SP42MATLAB TO AVOID CATASTROPHIC CRASHES!
%
% For more info on how sp42matlab works, see comments in sp42matlab.c file.


if(nargin<2) 
    print2screen=1;
end


if isa(filename, 'char')
    
    fid = fopen(filename);

    if(fid==-1) 
        error('The sp4 file specified does not exist');
    else
        fclose(fid);
        [data, printout] = sp42matlab(filename);
        if(print2screen) fprintf(printout); end
    end
    
else
    error('Filename must be a character string');
end
     
    
    