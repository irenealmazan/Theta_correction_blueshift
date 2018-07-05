function data = loadmda(filename, detector, print2screen, plotflag)

% data = loadmda(filename, detector, print2screen, plotflag)
%
% Stephan Hruszkewycz 1-12-09
%
% This function serves as a gateway to the Matlab executable file compiled
% from 'mda2matlab.c'.  Without the proper numbers input and output
% variables, running mda2matlab will crash Matlab causing it to exit with no
% warning.  To minimize such occurences, this function tries to catch errors
% and return gracefully to the Matlab prompt rather than calling a
% catastrophic mda2matlab command.  
%
% 'mda2matlab.c' is meant to be compiled as a Matlab executable so that APS 
% .mda files can be imported into a Matlab work environment as a matrix.  
% This function takes input arguments from the Matlab command line and 
% returns its output variables directly into Matlab.  This function can
% handle line and area scans where fluorescence spectra are recorded at 
% every point.  These spectra are identified as such because they contain
% 2000 points, much more than a normal line or area scan would contain. 
% These fluorescence scans are not loaded or returned to the Matlab 
% environment at this point in development.  
% 
% Input variables for loadmda:
% 
%   filename - Matlab string containing the name of the file.  The name 
%               must contain the entire path of the file if working outside 
%               the data folder.
%   detector - mda2matlab extracts ONE detector channel every time it is 
%               invoked.  The 'detector' input variable can be either the 
%               process varible channel number or a string containing the
%               process varible detector name, e.g.: '26idcXMAP:mca8.R0'.
%   print2screen - 'y' or 'n' flag that determines whether the screen 
%               readout from mda2matlab is displayed in Matlab. Suppressing
%               output is useful when loading large batches of .mda files.
%   plotflag - 'y' or 'n' flag that determines whether loaded data will be
%               automatically plotted in a Matlab figure.  Area scans are
%               plotted using 'surf' and can be rotated to show topology.
%
% Output variables:
% 
% 	data - This is a Matlab matrix containing the positioner and detector
%           data from the .mda file.  Line scans will return a two-
%           columned Matlab matrix where data(:,1) is the positioner and
%           data(:,2) is the detector readout.  Area scans are returned as
%           a 3D Matlab matrix where data(:,:,1) is the detector 
%           readout, data(:,:,2) is the corresponding values of the outer
%           loop positioner, and data(:,:,3) contains the positions of the
%           inner loop positioner.
% 			
% Compilation of mda2matlab:
% 	
% 	This C file is to be compiled (but not linked) using gcc (the free GNU
% 	compiler collection).  It is to be linked using the 'mex' compiler 
% 	that resides in '$MATLABDIRECTORY/bin/'.  This code calls both
% 	'mex.h' and 'mda-load.h' and utilizes functions from the mex and 
% 	mda-load libraries.  The 'mex.h' library resides in 
% 	'$MATLABDIRECTORY/extern/include/' and the mda-load library developed by
% 	Dohn Arms is avaiable at http://sector7.xor.aps.anl.gov/mda/ . 
% 	
% 	Compilation succeded using gcc 4.0.1 on Mac Dual Intel OS 10.5.5 
% 	running Matlab version 7.6.0.324 (R2008a)
% 	
% 	gcc mda2matlab.c -I/Applications/MATLAB_R2008a/extern/include 
% 		-I MDAUTILSDIRECTORY/ -c
% 
% 	/Applications/MATLAB_R2008a/bin/mex CC=gcc LD=gcc 
% 		-L MDAUTILSDIRECTORY/libmda-load.a -lm -output mda2matlab
% 		mda2matlab.o
% 		
% WARNING: Matlab executables compiled from external source code theoretically
% support exception handling that should return the function to the Matlab 
% prompt, however in practice these catches rarely work.  Especially in the 
% case of passing incorrect numbers of input or output variables, Matlab shuts
% down without warning and all workspace variables are lost.  This terribly 
% frustrating situation can be partially avoided by using the Matlab file 
% 'loadmda.m'.  This function serves as a simple calling function for 
% mda2matlab that tries to catch errors and return gracefully to the Matlab
% command line rather than calling a doomed mda2matlab command.  
% 
% IT IS HIGHLY RECOMMENDED THAT USERS CALL 'LOADMDA.M' WITHIN MATLAB RATHER 
% THAN DIRECTLY CALLING MDA2MATLAB TO AVOID CATASTROPHIC CRASHES!
%
% For more info on mda2matlab, see header in 'mda2matlab.c'.

fid = fopen(filename);

if(fid==-1) 
    error('The MDA file specified does not exist');
end

if(nargin<2) 
    detector=0; 
    print2screen=1;
    plotflag=1;
end

if(nargin<3) 
    print2screen=1;
    plotflag=1;
end

if(nargin<4) 
    plotflag=1;
end

if(fid~=-1)
    
    fclose(fid); %it will be re-opened by mda2matlab, here we only need 
            %to see if it exists.
    
    [data, printout]=mda2matlab(filename, detector);

    %print MDA info to screen if requested by user
    if(print2screen==1 || print2screen(1) =='y') 
        fprintf(printout);
    end


    %plot MDA results if requested by user
    if(plotflag==1 || plotflag(1) == 'y')

        if ndims(data)==2 & max(size(data))>1
            plot(data(:,1), data(:,2));
            
            %grab x label info out of printout
            posind = findstr(printout, 'positioners:');
            posind = posind+20;
            xlab = printout(posind);
            i=posind+1;
            while(printout(i) ~= ',')
                xlab = [xlab printout(i)];
                i=i+1;
            end
            xlabel(xlab);
            
            %grab y label info out of printout
            detind = findstr(printout, 'found detector:');
            detind = detind+6;
            ylab = printout(detind);
            i=detind+1;
            while(printout(i) ~= ',')
                ylab = [ylab printout(i)];
                i=i+1;
            end
            ylabel(ylab);
            
        end

        if ndims(data)==3 & max(size(data))>1

            x=data(1,:,3); 
            y=data(:,1,2); 
            
            %h=surf(x,y,data(:,:,1), 'EdgeColor', 'none');
            %view([0 0 1]);
            %lighting phong;
            %axis tight;
            
            imagesc(x,y,data(:,:,1));
            set(gca, 'Ydir', 'normal');
            
            %grab x label info out of printout
            posind2 = findstr(printout, 'positioners:');
            posind = posind2(1)+20;
            xlab = printout(posind);
            i=posind+1;
            while(printout(i) ~= ',')
                xlab = [xlab printout(i)];
                i=i+1;
            end
            xlabel(xlab,'Interpreter', 'none');
            
            %grab y label info out of printout
            posind = posind2(2)+22;
            ylab = printout(posind);
            i=posind+1;
            while(printout(i) ~= ',')
                ylab = [ylab printout(i)];
                i=i+1;
            end
            ylabel(ylab, 'Interpreter', 'none');
            
            colorbar
                 
        end

        axis square
        set(gcf, 'PaperPosition', [.25 6.75 4 4]);

        if isa(detector, 'char')
            title([filename ' detector ' detector]);
        else
            title([filename ' detector number ' num2str(detector)], ... 
                'Interpreter', 'none');
        end

    end %end if plotflag
    
    
    if data==0 
        fprintf('\n!!!!!!!WARNING!!!!!!!\n');
        fprintf('Specified detector NOT FOUND.\n');
        fprintf('Please choose a valid detector from the list above\n');
        fprintf('No data were loaded!\n');
    end
    
end


