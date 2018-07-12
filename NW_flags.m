% This scripts contains all the flags to operate NW_masterscript_BCDI.

noiseflag = 0; 
if(noiseflag) 
    display('ADDING NOISE')
else
    display('NO NOISE')
end

addNWstrain = 1; 
if(addNWstrain) 
    display('ADDING STRAIN FIELD')
else
    display('NO STRAIN')
end

newSample = 0; 
if(newSample) 
    display('MAKING A NEW SAMPLE')
else
    display('USING AN ALREADY EXISTING SAMPLE')
    
end

whichSample = 'Random'; %swith to 'Random" to use Sid's code and to 'Hexagone' for an hexagonal sample
display(['MAKING A ' whichSample ' SAMPLE'])

addNWsf = 0; 
if(addNWsf) display('ADDING STACKING FAULTS');end

addAngJitter = 3;
switch addAngJitter
    case 0
        display('NO ANGULAR JITTER')
    case 1
        display('ADDING CONSTANT ANGULAR JITTER')
    case 2
        display('ADDING NON REGULAR ANGULAR JITTER')    
    case 3
        display('ADDING AN ALREADY GENERATED ANGULAR JITTER')    
end




plotdqshift = 0; 
if(plotdqshift) 
    display('PLOTTING DQ ROCKING CURVE'); 
end

plotResults = 1;
if(plotResults)
   display('PLOTTING RESULTS') 
end

flagDebug = 0;
if (flagDebug)
    display('DEBUG MODE')
end

usesimI = 1; 
if(usesimI) 
    display('USING SIMULATED DATA')
else
    display('USING REAL DATA')
end


flagERHIOinitial = 0;
switch flagERHIOinitial
    case 0
        display('USE ER_HIO TO REFINE THE SUPPORT')
    case 1
        display('USE ER_HIO TO CREATE INITIAL GUESS FROM SCRATCH')       
end

initialGuess = 1;
switch initialGuess 
    case 0
         display('USING TRUE OBJECT');
    case 1
        display('USING ER/HIO ROUTINE TO FIND  GOOD INITIAL GUESS')
    case 2
       display('USING RANDOM INITIAL GUESS');
    case 3
        display('USING AN ALREADY GENERATED INITIAL GUESS')
end

smoothSupportFlag = 3;
switch (smoothSupportFlag)
    case 0
        display('USING A SMOOTH SUPPORT')
    case 1
        display('USING EXACT SUPPORT')
    case 2
        display('USING ER_HIO SUPPORT')
    case 3
        display('USING AN ALREADY CREATED ER_HIO SUPPORT')
end

flagContinue = 1;