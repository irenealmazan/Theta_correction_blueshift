global d2_bragg X Y Z ki_o kf_o

warning off;


addpath(genpath('/Users/ialmazn/Box Sync/Nanowire_ptychography/NSLS II/NSLS II March 2017/Analysis_end_of_beamtime'));
addpath(genpath('/Users/ialmazn/Documents/MATLAB/ptycho/m_scripts/'))
addpath(genpath('./calc_functions'));
addpath(genpath('./display_functions'));

%% Flags:
noiseflag = 1; 
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
if(newSample) display('MAKING A NEW SAMPLE');end
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

smoothSupportFlag = 1;
switch (smoothSupportFlag)
    case 0
        display('USING A SMOOTH SUPPORT')
    case 1
        display('USING EXACT SUPPORT')
    case 2
        display('USING THRESHOLD SUPPORT')
end


plotdqshift = 0; 
if(plotdqshift) display('PLOTTING DQ ROCKING CURVE'); end

usesimI = 1; 
if(usesimI) 
    display('USING SIMULATED DATA')
else
    display('USING REAL DATA')
end

initialGuess = 0;
switch initialGuess 
    case 0
         display('USING TRUE OBJECT');
    case 1
        display('USING ER/HIO ROUTINE TO FIND  GOOD INITIAL GUESS');
    case 2
       display('USING RANDOM INITIAL GUESS');
end

flagContinue = 0;
if flagContinue == 0
    display('STARTING A NEW PHASE RETRIEVAL OPERATION');
    
    %% initialize parameters relative to experimental setup and sample
    [corewidth,NW_length,NW_phi,a_latparam,c_latparam,...
        mplane_spacing,aplane_spacing,q_mplane,q_aplane,q_cplane,...
        zpdiam,outerzone,bsdiam,binaryprobe_flag,...
        meshdata,cutrad,edgepad,mncntrate,Niter_rho, Niter_pos,...
        Niter_theta,freq_pos,freq_restart,tau_backtrack_rho,beta_ini_rho,...
        counter_max_rho,tau_backtrack_theta,beta_ini_theta,counter_max_theta,ERflag,percent] = ...
        InitializeFunctions.NW_experimental_phretrieval_parameters();
    
    
    %% Scattering condition:
    scatgeo = 2110; %for strain image
    %scatgeo = 1010; %for stacking-fault image
    switch scatgeo
        case 1010 %SF
            [pixsize,lam,Npix,detdist,d2_bragg,depth,defocus,th,del,gam,...
                thscanvals,alphavals,phivals,...
                delta_thscanvals] = InitializeFunctions.NW_scatgeo_1010();
            %NW_scatgeo_1010;
        case 2110 %strain
            [pixsize,lam,Npix,detdist,d2_bragg,depth,defocus,th,del,gam,...
                thscanvals,alphavals,phivals,...
                delta_thscanvals] = InitializeFunctions.NW_scatgeo_2110();
            %NW_scatgeo_2110;
    end
    
    
    %% Create sample
    NW_diff_vectors_BCDI; % does both the vectors ki and kf and creates the object
    
    if newSample
        switch whichSample
            case 'Random'
                MakeSample;
            case 'Hexagone'
                NW_make_InGaAs_nocoreshell_BCDI;
        end
        
        
        NW = img;
        
    else
        
        load('Original_Sample');
    end
    
    if(addNWstrain)
        NW  = img;
    else
        NW  = abs(img);
    end
    
    
    NW_plot_diff_vectors_sample_BCDI;
    probe = ones(size(X));
    
    %% Calculate diffraction patterns
    NW_calc_dp_BCDI;
    NW_add_dp_noise;
    
    %%
    if initialGuess == 1
        ER_HIO;
    end
    
else
    display('CONTINUING A PHASE RETRIEVAL OPERATION');
end


%% Phase retrieval algorithm
NW_ph_retrieval_BCDI;
