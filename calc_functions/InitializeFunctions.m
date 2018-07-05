classdef InitializeFunctions
    % This library contains all the functions initializing the parameters
    % to be used in the different scripts
    properties(Constant)
    end
    
    
    methods(Static)
        
        function [corewidth,NW_length,NW_phi,a_latparam,c_latparam,...
                mplane_spacing,aplane_spacing,q_mplane,q_aplane,q_cplane,...
                zpdiam,outerzone,bsdiam,binaryprobe_flag,...
                meshdata,cutrad,edgepad,mncntrate,Niter_rho, Niter_pos,...
                Niter_theta,freq_pos,freq_store,freq_restart,tau_backtrack_rho,beta_ini_rho,...
                counter_max_rho,tau_backtrack_theta,beta_ini_theta,counter_max_theta,ERflag,percent_jitter] = NW_experimental_phretrieval_parameters()
            
            % In this script we initialize the values of the experimental set-up and
            % the details concerning the sample:
            
            
            
            %% sample details
            
            %NW details
            corewidth= .2; %edge-to-edge distance %NW 22
            %corewidth = 2*cosd(50)* .087; %take estimate from SEM image
            NW_length = 2;
            NW_phi=0; %misalignment of NW from vertical
            
            % real space
            a_latparam = 4.24e-4; %for InGaAs, 14.3% Ga fraction
            c_latparam = 6.93e-4;
            mplane_spacing = sind(60)*a_latparam;
            aplane_spacing = a_latparam/2;
            
            
            
            % reciprocal space
            q_mplane = 2*pi/(mplane_spacing*1e4);
            q_aplane = 2*pi/aplane_spacing *1e-4;
            q_cplane = 2*pi/c_latparam * 1e-4;
            %th_mplane110 = asind(lam/(2*mplane_spacing));
            %th_mplane220 = asind(2*lam/(2*mplane_spacing));
            %th_aplane112 = asind(lam/(2*a_latparam/2));
            
            
            %% ZP details
            zpdiam = 240;
            outerzone = .04;
            bsdiam = 50;
            binaryprobe_flag =0; %use binary probe defined at FWHM of probe
            
            %cutrad = .12;
            meshdata = 0;
            cutrad = .06;
            edgepad = 1.25;
            mncntrate = 1e3;
            
            %% Phase retrieval parameters:
            
            % Iteration parameters:
            Niter_rho = 8000;
            Niter_pos = 1;
            Niter_theta = 1;
            freq_pos = 1;
            %freq_rho = 10;
            freq_store = 100;
            
            % Beta adaptative step parameters and conjugated gradient restart parameter:
            freq_restart = 20;
            
            tau_backtrack_rho = 0.1;
            beta_ini_rho = 0.5;
            counter_max_rho = 10;
            
            tau_backtrack_theta = 0.5;
            beta_ini_theta = 0.1;
            counter_max_theta = 15;
            
            percent_jitter = 10;
            
            ERseries = [ones(Niter_rho,1);zeros(0,1)];
            ERflag = repmat(ERseries,round(Niter_rho/numel(ERseries)),1);
            
            
            
        end
        
        
        function [pixsize,lam,Npix,detdist,d2_bragg,depth,defocus,th,del,gam,...
                thscanvals,alphavals,phivals,...
                delta_thscanvals] = NW_scatgeo_2110()
            
            % This function contains the value of the experimental parameterse
            % when the [2110] reflection was measured.
            
            % experimental setup details
            pixsize = 55; %microns, for Merlin
            lam = etolambda(10400)*1e-4;
            
            Npix = 128;
            detdist = 0.529e6; % in micrometers
            d2_bragg = detdist * lam /(Npix*pixsize);
            depth = 128;%60;
            defocus = 0;
            
            % expected values of the motors for m-plane
            th = 73.3;
            del = -32.6; %in plane
            gam = 0; %out of plane
            
            
            % rocking curve scans
            if mod(depth,2)
                thscanvals =  [72.8:1.0/depth:73.3-1.0/depth 73.3 73.3+1.0/depth:1.0/depth:73.8];
            else
                thscanvals =  [72.8:1.0/depth:73.8-1.0/depth];
            end
            
            alphavals = zeros(numel(thscanvals),1);
            phivals = zeros(numel(thscanvals),1);
            delta_thscanvals = thscanvals-th;%[thscanvals-th;alphavals';phivals']';
            
            %usesimI = 1;
            
            
        end
        
        function [pixsize,lam,Npix,detdist,d2_bragg,depth,defocus,th,del,gam,...
                thscanvals,alphavals,phivals,...
                delta_thscanvals] = NW_scatgeo_1010()
            
            
            % This script contains the value of the experimental parameterse
            % when the [1010] reflection was measured.
            
            % experimental setup details
            pixsize = 55; %microns, for Merlin
            lam = etolambda(10400)*1e-4;
            
            Npix = 516;
            detdist = 0.33e6; % in micrometers
            d2_bragg = detdist * lam /(Npix*pixsize);
            depth = 100;
            defocus = 0;
            
            % expected values of the motors for m-plane
            th = -9.3;
            del = -18.6; %in plane
            gam = 0; %out of plane
            
            % rocking curve scans
            thscanvals =  [-9.42:.02:-9.13];
            alphavals = zeros(numel(thscanvals),1);
            phivals = zeros(numel(thscanvals),1);
            delta_thscanvals = thscanvals-th;%[thscanvals-th;alphavals';phivals']';
           
          
        end
        
    end   
        
end
    

