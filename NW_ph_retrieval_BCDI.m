% Phase retrieval algorithm which combines the retrieval of the object and
% the angles simultaneously. The variation of the object and the angles is
% based on the gradient of the error metric and how far we correct in the
% given direction is determined by a linear search algorithm 

display(['set # rho iterations to ' num2str(Niter_rho) ' temp, position frequency to ' num2str(Niter_pos) 'per rho iteration'])

midsl = round(depth/2); % index of the section of the object represented in fig. 5
printind = round( [10:10:100]*(numel(data_exp)/100));


% Is this the continuation of a phase retrieval process or does this start
% from scracth?
if flagContinue == 0
    cnt_ntheta = 1;
    cnt_store = 1;
    
   
    
    % initial list of angles:
    angles_list = zeros(numel(data_exp),1);
    for ii = 1:numel(data_exp)
        angles_list(ii) = data_exp(ii).dth_iter;
        
        data_exp(ii).theta_iter_ini.theta = angles_list(ii);
        data_exp(ii).theta_iter_ini.grad_final = 0;
        data_exp(ii).theta_iter_ini.beta = beta_ini_theta;
        data_exp(ii).theta_iter_ini.dth_new_iter = angles_list(ii);
        data_exp(ii).theta_iter_ini.dqshift(:) = data_exp(ii).dqshift;
        
        
    end
    
     % support:
     switch smoothSupportFlag
         case 0
             support_ini = abs(NW);
             [support_smooth] = Phretrieval_functions.smooth_support(support_ini,X,Y,Z);
             support = abs(support_smooth/max(support_smooth(:)));
         case 1
             support = abs(NW);
         case 2
             support_3DFT = support_new;
             support_2DFT = DiffractionPatterns.From3DFT_to_2DFT(support_3DFT,angles_list,probe,ki_o,kf_o,X,Y,Z);
             support = (support_2DFT > 0.1*max(support_2DFT(:)));
         case 3             
             support_3DFT = support_new;
             support_2DFT = DiffractionPatterns.From3DFT_to_2DFT(support_3DFT,angles_list,probe,ki_o,kf_o,X,Y,Z);
             support = (support_2DFT > 0.6*max(support_2DFT(:)));
     end
     
    
    
    % initial guess and initial error:
    switch initialGuess
        case 0
            rho_ini = NW.*support;
            %[scale_fact] = Phretrieval_functions.ini_guess_scalefactor(probe, rho_ini,angles_list, data_exp,ki_o,kf_o,X,Y,Z);
            scale_fact = sqrt(mncntrate/mn);
            rho = rho_ini.*scale_fact ;
        case 1
            rho_3DFT= (ifftn(newobj.dp));
            rho_ini = DiffractionPatterns.From3DFT_to_2DFT(rho_3DFT,angles_list,probe,ki_o,kf_o,X,Y,Z);
            rho = rho_ini.*support;
        case 2
            
            %rho_ini = rand(Npix,Npix,depth).* exp(i*2*pi*rand(Npix,Npix,depth)).* support;
            rho_ini = ones(Npix,Npix,depth).* exp(i*2*pi*rand(Npix,Npix,depth)).* support;
            [scale_fact] = Phretrieval_functions.ini_guess_scalefactor(probe, rho_ini,angles_list, data_exp,ki_o,kf_o,X,Y,Z);
            rho = rho_ini.*scale_fact ;
         case 3
            rho_3DFT= (ifftn(newobj.dp));
            rho_ini = DiffractionPatterns.From3DFT_to_2DFT(rho_3DFT,angles_list,probe,ki_o,kf_o,X,Y,Z);
            rho = rho_ini.*support;            
    end
  
   
     [err_ini] = DiffractionPatterns.calc_error_multiangle(probe, rho.*support, data_exp,angles_list,ki_o,kf_o,X,Y,Z);
    
    % initial value of the gradient in rho and theta, assumed to be zero
    norm_grad_rho = zeros(Niter_rho,1);
    beta_rho = zeros(Niter_rho,1);
    norm_grad_theta = zeros(round(Niter_rho/freq_pos),1);
    beta_theta = zeros(round(Niter_rho/freq_pos),1);
    
    errlist = [err_ini];   
    fprintf('initial  error: %4.4d \n',errlist);
    
    % number of iterations
    nrho_vect = [1:Niter_rho];
    
    
        % Gradient of the initial guess: 
     [gPIEiter] = 0;
     direction_rho = 0;
    
else
        
    nrho_vect = [nrho:Niter_rho];
    
     % initial value of the gradient in rho and theta, assumed to be zero
    norm_grad_rho = [norm_grad_rho;zeros(Niter_rho,1)];
    beta_rho = [beta_rho;zeros(Niter_rho,1)];
    norm_grad_theta = [norm_grad_theta;zeros(round(Niter_rho/freq_pos),1)];
    beta_theta = [beta_theta;zeros(round(Niter_rho/freq_pos),1)];

    % remultiply rho by support: 1) if support remained unchanged, then
    % this won't make any change. 2) if the support has changed
    % (shrink-wrap) then it should provide a sort of new rho_ini
    rho = rho.*support;
end

if plotResults
    DisplayResults.show_rho_theta_update(5,errlist,rho.*support,midsl,angles_list,delta_thscanvals'+dth_disp,norm_grad_rho(1),beta_rho(1),norm_grad_theta(1),beta_theta(1),'Ini');
end


%% Iterative engine:

for nrho = nrho_vect

    tic;
   err=0;
   fprintf('PIE iter %i: ', nrho);

    %RHO ITERATIONS

    if(1)
           
  
        [rho,beta_rho(nrho),norm_grad_rho(nrho),gPIEiter,direction_rho] = Phretrieval_functions.rho_update(probe, rho,gPIEiter,direction_rho,angles_list,support, nrho, data_exp,depth,errlist(end),freq_restart,tau_backtrack_rho,beta_ini_rho,counter_max_rho,ki_o,kf_o,X,Y,Z,ERflag(nrho));
        
        [err] = DiffractionPatterns.calc_error_multiangle(probe, rho, data_exp,angles_list,ki_o,kf_o,X,Y,Z);
        fprintf('     error: %4.4d     norm_grad_rho: %4.4d \n', err,norm_grad_rho(nrho));
        errlist = [errlist err];
        
        % store the current reconstruction:
        rho_store(nrho).rho_square = rho(:,:,midsl);
        rho_store(nrho).rho_hex = squeeze(rho(100,:,:));
        rho_store(nrho).beta_rho = beta_rho(nrho);
        
        if plotResults
            DisplayResults.show_rho_theta_update(5,errlist,rho,midsl,angles_list,delta_thscanvals'+dth_disp,norm_grad_rho(1:nrho),beta_rho(1:nrho),norm_grad_theta(1:cnt_ntheta-1),beta_theta(1:cnt_ntheta-1),'rho');
        end
    end

     % THETA ANNEALING
    %%{
    tic;
    if mod(nrho,freq_pos) == 0

        [angles_list,dq_shift,grad_final_theta,norm_grad_theta(cnt_ntheta),beta_theta(cnt_ntheta)] = Phretrieval_functions.theta_update(probe, rho,angles_list,data_exp,Niter_theta,errlist(end),tau_backtrack_theta,beta_ini_theta,counter_max_theta,ki_o,kf_o,X,Y,Z,flagDebug);

        % store the updated theta list
        for ii = 1:numel(data_exp)%index_to_distort%
           data_exp(ii).dqshift(:) = dq_shift(ii,:); 
           data_exp(ii).dth_iter = angles_list(ii);
           data_exp(ii).theta_iter(cnt_ntheta).grad_final = grad_final_theta(ii);
           data_exp(ii).theta_iter(cnt_ntheta).beta = beta_theta(cnt_ntheta);
           data_exp(ii).theta_iter(cnt_ntheta).dth_new_iter = angles_list(ii);
           data_exp(ii).theta_iter(cnt_ntheta).dqshift(:) = dq_shift(ii,:);
        end

        [err] = DiffractionPatterns.calc_error_multiangle(probe, rho, data_exp,angles_list,ki_o,kf_o,X,Y,Z);
        fprintf('    error: %4.4d     norm_grad_theta: %4.4d \n', err,norm_grad_theta(cnt_ntheta));
        errlist = [errlist err];
       
        % plot
        if plotResults
            DisplayResults.show_rho_theta_update(5,errlist,rho,midsl,angles_list,delta_thscanvals'+dth_disp,norm_grad_rho(1:nrho),beta_rho(1:nrho),norm_grad_theta(1:cnt_ntheta),beta_theta(1:cnt_ntheta),'theta');
        end

        cnt_ntheta = cnt_ntheta + 1;
    end
    
    if mod(nrho,freq_store) == 0
        save('../results_files/results.mat');
        display(['saving at iteration ' num2str(nrho)])
    end

    %}

end

   



