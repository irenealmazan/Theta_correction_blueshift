classdef Phretrieval_functions
    % This library contains all the functions which allow us to rotate the
    % detector and the sample
    properties(Constant)
    end
    
    
    methods(Static)
        
        function [original_object,sup,dp] = prepare_data_ER_HIO(rho,data_exp)
            % this function prepares the data for test 4 conditions: 2DFT dp
            % from original object, no noise and shapr support
            
            original_object = rho;
            sup = zeros(size(rho));%abs(original_object);
           
            range = [-20:20]+65;
            
            sup(range,range,range) = ones(numel(range),numel(range),numel(range));
            
            
            for ii=1:numel(data_exp)
                dp(:,:,ii) = data_exp(ii).I;
            end
        end

         
        function [support] = make_support(corners,facet_spacing,edgepad)
            % This function calculates the dqshift connecting the diffraction
            % patterns at th = thBragg and at th = thBragg + dth
            
            
            support_edge = facet_spacing * edgepad;
            
            %make parallel planes parallel to each pair of facets
            v1 = [corners(2,:) - corners(8,:)]; v1 = v1/norm(v1);
            v2 = [corners(1,:) - corners(8,:)]; v2 = v2/norm(v2);
            v3 = cross(v1,v2); v3 = v3/norm(v3);
            T1 = v3(1)*X + v3(2)*Y + v3(3)*Z;
            T1 = (T1>-support_edge/2 & T1<support_edge/2); %two parallel lines
            
            v1 = [corners(2,:) - corners(8,:)]; v1 = v1/norm(v1);
            v2 = [corners(3,:) - corners(8,:)]; v2 = v2/norm(v2);
            v3 = cross(v1,v2); v3 = v3/norm(v3);
            T2 = v3(1)*X + v3(2)*Y + v3(3)*Z;
            T2 = (T2>-support_edge/2 & T2<support_edge/2); %two parallel lines
            
            v1 = [corners(4,:) - corners(9,:)]; v1 = v1/norm(v1);
            v2 = [corners(3,:) - corners(9,:)]; v2 = v2/norm(v2);
            v3 = cross(v1,v2); v3 = v3/norm(v3);
            T3 = v3(1)*X + v3(2)*Y + v3(3)*Z;
            T3 = (T3>-support_edge/2 & T3<support_edge/2); %two parallel lines
            
            v3 = [0 1 0];
            T4 = v3(1)*X + v3(2)*Y + v3(3)*Z;
            T4 = (T4>-support_edge/2 & T4<support_edge/2); %two parallel lines
            
            
            support = T1&T2&T3&T4;
            
            %support = T3;
            %support = abs(NW);
            
            
        end
       
         function [new_support] = shrink_wrap_support(mod_object,threshold,X,Y,Z)
            % this function updates the support by convolving the updated rho
            % with a gaussian function and then selecting the points which have
            % an intensity above a threshold
           
            
           mod_object_blur = Phretrieval_functions.smooth_support(mod_object,X,Y,Z);
           new_support = (mod_object_blur > threshold*max(mod_object_blur(:)));
            
        end
        
        function [support_new] = smooth_support(support,X,Y,Z)
            % this function updates the support by convolving the updated rho
            % with a gaussian function and then selecting the points which have
            % an intensity above a threshold
           
            
            [Npix_x,Npix_y,Npix_z] = size(X);
            
            threshold = 0.20; % 20% of the signal
            
            support_fft = fftshift(fftn(fftshift(support)));
            
            FWHM_x = 1e6;
            FWHM_y = 1e6;
            FWHM_z = 1e6;
            
            gauss_to_convolve = exp(-FWHM_x.*X.^2/Npix_x-FWHM_y.*Y.^2/Npix_y-FWHM_z.*Z.^2/Npix_z);
            
            gauss_fft = fftshift(fftn(fftshift(gauss_to_convolve)));
            
            support_fft_gauss = support_fft.* gauss_fft;
            
            support_new = fftshift(ifftn(fftshift(support_fft_gauss)));
            
            %max_of_int = max(max(max(support_conv_gauss)));
            
            %nonzeros_elements = find(support_conv_gauss> threshold*max_of_int);
            
            %new_support = zeros(Npix_x,Npix_y,Npix_z);
            %new_support(nonzeros_elements) = 1;
            
        end
        
        function [scale_fact] = ini_guess_scalefactor(probe, rho,angle_list_ini, data_exp,ki,kf,X,Y,Z)
            % This function estimates the scaling factor by which the initial guess
            % has to be calculated in order to establish good initial conditions
            % for phase retrieval.
            
            
            [dq_shift] = DiffractionPatterns.calc_dqshift_for_given_th(angle_list_ini,ki,kf,kf-ki);
            
            % modulus square of the wavefront
            [Psi_mod] = DiffractionPatterns.calc_dp( dq_shift,probe,rho,X,Y,Z);
            
            denom = 0;
            numerator = 0;
            for jj = 1:numel(Psi_mod)
               denom = denom + sum(sum(Psi_mod(jj).Psi_mod));
               numerator = numerator + sum(sum(sqrt(Psi_mod(jj).Psi_mod.*data_exp(jj).I)));
            end
            
            scale_fact = numerator/denom;
            
        end
                
        function [rho_new,beta_rho,norm_gradient,gPIEiter,direction_rho] = rho_update(probe, rho,gPIEiter_old,direction_rho_old,angles_list,support,niter, data_exp,depth,err_0,freq_restart,tau_backtrack,beta_ini,counter_max,ki,kf,X,Y,Z,ERflag)
            % this functions updates rho
            
            % gradient calculation
            [gPIEiter] = GeneralGradient.calc_grad_multiangle(probe, rho,angles_list,data_exp,ki,kf,X,Y,Z);
            
             % calculate the norm of the gradient:
            norm_gradient = gPIEiter(:)'*gPIEiter(:);   
            
            % conjugate parameter calculation
            if niter == 1 || mod(niter,freq_restart) == 0               
                gamma = 0;
            else
               norm_gradient_old = gPIEiter_old(:)'*gPIEiter_old(:);
               gamma = norm_gradient/norm_gradient_old ;
            end
            
            % direction of the scaled steepest descent:
            D = 1/(max(max(max(abs(probe).^2))));
            
            
            %       direction_rho for conjugated gradient ;
            direction_rho = (D/depth)*(-gPIEiter+gamma*direction_rho_old).*support;
            
            % calculate the adaptative step length
            [beta_rho] = GeneralGradient.calc_beta_adaptative_step(probe, rho,angles_list,data_exp,gPIEiter,err_0,direction_rho,tau_backtrack,beta_ini,counter_max,'rho',ki,kf,X,Y,Z);
            
           
            
            % update the object:
             rho_new = rho + beta_rho * direction_rho;
            
             
             
            %{
              % calculate the modulus projection:
            %Pmrho = rho + beta_rho * direction_rho;
             
            % apply constraint of support in real space
%            if ERflag ==1
% 
%                rho_new = support.*Pmrho;
%            else
%                %HIO:
%                rho_new = support.*Pmrho + (1-support).*(rho-0.9*Pmrho);
% 
%            end
             
             %}
            
            
        end
        
        function [dth_new,dq_shift, grad_final_theta,norm_grad_theta,beta] = theta_update(probe, rho,angles_list,data_exp,Niter_theta,error_0,tau_backtrack,beta_ini,counter_max,ki,kf,X,Y,Z)
            %%% this function calculates the gradient of the error metric with
            %%% respect to the position of the angles analytically, and
            %%% the correction theta  step
                        
            qbragg = kf - ki;
                     
            % initialize varibles:
            dth_new = angles_list;
            
            dq_shift = zeros(numel(data_exp),3);
            
            grad_final_theta = zeros(Niter_theta,numel(data_exp));
            
            for ntheta = 1:Niter_theta
                                
                [grad_final_theta(ntheta,:)] = GeneralGradient.calc_grad_theta(probe, rho, data_exp, dth_new,ki,kf,X,Y,Z);
                
                
                for ii = 1:numel(data_exp)%
                    grad_final_theta(ntheta,ii) = grad_final_theta(ntheta,ii)./sum(sum(data_exp(ii).I));
                    display(['dth_new = ' num2str(dth_new(ii)) 'gradient = ' num2str(grad_final_theta(ntheta,ii))]);
                end
                
                % define the direction of the descent:
                direction = -squeeze(grad_final_theta(ntheta,:)');
                
                % calculate an adaptative step size:
                [beta] = GeneralGradient.calc_beta_adaptative_step(probe, rho, dth_new,data_exp,grad_final_theta,error_0,direction,tau_backtrack,beta_ini,counter_max,'theta',ki,kf,X,Y,Z);

                % corrected theta :
                dth_new = dth_new + beta* direction;
                
                % corrected dqshift:
                [dq_shift] = DiffractionPatterns.calc_dqshift_for_given_th(dth_new,ki,kf,qbragg);
                
                % norm of the gradient:
                norm_grad_theta = grad_final_theta(ntheta,:)*grad_final_theta(ntheta,:)';
                
            end
                                    
            
            %%% TEST OF  GRADIENT
            %{
            for jj= dthsearchind(2:end-1)
                grad_manual = test_grad_theta_manually(jj,thscan,fitQerr,data_exp(ii),bmtemp,rho,grad_final_theta(ii).grad(jj));

            end
        
            %}
            
            
        end
        
        function [dth_disp] = generate_angular_jitter(delta_thscanvals,index_to_distort,percentage)
           % this function calculates a distribution of angular gitters
           % between 0 and percentage*diff(delta_thsanvals)
            
           delta_theta = unique(diff(delta_thscanvals));
           delta_theta_ini = delta_theta; % corresponds to 100 % jitter
           
           dth_disp = zeros(numel(delta_thscanvals),1);

           
           for jj =1:numel(index_to_distort)
               dth_disp(index_to_distort(jj)) = (delta_theta_ini*percentage/100)*(-1+2*rand);                
           end
            
        end
    end
        
   
        
        
end
    

