
clear data_sim data_exp

% initialize
rock_curve = zeros(numel(delta_thscanvals),1);
dq_shift_nominal = zeros(numel(delta_thscanvals),3);
dq_shift_real = zeros(numel(delta_thscanvals),3);
mxI = zeros(size(delta_thscanvals));
im_sum_sim = zeros(Npix);

% Angular gitter:
 index_to_distort = [1:numel(delta_thscanvals)];
switch addAngJitter
    case 0
        [dth_disp] =  Phretrieval_functions.generate_angular_jitter(delta_thscanvals,index_to_distort,percent);
    case 1
        dth_disp = zeros(numel(delta_thscanvals),1);
        dth_disp(index_to_distort) = [1e-3];
    case 2
        [dth_disp] = Phretrieval_functions.generate_angular_jitter(delta_thscanvals,index_to_distort,percent); 
    case 3
        load('Original_jitter');
end


 [dq_shift_nominal] = DiffractionPatterns.calc_dqshift_for_given_th(delta_thscanvals,ki_o,kf_o,kf_o-ki_o);
     
 [dq_shift_real] = DiffractionPatterns.calc_dqshift_for_given_th(delta_thscanvals + dth_disp',ki_o,kf_o,kf_o-ki_o);
    
 [ simI,rock_curve,Proj_vol,FT_Proj_vol,Qterm] = DiffractionPatterns.calc_dp(dq_shift_real,probe,NW,X,Y,Z);

 if plotResults
     figure(27);clf;
 end
 
for ii = 1:numel(delta_thscanvals)
    
     mxI(ii) = max(simI(ii).I(:));
     
     im_sum_sim = im_sum_sim + simI(ii).I;    
    
     % store angles, dqshifts and diffraction patterns in structure
    data_exp(ii).dth_real = delta_thscanvals(ii)+dth_disp(ii);
    data_exp(ii).dth_nominal = delta_thscanvals(ii);
    data_exp(ii).dth_iter = delta_thscanvals(ii);
    data_exp(ii).dth_disp = dth_disp(ii);
    data_exp(ii).dshift_nominal = dq_shift_nominal;
    data_exp(ii).dqshift_real = dq_shift_real(ii,:);
    data_exp(ii).dqshift = dq_shift_real(ii,:); % initial value of dq
    data_exp(ii).simI = simI(ii).I;
    data_exp(ii).rock = rock_curve(ii);
    
    % plot:
    if plotResults
        if(mod(ii,1)==0)
            display(['simulating dp, ' num2str(ii) ' of ' num2str(numel(data_exp))]);
            subplot(121);
            imagecomp(squeeze(Proj_vol(ii).Psij));
            %imagecomp(NW(:,:,ii));
            title('Object')
            colorbar;
            axis image;
            
            subplot(122);
            imagesc(simI(ii).I);
            axis image;
            colorbar;
            title([' ii = ' num2str(ii)]);
            
            pause(.5);
            
            drawnow;
        end
    end

   
end

middpind = round(numel(data_exp)/2);

display(['simulating dp, ' num2str(middpind) ' of ' num2str(numel(data_exp))]);

if plotResults
    subplot(221);
    imagecomp(squeeze(Proj_vol(middpind).Psij));
    axis image;
    
    subplot(221);
    imagesc(simI(middpind).I);
    axis image;
    title([' ii = ' num2str(middpind)]);
    
    
    drawnow;
end
