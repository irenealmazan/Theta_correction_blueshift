% In this script the figures for the paper are generated:

%%%%% Figure 1 %%%%%%%%%%%%%%


fig_num = 1;
slice_array = [1:floor(round(numel(angles_list)/4)):numel(angles_list)];
[Psi_mod_matrix] = FiguresForPaper.display_slice_dp(NW,probe,delta_thscanvals,slice_array,ki_o,kf_o,d2_bragg,th,X,Y,Z,fig_num);

figure(2);
clf;

for jj = 1:numel(slice_array)
    %subplot(1,numel(slice_array),jj);
    figure;
    imagesc(Psi_mod_matrix(:,:,slice_array(jj)));
    axis image;
    colorbar;
    title(num2str(slice_array(jj)));
end

fig_num = 3;
DisplayFunctions.display_diff_geom(NW,ki,kf,qbragg,fig_num,X,Y,Z);




%%%%%%%%%%%% Figure 1: low panel %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





%%%%%%%%% calculation of the errors in the real space:

err_final = FiguresForPaper.calculate_error_realspace(NW,flipdim(flipdim(rho,3),2),[65],'3',16);
err_ERHIO = FiguresForPaper.calculate_error_realspace(NW,flipdim(flipdim(rho_ini,3),2),[65],'3',17);
