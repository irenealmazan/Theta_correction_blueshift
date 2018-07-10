% this scripts only tests the optimal support

[err_support_01e6_20percent,support_01e6_20percent] = Phretrieval_functions.optimize_support(rho,[0.2],[1 1 1]*0.1e6,probe,data_exp,angles_list,ki_o,kf_o,X,Y,Z);
[err_support_05e6_20percent,support_05e6_20percent] = Phretrieval_functions.optimize_support(rho,[0.2],[1 1 1]*0.5e6,probe,data_exp,angles_list,ki_o,kf_o,X,Y,Z);
[err_support_1e6_20percent,support_1e6_20percent] = Phretrieval_functions.optimize_support(rho,[0.2],[1 1 1]*1e6,probe,data_exp,angles_list,ki_o,kf_o,X,Y,Z);
[err_support_2e6_20percent,support_2e6_20percent] = Phretrieval_functions.optimize_support(rho,[0.2],[1 1 1]*2e6,probe,data_exp,angles_list,ki_o,kf_o,X,Y,Z);

%%{
figure(7);
clf;

subplot(151);
di(abs(NW),-.9,'g',X,Y,Z);
title('True object');


support_flip_01e6_20percent = flipdim(flipud(fliplr(support_01e6_20percent)),3);
subplot(152);
di(support_flip_01e6_20percent,-.9,'g',X,Y,Z);
title(['FWHWM = 0.1e6, 20 % error = ' num2str(err_support_01e6_20percent)]);


support_flip_05e6_20percent = flipdim(flipud(fliplr(support_05e6_20percent)),3);
subplot(153);
di(support_flip_05e6_20percent,-.9,'g',X,Y,Z);
title(['FWHWM = 0.5e6, 20 % error = ' num2str(err_support_05e6_20percent)]);


support_flip_1e6_20percent = flipdim(flipud(fliplr(support_1e6_20percent)),3);
subplot(154);
di(support_flip_1e6_20percent,-.9,'g',X,Y,Z);
title({'FWHWM = 1e6, 20 % error = ' num2str(err_support_1e6_20percent)});

support_flip_2e6_20percent = flipdim(flipud(fliplr(support_2e6_20percent)),3);
subplot(155);
di(support_flip_2e6_20percent,-.9,'g',X,Y,Z);
title(['FWHWM = 2e6, 20 % error = ' num2str(err_support_2e6_20percent)]);
%}

figure(6);
clf;

subplot(151);
di(abs(NW),-.9,'g',X,Y,Z);
title('True object');


support_flip_01e6 = flipdim(flipud(fliplr(support_01e6)),3);
subplot(152);
di(support_flip_01e6,-.9,'g',X,Y,Z);
title(['FWHWM = 0.1e6, 10 % error = ' num2str(err_support_01e6)]);


support_flip_05e6 = flipdim(flipud(fliplr(support_05e6)),3);
subplot(153);
di(support_flip_05e6,-.9,'g',X,Y,Z);
title(['FWHWM = 0.5e6, 10 % error = ' num2str(err_support_05e6)]);


support_flip_1e6 = flipdim(flipud(fliplr(support_1e6)),3);
subplot(154);
di(support_flip_1e6,-.9,'g',X,Y,Z);
title({'FWHWM = 1e6, 10 % error = ' num2str(err_support_1e6)});

support_flip_2e6 = flipdim(flipud(fliplr(support_2e6)),3);
subplot(155);
di(support_flip_2e6,-.9,'g',X,Y,Z);
title(['FWHWM = 2e6, 10 % error = ' num2str(err_support_2e6)]);


%{
figure(7);
subplot(151);
imagesc(abs(NW(:,:,65)));
axis image;
colorbar;
title('True object');


support_flip_01e6_20percent = flipdim((fliplr(support_01e6_20percent)),3);
subplot(152);
imagesc(abs(support_flip_01e6_20percent(:,:,65)));
axis image;
colorbar;
title(['FWHWM = 0.1e6, 20 % error = ' num2str(err_support_01e6_20percent)]);


support_flip_05e6_20percent = flipdim((fliplr(support_05e6_20percent)),3);
subplot(153);
imagesc(abs(support_flip_05e6_20percent(:,:,65)));
axis image;
colorbar;
title(['FWHWM = 0.5e6, 20 % error = ' num2str(err_support_05e6_20percent)]);


support_flip_1e6_20percent = flipdim((fliplr(support_1e6_20percent)),3);
subplot(154);
imagesc(abs(support_flip_1e6_20percent(:,:,65)));
axis image;
colorbar;
title(['FWHWM = 1e6, 20 % error = ' num2str(err_support_1e6_20percent)]);

support_flip_2e6_20percent = flipdim((fliplr(support_2e6_20percent)),3);
subplot(155);
imagesc(abs(support_flip_2e6_20percent(:,:,65)));
axis image;
colorbar;
title(['FWHWM = 2e6, 20 % error = ' num2str(err_support_2e6_20percent)]);
%}

%{
figure(6);
subplot(151);
imagesc(abs(NW(:,:,65)));
axis image;
colorbar;
title('True object');

support_flip_01e6 = flipdim((fliplr(support_01e6)),3);
subplot(152);
imagesc(abs(support_flip_01e6(:,:,65)));
axis image;
colorbar;
title(['FWHWM = 0.1e6 error = ' num2str(err_support_01e6)]);

support_flip_05e6 = flipdim((fliplr(support_05e6)),3);
subplot(153);
imagesc(abs(support_flip_05e6(:,:,65)));
axis image;
colorbar;
title(['FWHWM = 0.5e6 error = ' num2str(err_support_05e6)]);


support_flip_1e6 = flipdim((fliplr(support_1e6)),3);
subplot(154);
imagesc(abs(support_flip_1e6(:,:,65)));
axis image;
colorbar;
title(['FWHWM = 1e6 error = ' num2str(err_support_1e6)]);


support_flip_2e6 = flipdim((fliplr(support_2e6)),3);
subplot(155);
imagesc(abs(support_flip_2e6(:,:,65)));
axis image;
colorbar;
title(['FWHWM = 2e6 error = ' num2str(err_support_2e6)]);
%}