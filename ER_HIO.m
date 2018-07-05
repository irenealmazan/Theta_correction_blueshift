% This script combines ER and HIO perforemd with the 3DFT in order to solve
% the bcdi problem 
flagContinuation_ERHIO = 1;

% prepare data (matrix 3D):

%test_number = 43;

%[original_object,sup,sup_ini,dp] = TestErHioFunctions.test4p3(NW,data_exp,delta_thscanvals+dth_disp',ki_o,kf_o,probe,0.2,mncntrate,mn,X,Y,Z);

%{
figure(test_number);
clf;
for ii=1:size(dp,3)
    imagesc(dp(:,:,ii));
    axis image;
    colorbar;
    title(num2str(ii));
    pause(.2);
end
%}
er_iter = 850;
hio_iter = 200;

if flagContinuation_ERHIO == 0
   [~,dp] = DiffractionPatterns.From2DFT_to_3DFT(rho_ini,data_exp);

    [retrphase newobj] = erred3( sqrt(dp), sup, er_iter, 10, []);
else
    er_iter = 40;
    [FTPsi_ini,dp] = DiffractionPatterns.From2DFT_to_3DFT(rho_ini,data_exp);
    object_3DFT = (ifftn(fftshift(FTPsi_ini)));
    
    sup = abs(object_3DFT);
    sup_ini = sup;
    [retrphase newobj] = erred3( sqrt(dp), sup, 1, 10, []);
    finalobj = ifftn(newobj.dp);%fliplr(flipud((ifftn(newobj.dp))));
    error_to_compare(1) = DiffractionPatterns.calc_error_multiangle(probe, finalobj,data_exp ,delta_thscanvals,ki_o,kf_o,X,Y,Z) 

    
    newobj.object = fliplr(flipud((ifftn(fftshift(FTPsi_ini)))));
    
   
    
    
    for kk = 2:21
        [retrphase newobj] = erred3( sqrt(dp), sup, 1, 10, newobj);
        finalobj = ifftn(newobj.dp);%fliplr(flipud((ifftn(newobj.dp))));
        error_to_compare(kk) = DiffractionPatterns.calc_error_multiangle(probe, finalobj,data_exp ,delta_thscanvals,ki_o,kf_o,X,Y,Z) 
        display(num2str(kk))
    end
    
    for kk = [1:21]
        [retrphase newobj] = erred3( sqrt(dp), sup, 40, 10, newobj);
        finalobj = (ifftn(newobj.dp));
        error_to_compare(kk+20) = DiffractionPatterns.calc_error_multiangle(probe, finalobj,data_exp ,delta_thscanvals,ki_o,kf_o,X,Y,Z); 
        display(num2str(kk))
    end
    
     
end
[retrphase newobj] = hio3( sqrt(dp), sup_ini,hio_iter, 10, newobj, .7);

er_iter = 40;
for kk = 1:3
    [retrphase newobj] = erred3( sqrt(dp), sup, er_iter, 10, newobj);
   
end

original_object = rho_ini;
finalobj = (ifftn(newobj.dp));
%finalobj = (ifftn(newobj.dp)*exp(-1i*(0.26+1.56+2.7)));

DisplayResults.compare_two_objects(original_object,finalobj,'Original object','retrieved object',[50 80 50 80],[64],'1',10);
DisplayResults.compare_two_objects(original_object,finalobj,'Original object','retrieved object',[64 64],'12',11);

DisplayResults.compare_two_objects(original_object,finalobj,'Original object','retrieved object',[64],'2',12);
DisplayResults.compare_two_objects(original_object,finalobj,'Original object','retrieved object',[64 64],'13',13);

DisplayResults.compare_two_objects(original_object,finalobj,'Original object','retrieved object',[64],'3',14);
DisplayResults.compare_two_objects(original_object,finalobj,'Original object','retrieved object',[64 64],'23',15);

% for kk = 1:200
%     [retrphase newobj] = erred3( sqrt(dp), support_new, er_iter, 10, newobj);
%     % mod_object = abs(newobj.object);
%     %support_new = Phretrieval_functions.shrink_wrap_support(mod_object,0.1,X,Y,Z);
% end
%[retrphase newobj] = erred3( sqrt(dp), support_new, er_iter, 10, newobj);


save(['test_' num2str(test_number) ],'newobj','original_object','sup','sup_ini')

return;

%