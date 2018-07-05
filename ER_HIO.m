% This script combines ER and HIO perforemd with the 3DFT in order to solve
% the bcdi problem 
flagContinuation_ERHIO = 0;

% prepare data (matrix 3D):

[~,sup,sup_ini,dp] = Phretrieval_functions.prepare_data_ER_HIO(rho,data_exp);
errlist_ERHIO = [];

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




er_iter = 20;
hio_iter = 150;

[retrphase,newobj] = erred3( sqrt(dp), sup, er_iter, 10, []);
mod_object = abs(newobj.object);
support_new = Phretrieval_functions.shrink_wrap_support(mod_object,0.1,X,Y,Z);

for kk = 1:6
    [newobj,err] = erred3_and_calcerror_2DFT(dp,support_new,er_iter,update_iter,newobj,probe,data_exp,delta_thscanvals,ki_o,kf_o,X,Y,Z);
    
    mod_object = abs(newobj.object);
    support_new = Phretrieval_functions.shrink_wrap_support(mod_object,0.1,X,Y,Z);
end

[retrphase newobj] = hio3( sqrt(dp), support_new,hio_iter, 10, newobj, .7);

for kk = 1:3
    [retrphase newobj] = erred3( sqrt(dp), support_new, er_iter, 10, newobj);
     mod_object = abs(newobj.object);
    support_new = Phretrieval_functions.shrink_wrap_support(mod_object,0.1,X,Y,Z);
end
[retrphase newobj] = erred3( sqrt(dp), support_new, er_iter, 10, newobj);


finalobj = (ifftn(newobj.dp));



original_object = NW;
finalobj = (ifftn(newobj.dp));
%finalobj = (ifftn(newobj.dp)*exp(-1i*(0.26+1.56+2.7)));

DisplayResults.compare_two_objects(original_object,finalobj,'Original object','retrieved object',[50 80 50 80],[64],'1',10);
DisplayResults.compare_two_objects(original_object,finalobj,'Original object','retrieved object',[64 64],'12',11);

DisplayResults.compare_two_objects(original_object,finalobj,'Original object','retrieved object',[64],'2',12);
DisplayResults.compare_two_objects(original_object,finalobj,'Original object','retrieved object',[64 64],'13',13);

DisplayResults.compare_two_objects(original_object,finalobj,'Original object','retrieved object',[64],'3',14);
DisplayResults.compare_two_objects(original_object,finalobj,'Original object','retrieved object',[64 64],'23',15);

