% This script combines ER and HIO perforemd with the 3DFT in order to solve
% the bcdi problem 

% prepare data (matrix 3D):

[~,sup_ini,dp] = Phretrieval_functions.prepare_data_ER_HIO(NW,data_exp);
er_iter = 150;
hio_iter = 100;

[retrphase newobj] = erred3( sqrt(dp), sup_ini, er_iter, 10, []);


mod_object = abs(newobj.object);
support_new = Phretrieval_functions.shrink_wrap_support(mod_object,0.1,X,Y,Z);

er_iter = 60;
[retrphase newobj] = erred3( sqrt(dp), support_new, er_iter, 10, newobj);


er_iter = 20;
for kk = 1:3
    [retrphase newobj] = erred3( sqrt(dp), support_new, er_iter, 10, newobj);
    mod_object = abs(newobj.object);
    support_new = Phretrieval_functions.shrink_wrap_support(mod_object,0.1,X,Y,Z);
end

[retrphase newobj] = hio3( sqrt(dp), support_new,hio_iter, 10, newobj, .7);

er_iter = 60;
for kk = 1:3
    [retrphase newobj] = erred3( sqrt(dp), support_new, er_iter, 10, newobj);
     mod_object = abs(newobj.object);
    support_new = Phretrieval_functions.shrink_wrap_support(mod_object,0.1,X,Y,Z);
end
[retrphase newobj] = erred3( sqrt(dp), support_new, er_iter, 10, newobj);


finalobj = (ifftn(newobj.dp));

save('ER_HIO_initial_guess','newobj','support_new');

if plotResults
    original_object = NW*sqrt(mncntrate/mn);
    finalobj = ((ifftn((newobj.dp))));%fliplr(flipud((ifftn(newobj.dp))));
    %finalobj = (ifftn(newobj.dp)*exp(-1i*(0.26+1.56+2.7)));
    
    DisplayResults.compare_two_objects(original_object,finalobj,'Original object','retrieved object',[40 90 40 90],[64],'1',10);
    DisplayResults.compare_two_objects(original_object,finalobj,'Original object','retrieved object',[50 80],[64 64],'12',11);
    
    DisplayResults.compare_two_objects(original_object,finalobj,'Original object','retrieved object',[40 90 40 90],[64],'2',12);
    DisplayResults.compare_two_objects(original_object,finalobj,'Original object','retrieved object',[50 80],[64 64],'13',13);
    
    DisplayResults.compare_two_objects(original_object,finalobj,'Original object','retrieved object',[40 90 40 90],[64],'3',14);
    DisplayResults.compare_two_objects(original_object,finalobj,'Original object','retrieved object',[50 80],[64 64],'23',15);
    
end