function txt=click4ccd(empt, event_obj)

pos = get(event_obj, 'Position');
pixel=get(event_obj, 'DataIndex');
himage=get(event_obj, 'Target');

info = get(himage, 'UserData');
ccdnum=info.ccdnums(pixel(2),pixel(1));
mdanum=info.mdanum;

ccdnum=ccdnum-1;
hfig=gcf;

filename=['Images/' num2str(mdanum) '/scan_' num2str(mdanum) '_img_' num2str(ccdnum, '%6.6d') '.tif'];
%display(filename);
showccd(filename,0);

%filename=['Images/' num2str(mdanum) '/scan_' num2str(mdanum) '_img_' num2str(ccdnum, '%5.5d') '.hdf'];
%showccd_hdf(filename,0);

txt = { 
    ['X: ' num2str(pos(1))], ...
    ['Y: ' num2str(pos(2))], ...
    ['pixel: ' num2str(pixel) ', CCD: ' num2str(ccdnum)] ...
    };

%datacursormode off;
%pause(.5);
%datacursormode on;