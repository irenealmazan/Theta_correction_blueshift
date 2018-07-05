function txt=click4ccd_line(empt, event_obj)

pos = get(event_obj, 'Position');
pixel=get(event_obj, 'DataIndex');
himage=get(event_obj, 'Target');

info = get(himage, 'UserData');
ccdnum=info.ccdnums(pixel(1))-1;
mdanum=info.mdanum;

hfig=gcf;

filename=['Images/' num2str(mdanum) '/scan_' num2str(mdanum) '_img_' num2str(ccdnum, '%5.5d') '.tif'];
showccd(filename,0);

txt = { 
    ['X: ' num2str(pos(1))], ...
    ['Y: ' num2str(pos(2))], ...
    ['pixel: ' num2str(pixel)], ...
    ['CCD: ' num2str(ccdnum)] ...
    };

%datacursormode off;
%pause(.5);
%datacursormode on;