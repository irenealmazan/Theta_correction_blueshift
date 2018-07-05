function ccd=showccd_hdf(filename,logflag, chipsize, frac)
%display('DESKTOP VERSION');
if(nargin<3 || isempty(chipsize))
    %have to go get the image size from the HDF header
    info = hdfinfo(filename);
    chipsize(1) = info.Vgroup(2).Vgroup.SDS.Dims(1).Size;
    chipsize(2) = info.Vgroup(2).Vgroup.SDS.Dims(2).Size;
end

if(nargin<4) frac=1; end  

ccd =hdfread(filename, '/entry1/data/data', 'Index',  {[1 1],[1 1],chipsize});
ccd=double(ccd);

figure(10000);
clf
set(gcf, 'Name', 'CCD');
pos=get(gcf,'Position');
if(pos(3)>chipsize(2))
    sizex=pos(3);
else
    sizex=chipsize(2);
end
if(pos(4)>chipsize(1))
    sizey=pos(4);
else
    sizey=chipsize(1);
end
set(gcf, 'Position', [pos(1:2) sizex sizey]);
%set(gcf, 'Position', [pos(1:2) chipsize(2) chipsize(1)]);
%subplot('Position', [0 0 1 1]);
set(gcf, 'PaperPosition', [.25 6.75 3 3]);
set(gcf, 'Color','w');

if(nargin<2) logflag=0;end

if(~logflag)
    h=imagesc(ccd);axis equal tight off
    caxis([min(min(ccd)) max(max(ccd))*frac]);
else
    h=imagesc(log10(ccd)); axis equal tight off;
    caxis([min(min(log10(ccd))) max(max(log10(ccd)))*frac]);
end

title(filename,'Interpreter','none');

set(h, 'UserData', ccd);
%datacursormode on;
dcm_obj = datacursormode(gcf);
set(dcm_obj, 'DisplayStyle', 'window');
set(dcm_obj, 'UpdateFcn', @click4value);

%----------------------------------------------
function txt=click4value(empt, event_obj)

pos = get(event_obj, 'Position');
pixel=get(event_obj, 'DataIndex');
hfig=get(event_obj, 'Target');

values = get(hfig, 'UserData');
val=values(pixel(2),pixel(1));

txt = { 
    ['pixel: ' num2str(pixel)], ...
    ['val: ' num2str(val)] ...
    };