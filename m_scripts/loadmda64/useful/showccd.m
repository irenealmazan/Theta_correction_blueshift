function ccd=showccd(filename,logflag,frac)

if(nargin<3) frac=1; end

ccd = double(imread(filename));
sz=size(ccd);

figure(10000);
clf
set(gcf, 'Name', 'CCD');
pos=get(gcf,'Position');

set(gcf, 'Position', [pos(1:2) 800 800]);
axes('Position', [0 0 1 1]);
set(gcf, 'PaperPosition', [.25 6.75 3 3]);
set(gcf, 'Color','w');
axis off

if(nargin<2) logflag=0;end

%logflag=1; frac=.3;

if(~logflag)
    % uncomment following line to show entire image (old functionality)
        h=imagesc(ccd);axis equal tight off
        % plot a smaller ROI and a line scan of sum over y
        %subplot(2,1,1); h=imagesc(ccd(320:650, 320:650)); axis equal tight off
        %subplot(2,1,2); plot(sum(ccd(320:650,320:650))); axis square tight 
    %caxis([min(min(ccd)) max(max(ccd))*frac]);
    ca=caxis;
    %caxis([ca(1) ca(1)+frac*(abs(ca(2)-ca(1)))]);
    %caxis([730 900]);
else
    h=imagesc(log10(ccd)); axis equal tight off;
    %caxis([min(min(log10(ccd))) max(max(log10(ccd)))*frac]);
    ca=caxis;
    %caxis([ca(1) ca(1)+frac*(abs(ca(2)-ca(1)))]);
    %caxis([600 1900]);
end

text(50,50,[filename ',   max val: ' num2str(max(ccd(:)))],'Interpreter','none', 'color','w');
if (max(ccd(:)==65535)) text(50, 75, 'SATURATED', 'color','r'); end

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