function [minval,maxval] = tiled(base,iw,ih,increm,ca,ax,logflag)
% Base = first image
% iw = number of images horizontally
% ih = number of images vertically
% increm = number of images per point, if negative only first
%          image will be displayed.  If positive, images will
%          be averaged.
%

if nargin<4
    increm = 1;
end
if nargin<5 
    ca =[];
end
if nargin<6
    ax = [1 1028 1 1028];
end
if nargin<7
    logflag=0;
end

%thisfigure = figure;
clf

hor = 1/iw;
ver = 1/ih;

%ver = 1/iw;
%hor = 1/ih;

counter=-1;
minval = 100000;
maxval = 0;
frame=1;

chipsize=[1300 1340];

for i=1:ih
    for j=1:iw
        subplot('position',[(j-1)/iw,(i-1)/ih,hor,ver]);
%        subplot('position',[(i-1)/ih,(j-1)/iw,ver,hor]);
        if increm < 0 
            counter=counter-increm;
            ccd = double(hdfread(['image' num2str(base+counter, '%1.5d') '.hdf'], '/entry1/data/data', 'Index', {[1 1],[1 1],chipsize}));
        else
            counter = counter + 1;
            ccd = double(hdfread(['image' num2str(base+counter, '%1.5d') '.hdf'], '/entry1/data/data', 'Index', {[1 1],[1 1],chipsize}));
            for k=2:increm
                counter = counter + 1;
                ccd = ccd + double(hdfread(['image' num2str(base+counter, '%1.5d') '.hdf'], '/entry1/data/data', 'Index', {[1 1],[1 1],chipsize}));
            end
            %ccd = ccd/increm;
            %lmax = max(max(ccd));
            %lmin = min(min(ccd));
        end
        
        if logflag
            imagesc(log10(ccd));
        else
            imagesc(ccd);
        end
        
        %set(gca, 'DataAspect', 
        axis off;
        %if ~isempty(ca)
        %    caxis(ca); 
        %end
        if ~isempty(ax)
            axis(ax); 
        end
        
        %F(frame) = getframe;
        %frame = frame+1;
       
        % axis equal;
    end
end

