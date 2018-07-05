function [F] = tiled(base,iw,ih,increm,ca,ax)
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
    ca ='auto';
end
if nargin<6
    ax = [350 600 100 925];
end

%thisfigure = figure;
clf

hor = .98/iw;
ver = .98/ih;

counter=-1;
minval = 100000;
maxval = 0;
frame=1;

for i=1:ih
    for j=1:iw
        %subplot('position',[(j-1)/iw,(i-1)/ih,hor,ver]);
        if increm < 0 
            counter=counter-increm;
            ccd = hdfread(['image' num2str(base+counter, '%1.5d') '.hdf'], '/entry1/data/data', 'Index', {[1 1],[1 1],[1024 1024]});
        else
            counter = counter + 1;
            ccd = hdfread(['image' num2str(base+counter, '%1.5d') '.hdf'], '/entry1/data/data', 'Index', {[1 1],[1 1],[1024 1024]});
            for k=2:increm
                counter = counter + 1;
                ccd = ccd + hdfread(['image' num2str(base+counter, '%1.5d') '.hdf'], '/entry1/data/data', 'Index', {[1 1],[1 1],[1024 1024]});
            end
            ccd = ccd/increm;
            lmax = max(max(ccd));
            lmin = min(min(ccd));
            if lmax > maxval
                maxval = lmax;
            end
            if lmin < minval
                minval = lmin;
            end
        end
        
        imagesc(ccd); 
        axis off;
        caxis(ca); 
        axis(ax); 
        
        F(frame) = getframe;
        frame = frame+1;
       
        pause(.2)
    end
end

minval;
maxval;

