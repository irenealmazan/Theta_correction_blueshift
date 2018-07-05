function blurredobj= convolvegauss(image, sigma)

%SH 11-10-08

blurredobj = zeros(size(image));

[maxrows maxcols maxz]= size(image);

obj3d = 0;
if(maxz>1) obj3d =1;end

x=[1:maxcols];
y=[1:maxrows];
z=[1:maxz];

if(obj3d)
    [xmesh, ymesh, zmesh] = meshgrid(x,y,z);
else
    [xmesh, ymesh] = meshgrid(x,y);
end

% h=waitbar(0, 'waiting');
% counter =1;
% 
% for i=1:maxrows
%     for j=1:maxcols
%         for k=1:maxz
%         
%             if image(i,j,k)>0
%                 
%                 mux=j;
%                 muy=i;
%                 if(obj3d) muz = k; end
% 
%                 g2dx =(1/(sigma * sqrt(2*pi)) * exp(-(xmesh-mux).^2 ./(2*sigma^2)));
%                 g2dy = (1/(sigma * sqrt(2*pi)) * exp(-(ymesh-muy).^2 ./(2*sigma^2)));
% 
%                 if(obj3d)
%                     g2dz = (1/(sigma * sqrt(2*pi)) * exp(-(zmesh-muz).^2 ./(2*sigma^2)));
%                 else
%                     g2dz = ones(size(g2dx));
%                 end
% 
%                 %imagesc(g2dx(:,:,1));pause;
%                 %imagesc(g2dy(:,:,1));pause;
%                 %imagesc(g2dz(:,:,1));pause;
% 
%     %             clf
%     %             displayisosurf(g2dx.*g2dy.*g2dz, .001, 'y');
%     %             axis([0 20 0 20 0 20]);pause(.01);
% 
%     %             for(l=1:maxz) 
%     %                 scale = max(max(max(g2dx.*g2dy.*g2dz)));
%     %                 imagesc(g2dx(:,:,l).*g2dy(:,:,l).*g2dz(:,:,l)); 
%     %                 caxis([0 scale]);
%     %                 colorbar; pause(.1); 
%     %             end
% 
%                 blurredobj = blurredobj + image(i,j,k).*g2dx.*g2dy.*g2dz;
%             end
%             
%         end
%     end
% 
%     waitbar(i/maxrows);
% 
% end 
% 
% close(h);
% 
% clf;displayisosurf(blurredobj,  .001, 'y');axis([0 20 0 20 0 20]);

mux = round(maxcols/2)+1;
muy = round(maxrows/2)+1;
muz = round(maxz/2)+1;


%normal gaussian distribution, ie integral is always 1 and integrated
%intensity of the object is preserved

g2dx =(1/(sigma * sqrt(2*pi)) * exp(-(xmesh-mux).^2 ./(2*sigma^2)));
g2dy = (1/(sigma * sqrt(2*pi)) * exp(-(ymesh-muy).^2 ./(2*sigma^2)));

if(obj3d)
    g2dz = (1/(sigma * sqrt(2*pi)) * exp(-(zmesh-muz).^2 ./(2*sigma^2)));
else
    g2dz = ones(size(g2dx));
end


% gauss with peak amplitued = 1.  integrated intensity of object increases

% g2dx = exp(-.5*((xmesh-mux)/sigma) .^2);
% g2dy = exp(-.5*((ymesh-muy)/sigma) .^2);
% 
% if(obj3d)
%     g2dz = exp(-.5*((zmesh-muz)/sigma) .^2);
% else
%     g2dz = ones(size(g2dx));
% end


g3d = fftshift(g2dx.*g2dy.*g2dz);

g3dfft = fftn(g3d);
imagefft = fftn(image);

blurredobj = g3dfft.*imagefft;

blurredobj = ifftn(blurredobj);


