function h=imagecomp(im, xax,yax, cax, scale)

if nargin<5 || isempty(scale)
    scale = 1;
else
    if numel(scale)==1
        imtemp = abs(im);
        imtemp(find(imtemp>scale*max(imtemp(:)))) = scale*max(imtemp(:));
        im = imtemp.*exp(i*angle(im));
    else
        imtemp = abs(im);
        imtemp(find(imtemp<scale(1))) = scale(1);
        imtemp(find(imtemp>scale(2))) = scale(2);
        im = imtemp.*exp(i*angle(im));
    end
end

if nargin<2 || isempty(xax)

    h=imagesc(angle(im), ...
        'AlphaData',abs(im),...
        'AlphaDataMapping', 'scaled');

else

    h=imagesc(xax, yax, angle(im), ...
        'AlphaData',abs(im),...
        'AlphaDataMapping', 'scaled');
end
set(gca, 'color','k');
set(gcf,'Inverthardcopy','off');

if nargin<4 || isempty(cax)
    caxis([-pi pi]);
else
    caxis(cax);
end