function h=linecomp(xax, ydata)


hold on;
colormap jet;
%colormap hsv;
cm=colormap;

%for ii=1:numel(ydata)
%    h=plot(xax(ii), abs(ydata(ii)), 'o');
%    ph=angle(ydata(ii));
%    ind = round((ph+pi)/(2*pi)*63 + 1);
%    set(h,'markeredgecolor',cm(ind,:));
%end

h=[];

for ii=1:numel(ydata)-1
    
    ph=angle(ydata(ii));
    ind = round((ph+pi)/(2*pi)*63 + 1);
    
    %h=plot( xax(ii:ii+1), abs(ydata(ii:ii+1)));
    %set(h, 'color',cm(ind,:), 'linewidth', 1.5);
    
    htemp=area(xax(ii:ii+1), abs(ydata(ii:ii+1)));
    set(htemp, 'edgecolor', 'none', 'facecolor', cm(ind,:));
    h=[h htemp];
    
end
hold off;