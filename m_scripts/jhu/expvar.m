function [v mean1 vstat vlines] = expvar(image2d)

%SH 5/19/06
%This function returns the normalized variance of the image to produce a
%V(k) plot.  V = -1 + <I2> / <I>2 
%The image is also plotted in a matlab window

%display the image
%pcolor(double(image2d));
%set(findobj('Type','surface'), 'EdgeColor', 'none');
%colormap gray;
%axis equal;

%create a huge 1 column matrix with all the data
s = size(image2d);
image2d = double(image2d);

for i=1:s(1)
    line = image2d(i,:);
    ml(i) = mean(line);
    stdl = std(line);
    vl(i) = (stdl/ml(i))^2;
end
vlines = mean(vl);
mean(ml);

for(i=1:s(2))
    
    lowindex=s(1)*(i-1)+1;
    highindex = s(1)*i;
    image1d(lowindex:highindex, 1) = image2d(:,i);

end

image1d = double(image1d);
stdev = std(image1d);
mean1 = mean(image1d);

%diff = 128 - mean1; % this is to bring all mean intensities to the middle
%image1d = image1d + diff;


size(image1d);
hist(double(image1d),256);
mean1 = mean(image1d);
mean2 = mean(image1d.*image1d);
v = mean2/(mean1^2) -1;
vstat = (std(image1d)/mean1)^2;