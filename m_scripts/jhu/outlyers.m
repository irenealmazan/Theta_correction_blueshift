function imagenew = outlyers(image, stds)

%SH Aug 07 IMR japan.
%this function gets rid of outlying intensity points in a digital image
%usually from dust showing up during the scanning process.
% image - after using the commad imread
% stds - the number of standard deviations you want to use to determine
%   what an oulyer is

image = double(image);

s = size(image);

for(i=1:s(2)) 
    
    lowindex=s(1)*(i-1)+1;
    highindex = s(1)*i;
    image1d(lowindex:highindex, 1) = image(:,i);

end

stdev = std(image1d);
mean1 = mean(image1d);

%find the outlyers, which are probably particles of dust

for(j=1:s(2)) %from first column to last column
    
    foundhigh = find(image(:,j) > stds*stdev + mean1);
    foundlow = find(image(:,j) < mean1 - stds*stdev);
    [length(foundhigh) length(foundlow)];

    for i=1:length(foundhigh)
        image(foundhigh(i),j) = mean1;
    end

    for i=1:length(foundlow)
        image(foundlow(i), j) = mean1;
    end
end
imagenew = image;