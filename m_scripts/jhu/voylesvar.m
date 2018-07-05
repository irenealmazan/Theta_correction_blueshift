function [v mean_real]=voylesvar(set, norm)

%if it is a 1d vector matrix
if size(set,1)==1 | size(set,2)==1
    
    mean1 = mean(set);
    mean_real = mean1;

    %diff = 128 - mean1; % this is to bring all mean intensities to the middle
    %image1d = image1d + diff;
    
    mean2 = mean(set.*set);
    v = (mean2/(mean1^2)) -1;

elseif size(set,1)>1 & size(set,2)>1
    
    mean1 = mean(mean(set));
    mean_real = mean1;
    
    mean2 = mean(mean(set.*set));
    
    v = (mean2/(mean1^2)) -1;
    
end