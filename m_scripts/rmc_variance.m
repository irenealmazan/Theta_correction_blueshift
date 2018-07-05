function [v,m,immat]= rmc_variance(fileroot, numbers)

h=waitbarpos(0, ['loading images from ' fileroot]);

lettermat=['AA' 'AB' 'AC' 'AD' 'AE' 'AF' 'AG' 'AH' 'AI' 'AJ' 'AK' 'AL' 'AM'];
lettermat=[ lettermat 'AN' 'AO' 'AP' 'AQ' 'AR' 'AS' 'AT' 'AU' 'AV' 'AW' 'AX' 'AZ'];

for i=1:min([length(numbers) length(lettermat)])
    filename = [fileroot lettermat((i-1)*2+1:(i-1)*2+2) '.im'];
    %display(filename);
    im = load(filename);
    immat(:,:,i)=im;
    dims = size(im);
    im = reshape(im, dims(1)*dims(2),1);
    
    m(i,1) = mean(im);
    %v(i,1) = var(im);
    v(i,1) = mean(im.^2)/((mean(im))^2) -1;
    
    waitbar(i/length(numbers));
end

close(h)