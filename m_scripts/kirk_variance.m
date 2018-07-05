function [v,m]= kirk_variance(folder, numbers)

h=waitbarpos(0, ['loading images from ' folder]);

for i=1:length(numbers)
    filename = [folder '/' num2str(numbers(i), '%5.2f') '.im'];
    %display(filename);
    im = load(filename);
    dims = size(im);
    im = reshape(im, dims(1)*dims(2),1);
    
    m(i,1) = mean(im);
    %v(i,1) = var(im);
    v(i,1) = mean(im.^2)/((mean(im))^2) -1;
    
    waitbar(i/length(numbers));
end

close(h)