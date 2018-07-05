function [X, map] = aRGB2IND(RGB)

% writted by Nicolae CINDEA

m = size(RGB, 1);
n = size(RGB, 2);
X = zeros(m, n);
%map(1,:) = RGB(1, 1, :)./255;
map(1,:) = reshape(double(RGB(1,1,:)),1,3);

for i = 1:m
    for j = 1:n
        %RGBij = double(reshape(RGB(i,j,:), 1, 3)./255);
        %display(num2str(RGBij));
        
        RGBij=reshape(double(RGB(i,j,:)),1,3);
        %display(num2str(RGBij));
        
        isNotFound = true;
        k = 0;
        while isNotFound && k < size(map, 1)
            k = k + 1;
            if map(k,:) ==  RGBij
                isNotFound = false;
                %display(num2str(RGBij));
            end
        end
        if isNotFound
            map = [map; RGBij];
        end
        X(i,j) = double(k);
    end
end
map = double(map)/255;