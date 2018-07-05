function displaypilatus(filename, ca)

if nargsin <2 
    caflag = 0;
end

[stack, im]=tiffread2(filename); 
im=double(stack.data); 
imagesc(im(1:160,:));
colorbar; axis equal tight;

if ca caxis(ca); end
