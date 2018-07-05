function [imout] = circshift_interp(xax, yax, im, newxcen, newycen)

dx=xax(2)-xax(1);
dy=yax(2)-yax(1);
szsm=size(im);

xshift = (-xax(1) - dx*szsm(2)/2 + newxcen)/dx;
yshift = (-yax(1) - dy*szsm(1)/2 + newycen)/dy;
mx = mod(xshift,1);
my = mod(yshift,1);

im=circshift(im, [floor(yshift)+1 floor(xshift)+1]);
im=conv2(im,[1-mx mx],'same'); 
imout=(conv2(im',[1-my my],'same'))';




