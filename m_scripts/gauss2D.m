function f= gauss2D(X,Y,lam)

%SH 11-22-09
% f= gauss(X,Y,lam)
%
% X,Y, are two dimensional meshgrid type arrays that are the same size
%
% lam(1)= X0
% lam(2)= Amplitude
% lam(3)= sigma

f= lam(2) * exp(-0.5*((sqrt(X(:,:,1).^2+Y(:,:,1).^2)-lam(1))/lam(3)).^2);
