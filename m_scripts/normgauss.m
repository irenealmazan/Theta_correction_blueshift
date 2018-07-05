function ngauss = normgauss(x, sigma, mu)

%SH 11-10-08
%returns a normal gaussian distribution of width sigma and offset mu

ngauss = 1/(sigma * sqrt(2*pi)) * exp( -(x-mu).^2 ./(2*sigma^2));