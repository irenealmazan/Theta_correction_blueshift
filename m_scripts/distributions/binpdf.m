function y = binpdf(x,n,p)
% BINOPDF Binomial probability density function.
%	Y = BINOPDF(X,N,P) returns the binomial probability density 
%	function with parameters N and P at the values in X.
%	Note that the density function is zero unless X is an integer.
%
%	The size of Y is the common size of the input arguments. A scalar input  
%	functions as a constant matrix of the same size as the other inputs.	 

%	Reference:
%	   [1]  M. Abramowitz and I. A. Stegun, "Handbook of Mathematical
%	   Functions", Government Printing Office, 1964, 26.1.20.

%	Copyright (c) 1993 by The MathWorks, Inc.
%	$Revision: 1.1 $  $Date: 1993/05/24 18:53:34 $


if nargin < 3, 
    error('Requires three input arguments');
end

[errorcode x n p] = distchck(3,x,n,p);

if errorcode > 0
    error('The arguments must be the same size or be scalars.');
end

% Initialize Y to zero.
y = zeros(size(x));
 
% Binomial distribution is defined on positive integers less than N.
k = find(x >= 0  &  x == round(x)  &  x <= n);
if any(k),
    nk = round(exp(gammaln(n(k) + 1) - gammaln(x(k) + 1) - gammaln(n(k) - x(k) + 1)));
 y(k) = nk .* p(k) .^x(k) .* (1 - p(k)) .^ (n(k) - x(k));
end

k1 = find(n < 0 | p < 0 | p > 1 | round(n) ~= n); 
if any(k1)
    y(k1) = NaN * ones(size(k1)); 
end
