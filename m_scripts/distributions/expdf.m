function y = expdf(x,lambda)
%EXPPDF	Exponential probability density function.
%	Y = EXPPDF(X,LAMBDA) returns the exponential probability density 
%	function with parameter LAMBDA at the values in X.
%
%	The size of Y is the common size of X and LAMBDA. A scalar input   
%	functions as a constant matrix of the same size as the other input.	 

%	Reference:
%	   [1]  M. Abramowitz and I. A. Stegun, "Handbook of Mathematical
%	   Functions", Government Printing Office, 1964, 26.1.28.

%	Copyright (c) 1993 by The MathWorks, Inc.
%	$Revision: 1.1 $  $Date: 1993/05/24 18:54:14 $

if nargin <  1, 
    error('Requires at least one input argument.');
end

% Set the default mean to 1.
if nargin < 2,
    lambda = 1;
end

[errorcode x lambda] = distchck(2,x,lambda);

if errorcode > 0
    error('The arguments must be the same size or be scalars.');
end

% Initialize Y to zero.
y=zeros(size(x));

% Return NaN if LAMBDA is not positive.
k1 = find(lambda <= 0);
if any(k1) 
    y(k1) = NaN * ones(size(k1));
end

k=find(lambda > 0 & x >= 0);
if any(k),
    y(k) = exp(-x(k) ./ lambda(k)) ./ lambda(k);
end
