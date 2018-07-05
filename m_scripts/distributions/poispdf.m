function y = poisspdf(x,lambda)
%POISSPDF Poisson probability density function.
%   Y = POISSPDF(X,LAMBDA) returns the Poisson probability density 
%   function with parameter LAMBDA at the values in X.
%
%   The size of Y is the common size of X and LAMBDA. A scalar input   
%   functions as a constant matrix of the same size as the other input.    
%
%   Note that the density function is zero unless X is an integer.
%
%   See also POISSCDF, POISSFIT, POISSINV, POISSRND, POISSTAT, PDF.

%   References:
%      [1]  M. Abramowitz and I. A. Stegun, "Handbook of Mathematical
%      Functions", Government Printing Office, 1964, 26.1.22.

%   Copyright 1993-2004 The MathWorks, Inc. 
%   $Revision: 2.13.2.5 $  $Date: 2004/12/06 16:38:01 $

if nargin <  2, 
    error('stats:poisspdf:TooFewInputs','Requires two input arguments.'); 
end

[errorcode x lambda] = distchck(2,x,lambda);

if errorcode > 0
    error('stats:poisspdf:InputSizeMismatch',...
          'Requires non-scalar arguments to match in size.');
end

if isa(x,'single') || isa(lambda,'single')
   y = zeros(size(x),'single');
else
   y = zeros(size(x));
end

if ~isfloat(x)
   x = double(x);
end

if (length(y) == 0), return; end
y(lambda < 0) = NaN;

k = (x >= 0 & x == round(x) & lambda > 0);
% using exp(gammaln) instead of gamma can avoid possible overflow.
if (any(k(:)))
   y(k) = exp(-lambda(k) + x(k) .* log(lambda(k)) - gammaln(x(k) + 1));
end

y(x==0 & lambda==0) = 1;


function [errorcode,varargout] = distchck(nparms,varargin)
%DISTCHCK Checks the argument list for the probability functions.

%   Copyright 1993-2004 The MathWorks, Inc. 
%   $Revision: 2.12.4.2 $  $Date: 2004/01/24 09:33:22 $

errorcode = 0;
varargout = varargin;

if nparms == 1
    return;
end

% Get size of each input, check for scalars, copy to output
isscalar = (cellfun('prodofsize',varargin) == 1);

% Done if all inputs are scalars.  Otherwise fetch their common size.
if (all(isscalar)), return; end

n = nparms;

for j=1:n
   sz{j} = size(varargin{j});
end
t = sz(~isscalar);
size1 = t{1};

% Scalars receive this size.  Other arrays must have the proper size.
for j=1:n
   sizej = sz{j};
   if (isscalar(j))
      vj = varargin{j};
      if isnumeric(vj)
         t = zeros(size1,class(vj));
      else
         t = zeros(size1);
      end
      t(:) = varargin{j};
      varargout{j} = t;
   elseif (~isequal(sizej,size1))
      errorcode = 1;
      return;
   end
end

