function y = gapdf(x,a,b)
%GAMPDF	Gamma probability density function.
%	Y = GAMPDF(X,A,B) returns the gamma probability density function 
%	with parameters A and B, at the values in X.
%
%	The size of Y is the common size of the input arguments. A scalar input  
%	functions as a constant matrix of the same size as the other inputs.	 
%
%	Some references refer to the gamma distribution with a single
%	parameter. This corresponds to the default of B = 1.

%	References:
%	   [1]  L. Devroye, "Non-Uniform Random Variate Generation", 
%	   Springer-Verlag, 1986, pages 401-402.

%	Copyright (c) 1993 by The MathWorks, Inc.
%	$Revision: 1.2 $  $Date: 1993/09/10 22:00:05 $

if nargin < 3, 
    b = 1; 
end

if nargin < 2, 
    error('Requires at least two input arguments'); 
end

[errorcode x a b] = distchck(3,x,a,b);

if errorcode > 0
    error('The arguments must be the same size or be scalars.');
end

% Initialize Y to zero.
y = zeros(size(x));

%   Return NaN if the arguments are outside their respective limits.
k1 = find(a <= 0 | b <= 0);     
if any(k1)
    y(k1) = NaN * ones(size(k1));
end

k=find(x > 0 & ~(a <= 0 | b <= 0));
if any(k)
    y(k) = (a(k) - 1) .* log(x(k)) - (x(k) ./ b(k)) - gammaln(a(k)) - a(k) .* log(b(k));
    y(k) = exp(y(k));
end
k1 = find(x == 0 & a < 1);
if any(k1)
  y(k1) = Inf*ones(size(k1));
end


function [errorcode,out1,out2,out3,out4] = distchck(nparms,arg1,arg2,arg3,arg4)
%DISTCHCK Checks the argument list for the probability functions.

%	B.A. Jones  1-22-93
%	Copyright (c) 1993 by The MathWorks, Inc.
%	$Revision: 1.1 $  $Date: 1993/05/24 18:54:03 $

errorcode = 0;

if nparms == 1
    out1 = arg1;
    return;
end
    
if nparms == 2
    [r1 c1] = size(arg1);
    [r2 c2] = size(arg2);
    scalararg1 = (prod(size(arg1)) == 1);
    scalararg2 = (prod(size(arg2)) == 1);
    if ~scalararg1 & ~scalararg2
        if r1 ~= r2 | c1 ~= c2
            errorcode = 1;
            return;         
        end     
    end
    if scalararg1
        out1 = arg1(ones(r2,1),ones(c2,1));
    else
        out1 = arg1;
    end
    if scalararg2
        out2 = arg2(ones(r1,1),ones(c1,1));
    else
        out2 = arg2;
    end
end
    
if nparms == 3
    [r1 c1] = size(arg1);
    [r2 c2] = size(arg2);
    [r3 c3] = size(arg3);
    scalararg1 = (prod(size(arg1)) == 1);
    scalararg2 = (prod(size(arg2)) == 1);
    scalararg3 = (prod(size(arg3)) == 1);

    if ~scalararg1 & ~scalararg2
        if r1 ~= r2 | c1 ~= c2
            errorcode = 1;
            return;         
        end
    end

    if ~scalararg1 & ~scalararg3
        if r1 ~= r3 | c1 ~= c3
            errorcode = 1;
            return;                 
        end
    end

    if ~scalararg3 & ~scalararg2
        if r3 ~= r2 | c3 ~= c2
            errorcode = 1;
            return;         
        end
    end

    if ~scalararg1
        [rows columns] = size(arg1);
    elseif ~scalararg2
    [rows columns] = size(arg2);
    else
        [rows columns] = size(arg3);
    end
    out1 = arg1(ones(rows,1),ones(columns,1));
    out2 = arg2(ones(rows,1),ones(columns,1));
    out3 = arg3(ones(rows,1),ones(columns,1));
    out4 =[];
end

if nparms == 4
    [r1 c1] = size(arg1);
    [r2 c2] = size(arg2);
    [r3 c3] = size(arg3);
    [r4 c4] = size(arg4);
    scalararg1 = (prod(size(arg1)) == 1);
    scalararg2 = (prod(size(arg2)) == 1);
    scalararg3 = (prod(size(arg3)) == 1);
    scalararg4 = (prod(size(arg4)) == 1);

    if ~scalararg1 & ~scalararg2
        if r1 ~= r2 | c1 ~= c2
            errorcode = 1;
            return;         
        end
    end

    if ~scalararg1 & ~scalararg3
        if r1 ~= r3 | c1 ~= c3
            errorcode = 1;
            return;                 
        end
    end

    if ~scalararg1 & ~scalararg4
        if r1 ~= r4 | c1 ~= c4
            errorcode = 1;
            return;                 
        end
    end

    if ~scalararg3 & ~scalararg2
        if r3 ~= r2 | c3 ~= c2
            errorcode = 1;
            return;         
        end
    end

    if ~scalararg4 & ~scalararg2
        if r4 ~= r2 | c4 ~= c2
            errorcode = 1;
            return;         
        end
    end

    if ~scalararg3 & ~scalararg4
        if r3 ~= r4 | c3 ~= c4
            errorcode = 1;
            return;         
        end
    end


    if ~scalararg1
        [rows columns] = size(arg1);
    elseif ~scalararg2
    [rows columns] = size(arg2);
    elseif ~scalararg3
        [rows columns] = size(arg3);
    else
        [rows columns] = size(arg4);
    end
    out1 = arg1(ones(rows,1),ones(columns,1));
    out2 = arg2(ones(rows,1),ones(columns,1));
    out3 = arg3(ones(rows,1),ones(columns,1));
    out4 = arg4(ones(rows,1),ones(columns,1));
end

