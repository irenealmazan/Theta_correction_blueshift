function r = garnd(a,b,m,n);
%GAMRND	Random matrices from gamma distribution.
%	R = GAMRND(A,B) returns a matrix of random numbers chosen   
%	from the gamma distribution with parameters A and B.
%	The size of R is the common size of A and B if both are matrices.
%	If either parameter is a scalar, the size of R is the size of the other
%	parameter. Alternatively, R = GAMRND(A,B,M,N) returns an M by N matrix. 
%
%	Some references refer to the gamma distribution
%	with a single parameter. This corresponds to GAMRND
%	with B = 1. (See Devroye, pages 401-402.)

%	GAMRND uses a rejection or an inversion method depending on the
%	value of A. 

%	References:
%	   [1]  L. Devroye, "Non-Uniform Random Variate Generation", 
%	   Springer-Verlag, 1986

%	B.A. Jones 2-1-93
%	Copyright (c) 1993 by The MathWorks, Inc.
%	$Revision: 1.3 $  $Date: 1993/09/24 20:05:40 $

if nargin < 2, 
    error('Requires at least two input arguments.'); 
end


if nargin == 2
    [errorcode rows columns] = rndcheck(2,2,a,b);
end

if nargin == 3
    [errorcode rows columns] = rndcheck(3,2,a,b,m);
end

if nargin == 4
    [errorcode rows columns] = rndcheck(4,2,a,b,m,n);
end

if errorcode > 0
    error('Size information is inconsistent.');
end

% Initialize R to zero.
r = zeros(rows,columns);

scalara = (prod(size(a)) == 1);
if scalara 
    a = a(ones(rows,1),ones(columns,1));
end

scalarb = (prod(size(b)) == 1);
if scalarb 
    b = b(ones(rows,1),ones(columns,1));
end

% If a == 1, then gamma is exponential. (Devroye, page 405).
k = find(a == 1);
if any(k)
        r(k) = -b(k) .* log(rand(size(k)));
end 

c = zeros(rows,columns);
d = zeros(rows,columns);

k = find(a < 1 & a > 0);
% (Devroye, page 418 Johnk's generator)
if any(k)
    c(k) = 1 ./ a(k);
    d(k) = 1 ./ (1 - a(k));
    accept = k;
    while(length(accept)>0),
  u = rand(size(accept));
        v = rand(size(accept));
        x = u .^ c(accept);
        y = v .^ d(accept);
        k1 = find(x + y <= 1); 
        e = -log(rand(size(k1))); 
        r(accept(k1)) = e .* x(k1) ./ (x(k1) + y(k1));
        accept(k1) = [];
    end
    r(k) = r(k) .* b(k);
end

% Use a rejection method for A > 1.
k = find(a > 1);
% (Devroye, page 410 Best's algorithm)
bb = zeros(max(k),1);
c  = bb;
if any(k)
    bb(k) = a(k) - 1;
    c(k) = 3 * a(k) - 3/4;
    accept = k;
    while(length(accept)>0);
	   [m, n] = size(c(accept));
        u = rand(m,n);
        v = rand(m,n);
        w = u .* (1 - u);
        y = sqrt(c(accept) ./ w) .* (u - 0.5);
        x = bb(accept) + y;
        z = 64 * w .^ 3 .* v .^ 2;
        k1 = find(x >= 0);
        rhs = 2 * (bb(k(k1)) .* log(x(k1) ./ bb(k(k1))) - y(k1));
        k2 = find(z(k1) <= 1 - 2 * y(k1) .^2 ./ x(k1));
        r(accept(k1(k2))) = x(k1(k2)); 
        k3 = find(log(z(k1)) <= rhs & z(k1) > 1 - 2 * y(k1) .^2 ./ x(k1));
        r(accept(k1(k3))) = x(k1(k3));
        k1 = k1(:);
        omit = [k1(k2(:));k1(k3(:))];
        accept(omit) = [];
    end
        r(k) = r(k) .* b(k);
end

% Return NaN if b is not positive.
if any(any(b <= 0));
    if prod(size(b) == 1)
        r = NaN * ones(rows,columns);
    else
        k = find(b <= 0);
        r(k) = NaN * ones(size(k));
    end
end

% Return NaN if a is not positive.
if any(any(a <= 0));
    if prod(size(a) == 1)
        r = NaN * ones(rows,columns);
    else
        k = find(a <= 0);
        r(k) = NaN * ones(size(k));
    end
end

function [errorcode, rows, columns] = rndcheck(nargs,nparms,arg1,arg2,arg3,arg4,arg5)
%RNDCHECK error checks the argument list for the random number generators.

%   B.A. Jones  1-22-93
%   Copyright (c) 1993-98 by The MathWorks, Inc.
%   $Revision: 2.5 $  $Date: 1997/11/29 01:46:40 $

sizeinfo = nargs - nparms;
errorcode = 0;

if nparms == 3
    [r1 c1] = size(arg1);
    [r2 c2] = size(arg2);
    [r3 c3] = size(arg3);
end

if nparms == 2
    [r1 c1] = size(arg1);
    [r2 c2] = size(arg2);
end 

if sizeinfo == 0        
    if nparms == 1
        [rows columns] = size(arg1);
    end
    
    if nparms == 2
        scalararg1 = (prod(size(arg1)) == 1);
        scalararg2 = (prod(size(arg2)) == 1);
        if ~scalararg1 & ~scalararg2
            if r1 ~= r2 | c1 ~= c2
                errorcode = 1;
                return;         
            end
        end
        if ~scalararg1
            [rows columns] = size(arg1);
        elseif ~scalararg2
            [rows columns] = size(arg2);
        else
            [rows columns] = size(arg1);
        end
    end
    
    if nparms == 3
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
    end 
end

if sizeinfo == 1
    scalararg1 = (prod(size(arg1)) == 1);
    if nparms == 1
        if prod(size(arg2)) ~= 2
            errorcode = 2;
            return;
        end
        if  ~scalararg1 & arg2 ~= size(arg1)
            errorcode = 3;
            return;
        end
        if (arg2(1) < 0 | arg2(2) < 0 | arg2(1) ~= round(arg2(1)) | arg2(2) ~= round(arg2(2))),
            errorcode = 4;
            return;
        end 
        rows    = arg2(1);
        columns = arg2(2);
    end
    
    if nparms == 2
        if prod(size(arg3)) ~= 2
            errorcode = 2;
            return;
        end
        scalararg2 = (prod(size(arg2)) == 1);
        if ~scalararg1 & ~scalararg2
            if r1 ~= r2 | c1 ~= c2
                errorcode = 1;
                return;         
            end
        end
        if (arg3(1) < 0 | arg3(2) < 0 | arg3(1) ~= round(arg3(1)) | arg3(2) ~= round(arg3(2))),
            errorcode = 4;
            return;
        end 
        if ~scalararg1
            if any(arg3 ~= size(arg1))
                errorcode = 3;
                return;
            end
            [rows columns] = size(arg1);
        elseif ~scalararg2
            if any(arg3 ~= size(arg2))
                errorcode = 3;
                return;
            end
            [rows columns] = size(arg2);
        else
            rows    = arg3(1);
            columns = arg3(2);
        end
    end
    
    if nparms == 3
        if prod(size(arg4)) ~= 2
            errorcode = 2;
            return;
        end
        scalararg1 = (prod(size(arg1)) == 1);
        scalararg2 = (prod(size(arg2)) == 1);
        scalararg3 = (prod(size(arg3)) == 1);

        if (arg4(1) < 0 | arg4(2) < 0 | arg4(1) ~= round(arg4(1)) | arg4(2) ~= round(arg4(2))),
            errorcode = 4;
            return;
        end 

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
            if any(arg4 ~= size(arg1))
                errorcode = 3;
                return;
            end
            [rows columns] = size(arg1);
        elseif ~scalararg2
            if any(arg4 ~= size(arg2))
                errorcode = 3;
                return;
            end
            [rows columns] = size(arg2);
        elseif ~scalararg3
            if any(arg4 ~= size(arg3))
                errorcode = 3;
                return;
            end
            [rows columns] = size(arg3);
        else
            rows    = arg4(1);
            columns = arg4(2);
        end
    end 
end

if sizeinfo == 2
    if nparms == 1
        scalararg1 = (prod(size(arg1)) == 1);
        if ~scalararg1
            [rows columns] = size(arg1);
            if rows ~= arg2 | columns ~= arg3 
                errorcode = 3;
                return;
            end
        end
    if (arg2 < 0 | arg3 < 0 | arg2 ~= round(arg2) | arg3 ~= round(arg3)),
        errorcode = 4;
        return;
    end 
        rows = arg2;
        columns = arg3;
    end
    
    if nparms == 2
        scalararg1 = (prod(size(arg1)) == 1);
        scalararg2 = (prod(size(arg2)) == 1);
        if ~scalararg1 & ~scalararg2
            if r1 ~= r2 | c1 ~= c2
                errorcode = 1;
                return;         
            end
        end
        if ~scalararg1
            [rows columns] = size(arg1);
            if rows ~= arg3 | columns ~= arg4 
                errorcode = 3;
                return;
            end     
        elseif ~scalararg2
            [rows columns] = size(arg2);
            if rows ~= arg3 | columns ~= arg4 
                errorcode = 3;
                return;
            end     
        else
            if (arg3 < 0 | arg4 < 0 | arg3 ~= round(arg3) | arg4 ~= round(arg4)),
                errorcode = 4;
                return;
            end 
            rows = arg3;
            columns = arg4;
        end
    end
    
    if nparms == 3
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
            if rows ~= arg4 | columns ~= arg5 
                errorcode = 3;
                return;
            end     
        elseif ~scalararg2
            [rows columns] = size(arg2);
            if rows ~= arg4 | columns ~= arg5 
                errorcode = 3;
                return;
            end
        elseif ~scalararg3
            [rows columns] = size(arg3);
            if rows ~= arg4 | columns ~= arg5 
                errorcode = 3;
                return;
            end     
        else
            if (arg4 < 0 | arg5 < 0 | arg4 ~= round(arg4) | arg5 ~= round(arg5)),
                errorcode = 4;
                return;
            end 
            rows    = arg4;
            columns = arg5;
        end
    end 
end




