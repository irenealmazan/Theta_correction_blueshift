%
% 2D Ising
% ising2(n,m,b,q)
%    n - size of domain.  Either an integer, in which domain is n x n
%                         Or a two-vector [n1 n2]
%    m - number of iterations, positive integer
%    b - coupling parameter / temperature
%                         b = 0 -> Very high temp / low coupling
%                         b > 1 -> Very low temp / high coupling
%                         b = 0.440687  critical point
%    q - number of independent runs to concatenate (default q=1)
%

%
% P. Fieguth
% University of Waterloo
% Oct, 2002
%

%from: http://ocho.uwaterloo.ca/~pfieguth/Software/Ising/ising2.m

function s = ising2(n,m,b,q)

n=round(n); m=round(m);
if (nargin < 4), q=1; else, q=round(q); end;
p=n; if (length(n)>1), p=n(1); n=n(2); end;

s = zeros(p*q,n);

t = zeros(p+2,n+2); for i=2:(p+1), for j=(2+rem(i,2)):2:(n+1), t(i,j)=1; end; end; e = find(t);
t = zeros(p+2,n+2); for i=2:(p+1), for j=(3-rem(i,2)):2:(n+1), t(i,j)=1; end; end; o = find(t);
nr = length(e);
ex = exp(-b*[-4:4]);

for qi=1:q,
  t = sign(randn(p+2,n+2));

  for j=1:m,
    el = t(e-1) + t(e+1) + t(e-p-2) + t(e+p+2);
    t(e) = 1;
    eel = ex(5+el);
    t(e(find(rand(1,nr)<(eel./(eel+1./eel))))) = -1;
    t(1,:)=t(p+1,:); t(p+2,:)=t(2,:); t(:,1)=t(:,n+1); t(:,n+2)=t(:,2);
    
    el = t(o-1) + t(o+1) + t(o-p-2) + t(o+p+2);
    t(o) = 1;
    eel = ex(5+el);
    t(o(find(rand(1,nr)<(eel./(eel+1./eel))))) = -1;
    t(1,:)=t(p+1,:); t(p+2,:)=t(2,:); t(:,1)=t(:,n+1); t(:,n+2)=t(:,2);
  end;
  s((1:p)+(qi-1)*p,:) = t(2:(p+1),2:(n+1));
end;