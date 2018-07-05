function F = fts(x)

n=length(x);
F=zeros(size(x));
omeg = exp(-2*pi*i/n);

for p = 0:n-1 %hit up each frequency component    
    for j = 0:n-1  %each frequency is the weighted sum of all x's
        
        F(p+1) = F(p+1) + omeg^(j*p)*x(j+1);
        
    end
end
        