function x = ifts(F)

n=length(F);
x=zeros(size(F));
omeg = exp(2*pi*i/n);

for j = 0:n-1 %hit up each real space value    
    for p = 0:n-1  %sum up frequency matrix
        
        x(j+1) = x(j+1) + omeg^(j*p)*F(p+1);
        
    end
end

x=x/n;