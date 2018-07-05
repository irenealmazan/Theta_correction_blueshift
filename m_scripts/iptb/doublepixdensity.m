function out = doublepixdensity(in)

out = zeros(size(in)*2);

for i=1:size(in,2) %x direction
    
    indi = (i-1)*2 +1;
    
    for j =1:size(in,1) %y direction
        
        indj = (j-1)*2+1;
        
        out(indj:indj+1, indi:indi+1) = ones(2,2)*in(j,i);
    end
end