function [weight partialweight dsf] = fuossweight(gammas, c, qm, ffactorsE1, ffactorsE2, fpE1, fpE2, qiq, anomalel)

prefactor = c(anomalel) * (fpE1(anomalel) - fpE2(anomalel));
if (prefactor < 0)
    prefactor = -prefactor;
end

denom = zeros(1,length(ffactorsE1));
denom = denom';

for(i=1:length(c))
    denom = denom + (gammas(i) *(ffactorsE1(:,i) + ffactorsE2(:,i)));
end

weight = prefactor * denom;
    
for(i=1:length(c))
    partialweight(:,i) = c(anomalel) * (ffactorsE1(:,anomalel).*conj(ffactorsE1(:,i)) - ffactorsE2(:,anomalel).*conj(ffactorsE2(:,i)));
    if (i~=anomalel) 
        partialweight(:,i) = partialweight(:,i)*2;
    end
    partialweight(:,i) = partialweight(:,i) ./ weight;
    partialweight2(:,i) = spline(qm, partialweight(:,i), qiq(:,1));
end

dsf= partialweight2(:,1).*qiq(:,3) + partialweight2(:,2).*qiq(:,4) + partialweight2(:,3).*qiq(:,6);
dsf = dsf./qiq(:,1);[