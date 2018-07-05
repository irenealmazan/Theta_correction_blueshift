function support = make3dsupport(dp, bounds, offset)

support = zeros(size(dp));
center = round(size(dp)/2)
center = center+offset;

bounds = round(bounds/2)

for i= center(1)-bounds(1): center(1)+bounds(1)
    for j = center(2)-bounds(2) : center(2)+bounds(2)
        for k = center(3)-bounds(3) : center(3) + bounds(3)
            
            support(i,j,k) = 1;
            
        end
    end
end
