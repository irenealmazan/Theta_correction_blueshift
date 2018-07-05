function sw = calcspherewave(ksi,eta,d0,A0,lambda)

sw = sqrt(ksi.^2 + eta.^2 + d0^2); %distance from origin to an arbitrary point on the xy plane
sw = sw-d0; %this is the path length difference traveled between the two rays
sw = A0 .* exp (i*2*pi*sw/lambda); %this is the exact formulation
