function C=correlogram(linearpix)

linearpix=linearpix(:);

numpix= length(linearpix);
dphi = 360/numpix;

linearpix = vertcat(linearpix, linearpix);

C=zeros(numpix,1);

for smallphi=0:numpix %goes through small phi
    for bigphi=1:numpix %goes through big phi
        term1 = linearpix(bigphi);
        term2 = linearpix(bigphi+smallphi);
        term12= term1*term2;
        C(bigphi) = C(bigphi)+term12*dphi;
    end
end