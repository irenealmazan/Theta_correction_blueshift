function cohstats(imgs)
% count single and double hits
global fina;
d=size(imgs);
if(numel(d)==2)
    d=[d 1];
end

%for raw images, not photonized
%{
w=find(imgs(:)>300);
nophot2=size(w,1);  %counts all pixels w at least one photon in them
nophot1=sum(imgs(w))/780.; %counts all photons that hit the detector
nopixels=d(1)*d(2)*d(3);
prob1=nophot1/nopixels;
prob2=nophot2/nopixels;
sing=size(find((imgs(:)>300)&(imgs(:)<=850))) ;
doub=size(find((imgs(:)>850)&(imgs(:)<=1600.)));
%}

%for photonized images (pixels in image represent integer photon events)
%%{
nophot1 = length(find(imgs>0)); %counts all pixels w at least one photon in them
nophot2 = sum(sum(sum(imgs)));  %counts all photons that hit the detector
nopixels=numel(imgs);

%take two different estimates of the poisson expectation value lambda
prob1=nophot1/nopixels; %fraction of pixels that detected at least one photon
prob2=nophot2/nopixels; %total counted photons / total num pixels
pois1=(prob1)^2*nopixels/2;
pois2=(prob2)^2*nopixels/2;

sing= length(find(imgs==1));
doub= length(find(imgs==2));
p1=sing(1)/nopixels;  %probability of single hits
p2 = doub(1)/nopixels;%probability of a double hit

M=1/(2*p2/(p1^2)-1.);
%}

%playing around with different stats 
%{
nopixels=numel(imgs);
sing= length(find(imgs==1));
doub= length(find(imgs==2));
tot = sing+doub;

avephotppix = tot/nopixels;  %use this as expectation value for poisson
%}



fprintf('%s:\n\n  number singles %.0f doubles %.0f  Poisson  %0.f,%.0f\n', ...
    fina, sing, doub, pois1, pois2);
fprintf('  probs %.3g, %.3g. nos %.0f %.0f, No. of pixels %d\n', ...
    prob1, prob2, nophot1,nophot2,nopixels); 
fprintf('  Possion(n=2) %.3g,%.3g Poisson(n=3) %.3g %.3g\n', ...
    poiss(2,prob1),poiss(2,prob2),poiss(3,prob1),poiss(3,prob2));
fprintf('M = %.3f P1 = %.3g  P2 = %.3g\n\n',M,p1,p2);

