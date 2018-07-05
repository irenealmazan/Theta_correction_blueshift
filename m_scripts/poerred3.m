function [retrphase newobj] = poerred3(dp, support, num, savenum, prevobj)

dp = abs(dp);

A=zeros(size(dp));
B=zeros(size(dp));
C=zeros(size(dp));
D=zeros(size(dp));

ssqexpamp = sum(sum(sum(abs(dp.^2))));

if isempty( findobj('Name', 'phasing'))
    hfig = gca;
else
    hfig = findobj('Name', 'phasing');
end

if isempty( findobj('Name', 'chi'))
    hchi = gca+1;
else
    hchi = findobj('Name', 'chi');
end

figure(hchi);clf;

%initial guess as to the phase of the DP comes from a FT of the support
if isempty(prevobj) A = fftn(support); end

if ~isempty(prevobj) A = fftn(prevobj.object); end
%imagesc(fftshift(log(abs((A)))));pause%
%imagesc(fftshift(angle((A))));pause

h=waitbarpos(0, 'error reduction phase iterations');

chi =[];
chilen = 0;
if ~isempty(prevobj) 
    chi = prevobj.chi; 
    chilen = length(chi);
end

outsup = find(abs(support) == 0);
insup = find(abs(support) ~= 0);

ii= sqrt(-1);

for i=1:num
    
    %we apply the phase guess to the measured dp amplitudes
    B = abs(dp) .* exp(ii * angle(A));
    
    %length(find(angle(A) ~= angle(B)))
    
    %we bring the new guess to real space via IFT
    C = ifftn(B);

    %we knock off all intensity in real space that falls outside the bounds
    %of the support
    
    D=C;
    temp = find(real(D)<0);
    display(['number of real values < 0 ' num2str(length(temp))]);
    temp = intersect(outsup, temp);
    display(['intersect temp and outsup' num2str(length(temp))]);
    D(temp) = abs(D(temp)) * exp(ii * -pi);
    temp = find(real(D)>0);
    display(['number of real values > 0 ' num2str(length(temp))]);
    temp = intersect(outsup, temp);
    display(['intersect temp and outsup' num2str(length(temp))]);
    D(temp) = abs(D(temp));

    
    if(mod(i,savenum)==0)
        displayphasing(A, B, C, D, support, 3, hfig, hchi); 
    end
    
    %return to dp space and apply the 
    A = fftn(D);
    
    %find a chi fit.
    chi(chilen+i,1) = sum(sum(sum( (abs(A)-abs(dp)).^2))) / ssqexpamp ;
    plot([1:length(chi)], log10(chi));
    title(['chi fit is ' num2str(chi(i,1))]);
    ylabel('log chi');
    
    waitbar(i/num);
end

close(h);

%imagesc(angle(A));

retrphase = abs(dp).*exp(sqrt(-1)*angle(A));

%find the minimum number of points to describe the support
[ind1 ind2 ind3] = ind2sub(size(support), find(support==1));
suphull = [ind2 ind1 ind3];
K = convhulln(suphull);
suphull = suphull(unique(K),:);

if isempty(prevobj)
    histnum=1;
    newobj.history(histnum).alg = 'erred';
    newobj.history(histnum).sup = inputname(2);
    newobj.history(histnum).suphull = suphull;
    newobj.history(histnum).dp = inputname(1);
    newobj.history(histnum).iter = i;
    newobj.chi = chi;
end

if ~isempty(prevobj)
    
    newobj = prevobj;
    histnum = length(prevobj.history);
    histnum= histnum+1;
    newobj.history(histnum).alg = 'erred';
    newobj.history(histnum).sup = inputname(2);
    newobj.history(histnum).suphull = suphull;
    newobj.history(histnum).dp = inputname(1);
    newobj.history(histnum).iter = i;
    newobj.chi = chi; 
end

newobj.object = D;
newobj.dp = fftshift(A);
totiter = 0;
for i=1:histnum totiter = totiter + newobj.history(i).iter; end
newobj.totiter = totiter;

