function [retrphase newobj] = pchio3(dp, support, num, savenum, prevobj, beta, phirange)

dp = abs(dp);

%initial guess as to the phase of the DP comes from a FT of the support
if isempty(prevobj) A = fftn(support); end

if ~isempty(prevobj) A = fftn(prevobj.object); end
%imagesc(fftshift(log(abs((A)))));pause%
%imagesc(fftshift(angle((A))));pause

h=waitbarpos(0, 'hio phase iterations');

if isempty( findobj('Name', 'phasing'))
    hfig = figure;
    set(gcf, 'Name','phasing');
else
    hfig = findobj('Name', 'phasing');
end

if isempty( findobj('Name', 'chi'))
    hchi = figure;
    set(gcf, 'Name','chi');
else
    hchi = findobj('Name', 'chi');
end

figure(hchi);clf;

ssqexpamp = sum(sum(sum(abs(dp))));

if isempty(prevobj) Dold= abs(support); end
if ~isempty(prevobj) Dold = prevobj.object; end

outsup = find(abs(support) == 0);
insup = find(abs(support) ~= 0);

chi=[];
chilen=0;
if ~isempty(prevobj) 
    chi = prevobj.chi; 
    chilen = length(chi);
end

ii=sqrt(-1);

for i=1:num
    
    
    %we apply the phase guess to the measured dp amplitudes
    B = abs(dp) .* exp(ii * angle(A));
    
    %we bring the new guess to real space via IFT
    C = ifftn(B);
    
    if(i>1) Dold = D; end
    
    D=C;
    D(outsup) = Dold(outsup) - beta * C(outsup);
    Din = D(insup);
    indexDinbad = find(angle(Din)<phirange(1) | angle(Din)>phirange(2));
    D(indexDinbad) = Dold(indexDinbad) - beta * C(indexDinbad);
    
    if(mod(i,savenum)==0)        
        displayphasing(A, B, C, D, support, 3, hfig, hchi);
    end
    
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
    newobj.history(histnum).alg = 'hio';
    newobj.history(histnum).sup = inputname(2);
    newobj.history(histnum).suphull = suphull;
    newobj.history(histnum).dp = inputname(1);
    newobj.history(histnum).beta = beta;
    newobj.history(histnum).iter = i;
    newobj.chi = chi;
end

if ~isempty(prevobj)
    
    newobj = prevobj;
    histnum = length(prevobj.history);
    histnum= histnum+1;
    newobj.history(histnum).alg = 'hio';
    newobj.history(histnum).sup = inputname(2);
    newobj.history(histnum).suphull = suphull;
    newobj.history(histnum).dp = inputname(1);
    newobj.history(histnum).beta = beta;
    newobj.history(histnum).iter = i;
    newobj.chi = chi; 
end

newobj.object = D;
newobj.dp = fftshift(A);
totiter = 0;
for i=1:histnum totiter = totiter + newobj.history(i).iter; end
newobj.totiter = totiter;
