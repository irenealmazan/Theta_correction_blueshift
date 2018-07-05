function [retrphase newobj] = hio(dp, support, num, savenum, prevobj, beta, phirange)

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

dp = abs(dp);
ssqexpamp = sum(sum(sum(abs(dp.^2))));

imagesc(real(support));pause(1)

%initial guess as to the phase of the DP comes from a FT of the support
if isempty(prevobj) 
    %A = fft2(support); 
    A= rand(size(support)) * 2*pi -pi;
end

if nargin < 7 | isempty(phirange) phirange=[-pi pi]; end

if ~isempty(prevobj) A = fft2(prevobj.object); end
%imagesc(fftshift(log(abs((A)))));pause%
%imagesc(fftshift(angle((A))));pause

h=waitbarpos(0, 'phase iterations');

chi =[];
chilen = 0;
if ~isempty(prevobj) 
    chi = prevobj.chi; 
    chilen = length(chi);
end

if isempty(prevobj) Dold= abs(support); end
if ~isempty(prevobj) Dold = prevobj.object; end

outsup = find(abs(support) == 0);
insup = find(abs(support) ~= 0);

for i=1:num
    
    
    %we apply the phase guess to the measured dp amplitudes
    B = abs(dp) .* exp(sqrt(-1) * angle(A));
    
    %we bring the new guess to real space via IFT
    C = ifft2(B);

    %we knock off all intensity in real space that falls outside the bounds
    %of the support
    
    if(i>1) Dold = D; end
    
    D=C;
    %D(insup) = C(insup);
    D(outsup) = Dold(outsup) - beta * C(outsup);
    Din = D(insup);
    indexDinbad = find(angle(Din)<phirange(1) | angle(Din)>phirange(2));
    D(indexDinbad) = Dold(indexDinbad) - beta * C(indexDinbad);
        
    if(mod(i,savenum)==0)
         
        figure(hfig);
        
        subplot(4,2,1);
        imagesc(fftshift(log(abs(A)+1)));
        axis equal tight
        subplot(4,2,2);
        imagesc(fftshift(angle(A)));
        axis equal tight

        subplot(4,2,3);
        imagesc(fftshift(log(abs(B)+1)))
        axis equal tight
        subplot(4,2,4);
        imagesc(fftshift(angle(B)));
        axis equal tight

        subplot(4,2,5);
        imagesc(abs(C));
        axis equal tight; zoom(3);
        subplot(4,2,6);
        imagesc(angle(C));
        axis equal tight ; zoom(3);

        subplot(4,2,7);
        imagesc(abs(D));
        axis equal tight; zoom(3);
        subplot(4,2,8);
        imagesc(angle(D)); 
        axis equal tight;zoom(3);
        
        pause(.01)

    end
    
    %return to dp space and apply the 
    A = fft2(D);
    
    %find a chi fit.
    chi(chilen+i,1) = sum(sum(sum( (abs(A)-abs(dp)).^2))) / ssqexpamp ;
    figure(hchi);
    plot([1:length(chi)], log10(chi));
    title(['chi fit is ' num2str(chi(chilen+i,1)) ' HIO']);
    ylabel('log chi');

    waitbar(i/num,h);
    
end

close(h);

%imagesc(angle(A));

retrphase = abs(dp).*exp(sqrt(-1)*angle(A));

if isempty(prevobj)
    histnum=1;
    newobj.history(histnum).alg = 'hio';
    newobj.history(histnum).sup = inputname(2);
    newobj.history(histnum).suphull = support;
    newobj.history(histnum).dp = inputname(1);
    newobj.history(histnum).beta = beta;
    newobj.history(histnum).iter = i;
    newobj.chi = chi;
    newobj.history(histnum).sigX = NaN;
    newobj.history(histnum).sigY = NaN;
end

if ~isempty(prevobj)
    
    newobj = prevobj;
    histnum = length(prevobj.history);
    histnum= histnum+1;
    newobj.history(histnum).alg = 'hio';
    newobj.history(histnum).sup = inputname(2);
    newobj.history(histnum).suphull = support;
    newobj.history(histnum).dp = inputname(1);
    newobj.history(histnum).beta = beta;
    newobj.history(histnum).iter = i;
    newobj.chi = chi; 
    newobj.history(histnum).sigX = NaN;
    newobj.history(histnum).sigY = NaN;
end

newobj.object = D;
newobj.dp = fftshift(A);
totiter = 0;
for i=1:histnum totiter = totiter + newobj.history(i).iter; end
newobj.totiter = totiter;
newobj.sigX = NaN;
newobj.sigY = NaN;
