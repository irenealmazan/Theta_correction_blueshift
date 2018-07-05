function [dsf,wt] = dsfstanford(k,inte1,inte2,s,eoi,e1,e2,ff,fp1,fp2,fpp1,fpp2,gamma,adasf1,adasf2)
% DSF    Differential structure factor (weighted method)
%
% Usage: [dsf,wt] = dsf(k,inte1,inte2,s,eoi,e1,e2,ff,fp1,fp2,fpp1,fpp2{,gamma,adasf1,adasf2})
%
% Calculates the differential structure factor, defined
% as dsf = (inte1-inte2-(coh1-coh2))./(weighting factor), for use
% the differential anomalous scattering technique. The weighting factor
% approximately removes any angular dependence and normalizes the DSF 
% to a per atom basis. Note that we assume that differences in f'' 
% between the two energies are small. Like the total structure factor, 
% the differential structure factor (DSF) -> -1 as k -> 0 and DSF -> 0 
% as k -> infinity. Hope believes that gamma can be changed accordingly 
% to give the best satisfaction of these boundary conditions.
%
% Input:
%       k   =   vector containing scattering vector magnitudes
%       inte1, inte2 =  corrected, normalized experimental
%               elastic intensities at two photon energies (see WAXS.M)
%       s       = sample composition string, e.g. '2H1O'
%       eoi     = symbol or atomic number of element of interest
%       ff       = atomic form factor f(k) for elastic x-ray scattering. 
%           If ff=[], a parameterization will be used; which one depends 
%           on the global FF_FLAG - see the documentation. If specified, 
%           ff is a two-dimensional array with one column for each  element
%           (in the same order as in  s above), with the rows corresponding
%           to the values of k.
%       e1, e2  = the two x-ray energies (eV)
%       fp1,fp2 = anomalous scattering factors, f', at the two energies.
%           If you want to specify the f' explicitly, then make fp a row
%           vector, with one entry for each element (in the same order,
%           as in s above. If empty, a parameterization will be used
%           to calculate the f' values. See CALSCAT.M for details.
%       fpp1,fpp2 = anomalous scattering factors f'' (see above)
%       gamma   = estimated fraction of like near-neighbors for
%                 the nearest-neighbor shell of the element of interest.
%                 Set gamma=[] to force DSF to use the difference of the
%                 average scattering factors for weighting, same as
%                 the old version of DSF.M
%       adasf1,adasf2 = angular dependence of the anomalous scattering factors
%                 (optional; see CALSCAT.M for description)  
% Output:
%       dsf     = the differential structure factor around atoms of type A
%       wt      = the weighting function
%

% Originally by Hope Ishii (1/26/98 and 9/27/00)
% 29-Jun-05 Todd Hufnagel Completely rewritten; Hope's version was
% hard-wired for Mo-Ge. This version is much more general and flexible.
% Note well, though, that I have not tested it thoroughly, so use with
% caution.

% References:   Hope Ishii, Ph.D. dissertation, Stanford (2002)
%               P. H. Fuoss et.al., PRL 46, 1537 (1981)
%               P. H. Fuoss, Ph.D. dissertation, Stanford (1980)

global ANOMAL_FLAG FF_FLAG

if isempty(FF_FLAG)
	FF_FLAG=0; % Use parameterization of Waasmaier & Kirfel for form factor by default
end

if isempty(ANOMAL_FLAG)
	ANOMAL_FLAG=0; % Use parameterization of Henke for anomalous scattering factors by default
end

% if not specified, the ASF are assumed to be independent of angle
if nargin<15
	adasf2=[];
end
if nargin<14
    adasf1=[];
end

% If not given, assign weighting function to make the definition of DSF
% used here match the old DSF.M (from Ritvaa Serimaa)
if nargin<13
    gamma=[];
end

% Do a little error and reality checking
if ~isempty(gamma) & ((gamma<0)|(gamma>1))
	error('Gamma out of range - must be between 0 and 1 (or empty)')
end

% figure out our composition
[z,w] = strtoz(s);
nz = length(z);
w = w./(sum(w));    % atomic fractions

% make sure the inputs are column vector
k=k(:);
inte1=inte1(:);
inte2=inte2(:);

% determine which edge we're working on, and do a little error checking
eoi = element(eoi);
j = find(z==eoi);
if isempty(j)
    error('Element of interest is not in composition string, s!')
end

% Now we get down to business.
% If f', f'' are not given, look them up
if isempty(fp1)
    if isempty(fpp1)
        [fp1,fpp1] = fprime(e1,s);
    else
        fp1 = fprime(e2,s);
    end
end
if isempty(fp2)
    if isempty(fpp2)
        [fp2,fpp2] = fprime(e2,s);
    else
        fp2 = fprime(e2,s);
    end
end
if isempty(fpp1)
    [junk,fpp1]=fprime(e1,s);
end
if isempty(fpp2)
    [junk,fpp2]=fprime(e2,s);
end

coh1 = zeros(size(k));
coh2 = zeros(size(k));
ave1 = zeros(size(k));
ave2 = zeros(size(k));

for l = 1:nz,       % For each element
    % Calculate the complex scattering factor
    % Do angular correction on fp and fpp, as necessary.
    [fpo1,fppo1]= localanglecorr(k,e1,fp1(l), fpp1(l), z(l), adasf1);
    [fpo2,fppo2]= localanglecorr(k,e2,fp2(l), fpp2(l), z(l), adasf2);

    if isempty(ff)  % use parameterization for form factor
        f1(:,l) = formfact(z(l),k) + fpo1 + i.*fppo1;
        f2(:,l) = formfact(z(l),k) + fpo2 + i.*fppo2; 
    else            % user-supplied form factor
        f1(:,l) = ff(:,l) + fp1 + i.*fpp1;	
        f2(:,l) = ff(:,l) + fp2 + i.*fpp2;
    end
    
    coh1   = coh1 + w(l).*f1(:,l).*conj(f1(:,l));	    		% Calculate the coherent scattering
    coh2   = coh2 + w(l).*f2(:,l).*conj(f2(:,l));
    ave1   = ave1 + w(l).*f1(:,l);
    ave2   = ave2 + w(l).*f2(:,l); 
end;

% Ritvaa's version of DSF used the difference in the average scattering
% factor for weighting. The Stanford crew use something a bit more 
% complex (see references above). Here, we let the user pick.
if ~isempty(gamma) % use the Stanford weighting scheme
    disp('Using new (Stanford) weighting scheme')
    % Calculate the prefactor on the weighting funcntion
    prefactor=w(j).*real(fp2(j)+fpp2(j)-fp1(j)-fpp1(j));

    % Estimate the fraction of atoms of each type in the 1st shell around A
    g=zeros(size(z));
    for l=1:nz
        if l==j % this is the element of interest
            g(l) = gamma;
        else    % this is one of the other elements
            g(l) = (1-gamma).*(1-w(j));
    end
end

% Now calculate the weighting function. Note that the following are array 
% calculations, performed for the whole k range at once.
    wt=zeros(size(k));
    for l = 1:nz,			                		% For each element
 		    % Do angular correction on fp and fpp, as necessary.
            [fpo1,fppo1]= localanglecorr(k,e1,fp1(l), fpp1(l), z(l), adasf1);
            [fpo2,fppo2]= localanglecorr(k,e2,fp2(l), fpp2(l), z(l), adasf2);
		    
            if isempty(ff)
			    f1 = formfact(z(l),k) + fpo1 + i*fppo1;% use either C&M or W&K for form factor
			    f2 = formfact(z(l),k) + fpo2 + i*fppo2;% use either C&M or W&K for form factor
		    else					
			    f1= ff(:,l) + fpo1 + i*fppo1;		% User-supplied form factor
			    f2= ff(:,l) + fpo2 + i*fppo2;		% User-supplied form factor
            end
          	wt = wt + g(l) .* real(f1 + f2); % Calculate the weighting
    end
wt=prefactor.*wt(:); 

else    % use Ritvaa's weighting scheme
    disp('Using old-style (avef) weighting scheme.')
    wt = abs(ave1.^2) - abs(ave2.^2);
end

dsf = ((inte2-coh2)-(inte1-coh1))./wt;
if isempty(gamma)   % for Ritvaa's scheme, make sure dsf comes out positive
    dsf=-dsf;
end

%-----------------------------------------
function [fpo,fppo]= localanglecorr(k,energy,fpin,fppin,zed,adasf)
% [fpo,fppo]= localanglecorr(fpin,fppin,zed,adasf)
% Note: This function is also incorporated into DSF.M, so if
% any changes are made that affect the angle-dependent ASF input/output
% parameters of CALSCAT, they should be made in both places.
if isempty(adasf) | adasf(1)~=zed,
	fpo = fpin;
	fppo= fppin;
else
	ksq= sind(ktotheta(k,etolambda(energy))).^2;
	fpo = fpin  *(1 + adasf(2)*ksq + adasf(3)*ksq.^2);
	fppo= fppin *(1 + adasf(4)*ksq + adasf(5)*ksq.^2);
end
