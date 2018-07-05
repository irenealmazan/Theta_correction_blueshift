function [fp,fpp,den]=xastoasf(energy,xas,atom,edgename,erange,nfit)
% XASTOASF Calculate f' and f'' from x-ray absorption data
%
% Usage: [fp,fpp,den]=xastoasf(energy,xas,atom{,edgename,erange,nfit})
%
% Calculates the anomalous scattering factors f' and f'' from x-ray
% absorption data in two steps. First, FPPFIT is called to calculate f''
% from the absorption data. Then, f' is calculated from f'' by a call to
% KRAMKRON. See FPPFIT and KRAMKROM for details on the calculations.
%
% Input:
%   energy   = vector of x-ray energies (eV)
%   xas      = vector of x-ray absorption data (ln(I0/I))
%   atom     = atomic number or name of element in question
%   edgename = absorption edge ('K' or 'LIII')
%              (optional; the default is 'K')
%   erange   = [emin emax] range of energies near edge to exclude
%              from fit to determine f'' (optional; default is from
%              50 eV below to 100 eV above the edge)
%   nfit     = order of background fit for f'' (optional; default
%              is 2, allowable values are 1-6.
%
% Output:
%  fp, fpp   = real and imaginary parts of the anomalous scattering
%              factor (electron units)
%  den       = number density * thickness for this element (1/Å^2)
%
% See also: FPPFIT, KRAMKRON

% 16-May-03 Todd Hufnagel (hufnagel@jhu.edu), Johns Hopkins University

% Requires FPPFIT, KRAMKRON

% Call FPPFIT and KRAMKRON to do the work; all necesary error
% checking is done there.
switch nargin
    case 3
        [fpp,den]=fppfit(energy,xas,atom);
    case 4
        [fpp,den]=fppfit(energy,xas,atom,edgename);
    case 5
        [fpp,den]=fppfit(energy,xas,atom,edgename,erange);
    case 6
        [fpp,den]=fppfit(energy,xas,atom,edgename,erange,nfit);
    otherwise
        error('Must provide at least the energy scale, absorption data, and name of element.')
end
fp=kramkron(energy,fpp,atom);