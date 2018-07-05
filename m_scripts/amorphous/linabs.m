function [mu,murho] = linabs(energy,s,rho);% LINABS Linear absorption coefficient%% Usage: [mu,murho] = linabs(energy,s,rho)%% Determines the linear absorption coefficient and% mass absorption coefficient for a sample of% known composition and density.%% Input:%   energy = x-ray energy (eV)%   s      = formula unit string, e.g. '2H1O' for water%   rho    = density (gm/cm^3)%% Output:%   mu     = linear absorption coefficient (1/cm)%   murho  = mass absorption coefficient (cm^2/g)%% Todd Hufnagel (hufnagel@jhu.edu) 4-30-00% Requires: AVGATOMmurho=avgatom(energy,s);mu=murho.*rho;