function [pos,atden,nelem,natoms,atrad]=readconfig(inputfile)
%
% READCONFIG - read data and atom positions from RMC configuration file
%
%  [pos,atden,nelem,natoms,atrad]=readconfig(inputfile)
%
% Input:
%   inputfile = configuration file name (string)
%
% Output:
%   pos = array (natoms x 3) of atomic positions (box units)
%   atden = atomic density (atoms/Å^3)
%   nelem = number of chemical elements
%   natoms = vector (nelem long) with number of atoms of each element
%   atdad = vector (nelem long) with radius of each atom type
%

% Todd Hufnagel (hufnagel@jhu.edu) 04-Jun-2007
% Stephan Hruszkewyca, entirely rewritten 6-13-08
fid1=fopen(inputfile, 'r');
atden = str2num(fgetl(fid1));
nelem = str2num(fgetl(fid1));

for i=1:nelem
    row = str2num(fgetl(fid1));
    natoms(i) = row(1);
    atrad(i)= row(2);
end

pos = fscanf(fid1, '%f %f %f', [3 inf]);
pos = pos';

fclose(fid1);
totatoms = sum(natoms);

% Summary report
fprintf('\r')
fprintf('%d atoms read from file %s\r',totatoms,inputfile)
fprintf('\r')
fprintf('Atomic density (atoms/≈^3): %6.4f\r',atden)
fprintf('\r')
fprintf('Element   Atoms   Fraction   Radius (≈)\r')
fprintf('-------   -----   --------   ----------\r')
for i=1:nelem
    fprintf('   %d      %5d     %4.2f        %4.2f\r',...
        i,natoms(i),natoms(i)./totatoms,atrad(i))
end
fprintf('\r') 
