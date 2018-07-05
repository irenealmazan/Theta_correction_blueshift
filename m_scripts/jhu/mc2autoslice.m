function mc2autoslice(mcfilename, atslfilename, atomicnums);

%This function takes standard rmc configuration files and turns them into
%configuration files that can be used with Kirkland's TEM simulation
%packages "autoslice", from Advanced Computing in Electron
%Microscopy, Earl J Kirkland 1998 Plenum Press
%
%Input
%mcfilename - string containing the rmc file name
%atslfilename - name to be given to autoslice-compatible file
%atomicnums - atomic numbers in atomic order, eg: [Z1 Z2 Z3...]

%first read the config file
fid1=fopen(mcfilename, 'r');
rho = str2num(fgetl(fid1));
atomtypes = str2num(fgetl(fid1));
%numdifat=atomtypes;

for i=1:atomtypes
    row = str2num(fgetl(fid1));
    numdifat(i) = row(1);
    radii(i)= row(2);
end

xyz = fscanf(fid1, '%f %f %f', [3 inf]);
xyz = xyz';

xyzrotate = xyz;
xyzrotate(:,1) = xyz(:,2);
xyzrotate(:,2) = -xyz(:,1);

fclose(fid1);

halfedge = .5*(sum(numdifat)/rho)^(1/3);

fid2 = fopen(atslfilename, 'w');

fprintf(fid2, ['converted from monte carlo file: ' mcfilename '\n']);
fprintf(fid2, [num2str(halfedge*2) '  ' num2str(halfedge*2) ...
    '  ' num2str(halfedge*2) '\n']);

counter=1;
for i=1:atomtypes
    for j=1:numdifat(i)
        fprintf(fid2, '%i  %f  %f  %f  1.0  0\n', atomicnums(i), ... 
            (xyzrotate(counter,:)+1)*halfedge);
        counter = counter +1;
    end
end

fprintf(fid2, '-1\n\n\n');
fclose(fid2);