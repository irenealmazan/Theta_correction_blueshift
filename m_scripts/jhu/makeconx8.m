function makeconx8(confilein, confileout)

%first read the config file
fid1=fopen(confilein, 'r');
rho = str2num(fgetl(fid1));
atomtypes = str2num(fgetl(fid1));
%numdifat=atomtypes;

for i=1:atomtypes
    row = str2num(fgetl(fid1));
    numdifat(i) = row(1);
    radii(i)= row(2);
end

fid2=fopen(confileout, 'w');
fprintf(fid2, '%f\n',rho);
fprintf(fid2, '%i\n',atomtypes);
for i=1:atomtypes
    fprintf(fid2, '%i  %f\n', numdifat(i)*8, radii(i));
end

for i=1:sum(numdifat)
    %get the original position
    pos = fscanf(fid1, '%f %f %f', [3 1]);
    
    %rewrite the original position to 8 new coordinates in the box
    pos = pos/2+.5; %the +x +y +z quadrant
    
    %+z
    fprintf(fid2, '%f %f %f\n', [pos(1) pos(2) pos(3)]);
    fprintf(fid2, '%f %f %f\n', [-pos(1) pos(2) pos(3)]);
    fprintf(fid2, '%f %f %f\n', [pos(1) -pos(2) pos(3)]);
    fprintf(fid2, '%f %f %f\n', [-pos(1) -pos(2) pos(3)]);

    %-z
    fprintf(fid2, '%f %f %f\n', [pos(1) pos(2) -pos(3)]);
    fprintf(fid2, '%f %f %f\n', [-pos(1) pos(2) -pos(3)]);
    fprintf(fid2, '%f %f %f\n', [pos(1) -pos(2) -pos(3)]);
    fprintf(fid2, '%f %f %f\n', [-pos(1) -pos(2) -pos(3)]);
    
end

%xyz = fscanf(fid1, '%f %f %f', [3 inf]);
%xyz = xyz';

fclose(fid1);
fclose(fid2);