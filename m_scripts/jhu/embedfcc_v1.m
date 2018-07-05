function embedfcc(config, fracfcc, numxtal, radius, filename_xtals, ...
        filename_newcon, fwhm, quadfrac)

%This function takes a standard rmc configuration file, scoops out
%spherical sections from that configuration and replaces them with a
%randomly oriented fcc crystal with a prescribed lattice parameter.  The
%atoms on the lattice sites are randomly assigned from the list of atoms
%originally specified in the con file.  Output files are generated that
%contain the new configuration with rearranged atoms (filename_newcon), and
%a file is made that contains only the positions of atoms that were
%altered.
%
%Input:
%config - con file with original positions to be altered
%fracfcc - defines fcc lattice parameter as (average atomic radius)*fracfcc
%numxtals - the number of random centers onto which an fcc grain will be
%           put
%radius - the radius of the fcc spherical grains, in angstroms
%filename_xtals - file containing only the atoms in the grains
%filename_newcon - file containg the new configuration, with adjusted
%           density and atom numbers
%
%SH Feb 08

%first read the config file
fid1=fopen(config, 'r');
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

fclose(fid1);

xyz_original = xyz;

halfedge = .5 * (sum(numdifat)/rho)^(1/3);
radius = radius / halfedge;
fwhm = fwhm / halfedge;

averad = (sum(radii)/3) * fracfcc / halfedge;
avelattice = 4*averad/sqrt(2);

lowindex(1) = 0;
highindex(1) = numdifat(1);
for i=2:length(numdifat)
    lowindex(i) = lowindex(i-1) + numdifat(i-1);
    highindex(i) = highindex(i-1) + numdifat(i);
end

numuc = round(2* radius/avelattice);
fcc = genfcc(numuc);
fcc = fcc * avelattice ;
fcc = fcc - radius; %center it around 000
fcc = fccsphere(fcc, radius, quadfrac);

display(['crystals containing about ' num2str(length(fcc)) ' atoms will be inserted.']);
display(['they will be about ' num2str(2*halfedge*radius/10) ' nm in diameter.']);

%imbed the fcc in a random spot

for i=1:numxtal
    %display(['crystal ' num2str(i)]);
    
    numuc = round(2* radius/avelattice);
    fcc = genfcc(numuc);
    fcc = fcc * avelattice ;
    fcc = fcc - radius; %center it around 000
    fcc = fccsphere(fcc, radius, quadfrac);
    if(fwhm~=0) fcc = fccgauss(fcc, fwhm); end
    
    x=rand*2-1;
    y=rand*2-1;
    z=rand*2-1;
    
    if x<1-radius & x>-1+radius %we are fully inside the box
        xok = find(xyz(:,1)>x-radius & xyz(:,1)<x+radius);
    elseif x> 1-radius %too far to the positive, need some negative x's to fill in
        diff = radius - (1-x);
        xok = find((xyz(:,1)>x-radius & xyz(:,1)<1) | (xyz(:,1)>-1 & xyz(:,1)< -1+diff));
    elseif x< -1+radius % too far into the negative, need some positives to fill in
        diff = radius -(x+1);
        xok = find((xyz(:,1)> -1 & xyz(:,1)< x+radius) |( xyz(:,1)<1 & xyz(:,1) >1-diff));
    end
    
    if y<1-radius & y>-1+radius %we are fully inside the box
        yok = find(xyz(:,2)>y-radius & xyz(:,2)<y+radius);
    elseif y> 1-radius %too far to the positive, need some negative x's to fill in
        diff = radius - (1-y);
        yok = find((xyz(:,2)>y-radius & xyz(:,2)<1) | (xyz(:,2)>-1 & xyz(:,2)< -1+diff));
    elseif y< -1+radius % too far into the negative, need some positives to fill in
        diff = radius -(y+1);
        yok = find((xyz(:,2)> -1 & xyz(:,2)< y+radius) |( xyz(:,2)<1 & xyz(:,2) >1-diff));
    end

    if z<1-radius & z>-1+radius %we are fully inside the box
        zok = find(xyz(:,3)>z-radius & xyz(:,3)<z+radius);
    elseif z> 1-radius %too far to the positive, need some negative x's to fill in
        diff = radius - (1-z);
        zok = find((xyz(:,3)>z-radius & xyz(:,3)<1) | (xyz(:,3)>-1 & xyz(:,3)< -1+diff));
    elseif z< -1+radius % too far into the negative, need some positives to fill in
        diff = radius -(z+1);
        zok = find((xyz(:,3)> -1 & xyz(:,3)< z+radius) |( xyz(:,3)<1 & xyz(:,3) >1-diff));
    end
        
    ok = intersect(xok, yok);
    ok = intersect(ok,  zok);
    
    %plot3(xyz(ok,1),xyz(ok,2),xyz(ok,3), '.'); axis equal; pause
    
    counter=1;
    
    for j=1:length(ok)
        delx = abs(xyz(ok(j),1) - x);
        if delx >1
            delx = delx-2;
        end 
        
        dely = abs(xyz(ok(j),2) -y);
        if dely >1
            dely = dely-2;
        end
        
        delz = abs(xyz(ok(j),3) -z);
        if delz >1
            delz = delz-2;
        end
            
        dist = sqrt(delx^2 + dely^2 + delz^2);
        if dist <= radius
            oksphere(counter,:) = xyz(ok(j),:);
            oksphereindex(counter) = ok(j);
            %display(num2str(oksphere(counter,:)));
            counter = counter +1;
        end
    end
    
    %plot3(xyz(ok,1),xyz(ok,2),xyz(ok,3), '.',oksphere(:,1),oksphere(:,2),oksphere(:,3),'ro');axis equal; pause;
    
    if(~isempty(oksphereindex))
        order = randperm(length(oksphereindex));
        oksphereindex_o=oksphereindex;

        for i=1:length(oksphereindex)
            oksphereindex(i) = oksphereindex_o(order(i));
        end
    end
    
    fcc(:,1) = fcc(:,1)+x;
    fcc(:,2) = fcc(:,2)+y;
    fcc(:,3) = fcc(:,3)+z;
    
    for i=1:length(fcc)
        if fcc(i,1)<=-1
            fcc(i,1)=fcc(i,1)+2;
        end
        if fcc(i,1)>1
            fcc(i,1)=fcc(i,1)-2;
        end
        
        if fcc(i,2) <= -1
            fcc(i,2) = fcc(i,2)+2;
        end
        if fcc(i,2) > 1
            fcc(i,2) = fcc(i,2)-2;
        end
        
        if fcc(i,3) <= -1
            fcc(i,3) = fcc(i,3)+2;
        end
        if fcc(i,3) > 1
            fcc(i,3) = fcc(i,3)-2;
        end
    end
    
    %display(['size fcc: ' num2str(size(fcc)) ', size oksphere: ' num2str(size(oksphere))]);
    
    %oksphereindex
    %oksphere
    %fcc
    
    if size(fcc,1) > size(oksphere,1)
        difference = size(fcc,1) - size(oksphere,1);
        
        for j=1:size(oksphere,1) %substitute these for existing atom positions
            xyz(oksphereindex(j),:) = fcc(j,:);
        end
        
        numat = length(xyz);
        for j=1:difference
            %randindex = ceil(rand*(difference+size(oksphere,1))); %insert a value between this and the next xyz index
            randindex = ceil(rand*numat);
            %display([' adding an atom ' num2str(randindex) ' coord: ' num2str(fcc(j+size(oksphere,1),:))]);
            xyz = addrow(xyz, randindex, fcc(j+size(oksphere,1),:) );
            xyz_original = addrow(xyz_original, randindex, [-30 -30 -30]);
            
            oncethrough =0;
            for k=1:atomtypes
                if randindex <= highindex(k)
                    highindex(k) = highindex(k) +1;
                    if(~oncethrough) numdifat(k) = numdifat(k)+1; end
                    oncethrough =1;
                end
            end
            
        end
    end
    
    if size(fcc,1) < size(oksphere,1) %need to take values out of xyz
        difference = size(oksphere,1) - size(fcc,1);  
        numat = length(xyz);
        
        for j=1:size(fcc,1) %substitute these for existing atom positions
            xyz(oksphereindex(j),:) = fcc(j,:);
            %xyz(oksphereindex(j),:) = [0 0 0];
        end
        
        for j=1:difference
            %display([' removing an atom ' num2str(oksphereindex(j+size(fcc,1)) ) ]);
            xyz(oksphereindex(j+size(fcc,1)),:)=[-1.2 -1.2 -1.2];
            %display([' ' num2str(xyz(oksphereindex(j+size(fcc,1)),:))]);
            %xyz = removerow(xyz, oksphereindex(size(fcc,1)+j)); %removes this atom from the list
             
        end

    end

    %clear ok
    %clear fcc
    %clear oksphere
    %clear oksphereindex
    
    ok=[];
    fcc=[];
    oksphere=[];
    oksphereindex=[];
    
    %plot3(xyz(:,1),xyz(:,2),xyz(:,3),'r.');
    %plot3(oksphere(:,1),oksphere(:,2),oksphere(:,3),'g.');
    %plot3(fcc(:,1),fcc(:,2),fcc(:,3),'.');
 
    %axis equal
end

%hold off;

xtals = xyz(find(xyz(:,1)~= xyz_original(:,1)), :);
xtals = xtals(find(xtals(:,1)~= -1.2), :);  %given in case of more oksphere than fcc

byebye = find(xyz(:,1)<-1); %these are the positions which are being removed
byebye1 = byebye;
hello= find(xyz_original(:,1) <-1); 
display(['length byebye: ' num2str(length(byebye)) ', length hello: ' num2str(length(hello)) ', size orig xyx: ' num2str(size(xyz_original))]);

%xyz_original
highindex_mod = highindex;

for i=1:length(xyz)
     out = xyz(i,1);
     if out<-1
         oncethrough = 0;
         for j=1:atomtypes
             if i <= highindex(j)
                 highindex_mod(j) = highindex_mod(j)-1;
                 if(~oncethrough) numdifat(j) = numdifat(j) -1;end
                 oncethrough=1;
             end
         end
     end
end
        

while (~isempty(byebye))
    xyz = removerow(xyz, byebye(1));
    byebye = find(xyz(:,1)<-1);
end

final = xyz;

slice=find(final(:,1)<0.2 & final(:,1)>-0.2);
plot3(final(slice,1), final(slice,2), final(slice,3),'.');axis equal


savextl(xtals, filename_xtals, halfedge, numdifat, radii);

save(xyz, filename_newcon, halfedge, numdifat, radii);


%--------------------------------------------------------------------

function fcc = genfcc(numuc)

%do an even layer with the face centers in plane
evenlayer = 't';
atom = 1;
for i=0:numuc %goes up the z axis
    for j=0:numuc %x coord
        for k=0:numuc %y coord
            fcc(atom,1)=j;
            fcc(atom,2)=k;
            fcc(atom,3)=i;
            atom=atom+1;
        end
    end

    for j=0:numuc %x coord
        for k=0:numuc %y coord
            fcc(atom,1)=j+.5;
            fcc(atom,2)=k+.5;
            fcc(atom,3)=i;
            atom=atom+1;
        end
    end
    
    for j=0:numuc %x coord
            for k=0:numuc %y coord
                fcc(atom,1)=j;
                fcc(atom,2)=k+.5;
                fcc(atom,3)=i+.5;
                atom=atom+1;
            end
    end
    
    for j=0:numuc %x coord
            for k=0:numuc %y coord
                fcc(atom,1)=j+.5;
                fcc(atom,2)=k;
                fcc(atom,3)=i+.5;
                atom=atom+1;
            end
    end
end

%randomize the order of the listed positions.
order = randperm(length(fcc));
fcc_o=fcc;
for i=1:length(fcc)
    fcc(i,:) = fcc_o(order(i),:);
end
    

%---------------------------------------------------

function fcc = fccgauss(fcc_o, fwhm)

for i=1:length(fcc_o)
    displace = randn * fwhm/2;
    fcc(i,1) = fcc_o(i,1) + displace; 
    displace = randn * fwhm/2;
    fcc(i,2) = fcc_o(i,2) + displace; 
    displace = randn * fwhm/2;
    fcc(i,3) = fcc_o(i,3) + displace; 
end


%-----------------------------------------------------

function fcc = fccsphere(fcc_o, radius, quadfrac)

%fcc_o comes as cube of edge length 2*radius centered at (000)
%here we have to pick out the atoms that are in a sphere radius
counter=1;
for i=1:length(fcc_o)
    
    %first impose the quadriatic strain (if requested)
    if(quadfrac)
        fcc_o(i,1) = fcc_o(i,1)+ quadfrac * fcc_i(i,1)^2;
        fcc_o(i,2) = fcc_o(i,2)+ quadfrac * fcc_i(i,2)^2;
        fcc_o(i,3) = fcc_o(i,3)+ quadfrac * fcc_i(i,3)^2;        
    end
    
    xsq = fcc_o(i,1)^2;
    ysq = fcc_o(i,2)^2;
    zsq = fcc_o(i,3)^2;
    dist = sqrt(xsq+ysq+zsq);
    if(dist<=radius)
        fcc(counter,:) = fcc_o(i,:);
        counter=counter+1;
    end
end

%three degrees of freedom for x axis
nxi = rand-.5;
nxj = rand-.5;
nxk = rand-.5;

%two degrees for y
nyi = rand-.5;
nyj = rand-.5;
nyk = (-nxi*nyi - nxj*nyj)/nxk;

%one degree for z
nzi = rand-.5;
A = [nyj nyk; nxj nxk];
B = [-nzi*nyi -nzi*nxi]';
temp = A^-1 * B;
nzj = temp(1);
nzk = temp(2);

nx = [nxi nxj nxk];
nx = nx/norm(nx);

ny = [nyi nyj nyk];
ny = ny/norm(ny);

nz = [nzi nzj nzk];
%magnz = sqrt(nzi^2+nzj^2+nzk^2);
%nz = nz/magnz;

dot(nx,ny);
dot(ny,nz);
dot(nx,nz);

%plot3([0 nxi],[0 nxj],[0 nxk],'b');
%plot3([0 nyi],[0 nyj],[0 nyk],'g');
%plot3([0 nzi],[0 nzj],[0 nzk],'c');

nz = cross( nx, ny);
nz = nz/norm(nz);

%plot3([0 nzi],[0 nzj],[0 nzk],'r--');

%hold off


for i=1:length(fcc)
    
    atvec=fcc(i,:);
    
    fcc(i,1) = dot(atvec, nx);
    fcc(i,2) = dot(atvec, ny);
    fcc(i,3) = dot(atvec, nz);
    
end


%------------------------------------------------

function xyz = removerow(xyz_o, index)

num = size(xyz_o,1);
xyz_rest= xyz_o(index+1:num, :);
xyz=xyz_o(1:index-1, :);
xyz(index:num-1, :) = xyz_rest;

%-------------------------------------------------

function xyz = addrow(xyz_o, index, row)
%puts a new row at index+1
num = size(xyz_o,1);
xyz_rest= xyz_o(index+1:num, :);
xyz(1:index,:) = xyz_o(1:index,:);
xyz(index+1,:)= row;
xyz(index+2:num+1,:) = xyz_rest;


%-----------------------------------------------

function save(pos, filename, halfedge, numdifat, radius)

% makes a config file
fid = fopen(filename, 'wt');

rho = length(pos) / (halfedge*2)^3;
fprintf(fid, '%6.6f\n', rho);
fprintf(fid, '%6.0d\n', length(numdifat));

for i=1:length(numdifat)
    fprintf(fid, '%6.0d %6.6f\n', numdifat(i), radius(i));
end

for i=1:length(pos)
    if pos(i,1) >-1
        fprintf(fid, '%6.9f %6.9f %6.9f\n', pos(i,1), pos(i,2), pos(i,3));
    end
end

fclose(fid);

%-----------------------------------------------

function savextl(pos, filename, halfedge, numdifat, radius)
% makes a config file
fid = fopen(filename, 'wt');

rho = length(pos) / (halfedge*2)^3;
fprintf(fid, '%6.6f\n', rho);
fprintf(fid, '%6.0d\n', length(rho));

fprintf(fid, '%6.0d %6.6f\n', length(pos), mean(radius));

for i=1:length(pos)
    fprintf(fid, '%6.9f %6.9f %6.9f\n', pos(i,1), pos(i,2), pos(i,3));
end

fclose(fid);

