function embedmro(config, sro_fori, sro_frac, mro_fori ,mro_frac, numxtal, filename_xtals, filename_newcon)


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

averadius = sum(radii)/(mean(radii)*halfedge*2);
avelattice = 4*averadius/sqrt(2);
radius = mro_frac*6*averadius + sro_frac*2*averadius;

lowindex(1) = 0;
highindex(1) = numdifat(1);
for i=2:length(numdifat)
    lowindex(i) = lowindex(i-1) + numdifat(i-1);
    highindex(i) = highindex(i-1) + numdifat(i);
end

display(['mro crystals containing 169 atoms will be inserted.']);
display(['they will be about ' num2str(radius*2*halfedge/10) ' nm in diameter.']);

%imbed the fcc in a random spot
sro_xtal(13,3,13)=0;
phi = (1+sqrt(5))/2; %used in the icosahedral distances calculations

for i=1:numxtal
    %display(['crystal ' num2str(i)]);
    
    %make the short range clusters
    for j=1:13
        if(sro_fori(1) == 'f')
            sro_xtal(:,:,j) = genfcc;
            sro_xtal(:,:,j) = sro_xtal(:,:,j) * avelattice *sro_frac;
            sro_xtal(:,:,j) = rotate(sro_xtal(:,:,j));
            %sro_xtal(:,:,j) = broaden(sro_xtal(:,:,j));
        end

        if(sro_fori(1) == 'i')
            sro_xtal(:,:,j) = genicos;
            sro_xtal(:,:,j) = sro_xtal(:,:,j) * sqrt((2*averadius*sro_frac)^2/(phi^2+1));
            sro_xtal(:,:,j) = rotate(sro_xtal(:,:,j));
            %sro_xtal(:,:,j) = broaden(sro_xtal(:,:,j));
        end
        %plot3(sro_xtal(:,1,j),sro_xtal(:,2,j), sro_xtal(:,3,j),'ro');axis equal;pause;
    end
    
    if(mro_fori(1) == 'f')
        mrolattice = sqrt(2)*6*averadius * mro_frac;
        mro_xtal = genfcc;
        mro_xtal = mro_xtal * mrolattice;
        mro_xtal = rotate(mro_xtal);
        %mro_xtal = broaden(mro_xtal);
    end
    
    if(mro_fori(1) == 'i')
        mro_xtal = genicos;
        mro_xtal = mro_xtal * sqrt((6*averadius * mro_frac)^2/(phi^2+1));
        mro_xtal = rotate(mro_xtal);
        %mro_xtal = broaden(mro_xtal);
    end
    
    
    for j=1:length(mro_xtal)
        neworigin = mro_xtal(j,:);
        for k=1:length(sro_xtal)
            xtal((j-1)*13 + k, :) = sro_xtal(k,:,j) + neworigin;
        end
    end
    
    %plot3(xtal(:,1),xtal(:,2),xtal(:,3),'ro');axis equal; pause;
    
    x=rand*2-1;
    y=rand*2-1;
    z=rand*2-1;
    
%     x = 0;
%     y=0;
%     z=0;
    
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
    
    order = randperm(length(oksphereindex));
    oksphereindex_o=oksphereindex;
    
    for i=1:length(oksphereindex)
        oksphereindex(i) = oksphereindex_o(order(i));
    end
    
    xtal(:,1) = xtal(:,1)+x;
    xtal(:,2) = xtal(:,2)+y;
    xtal(:,3) = xtal(:,3)+z;
    
    for i=1:length(xtal)
        if xtal(i,1)<=-1
            xtal(i,1)=xtal(i,1)+2;
        end
        if xtal(i,1)>1
            xtal(i,1)=xtal(i,1)-2;
        end
        
        if xtal(i,2) <= -1
            xtal(i,2) = xtal(i,2)+2;
        end
        if xtal(i,2) > 1
            xtal(i,2) = xtal(i,2)-2;
        end
        
        if xtal(i,3) <= -1
            xtal(i,3) = xtal(i,3)+2;
        end
        if xtal(i,3) > 1
            xtal(i,3) = xtal(i,3)-2;
        end
    end
    
    %display(['size fcc: ' num2str(size(fcc)) ', size oksphere: ' num2str(size(oksphere))]);
    
    %oksphereindex
    %oksphere
    %fcc
    
    if size(xtal,1) > size(oksphere,1)
        difference = size(xtal,1) - size(oksphere,1);
        
        for j=1:size(oksphere,1) %substitute these for existing atom positions
            xyz(oksphereindex(j),:) = xtal(j,:);
        end
        
        numat = length(xyz);
        for j=1:difference
            %randindex = ceil(rand*(difference+size(oksphere,1))); %insert a value between this and the next xyz index
            randindex = ceil(rand*numat);
            %display([' adding an atom ' num2str(randindex) ' coord: ' num2str(fcc(j+size(oksphere,1),:))]);
            xyz = addrow(xyz, randindex, xtal(j+size(oksphere,1),:) );
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
    
    if size(xtal,1) < size(oksphere,1) %need to take values out of xyz
        difference = size(oksphere,1) - size(xtal,1);  
        numat = length(xyz);
        
        for j=1:size(xtal,1) %substitute these for existing atom positions
            xyz(oksphereindex(j),:) = xtal(j,:);
            %xyz(oksphereindex(j),:) = [0 0 0];
        end
        
        for j=1:difference
            %display([' removing an atom ' num2str(oksphereindex(j+size(fcc,1)) ) ]);
            xyz(oksphereindex(j+size(xtal,1)),:)=[-1.2 -1.2 -1.2];
            %display([' ' num2str(xyz(oksphereindex(j+size(fcc,1)),:))]);
            %xyz = removerow(xyz, oksphereindex(size(fcc,1)+j)); %removes this atom from the list
             
        end

    end

    clear ok
    clear xtal
    clear oksphere
    clear oksphereindex

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

function fcc = genfcc

fcc(1,:) = [0 0 0];

fcc(2,:) = [.5 .5 0];
fcc(3,:) = [-.5 -.5 0];
fcc(4,:) = [0 .5 .5];
fcc(5,:) = [-.5 0 .5];
fcc(6,:) = [.5 0 -.5];
fcc(7,:) = [0 -.5 -.5];

fcc(8,:) = [.5 -.5 0 ];
fcc(9,:) = [.5 0 .5];
fcc(10,:) = [0, -.5 .5];

fcc(11,:) = [-.5 .5 0];
fcc(12,:) = [-.5 0 -.5];
fcc(13,:) = [0 .5 -.5];

%randomize the order of the listed positions.
order = randperm(length(fcc));
fcc_o=fcc;
for i=1:length(fcc)
    fcc(i,:) = fcc_o(order(i),:);
end
    

%---------------------------------------------------

function icos = genicos

phi = (1+sqrt(5))/2;

icos(1,:) = [0 1 phi];
icos(2,:) = [0 -1 phi];
icos(3,:) = [0 1 -phi];
icos(4,:) = [0 -1 -phi];
icos(5,:) = [1 phi 0];
icos(6,:) = [-1 phi 0];
icos(7,:) = [1 -phi 0];
icos(8,:) = [-1 -phi 0];
icos(9,:) = [phi 0 1];
icos(10,:) = [phi 0 -1];
icos(11,:) = [-phi 0 1];
icos(12,:) = [-phi 0 -1];
icos(13,:) = [0 0 0];

order = randperm(13);
icos_o = icos;
for i=1:13
    icos(i,:) = icos_o(order(i),:);
end

%----------------------------------------------------

function fcc = fccsphere(fcc_o, radius)

counter=1;
for i=1:length(fcc_o)
    xsq = fcc_o(i,1)^2;
    ysq = fcc_o(i,2)^2;
    zsq = fcc_o(i,3)^2;
    dist = sqrt(xsq+ysq+zsq);
    if(dist<=radius)
        fcc(counter,:) = fcc_o(i,:);
        counter=counter+1;
    end
end

nxi = rand-.5;
nxj = rand-.5;
nxk = rand-.5;

nyi = rand-.5;
nyj = rand-.5;
nyk = (-nxi*nyi - nxj*nyj)/nxk;

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


%---------------------------------------------------

function xtal = rotate(xtal)

%pick a random number from -.5 to .5 for x
nxi = rand-.5;
nxj = rand-.5;
nxk = rand-.5;

%pick two random numbers for y direction, but one is restrained
nyi = rand-.5;
nyj = rand-.5;
nyk = (-nxi*nyi - nxj*nyj)/nxk;

% only one degree of freedom for the z axis at this point
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


for i=1:length(xtal)
    
    atvec=xtal(i,:);
    
    xtal(i,1) = dot(atvec, nx);
    xtal(i,2) = dot(atvec, ny);
    xtal(i,3) = dot(atvec, nz);
    
end


%-----------------------------------------------
function xtal = broaden(xtal)


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

