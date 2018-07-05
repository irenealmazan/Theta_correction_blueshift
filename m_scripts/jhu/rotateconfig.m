function [config numdifat]=rotateconfig(configfile,configfileout,periodiccube,nx,ny,nz)

%SH 6-14-08
%
%function [config ]=rotateconfig(config,{periodiccube},{nx,ny,nz})
%
%This fucntion most generally takes any configuration of [x y z]
%coordinates and rotates them arbitrarily in space.  The rotation can be
%random or given by the user.  It was designed to work with RMC
%confgurations that are centered at zero and have boundaries at [-1 1] in
%xyz.  In this case, if you want another such cube to be output with
%rotated atomic positions, then with periodiccube, the program will 
%generate this by using the periodic boundary conditions.  Arbitrary
%rotation, however, does not ensure that all atoms will end up in the new
%cube.  Some will be outside and some that are included will have periodic
%images of themselves cast in the new cube.  To ensure that all atoms are
%in the new cube and that no double-atoms are included, specify the
%rotation manually and limit it to orthogonal rotations e.g.:
% x face up:   nx=[0 0 1], ny=[0 1 0], nz=[-1 0 0]
% y face up:   nx=[1 0 0], ny=[0 0 1], nz=[0 -1 0]
%input { } indicate an optional input:
%   config - matrix of xyz postions, sized N x 3
%   {periodcube} - 0 for simple rotation of config positions onto new axis
%                1 for rotation of cube with periodic buffer.  this will 
%                give a cube with sides [-1 1] 
%   {nx, ny, nz} - each of the form [j k l]. expressed in the original 
%                coordinate system, nx 

if nargin<6
    randmode=1;
else
    randmode=0;
end

if nargin<3
    periodiccube = 0;
end

if nargin<2
    savemode=0;
else
    savemode=1;
end

%first open the config file
fid1=fopen(configfile, 'r');
rho = str2num(fgetl(fid1));
atomtypes = str2num(fgetl(fid1));
%numdifat=atomtypes;

for i=1:atomtypes
    row = str2num(fgetl(fid1));
    numdifat(i) = row(1);
    radii(i)= row(2);
end

config = fscanf(fid1, '%f %f %f', [3 inf]);
config = config';

fclose(fid1);

%need to add a fourth row to config that has the atom type

totatoms = sum(numdifat);

counter =1;
flag=numdifat(counter);
for i=1:totatoms
    if i>flag
        counter=counter+1;
        flag=flag+numdifat(counter);
    end
    config(i,4) = counter;
end


if(randmode)

    %these are the components of the new x axis projected onto the existing
    %axes corresponding to a random distance along the new axis
    %three degrees of freedom for x axis
    nxi = rand-.5;
    nxj = rand-.5;
    nxk = rand-.5;

    %once two components are known for y, the third is also known
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

    nz = cross( nx, ny);
    nz = nz/norm(nz);
end

%just to make sure the vecors the user put in are unit vectors
nx = nx/norm(nx);
ny = ny/norm(ny);
nz = nz/norm(nz);

%make sure that the cube comes back out as a cube with periodic boundary
%conditions
%comment this out if you want your original shape to be intact without
%adding atoms later

if(periodiccube)
    %distance of new [111] vector from origin
    corner_o = [1 1 1];
    corner(1) = dot(corner_o, nx);
    corner(2) = dot(corner_o, ny);
    corner(3) = dot(corner_o, nz);

    %this is the maximum amount that the cube "sticks out" at its maximum point (the
    %corner) from the original cube.  we need to layer the periodic atoms up to
    %this thickness on all sides to be able to carve out a cube 
    diff=norm(corner)-1;  %this is the same number for all orientations

    config_o=config;
    config= imposeperiodicity(config,diff);
end
    
for i=1:length(config)
    
    atvec=config(i,1:3);
    
    config(i,1) = dot(atvec, nx);
    config(i,2) = dot(atvec, ny);
    config(i,3) = dot(atvec, nz);
    
end

if(periodiccube)
    %finally, we throw out atoms that are not in the [-1 +1] range on the
    %original cartesian coordinate system
    config = config(find(config(:,1)>=-1 & config(:,1)<1),:);
    config = config(find(config(:,2)>=-1 & config(:,2)<1),:);
    config = config(find(config(:,3)>=-1 & config(:,3)<1),:);
end

%put the atoms back in atomic order

config_final=[];
for i=1:length(numdifat)
    temp = find(config(:,4)==i);
    numdifat(i) = length(temp);
    config_final = [config_final; config(temp,:)];
end

config=config_final;

if(savemode)
    headerstring =[num2str(rho) '\n' num2str(atomtypes)];
    for i=1:atomtypes
        headerstring = [headerstring '\n' num2str(numdifat(i)) '  ' num2str(radii(i))];
    end
    display(headerstring);
    saveigor(config(:,1:3), configfileout, headerstring);
end
%--------------------------------------------------------------------------

function pos2=imposeperiodicity(pos, t)
%
% Impose periodic boundary conditions on RMC configuration for Voronoi
% near-neighbor calculations. Do this by replicating a "shell" about
% two atoms thick around the entire configuration box.

% Start with the box itself:
pos2=pos;

% Now append a slice of this thickness to each of the six faces of the
% box. Do each face in turn, starting with the face at x=+L (again assuming
% that the total size of the box is 2L). First find the atoms within t
% from x=-L,since these are the ones that will be replicated:
j=find(pos(:,1)<=-(1-t));
pos2=[pos2; pos(j,1)+2 pos(j,2) pos(j,3) pos(j,4)];


% Now do face near x=-L
j=find(pos(:,1)>=1-t);
pos2=[pos2; pos(j,1)-2 pos(j,2) pos(j,3) pos(j,4)];


% Next, do the faces near y=+L and -L:
j=find(pos(:,2)<=-(1-t));
pos2=[pos2; pos(j,1) pos(j,2)+2 pos(j,3) pos(j,4)];

j=find(pos(:,2)>=1-t);
pos2=[pos2; pos(j,1) pos(j,2)-2 pos(j,3) pos(j,4)];


% Now the faces near z=+L and -L:
j=find(pos(:,3)<=-(1-t));
pos2=[pos2;  pos(j,1) pos(j,2) pos(j,3)+2 pos(j,4)];

j=find(pos(:,3)>=1-t);
pos2=[pos2; pos(j,1) pos(j,2) pos(j,3)-2 pos(j,4)];


% That does it for the faces, but there still are the regions near the
% edges and the corners. Start with the edge near x,y=L:
j=find((pos(:,1)<=-(1-t))&(pos(:,2)<=-(1-t)));
pos2=[pos2; pos(j,1)+2 pos(j,2)+2 pos(j,3) pos(j,4)];


% Now the remaining x,y pairs:
j=find((pos(:,1)>=1-t)&(pos(:,2)>=1-t));
pos2=[pos2; pos(j,1)-2 pos(j,2)-2 pos(j,3) pos(j,4)];

j=find((pos(:,1)<=-(1-t))&(pos(:,2)>=1-t));
pos2=[pos2; pos(j,1)+2 pos(j,2)-2 pos(j,3) pos(j,4)];

j=find((pos(:,1)>=1-t)&(pos(:,2)<=-(1-t)));
pos2=[pos2; pos(j,1)-2 pos(j,2)+2 pos(j,3) pos(j,4)];


% Now do all of the x,z pairs
j=find((pos(:,1)<=-(1-t))&(pos(:,3)<=-(1-t)));
pos2=[pos2; pos(j,1)+2 pos(j,2) pos(j,3)+2 pos(j,4)];

j=find((pos(:,1)>=1-t)&(pos(:,3)>=1-t));
pos2=[pos2; pos(j,1)-2 pos(j,2) pos(j,3)-2 pos(j,4)];

j=find((pos(:,1)<=-(1-t))&(pos(:,3)>=1-t));
pos2=[pos2; pos(j,1)+2 pos(j,2) pos(j,3)-2 pos(j,4)];

j=find((pos(:,1)>=1-t)&(pos(:,3)<=-(1-t)));
pos2=[pos2; pos(j,1)-2 pos(j,2) pos(j,3)+2 pos(j,4)];


% Finally, do all of the y,z pairs
j=find((pos(:,2)<=-(1-t))&(pos(:,3)<=-(1-t)));
pos2=[pos2; pos(j,1) pos(j,2)+2 pos(j,3)+2 pos(j,4)];

j=find((pos(:,2)>=1-t)&(pos(:,3)>=1-t));
pos2=[pos2; pos(j,1) pos(j,2)-2 pos(j,3)-2 pos(j,4)];

j=find((pos(:,2)<=-(1-t))&(pos(:,3)>=1-t));
pos2=[pos2; pos(j,1) pos(j,2)+2 pos(j,3)-2 pos(j,4)];

j=find((pos(:,2)>=1-t)&(pos(:,3)<=-(1-t)));
pos2=[pos2; pos(j,1) pos(j,2)-2 pos(j,3)+2 pos(j,4)];


% Last (and least) do all of the corners:
j=find((pos(:,1)<=-(1-t))&(pos(:,2)<=-(1-t))&(pos(:,3)<=-(1-t)));
pos2=[pos2; pos(j,1)+2 pos(j,2)+2 pos(j,3)+2 pos(j,4)];

j=find((pos(:,1)>=1-t)&(pos(:,2)<=-(1-t))&(pos(:,3)<=-(1-t)));
pos2=[pos2; pos(j,1)-2 pos(j,2)+2 pos(j,3)+2 pos(j,4)];

j=find((pos(:,1)<=-(1-t))&(pos(:,2)>=1-t)&(pos(:,3)<=-(1-t)));
pos2=[pos2; pos(j,1)+2 pos(j,2)-2 pos(j,3)+2 pos(j,4)];

j=find((pos(:,1)<=-(1-t))&(pos(:,2)<=-(1-t))&(pos(:,3)>=1-t));
pos2=[pos2; pos(j,1)+2 pos(j,2)+2 pos(j,3)-2 pos(j,4)];

j=find((pos(:,1)>=1-t)&(pos(:,2)>=1-t)&(pos(:,3)<=-(1-t)));
pos2=[pos2; pos(j,1)-2 pos(j,2)-2 pos(j,3)+2 pos(j,4)];

j=find((pos(:,1)>=1-t)&(pos(:,2)<=-(1-t))&(pos(:,3)>=1-t));
pos2=[pos2; pos(j,1)-2 pos(j,2)+2 pos(j,3)-2 pos(j,4)];

j=find((pos(:,1)<=-(1-t))&(pos(:,2)>=1-t)&(pos(:,3)>=1-t));
pos2=[pos2; pos(j,1)+2 pos(j,2)-2 pos(j,3)-2 pos(j,4)];

j=find((pos(:,1)>=1-t)&(pos(:,2)>=1-t)&(pos(:,3)>=1-t));
pos2=[pos2; pos(j,1)-2 pos(j,2)-2 pos(j,3)-2 pos(j,4)];

%end %function imposeperiodicity
