function [atom,vertices,c,atden,nelem,bh]=voronoinn(inputfile,cutoff,minarea,r)
%
% VORONOINN - Voronoi construction for RMC configuration file
%
%   [atom,vertices,atden,nelem]=voronoinn(inputfile,cutoff,minarea)
%
% Reads in an RMC configuration file and performs a Voronoi near-neighbor
% construction on it. Assumes boxlength=2 (-1,1 convention).
%
% Input:
%   inputfile = configuration file name (string)
%   cutoff = distance cutoff for neighbors (Å) (optional; default is
%       five times the hard-sphere radius of the largest atom type)
%   minarea = minimum area (as a fraction of total area of Voronoi
%       polyhedron) that a face can have and still be counted (optional;
%       defaul=0.01 or 1%)
%
% Output:
%   atom =  Structure array (one element for each atom) with information
%           about the atom and its neighbors, in these fields:
%           'pos' = position [x y z]
%           'type' = chemical identity (element number)
%           'vol' = Voronoi cell volume (in box units)
%           'nn' = number of near neighbors (of all types)
%           'vertices' = indices of the vertices of the Voronoi polyhedron
%               around this atom
%           'neighbor' = structure array with information about each
%               neighboring atom with these fields:
%               'id' = index of the neighbor atom
%               'type' = type of the neighbor (element number)
%               'nedges' = number of edges of the shared face of polyhedron
%               'vertices' = indices of vertices common between central
%                   atom and this neighbor
%               'dpos' = relative position of neighbor [dx dy dz] (box
%                   units)
%               'dist' = distance of neighbor (box units)
%   vertices = array of positions [x y z] of Voronoi vertices
%   atden = atomic density (from config file, atoms/Å^3)
%   nelem = number of chemical elements
%

% Requires: READCONFIG,GETNUMFROMSTRING
% Reference: Allen and Tildesley, "Computer Simulation of Liquids" for
% a discussion of the linked list structure for atoms in each subcell
%
% Todd Hufnagel (hufnagel@jhu.edu)
% 21-Jun-07 First working version (based on calcvoronoi4.m)
% 22-Jun-07 Added removal of neighbors for whom the shared face of the
%   Vornoi polyhedron was very small.


%--------------------------------------------------------------------------
% Start by reading in the configuration file

[pos,atden,nelem,natomselem,atrad]=readconfig(inputfile);

%
%--------------------------------------------------------------------------
%
% Setting up:
%
% First, calculate the real size of the simulation box, and determine
% the cutoff for near-neighbor calculations. If it's not specified in
% the function call, the cutoff defaults to 3 times the radius of the
% largest atom.

natoms=length(pos);
halfedge = ((1/atden)*natoms)^(1/3) /2;
rspacing= r(2)-r(1);
lreal=(natoms/(8*atden)).^(1/3);    % real complete box length (Å)
fprintf('%d atoms in model with box %4.2f Å on a side.\r',natoms,2*lreal)

if nargin<2
    cutoff=5.*max(atrad);
    fprintf('Neighbor cutoff not specified; using %4.2f Å ',cutoff)
else
    fprintf('Neighbor cutoff: %4.2f Å ',cutoff)
end
cutoff=cutoff./lreal;
fprintf('(%4.2f in box units.)\r',cutoff)

if nargin<3
    minarea=0.01;
end
fprintf('Minimum fractional area for near-neighbor faces: %5.3f\r',...
        minarea)

% Below, we will set up a linked list to store information about which
% atoms are in which subcell. Here, the question is how many cells should 
% we use? There is no point in having cells of size greater than the
%cutoff distance:

m=floor(2./cutoff);

% But there is no sense in having fewer than fewer than 4 subcells per
% side of the simulation box. Since this is a pathological case, we
% disallow it rather than write separate code to deal with it.

if m<4 
    error('Cutoff must be less than %4.2f Å (%4.2f in box units) to implement cells.'...
        ,lreal./2,0.5)
end
ncell=m.^3;
fprintf('Using %d cells (m=%d).\r',ncell,m);

%
%--------------------------------------------------------------------------
%
% Boundary conditions and calculation of Voronoi polyhedra
%
% Periodic boundary conditions are assumed in RMC, but voronoin (the 
% Matlab function) doesn't know anything about that. So we need to make
% sure that the atoms near the faces of the box have neighbors on all 
% sides. We do this by artifically extending the box, using periodicity.
%
pos=imposeperiodicity(pos);
addedatoms=length(pos)-natoms;
fprintf('%d atoms added to provide periodic boundary conditions (%d total).\r',...
    addedatoms,natoms+addedatoms)

% Now determine the Voronoi cells. v is a numv x 3 array giving the
% coordinates of the vertices (where numv is the number of vertices). 
% c is a vector cell array where each element contains the indices of
% the vertices for that cell. (See 'help voronoin' for more information.)

[v,c]=voronoin(pos);
vertices=v; % just make this assignment to give output easy name

% Calculate the volumes of the Voronoi polyhedra. We also check to make
% sure that the polyhedra are valid in the sense of not having any
% vertices at infinity. (This should never happen, since we impose
% periodic boundary conditions, but you can't be too careful.)

vol=zeros(natoms,1);
for i=1:natoms
    if ~isinf(v(c{i},1))
        x = v(c{i},:);       % View ith Voronoi cell.
        [k,cellv] = convhulln(x);
        vol(i)=cellv;
    else
        warning(['Atom ',i,' has invalid Voronoi cell!'])
        vol(i)=nan;
    end
end

% Plot a histogram of the Voronoi cell volumes for inspection;
% this is a useful reality check

[wowy,wowx]=hist(vol,20);
wowy=wowy./natoms;
bar(wowx,wowy)
title('Distribution of Voronoi cell volumes')

% Report the average cell volume, and warn user if different from model
avgcellvol=sum(vol)./natoms;
fprintf('Average Voronoi cell volume: %6.4f\r',avgcellvol)

atvol=8./natoms;
if abs((avgcellvol-atvol)./atvol)>0.001
    warning('Atomic volume from model is %6.4f\r',atvol);
end
input('Paused; hit RETURN to continue.')

%
%--------------------------------------------------------------------------
%
% Determination of relevant information about each atom, including all
% of its near neighbors
%
% Start by preallocating a structure to hold all of the information about
% each atom. Information includes:
%   pos - position [x y z] (box units)
%   type - chemical identity (element number)
%   vol - Voronoi cell volume (in box units)
%   nn - number of neighbors (of all types)
%   vertices - indices of the Voronoi vertices associated with this atom
%   neighbor - structure containing information about neighboring atoms
%       (id, type, number of edges shared, indices of vertices shared,
%       area of shared face, relative position (dpos), and distance)
%
atom = repmat(struct('id',0,'pos',[0 0 0],'type',0,'vol',0,'nn',0,...
    'vertices',zeros(40,1),'vorind',[0 0 0 0 0 0 0],'neighbor', ... 
    struct('id',0 ,...
    'type',0,'nedges',[],'vertices',zeros(10,1),...
    'area',0,'dpos',zeros(1,3),'dist',0)),natoms,1);

% Now we need to figure out all of the above for each atom. Position,
% Type and volume are easy, since we pretty much already know them:
%
% Assign proper type (i.e. chemical element) to each atom

id=1;
for elem=1:nelem
    for i=id:id+natomselem(elem)-1
        atom(i).type=elem;
    end
    id=id+natomselem(elem);
end

% Assign volume (already calculated above), position, and vertices
for i=1:natoms
    atom(i).id = i;
    atom(i).vol=vol(i);
    atom(i).pos=pos(i,:);
    nvert=length(c{i});
    atom(i).vertices(1:nvert)=c{i};
end

% We set up a linked list to describe which atoms are in each cell.
% The array head (length ncell) is the index of the first atom in 
% each cell, while list (length natoms) is the linked list itself.
% (See Allen and Tildesley, p 151 for more details.)

head=zeros(1,ncell);    % entry point into linked list for each cell
list=zeros(1,natoms);   % linked list
celli=m./2;
% Assign atoms to subcells using linked list:
for i=1:natoms
    icell=1+fix((pos(i,1)+1).*celli) ...
        +fix((pos(i,2)+1).*celli).*m ...
        +fix((pos(i,3)+1).*celli).*m.*m;
    list(i)=head(icell);
    head(icell)=i;
end

%
% Now figure out what the cell neighbors of each cell are
% cellneighbor is of size (ncell,27) and gives the indices of
% each of the 27 cell neighbors. Note that "cell neighhbors" includes
% the cell itself (so that when we index over cells to find atomic
% neighbors, we find potential neighbors within the cell itself)

cellneighbor=zeros(ncell,27);
for icell=1:ncell
        [ix,iy,iz] = ind2sub([m m m],icell);
        cn=[ix-1 iy-1 iz-1;
            ix   iy-1 iz-1;
            ix+1 iy-1 iz-1;
            ix-1 iy   iz-1;
            ix   iy   iz-1;
            ix+1 iy   iz-1;
            ix-1 iy+1 iz-1;
            ix   iy+1 iz-1;
            ix+1 iy+1 iz-1;
            ix-1 iy-1 iz;
            ix   iy-1 iz;
            ix+1 iy-1 iz;
            ix-1 iy   iz;
            ix   iy   iz;
            ix+1 iy   iz;
            ix-1 iy+1 iz;
            ix   iy+1 iz;
            ix+1 iy+1 iz;
            ix-1 iy-1 iz+1;
            ix   iy-1 iz+1;
            ix+1 iy-1 iz+1;
            ix-1 iy   iz+1;
            ix   iy   iz+1;
            ix+1 iy   iz+1;
            ix-1 iy+1 iz+1;
            ix   iy+1 iz+1;
            ix+1 iy+1 iz+1];
        cn(cn==m+1)=1; % impose periodic boundaries
        cn(cn==0)=m;
           
        for j=1:27
            cellneighbor(icell,j)=sub2ind([m m m],cn(j,1),cn(j,2),cn(j,3));
        end
end

%
%--------------------------------------------------------------------------
%
% Figure out neighbors for each atom
%
fprintf('\rDetermining atomic near-neighbors...\n')
%h = waitbar(0,'Determining atomic neighbors');
t = cputime;

bh=zeros(length(r),1);

for iatom=1:natoms           % index over all cells

    if mod(iatom,25)==0
        display(['now processing origin atom number ' num2str(iatom) ' of ' num2str(natoms)]);
    end
    atom(iatom).id=iatom;
    clear neighborlist
    clear extendednn
    counter =1;
    %iatom
    nn=0;
    for jatom=1:natoms      % index over all local cells
        %display([num2str([iatom jatom])])

        %first just record the atoms that are in the 27 cell
        %structure

        atom(iatom).extendednn(counter)=jatom;
        counter=counter+1;

        % The Voronoi polyhedra of neighboring atoms have
        % some vertices in common:

        %SH 7-9-07
        %need to account for boundary conditions. 
        %temporarily move vertices of jth atom to their periodic
        %position by adding or subtracting 2 from current values in
        %x y or z as needed

        delx = atom(iatom).pos(1) - atom(jatom).pos(1);
        dely = atom(iatom).pos(2) - atom(jatom).pos(2);
        delz = atom(iatom).pos(3) - atom(jatom).pos(3);

        jvertices = v(c{jatom},:); %this will hold the temporary vertex positions
        ivertices = v(c{iatom},:);

        %check bc in x direction
        if delx >= 1
            jvertices(:,1) = jvertices(:,1)+2;
            delx = 2-delx;
        end
        if delx < -1
            jvertices(:,1) = jvertices(:,1)-2;
            delx = 2+delx;
        end

        %check bc in y direction
        if dely >= 1
            jvertices(:,2) = jvertices(:,2)+2;
            dely = 2-dely;
        end
        if dely < -1
            jvertices(:,2) = jvertices(:,2)-2;
            dely = 2+dely;
        end

        %check bc in z direction
        if delz >= 1
            jvertices(:,3) = jvertices(:,3)+2;
            delz= 2-delz;
        end
        if delz < -1
            jvertices(:,3) = jvertices(:,3)-2;
            delz = 2+delz;
        end
        
        dist = sqrt(delx^2 + dely^2 + delz^2) * halfedge;
        
        bin = round( dist / rspacing) +1;
    
        if bin<=length(r) & dist>0
            bh(bin) = bh(bin)+1;
        end
        
        sharedvert=ismember(single(ivertices), single(jvertices));
        neighbors = find(sharedvert(:,1)==1 & sharedvert(:,2)==1 ...
            & sharedvert(:,3)==1);
        sharedvert = sharedvert(:,1); %to fit assignment in atom.neighbor.vertices

        %end SH
        
        if ~isempty(neighbors) & (iatom ~= jatom)  % Do we share vertices?
                nn=nn+1;
                neighborlist(nn)=jatom;
                atom(iatom).neighbor(nn).id=jatom;
                atom(iatom).neighbor(nn).type=atom(jatom).type;
                atom(iatom).neighbor(nn).nedges=length(neighbors);
                atom(iatom).neighbor(nn).vertices=c{iatom}(sharedvert);
                dx=pos(jatom,1)-pos(iatom,1);
                dy=pos(jatom,2)-pos(iatom,2);
                dz=pos(jatom,3)-pos(iatom,3);
                atom(iatom).neighbor(nn).dpos(1)=dx;
                atom(iatom).neighbor(nn).dpos(2)=dy;
                atom(iatom).neighbor(nn).dpos(3)=dz;
                dist=sqrt(dx.^2+dy.^2+dz.^2);
                atom(iatom).neighbor(nn).dist=dist;
        end


    end

    if nn>0
        atom(iatom).neighborlist = neighborlist;
        atom(iatom).nn=nn;
    end
    
%waitbar(icell/ncell)
end
%close(h)



fprintf('that took %5.1f s.\r',cputime-t)

% Calculate area of near-neighbor faces (if we're checking this)
if minarea>0
    fprintf('\rEliminating neighbors with small shared faces...')
%    h = waitbar(0,'Eliminating neighbors with small shared faces');
    t = cputime;
    for iatom=1:natoms % Loop over all atoms to calculate face areas
%        fprintf('Atom: %d\r',iatom)
%        fprintf('Areas: ');

        nn=atom(iatom).nn;
        facearea=zeros(nn,1);
        for iface=1:nn % Loop over all faces (i.e. neighbors)
% I found the following little bit of clever code for calculating the area
% of an arbitray polygon on the internet, but didn't note where. It makes
% use of the built-in function POLYAREA, but POLYAREA only knows about
% points in (x,y). The clever bit is an easy way to rotate our
% polygon in (x,y,z) into the (x,y) plane so that POLYAREA can work on it.
            vec1=v(atom(iatom).neighbor(iface).vertices(1),:)-...
                v(atom(iatom).neighbor(iface).vertices(2),:);
            vec2=v(atom(iatom).neighbor(iface).vertices(1),:)-...
                v(atom(iatom).neighbor(iface).vertices(3),:);
            fnorm=cross(vec1,vec2);
            p=v(atom(iatom).neighbor(iface).vertices,:);
 %           size(fnorm)
 %           size(p)
            xyproj=p*null(fnorm);
            k=convhull(xyproj(:,1),xyproj(:,2));
            facearea(iface)=polyarea(xyproj(k,1),xyproj(k,2));
            atom(iatom).neighbor(iface).area=facearea(iface);
 %       fprintf('%6.4f, ',facearea(iface))
        end % iface
        totarea=sum(facearea);

        minfacearea=minarea*totarea;

        iface=1;
 %       [iatom length(atom(iatom).neighbor(iface).area)]

        while iface<=length(atom(iatom).neighbor)
  %      fprintf('%d ',iface)
%        for iface=1:nn  % Now need to check to see if faces are too small
%            iface
            if atom(iatom).neighbor(iface).area<minfacearea
%                fprintf('\rAtom %d face %d has area %6.4f (tot=%6.4f)\r',...
%                    iatom,iface,atom(iatom).neighbor(iface).area,totarea)
                atom(iatom).neighbor(iface)=[];
                atom(iatom).nn=atom(iatom).nn-1;
               
            else
                iface=iface+1;
            end
        end

        
% The following was a double-check of the face removal, but it
% doesn't seem to be needed.
        %        iface=1;
%        while iface<=length(atom(iatom).neighbor)
%            [iface length(atom(iatom).neighbor(iface).area)]

%        for iface=1:atom(iatom).nn  % Now need to check to see if faces are too small
%            if atom(iatom).neighbor(iface).area<minfacearea
%                fprintf('Huh? Atom %d, face %d has area %6.4f (tot=%6.4f)\r',...
%                    iatom,iface,atom(iatom).neighbor(iface).area,totarea)
%            end
%            iface=iface+1;
%        end

        %       fprintf('\rTotal area: %6.4f\r',totarea)
 %       pause
 %   waitbar(iatom/natoms)
    end %iatom


    %close(h)
fprintf('that took %5.1f s.\r',cputime-t)

end % minarea


% Look at some simple statistics about the number of near neighbors
% as a reality check:
nn=zeros(natoms,1);
for i=1:natoms
    nn(i)=length(atom(i).neighbor);
%    fprintf('Atom: %d has %d neighbors\r',i,nn(i))
end

[wowy,wowx]=hist(nn,1:20);
wowy=wowy./natoms;
%bar(wowx,wowy)
title('Distribution of number of neighbors')

avgnn=sum(nn)./natoms;
fprintf('Average number of neighbors: %5.1f\r',avgnn)
%input('Paused; hit RETURN to continue.')



%end % function voronoinn


%--------------------------------------------------------------------------

function pos2=imposeperiodicity(pos)
%
% Impose periodic boundary conditions on RMC configuration for Voronoi
% near-neighbor calculations. Do this by replicating a "shell" about
% two atoms thick around the entire configuration box.

% Start with the box itself:
pos2=pos;

% We don't want to replicate the entire box on each side, which would be
% overkill. Instead, we just use a slice from each side that is about
% two atoms thick. That should be enough to allow reliable calculation
% of the near neigbors near the faces. So we want a slice of thickness
% t=2*atvol^(1/3) where atvol=average atomic volume
natoms=length(pos);
atvol=8./natoms;    % note assumes box side length = 2
t=2.*atvol.^(1/3);

% Now append a slice of this thickness to each of the six faces of the
% box. Do each face in turn, starting with the face at x=+L (again assuming
% that the total size of the box is 2L). First find the atoms within t
% from x=-L,since these are the ones that will be replicated:
j=find(pos(:,1)<=-(1-t));
pos2=[pos2; pos(j,1)+2 pos(j,2) pos(j,3)];


% Now do face near x=-L
j=find(pos(:,1)>=1-t);
pos2=[pos2; pos(j,1)-2 pos(j,2) pos(j,3)];


% Next, do the faces near y=+L and -L:
j=find(pos(:,2)<=-(1-t));
pos2=[pos2; pos(j,1) pos(j,2)+2 pos(j,3)];

j=find(pos(:,2)>=1-t);
pos2=[pos2; pos(j,1) pos(j,2)-2 pos(j,3)];


% Now the faces near z=+L and -L:
j=find(pos(:,3)<=-(1-t));
pos2=[pos2;  pos(j,1) pos(j,2) pos(j,3)+2];

j=find(pos(:,3)>=1-t);
pos2=[pos2; pos(j,1) pos(j,2) pos(j,3)-2];


% That does it for the faces, but there still are the regions near the
% edges and the corners. Start with the edge near x,y=L:
j=find((pos(:,1)<=-(1-t))&(pos(:,2)<=-(1-t)));
pos2=[pos2; pos(j,1)+2 pos(j,2)+2 pos(j,3)];


% Now the remaining x,y pairs:
j=find((pos(:,1)>=1-t)&(pos(:,2)>=1-t));
pos2=[pos2; pos(j,1)-2 pos(j,2)-2 pos(j,3)];

j=find((pos(:,1)<=-(1-t))&(pos(:,2)>=1-t));
pos2=[pos2; pos(j,1)+2 pos(j,2)-2 pos(j,3)];

j=find((pos(:,1)>=1-t)&(pos(:,2)<=-(1-t)));
pos2=[pos2; pos(j,1)-2 pos(j,2)+2 pos(j,3)];


% Now do all of the x,z pairs
j=find((pos(:,1)<=-(1-t))&(pos(:,3)<=-(1-t)));
pos2=[pos2; pos(j,1)+2 pos(j,2) pos(j,3)+2];

j=find((pos(:,1)>=1-t)&(pos(:,3)>=1-t));
pos2=[pos2; pos(j,1)-2 pos(j,2) pos(j,3)-2];

j=find((pos(:,1)<=-(1-t))&(pos(:,3)>=1-t));
pos2=[pos2; pos(j,1)+2 pos(j,2) pos(j,3)-2];

j=find((pos(:,1)>=1-t)&(pos(:,3)<=-(1-t)));
pos2=[pos2; pos(j,1)-2 pos(j,2) pos(j,3)+2];


% Finally, do all of the y,z pairs
j=find((pos(:,2)<=-(1-t))&(pos(:,3)<=-(1-t)));
pos2=[pos2; pos(j,1) pos(j,2)+2 pos(j,3)+2];

j=find((pos(:,2)>=1-t)&(pos(:,3)>=1-t));
pos2=[pos2; pos(j,1) pos(j,2)-2 pos(j,3)-2];

j=find((pos(:,2)<=-(1-t))&(pos(:,3)>=1-t));
pos2=[pos2; pos(j,1) pos(j,2)+2 pos(j,3)-2];

j=find((pos(:,2)>=1-t)&(pos(:,3)<=-(1-t)));
pos2=[pos2; pos(j,1) pos(j,2)-2 pos(j,3)+2];


% Last (and least) do all of the corners:
j=find((pos(:,1)<=-(1-t))&(pos(:,2)<=-(1-t))&(pos(:,3)<=-(1-t)));
pos2=[pos2; pos(j,1)+2 pos(j,2)+2 pos(j,3)+2];

j=find((pos(:,1)>=1-t)&(pos(:,2)<=-(1-t))&(pos(:,3)<=-(1-t)));
pos2=[pos2; pos(j,1)-2 pos(j,2)+2 pos(j,3)+2];

j=find((pos(:,1)<=-(1-t))&(pos(:,2)>=1-t)&(pos(:,3)<=-(1-t)));
pos2=[pos2; pos(j,1)+2 pos(j,2)-2 pos(j,3)+2];

j=find((pos(:,1)<=-(1-t))&(pos(:,2)<=-(1-t))&(pos(:,3)>=1-t));
pos2=[pos2; pos(j,1)+2 pos(j,2)+2 pos(j,3)-2];

j=find((pos(:,1)>=1-t)&(pos(:,2)>=1-t)&(pos(:,3)<=-(1-t)));
pos2=[pos2; pos(j,1)-2 pos(j,2)-2 pos(j,3)+2];

j=find((pos(:,1)>=1-t)&(pos(:,2)<=-(1-t))&(pos(:,3)>=1-t));
pos2=[pos2; pos(j,1)-2 pos(j,2)+2 pos(j,3)-2];

j=find((pos(:,1)<=-(1-t))&(pos(:,2)>=1-t)&(pos(:,3)>=1-t));
pos2=[pos2; pos(j,1)+2 pos(j,2)-2 pos(j,3)-2];

j=find((pos(:,1)>=1-t)&(pos(:,2)>=1-t)&(pos(:,3)>=1-t));
pos2=[pos2; pos(j,1)-2 pos(j,2)-2 pos(j,3)-2];

%end %function imposeperiodicity
%
%--------------------------------------------------------------------------


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

% Requires: GETNUMFROMSTRING

% Todd Hufnagel (hufnagel@jhu.edu) 04-Jun-2007

% ch is going to contain all of the string data - make it global to avoid
% having to pass it back and forth
global ch

% open the file, read in the data, deal with CR/LF, and make the final
% character array
fid=fopen(inputfile,'r');
[dch,count]=fread(fid,inf,'uchar');
fclose(fid);

% Old CR/LF code from readspec.m. I have no idea whether this works
% or not for machines other than Macintosh.
chdat=[dch(:)',char(13)];
% replace lf's with spaces or cr's, depending on platform
id10= find(chdat==char(10));
comp= computer;
if strcmp(comp(1:3),'PCW')|strcmp(comp(1:3),'VAX')|strcmp(comp(1:3),'ALP'),
	chdat(id10)= char(' '*ones(1,length(id10)));
else
	chdat(id10)= char(13*ones(1,length(id10)));
end

ch=char(dch)';

% startline is an array with all of the CRs, which is how we distinguish
% the individual lines. The variable 'line' keeps track of where we are
% in the file
lineno=1;
startline=find(dch==char(10))+1;
startline=[1 startline'];

% read in the header information
atden=readnumline(lineno,startline);
lineno=lineno+1;
nelem=readnumline(lineno,startline);
lineno=lineno+1;
for i=1:nelem
    nums=readnumline(lineno,startline);
    natoms(i)=nums(1); 
    atrad(i)=nums(2);
    lineno=lineno+1;
end

totatoms=sum(natoms);
% read in all atom positions
for i=1:totatoms
    pos(i,:)=readnumline(lineno,startline);
    lineno=lineno+1;
end

% Summary report
fprintf('\r')
fprintf('%d atoms read from file %s\r',totatoms,inputfile)
fprintf('\r')
fprintf('Atomic density (atoms/Å^3): %6.4f\r',atden)
fprintf('\r')
fprintf('Element   Atoms   Fraction   Radius (Å)\r')
fprintf('-------   -----   --------   ----------\r')
for i=1:nelem
    fprintf('   %d      %5d     %4.2f        %4.2f\r',...
        i,natoms(i),natoms(i)./totatoms,atrad(i))
end
fprintf('\r') 

%end % function readconfig

%-------------------------------------------------------------------------
% Function readnumline - read numbers from a line
%
function numbers=readnumline(lineno,startline)

global ch
string=ch(startline(lineno):startline(lineno+1)-1);
numbers=getnumfromstring(string);
%end % function extractnum

    
function numbers=getnumfromstring(inputstring)
%
% GETNUMFROMSTRING - extract numbers from string with other text
%
%   numbers=getnumfromstring(inputstring)
%
% Given an input string with artibtrary text, GETNUMFROMSTRING will
% extract any numbers (floating point or otherwise) contained in the 
% string. The only assumption is that the numbers are deliminated
% by white space (spaces and/or tabs, for instance). Numbers delineated
% by other text will return NaN.
%
% Input:
%       inputstring = arbitrary text string with any number of numbers
%                       (assumed to be delineated by whitespace)
% Output:
%       numbers = array containing all numbers in inputstring
%

% Todd Hufnagel (hufnagel@jhu.edu), 4-Jun-07

% First strip out any insignificant white space (at beginning or end)
inputstring=strtrim(inputstring);

% locate all digits
digits=isstrprop(inputstring,'digit');

% local all remaining white space
wspace=isstrprop(inputstring,'wspace');

% pos1 defines the start of the section of the string to be parsed
pos1=1; 

% numnum is the number of numbers we've found
numnum=0;

% alldone is a flag to indicate we've reached the end of the string
alldone=0;

while pos1<=length(inputstring) && ~alldone
    foundstart=0;           % have we found a starting point for a number?

    while ~foundstart && ~alldone
        if digits(pos1)     % starting a new number with a digit
            foundstart=1;
        elseif (inputstring(pos1)=='-')||(inputstring(pos1)=='+')
            if length(inputstring)>pos1
                if digits(pos1+1)
                    foundstart=1;   % starting number with unary + or -
                end
            else
                alldone=1;          % + or - at end of string
            end
        else
            pos1=pos1+1;
        end

        if length(inputstring)<pos1
                alldone=1;
        end

    end

    if ~alldone             % i.e. if we're not at the end of the string
        pos2=1;             % pos2 is the end of the section to be parsed
                            % and is relative to pos1

        % the end of the section to be parsed is where we find
        % either the next white space or the end of the string
        while length(inputstring)>(pos1+pos2)&&~wspace(pos1+pos2)
            pos2=pos2+1;
        end

        if length(inputstring)<pos1+pos2 || wspace(pos1+pos2)
            pos2=pos2-1;
        end

        % section is a string containing the portion of inputstring to
        % be parsed into a number
        section=inputstring(pos1:pos1+pos2);

        % Luckily for us, Matlab has code to do the parsing
        numnum=numnum+1;    % We've identified another number
        numbers(numnum)=str2double(section);

        % Increment to the next position in the string and start over
        pos1=pos1+pos2+1;
    end
end


