function [r,histogram,zpartial]=g2r(colpos,z,nelem,deltar,colrad,q)
% Calculate partial 2D histograms
% Inputs:
% 	colpos = (x,y) positions plus atomic number of each atom in the column
%	nelem = number of elements (total)
%	z = array of atomic numbers
%	deltar = bin size (Angstroms)
%	colrad = radius of column

npartials=0;
for i=1:nelem
	for j=i:nelem
		npartials=npartials+1;
			zpartial(npartials,1)=z(i);
			zpartial(npartials,2)=z(j);
	end
end
zpartial=sort(zpartial,2,'descend');

if npartials~=(nelem).*(nelem+1)./2
		error('You screwed something up!')
end

maxbin=fix(2.*colrad./deltar)+1;
r=[0:deltar:2*colrad];
histogram=zeros(maxbin,npartials);

if length(r)~=length(histogram);
		error('You screwed something up!')
end

natomscol=length(colpos);

for i=1:natomscol-1
	for j=i+1:natomscol
		
        rijsq=(colpos(i,1)-colpos(j,1)).*(colpos(i,1)-colpos(j,1))+(colpos(i,2)-colpos(j,2)).*(colpos(i,2)-colpos(j,2));
		rij=sqrt(rijsq);
        bin=round(rij./deltar)+1;
		zs=[colpos(i,3) colpos(j,3)];
		zs=sort(zs,2,'descend');
		
		partial=find(ismember(zpartial,zs,'rows')==1);
		
        rmri=sqrt(colpos(i,1).*colpos(i,1)+colpos(i,2).*colpos(i,2));
        rmrj=sqrt(colpos(j,1).*colpos(j,1)+colpos(j,2).*colpos(j,2));

        airyargi=2*pi.*rmri.*q;
        airyargj=2*pi.*rmrj.*q;

        ai=besselj(1,airyargi)./airyargi;
        aj=besselj(1,airyargj)./airyargj;
        
        
		if bin<=maxbin
			%histogram(bin,partial)=histogram(bin,partial)+2;
            histogram(bin,partial)=histogram(bin,partial)+2*ai*aj;
		end
	end
end

