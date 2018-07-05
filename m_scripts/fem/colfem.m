function [intensity,vk,imean, histogram, colpos]=colfem(config,z,k, colrad,deltar);
% 
% Calculation of fluctuation electron microscopy signal, by calculation
% of dark field images using the method of Dash et al
% (J Phys Cond Matter 15, S2425 (2003))
% Hollow cone illumination (Eqn 4 in Dash)
%
% Inputs:
%	config = file name of configuration (string)
%	z = array of atomic numbers of atoms
%	k = scattering vector magnitude (array)
%	q = radius of objective aperture (1/A)
%
% Outputs:
%	intensity = array of one image for each point in k
%	vk = calculated variance for each point in k
%	imean = mean intensity for each point in k
	
% 28-Mar-2009 Todd Hufnagel (hufnagel@jhu.edu)

% Read in the configuration file, and figure out what the atomic
% number of each atom is. Result is posz = [x,y,z,atomicnumber]
% with one row for each atom in the model
    [pos,atden,nelem,natoms,atrad]=readconfig(config);
    natomstot=sum(natoms);
    posz=zeros(length(pos),4);
    posz(:,1:3)=pos;
    atomn=0;
    for i=1:nelem
        for j=1:natoms(i)
            atomn=atomn+1;
            posz(atomn,4)=z(i);
        end
    end

% CONVERT FROM BOX UNITS TO ANGSTROMS
	l=(natomstot./atden).^(1/3);
	posz(:,1:3)=posz(:,1:3).*l/2;
	
% A few dimensions:
	q=0.61./colrad;	% Radius of the column (Airy disk)
    nstepsedge=fix(l./(colrad)); % How many steps (in real space) to take
    
    nstepsedge=2;
    
	stepsize=l./nstepsedge; % Size of the steps
	npts=nstepsedge.*nstepsedge; % Total number of points in the image
	
	nptsk=length(k);
    intensity=zeros(nstepsedge,nstepsedge,nptsk);
    vk=zeros(nptsk,1);

%DO THE CALC USING THE FULL J0 EQUATION    
%For each scattering vector, calculate an image and return it along
%with the mean intensity and v(k)
 	counter=1;
    for i=1:nptsk
        %figure(counter); counter=counter+1;
        intensity(:,:,i)=calcimage(posz,colrad,l,nstepsedge,stepsize,q,k(i));
        histogram =0;
        colpos =0;
    end
 
%DO IT USING GR APPROXIMATION
    %[intensity histogram colpos]=calcimagegr(posz,colrad,l,nstepsedge,stepsize,q,k, ...
    %    z, nelem, deltar);
    
    for i=1:nptsk
  		imean(i)=mean(mean(intensity(:,:,i)));
 		avgisquared=sum(sum(intensity(:,:,i).^2))./numel(intensity(:,:,i));
        vk(i)=avgisquared/imean(i)^2-1;
     end

    
%		deltar=0.05;	% bin width in angstroms
%		[r,histogram,zpartial]=g2r(colpos,z,nelem,deltar,colrad);
%		subplot(2,5,1);plot(r,histogram(:,1));
%		subplot(2,5,2);plot(r,histogram(:,2));
%		subplot(2,5,3);plot(r,histogram(:,3));
%		subplot(2,5,4);plot(r,histogram(:,4));
%		subplot(2,5,5);plot(r,histogram(:,5));
%		subplot(2,5,6);plot(r,histogram(:,6));
%		subplot(2,5,7);plot(r,histogram(:,7));
%		subplot(2,5,8);plot(r,histogram(:,8));
%		subplot(2,5,9);plot(r,histogram(:,9));
%		subplot(2,5,10);plot(r,histogram(:,10));
			
%		sizegr=size(histogram);
%		for i=1:sizegr(2)
%			plot(r,histogram(:,i));
%			pause
%		end
		
		