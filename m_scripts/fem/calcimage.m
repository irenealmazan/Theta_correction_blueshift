function intensity=calcintensity(posz,colrad,l,nstepsedge,stepsize,q,k)
%h = waitbarpos(0,'Calculating intensity',0,0);

intensity=zeros(nstepsedge);
for ix=1:nstepsedge
	for iy=1:nstepsedge	
		x=-l/2+(ix-1)*stepsize + .5*stepsize;
		y=-l/2+(iy-1)*stepsize + .5*stepsize;
	
		[colpos,ncol]=columnpos(posz,x,y,colrad,l);

		for iat=1:ncol
			for jat=1:ncol
				rmri=sqrt(colpos(iat,1).*colpos(iat,1)+colpos(iat,2).*colpos(iat,2));
				rmrj=sqrt(colpos(jat,1).*colpos(jat,1)+colpos(jat,2).*colpos(jat,2));
				airyargi=2*pi.*rmri.*q;
				airyargj=2*pi.*rmrj.*q;
				ai=besselj(1,airyargi)./airyargi;
				aj=besselj(1,airyargj)./airyargj;
				rij=sqrt((colpos(iat,1)-colpos(jat,1))*(colpos(iat,1)-colpos(jat,1))+(colpos(iat,2)-colpos(jat,2))*(colpos(iat,2)-colpos(jat,2)));
				
                if(rij > 0)
                
                    %rij=round(rij/.05) * .05;

                    j0=besselj(0,2*pi*k.*rij);

                    %REAL ALGORITHM
                    intensity(iy,ix)=intensity(iy,ix)+ai.*aj.*j0;

                    %TEST CODE
                    %intensity(iy,ix)=intensity(iy,ix)+rij*j0;
                end
                
            end
		end

        display(num2str(intensity(iy,ix)));
        %intensity(iy,ix) = ncol;
    end
%waitbarpos(ix/nstepsedge,h)
end

display(num2str(intensity));
%close(h)