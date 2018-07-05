function [intensity histogram colpos] = calcintensity(posz,colrad,l,nstepsedge,stepsize,q, ... 
    k,z,nelem,deltar)

h = waitbarpos(0,'Calculating intensity',0,0);

intensity=zeros(nstepsedge, nstepsedge, length(k));
counter=1;
r=0;

ffall = load('ffactors_e_kirk_noheader.txt');
ff(:,1) = spline(ffall(:,1),ffall(:,2), k);
ff(:,2) = spline(ffall(:,1),ffall(:,3), k);
ff(:,3) = spline(ffall(:,1),ffall(:,4), k);

for ix=1:nstepsedge
	for iy=1:nstepsedge	
        
		x=-l/2+(ix-1)*stepsize + .5*stepsize;
		y=-l/2+(iy-1)*stepsize + .5*stepsize;
	
        figure(counter); counter=counter+1;
		[colpos,ncol]=columnpos(posz,x,y,colrad,l);

        [r,histogram,zpartial]=g2r(colpos,z,nelem,deltar,colrad,q);
        plot(r, histogram);pause(1);
        title(['pixel ' num2str(ix) ' ' num2str(iy)]);pause(.1);
        
        for ik=1: length(k)
            
            j0=besselj(0, 2*pi*k(ik).*r)';
            plot(r, j0); pause(.5);
            
            %intensity(iy, ix, ik) = sum( j0.* sum(histogram,2) );
            %intensity(iy, ix, ik) = ncol;
            %intensity(iy, ix, ik) = sum(sum(histogram));
        
            intensity(iy, ix, ik) = 0;
            for il=1:size(zpartial,1)               
                
                %ffe1 = escatfa(zpartial(il,1), .5*k(ik), 's');
                %ffe2 = escatfa(zpartial(il,2), .5*k(ik), 's');
                
                if(zpartial(il,1) == 46) flag1=1; end
                if(zpartial(il,2) == 46) flag2=1; end
                if(zpartial(il,1) == 28) flag1=2; end
                if(zpartial(il,2) == 28) flag2=2; end
                if(zpartial(il,1) == 15) flag1=3; end
                if(zpartial(il,2) == 15) flag2=3; end
                ffe1 = ff(ik,flag1);
                ffe2 = ff(ik,flag2);
                
                %intensity(iy, ix, ik) = intensity(iy, ix, ik) + ... 
                %    sum(histogram(:,il).*j0) * ffe1 * ffe2;
                
                %intensity(iy, ix, ik) = intensity(iy, ix, ik) + ffe1 * ffe2;

                intensity(iy, ix, ik) = intensity(iy, ix, ik) + ... 
                    sum(histogram(:,il).*j0);
            end
            
        end
        
    waitbar(counter/(nstepsedge^2),h)
    end
end
close(h)