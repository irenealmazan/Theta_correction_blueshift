function [bkgnd, pixels] = hcdf_bkgnd(apdiam, kvec, nat, edge, ffactors_e, pix)

rad = apdiam/2;
Q=.61/rad;

x=rand(nat,1)*edge;
y=rand(nat,1)*edge;

pixspacing= edge/pix;
paxis = [1:pix]* pixspacing - .5*pixspacing;
 
bkgnd= zeros(length(kvec),1);
pixels = zeros(pix,pix,length(kvec)); 

colvol = 0;

plot(x,y,'o');axis equal; hold
pause;

for l=1:pix
    
    xpix = paxis(l);
    display(['doing row ' num2str(l)]);

    for m=1:pix
        ypix = paxis(m);
        
        %find which atoms are in the circle
        for i=1:nat
            xd = abs(x(i)-xpix);
            if xd >edge/2 xd = abs(xd-edge);end
             
            yd = abs(y(i)-ypix);
            if yd >edge/2 yd = abs(yd-edge);end
                
            dist2 = xd^2 + yd^2;
            if dist2 <= rad^2 
                column(counter,:) = [x(i)
        
        for i=1:colvol

            r1 = rand * rad;
            appfunc1 = calc_appfunc(Q,r1);
            at1 = rand;
            if at1 <=.4 at1 = 1; 
            elseif at1 <=.8 && at1 >.4 at1 = 2; 
            elseif at1 >.8 at1 = 3; end

            for j=1:colvol

                r2 = rand *rad;
                appfunc2 = calc_appfunc(Q,r2);

                at2 = rand;
                if at2 <=.4 at2 = 1; 
                elseif at2 <=.8 && at2 >.4 at2 = 2; 
                elseif at2 >.8 at2 = 3; end

                theta = rand * 2*pi;

                r = sqrt((r1-r2*cos(theta))^2 + (r2*sin(theta))^2);

                for k=1:length(kvec)

                    costerm = calc_costerm(kvec(k), r);
                    pixels(l,m,k) = pixels(l,k) + appfunc1*ffactors_e(k,at1) ...
                        *appfunc2*ffactors_e(k,at2) * costerm;
                end %for k

            end %for j
        end %for i 
    end %for m
end %for l
                     
bkgnd = mean(pixels,1);


function appfunc = calc_appfunc(Q,sig)
appfunc = besselj(1, sig*2*pi*Q) / (sig*2*pi*Q);

function costerm = calc_costerm(scatvec, sig)
costerm = besselj(0, 2*pi*scatvec*sig);