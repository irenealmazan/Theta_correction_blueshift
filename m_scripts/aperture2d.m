function tx = aperture2d(ksi,eta,apsize, apshape, offset) 

tx = zeros(size(ksi)); %transmission function

if(ismember('c', apshape))
    %circle
    %display('imposing circular aperture');
    radius = apsize(1)/2;
    radiussq = radius*radius;
    for(ii=1:numel(tx))
        rad = (offset(1)-ksi(ii))^2 + (offset(2)-eta(ii))^2;
        if rad <= radiussq
            tx(ii) =1;
        end
    end
end

if(ismember('r', apshape))
    %rectangle
    %display('imposing rectangual aperture');
    xwidth = apsize(1);
    ywidth = apsize(2);
    for(ii=1:numel(tx))
        if abs(ksi(ii)-offset(1))<=xwidth/2 && abs(eta(ii)-offset(2))<=ywidth/2
            tx(ii) =1;
        end
    end
end

if(ismember('z', apshape))
    
    %make a zone plate
    delr = apsize(3); %width of the outermost zone
    rmax= apsize(1)/2;  %radius of outermost zone
    rmaxsq = rmax*rmax;
        
    for(ii=1:numel(tx))
        radsq = (offset(1)-ksi(ii))^2 + (offset(2)-eta(ii))^2;
        order = floor( radsq/(2*rmax*delr) );
        if (mod(order,2) && radsq <=rmaxsq)
            tx(ii) = 1;
        end
    end
end

if(ismember('b', apshape)) %put in a beamstop
    %display('putting in a circular beamstop');
    radius = apsize(2)/2;
    radiussq = radius*radius;
    for(ii=1:numel(tx))
        rad = (offset(1)-ksi(ii))^2 + (offset(2)-eta(ii))^2;
        if rad <= radiussq
            tx(ii) =0;
        end
    end
end
