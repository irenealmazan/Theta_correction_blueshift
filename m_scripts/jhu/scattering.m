function I = scattering(filename, qm, fpd)

a = load(filename);
natoms= size(a,1);

halfedge = (natoms/.073)^(1/3) /2

x= a(:,1);
y= a(:,2);
z= a(:,3);

numq = size(qm,1);
scattering(1:numq, 1) = zeros;
n=1;

for(k=1:natoms) %origin atom
   
    xo = x(k);
    yo = y(k);
    zo = z(k);
    
    for(l=1:natoms) %reference atom
        
        xr = x(l);
        yr = y(l);
        zr = z(l);
        
        r= sqrt((xo-xr)^2 + (yo-yr)^2 + (zo-zr)^2) * halfedge;
        
        if(r>0)
            
            qr = qm .* r;
            sinterm = (sin(qr)./qr);
            scattering(:,n) = fpd.*fpd.*sinterm;
            n=n+1;
            
        end
    end
end

plot(qm, sum(scattering,2)/natoms, qm, fpd.*conj(fpd));
I=scattering;