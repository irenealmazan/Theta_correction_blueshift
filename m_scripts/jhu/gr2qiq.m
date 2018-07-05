function [q qiq] = gr2qiq(r, gr, rho, halfedge, fftnum);

%SH 12-13-05 takes pair correlation function, fft's it and gives back
%q*i(q)


delr = r(2)- r(1);

index=find(r>=halfedge);
index=min(index);

gr(index:max(size(r)))=1;

%extend r, gr to reach fftnum
if(fftnum)
    gr(index:fftnum)=1;
    r=[0:delr:(fftnum-1)*delr]';
end

rgr=r.*(gr-1);
rgr(1)=0;

Fy= -imag(fft(rgr, max(size(r)) ));

delq = 2*pi/(delr*max(size(r)));
q=([0:1:max(size(r))/2-1])*delq';
qiq=4*pi*rho*Fy*delr;
qiq=qiq(1:max(size(r))/2,1);

