function [r, gr] = qiq2gr(q, qiq, rho, fftnum, expqmax);

%SH 12-13-05 takes pair correlation function, fft's it and gives back
%q*i(q)


delq = q(2)- q(1);

index=length(q);

%extend r, gr to reach fftnum
if(fftnum)
    qiq(index+1:fftnum)=0;
    q=[0:delq:(fftnum-1)*delq]';
end

zer = find(q>expqmax);
qiq(zer) = 0;

%plot(q,qiq);

fy= -imag(fft(qiq, max(size(q)) ));

delr = 2*pi/(delq*max(size(q)));
r=([0:1:max(size(q))/2-1])*delr';
size(r)

%qiq=4*pi*rho*Fy*delr;
%qiq=qiq(1:max(size(r))/2,1);

rgr=fy(1:max(size(q))/2,1)';
rgr = rgr*delq ./ (2*pi*pi*rho*r);
rgr(1)=0;

gr = rgr./r +1;