function [linearpix2, C, Cfilter] = linpix(linearpix)

%input only one row of linarpix at a time


l=length(linearpix);

linearpix2 = zeros(l,1);

%get rid of any zeros that may be attached to the end of the input array
linearpixtemp = linearpix; 
ltemp = length(find(linearpixtemp~=0));
linearpixtemp = linearpixtemp(1:ltemp);

dphi = (2*pi)/ltemp;

meantemp = sum(linearpixtemp)/length(linearpixtemp);
linearpixtemp = linearpixtemp - meantemp;
F=fft(linearpixtemp);
F(1:8)=zeros;
lengthF= length(F);
f=ifft(F(1: round(lengthF/2)));

linearpix2(1:length(f),1) = real(f) + meantemp;

%plot([1:ltemp], linearpixtemp+meantemp, [1:length(f)]*2, linearpix2(1:length(f), 1));pause;

%now do the integral

C = zeros(ltemp, 1);


for i = 1:ltemp  % do each different small phi
    
    smphitopix = i;  %pixel phase shift
    
    for j= 1:ltemp  %go the length of the array to do the sum for this small phi
        
        phi = j;
        shiftindex = phi+smphitopix;
        if shiftindex > ltemp %greater than lower case L
            shiftindex = shiftindex -ltemp;
        end
        
        tempsum = linearpix(round(phi)) * linearpix(shiftindex) * dphi;

        C(i) = C(i) + tempsum;
    end
end


ltemp = length(find(linearpix2~=0));
pix2 = linearpix2(1:ltemp);
dphi = (2*pi)/ltemp;

%plot([1:ltemp], pix2); pause;

Cfilter = zeros(ltemp,1);

for i = 1:ltemp  % do each different small phi
    
    smphitopix = i;  %pixel phase shift
    
    for j= 1:ltemp  %go the length of the array to do the sum for this small phi
        
        phi = j;
        shiftindex = phi+smphitopix;
        if shiftindex > ltemp %greater than lower case L
            shiftindex = shiftindex -ltemp;
        end
        
        tempsum = linearpix2(round(phi)) * linearpix2(shiftindex) * dphi;

        Cfilter(i) = Cfilter(i) + tempsum;
    end
end

plot([1:ltemp], Cfilter)