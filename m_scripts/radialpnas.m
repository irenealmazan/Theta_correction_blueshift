function [v mean vaxis center peakdist linearpix] = radialvar(imagename, frac, ...
    ptsbefpeak, peakmult, linewidth, circlediv, center, peakdistin)

%SH August 06 at IMR
%output v - variance along the k axis that is defined bdy the resolution
%   mean - mean of the intensity along the radial loop
%   vaxis - the k axis
%input imagename - file name string, eg diffpat.bmp
%   frac - the fraction of random points that will be used as the sample
%       population for determining variance
%   resolution - the spacing of the q axis points where the bins occur in
%       units of pixels
%   limit - how many pixels away from the center will the calculations stop
%   center - [x y] vector pointing to the center of the circle.  if
%       unknownn, this value should be left as []


if isa(imagename,'char')
    image = double(imread(imagename)); 
else
    image = double(imagename);
end

displayimage=image;
sizeim=size(image)

colormap(gray);
imagesc(image);
axis equal;

thetalow=[];
thetahigh=[];


if isempty(center) %get the center by hand
    display('double click the center, then press space');
    pause;
    cu_pos = get(gca, 'CurrentPoint');
    cx = round(cu_pos(1,1));
    cy = round(cu_pos(1,2));   
    
    displayimage(cy,cx)=0;
    
    %show the close up of the dp
    colormap(gray);
    dist = 300;
    
    displayimage = drawcircle(displayimage, cx, cy, 200, 65536);
    displayimage = drawcircle(displayimage, cx, cy, 170, 65536);
    displayimage = drawcircle(displayimage, cx, cy, 40, 0);
    imagesc(displayimage); axis equal; axis([cx-dist cx+dist cy-dist cy+dist]);
    
    flag ='q';
    while(flag~='0') %this will loop until an appropriate center is found
    
        flag = input(['if center (' num2str(cx) ',' num2str(cy) ') is good, enter 0, to shift enter r/l/u/d: '],'s');
        if flag ~='0' 
            
            displayimage=image;
            [cx cy]= shift(flag, cx,cy,displayimage);
            displayimage = drawcircle(displayimage,cx, cy, 170, 65535);
            displayimage = drawcircle(displayimage, cx, cy, 200, 65535);
            displayimage = drawcircle(displayimage, cx, cy, 40, 0);
            displayimage = drawcircle(displayimage, cx, cy, 5, 0);
            displayimage(cy, cx)=0;
   
            imagesc(displayimage); axis equal; axis([cx-dist cx+dist cy-dist cy+dist]);
        end
    end
    displayimage=image; %clear the guiding circle from the screen
end


if ~isempty(center)
    cx=round(center(1));
    cy=round(center(2));
end

center(1)=cx;
center(2)=cy;

%pick out the first main peak    
xslice=[1:sizeim(1)];
Islice=image(:,cx);
xslice=xslice(:);
Islice=Islice(:);

if isempty(peakdistin)

    %find the center of the y direction
    for i=1:2
        plot(xslice, Islice, [cy cy],[0 65535]);
        flag=1;
        hc=[];
        hc1=[];
        while flag
            display('double click on left bounds, then press space');
            delete(hc);
            pause;
            cu_pos = get(gca, 'CurrentPoint');
            %peakdist=abs(cu_pos(1,1)-cx)
            leftbound=cu_pos(1,1);
            hc=text(cu_pos(1,1),cu_pos(1,2), '*','Color','r');
            display('double click on right bounds, then press space');
            delete(hc1);
            pause;
            cu_pos = get(gca, 'CurrentPoint');
            rightbound=cu_pos(1,1);
            hc1=text(cu_pos(1,1),cu_pos(1,2), '*','Color','r');
            flag=input('is x position of peak ok? y=0: ');
        end

        leftbound=round(leftbound);
        rightbound=round(rightbound);
        [area,chan_max,fity]=peakarea(xslice(leftbound:rightbound), ... 
            Islice(leftbound:rightbound),1);
        plot(xslice(leftbound-50:rightbound+50),Islice(leftbound-50:rightbound+50), ... 
            xslice(leftbound:rightbound),fity); pause(1);
        peakdist2(i)=abs(chan_max - cy);
        display(['distance from peak: ' num2str(peakdist2(i))]);
    end

    peakdist=(peakdist2(1)+peakdist2(2))/2;
    display(['ave dist to peak: ' num2str(peakdist)]);


    %find the center of the x direction

    yslice=[1:sizeim(2)];
    Iyslice=image(cy,:);
    yslice=yslice(:);
    Iyslice=Iyslice(:);

    for i=1:2
        plot(yslice, Iyslice, [cx cx],[0 65535]);
        flag=1;
        hc=[];
        hc1=[];
        while flag
            display('double click on left bounds, then press space');
            delete(hc);
            pause;
            cu_pos = get(gca, 'CurrentPoint');
            %peakdist=abs(cu_pos(1,1)-cx)
            leftbound=cu_pos(1,1);
            hc=text(cu_pos(1,1),cu_pos(1,2), '*','Color','r');
            display('double click on right bounds, then press space');
            delete(hc1);
            pause;
            cu_pos = get(gca, 'CurrentPoint');
            rightbound=cu_pos(1,1);
            hc1=text(cu_pos(1,1),cu_pos(1,2), '*','Color','r');
            flag=input('is x position of peak ok? y=0: ');
        end

        leftbound=round(leftbound);
        rightbound=round(rightbound);
        [area,chan_max,fity]=peakarea(yslice(leftbound:rightbound), ... 
            Iyslice(leftbound:rightbound),1);
        plot(yslice(leftbound-50:rightbound+50),Iyslice(leftbound-50:rightbound+50), ... 
            yslice(leftbound:rightbound),fity); pause(1);
        peakdist3(i)=abs(chan_max - cx);
        display(['distance from peak: ' num2str(peakdist3(i))]);
    end

    peakdist=(peakdist2(1)+peakdist2(2))/2;
    peakdistx=(peakdist3(1)+peakdist3(2))/2;
    display(['ave dist to peak in y: ' num2str(peakdist)]);
    display(['ave dist to peak in x: ' num2str(peakdistx)]);

    plot(yslice-cx,Iyslice,xslice-cy,Islice); pause

end

if ~isempty(peakdistin) peakdist=peakdistin; end

vspacing=(peakdist/ptsbefpeak);
limit=peakmult*ptsbefpeak;
vaxis = [vspacing:vspacing:limit*vspacing] + cx;


displayimage = drawcircle(displayimage, cx, cy, 5, 0);
imagesc(displayimage); axis equal;

%figure out which angular ranges are to be excluded

flag=input('exclude angular regions? 1=y ');
counter=1;

while(flag==1)
    
    sliceok=0;
    while(sliceok~=1) %make sure the user likes the slice
        
        %also can use ginput for mouse points
        
        display ('double click the lower bound of the excluded pie, then press space')
        pause;
        cu_pos = get(gca, 'CurrentPoint');
        lowerx = round(cu_pos(1,1));
        lowery = round(cu_pos(1,2));   
        [hlowline thetalow(counter)] = drawline(cx,cy,lowerx,lowery,limit,'c');

        display('double click the upper bound of the excluded pie, then press space')
        pause;
        cu_pos = get(gca, 'CurrentPoint');
        lowerx = round(cu_pos(1,1));
        lowery = round(cu_pos(1,2)); 
        [hhighline thetahigh(counter)] = drawline(cx,cy,lowerx,lowery,limit,'c');
    
        sliceok=input('is this slice ok? y=1 ');
        
        if sliceok~=1 %delete the handles
            delete(hlowline);
            delete(hhighline);
            pause(.1);
        end
        flag=input('define another slice? y=1 ');
    end
        
    counter=counter+1; %number of angular exclusions
end

% make axis of points where circles will be collected  
lengthv = length(vaxis);
text(vaxis,cy*ones(1,lengthv), '+','Color','r');

pause(.1);

iny = cy;
outy = cy;

color1=65536;
color2=30000;

displayimage=image;
maxvmag = max(vaxis) - cx;
linearpix(1).inten=0;
linearpix(1).radius=0;
linearpix(1).radax=0;

for i=2:lengthv
    temp=color1;
    color1=color2;
    color2=temp;
    [v(i,:) mean(i,:) displayimage linearpixtemp]=circlevar(cx,cy,vaxis(i),cy,displayimage, ... 
        image,frac,temp,thetalow,thetahigh,linewidth,circlediv);
    
    %linearpix(1:length(linearpixtemp), i) = linearpixtemp;
    linearpix(i).inten = linearpixtemp;
    
%Routine for filtering out the low frequency variations
    meantemp = sum(linearpixtemp)/length(linearpixtemp);
    linearpixtemp = linearpixtemp - meantemp;
    F=fft(linearpixtemp);
    F(1:8)=zeros;
    lengthF= length(F);
    f=ifft(F(1: round(lengthF/2)));
        
    %linearpix(1:length(f), i+length(vaxis)) = f + meantemp;
        
    linearpix(i).inten_ft = f+meantemp;
    linearpix(i).radius = vaxis(i);
    linearpix(i).radax= [1:length(linearpix(i).inten)]*360/length(linearpix(i).inten);
    linearpix(i).radax_ft= [1:length(linearpix(i).inten_ft)]*360/length(linearpix(i).inten_ft);

    imagesc(displayimage);axis equal;
    pause(.01); %allows the stripe to be displayed
end
    

%now do the Wochner PNAS analysis
for i=1:length(linearpix) %for every different Q
    
    %Use the FT'ed numbers
    %numpix=length(linearpix(i).radax_ft);
    %inten = linearpix(i).inten_ft;
    
    %or use the normal numbers (recommended)
    numpix=length(linearpix(i).radax);
    inten = linearpix(i).inten;

    deltarange = [1:numpix]; %increment delta by single pixels
    C=[];
    
    for del=deltarange 
        
        IIdelsum=0;
        IIsum=0;
        Isum=0;
        
        for k=1:numpix
            I=inten(k);
            offsetpix=k+(del-1);
            if offsetpix>numpix 
                offsetpix=offsetpix-numpix;
            end
            Idel=inten(offsetpix);
            
            if I>0 && Idel>0
                IIdelsum= IIdelsum + I*Idel;
                IIsum= IIsum + I*I;
                Isum = Isum + I;
            end
            
        end
        %C(del) = (IIdelsum - IIsum) / (IIsum);
        C(del) = (IIdelsum - Isum^2) / (Isum^2);

    end
    linearpix(i).C=C;
end
        


scaling=max(v)/max(mean);

vaxis=vaxis-cx;
vaxis=vaxis*(.4522/(vspacing*ptsbefpeak));

%figure;plot(vaxis,v,vaxis,mean*scaling,'--')
%if isa(imagename, 'char')
%    title(imagename);
%end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [v mean_real displayimage linearpix]= circlevar(cx,cy, ... 
    inx,iny,displayimage,image,frac,color,thetalow,thetahigh, ...
    linewidth,circlediv);

clear image1d
image1d = [];

%change thetalow and thetahigh to radians
thetalow=thetalow*pi/180;
thetahigh=thetahigh*pi/180;

dist = sqrt( (cx - inx)^2 + (cy -iny) ^2); %nominal radius of the circle

maxarcpix=round(dist*pi*2)*linewidth; %this is the nominal number of points on the circle line


%this makes linearpix an array of intensities in order 
pixtocollect=round(maxarcpix);
theta=0;
div=(pi*2)/pixtocollect;
for i=1:pixtocollect
    theta=theta+div;
    thetapass=1;
    
    xpx= round(dist*cos(theta))+cx;
    ypx= round(dist*sin(theta))+cy;
    
    for j=1:length(thetalow)
        if theta>=thetalow(j) && theta<=thetahigh(j)
            thetapass=0;
        end
    end
            
    if thetapass 
        linearpix(i)=image(ypx,xpx);
    else
        linearpix(i) = 0;
    end
end
%plot([1:length(linearpix)],linearpix);
%pause;


%pick out points along arc at random angles
for i=1:frac*maxarcpix
    
    thetapass=0;
    while thetapass~=1  %figure out if angle is in acceptable range
        theta=rand*2*pi;
        bin=ceil(theta/(2*pi/circlediv));
        thetapass=1;  %assume that it is outside of forbidden ranges
        for j=1:length(thetalow)
            if theta>=thetalow(j) && theta<=thetahigh(j)
                thetapass=0;
            end
        end
    end
    
    if linewidth>1 
        randdist=(dist-linewidth/2) + linewidth*rand;
    else
        randdist=dist;
    end
    
    xpx= round(randdist*cos(theta))+cx;
    ypx= round(randdist*sin(theta))+cy;
    
    if xpx>0 && ypx>0 && xpx <=size(displayimage,2) && ypx <=size(displayimage,1)
        image1d(bin,i) = image(ypx, xpx);
        displayimage(ypx,xpx)=color;
    end
    
end

[v mean_real]=calcvar(image1d);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [v mean_real] = calcvar(image1d)

%here calcvar is a matrix whose rows correspond to the intensities gathered
%in a bin of the circle

%plot([1:length(image1d)], image1d);axis([1 length(image1d) 0 65535]);pause

for k=1:size(image1d,1) %go through each bin
    
    clear image1dtemp;
    counter=1;
    for j=1:size(image1d,2)
        if image1d(k,j)>0
            image1dtemp(counter)=image1d(k,j);
            counter = counter+1;
        end
    end
    
    mean1(k) = mean(image1dtemp);
    mean_real(k) = mean1(k);

    stdstat=std(image1dtemp);

    for i=1:size(image1dtemp)
        if image1dtemp(i) > mean1(k)+3*stdstat
            image1dtemp(i)=mean1(k);
        end
        if image1dtemp(i) < mean1(k)-3*stdstat
            image1dtemp(i)=mean1(k);
        end
    end

    %plot([1:length(image1d)], image1d);axis([1 length(image1d) 0 65535]);pause(.1);

    mean2(k) = mean(image1dtemp.*image1dtemp);
    v(k) = (mean2(k)/(mean1(k)^2)) -1;
end
%v=var(image1d); %statistics definition



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [newcx newcy]= shift( flag, oldcx, oldcy, displayimage)

%graphics funciton to move circle that is dislplayed on screen either
%right, left, up, or down.

if flag=='r' %shift the center to the right
    newcx = oldcx + 1;
    newcy = oldcy;
end

if flag=='l' %shift the center to the left
    newcx = oldcx - 1;
    newcy = oldcy;
end

if flag=='u' %shift the center up
    newcx = oldcx;
    newcy = oldcy+1;
end

if flag=='d' %shift the center down
    newcx = oldcx;
    newcy = oldcy-1;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%function ringhandles = drawcircle(displayimage, cx, cy, radius, color)
function displayimage = drawcircle(displayimage, cx, cy, radius, color)

div_x = 360;
increment = (1/div_x) * 2*radius;
displayimage(cy,cx)=0;

for i=1:div_x
    x = round(cx - radius + i*increment);
    y= round(sqrt(radius^2 - (x-cx)^2) + cy);
    displayimage(y,x)=color;
    displayimage(cy-(y-cy), x)=color;
    
    %ringhandles(2*i -1) = text(x,y,'+','Color',color);
    %ringhandles(2*i) = text(x,cy-(y-cy),'+','Color',color);
end
% draw the point at the very left
x = round(cx - radius);
y= round(sqrt(radius^2 - (x-cx)^2) + cy);
   
%ringhandles(div_x*2+1) = text(x, y, '+', 'Color', color);
displayimage(y,x)=color;
ringhandles=[];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ringhandles theta] = drawline(cx,cy,ptx,pty,length,color)

%figure out the angle
x=ptx-cx;
y=pty-cy;
radius = sqrt( (cx - ptx)^2 + (cy -pty) ^2);
theta=acosd(x/radius);
if y<0
    theta=theta + 2*(180-theta);
end

%draw in the points
numpts=20;
spacing=length/numpts;
ringhandles(1)=text(cx,cy,'+','Color',color);
ringhandles(2)=text(ptx,pty,'+','Color','r');

for i=1:numpts
    rad=i*spacing;
    x1=rad*cosd(theta)+cx;
    y1=rad*sind(theta)+cy;
    ringhandles(i+2)= text(x1,y1,'+','Color',color);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [area,chan_max,fity]=peakarea(xarray,yarray,shape) 
%  PEAKAREA - Area of a peak fit to data
%
%  Input:   xpos  = array of x-data (angle, for instance)
%           ypos  = array of y-data (intensity, for instance)
%			shape = shape parameter (0=Gaussian, 1=Lorentzian, 2=Voigt [default])
%
%  Output:  area  = area under peak

%  functions called:
%  leasqr  (the actual fitting routine)
%  pk_qvt   (the function)
%  pk_qvtdf  (the jacobian of the function)
% 21-APR-2004 Gutted FITVOIGT to produce this function for use in STITCH.
% Here are the FITVOIGT comments:
%  Created: 29-APR-1997 
%  Author :  Sean M. Brennan (Bren@SLAC.stanford.edu)
%  Modifications:
%  14-JUL-99 Todd Hufnagel, JHU (hufnagel@jhu.edu)
%    1) Make default plot linear scale
%    2) Make default 1/sqrt(y) weighting
%    3) Print out all fit results to command window
%    4) Plot correct final result
%    5) Cleaned up help comments
%  10-MAR-02 Todd Hufnagel
%    1) Also return area under peak
%

global verbose
% tell leasqr what you want:  
verbose(1)= 0;  % Don't print output
verbose(2)= 0;  % Don't plot intermediate results
verbose(3)= 0;  % Linear scale (log scale would be 0)
verbose(4)= 0;  % Don't stop prematurely

xarray= xarray(:);
yarray= yarray(:);
npts= length(xarray);
[peak,index]= max(yarray);   % find peak of yarray
xpos= xarray(index);         % that's the xval, too.

% Now find speak width, sigma
delta_y=yarray(2:npts)-yarray(1:npts-1);
delta_x=xarray(2:npts)-xarray(1:npts-1);
int_diff=delta_y./delta_x;
[dum l_peak]=max(int_diff(1:index-1));
[dum r_peak]=min(int_diff(index:npts-1));
r_peak=r_peak+index;
sigma=xarray(r_peak)-xarray(l_peak);

lambda= [xpos peak sigma shape 0 0];

constrain=0.001*ones(6,1);  % only 3 parameters, fit them all
constrain(5)=0;constrain(6)=0; % no background
if (shape==0)|(shape==1)    % Gaussian or Lorentzian
    constrain(4)=0;
else
    shape=0.5;  % Guess 50% Gaussian
end    

% Weighting. For equal weights, make this wt= ones(npts,1);
wt= ones(size(xarray));
%plo
fity=leasqr(xarray,yarray,lambda,'pk_qvt',0.001,50,wt,constrain,'pk_qvtdf');


%SH 12-3-04 spit out location of maximum
splxarray=min(xarray):0.1:max(xarray);
splfity=spline(xarray,fity,splxarray);
[value_max,index_max]=max(splfity);
chan_max=splxarray(index_max);
%SH end corrections


area= trapz(xarray,fity);



