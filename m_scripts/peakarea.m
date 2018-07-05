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
