function y=deadfunc(x,p)
% DEADFUNC - Detector deadtime function
%
% x = array of monitor count rates (cps)
% p(1) = deadtime (microseconds), p(2)=scale factor
%

% TCH 7-1-94, originally by Paul Zschack

scale=p(2);
deadtime=p(1);

y = scale.*x.*exp(-scale.*deadtime.*x);
