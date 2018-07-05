function [del th phi gam] = qtoA(hvec, alpha, lambda, dspacing)
% stolen from my my caxis code geo_zaxis.c. only handling the constant
% alpha mode part right now...

om = (2*pi/dspacing)*eye(3);    % need this to be input
sigma = 0;                   % need this to be input
tau = 0;                     % need this to be input
alpha = alpha*pi/180;            % need this to be input

hphi = om*hvec';
q = norm(hphi)*lambda/(4*pi);
% normalized hphi used in the loop.
hphi = hphi/norm(hphi);

tth = 2*asin(q);
display(['    2 theta is ' num2str(tth *180/pi) ' deg, ' num2str(tth) ' rad']);

%line 581 in geo_zaxis.c
az = [0; 0; 0;];
az(1) = cos(tau)*sin(sigma);
az(2) =-sin(tau)*sin(sigma);
az(3) = cos(sigma);

npar = hphi'*az %line 614 in c code
argbeta = 2*q*npar - sin(alpha) %line 636 in .c
beta = asin(2*q*npar - sin(alpha)) %line 646 in .c

%estimate eta and del using alpha and beta
eta = alpha; %line 688 in .c
del = beta; %line 694 in .c

% iterative loop for finding nu, chi, eta and del
% this loop starts on line 708
for (j = 1:1000)
   %calculate nu
   argnu = (1-2*q^2)/cos(del);  %line 722
   nu = acos((1-2*q^2)/cos(del)); % line 728
    
   %calculate chi, lines 739 thru 742 in .c code
   t1 = cos(del)*sin(nu);
   t2 = cos(nu)*cos(del)*cos(eta) + sin(del)*sin(eta) - cos(eta);
   t3 = hphi(2)*t1 - hphi(1)*t2;
   t4 = hphi(1)*t1 + hphi(2)*t2;
   chi = atan2(t3,t4);

   %calc eta, lines 780 thru 794
   eta = asin((sin(alpha) + cos(eta)*(cos(chi)*az(2) - sin(chi)*az(1)))/az(3));

   %calc del, line 820
   del = asin((2*q*hphi(3) - 2*q^2*sin(eta))/cos(eta));

end

disp('    DEL       THETA        PHI         GAM')    
disp([nu*180/pi chi*180/pi eta*180/pi del*180/pi])

th= chi*180/pi;
phi=eta*180/pi;
gam=del*180/pi;
del = nu*180/pi;

end