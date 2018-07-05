% This script sets the experimental geometry: the direction of the
% incoming, outgoing beam and the momentum transfer; calculates the values
% of the momentum transfer for each value of the rocking curve; it
% calculates the corresponding spatial mesh in the real space


% intial directions and magnitude of incoming and outgoing wavectors:
kf = [0 0 1];
ki = [0 0 1];
kmag = 2*pi/lam;

[Ry,Rx] = RotationMatrix.detector(-del,gam);

% rotate incoming beam
ki = (Ry * ki.').';
ki = (Rx * ki.').';

% calculate momentum transfer:
qbragg = kf-ki;

% store the Bragg wavevectors and the Bragg momentum transfer in global variables:
ki_o = ki;
kf_o = kf;
qbragg_o = qbragg;

% calculate the momentum transfer shift for each angle of the rocking
% curve:
dqlist = zeros(numel(delta_thscanvals), 3);
ki_list = zeros(numel(delta_thscanvals), 3);
kf_list = zeros(numel(delta_thscanvals), 3);

[dqlist,ki_list,kf_list] = DiffractionPatterns.calc_dqshift_for_given_th(delta_thscanvals,ki_o,kf_o,qbragg_o);

%%
%scale all results by kmag:
ki_o = ki_o*kmag; 
kf_o = kf_o*kmag; 
qbragg_o = qbragg*kmag;
dqlist_o = dqlist * kmag;
ki_o_list = ki_list * kmag;
kf_o_list = ki_list * kmag;

% do the meshgrid:

th_pixel_size = 2*pi/abs((dqlist_o(end,3) - dqlist_o(1,3)));

[X,Y,Z] = meshgrid([-Npix/2:Npix/2-1]*d2_bragg, ...
                    [-Npix/2:Npix/2-1]*d2_bragg,...
                    [-depth/2:depth/2-1]*th_pixel_size);
                
                              
 
%}
