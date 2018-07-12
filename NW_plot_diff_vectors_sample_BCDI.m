% In this script we plot the incoming, outgoing wavevectors, the
% corresponding momentum transfer and the shifts in momentum transfer
% associated to the rocking curve measurement

% plot the frame:
%%{
[X_square,Y_square,Z_square] = meshgrid([-Npix/2 Npix/2-1]*d2_bragg,[-Npix/2 Npix/2-1]*d2_bragg,[-depth/2 depth/2-1]*d2_bragg);
X_square_toplot = X_square(:);
Y_square_toplot = Y_square(:);
Z_square_toplot = Z_square(:);

% display the geometry:

fig_num = 3;
figure(fig_num);
DisplayFunctions.display_diff_geom(NW,ki,kf,qbragg,fig_num,X,Y,Z);

