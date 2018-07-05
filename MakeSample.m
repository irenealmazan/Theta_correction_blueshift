N = 25; % number of Delaunay mesh points
arr = [ size(X)]; % desired array size

% generating random distribution of polar coordinates, normally distributed
% in the radial direction and uniform in azimuthal and polar directions.
cosTheta = -1 + 2*rand( 1, N );
sinTheta = sin( acos( cosTheta ) );
phi = pi * ( -1 + 2*rand( 1, N ) );
r = min( arr )/5.5 + 0.5.*randn( 1, N );

seed_cosTheta = (cosTheta+1)/2;
seed_phi = (phi/pi+1)/2;
seed_r = (r-min(arr)/5.5)*2;

% 'pts' contains the Delaunay mesh points
pts = [ ...
    r .* sinTheta .* cos( phi ) ; ...
    r .* sinTheta .* sin( phi ) ; ...
    r .* cosTheta .* ones( size( phi ) ) ...
];   

% rotating mesh points by a random rotation
[ R, ~, ~ ] = svd( rand( 3 ) );
pts = ( R * pts )';

% adding central point at origin.
pts = [ [ 0 0 0 ] ; pts ];
pts = pts + repmat( round(arr/2), N+1, 1 );

% generating Voronoi cells of Delaunay mesh points and 
% picking out the central one as synthetic grain.
[ img0, fullMap ] = SampleGenerator.getVoronoiCell(arr, pts);
[ img ] = SampleGenerator.generateRandomPhase( img0, 2,pi/2  );

[rock_curve_3D,intens] = DiffractionPatterns.calc_rock_curve_3DFT(img,addNWstrain,mncntrate);
[img_comp] = DiffractionPatterns.calc_compatible_rho(img,addNWstrain);

DisplayResults.compare_two_objects(img0,img_comp,'Orignal object','Compatible object',[65],'3',400);

%{
% subplot( 1, 2, 1 );
figure; 
% patch( isosurface( smooth3( img, 'gaussian', 7 ) ), 'EdgeColor', 'none', 'FaceColor', 'y' ); 
isosurface( ...
    X,Y,Z, ...
    smooth3( abs( rho ), 'gaussian', 9 ), 0.5, ...
    smooth3( angle( rho ), 'gaussian', 9 ) ); 
axis square; axis tight; axis off; colorbar; lighting gouraud;
% axis image; 
% axis( [ 0 arr(1) 0 arr(2) 0 arr(3) ] );
camlight;
lighting gouraud;
% axis off;

figure;
isosurface( ...
    X,Y,Z, ...
    smooth3( intens, 'gaussian', 15 ));
axis square; axis tight; axis off; colorbar; lighting gouraud;
% axis image; 
% axis( [ 0 arr(1) 0 arr(2) 0 arr(3) ] );
camlight;
lighting gouraud;
% axis off;


% subplot( 1, 2, 2 );
figure; 
imagesc( abs( fk(:,:,35) ) );
axis image; 
colormap( 'hot' );
colorbar; 


% clearvars -except img intens
% save( '/home/smaddali/ANL/simulatedCrystals/crystal_1.mat' );
%}