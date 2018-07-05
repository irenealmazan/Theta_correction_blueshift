scrsz = get(0,'ScreenSize');
load_mat_file = 0;
load_image = 1;
%--------------------------------------------------------------------------
if (load_mat_file == 1)
    %Load matlab binary file containing simulated domain pattern
    newData1 = load('-mat', '/home/ash/MATLAB/Work/PIE/domain2048x2048.mat');

    vars = fieldnames(newData1);
    for i = 1:length(vars)
        assignin('base', vars{i}, newData1.(vars{i}));
    end
%--------------------------------------------------------------------------
%Load an image to reconstruct:
elseif (load_image == 1)
    %U = rgb2hsv(imread('btvn.bmp'));
    U = double(imread('dom2.bmp'));
    U=U(:,:,1);
    %U = 1*(0+U ./ max(max(U)))/1;

    x=1:1:size(U,2);
    y=1:1:size(U,1);
    um_per_pixel_x = 10/size(U,2);
    um_per_pixel_y = 10/size(U,1);
end
%--------------------------------------------------------------------------
%must have slice_size <= qx_space_size_pixels, qy_space_size_pixels
slice_size = 1:512;
qx_space_size_pixels = 512;
qy_space_size_pixels = 512;
% qx_space_size_pixels = 768;
% qy_space_size_pixels = 768;
%576
%640

U=U(slice_size,slice_size);
x=x(slice_size);
y=y(slice_size);
clear('L','Q','ans','clims','dt','epsilon','g','i','k','l')
clear('mm','n','newData1','q','t','vars','u','N')
%--------------------------------------------------------------------------
%size of domain pattern array (pixels)
Nx = size(U,2);
Ny = size(U,1);
%--------------------------------------------------------------------------
%beam properties. need to go somewhat smaller than you think if blurring 
%the beam using fspecial() filters.
if (load_mat_file == 1)
    beam_diameter_um_x = 10; 
    beam_diameter_um_y = 10;
elseif (load_image == 1)
    beam_diameter_um_x = 5; 
    beam_diameter_um_y = 5;
end
beam_diameter_pixels_x = ceil(beam_diameter_um_x / um_per_pixel_x);
beam_diameter_pixels_y = ceil(beam_diameter_um_y / um_per_pixel_y);
%--------------------------------------------------------------------------
%create qx and qy (reciprocal space) vectors using x and y:
qx_max = 2*pi / um_per_pixel_x;
qy_max = 2*pi / um_per_pixel_y;
%qx = -qx_max : 2*qx_max/(size(x,2)-1) : qx_max;
%qy = -qy_max : 2*qy_max/(size(y,2)-1) : qy_max;
qx = -qx_max : 2*qx_max/(qx_space_size_pixels-1) : qx_max;
qy = -qy_max : 2*qy_max/(qy_space_size_pixels-1) : qy_max;
%--------------------------------------------------------------------------
%overlap parameter:
overlap_x = 0.6;
overlap_y = 0.6;
%--------------------------------------------------------------------------
%number of rows and columns
n_row = 3;
n_col = 3;
%--------------------------------------------------------------------------
%where to read/write speckle patterns from/to
data_path = '~/m_scripts/ash/PIE/Speckle/';
%--------------------------------------------------------------------------