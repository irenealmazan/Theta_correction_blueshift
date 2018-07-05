%{
HIO
CONSTRAINTS:    Fourier modulus
                positive electron density (NOT IF WE'VE A COMPLEX OBJECT)
                object must be inside finite region     
%}

scrsz = get(0,'ScreenSize');

%obj_abs = imread('shapes2.bmp');
obj_abs = imread('dom9_256.bmp');
obj_abs = double(obj_abs(:,:,3));
obj_abs = obj_abs / max(max(obj_abs));

%obj_phase = imread('shapes3.bmp');
obj_phase = imread('dom10_256.bmp');
obj_phase = double(obj_phase(:,:,2));
obj_phase = 2*pi*(0.4*obj_phase/max(max(obj_phase)));

%Plot the object abs val and phase
figure
set(gcf,'Position',[scrsz(1)*100 scrsz(2)*100 scrsz(3)*1 scrsz(4)*1])
subplot(121)
imagesc(obj_abs)
subplot(122)
imagesc(obj_phase)
colormap gray
pause

size_object_x = size(obj_abs,2);
size_object_y = size(obj_abs,1);

%{
pinhole = imread('pinhole2_256.bmp');
pinhole  = double(pinhole (:,:,1));
beam = fftshift(abs(fft2(pinhole)));
%beam = beam / max(max(beam));
beam = log10(1+beam); 
beam = beam / max(max(beam));
%figure
%imagesc(beam)
%break
%}

%%{
%gaussian beam:
[X,Y] = meshgrid(1:1:size_object_y, 1:1:size_object_x);
beam = exp(-((X-256).^2 + (Y-256).^2)/128^2);
%}

%{
beam = ones(size_object_y,size_object_x);
beam(1:80,:)=0;
beam(:,1:80)=0;
beam(:,420:size_object_y)=0;
beam(420:size_object_y,:)=0;
%}

%{
%H = fspecial('disk',30);
H = fspecial('motion',20,45);
beam = imfilter(imfilter(beam,H,'replicate'),H,'replicate');
%}

%use mouse pointer to select support (look up roipoly in help files)
figure
imagesc(beam .* obj_abs)
colormap gray
[x, y, region_selected, xi, yi] = roipoly;
close all
support = region_selected; 
clear('x','y','xi','yi')

%HIO beta parameter. 
beta = 0.75;

%This is the size we want in fourier space, oversampled by os_factor.
os_factor = 3;
size_x = ceil(size_object_x * os_factor);
size_y = ceil(size_object_y * os_factor);

padding_x = (size_x - size_object_x);
padding_y = (size_y - size_object_y);

obj_complex = padarray(beam .* obj_abs .* exp(i .* obj_phase),[padding_y padding_x],'post');
support     = padarray(support,[padding_y padding_x],'post');

inside_of_finite_support    = support;
outside_of_finite_support   = not(support);

%diffraction pattern:
obj_speckle = abs(fft2(obj_complex, size_y, size_x));

g_k         = 0*obj_speckle + 1;
phase       = 0*obj_speckle + 0;
G_k_prime   = obj_speckle;

figure(223)
set(223,'Position',[scrsz(1)*100 scrsz(2)*100 scrsz(3)*1 scrsz(4)*1])
colormap hsv

clear('os_factor','phase', 'obj_abs', 'obj_phase', 'support', 'obj_complex')

%Ross suggestion: force phase to be within some range as a constraint

%Other idea: force phase --> 0 outside of finite support

%Other idea: pick a guessed large support, run some iterations, use the
%approximate reconstructed image as a refinement on the support, run more
%iterations, use the better approximate reconstructed image as a refinement
%on the support, run more iterations, etc.

for kk=1:4
for jj=1:400
    
    g_k_prime_temp = ifft2(G_k_prime,size_y,size_x);

    %g_k_prime_outside = (g_k - beta * g_k_prime_temp) .* outside_of_finite_support;
    %g_k_prime_inside = g_k_prime_temp .* inside_of_finite_support;
    %g_k = g_k_prime_inside + g_k_prime_outside; 
    
    g_k = g_k_prime_temp .* inside_of_finite_support + ...
        (g_k - beta * g_k_prime_temp) .* outside_of_finite_support; 
    
    %{
    %try setting real space phase to zero if outside support:
    
    angle_g_k = angle(g_k) .* inside_of_finite_support; 
    g_k = abs(g_k) .* exp(i .* angle_g_k );
    
    %result: seems to just stagnate, worse than before. try tighter
    %support?
    %}
        
    G_k_prime = obj_speckle .* exp(i*angle(fft2(g_k, size_y, size_x)));
    
    subplot(1,2,1)
    imagesc(abs(g_k(1:size_object_y,1:size_object_x)))
    colorbar
    daspect([1 1 1])
    subplot(1,2,2)
    imagesc(angle(g_k(1:size_object_y,1:size_object_x)))
    colorbar
    daspect([1 1 1])
    title(strcat('iteration #',num2str(jj,'%d')))
    pause(0.1)
    
end

figure
set(gcf,'Position',[scrsz(1)*100 scrsz(2)*100 scrsz(3)*1 scrsz(4)*1])
colormap hsv

subplot(121)
imagesc(abs(g_k(1:size_object_y,1:size_object_x)))
daspect([1 1 1])
subplot(122)
imagesc(angle(g_k(1:size_object_y,1:size_object_x)))
daspect([1 1 1])

[x, y, region_selected, xi, yi] = roipoly;
close all
support = region_selected; 
clear('x','y','xi','yi')

support     = padarray(support,[padding_y padding_x],'post');

inside_of_finite_support    = support;
outside_of_finite_support   = not(support);

%g_k         = 0*obj_speckle + 1;
%phase       = 0*obj_speckle + 0;
%G_k_prime   = obj_speckle;

figure(223)
set(223,'Position',[scrsz(1)*100 scrsz(2)*100 scrsz(3)*1 scrsz(4)*1])
colormap hsv

clear('os_factor','phase', 'obj_abs', 'obj_phase', 'support', 'obj_complex')

end

%how to update support: use previous reconstruction absval, if we're below
%some low intensity value, set the pixel to zero, but if we're above this
%value, set the pixel to one?

%I = imread('circuit.tif');
%imcontour(I,3)

%I = imread('circuit.tif');
%[C,h] = imcontour(I,3);

%{
figure; 
for ii=1:256 
    plot(abs(g_k(ii,1:256))); 
    ylim([0 1.5]); 
    pause(0.1); 
end
%}


%{
%also play around with colormapeditor
colormap jet
pause(0.5)
colormap hot
pause(0.5)
colormap cool
pause(0.5)
colormap spring
pause(0.5)
colormap summer
pause(0.5)
colormap autumn
pause(0.5)
colormap winter
pause(0.5)
colormap gray
pause(0.5)
colormap bone
pause(0.5)
colormap copper
pause(0.5)
colormap pink
pause(0.5)
colormap lines
pause(0.5)
%}