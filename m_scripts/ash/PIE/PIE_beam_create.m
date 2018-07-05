ROI = zeros(beam_diameter_pixels_y,beam_diameter_pixels_x);
xpix = round((beam_diameter_pixels_x + 1) / 2);
ypix = round((beam_diameter_pixels_y + 1) / 2);

for ii = 1:beam_diameter_pixels_x
    for jj = 1:beam_diameter_pixels_y
        if ((xpix-ii)/ceil(beam_diameter_pixels_x/2))^2 + ((ypix-jj)/ceil(beam_diameter_pixels_y/2))^2 < 1
            ROI(jj,ii) = 1;
        end
    end
end

beam = padarray(ROI,[floor((size(U,1)-beam_diameter_pixels_y)/2) floor((size(U,2)-beam_diameter_pixels_x)/2)]);
if(size(beam,1) < size(U,1))
    beam = padarray(beam,[1 0],'post');
end
if(size(beam,2) < size(U,2))
    beam = padarray(beam,[0 1],'post');
end

new_width_x = fwhm(1:size(beam,2),beam(round(size(beam,2)/2),:));
new_width_y = fwhm(1:size(beam,1),beam(round(size(beam,1)/2),:));
delta_x = floor(beam_diameter_pixels_x*(1-overlap_x));
delta_y = floor(beam_diameter_pixels_y*(1-overlap_y));

%randoffset_y_pixels = round(random('unif',-new_width_y*0.2,new_width_y*0.2,[n_row n_col]));
%randoffset_y_um = randoffset_y_pixels * um_per_pixel_y;
%randoffset_x_pixels = round(random('unif',-new_width_x*0.2,new_width_x*0.2,[n_row n_col]));
%randoffset_x_um = randoffset_x_pixels * um_per_pixel_x;

%H = fspecial('gaussian',20,20);
%beam = imfilter(beam,H,'replicate');

%in pixels on U array
beam_start_y = -100;
beam_start_x = -100;

true_beam(1:n_row,1:n_col) = struct('A',zeros(Ny,Nx));

%beam_scaled = round(3E9 * beam / sum(sum(beam)));
beam_scaled = 3E4 * beam / sum(sum(beam));

for l = 1:n_row
    for m = 1:n_col
        
        %%{
        true_beam(l,m).A = circshift(beam, [(beam_start_y + delta_y*(l-1)) (beam_start_x + delta_x*(m-1))]); 
        true_beam_temp = true_beam(l,m).A;
        save(strcat(data_path,strcat(strcat('true_beam',strcat(num2str(l,'_%1.1d'),num2str(m,'_%1.1d'))),'.mat')),'true_beam_temp','-mat')
        fprintf(strcat(strcat(strcat('Beam',strcat(num2str(l,'_%d'), num2str(m,'_%d'))), ' has been generated.'),'\n'))
        %}
        
        %{
        %true_beam(l,m).A = circshift(beam_scaled, [ (beam_start_y + delta_y*(l-1) + 0*randoffset_y_pixels(l,m)) ...
                                                    %(beam_start_x + delta_x*(m-1) + 0*randoffset_x_pixels(l,m))]); 
                                                    
        true_beam(l,m).A = circshift(beam_scaled, [(beam_start_y + delta_y*(l-1)) (beam_start_x + delta_x*(m-1))]); 
                                                    
        true_beam_temp = true_beam(l,m).A;
        save(strcat(data_path,strcat(strcat('true_beam',strcat(num2str(l,'_%1.1d'),num2str(m,'_%1.1d'))),'.mat')),'true_beam_temp','-mat')
        fprintf(strcat(strcat(strcat('Beam',strcat(num2str(l,'_%d'), num2str(m,'_%d'))), ' has been generated.'),'\n')) 
        %}
        
        %imagesc(x,y,true_beam(l,m).A .* U_Fe_Gd_Attenuation)
        %imagesc(x,y, U_Fe_Gd_Attenuation)
        %pause
        
        %guessed_beam = circshift(beam_final, [(beam_start_y + delta_vertical*(l-1)) (beam_start_x + delta_horiz*(m-1))]); 
        %save(strcat(data_path,strcat(strcat('guessed_beam',strcat(num2str(l,'_%1.1d'),num2str(m,'_%1.1d'))),'.mat')),'guessed_beam','-mat')
        
    end
end

clear('H','xpix','ypix','ii','jj','ROI')
clear('overlap_x','overlap_y','beam_diameter_um_x','beam_diameter_um_y')
%clear('um_per_pixel_x','um_per_pixel_y','beam_diameter_pixels_x','beam_diameter_pixels_y')