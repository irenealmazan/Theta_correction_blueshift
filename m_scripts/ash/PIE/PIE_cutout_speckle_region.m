l=1;m=1;
colormap gray
imagesc(log(1+fftshift(speckle(l,m).A)))
[x_cut, y_cut, region_selected, xi_cut, yi_cut] = roipoly;

%region_cut_out = not(region_selected);
region_cut_out = region_selected;

%speckle_cut_out(1:n_row,1:n_col) = struct('A',zeros(qy_space_size_pixels,qx_space_size_pixels));

for l = 1:n_row
    for m = 1:n_col
        
        %speckle_cut_out(l,m).A = fftshift(region_cut_out .* fftshift(speckle(l,m).A));
        speckle(l,m).A = fftshift(region_cut_out .* fftshift(speckle(l,m).A)); %#ok<AGROW>
        
        imagesc(fftshift(speckle(l,m).A),[0 2E5])
        daspect([1 1 1])
        h = gcf;
        saveas(h,strcat(data_path,strcat(strcat('speckle_cutout',strcat(num2str(l,'_%1.1d'),num2str(m,'_%1.1d'))),'.jpg')))
        pause(0.1)
    

    end
end

