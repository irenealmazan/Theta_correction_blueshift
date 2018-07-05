%figure(99)
%set(99,'Position',[scrsz(1)*100 scrsz(2)*100 scrsz(3)*1 scrsz(4)*1])

speckle(1:n_row,1:n_col) = struct('A',zeros(qy_space_size_pixels, qx_space_size_pixels));

%figure(98)
%set(98,'Position',[scrsz(1)*100 scrsz(2)*100 scrsz(3)*1 scrsz(4)*1])

for l = 1:n_row
    for m = 1:n_col
        
        if (load_mat_file == 1)
            speckle(l,m).A = abs(fft2(true_beam(l,m).A .* U_Fe_Gd_Attenuation, qy_space_size_pixels, qx_space_size_pixels)) .^ 2; 
            %speckle(l,m).A = abs(fft2(true_beam(l,m).A .* U , q_space_size, q_space_size)) .^ 2;
            
            subplot(121)
            imagesc(x,y,true_beam(l,m).A .* U_Fe_Gd_Attenuation)
            daspect([1 1 1])
            subplot(122)
            imagesc(qx,qy,fftshift(speckle(l,m).A),[0 5E3])
            daspect([1 1 1])
            colormap gray
            pause
            
            %{
            subplot(121)
            imagesc(x,y,true_beam(l,m).A .* U)
            daspect([1 1 1])
            subplot(122)
            imagesc(qx,qy,speckle(l,m).A,[0 7E3])
            daspect([1 1 1])
            colormap gray
            pause
            %}
            
        elseif (load_image == 1)
            speckle(l,m).A = abs(fft2(true_beam(l,m).A .* U , qy_space_size_pixels, qx_space_size_pixels)) .^ 2;
            
            subplot(121)
            imagesc(x,y,true_beam(l,m).A .* U)
            daspect([1 1 1])
            subplot(122)
            imagesc(qx,qy,fftshift(speckle(l,m).A),[0 7E3])
            daspect([1 1 1])
            colormap gray
            pause(.2)
        end
        
        
        speckle_temp = speckle(l,m).A;
  
        
        save(strcat(data_path, strcat(strcat('speckle', strcat(num2str(l,'_%1.1d'), num2str(m,'_%1.1d'))), '.mat')), 'speckle_temp', '-mat')
        fprintf(strcat(strcat(strcat('Speckle',strcat(num2str(l,'_%d'), num2str(m,'_%d'))), ' has been generated.'),'\n'))

    end
end

clear('l','m','beam_start_y','beam_start_x','create_speckle','load_speckle')

