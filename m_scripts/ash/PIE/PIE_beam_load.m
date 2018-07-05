%preallocate the structure
true_beam(1:n_row,1:n_col) = struct('A',zeros(Ny,Nx)); 

for l = 1:n_row
    for m = 1:n_col
        true_beam_temp = load(strcat(data_path,strcat(strcat('true_beam',strcat(num2str(l,'_%1.1d'),num2str(m,'_%1.1d'))),'.mat')));
        true_beam(l,m).A = true_beam_temp.true_beam_temp;
        fprintf(strcat(strcat(strcat('Beam',strcat(num2str(l,'_%d'), num2str(m,'_%d'))), ' has been loaded.'),'\n'))
        
        %check to see if real space pixel dimensions are the same as recip space:
        if (size(true_beam(l,m).A,1) < qy_space_size_pixels)

            true_beam(l,m).A = padarray(true_beam(l,m).A,[qy_space_size_pixels - size(true_beam(l,m).A,1) 0],'post');

        end

        if (size(true_beam(l,m).A,2) < qx_space_size_pixels)

            true_beam(l,m).A = padarray(true_beam(l,m).A,[0 qx_space_size_pixels - size(true_beam(l,m).A,2)],'post');
            
        end

       
    end
end

clear('true_beam_temp','beam_diameter_um_x','beam_diameter_um_y')
clear('overlap_x','overlap_y','load_beam','create_beam','l','m')



