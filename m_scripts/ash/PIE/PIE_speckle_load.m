%preallocate the structure
speckle(1:n_row,1:n_col) = struct('A',zeros(qy_space_size_pixels,qx_space_size_pixels));

for l = 1:n_row
    for m = 1:n_col
        speckle_temp = load(strcat(data_path,strcat(strcat('speckle',strcat(num2str(l,'_%1.1d'),num2str(m,'_%1.1d'))),'.mat')));
        speckle(l,m).A = speckle_temp.speckle_temp;
        fprintf(strcat(strcat(strcat('Speckle',strcat(num2str(l,'_%d'), num2str(m,'_%d'))), ' has been loaded.'),'\n'))
    end
end

clear('speckle_temp','create_speckle','load_speckle','l','m')
