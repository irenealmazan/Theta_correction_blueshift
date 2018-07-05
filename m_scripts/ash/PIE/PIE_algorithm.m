beta = 0.75;
alpha = 1E-4;
N_iterations = 30;
%--------------------------------------------------------------------------
update_function(1:n_row,1:n_col) = struct('A',zeros(Ny,Nx)); 
area_covered_beam = 0;

for l = 1:n_row
    for m = 1:n_col
        
        %visual aid to see beam coverage on the object array:
        area_covered_beam = sign(area_covered_beam + true_beam(l,m).A);
        
        %calculate update function:
        update_function(l,m).A = (abs(true_beam(l,m).A) .* conj(true_beam(l,m).A)) ...
                ./ (max(max(abs(true_beam(l,m).A))) .* (alpha + abs(true_beam(l,m).A) .^ 2));

    end
end
%--------------------------------------------------------------------------
%view_slice = 285:715; %for 1024x1024
%view_slice = 155:595; %for 768x768
%view_slice = 195:535; %for 768x768
%--------------------------------------------------------------------------
%starting guess in real space
rho1 = ones(qy_space_size_pixels,qx_space_size_pixels);

%figure(99)
%set(99,'Position',[scrsz(1)*100 scrsz(2)*100 scrsz(3)*1 scrsz(4)*1])
view_slice = 1:512; %for 768x768

if (cutout_speckle == 0)
    for n = 1:N_iterations
        tic
        for l = 1:n_row
            for m = 1:n_col
                
                %multilply the beam profile to the current object guess to
                %obtain an illumination + sample density that will be FTed
                %into the far field 
                rhoiterate = true_beam(l,m).A .* rho1;
                
                %take this illumination*sample condition and FFT to the
                %detector field.  this will give an estimate of the
                %diffraction pattern and phases we were to observe by
                %hitting the sample in this location with the beam
                F2 = fft2(rhoiterate, qy_space_size_pixels, qx_space_size_pixels);

                %now we modify the diffraction pattern we just got because
                %we know the diffracted amplitudes eminating from this
                %position in the sample from experiment
                phases = angle(F2);
                
                %our new guess of what the object should look like in real
                %space is the IFFT of the reciprocal space guess that
                %includes experimental amplitudes
                rho2 = ifft2(sqrt(speckle(l,m).A) .* exp(i*phases));
                
                %our new guess for this object is our original guess for
                %this iteration + a weighted difference between the
                %object*beam(x,y) estimate before and after applying
                %reciprocal space amplitude constraints. 
                rhotemp=rho1;
                rho1 = rhotemp + update_function(l,m).A .* beta .* (rho2 - rhoiterate);
                pause;
            end
        end
        toc

        subplot(121)
        imagesc(x(view_slice),y(view_slice),abs(rho1(view_slice,view_slice))) 
        daspect([1 1 1])
        title(['current guess for iteration ' num2str(n,'%d')])
        subplot(122)
        %imagesc(x(view_slice),y(view_slice),area_covered_beam(view_slice,view_slice) .* U_Fe_Gd_Attenuation(view_slice,view_slice)) 
        imagesc(x(view_slice),y(view_slice),area_covered_beam(view_slice,view_slice) .* U(view_slice,view_slice)) 
        title('starting object');
        daspect([1 1 1])
        colormap gray
        pause(0.1)

    end
else
    for n = 1:N_iterations
        tic
        for l = 1:n_row
            for m = 1:n_col

                rhoiterate = true_beam(l,m).A .* rho1;
                F2 = fft2(rhoiterate, qy_space_size_pixels, qx_space_size_pixels);

                phases = angle(F2);
                rho2 = ifft2(sqrt(speckle(l,m).A) .* exp(i*phases));
                rho1 = rho1 + update_function(l,m).A .* beta .* (rho2 - rhoiterate);

            end
        end
        toc

        subplot(121)
        imagesc(x(view_slice),y(view_slice),abs(rho1(view_slice,view_slice))) 
        daspect([1 1 1])
        title(num2str(n,'%d'))
        subplot(122)
        %imagesc(x(view_slice),y(view_slice),area_covered_beam(view_slice,view_slice) .* U_Fe_Gd_Attenuation(view_slice,view_slice)) 
        imagesc(x(view_slice),y(view_slice),area_covered_beam(view_slice,view_slice) .* U(view_slice,view_slice)) 
        daspect([1 1 1])
        colormap gray
        pause(0.1)

    end
end

h = gcf;
saveas(h,'~/m_scripts/ash/PIE/Speckle/result.jpg')
%--------------------------------------------------------------------------
