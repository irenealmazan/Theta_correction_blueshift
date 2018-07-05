%--------------------------------------------------------------------------
%initialize:
PIE_init;
%--------------------------------------------------------------------------
%what if we cut out a part of the speckle pattern (beam stop)? set this to 1 to study this
cutout_speckle = 0;

%create or load illumination function/speckle patterns:
create_beam_speckle = 1;

domain_pattern_scaling; %account for absorption of Gd and Fe:
if (create_beam_speckle == 1)
    PIE_beam_create;
    PIE_speckle_create;
elseif (create_beam_speckle ~= 1)
    PIE_beam_load;
    PIE_speckle_load;
    if (cutout_speckle == 1)
        PIE_cutout_speckle_region
    end
end

PIE_algorithm;
