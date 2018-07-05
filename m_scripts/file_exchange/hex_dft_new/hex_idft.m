% This function computes the IDFT for hexagonally sampled data
% The input is in  the form of a  DFT vector. Please note that a hexagonally 
% sampled data can be represented easily in the form of a vector
% and such a representation has a property that neighboring elements
% in the data are neighboring in the vector as well
% Input :
%------------
% dft_vec : Inpute dft Vector
% Output :
% out : Output signal vector
% Author - Ujjwal
% IIIT Hyderabad
% Reference : "Hexagonal Image Processing by Lee Middleton and Jayanthi
% Sivaswamy
% E-mail ID - meetukme@gmail.com
%-------------------------------------------------------------------------
function out = hex_idft(dft_vec)
global spectral_mapping;
global spatial_mapping;
len = length(dft_vec); % Store the length of the DFT Vector
out = zeros(1,len); % Initialize the output vector
lambda = log(length(dft_vec))/log(7); % Compute the Lambda Value
N_lambda_inv = inv([3 -2 ; 2 1]^lambda); % Store it for future reference
for counter = 1 : len
    % Next 2 lines implement eq 3.55
     temp = arrayfun(@(x) dft_vec(x) .* exp(2 * pi * 1i * transpose(spectral_mapping{counter}) * N_lambda_inv * spatial_mapping{x}) ,1:len,'UniformOutput',0);
     out(counter) = sum(cat(3,temp{:}),3); % Store the result
end
N = [3 -2; 2 1];
out = 1/(7^(lambda-1)) * real(out / abs(det(N))); % Final Scaling needed 