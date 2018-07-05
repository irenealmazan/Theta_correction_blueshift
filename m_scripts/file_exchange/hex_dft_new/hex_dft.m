% This function computes the DFT for hexagonally sampled data
% The input is in  the form of a vector. Please note that a hexagonally 
% sampled data can be represented easily in the form of a vector
% and such a representation has a property that neighboring elements
% in the data are neighboring in the vector as well
% Input :
%------------
% signal_vec : Inpute Signal Vector
% Output :
% out : Output DFT vector
% Author - Ujjwal
% IIIT Hyderabad
% Reference : "Hexagonal Image Processing by Lee Middleton and Jayanthi
% Sivaswamy
% E-mail ID - meetukme@gmail.com
%-------------------------------------------------------------------------
function out = hex_dft(signal_vec)
global N;
global C;
global spectral_mapping;
global spatial_mapping;
len = length(signal_vec); % Length of the Signal
lambda = log(len)/log(7); % Compute the Lambda Value
N = arrayfun(@(x) [3 -2; 2 1]^x, 0:(lambda-1),'UniformOutput',0); % Store the periodicity matrices
N_lambda_inv = inv([3 -2 ; 2 1]^lambda); % for quick reference
C = cat(2,{[0;0;0]},arrayfun(@(x) [0 0 -1; -1 0 0; 0 -1 0]^(7-x) * [1;0;-1],1 : 6,'UniformOutput',0)); % Store the Her's mapping function
address_hex = arrayfun(@(x) str2double(dec2base(x,7)),0:(len-1)); % Store the aggregate hexagonal addresses
spatial_mapping = arrayfun(@(x) map_spatial(address_hex(x)),1:length(address_hex),'UniformOutput',0); % Store the spatial mapping function outputs
spectral_mapping = arrayfun(@(x) map_spectral(address_hex(x)),1:length(address_hex),'UniformOutput',0); % Store the spectral mapping function outputs
out = zeros(1,len); % Initialize the output vector
for counter = 1 : len
    % Next 2 lines implement Eq 3.54
    temp = arrayfun(@(x) signal_vec(x) .* exp(-2 * pi * 1i * transpose(spectral_mapping{counter}) * N_lambda_inv * spatial_mapping{x}) ,1:len,'UniformOutput',0);
    out(counter) = sum(cat(3,temp{:}),3); % Store in the output vector
end
    