function [density_mat] = electron_density(NEGF_result)
%ELECTRON_DENSITY Summary of this function goes here
%   Detailed explanation goes here
density_mat = zeros(size(NEGF_result.sample.units));
width = NEGF_result.sample.width;
length = NEGF_result.sample.length;
density_values = real(diag(NEGF_result.Gn));

for j = 1:length
    density_mat(:,j) = density_values(((j-1)*width+1):(j*width));
end

