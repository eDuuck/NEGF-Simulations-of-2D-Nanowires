function [remapped_data] = NEGF_result_remap(NEGF_result, data)
%ELECTRON_DENSITY Summary of this function goes here
%   Detailed explanation goes here
if nargin < 2
    data = "electrons";
end

remapped_data = zeros(size(NEGF_result.sample.units));
width = NEGF_result.sample.width;
length = NEGF_result.sample.length;
switch(data)
    case "electrons"
        data_values = real(diag(NEGF_result.Gn));
    case "fermi"
        Gn = real(diag(NEGF_result.getGn()));
        A = diag(real(1i*(NEGF_result.G - NEGF_result.G')));
        data_values = Gn ./A;
    case "spectral_function"
        data_values = diag(real(1i*(NEGF_result.G - NEGF_result.G')));
    otherwise 
        error("Not a supported sort of data to remap.");
end


    
for j = 1:length
    remapped_data(:,j) = data_values(((j-1)*width+1):(j*width));
end

