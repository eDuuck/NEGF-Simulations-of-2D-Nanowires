function [gamma] = gammCalc(inputMatrix)
%gammCalc Calculates the imaginary values of a vector
%   Detailed explanation goes here
gamma = zeros(size(inputMatrix));
for R = 1:size(inputMatrix,3)
    gamma(:,:,R) = real(1i * (inputMatrix(:,:,R) - inputMatrix(:,:,R)'));
end

