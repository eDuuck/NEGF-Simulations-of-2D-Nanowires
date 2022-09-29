function [T] = NEGF_transmission(NEGF_result)
%NEGF_TRANSMISSION Summary of this function goes here
%   Detailed explanation goes here
G = NEGF_result.getG();
gamma1 = NEGF_result.getGamma(1);
gamma2 = NEGF_result.getGamma(2);
T = real(trace(gamma1*G*gamma2*G'));
end

