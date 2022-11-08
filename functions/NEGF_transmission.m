function [T] = NEGF_transmission(NEGF_result,con)
%NEGF_TRANSMISSION calculates the transmission through a sample related to
%a NEGF_result.
%   NEGF_TRANSMISSION(NEGF_result) calculates the transmission for an
%   incoming wave at contact 1.
%
%   NEGF_TRANSMISSION(NEGF_result,con) calculates the transmission for an
%   incoming wave at contact 'con'.

if nargin < 2
    con = 1;
end
G = NEGF_result.getG();
A =1i*(G - G');
sigIn = NEGF_result.getSigmaIn();
sigInSum = zeros(size(G));
for j = 1:length(sigIn)
    sigInSum = sigInSum + sigIn{j};
end
Gn = G *(sigInSum + NEGF_result.getSigma0In()) * G';
T = real(trace(sigIn{con}*A- NEGF_result.getGamma(con)*Gn));
end