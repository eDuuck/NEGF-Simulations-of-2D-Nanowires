function [Curr] = NEGF_current(NEGF_results, E)
%NEGF_resistance Summary of this function goes here
%   Detailed explanation goes here
e_charge = 1.602176634E-19;
h = 6.62607015e-34;
T = zeros(1,length(E));
for i = 1:length(E)
    T(i) = NEGF_transmission(NEGF_results{i});
end
Curr = e_charge/h*T*triu([0,diff(E)].'*ones(1,length(E)));
end

