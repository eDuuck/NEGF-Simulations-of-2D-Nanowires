function [GDiag] = FIND(E,H,Sigma)
    %FIND Summary of this function goes here
    %   Detailed explanation goes here
    Ginv = E*eye(size(H))-H-Sigma;
end

