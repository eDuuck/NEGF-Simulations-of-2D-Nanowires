function [M] = decompress(N,method)
%DECOMPRESS Summary of this function goes here
%   Detailed explanation goes here
if method{1} == 'block'
    dims = N(1:2);
    M = zeros(N(1),N(2));
    i = 1;
    for j = 1:((length(N)-2)/2)
        cNum = N(2*j + 1);
        for k = 1:N(2*j + 2)
            if method{2} == 1
                M(mod(floor((i-1)/dims(2)), dims(1))+1, mod(i-1,dims(2))+1) = cNum;
            else
                M(mod(i-1,dims(1))+1,mod(floor((i-1)/dims(1)),dims(2))+1) = cNum;
            end
            i = i+1;
        end
    end
end

