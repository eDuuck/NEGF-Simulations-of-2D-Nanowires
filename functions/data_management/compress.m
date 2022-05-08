function [N,method,result] = compress(M,method)
%COMPRESS Summary of this function goes here
%   Detailed explanation goes here
if nargin < 2
    method = {'block',3};
end

if method{1} == 'block'
    N_size = 2; %Minimum size, first 2 digits width and height of M.
    dims = size(M);
    N = dims;
    j = 1;
    if method{2} == 1 || method{2} == 2
        while j <= numel(M)
            if method{2} == 1
                cNum = M(mod(floor((j-1)/dims(2)),dims(1))+1, mod(j-1,dims(2))+1);
            else
                cNum = M(mod(j-1,dims(1))+1,mod(floor((j-1)/dims(1)),dims(2))+1);
            end
            k = j;
            compNum = cNum;
            while compNum == cNum
                k = k+1;
                if k > numel(M)
                    break;
                end
                if method{2} == 1
                    compNum = M(mod(floor((k-1)/dims(2)),dims(1))+1, mod(k-1,dims(2))+1);
                else
                    compNum = M(mod(k-1,dims(1))+1,mod(floor((k-1)/dims(1)),dims(2))+1);
                end
            end
            N = [N, cNum, k-j];
            j = k;
        end
    else
        method{2} = 1;
        N1 = compress(M,method);
        method{2} = 2;
        N2 = compress(M,method);
        if length(N1) < length(N2)
            N = N1;
            method{2} = 1;
        else
            N = N2;
        end
    end
    result = true;
    if sum(M~=0,'all')*2 < length(N)
        N = sparse(M);
        result = false;
    end
end
    
end

