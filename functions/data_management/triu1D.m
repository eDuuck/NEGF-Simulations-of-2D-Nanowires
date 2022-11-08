function N = triu1D(M,side)
%TRIU1D will return a 1 dimensional vector representing the top triangular
%matrix of M. If a 1D-vector is fed into TRIU1D, a restored 2D matrix will
%be returned.

if min(size(M)) ~= 1
    len = size(M,2);
    N = zeros((len^2+len)/2,1);
    for i = 1:len
        startIndex = 1 + (len-(i-2)+len)*(i-1)/2;
        endIndex = (len-(i-1)+len)*(i)/2;
        N(startIndex:endIndex) = M(i,i:end);
    end
else
    if ~exist("side","var")
        side = sqrt(0.25+2*length(M))-0.5;
    end
    N = zeros(side);
    for i = 1:side
        startIndex = 1 + (side-(i-2)+side)*(i-1)/2;
        endIndex = (side-(i-1)+side)*(i)/2;
        N(i,i:end) = M(startIndex:endIndex);
    end
end
end