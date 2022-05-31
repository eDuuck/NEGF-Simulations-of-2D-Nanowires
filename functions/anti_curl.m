function [A] = anti_curl(B,c)
%ANTI_CURL Summary of this function goes here
%   Detailed explanation goes here
[width, length] = size(B);
A = zeros(width, length,2);

if isequal(0*B,B)
    return;
end
if nargin < 2
    c = 0;
end

dAydx = ones(width, length)*c;
dAxdy = dAydx - B;
A(:,:,1) = tril(ones(width))*dAxdy;
A(:,:,2) = dAydx*triu(ones(length));
A(:,:,1) = A(:,:,1) - mean(A(:,:,1),'all');
A(:,:,2) = A(:,:,2) - mean(A(:,:,2),'all');
end

