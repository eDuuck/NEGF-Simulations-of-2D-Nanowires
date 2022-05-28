function [A] = anti_curl(B,c)
%ANTI_CURL Summary of this function goes here
%   Detailed explanation goes here
if nargin < 2
    c = 0;
end
[width, length] = size(B(:,:,1));
dAydx = ones(width, length)*c;
dAxdy = dAydx - B(:,:,3);
A = zeros(width, length,2);
A(:,:,1) = tril(ones(width))*dAxdy;
A(:,:,2) = dAydx*triu(ones(length));
A(:,:,1) = A(:,:,1) - mean(A(:,:,1),'all');
A(:,:,2) = A(:,:,2) - mean(A(:,:,2),'all');
end

