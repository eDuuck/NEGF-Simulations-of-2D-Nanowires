function [result] = lin_discretize(A,byteSize)
%LIN_DISCRETIZE Returns the matrix B that is a linear disctetization of the matrix
%A with the resolution of bytesize. The second returned matrix RANGE is
%descriptive and necesarry to restore A.
%   Detailed explanation goes here
A = full(A);
realVal = [min(real(A),[],'all'),max(real(A),[],'all')];
imagVal = [min(imag(A),[],'all'),max(imag(A),[],'all')];
range = [realVal;imagVal];
if ~exist('byteSize','var')
    byteSize = 1;
end
switch byteSize
    case 1
        B = uint8((real(A)-realVal(1))*255 /diff(realVal) + ...
                  (imag(A)-imagVal(1))*255i/diff(imagVal));
    case 2
        B = uint16((real(A)-realVal(1))*65535 /diff(realVal)+ ...
                   (imag(A)-imagVal(1))*65535i/diff(imagVal));
    case 4
        B = uint32((real(A)-realVal(1))*(2^32-1)/diff(realVal) + ...
                1i*(imag(A)-imagVal(1))*(2^32-1)/diff(imagVal));
    case 8
        B = uint64((real(A)-realVal(1))*(2^64-1)/diff(realVal) + ...
                1i*(imag(A)-imagVal(1))*(2^64-1)/diff(imagVal));
    case deafult
        error('byteSize needs to be of size{1,2,4,8}')
end
result = struct('matrix',B,'range',range,'byteSize',byteSize);
end

