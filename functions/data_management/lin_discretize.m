function [result] = lin_discretize(A,byteSize)
%LIN_DISCRETIZE Returns the matrix B that is a linear disctetization of the matrix
%A with the resolution of bytesize. The second returned matrix RANGE is
%descriptive and necesarry to restore A.
%   LIN_DISCRETIZE is doing the same thing as disc_mat.m but does only
%   feature linear discrete space for the returned structure.

A = full(A);
realVal = [min(real(A),[],'all'),max(real(A),[],'all')];
imagVal = [min(imag(A),[],'all'),max(imag(A),[],'all')];
range = [realVal;imagVal];
if ~exist('byteSize','var')
    byteSize = 1;
end
if isequal(imagVal,[0 0])
    transFormMat = (real(A)-realVal(1))*2^(8*byteSize) /diff(realVal);
elseif isequal(realVal,[0 0])
    transFormMat = (imag(A)-imagVal(1))*2^(8*byteSize)*1i/diff(imagVal);
else
    transFormMat = (real(A)-realVal(1))*2^(8*byteSize) /diff(realVal) + ...
                  (imag(A)-imagVal(1))*2^(8*byteSize)*1i/diff(imagVal);
end

switch byteSize
    case 1
        B = uint8(transFormMat);
    case 2
        B = uint16(transFormMat);
    case 4
        B = uint32(transFormMat);
    case 8
        B = uint64(transFormMat);
    case deafult
        error('byteSize needs to be of size{1,2,4,8}')
end
result = struct('matrix',B,'range',range,'byteSize',byteSize);
end

