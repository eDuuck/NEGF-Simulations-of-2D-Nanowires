function [N,method,result] = compress(M,method)
%COMPRESS Summary of this function goes here
%   Detailed explanation goes here
if ~exist('method','var')
    method = {'block',3};
elseif ~iscell(method)
    method = {method};
end

switch(lower(method{1}))
    case 'block'
        [N,method,result] = block_compression(M,method);
    case 'qoi'
        [N,result] = QOI_compression(M,method);
    case 'gomp'
        fM = fft2(M);
        if length(method) > 2
            redMargin = method{3};
        else
            redMargin = 0;
        end
        fM = fM .* (abs(fM) > (max(abs(fM),[],'all')*redMargin));
        %fM = log(fM);
        [QOI_comp,result] = QOI_compression(triu1D(fM),method);
        N = struct("method","Gomp","height", [], "width", [], "QOI_result",QOI_comp);
        N.method = 'Gomp';
        [N.height,N.width] = size(M);
    otherwise
        error('Specified method of compression not supported')
end
end


function [N,method,result] = block_compression(M,method)
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

function [N,result] = QOI_compression(M,method)
[height,width] = size(M);
templength = ceil(2.5*numel(M));    %Needs to use non-linear discretization if fft is to be used.
if length(method) == 1
    M = lin_discretize(M);
    temp = zeros(1,templength,'uint8');
    old_values = zeros(2^(6),1,'uint8');
    shiftLen = 6;
    byteSize = 1;
else
    M = lin_discretize(M,method{2});
    byteSize = method{2};
    switch(byteSize)
        case 1
            old_values = zeros(2^(6),1,'uint8'); %Standard case
            temp = zeros(1,templength,'uint8');
        case 2
            old_values = zeros(2^(6),1,'uint16'); %Extra accurate
            temp = zeros(1,templength,'uint8');
        case 4
            old_values = zeros(2^(6),1,'uint32'); %Special case accuracy
            temp = zeros(1,templength,'uint32');
        case 8
            old_values = zeros(2^(6),1,'uint64'); %Just stupid.
            temp = zeros(1,templength,'uint64');
    end
    shiftLen = method{2}*8 - 2;
end
A = double(reshape(M.matrix,1,[]));

index_position = @(x) mod(3*real(x)+5*imag(x),64)+1;


temp(1) = bitshift(0b10, 6);      %First cell value is always a new value.
temp(2) = real(A(1));
temp(3) = imag(A(1));
old_values(index_position(A(1))) = A(1);
tempIndex = 4;
curIndex = 2;


while curIndex <= length(A)
    if A(curIndex) == A(curIndex-1)   %Repeating values
        runlength = 1;
        compIndex = curIndex-1;
        while runlength < 64
            curIndex = curIndex+1;
            if curIndex == 114
                    debug
            end
            if curIndex > length(A) %Algorithm has reached end of data.
                break
            elseif A(curIndex) == A(compIndex) %Still repeating pattern
                runlength = runlength + 1;
            else
                break;
            end
        end
        temp(tempIndex) = bitshift(0b11, shiftLen) + runlength-1;
        tempIndex = tempIndex + 1;

    elseif A(curIndex) == old_values(index_position(A(curIndex))) %Value has been seen earlier.
        temp(tempIndex) = index_position(A(curIndex))-1; %No bitshift needed as bitshift(0b00, shiftLen) = 0;
        curIndex = curIndex + 1;
        tempIndex = tempIndex + 1;

    elseif QOI_diff_range(A(curIndex),A(curIndex-1),byteSize) %Small change from previous value.
        temp(tempIndex) = bitshift(0b01, shiftLen) + QOI_diff_str(A(curIndex),A(curIndex-1),byteSize);
        old_values(index_position(A(curIndex))) = A(curIndex);
        curIndex = curIndex + 1;
        tempIndex = tempIndex + 1;

    else
        blockStart = curIndex;
        runlength = 1;
        old_values(index_position(A(curIndex))) = A(curIndex);
        while runlength < 64
            curIndex = curIndex+1;
            if curIndex > length(A) %Algorithm has reached end of data.
                break
            elseif ~(A(curIndex) == A(curIndex-1) || ...
                    A(curIndex) == old_values(index_position(A(curIndex))) || ...
                    QOI_diff_range(A(curIndex),A(curIndex-1),byteSize))

                runlength = runlength + 1;
                old_values(index_position(A(curIndex))) = A(curIndex);
            else
                curIndex = curIndex-1; %Necessary spaghetti code.
                break
            end
        end
        curIndex = curIndex + 1;
        temp(tempIndex) = bitshift(0b10, shiftLen) + runlength-1;
        tempIndex = tempIndex + 1;
        for i = 0:runlength-1
            temp(tempIndex) = real(A(blockStart + i));
            tempIndex = tempIndex + 1;
            temp(tempIndex) = imag(A(blockStart + i));
            tempIndex = tempIndex + 1;
        end
    end
end
range = M.range;
N = struct('heigth',height,'width',width,...
    'range',range,'byteSize',M.byteSize,'comp_data',uint8(temp(1:tempIndex-1)));
N.method = 'QOI';
result = true;
if tempIndex-1 > 2.5*length(A)
    result = false;
end
end

function bool = QOI_diff_range(A,B,byteSize)
modRange = bitshift(1,8);
realDiff = mod(real(A)-real(B)+4,modRange);
imagDiff = mod(imag(A)-imag(B)+4,modRange);

bool = imagDiff < byteSize && realDiff < byteSize;
end

function str = QOI_diff_str(A,B,byteSize)
diff = uint8(A-B);
modRange = bitshift(1,8*byteSize);

str = bitshift(mod(real(diff)+4,modRange),3) + mod(imag(diff)+4,modRange);
end
