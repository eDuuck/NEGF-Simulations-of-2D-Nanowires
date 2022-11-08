function [N,method,result] = compress(M,method)
%COMPRESS compresses a matrix M according to specified method and returns
%the data packaged in N.
%   Except returning the compressed matrix, the method chosen is also
%   returned. This is important for the block method which also needs to
%   specify the direction of compression which will be picked as the
%   optimal direction is not mentioned. Moreover, COMPRESS will also return
%   a result if the compression was successful or not. (Should always be now).
% 
%   Supported methods are:
%
%   "QOI" : Based on the QOI image compression. This also reduces the data
%   bytesize to 8-bit.
%   
%   "8-bit" : Simply reduces the matrix to a 8-bit format instead of the
%   64-bit format used in MATLAB. The fastest compression and uncompression.
%
%   "gomp" : (g-compression) Works if the matrix is symmetric along the
%   diagonal. This then only saves the one side and then performes QOI on
%   the remaining data.
%
%   "block" : A simple lossless method that should only be used for
%   matrices that are very homogeneous.
%
%   
%   Support for 16-bit and 32-bit structures are not supported for QOI or
%   GOMP currently. The code needs to be improved for these scenarios if
%   8-bit is not sufficient.


if ~exist('method','var')
    method = {'QOI'};
elseif ~iscell(method)
    method = {method};
end

switch(lower(method{1}))
    case 'block'
        [N,method,result] = block_compression(M,method);
    case 'qoi'
        [N,result] = QOI_compression(M,method);
    case 'gomp'
        %Triu1D converts the top triangular matrix of M to a 1D vector.
        [QOI_comp,result] = QOI_compression(triu1D(M),method);
        N = struct("method","Gomp","height", [], "width", [], "QOI_result",QOI_comp);
        N.method = 'Gomp';
        [N.height,N.width] = size(M);
    case '8-bit'
        N = disc_mat(M,"linear");
        result = true;
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
templength = ceil(2.5*numel(M));
if length(method) == 1
    M = disc_mat(M,"linear");
    temp = zeros(1,templength,'uint8');
    old_values = zeros(2^(6),1);
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
            old_values = zeros(2^(6),1,'uint16'); %Support for these aren't done yet.
            temp = zeros(1,templength,'uint8');
        case 4
            old_values = zeros(2^(6),1,'uint32'); 
            temp = zeros(1,templength,'uint32');
        case 8
            old_values = zeros(2^(6),1,'uint64'); 
            temp = zeros(1,templength,'uint64');
    end
    shiftLen = method{2}*8 - 2;
end
modRange = 2^(byteSize*8);

data = double(reshape(M.matrix,1,[]));

temp(1) = bitshift(0b10, 6);      %First cell value is always a new value.
temp(2) = real(data(1));
temp(3) = imag(data(1));
indPos = mod(3*real(data(1))+5*imag(data(1)),64)+1;
old_values(indPos) = data(1);
tempIndex = 4;
curIndex = 2;

while curIndex <= length(data)
    curData = data(curIndex);
    indPos = mod(3*real(curData)+5*imag(curData),64)+1;
    lastData = data(curIndex-1);
    diff = curData - lastData;

    dataDiff = mod(real(diff)+4,modRange) < 8 && ...
               mod(imag(diff)+4,modRange) < 8;

    if curData == lastData %Repeating values
        runlength = 1;
        compIndex = curIndex-1;
        while runlength < 64
            curIndex = curIndex+1;
            if curIndex == 114
                    debug
            end
            if curIndex > length(data) %Algorithm has reached end of data.
                break
            elseif data(curIndex) == data(compIndex) %Still repeating pattern
                runlength = runlength + 1;
            else
                break;
            end
        end
        temp(tempIndex) = bitshift(0b11, shiftLen) + runlength-1;
        tempIndex = tempIndex + 1;

    elseif curData == old_values(indPos) %Value has been seen earlier.
        temp(tempIndex) = indPos-1; %No bitshift needed as bitshift(0b00, shiftLen) = 0;
        curIndex = curIndex + 1;
        tempIndex = tempIndex + 1;

    elseif dataDiff %Small change from previous value.
        QOI_diff_str = mod(real(diff)+4,modRange)*2^3 + mod(imag(diff)+4,modRange);
        temp(tempIndex) = 1*2^shiftLen + QOI_diff_str; %0b01
        old_values(indPos) = curData;
        curIndex = curIndex + 1;
        tempIndex = tempIndex + 1;

    else
        blockStart = curIndex;
        runlength = 1;
        old_values(indPos) = curData;
        while runlength < 64
            lastData = data(curIndex);
            curIndex = curIndex+1;
            if curIndex > length(data) %Algorithm has reached end of data.
                break
            end
            curData = data(curIndex);
            indPos = mod(3*real(curData)+5*imag(curData),64)+1;
            dataDiff = mod(real(curData)-real(lastData)+4,modRange) < 8 && ...
               mod(imag(curData)-imag(lastData)+4,modRange) < 8;

            if curData ~= lastData && ...
                    curData ~= old_values(indPos)&& ...
                    ~dataDiff

                runlength = runlength + 1;
                old_values(indPos) = curData;
            else
                curIndex = curIndex-1; %Necessary spaghetti code.
                break
            end
        end
        curIndex = curIndex + 1;
        temp(tempIndex) = 2*2^shiftLen + runlength-1; %0b10
        tempIndex = tempIndex + 1;
        for i = 0:runlength-1
            temp(tempIndex) = real(data(blockStart + i));
            tempIndex = tempIndex + 1;
            temp(tempIndex) = imag(data(blockStart + i));
            tempIndex = tempIndex + 1;
        end
    end
end
range = M.range;
N = struct('heigth',height,'width',width,...
    'range',range,'byteSize',M.byteSize,'comp_data',uint8(temp(1:tempIndex-1)));
N.method = 'QOI';
N.discMethod = M.method;
result = true;
if tempIndex-1 > 2.5*length(data)
    result = false;
end
end