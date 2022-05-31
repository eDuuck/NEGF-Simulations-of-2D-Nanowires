function [N,method,result] = compress(M,method)
%COMPRESS Summary of this function goes here
%   Detailed explanation goes here
if ~exist('method','var')
    method = {'block',3};
elseif ~iscell(method)
    method = {method};
end

switch(method{1})
    case 'block'
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
        return;
    case 'QOI'
        [N,result] = QOI_compression(M,method);
    otherwise
        error('Specified method of compression not supported')
end
end

function [N,result] = QOI_compression(M,method)
    [height,width] = size(M);
        if length(method) == 1
            M = lin_discretize(M);
            old_values = uint8(zeros(2^(6),1));
        else
            M = lin_discretize(M,method{2});
            switch(method{2})
                case 1
                    old_values = uint8(zeros(2^(6),1)); %Standard case
                case 2
                    old_values = uint16(zeros(2^(6),1)); %Extra accurate
                case 4
                    old_values = uint32(zeros(2^(6),1)); %Special case accuracy
                case 8
                    old_values = uint64(zeros(2^(6),1)); %Just stupid.
            end
        end
        temp = uint8(zeros(1,ceil(2.5*numel(M.matrix)))); %If length(temp) > length(M), compression has failed.
        A = double(reshape(M.matrix,1,[]));

        index_position = @(x) mod(3*real(double(x))+5*imag(double(x)),64)+1;
        

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
                    if curIndex > length(A) %Algorithm has reached end of data.
                        break
                    end
                    if A(curIndex) == A(compIndex) %Still repeating pattern
                        runlength = runlength + 1;
                    else
                        break;
                    end
                end
                temp(tempIndex) = bitshift(0b11, 6) + runlength-1;
                tempIndex = tempIndex + 1;

            elseif A(curIndex) == old_values(index_position(A(curIndex))) %Value has been seen earlier.
                temp(tempIndex) = bitshift(0b00, 6) + index_position(A(curIndex));
                curIndex = curIndex + 1;
                tempIndex = tempIndex + 1;

            elseif QOI_diff_range(A(curIndex),A(curIndex-1)) %Small change from previous value.
                temp(tempIndex) = bitshift(0b01, 6) + QOI_diff_str(A(curIndex),A(curIndex-1));
                old_values(index_position(A(curIndex))) = A(curIndex);
                curIndex = curIndex + 1;
                tempIndex = tempIndex + 1;
            else
                blockStart = curIndex;
                runlength = 1;
                while runlength < 64
                    curIndex = curIndex+1;
                    if curIndex > length(A) %Algorithm has reached end of data.
                        break
                    end
                    if ~(A(curIndex) == A(curIndex-1) || ...
                       A(curIndex) == old_values(index_position(A(curIndex))) || ...
                       QOI_diff_range(A(curIndex),A(curIndex-1)))

                        runlength = runlength + 1;
                        old_values(index_position(A(curIndex))) = A(curIndex);
                    else
                        break;
                    end
                end
                temp(tempIndex) = bitshift(0b10, 6) + runlength-1;
                tempIndex = tempIndex + 1;
                for i = 0:runlength-1
                    temp(tempIndex) = real(A(blockStart + i));
                    tempIndex = tempIndex + 1;
                    temp(tempIndex) = imag(A(blockStart + i));
                    tempIndex = tempIndex + 1;
                end
            end
        end
        
        N = struct('height',height,'width',width,...
            'range',M.range','byteSize',M.byteSize,'comp_data',uint8(temp(1:tempIndex-1)));
        result = true;
        if tempIndex-1 > length(A)
            result = false;
        end
end

function bool = QOI_diff_range(A,B,byteSize)
    if ~exist('byteSize','var')
        byteSize = 1;
    end
    switch byteSize
        case 1
            bool = 0 <= mod(imag(A-B)+4,255) & mod(imag(A-B)+4,255) <= 7 & ...
                   0 <= mod(real(A-B)+4,255) & mod(real(A-B)+4,255) <= 7;
    end
end

function str = QOI_diff_str(A,B,byteSize)
    if ~exist('byteSize','var')
        byteSize = 1;
    end
    switch byteSize
        case 1
            str = bitshift(mod(real(A-B)+4,255),3) + mod(imag(A-B)+4,255);
    end
end
