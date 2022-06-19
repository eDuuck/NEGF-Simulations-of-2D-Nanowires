function [M] = decompress(N,method,debugstuff)
%DECOMPRESS Summary of this function goes here
%   Detailed explanation goes here
if ~exist("method","var")
    method = {N.method};
end
if ~exist("debugstuff","var")
    debugstuff = [];
end
if ~iscell(method)
    method = {method};
end
switch lower(method{1})
    case 'block'
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
    case 'qoi'
        M = QOI_decompress(N,debugstuff);
        
    case 'gomp'
        A = QOI_decompress(N.QOI_result);
        ftNu = triu1D(A,N.width);
        ftN = ftNu + ftNu.' - diag(diag(ftNu));
        M = ifft2(ftN);
end
end

function M = QOI_decompress(N,debugstuff)
        data = N.comp_data;
        index_position = @(x) mod(3*real(double(x))+5*imag(double(x)),64)+1;
        old_values = uint8(zeros(2^(6),1));
        A = uint8(zeros(N.width*N.heigth,1));
        index = 1;
        Aindex = 1;
        while index <= length(data)
%             if Aindex == 5581
%                 disp('lmao')
%             end

            header = data(index);
            index = index + 1;
            switch bitshift(header,-6)
                case 0b10 %New values
                    runlength = bitand(header,0x3F);
                    for j = 0:runlength
                        if index+1 > length(data)
                            break
                        end
                        A(Aindex) = double(data(index)) + double(data(index+1))*1i;
                        old_values(index_position(A(Aindex))) = A(Aindex);
                        index = index + 2;
                        Aindex = Aindex + 1;
                    end
                case 0b11   %Repeating values
                    repVal = A(Aindex - 1);
                    replength = double(bitand(header,0x3F))+1;
                    A(Aindex:(Aindex+replength-1)) = repVal;
                    Aindex = Aindex + replength;
                case 0b01   %Close value
                    lastVal = double(A(Aindex - 1));
                    A(Aindex) = mod(real(lastVal) +...
                        double(bitshift(bitand(header,0b111000),-3)) - 4,256);
                    A(Aindex) = double(A(Aindex)) + 1i * mod(imag(lastVal) +...
                         double(bitand(header,0b111)) - 4,256);
                    old_values(index_position(A(Aindex))) = A(Aindex);
                    Aindex = Aindex + 1;
                case 0b00   %Old value
                    A(Aindex) = old_values(bitand(header,0x3F)+1);
                    Aindex = Aindex + 1;
            end 
        end
        temp = struct("matrix",reshape(A,[N.heigth,N.width]), ...
            "byteSize",1,"range",N.range);
        M = contin_mat(temp);
end


function plotDebugStuff(A,N,debugstuff,range)
temp = struct("matrix",reshape(A(1:(N.heigth*N.width)),[N.heigth,N.width]), ...
    "byteSize",1,"range",N.range);
M = contin_mat(temp);
subplot(1,3,1)
imagesc(abs(debugstuff))
title("Original")
subplot(1,3,2)
imagesc(abs(M))
title("Decompressed")
subplot(1,3,3)
imagesc(abs(debugstuff - M))
title("Diff")
if exist("range","var")
    if length(range) == 1
        point = find(A == 0,1);
        y = mod(point,N.heigth)+1;
        x = ceil(point/N.heigth);
        range = [x-5,x+5,y-5,y+5];
    end
    axis(range);
end
end