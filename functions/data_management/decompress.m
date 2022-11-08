function [M] = decompress(N,method)
%DECOMPRESS returns the full size matrix from a compressed datastructure
%returned from compress.m.
%   DECOMPRESS(N) uses the data structure N from compress.m to decompress
%   the data and returns a matrix M.
%
%   DECOMPRESS(N,method) only needs to be used if the method isn't
%   specified in the data structure N. This was used in earlier
%   implementations when N wasn't a structure but simply a matrix. This
%   shouldn't be needed to be used.


if ~exist("method","var")
    method = {N.method};
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
        M = QOI_decompress(N);
    case 'gomp'
        A = triu1D(QOI_decompress(N.QOI_result),N.width);
        M = A + A.' - diag(diag(A));
    case '8-bit'
        M = contin_mat(N);
end
end

function M = QOI_decompress(N)
        data = double(N.comp_data);
        old_values = zeros(2^(6),1);
        A = zeros(N.width*N.heigth,1);
        index = 1;
        Aindex = 1;
        while index <= length(data)
            header = data(index);
            index = index + 1;
            payLoad = mod(header,2^6);
            switch floor(header/2^6)
                case 2 %0b10 New values
                    runlength = payLoad;
                    for j = 0:runlength
                        if index+1 > length(data)
                            break
                        end
                        curData = complex(data(index),data(index+1));
                        A(Aindex) = curData;
                        indPos = mod(3*real(curData)+5*imag(curData),64)+1;
                        old_values(indPos) = curData;
                        index = index + 2;
                        Aindex = Aindex + 1;
                    end
                case 3  %0b11 Repeating values
                    repVal = A(Aindex - 1);
                    replength = payLoad;
                    A(Aindex:(Aindex+replength)) = repVal;
                    Aindex = Aindex + replength+1;
                case 1  %0b01 Close value
                    lastVal = A(Aindex - 1);
                    A(Aindex) = mod(real(lastVal) +...
                        floor(payLoad/2^3) - 4,256);
                    A(Aindex) = A(Aindex) + 1i * mod(imag(lastVal) +...
                         mod(header,2^3) - 4,256);
                    indPos = mod(3*real(A(Aindex))+5*imag(A(Aindex)),64)+1;
                    old_values(indPos) = A(Aindex);
                    Aindex = Aindex + 1;
                case 0   %Old value
                    A(Aindex) = old_values(payLoad+1);
                    Aindex = Aindex + 1;
            end 
        end
        temp = struct("matrix",reshape(A,[N.heigth,N.width]), ...
            "byteSize",1,"range",N.range,"method",N.discMethod);
        M = contin_mat(temp);
end



%This was only used for debugging. I'll keep it in if tempering with the
%algorithms needs to be done.
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