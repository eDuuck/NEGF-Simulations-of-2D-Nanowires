function [M] = decompress(N,method)
%DECOMPRESS Summary of this function goes here
%   Detailed explanation goes here
switch method{1}
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
    case 'QOI'
        data = N.comp_data;
        index_position = @(x) mod(3*real(double(x))+5*imag(double(x)),64)+1;
        old_values = uint8(zeros(2^(6),1));
        A = uint8(zeros(N.width*N.length));

        index = 1;
        Aindex = 1;
        while index <= length(data)
            header = data(index);
            switch bitshift(header,-6)
                case 0b10 %New values
                    index = index + 1;
                    for j = 1:bitand(header,0x3F)
                        A(Aindex) = data(index);
                        index = index + 1;
                        A(Aindex) = A(Aindex) + data(index)*1i;
                        old_values(index_position(A(Aindex))) = A(Aindex);
                        index = index + 1;
                        Aindex = Aindex + 1;
                    end
                case 0b11

                case 0b01

                case 0b00
            end
        end
end

