function [B] = contin_mat(discrete_rep)
%CON_MAT Summary of this function goes here
%   Detailed explanation goes here
A = discrete_rep.matrix;
range = discrete_rep.range;
switch(lower(discrete_rep.method))
    case "linear"
        switch discrete_rep.byteSize
            case 1
                B =double(real(A))*diff(range(1,:))/255 + range(1,1) + ...
                    1i*(double(imag(A))*diff(range(2,:))/255 + range(2,1));
            case 2
                B =double(real(A))*diff(range(1,:))/(2^16-1) + range(1,1) + ...
                    1i*(double(imag(A))*diff(range(2,:))/(2^16-1) + range(2,1));
            case 4
                B =double(real(A))*diff(range(1,:))/(2^32-1) + range(1,1) + ...
                    1i*(double(imag(A))*diff(range(2,:))/(2^32-1) + range(2,1));
            case 8
                B =double(real(A))*diff(range(1,:))/(2^64-1) + range(1,1) + ...
                    1i*(double(imag(A))*diff(range(2,:))/(2^64-1) + range(2,1));
        end
    case "non-linear"
        B = reshape(range(1,real(A)) + 1i*range(2,imag(A)),size(A));
end
end

