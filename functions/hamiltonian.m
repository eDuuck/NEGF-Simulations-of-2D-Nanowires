function [H] = hamiltonian(sample, B)
%HAMILTONIAN Returns the hamiltonian of a given sample.
%   Returns the hamiltonian matrix for a sample, a magnetic field described
%   by the matrix B can be applied to the sample. If the sample has a
%   size of [M N] then B needs to be of size [M N 3] for the different
%   coordinates of the vector field matching every point of the sample's
%   units.

wid = sample.width;
len = sample.length;
if nargin < 2
    B = zeros(wid,len,3);
end
if length(B) == 3
    B_field = ones(wid,len,3);
    B_field(:,:,1) = ones(wid,len) * B(1);
    B_field(:,:,2) = ones(wid,len) * B(2);
    B_field(:,:,3) = ones(wid,len) * B(3);
else
    B_field = B;
end

A = anti_curl(B_field);
if lower(sample.arch) == "rectangular"
    
    t = sample.t;
    units = sample.getUnits();
%     H_elements_max = (wid+length(t)*(2*wid-length(t)-1))*len + ... %Main column elements.
%         wid*length(t)*(2*len-length(t)-1);
%     
%     H = spalloc(sample.M, sample.M,H_elements_max);
    %     for j = 1:x
    %         H(((y*(j-1))+1):(y*j),((y*(j-1))+1):(y*j)) = ...
    %                        column(y,units(:,j),t); %Fills column blocks.
    %     end
    %
    %
    %     for j = 1:length(t)
    %         H = H + sparse(t(j)*(diag(ones((x-j)*y,1),y*j)+diag(ones((x-j)*y,1),-y*j)));
    %     end
    H = zeros(sample.M);
    for x = 1:len
        for y = 1:wid
            pos = (x-1)*wid+y;
            H(pos,pos) = units(y,x);
            for j = 1:length(t)
                if y+j-1 < wid
                    H(pos,pos+j) = t(j) * exp(1i * A(y,x,2));
                end
                if y-(j-1) > 1
                    H(pos,pos-j) = t(j) * exp(-1i * A(y,x,2));
                end
            end
        end
        startpos = (x-1)*wid + 1;
        for j = 1:length(t)
            if x + j - 1 < len
                H(startpos:(startpos+wid-1),(startpos+wid*j):(startpos+(j+1)*wid-1)) ...
                    = diag(t(j)*exp(1i * A(:,x,1)));
            end
            if x - (j - 1) > 1
                H(startpos:(startpos+wid-1),(startpos-j*wid):(startpos-wid*(j-1)-1)) ...
                    = diag(t(j)*exp(-1i * A(:,x,1)));
            end
        end
    end
end
end

function A = column(y,eps,t)
    A = diag(eps);
    for j = 1:length(t)
        A = A + t(j)*(diag(ones(y-j,1),j) + diag(ones(y-j,1),-j));
    end
end

