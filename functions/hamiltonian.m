function [H] = hamiltonian(sample, B)
%HAMILTONIAN Returns the hamiltonian of a given sample.
%   Returns the hamiltonian matrix for a sample, a magnetic field described
%   by the matrix B can be applied to the sample. If the sample has a
%   size of [M N] then B needs to be of size [M N 3] for the different
%   coordinates of the vector field matching every point of the sample's
%   units.

y = sample.width;
x = sample.length;
if nargin < 2
    B = 0;
end
if length(B) == 3
    B_field = ones(y,x,3);
    B_field(:,:,1) = ones(y,x) * B(1);
    B_field(:,:,2) = ones(y,x) * B(2);
    B_field(:,:,3) = ones(y,x) * B(3);
else
    B_field = B;
end

if lower(sample.arch) == "rectangular"

    t = sample.t;
    units = sample.getUnits();
    H_elements_max = (y+length(t)*(2*y-length(t)-1))*x + ... %Main column elements.
                    y*length(t)*(2*x-length(t)-1);
    
    H = spalloc(sample.M, sample.M,H_elements_max);
    for j = 1:x
        H(((y*(j-1))+1):(y*j),((y*(j-1))+1):(y*j)) = ...
                       column(y,units(:,j),t); %Fills column blocks.
    end
   
    
    for j = 1:length(t)
        H = H + sparse(t(j)*(diag(ones((x-j)*y,1),y*j)+diag(ones((x-j)*y,1),-y*j)));
    end
end
end

function A = column(y,eps,t)
    A = diag(eps);
    for j = 1:length(t)
        A = A + t(j)*(diag(ones(y-j,1),j) + diag(ones(y-j,1),-j));
    end
end

