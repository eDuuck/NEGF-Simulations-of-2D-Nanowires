function [H] = hamiltonian(sample, B)
%HAMILTONIAN Returns the hamiltonian of a given sample.
%   Returns the hamiltonian matrix for a sample, a magnetic field described
%   by the matrix B can be applied to the sample. If the sample has a
%   size of [M N] then B needs to be of size [M N 3] for the different
%   coordinates of the vector field matching every point of the sample's
%   units.
if lower(sample.arch) == "honeycomb"
    N = 2*(prod(sample.dim)+sum(sample.dim));
    H = zeros(N);
    
    if size(sample.eps) == 1
        H = eye(N)*sample.eps;
    end
elseif lower(sample.arch) == "rectangular"
    y = sample.width;
    x = sample.length;
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

