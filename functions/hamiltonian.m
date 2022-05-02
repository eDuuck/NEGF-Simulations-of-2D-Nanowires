function [H] = hamiltonian(sample)
%HAMILTONIAN Summary of this function goes here
%   Detailed explanation goes here
if lower(sample.arch) == "honeycomb"
    N = 2*(prod(sample.dim)+sum(sample.dim));
    H = zeros(N);
    
    if size(sample.eps) == 1
        H = eye(N)*sample.eps;
    end
elseif lower(sample.arch) == "rectangular"
    width = sample.width;
    length = sample.length;
    t = sample.conn;
    
    alpha = @(i) diag(ones(width,1)).*diag(sample.units(:,i)) + ...
    t*(diag(ones(width-1,1),1) + diag(ones(width-1,1),-1));
    
    beta = diag(ones(width,1))*t;
    
    H = sparse(sample.M, sample.M);
	H(1:width,1:width) = alpha(1);
    for j = 1:length-1
        H(((width*j)+1):(width*(j+1)),((width*j)+1):(width*(j+1))) = alpha(j+1);
        H(((width*j)+1):(width*(j+1)),((width*(j-1))+1):(width*j)) = beta;
        H(((width*(j-1))+1):(width*(j)),((width*j)+1):(width*(j+1))) = beta;
    end
end
end

