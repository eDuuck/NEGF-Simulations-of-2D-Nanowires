function [H,V] = twoDimmRed(struct, t)
%TWODIMMRED Returns a set of one-dimensional hamiltonians that can be used
%to solve a two dimensional problem.
%   Detailed explanation goes here
[w,l] = size(struct);

H = zeros(l,l,w);

cH = zeros(w,w,l);
V = zeros(w,w,l);

cH(:,:,1) = diag(struct(:,1)) + t*diag(ones(w-1,1),1) + ...
    t*diag(ones(w-1,1),-1);
[V(:,:,1),D] = eig(cH(:,:,1));
H(1,1,:) = diag(D);

for j = 2:l
    %Could expand this check to see if the same column appears more than
    %once, significant improvement in calculation time in bigger matrices.
    
    %Currently only comparing current row to last one.
    if(isequal(struct(:,j),struct(:,j-1)) && false)
        H(j,j,:) = H(j-1,j-1,:);
        V(:,:,j) = V(:,:,j-1);
    else
        cH(:,:,j) = diag(struct(:,j)) + t*diag(ones(w-1,1),1) + ...
            t*diag(ones(w-1,1),-1);        %Create the column matrix.
        [V(:,:,j),D] = eig(cH(:,:,j));     %Diagonalize the matrix.
        H(j,j,:) = diag(D);                %Insert diagonal elements into H.
    end
end

for j = 1:w
    H(:,:,j) = H(:,:,j) + t*diag(ones(l-1,1),1) + t*diag(ones(l-1,1),-1);
end

function c