function [A] = anti_curl(B,a,c)
%ANTI_CURL A very simple function that returns the magnetic vector 
%potential to the magnetic field B. Note that B should only contain the Z
%component. The value a specifies the seperation between lattice points in
%B. C can be used to specify the constant dAydx.
[width, length] = size(B);
A = zeros(width, length,2);

if isequal(0*B,B)
    return;
end
if nargin < 3
    c = 0;
end

dAydx = ones(width, length)*c;
dAxdy = dAydx - B*a;
A(:,:,1) = -tril(ones(width))*dAxdy;
A(:,:,2) = dAydx*triu(ones(length));
A(:,:,1) = A(:,:,1) - mean(A(:,:,1),'all');
A(:,:,2) = A(:,:,2) - mean(A(:,:,2),'all');
end

