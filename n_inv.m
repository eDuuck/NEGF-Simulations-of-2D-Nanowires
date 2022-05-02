function [X_inv] = n_inv(X, G, error, max_it,rate)
%N_INV Summary of this function goes here
%   Detailed explanation goes here
if nargin < 5
    rate = 0.1;
end
if nargin < 4
    max_it = 15;
end
if nargin < 3
    error = 1e-3;
end
if nargin < 2 || min(G == 0,[],'all')
    X_inv = X;
else
    X_inv = G;
end
for j = 1:max_it
    X_new = 2*X_inv - X_inv*X*X_inv;
    change = X_new - X_inv;
    X_inv = X_inv + rate*change;
    if norm(change,2) < error
        break;
    end
end
end

