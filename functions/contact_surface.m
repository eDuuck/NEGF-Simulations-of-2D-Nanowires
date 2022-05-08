function [SGF] = contact_surface(contact,E,acc,rate)
%CONTACT_SURFACE Returns the surface greens function for a given contact.
%   Detailed explanation goes here
if nargin < 3
    acc = mean(contact.SC)*1e-3;
end
if nargin < 4
    rate = 0.75;
end

SGF = zeros(numel(contact.SC));
alpha = contact.alpha;
beta = contact.beta;
I = eye(numel(contact.SC));

if contact.size == Inf
    iterations = 1000;
else
    iterations = contact.size;
    rate = 1;
end

for j = 1:iterations
    SGFnew = (E*I+1i*contact.eta*I-alpha-beta'*SGF*beta)^-1;
    change = SGFnew - SGF;
    %max_change(j) = max(abs(change),[],'all');
    SGF = SGF + rate*change;
    if (max(abs(change), [],'all') < acc) && (contact.size == Inf)
        break
    end
end
%SGF = sparse(SGF);
end


