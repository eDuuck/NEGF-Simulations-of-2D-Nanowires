function [SGF] = contact_surface(contact,E, rate, iterations, acc)
%CONTACT_SURFACE Returns the surface greens function for a given contact.
%   Detailed explanation goes here
if nargin < 4
    acc = mean(contact)*1e-3;
end
if nargin < 3
    iterations = 500;
end
if nargin < 2
    rate = 0.75;
end

SGF = zeros(numel(contact.SC));
alpha = contact.alpha;
beta = contact.beta;


for j = 1:iterations
    SGFnew = (E*I+1i*contact.eta*I-alpha-beta'*SGF*beta)^-1;
    change = SGFnew - SGF;
    %max_change(j) = max(abs(change),[],'all');
    gn = gn + rate(k)*change;
    if max(abs(change)) < acc
        break
    end
end

