function [SGF] = contact_surface(contact,E,rate,acc)
%CONTACT_SURFACE Returns the surface greens function for a given contact.
%   Detailed explanation goes here
if nargin < 4
    acc = min(abs(contact.SC),[],'all')*1e-3;
end
if nargin < 3
    rate = 0.5;
end

SGF = zeros(size(contact.alpha,1));
alpha = contact.alpha;
beta = contact.beta;
I = eye(size(contact.alpha,1));

if contact.basis_length == Inf
    iterations = 1000;
else
    iterations = contact.basis_length;
    rate = 1;
end
max_change = zeros(iterations,1);
for j = 1:iterations
    SGFnew = (E*I+1i*contact.eta*I-alpha-beta'*SGF*beta)^-1;
    change = SGFnew - SGF;
    SGF = SGF + rate*change;
    %max_change(j) = max(abs(change), [],'all');
    if (max(abs(change), [],'all') < acc) && (contact.basis_length == Inf)
        %plot(max_change);
        %hold on
        %plot(j,0,'rx');
        %hold off
        %disp(j)
        break
    end
end

SGF = contact.con_mat * SGF * contact.con_mat';
%SGF = sparse(SGF);
end


