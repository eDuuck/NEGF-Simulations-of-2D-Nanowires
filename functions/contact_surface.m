function [SGF,g0] = contact_surface(contact,E,B,r0,rate,acc)
%CONTACT_SURFACE Returns the surface greens function for a given contact.
%   Detailed explanation goes here
if nargin < 6
    acc = min(abs(contact.SC),[],'all')*1e-4;
end
if nargin < 5
    rate = 0.8;
end

left = contact.face == 1; %%Since it's facing right, it's coming from the left side.



q = 1.6022e-19; hbar = 1.0546e-34; qha = 1i*contact.a*q/hbar;

%SGF = zeros(size(contact.SC,1));

A = anti_curl(B,contact.a);

t_l = length(contact.tau);
M = numel(contact.SC);
temp = Sample(M,t_l,contact.SC(1),contact.tau);
alpha = hamiltonian(temp,B);
beta = zeros(t_l*M);
for j = 1:t_l
    beta = beta + diag(ones(M*(t_l-j+1),1) .* (contact.tau(t_l-j+1).*exp(qha*A(:,1))) ...
        ,M*(j-1));
end
I = eye(size(contact.SC,1));
%imagesc(imag(beta))
%pause(0.01);

if contact.basis_length == Inf
    iterations = 100000;
else
    iterations = contact.basis_length;
    rate = 1;
end
SGC = (E+1i*contact.eta)*I-alpha;
if ~isequal(r0,0)
    SGF = r0;
else
    SGF = (SGC)^-1; 
end
for j = 1:iterations
    if j == 500
        rate = 0.5;
    end
    if left
        SGFnew = (SGC-beta'*SGF*beta)^-1;
    else
        SGFnew = (SGC-beta*SGF*beta')^-1;
    end
    change = SGFnew - SGF;
    c_val = sum(abs(change),"all")/sum(abs(SGFnew)+abs(SGF),"all");
    SGF = SGF + rate*change;
    if c_val < acc && (contact.basis_length == Inf)
        break
    end
end
g0 = SGF;
if left
    SGF = beta'* SGF * beta;
else
    SGF = beta* SGF * beta';
end
% figure(1)
% subplot(3,3,1 + left*3);imagesc(abs(SGF));subplot(3,3,2 + left*3);imagesc(real(SGF));subplot(3,3,3 + left*3);imagesc(imag(SGF));
% if left
% subplot(3,3,1 + 6);imagesc(abs(g0));subplot(3,3,2 + 6);imagesc(real(g0));subplot(3,3,3 + 6);imagesc(imag(g0));
% end
% pause(0.001);
SGF = sparse(SGF);
end


