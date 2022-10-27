function [SGF,g0] = contact_surface(contact,NEGF_param,r0)
%CONTACT_SURFACE Returns the surface greens function for a given contact.
%   Detailed explanation goes here
E = NEGF_param.E(1);
B = ones(size(contact.SC))*NEGF_param.B(1);

scaling = 1/(E+contact.SC(1));
Es = E*scaling;
eta = contact.eta*scaling;

acc = NEGF_param.errorMarg;
if nargin < 5
    rate = 0.8;
end

left = contact.face == 1; %%Since it's facing right, it's coming from the left side.

if length(B) == 1
    B = ones(size(contact.SC))*B;
end

q = 1.6022e-19; hbar = 1.0546e-34; qha = 1i*contact.a*q/hbar;

A = anti_curl(B,contact.a);

t_l = length(contact.tau);
M = numel(contact.SC);
temp = Sample(M,t_l,contact.SC(1), contact.tau);  %Currently only uniform contacts are possible.
alpha = hamiltonian(temp,B);
beta = zeros(t_l*M);

for j = 1:t_l
    beta = beta + diag(ones(M*(t_l-j+1),1) .* (contact.tau(t_l-j+1).*exp(qha*A(:,1))) ...
        ,M*(j-1));
end

alpha = alpha * scaling;
beta = beta * scaling;

I = eye(size(contact.SC,1));

if contact.basis_length == Inf
    iterations = NEGF_param.con_it_lim;
else
    iterations = contact.basis_length;
    rate = 1;
end
SGC = (Es+1i*eta)*I-alpha;
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
        %SGFnew = pinv(SGC-beta'*SGF*beta);
    else
        SGFnew = (SGC-beta*SGF*beta')^-1;
        %SGFnew = pinv(SGC-beta*SGF*beta');
    end
    change = SGFnew - SGF;
    c_val = sqrt(sum(real(change).^2 + imag(change).^2,"all")/...
        (sum(real(SGF).^2 + imag(SGF).^2,"all")));
    SGF = SGF + rate*change;
    if c_val < acc && (contact.basis_length == Inf)
        break
    end
end
if j == iterations && (contact.basis_length ~= Inf)
    if NEGF_param.error_halt
        error("Did not converge in " + iterations + " iterations. C_val = "+ c_val);
    else
        disp("Did not converge in " + iterations + " iterations. C_val = "+ c_val);
    end
end
g0 = SGF;
if left
    SGF = beta'* SGF * beta;
else
    SGF = beta* SGF * beta';
end
SGF = sparse(SGF/scaling);
end


