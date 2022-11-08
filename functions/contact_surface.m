function [SGF,g0] = contact_surface(contact,NEGF_param,r0)
%CONTACT_SURFACE Returns the surface greens function for a given contact.
%   This is done iteratively with parameters from NEGF_param. r0 is an
%   initial guess for the given contacts green function. Without this guess
%   certain setups converge very slowely and thus it may be faster to
%   simulate a range with nearby values first.



%First import parameters from contact and NEGF_param.
E = NEGF_param.E(1);
B = ones(size(contact.SC))*NEGF_param.B(1);

scaling = 1/(E+contact.SC(1));%Not sure if scaling is needed.
E = E*scaling;
eta = contact.eta*scaling;  
acc = NEGF_param.errorMarg;
rate = NEGF_param.rate;
if contact.unit_length == Inf
    iterations = NEGF_param.con_it_lim;
else
    iterations = contact.unit_length;
    rate = 1;
end

if length(B) == 1 %This assumes a uniform magnetic field.
    B = ones(size(contact.SC))*B;
end


%Next step is to form the unit cells that are used to calculate SGF.
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


SGC = (E+1i*eta)*I-alpha; %Constant during the algorithm.
if ~isequal(r0,0)
    SGF = r0;
else
    SGF = (SGC)^-1; 
end
%Begin the converging algorithm.
for j = 1:iterations
%     if j == 500
%         rate = 0.5;
%     end
    if contact.face == 1
        SGFnew = (SGC-beta'*SGF*beta)^-1;
    else
        SGFnew = (SGC-beta*SGF*beta')^-1;
    end
    change = SGFnew - SGF;
    c_val = sqrt(sum(real(change).^2 + imag(change).^2,"all")/...
        (sum(real(SGF).^2 + imag(SGF).^2,"all")));
    SGF = SGF + rate*change;
    if c_val < acc && (contact.unit_length == Inf)
        break
    end
end

%If it didn't converge in iteration limit.
if j == iterations && (contact.unit_length ~= Inf)
    if NEGF_param.error_halt
        error("Did not converge in " + iterations + " iterations. C_val = "+ c_val);
    else
        disp("Did not converge in " + iterations + " iterations. C_val = "+ c_val);
    end
end
g0 = SGF;

%This final steps connects the SGF to the channel.
if contact.face == 1
    SGF = beta'* SGF * beta;
else
    SGF = beta* SGF * beta';
end
SGF = sparse(SGF/scaling);
end


