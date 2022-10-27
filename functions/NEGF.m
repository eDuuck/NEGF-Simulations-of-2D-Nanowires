function [result] = NEGF(NEGF_param)
%NEGF Summary of this function goes here
%   NEGF(sample,E,B,errorMarg,rate,it_lim,reduce,G0)
E = NEGF_param.E(1);
B = NEGF_param.B(1);
sample = NEGF_param.sample;

it_lim = NEGF_param.it_lim;
rate = NEGF_param.rate;
errorMarg = NEGF_param.errorMarg;


result = NEGF_result(sample,E,B);
H = hamiltonian(sample,B);

%Not sure if scaling is needed.
%scaling = 1/(rms(E)+rms(B)+min(abs(sample.units),[],'all') + min(abs(sample.t)));
scaling = 1;

[sigma,sigmaIn,s0] = sigma_from_sample(NEGF_param);
sigSum = sparse(sample.M,sample.M);
sigInSum = sparse(sample.M,sample.M);
for j = 1:length(sigma)
    sigSum = sigSum + sigma{j};
    sigInSum = sigInSum + sigmaIn{j};
end
H = H*scaling; E = E*scaling; B = B*scaling;
sigSum = sigSum * scaling; sigInSum = sigInSum * scaling;

EI = eye(sample.M) * E;


D = sample.D * scaling^2;
if ~isequal(D,D*0)
    if NEGF_param.g0 ~= 0
        sigma0 = NEGF_param.g0.getSigma0() * scaling;
    else
        G = (EI - H - sigSum)^-1;
        sigma0 = D .* G;
    end
    
    for j = 1:it_lim
        G =(EI - H - sigSum - sigma0)^-1;
        sigma0New = D .* G;
        change = sigma0New - sigma0;
        c_val = sqrt(sum(real(change).^2 + imag(change).^2,"all")/...
        (sum(real(sigma0).^2 + imag(sigma0).^2,"all")));
        if c_val < errorMarg
            break;
        end
        sigma0 = sigma0 + rate*change;
    end
    if j == it_lim && NEGF_param.error_halt
        error("Sigma0 did not converge in " + it_lim + ...
            " iterations. Error margin = " + c_val);
    end

    if NEGF_param.g0 ~= 0
        sigma0In = NEGF_param.g0.getSigma0In() * scaling;
    else
        Gn = G*sigInSum*G';
        sigma0In = D .* Gn;
    end

    for j = 1:it_lim
        Gn = G *(sigInSum + sigma0In) * G';
        Sigma0InNew = D .* Gn;
        change = Sigma0InNew - sigma0In;
        c_val = sqrt(sum(real(change).^2 + imag(change).^2,"all")/...
        (sum(real(sigma0In).^2 + imag(sigma0In).^2,"all")));
        if c_val < errorMarg
            break;
        end
        sigma0In = sigma0In + rate*change;
    end

    if j == it_lim && NEGF_param.error_halt
        error("Sigma0In did not converge in " + it_lim + ...
            " iterations. Error margin = " + c_val);
    end
else
    %G = (EI - H - sigSum)^-1;
    sigma0 = 0;
    sigma0In = 0;
end
fermiLevels = zeros(sample.nbrOfContacts,1);
for j = 1:sample.nbrOfContacts
    fermiLevels(j) = sample.contacts{j}.fermi;
end
%result.G = G;  %This isn't needed since G is fast to calculate once sigma0
%and surface green functions are extracted.
result.sigma = sigma;
%result.sigmaIn = sigmaIn; (Same as G)
result.sigma0 = sigma0/scaling;
result.sigma0In = sigma0In/scaling;
result.fermiLevels = fermiLevels;
result.s0 = s0;
if NEGF_param.compress
    if B == 0 %When no magnetic fields are present, matrices are symmetric. (It seems)
        result.compressionMethod = 'gomp'; 
    end
    result.compress(true);
end
end