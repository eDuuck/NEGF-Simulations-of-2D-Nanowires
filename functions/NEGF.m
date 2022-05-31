function [result] = NEGF(sample,E,B,errorMarg,rate,it_lim,reduce,G0)
%NEGF Summary of this function goes here
%   NEGF(sample,E,B,errorMarg,rate,it_lim,reduce,G0)

if nargin < 7
    reduce = false;
end
if nargin < 6
    it_lim = 100;
end
if nargin < 5
    rate = 0.5;
end
if nargin < 4
    errorMarg = 1e-6 * min(sample.getUnits,[],'all');
end
if nargin < 3
    B = 0;
end

result = NEGF_result(sample,E);
H = hamiltonian(sample,B);

[sigma,sigmaIn] = sigma_from_sample(sample,E);
sigSum = sparse(sample.M,sample.M);
sigInSum = sparse(sample.M,sample.M);
for j = 1:length(sigma)
    sigSum = sigSum + sigma{j};
    sigInSum = sigInSum + real(1i * (sigmaIn{j} - sigmaIn{j}'));
end
EI = speye(sample.M) * E;

if nargin == 8
    G = G0;
else
    G = (EI - H - sigSum)^-1;
end

D = sample.D;
if ~isequal(D,0)
    
    sigma0 = D .* G;
    
    for j = 1:it_lim
        G =(EI - H - sigSum - sigma0)^-1;
        Sigma0New = D .* G;
        change = Sigma0New - sigma0;
        if max(abs(change), [], 'all') < errorMarg
            
            %disp("sigma0 in " + j)
            break;
        end
        sigma0 = sigma0 + rate*change;
    end
    Gn = G*sigInSum*G';
    sigma0In = D .* Gn;
    for j = 1:it_lim
        Gn = G *(sigInSum + sigma0In) * G';
        Sigma0InNew = D .* Gn;
        change = Sigma0InNew - sigma0In;
        if max(abs(change), [], 'all') < errorMarg
            %disp("sigma0in in " + j)
            break;
        end
        sigma0In = sigma0In + rate*change;
    end
else
    G = (EI - H - sigSum)^-1;
    sigma0 = 0;
    sigma0In = 0;
    Gn = G*sigInSum*G';
end
if reduce
    %DoReduceStuff
end
result.G = G * 1;
result.Gn = Gn;
result.sigma = sigma;
result.sigmaIn = sigmaIn;
result.sigma0 = sigma0;
result.sigma0In = sigma0In;
result.reduced = reduce;
end