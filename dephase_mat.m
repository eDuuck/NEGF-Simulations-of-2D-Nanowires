function [G, Gn, Sigma0, Sigma0In] = dephase_mat(E,H,Sigma,Fermi,D,G0,errorMarg,rep_lim, g)
%CALC_SIGMA_DEPHASE Summary of this function goes here
%   Detailed explanation goes here
if nargin < 9
    g = 0.6;
end
if nargin < 8
    rep_lim = 1E2;
end
if nargin < 7
    errorMarg = 1e-6;
end

    

EI = E*eye(size(H,1));

sigSum = sum(Sigma,3);
SigmaIn = zeros(size(Sigma));
for k = 1:length(Fermi)
    SigmaIn(:,:,k) = Fermi(k)*real(1i * (Sigma(:,:,k) - Sigma(:,:,k)'));
end

if nargin < 6
    G =(EI - H - sigSum)^-1;
else
    G = G0;
end

if ~isequal(D,0)
    
    Sigma0 = D .* G;
    
    for j = 1:rep_lim
        G =(EI - H - sigSum - Sigma0)^-1;
        Sigma0New = D .* G;
        if max(abs(Sigma0New - Sigma0), [], 'all') < errorMarg
            break;
        end
        Sigma0 = Sigma0 + g*(Sigma0New - Sigma0);
    end
    if j == rep_lim
        disp("lmaoooo")
    end
    Gn=G*sum(SigmaIn,3)*G';
    Sigma0In = D .* Gn;
    for j = 1:rep_lim
        Gn = G *(sum(SigmaIn,3) + Sigma0In) * G';
        Sigma0InNew = D .* Gn;
        if max(abs(Sigma0InNew - Sigma0In), [], 'all') < errorMarg
            break;
        end
        Sigma0In = Sigma0In + g*(Sigma0InNew - Sigma0In);
    end
else
    Sigma0 = 0;
    Sigma0In = 0;
    Gn=G*sum(SigmaIn,3)*G';
end