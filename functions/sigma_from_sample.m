function [sigma,sigmaIn,g0_ret] = sigma_from_sample(NEGF_param)
%SIGMA_FROM_SAMPLE Summary of this function goes here
%   Detailed explanation goes here
sample = NEGF_param.sample;
r0 = NEGF_param.g0;

sigma = cell(1,sample.nbrOfContacts);
sigmaIn = sigma;
g0_ret = cell(sample.nbrOfContacts,1);
for j = 1:sample.nbrOfContacts
    contact = sample.contacts{j};
    con_width = size(contact.con_mat,1);
    con_length = 1;
    y = contact.pos(1); x = contact.pos(2);
    if r0 ~= 0
        g0 = r0.getS0{j};
    else
        g0 = 0;
    end
    
    [SGF,g0_ret{j}] = contact_surface(contact,NEGF_param,g0);
    
    Hpos = zeros(sample.length);
    Hpos(contact.pos(2),contact.pos(2)) = 1;
    sigma{j} = kron(Hpos,SGF);

%     sigma{j} = spalloc(sample.M, sample.M,numel(SGF));
%     for k = 1:con_length %This might be unnecessary convoluted. (Edit: It is.)
%         sIndx = (x+k-2)*sample.width + y;
%         eIndx = sIndx + con_width-1;
%         sigma{j}(sIndx:eIndx,sIndx:eIndx) = SGF;
%     end
    sigmaIn{j} = 1i*(sigma{j}-sigma{j}') * contact.fermi;
end

end

