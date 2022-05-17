function [sigma,sigmaIn] = sigma_from_sample(sample,E)
%CONTACTSIGMA Summary of this function goes here
%   Detailed explanation goes here
sigma = cell(1,sample.nbrOfContacts);
sigmaIn = sigma;
for j = 1:sample.nbrOfContacts
    contact = sample.contacts{j};
    con_width = size(contact.con_mat,1);
    con_length = 1;
    y = contact.pos(1); x = contact.pos(2);
    SGF = contact_surface(contact,E);
    sigma{j} = spalloc(sample.M, sample.M,numel(SGF));
    for k = 1:con_length %This might be unnecessary convoluted.
        sIndx = (x+k-2)*sample.width + y;
        eIndx = sIndx + con_width-1;
        sigma{j}(sIndx:eIndx,sIndx:eIndx) = SGF;
    end
    sigmaIn{j} = sigma{j} * contact.fermi;
end

end

