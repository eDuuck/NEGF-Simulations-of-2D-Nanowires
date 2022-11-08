function [sigma,sigmaIn,g0_ret] = sigma_from_sample(NEGF_param)
%SIGMA_FROM_SAMPLE returns the sigma and sigmaIn matrices for a sample
%specified in the NEGF_param.
%   The returned matrices are placed in cells. Thus to access the Sigma
%   matrix for the first contact, use sigma{1} from the returned cell.

sample = NEGF_param.sample;
r0 = NEGF_param.g0; %g0 is an initial guess for the contacts.

sigma = cell(1,sample.nbrOfContacts);
sigmaIn = cell(1,sample.nbrOfContacts);
g0_ret = cell(1,sample.nbrOfContacts);

for j = 1:sample.nbrOfContacts
    contact = sample.contacts{j};

    if r0 ~= 0
        g0 = r0.getS0(j);
    else
        g0 = 0;
    end

    [SGF,g0_ret{j}] = contact_surface(contact,NEGF_param,g0);
    
    %Place SGF in correct place of the hamiltonian.
    Hpos = zeros(sample.length);
    Hpos(contact.pos(2),contact.pos(2)) = 1;
    sigma{j} = kron(Hpos,SGF);

    sigmaIn{j} = 1i*(sigma{j}-sigma{j}') * contact.fermi;
end

end

