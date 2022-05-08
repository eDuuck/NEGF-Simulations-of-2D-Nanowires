clear
sa = Sample(20,100,1,linspace(-1,-0.01,3));
H = hamiltonian(sa);
%imagesc(abs(H));
sa.addContact(1.2*ones(5,1),-0.05)
%%
clear
for E = 0.9:0.01:2
    tic
tau = -0.05;
eta = 0.001;
%E = 1.01;
Ec = 1;
eps = Ec - 4*tau;
con = Contact(ones(9,1)*eps,tau,Inf);
con.eta = eta;
SGF = contact_surface(con,E,1e-10,0.7);
%contact_testing
figure(2)
imagesc(abs(SGF),[0 15]);
title(E);
xlabel(toc)
pause(0.05)
end

