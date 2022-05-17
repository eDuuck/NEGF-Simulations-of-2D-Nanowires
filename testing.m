clear
sa = Sample(10,20,1,linspace(-0.7,-0.01,3));
sa.append(1.2*ones(3,20),'u');
sa.compress();
H = hamiltonian(sa);
%imagesc(abs(H));
sa.addContact(1.2*ones(10,1),-0.05,[1,1]);
con_mat = [eye(4),zeros(4,6);zeros(1,10);
           0, 0, 0, 0, 0, 0, 1, 0, 0, 0];
sa.addContact(1.2*ones(10,1),-0.05,[1,1],con_mat)
con = sa.contacts{1};
SGF = contact_surface(con,1.1,0.5);
con = sa.contacts{2};
SGF2 = contact_surface(con,1.1,0.5);
subplot(1,2,1);
imagesc(abs(SGF));
subplot(1,2,2);
imagesc(abs(SGF2));

%%
clear
for E = 0.9:0.01:1.5
    tic
tau = -0.05;
eta = 0.001;
Ec = 1;
eps = Ec - 4*tau;
con = Contact(ones(10,1)*eps,tau,[1,1]);
SGF = contact_surface(con,E,0.5);
%contact_testing
figure(2)
imagesc(abs(SGF));
title(E);
xlabel(toc)
pause(0.05)
end

%%
clear
sa = Sample(10,20,1,linspace(-0.7,-0.01,3));
sa.append(1.2*ones(3,20),'u');
sa.compress();
H = hamiltonian(sa);
%imagesc(abs(H));
sa.addContact(1.2*ones(10,1),-0.05,[1,1]);
con_mat = [eye(4),zeros(4,6);zeros(1,10);
           0, 0, 0, 0, 0, 0, 1, 0, 0, 0];
sa.addContact(1.2*ones(10,1),-0.05,[1,5],con_mat,1)

[sigma,sigmaIn] = sigma_from_sample(sa,1.1);
imagesc(abs(sigmaIn{2}))




%% Test if NEGF is working correctly, compare to results from uniform2Dwire.

Ec = 1;
a = 1;
t0 = 0.05/a^2;
eps = Ec + 4*t0;
t = -t0;

sample = Sample(15,10,eps,t);
sample.addContact(ones(sample.width,1)*eps,t,[1,1]);
sample.addContact(ones(sample.width,1)*eps,t,[1,sample.length]);
sample.contacts{end}.fermi = 0;

res = NEGF(sample,1.1);
electron_den = electron_density(res);
figure(2);
imagesc(electron_den)
