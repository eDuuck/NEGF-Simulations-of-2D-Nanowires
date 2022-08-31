clear all
load('physical_constants.mat');

Ec = 1;
a = 25E-10;
effective_mass = 0.1;
t0 =  h_bar^2/(2*electron_mass*effective_mass* a^2 * e_charge);
eps = 4*t0;
t =-t0;

sample = Sample(25,1,eps,t,a);
%sample.append(ones(sajmple.width,2)*eps*1.1);
%sample.append(ones(sample.width,5)*eps);
%sample.append(ones(sample.width,2)*eps*1.1);
%sample.append(ones(sample.width,5)*eps);
sample.addContact(ones(sample.width,1)*eps,t,[1,1]);
sample.addContact(ones(sample.width,1)*eps,t,[1,sample.length]);
sample.contacts{end}.fermi = 0;

%sample.D = 1e-5;%eye(sample.M)*1e-4;
%sample.applyNoise(0.01,3);


E = t0;
B = 0:0.1:50;%linspace(0,10,100);

tic
result = NEGF_map(sample,E,B,0,300,1e-6,0.5,50,false);
toc
disp(getSize(result,2))