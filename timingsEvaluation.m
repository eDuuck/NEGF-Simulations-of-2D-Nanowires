clear all
Ec = 1;
a = 1;
t0 = 0.05/a^2;
eps = Ec + 4*t0;
t =-t0;

sample = Sample(15,10,eps,t);
sample.append(ones(sample.width,2)*eps*1.1);
sample.append(ones(sample.width,5)*eps);
sample.append(ones(sample.width,2)*eps*1.1);
sample.append(ones(sample.width,5)*eps);
sample.addContact(ones(sample.width,1)*eps,t,[1,1]);
sample.addContact(ones(sample.width,1)*eps,t,[1,sample.length]);
sample.contacts{end}.fermi = 0;

sample.D = 1e-4;
sample.applyNoise(0.01,3);

E = 1.2;
B = 0;

NEGF_result = NEGF(sample,E,B,1e-6,0.8,100);
subplot(1,3,1)
imagesc(abs(NEGF_result.G)); 
orgG = NEGF_result.G;
NEGF_result.reduce();

NEGF_result.reduce(0);
subplot(1,3,2)
imagesc(abs(NEGF_result.G));
subplot(1,3,3)
imagesc(abs(NEGF_result.G-orgG));
disp(mean(abs(NEGF_result.G-orgG),"all"));