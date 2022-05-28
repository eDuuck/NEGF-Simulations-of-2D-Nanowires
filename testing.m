clear
sa = Sample(10,20,1,linspace(-0.7,-0.01,3));
sa.append(1.2*ones(3,20),'u');
sa.compress();
H = hamiltonian(sa,[0,1,1]);
imagesc(imag(H));
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
delete('testing.mat')

Ec = 1;
a = 1;
t0 = 0.05/a^2;
eps = Ec + 4*t0;
t =-t0;

sample = Sample(15,10,eps,t);
% sample.append(ones(sample.width,2)*eps*1.1);
% sample.append(ones(sample.width,5)*eps);
% sample.append(ones(sample.width,2)*eps*1.1);
% sample.append(ones(sample.width,5)*eps);
sample.addContact(ones(sample.width,1)*eps,t,[1,1]);
sample.addContact(ones(sample.width,1)*eps,t,[1,sample.length]);
sample.contacts{end}.fermi = 0;


sample.D = eye(sample.M)*1e-4;
%sample.applyNoise(0.01,3);

% res = NEGF(sample,1.1);
% electron_den = NEGF_result_remap(res,'electrons');
% fermi_levels = NEGF_result_remap(res,'fermi');
% figure(2);
% subplot(1,2,1)
% imagesc(electron_den)
% subplot(1,2,2)
% imagesc(fermi_levels)

sim_points = 100;
results = struct('NEGF_results',{},'complete',zeros(1,sim_points));

filename = 'testing.mat';
if isfile(filename)
    load(filename);
else
    results = struct('complete',zeros(1,sim_points));
    results.NEGF_results = cell(1,sim_points);
end
E = linspace(1,1.4,sim_points);
B = [0,0,1];
k = 0;
disp("0/" + sim_points)
calc_timewo = zeros(1,sim_points);
for x = 1:sim_points
    tic
    if ~results.complete(x)
        if x > 1
            G0 = 0;%results.NEGF_results{x-1}.G;
        else
            G0 = 0;
        end
        results.NEGF_results{x} = NEGF(sample,E(x),B,1e-6,0.8,100,false,G0);
        results.complete(x) = 1;
        k = k + toc;
        calc_timewo(x) = toc;
    end
    disp(x+"/"+sim_points)
   if k > 300 || x == sim_points
       disp("Saving results...");
       k = 0;
       %save(filename,'results');
   end
end
temp = zeros(1,sim_points-1);
for i = 2:sim_points
    temp(i-1) = sqrt(sum(abs(results.NEGF_results{i}.G-results.NEGF_results{i-1}.G).^2,'all'));
end

T = zeros(1,sim_points);
%load('testing.mat')
for x = 1:sim_points
    gamma1 = real(1i*(results.NEGF_results{x}.sigma{1} - results.NEGF_results{x}.sigma{1}'));
    gamma2 = real(1i*(results.NEGF_results{x}.sigma{2} - results.NEGF_results{x}.sigma{2}'));
    T(x) = real(trace(gamma1*results.NEGF_results{x}.G*gamma2*results.NEGF_results{x}.G'));
end
figure(2)
plot(temp/max(temp))
hold on
plot(calc_time/max(calc_time))
plot(T)
hold off
%%

if true
    for x = 1:sim_points
        electrons = NEGF_result_remap(results.NEGF_results{x},'electrons');
        fermi_levels = NEGF_result_remap(results.NEGF_results{x},'fermi');
        %disp(max(fermi_levels,[],'all'));
        subplot(1,2,1)
        imagesc(electrons,[0,100])
        subplot(1,2,2)
        imagesc(fermi_levels,[0,1])
        pause(0.1);
    end
end
figure(2)
plot(T,E)