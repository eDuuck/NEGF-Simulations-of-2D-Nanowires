%% Generate data for the 1D-wire situation.
% Initial setup of sample.
clear
load('physical_constants.mat');

Ec = 1.12* eV;
a = 2E-10;
effective_mass = 1;
t0 =  h_bar^2/(2*electron_mass*effective_mass* a^2);
t =-t0;
eps = Ec - 2*t;

wid = 1; len = 20;  %Width 1; Length 40.

sample = Sample(wid,len,eps,t,a);
sample.addContact(ones(sample.width,1)*eps,t,[1,1]);
sample.addContact(ones(sample.width,1)*eps,t,[1,sample.length]);
sample.contacts{end}.fermi = 0;

E0 = Ec + 0 * t;
E_end = Ec - 4 * t;
sim_points = 100;

E = linspace(E0,E_end,sim_points);
B = 0;

tic

%----------------Plotting the simulations for a 1D uniform wire------------
figure(1);clf
hold on
D_vals = eV*eps*10.^[-inf,-3,-2];
for j = 1:length(D_vals)
    sample.D = D_vals(j) * diag(ones(1,sample.M));
    legends{j} = "D = " + D_vals(j);
    
    params = NEGF_param(sample,E,B);
    params.error_halt = false;
    result = NEGF_map(params);
    T = zeros(1,length(E));
    for i = 1:length(T)
        T(i) = NEGF_transmission(result.NEGF_result{i});
    end

    plot(E/eV,T,'Color',[j/length(D_vals),0,1-j/length(D_vals)])
end
hold off
xlabel("E [eV]"); ylabel("Transmission");
legend(legends); title("1D wire with no barrier.")



%----------------Adding double barrier to sample---------------------------
%sample.units(floor(sample.length/3)) = eps + 1*eV;
sample.units(floor(2*sample.length/4)) = eps + 1*eV;
D_vals = eV*eps*10.^[-inf,-3,-2,-1];

figure(2);clf
legends = cell(1,length(D_vals));
subplot(1,2,1)
hold on
for j = 1:length(D_vals)
    sample.D = D_vals(j) * diag(ones(1,sample.M));  %Phase + momentum relaxation.
    legends{j} = "D = " + D_vals(j);

    params = NEGF_param(sample,E,B);
    params.error_halt = false;
    result = NEGF_map(params);
    T = zeros(1,length(E));
    for i = 1:length(T)
        T(i) = NEGF_transmission(result.NEGF_result{i});
    end

    plot(E/eV,T,'Color',[j/length(D_vals),0,1-j/length(D_vals)])
end
hold off
axis([E(1)/eV, E(end)/eV, 0, 1]);
xlabel("E [eV]"); ylabel("Transmission");
legend(legends); title("1D wire with double barrier.")

legends = cell(1,length(D_vals));
subplot(1,2,2)
hold on
for j = 1:length(D_vals)
    sample.D = D_vals(j) *ones(sample.M);   %Only phase relaxation.
    legends{j} = "D = " + D_vals(j);

    params = NEGF_param(sample,E,B);
    params.error_halt = false;
    result = NEGF_map(params);
    T = zeros(1,length(E));
    for i = 1:length(T)
        T(i) = NEGF_transmission(result.NEGF_result{i});
    end

    plot(E/eV,T,'Color',[j/length(D_vals),0,1-j/length(D_vals)])
end
hold off
axis([E(1)/eV, E(end)/eV, 0, 1]);
xlabel("E [eV]"); ylabel("Transmission");
legend(legends); title("1D wire with double barrier.")

toc
%% Generate data for the 2D-wire situation.
% Initial setup of sample.
clear
load('physical_constants.mat');

Ec = 1.12* eV;
a = 2E-10;
effective_mass = 0.1;
t0 =  h_bar^2/(2*electron_mass*effective_mass* a^2);
t =-t0;
eps = Ec - 4*t;

wid = 10; len = 5;  %Width 1; Length 40.
sample = Sample(wid,len,eps,t,a);
sample.addContact(ones(sample.width,1)*eps,t,[1,1]);
sample.addContact(ones(sample.width,1)*eps,t,[1,sample.length]);
sample.contacts{end}.fermi = 0;

E0 = Ec;
E_end = Ec - 4 * t;
sim_points = 200;

E = linspace(E0,E_end,sim_points);
B = 0;
D_vals = eV*eps*10.^[-inf,-1];

figure(3);clf
legends = cell(1,length(D_vals));
hold on
for j = 1:length(D_vals)
    sample.D = D_vals(j) * diag(ones(1,sample.M));
    legends{j} = "D = " + D_vals(j);

    params = NEGF_param(sample,E,B);
    params.it_lim = 1000;
    params.error_halt = false;
    params.print = true;

    result = NEGF_map(params);
    T = zeros(1,length(E));
    for i = 1:length(T)
        T(i) = NEGF_transmission(result.NEGF_result{i});
    end

    plot(E/eV,T,'Color',[j/length(D_vals),0,1-j/length(D_vals)])
end
hold off
xlabel("E [eV]"); ylabel("Transmission");
legend(legends); title("2D wire with double barrier.")