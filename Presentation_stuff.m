%% 1D-tråd med dubbelbarriär
clear
load('physical_constants.mat');

Ec = 0* eV;
a = 5E-10;
effective_mass = 0.3;
t0 =  h_bar^2/(2*electron_mass*effective_mass* a^2);
t =-t0;
eps = Ec - 2*t;

E_start = Ec - 0.005 * t;
E_end = Ec - 2 * t;
sim_points = 200;

E = linspace(E_start,E_end,sim_points);
B = 0;

wid = 1; barr_len = 1;
barr_sep = 8; marg_len = 10;

sample = Sample(wid,marg_len,eps,t,a);
sample.append(ones(wid,barr_len)*eps*2);
sample.append(ones(wid,barr_sep)*eps);
sample.append(ones(wid,barr_len)*eps*2);
sample.append(ones(wid,marg_len)*eps);
%sample.append(ones(2,sample.length)*eps,'d')
sample.addContact(ones(sample.width,1)*eps,t,[1,1]);
sample.addContact(ones(sample.width,1)*eps,t,[1,sample.length]);
sample.contacts{end}.fermi = 0;

%sample.D = diag(ones(sample.M))*eps*eV*1E-2;
%sample.D = ones(sample.M)*eps*eV*5E-3;
%sample.D = diag(ones(sample.M))*eps*eV*1E-1;



params = NEGF_param(sample,E,B);
params.print = true;
params.error_halt = false;
params.rate = 0.90;
params.it_lim = 5000;
result = NEGF_map(params);

T = zeros(1,sim_points);
for i = 1:sim_points
    T(i) = NEGF_transmission(result.NEGF_result{i});
end
plot(E./eV,T)
axis([0,E(end)./eV,0,1])

%%
figure(1)
for i = 1:sim_points
    subplot(1,2,2);
    plot(NEGF_result_remap(result.NEGF_result{i},'fermi'));
    hold on
    plot(ones(2,1)*marg_len,[0,1],'color', [1,0,0,0.5], 'LineWidth',2);
    plot(ones(2,1)*(sample.length-marg_len),[0,1],'color', [1,0,0,0.5], 'LineWidth',2);
    hold off
    axis([1,sample.length,0,1])
    pause(0.05);
    subplot(1,2,1);
    plot(E,T);
    hold on
    plot(E(i),T(i),'rx','LineWidth',2);
    hold off
    axis([E(1),E(end),0,1])
end

%%
figure(1)
imagesc(sample.getUnits)
figure(2)
for i = 1:sim_points
    imagesc(NEGF_result_remap(result.NEGF_result{i},'fermi'))
    pause(0.05);
end











%% Transmittance data
clear
load('Transmittance_animation.mat'); 
val1 = 1;
val2 = 2;

data = amp_figures;

figure(4)
for i = 1:length(E)
    subplot(1,2,1)
    imagesc(contin_mat(amp_figures{val1}{i}))
    subplot(1,2,2)
    imagesc(contin_mat(amp_figures{val2}{i}))
    title("E = " + E(i))
    pause(0.1)
end

%% Hall effect data
clear
load('Hall_effect_animation.mat'); 
val1 = 1;
val2 = 4;

data = corr_figures;

figure(5)
for i = 1:length(B)
    subplot(1,2,1)
    imagesc(contin_mat(amp_figures{val1}{i}))
    subplot(1,2,2)
    imagesc(contin_mat(amp_figures{val2}{i}))
    title("B = " + B(i))
    pause(0.05)
end






%%
