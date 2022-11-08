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

wid = 1; len = 42;  %Width 1; Length 40.

sample = Sample(wid,len,eps,t,a);
sample.addContact(ones(sample.width,1)*eps,t,[1,1]);
sample.addContact(ones(sample.width,1)*eps,t,[1,sample.length]);
sample.contacts{end}.fermi = 0;

E_start = Ec + 0 * t;
E_end = Ec - 4 * t;
sim_points = 1000;

E = linspace(E_start,E_end,sim_points);
B = 0;

tic

%----------------Plotting the simulations for a 1D uniform wire------------
% figure(1);clf
% hold on
% D_vals = [-inf]%,-3,-2,-1];
% for j = 1:length(D_vals)
%     sample.D = eV*eps*10^D_vals(j) * diag(ones(1,sample.M));
%     if j == 1
%         legends{j} = "D_0 = 0";
%     else
%         legends{j} = "D_0 = \Gamma \cdot 1E"+ D_vals(j);
%     end
%     
%     params = NEGF_param(sample,E,B);
%     params.error_halt = false;
%     params.compress = false;
%     result = NEGF_map(params);
%     T = zeros(1,sim_points);
%     Curr = zeros(1,sim_points);
%     for i = 1:sim_points
%         T(i) = NEGF_transmission(result.NEGF_result{i});
%         Curr(i) = e_charge/h*T*ones(1,sim_points)'*(E_end-E_start)/sim_points;
%     end
% 
%     plot(E/eV,T,'Color',[j/length(D_vals),0,1-j/length(D_vals)])
% end
% toc
% hold off
% xlabel("E [eV]"); ylabel("Transmission");
% legend(legends); title("1D wire with no barrier.")


%----------------Adding double barrier to sample---------------------------
sample.units(floor(3*sample.length/7)) = eps + 1*eV;
sample.units(floor(4*sample.length/7)) = eps + 1*eV;
D_vals = [-inf,-2,-1.5];

it_lim = 10000;
error_halt = false;

starttic = tic;

figure(2);clf
tcl2 = tiledlayout(1,2);
figure(1);clf
tcl1 = tiledlayout(1,2);

phaseTot = 0;phaseComps = 0;
momTot = 0;momComps = 0;
while true
s = tic;
nexttile(2)
legends = cell(1,length(D_vals));
hold on
for j = 1:length(D_vals)
    sample.D = eV*eps*10^D_vals(j) * diag(ones(1,sample.M));  %Phase + momentum relaxation.
    if j == 1
        legends{j} = "D_0 = 0";
    else
        legends{j} = "D_0 = \Gamma \cdot 1E"+ D_vals(j);
    end
    params = NEGF_param(sample,E,B);
    params.it_lim = it_lim; params.error_halt = error_halt;
    result = NEGF_map(params);

    T = zeros(1,length(E));
    if j == 2
        figure(2)
        nexttile(2)
        hold on
        plot([-5,0,1:len,len+1,len+5], ...
            [1,1,NEGF_result_remap(result.NEGF_result{100},'fermi'),0,0],"linewidth", 1.5)
        plot(ones(1,2)*floor(3*sample.length/7),[-1,2],'Color',[1 0 0. 0.4],"linewidth", 3)
        plot(ones(1,2)*floor(4*sample.length/7),[-1,2],'Color',[1 0 0. 0.4],"linewidth", 3)
        axis([-5,len+5,-0.05,1.1]); set(gca, 'FontWeight', 'bold')
        title("Phase + momentum relaxati"); xlabel("Lattice point"); ylabel("Fermi Level");
        hold off
        grid
        figure(1)
    end
    for i = 1:length(T)
        res = result.NEGF_result{i};
        T(i) = NEGF_transmission(res);
    end
    plot(E/eV,T,'Color',[j/length(D_vals),0,1-j/length(D_vals)],'linewidth',1.5)
end
grid
hold off
axis([E(1)/eV, E(end)/eV, 0, 1]);
xlabel("E [eV]"); ylabel("Transmission");
legend(legends, 'Location', "southeast"); title("Phase + momentum relaxation")
set(gca, 'FontWeight', 'bold')

%toc(starttic)
momTot = momTot + toc(s); momComps = momComps + 1;
s = tic;
legends = cell(1,length(D_vals));
nexttile(1)
hold on
for j = 1:length(D_vals)
    sample.D = eV*eps*10^D_vals(j) * ones(sample.M);   %Only phase relaxation.
    if j == 1
        legends{j} = "D_0 = 0";
    else
        legends{j} = "D_0 = \Gamma \cdot 1E"+ D_vals(j);
    end

    params = NEGF_param(sample,E,B);
    params.it_lim = it_lim;
    params.error_halt = error_halt;
    result = NEGF_map(params);
    
    if j == 2
        figure(2)
        nexttile(1)
        hold on
        plot([-5,0,1:len,len+1,len+5], ...
            [1,1,NEGF_result_remap(result.NEGF_result{100},'fermi'),0,0],"linewidth", 1.5)
        plot(ones(1,2)*floor(3*sample.length/7),[-1,2],'Color',[1 0 0. 0.4],"linewidth", 3)
        plot(ones(1,2)*floor(4*sample.length/7),[-1,2],'Color',[1 0 0. 0.4],"linewidth", 3)
        axis([-5,len+5,-0.05,1.1]); set(gca, 'FontWeight', 'bold')
        title("Phase relaxation"); xlabel("Lattice point"); ylabel("Fermi Level"); 
        hold off
        grid
        figure(1)
    end

    T = zeros(1,length(E));
    for i = 1:length(T)
        T(i) = NEGF_transmission(result.NEGF_result{i});
    end

    plot(E/eV,T,'Color',[j/length(D_vals),0,1-j/length(D_vals)],'linewidth',1.5)
end
grid
hold off
axis([E(1)/eV, E(end)/eV, 0, 1]);
xlabel("E [eV]"); ylabel("Transmission");
set(gca, 'FontWeight', 'bold')
legend(legends, 'Location', "southeast"); title("Phase relaxation")
%toc(phaseTime)
phaseTot = phaseTot + toc(s); phaseComps = phaseComps + 1
break %This loop was just here for timing evaluations.
end
%% Generate data for the 2D-wire width variation.
% Initial setup of sample.
clear
load('physical_constants.mat');

Ec = 1.12* eV;
a = 10E-10;
effective_mass = 0.5;
t0 =  h_bar^2/(2*electron_mass*effective_mass* a^2);
t =-t0;
eps = Ec - 4*t;



E0 = Ec;
E_end = Ec - 4 * t;
sim_points = 200;

E = linspace(E0,E_end,sim_points);
B = 0;
D_vals = Ve*eps*10.^[-inf,-1];

figure(3);clf

legends = cell(1,length(D_vals));
hold on
widPoints = 60
T = zeros(1,widPoints);
for j = 1:widPoints
    wid = j; len = 5;  %Width 1; Length 40.
    sample = Sample(wid,len,eps,t,a);
    sample.addContact(ones(sample.width,1)*eps,t,[1,1]);
    sample.addContact(ones(sample.width,1)*eps,t,[1,sample.length]);
    sample.contacts{end}.fermi = 0;
    sample.D = 0;

    params = NEGF_param(sample,E(10),B);
    params.it_lim = 1000;
    params.error_halt = false;
    params.print = true;
    
    result = NEGF(params);

    if j == 17 || j == 40
        subplot(2,2,ceil(j/10))
        imagesc(NEGF_result_remap(result,'electrons'))
        title("Electron density, Width = " +j +" lattice points.");
    end
    disp(j + "/" + widPoints);
    T(j) = NEGF_transmission(result);
end
subplot(2,2,[1,3]);
hold on
plot(1:widPoints,T, 'linewidth',1.5)
plot([17 40],T([17 40]), 'rx', 'linewidth',1.5);
xlabel("Width [lattice points]"); ylabel("Transmission");
set(gca, 'FontWeight', 'bold')
grid

%% Generate data Quantum Hall Effect. - Amplitude variation

load('physical_constants.mat');


Ec = 0*eV; 
a = 50E-10; 
effective_mass = 0.1; 
t0 =  h_bar^2/(2*electron_mass*effective_mass* a^2); 
t =-t0;
eps = Ec + 4*t0;

wid = 22; len = 40;  

E = Ec + 10*t0/10;
B = 0:0.02:12;
B = 0:0.05:3;


figure(6);clf

hold on
noiseLevels = [0];% 0.01 0.025 0.05 0.1];
totTime = tic;
calcTimes = zeros(length(noiseLevels),1);



V = zeros(length(noiseLevels),length(B));

Rvar = 10;
R = zeros(length(noiseLevels),length(B),Rvar+1);


amp_data = struct;
amp_data.a = a;
amp_data.Ec = Ec;
amp_data.t0 = t0;
amp_data.eps = eps;
amp_data.eff_mass = effective_mass;
amp_data.noiseLevels = noiseLevels;
amp_data.res = cell(1,length(noiseLevels));        


for j = 1:length(noiseLevels)
    sample = Sample(wid,len,eps,t,a);
    sample.addContact(ones(wid,1)*eps,t,[1,1]);
    sample.addContact(ones(wid,1)*eps,t,[1,len]);
    sample.contacts{end}.fermi = 0;
    sample.contacts{end}.face = -1;
    sample.D = 0;%eV*eps*10^(-6) * ones(1,sample.M);
    sample.applyNoise(noiseLevels(j),wid/5);

    params = NEGF_param(sample,E,B);
    params.it_lim = 1000;
    params.error_halt = false;
    params.errorMarg = 1e-5;
    params.print = false;
    params.save_time = 600;
    
    startCalc = tic;
    
    disp("Finding sigma...")
    result = NEGF_map(params);
    amp_data.res{j} = result;

    disp("Solving R...")
    for k = 1:length(B)
        T(j,k) = NEGF_transmission(result.NEGF_result{k});
        fermi_map = NEGF_result_remap(result.NEGF_result{k},"fermi"); 

        %result.compress();

        v_range = 0;
        v_pon = ceil(len/2)-v_range:ceil(len/2)+v_range;
        V(j,k) = mean(fermi_map(1:(1+v_range),v_pon),"all")-mean(fermi_map(end-v_range:end,v_pon),"all");
        R(j,k,1) = V(j,k)/T(j,k);
        for l = 1:Rvar
            v_range = ceil(len/(2*Rvar)) * l;
            v_pon = ceil(len/2)-v_range+1:ceil(len/2)+v_range;
            V(j,k) = mean(fermi_map(1:min(floor(wid/8),v_range)+1,v_pon),"all")-mean(fermi_map(end-min(floor(wid/8),v_range):end,v_pon),"all");
            R(j,k,l+1) = V(j,k)/T(j,k);
        end
        disp(k + "/" + length(B));
    end
    calcTimes(j) = toc(startCalc);
    disp(j + "/" + length(noiseLevels));
    legends{j} = "A_n = " + noiseLevels(j);
%     plot(B,R(j,:),'Color',[j/length(noiseLevels),0,1-j/length(noiseLevels)],'linewidth',1.5);
%     plot(B,R1(j,:),'Color',[j/length(noiseLevels),0,1-j/length(noiseLevels)]*0.6,'linewidth',1.5);
    plot(B,R(j,:,end),'Color',[j/length(noiseLevels),0,1-j/length(noiseLevels)],'linewidth',1.5);
    pause(0.01);
end

xlabel("Amplitude of B-Field [T]"); ylabel("Hall Resistance [\Omega]")
legend(legends, 'location', 'northwest')
totTime = toc(totTime);
axis([B(1),B(end),-0.1,1*1.1]);
grid
set(gca, 'FontWeight', 'bold')
hold off

amp_data.B = B;
amp_data.E = E;
amp_data.R = R;
amp_data.T = T;
amp_data.V = V;
amp_data.totTime = totTime;
amp_data.calcTimes = calcTimes;

%% Generate data Quantum Hall Effect. - Corr_length variation

Ec = 0*eV;
a = 50E-10;
effective_mass = 0.1;
t0 =  h_bar^2/(2*electron_mass*effective_mass* a^2);
t =-t0;
eps = Ec + 4*t0;

wid = 22; len = 40;  %Width 1; Length 40.

E = Ec + 10*t0/10;
B = 0:0.02:12;


figure(7);clf

hold on
noiseSpacing = [wid*10 wid*5 wid wid/2 wid/5 wid/10 wid/20];
totTime = tic;
calcTimes = zeros(length(noiseSpacing),1);

T = zeros(length(noiseSpacing),length(B));
V = zeros(length(noiseSpacing),length(B));

Rvar = 10;
R_corr = zeros(length(noiseSpacing),length(B),Rvar+1);

corr_data = struct;
corr_data.a = a;
corr_data.Ec = Ec;
corr_data.t0 = t0;
corr_data.eps = eps;
corr_data.eff_mass = effective_mass;
corr_data.noiseSpacing = noiseSpacing;
corr_data.res = cell(1,length(noiseSpacing));      

for j = 1:length(noiseSpacing)
    sample = Sample(wid,len,eps,t,a);
    sample.addContact(ones(wid,1)*eps,t,[1,1]);
    sample.addContact(ones(wid,1)*eps,t,[1,len]);
    sample.contacts{end}.fermi = 0;
    sample.contacts{end}.face = -1;
    sample.D = 0;%eV*eps*10^(-6) * ones(1,sample.M);
    sample.applyNoise(0.025,noiseSpacing(j));

    params = NEGF_param(sample,E,B);
    params.it_lim = 1000;
    params.error_halt = false;
    params.errorMarg = 1e-5;
    params.print = false;
    params.save_time = 600;
    
    startCalc = tic;
    
    disp("Finding sigma...")
    result = NEGF_map(params);
    corr_data.res{j} = result;
    disp("Solving R...")
    for k = 1:length(B)
        T(j,k) = NEGF_transmission(result.NEGF_result{k});
        fermi_map = NEGF_result_remap(result.NEGF_result{k},"fermi"); 
        v_range = 0;
        v_pon = ceil(len/2)-v_range:ceil(len/2)+v_range;
        V(j,k) = mean(fermi_map(1:(1+v_range),v_pon),"all")-mean(fermi_map(end-v_range:end,v_pon),"all");
        R_corr(j,k,1) = V(j,k)/T(j,k);
        for l = 1:Rvar
            v_range = ceil(len/(2*Rvar)) * l;
            v_pon = ceil(len/2)-v_range+1:ceil(len/2)+v_range;
            V(j,k) = mean(fermi_map(1:min(floor(wid/8),v_range)+1,v_pon),"all")-mean(fermi_map(end-min(floor(wid/8),v_range):end,v_pon),"all");
            R_corr(j,k,l+1) = V(j,k)/T(j,k);
        end
        disp(k + "/" + length(B));
    end
    calcTimes(j) = toc(startCalc);
    if mod(j,floor(length(noiseSpacing)/100)) == 0
        disp(j + "/" + length(noiseSpacing));
    end
    legends{j} = "A_n = " + noiseSpacing(j);
    %plot(B,R(j,:),'Color',[j/length(noiseLevels),0,1-j/length(noiseLevels)],'linewidth',1.5);
    %plot(B,R1(j,:),'Color',[j/length(noiseLevels),0,1-j/length(noiseLevels)]*0.6,'linewidth',1.5);
    plot(B,R_corr(j,:,end),'Color',[j/length(noiseSpacing),0,1-j/length(noiseSpacing)],'linewidth',1.5);
    pause(0.01);
end
xlabel("Amplitude of B-Field [T]"); ylabel("Hall Resistance [\Omega]")
legend(legends, 'location', 'northwest')
totTime = toc(totTime);
axis([B(1),B(end),-0.1,1*1.1]);
grid
set(gca, 'FontWeight', 'bold')
hold off

corr_data.B = B;
corr_data.E = E;
corr_data.R = R_corr;
corr_data.T = T;
corr_data.V = V;
corr_data.totTime = totTime;
corr_data.calcTimes = calcTimes;

%% Generate data 2D wire transmission. - Amplitude variation
clear
load('physical_constants.mat');


Ec = 1.12*eV; 
a = 10E-10; 
effective_mass = 0.1; 
t0 =  h_bar^2/(2*electron_mass*effective_mass* a^2); 
t =-t0;
eps = Ec + 4*t0;

wid = 22; len = 40;  %Width 1; Length 40.

E = linspace(Ec,Ec + 2*t0, 200);
B = 0;


noiseLevels = [0 0.01 0.025 0.05 0.1];
totTime = tic;
calcTimes = zeros(length(noiseLevels),1);


amp_data = struct;
amp_data.a = a;
amp_data.Ec = Ec;
amp_data.t0 = t0;
amp_data.eps = eps;
amp_data.eff_mass = effective_mass;
amp_data.noiseLevels = noiseLevels;
amp_data.res = cell(1,length(noiseLevels));        


for j = 1:length(noiseLevels)
    sample = Sample(wid,len,eps,t,a);
    sample.addContact(ones(wid,1)*eps,t,[1,1]);
    sample.addContact(ones(wid,1)*eps,t,[1,len]);
    sample.contacts{end}.fermi = 0;
    sample.contacts{end}.face = -1;
    sample.D = 0;%eV*eps*10^(-6) * ones(1,sample.M);
    sample.applyNoise(noiseLevels(j),wid/5);

    params = NEGF_param(sample,E,B);
    params.it_lim = 1000;
    params.error_halt = false;
    params.errorMarg = 1e-5;
    params.print = false;
    params.save_time = 600;
    
    startCalc = tic;
    
    disp("Finding sigma...")
    result = NEGF_map(params);
    amp_data.res{j} = result;

    disp("Solving R...")
    for k = 1:length(E)
        T(j,k) = NEGF_transmission(result.NEGF_result{k});
        disp(k + "/" + length(E));
    end
    calcTimes(j) = toc(startCalc);
    disp(j + "/" + length(noiseLevels));
    legends{j} = "A_n = " + noiseLevels(j);
end

figure(8);
clf
hold on
for j = 1:length(noiseLevels)
    plot(E./eV,T(j,:),'Color',[j/length(noiseLevels),0,1-j/length(noiseLevels)],'linewidth',1.5);
end
hold off

xlabel("E [eV]"); ylabel("Transmittance [\Omega]")
legend(legends, 'location', 'northwest')
totTime = toc(totTime);
axis([E(1)/eV,E(end)/eV,-0.1,T(1,end)*1.1]);
grid
set(gca, 'FontWeight', 'bold')


amp_data.B = B;
amp_data.E = E;
amp_data.T = T;
amp_data.totTime = totTime;
amp_data.calcTimes = calcTimes;

%% Generate data 2D wire transmission. - Corr_length variation


Ec = 1.12*eV; 
a = 10E-10; 
effective_mass = 0.1; 
t0 =  h_bar^2/(2*electron_mass*effective_mass* a^2); 
t =-t0;
eps = Ec + 4*t0;

wid = 22; len = 40;  %Width 1; Length 40.

E = linspace(Ec,Ec + 2*t0, 200);
B = 0;


noiseSpacing = [wid*10 wid*5 wid  wid/10 wid/22];
totTime = tic;
calcTimes = zeros(length(noiseSpacing),1);

T = zeros(length(noiseSpacing),length(B));


corr_data = struct;
corr_data.a = a;
corr_data.Ec = Ec;
corr_data.t0 = t0;
corr_data.eps = eps;
corr_data.eff_mass = effective_mass;
corr_data.noiseSpacing = noiseSpacing;
corr_data.res = cell(1,length(noiseSpacing));      

for j = 1:length(noiseSpacing)
    sample = Sample(wid,len,eps,t,a);
    sample.addContact(ones(wid,1)*eps,t,[1,1]);
    sample.addContact(ones(wid,1)*eps,t,[1,len]);
    sample.contacts{end}.fermi = 0;
    sample.contacts{end}.face = -1;
    sample.D = 0;%eV*eps*10^(-6) * ones(1,sample.M);
    sample.applyNoise(0.025,noiseSpacing(j));

    params = NEGF_param(sample,E,B);
    params.it_lim = 1000;
    params.error_halt = false;
    params.errorMarg = 1e-5;
    params.print = false;
    params.save_time = 600;
    
    startCalc = tic;
    
    disp("Finding sigma...")
    result = NEGF_map(params);
    corr_data.res{j} = result;
    disp("Solving R...")
    for k = 1:length(E)
        T(j,k) = NEGF_transmission(result.NEGF_result{k});
        disp(k + "/" + length(E));
    end
    calcTimes(j) = toc(startCalc);
    disp(j + "/" + length(noiseSpacing));
    legends{j} = "A_l = " + noiseSpacing(j);
end

figure(9);
clf
hold on
for j = 1:length(noiseSpacing)
    plot(E./eV,T(j,:),'Color',[j/length(noiseSpacing),0,1-j/length(noiseSpacing)],'linewidth',1.5);
end
hold off
xlabel("E [eV]"); ylabel("Transmittance [\Omega]")
legend(legends, 'location', 'northwest')
axis([E(1)/eV,E(end)/eV,-0.1,T(1,end)*1.1]);
grid
set(gca, 'FontWeight', 'bold')


totTime = toc(totTime);
corr_data.B = B;
corr_data.E = E;
corr_data.T = T;
corr_data.totTime = totTime;
corr_data.calcTimes = calcTimes;
%% Compression rate
clear
load('physical_constants.mat');

Ec = 1.12*eV; 
a = 10E-10; 
effective_mass = 0.1; 
t0 =  h_bar^2/(2*electron_mass*effective_mass* a^2); 
t =-t0;
eps = Ec + 4*t0;

wid = 20; len = 10; %From 10 to 60.

E = linspace(Ec,Ec + 2*t0, 1000);
B = 0.1;

loops = 3;

noiseSpacing = wid;
noiseAmp = 0.025;
totTime = tic;

NEGF_times = zeros(loops,length(E));
comp_times = NEGF_times;
T_times = NEGF_times;
cT_times = NEGF_times;

T = zeros(loops,length(E));
cT = T;


comp_data = struct;
comp_data.a = a;
comp_data.Ec = Ec;
comp_data.t0 = t0;
comp_data.eps = eps;
comp_data.eff_mass = effective_mass;
comp_data.noiseSpacing = noiseSpacing;
comp_data.noiseAmp = noiseAmp;
comp_data.res = cell(loops,length(E));

comp_data.size = T;
comp_data.comp_size = T;
comp_data.comp_rate = T;
comp_data.comp_error = T;

sample = Sample(wid,len,eps,t,a);
sample.addContact(ones(wid,1)*eps,t,[1,1]);
sample.addContact(ones(wid,1)*eps,t,[1,len]);
sample.contacts{end}.fermi = 0;
sample.contacts{end}.face = -1;
sample.D = eV*eps*10^(-6) * ones(1,sample.M);
for j = 1:loops
    calcSam = copy(sample);
    calcSam.applyNoise(noiseAmp,noiseSpacing);

    for k = 1:length(E)
        params = NEGF_param(calcSam,E(k),B);
        params.it_lim = 1000;
        params.errorMarg = 1e-5;
        params.error_halt = false;
        if k > 1
            params.g0 = result;
        end

        startTime = tic;
        result = NEGF(params);
        result.compressionMethod = 'qoi';
        NEGF_times(j,k) = toc(startTime);
        
        startTime = tic;
        T(j,k) = NEGF_transmission(result);
        T_times(j,k) = toc(startTime);

        comp_data.size(j,k) = getSize(result,1);
        
        startTime = tic;
        result.compress();
        comp_times(j,k) = toc(startTime);
        
        startTime = tic;
        cT(j,k) = NEGF_transmission(result);
        cT_times(j,k) = toc(startTime);

        comp_data.comp_size(j,k) = getSize(result,1);
        comp_data.comp_error(j,k) = abs((cT(j,k)-T(j,k)));
        comp_data.comp_rate(j,k) = 1-comp_data.comp_size(j,k)/comp_data.size(j,k);

        comp_data.res{j,k} = result;
        disp(k)
    end
    
    disp("Done with calculations in loop " + j)
    
    sample.append(ones(wid,1)*eps*2)
    sample.append(ones(wid,9)*eps)
    sample.contacts{2}.pos = [1,sample.length];
    sample.D = eV*eps*10^(-6) * ones(1,sample.M);
end




totTime = toc(totTime);
comp_data.B = B;
comp_data.E = E;
comp_data.T = T;
comp_data.cT = cT;
comp_data.totTime = totTime;
comp_data.NEGF_times = NEGF_times; 
comp_data.comp_times = comp_times;
comp_data.T_times = T_times;
comp_data.cT_times = cT_times;


figure(10)
subplot(1,2,1)
histogram(comp_data.comp_error);
title("Transmission error with compression")
subplot(1,2,2)
histogram(comp_data.comp_rate);
title("Compression rate")

save("Compress_data_2","comp_data")


%% VWW (Very Wide Wire)
clear
load('physical_constants.mat');
Ec = 0* eV;
a = 1E-9;
effective_mass = 0.03; 
t0 =  h_bar^2/(2*electron_mass*effective_mass* a^2); 
t =-t0;
eps = Ec + 4*t0;

VWW_data = struct;

E = Ec + 0.1*t0;
B = 0:0.001:0.4;

wid = 500; len = 10;

sample = Sample(wid,len,eps,t,a);
sample.addContact(ones(wid,1)*eps,t,[1,1]);
sample.addContact(ones(wid,1)*eps,t,[1,len]);
sample.contacts{end}.fermi = 0;
sample.contacts{end}.face = -1;
sample.D = 0;

params = NEGF_param(sample,E,B);
params.error_halt = false;
params.errorMarg = 1E-4;
params.con_it_lim = 5000;
params.print = true;
params.compress = true;
startTic = tic;
result = NEGF_map(params);
conItTime = toc(startTic)


%%
VWW_data.E = E;
VWW_data.B = B;
VWW_data.result = result;
VWW_data.a = a;
VWW_data.Ec = Ec;
VWW_data.m = effective_mass;


startTic = tic;
for i = 1:106
    VWW_data.T(i) = NEGF_transmission(result.NEGF_result{i});
    disp(i)
end
transTime = toc(startTic)
save("VWW_data","VWW_data")