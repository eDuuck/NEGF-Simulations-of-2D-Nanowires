%% Load Quantum Hall Effect Data.
clear
load("Hall_effect_data.mat")
load("physical_constants.mat")
%% Amp-Variation
R = amp_data.R;
T = amp_data.T;

E = amp_data.E;
B = amp_data.B;

noiseLevels = amp_data.noiseLevels;
MesVar = length(noiseLevels);

IP = [512 554];

figure(1)
for i = 11%1:length(R(1,1,:))
    legends = cell(1,3);
    clf
    for j = 1:2:MesVar
        subplot(2,3,[1 4])
        hold on
        plot(B,R(j,:,i),'Color',[j/length(noiseLevels),0,1-j/length(noiseLevels), 0.9],'linewidth',1.5);
        hold off
        subplot(2,3,[2 5])
        hold on
        plot(B,T(j,:),'Color',[j/length(noiseLevels),0,1-j/length(noiseLevels)],'linewidth',1.5);
        hold off
        legends{(j+1)/2} = "\sigma = " + noiseLevels(j);
    end
    subplot(2,3,[1 4])
    hold on
    plot(B(IP),R(MesVar,IP,i),'ko','linewidth',1.5);
    hold off
    xlabel("B_z [T]"); ylabel("R_{xy} [R_K]")
    legend(legends, 'location', 'northwest')
    title("Quantum Hall Resistance")
    axis([B(1),B(end),-0.1,1*1.1]);
    grid; set(gca, 'FontWeight', 'bold')

    subplot(2,3,[2 5])
    hold on
    plot(B(IP),T(MesVar,IP),'ko','linewidth',1.5);
    hold off
    xlabel("B_z [T]"); ylabel("T")
    axis([B(1),B(end),-0.1,T(1,1)*1.1]);
    title("Transmission")
    grid; set(gca, 'FontWeight', 'bold')

    subplot(2,3,3)
    imagesc(NEGF_result_remap(amp_data.res{MesVar}.NEGF_result{IP(1)},'electrons'))
    title("B = " + B(IP(1)) + "T")
    subplot(2,3,6)
    imagesc(NEGF_result_remap(amp_data.res{MesVar}.NEGF_result{IP(2)},'electrons'))
    title("B = " + B(IP(2)) + "T")
    pause(1)
end


%% Corr-Variation
R = corr_data.R;
T = corr_data.T;

E = corr_data.E;
B = corr_data.B;

noiseSpacing = corr_data.noiseSpacing;
MesVar = length(noiseSpacing);

IP = [479 480];

figure(2)
for i = 11%1:length(R(1,1,:))
    printVals = [1,3,6];
    legends = cell(1,length(printVals));
    clf
    for j = 1:length(printVals)
        subplot(2,3,[1 4])
        hold on
        plot(B,R(printVals(j),:,i),'Color',[j/length(printVals),0,1-j/length(printVals), 0.9],'linewidth',1.5);
        hold off
        subplot(2,3,[2 5])
        hold on
        plot(B,T(printVals(j),:),'Color',[j/length(printVals),0,1-j/length(printVals)],'linewidth',1.5);
        hold off
        legends{j} = "C_l = " + noiseSpacing(printVals(j));
    end
    subplot(2,3,[1 4])
    hold on
    plot(B(IP),R(printVals(end),IP,i),'ko','linewidth',1.5);
    hold off
    xlabel("B_z [T]"); ylabel("R_{xy} [R_K]")
    legend(legends, 'location', 'northwest')
    title("Quantum Hall Resistance")
    axis([B(1),B(end),-0.1,1*1.1]);
    grid; set(gca, 'FontWeight', 'bold')

    subplot(2,3,[2 5])
    hold on
    plot(B(IP),T(printVals(end),IP),'ko','linewidth',1.5);
    hold off
    xlabel("B_z [T]"); ylabel("T")
    axis([B(1),B(end),-0.1,T(1,1)*1.1]);
    title("Transmission")
    grid; set(gca, 'FontWeight', 'bold')

    subplot(2,3,3)
    imagesc(NEGF_result_remap(corr_data.res{printVals(end)}.NEGF_result{IP(1)},'electrons'))
    title("B = " + B(IP(1)) + "T")
    subplot(2,3,6)
    imagesc(NEGF_result_remap(corr_data.res{printVals(end)}.NEGF_result{IP(2)},'electrons'))
    title("B = " + B(IP(2)) + "T")
    pause(1)
end

%% Generate transmittance data.
clear
load("Trasmittance_data.mat")  
load("physical_constants.mat")
%% Amp variation
T = amp_data.T;

E = amp_data.E;
B = amp_data.B;

noise = amp_data.noiseLevels;
MesVar = length(noise);

IP = 150;

figure(3)

printVals = [1,2,3,4,5];
legends = cell(1,length(printVals));
clf
subplot(2,4,[1 2 5 6])
hold on
for j = 1:length(printVals)
    plot(E./eV,T(j,:),'Color',[j/MesVar,0,1-j/MesVar],'linewidth',1.5);
    legends{j} = "\sigma = " + noise(printVals(j));
end
plot([1 1]*E(IP)./eV,[0 12],'r--','linewidth',1);
hold off
xlabel("E [eV]"); ylabel("T")
legend(legends, 'location', 'northwest')
title("Transmittance")
axis([E(1)./eV,E(end)./eV,-0.1,T(1,end)*1.1]);
grid; set(gca, 'FontWeight', 'bold')

dataType = 'electrons';
subplot(2,4,3)
imagesc(NEGF_result_remap(amp_data.res{printVals(1)}.NEGF_result{IP},dataType))
title("\sigma = " + noise(1))
subplot(2,4,4)
imagesc(NEGF_result_remap(amp_data.res{printVals(3)}.NEGF_result{IP},dataType))
title("\sigma = " + noise(3))
subplot(2,4,7)
imagesc(NEGF_result_remap(amp_data.res{printVals(4)}.NEGF_result{IP},dataType))
title("\sigma = " + noise(4))
subplot(2,4,8)
imagesc(NEGF_result_remap(amp_data.res{printVals(5)}.NEGF_result{IP},dataType))
title("\sigma = " + noise(5))

%% Corr variation
T = corr_data.T;

E = corr_data.E;
B = corr_data.B;

noise = corr_data.noiseSpacing;
MesVar = length(noise);

IP = 155;

figure(4)

printVals = [1,3,4,5];
legends = cell(1,length(printVals));
clf
subplot(2,4,[1 2 5 6])
hold on
for j = 1:length(printVals)
    plot(E./eV,T(printVals(j),:),'Color',[j/MesVar,0,1-j/MesVar],'linewidth',1.5);
    legends{j} = "C_l = " + noise(printVals(j));
end
plot([1 1]*E(IP)./eV,[0 12],'r--','linewidth',1);
hold off
xlabel("E [eV]"); ylabel("T")
legend(legends, 'location', 'northwest')
title("Transmittance")
axis([E(1)./eV,E(end)./eV,-0.1,T(1,end)*1.1]);
grid; set(gca, 'FontWeight', 'bold')

dataType = 'electrons';
subplot(2,4,3)
imagesc(NEGF_result_remap(amp_data.res{printVals(1)}.NEGF_result{IP},dataType))
title("C_l = " + noise(printVals(1)))
subplot(2,4,4)
imagesc(NEGF_result_remap(amp_data.res{printVals(2)}.NEGF_result{IP},dataType))
title("C_l = " + noise(printVals(2)))
subplot(2,4,7)
imagesc(NEGF_result_remap(amp_data.res{printVals(3)}.NEGF_result{IP},dataType))
title("C_l = " + noise(printVals(3)))
subplot(2,4,8)
imagesc(NEGF_result_remap(amp_data.res{printVals(4)}.NEGF_result{IP},dataType))
title("C_l = " + noise(printVals(4)))


%% Load compression data.
clear
load("Compress_data_2.mat")  
%%
figure(5)
subplot(1,3,1)
histogram(comp_data.comp_error);
title("Transmission error with compression")
ylabel("Counts"); xlabel("Error in transmission")
grid;  set(gca, 'FontWeight', 'bold')
subplot(1,3,2)
histogram(comp_data.comp_rate*100);
title("Compression rate")
xlabel("Data reduction [%]");ylabel("Counts")
grid;  set(gca, 'FontWeight', 'bold')
subplot(1,3,3)
timeIncrease = ((comp_data.cT_times+comp_data.comp_times ...
    +comp_data.NEGF_times)./(comp_data.T_times+comp_data.NEGF_times)-1)*100;
histogram(timeIncrease);
title("Time increase from compression")
xlabel("Time increase [%]");ylabel("Counts")
grid;  set(gca, 'FontWeight', 'bold')
