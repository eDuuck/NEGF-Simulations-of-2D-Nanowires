clear
load("Hall_effect_data.mat")
load("physical_constants.mat")
%% Hall effect data
len = length(amp_data.res{1}.completed);
amp_figures = cell(1,length(amp_data.res));
corr_figures = cell(1,length(corr_data.res));

f = waitbar(0);
for j = 1:length(amp_data.res)
    amp_figures{j} = cell(1,len);
    for i = 1:len
        mat = NEGF_result_remap(amp_data.res{j}.NEGF_result{i},'electrons');
        amp_figures{j}{i} = disc_mat(mat,'linear');
        waitbar(i/len,f,"Amp data: " + j + "/" + length(amp_data.res))
    end
end

for j = 1:length(amp_data.res)
    corr_figures{j} = cell(1,len);
    for i = 1:len
        mat = NEGF_result_remap(amp_data.res{j}.NEGF_result{i},'electrons');
        corr_figures{j}{i} = disc_mat(mat,'linear');
        waitbar(i/len,f,"Amp data: " + j + "/" + length(amp_data.res))
    end
end
E = amp_data.res{1}.E;
B = amp_data.res{1}.B;
save("Hall_effect_animation.mat","amp_figures","corr_figures","E","B");


%%
clear
load("Trasmittance_data.mat")
load("physical_constants.mat")
%% Amp data
len = length(amp_data.res{1}.completed);
amp_figures = cell(1,length(amp_data.res));
corr_figures = cell(1,length(corr_data.res));

f = waitbar(0);
for j = 1:length(amp_data.res)
    amp_figures{j} = cell(1,len);
    for i = 1:len
        mat = NEGF_result_remap(amp_data.res{j}.NEGF_result{i},'electrons');
        amp_figures{j}{i} = disc_mat(mat,'linear');
        waitbar(i/len,f,"Amp data: " + j + "/" + length(amp_data.res))
    end
end

for j = 1:length(amp_data.res)
    corr_figures{j} = cell(1,len);
    for i = 1:len
        mat = NEGF_result_remap(amp_data.res{j}.NEGF_result{i},'electrons');
        corr_figures{j}{i} = disc_mat(mat,'linear');
        waitbar(i/len,f,"Amp data: " + j + "/" + length(amp_data.res))
    end
end
E = amp_data.res{1}.E;
B = amp_data.res{1}.B;
save("Transmittance_animation.mat","amp_figures","corr_figures","E","B");


% %%
% clear
% load("Transmittance_animation.mat")
% for i = 1:length(E)
%     imagesc(contin_mat(corr_figures(:,:,i)))
%     title("E = " + E)
% end
% %%
% clear
% load("Transmittance_animation.mat")
% for i = 1:length(E)
%     imagesc(contin_mat(corr_figures(:,:,i)))
%     title("E = " + E)
% end