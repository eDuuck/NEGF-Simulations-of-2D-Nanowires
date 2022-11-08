%%
clear
for E = 0.9:0.01:1.5
    tic
tau = -0.05;
eta = 0.001;
Ec = 1;
eps = Ec - 4*tau;
con = Contact(ones(10,1)*eps,tau,[1,1],1);
SGF = contact_surface(con,E,0,0);
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
sample.applyNoise(0.01,3);

% res = NEGF(sample,1.1);
% electron_den = NEGF_result_remap(res,'electrons');
% fermi_levels = NEGF_result_remap(res,'fermi');
% figure(2);
% subplot(1,2,1)
% imagesc(electron_den)
% subplot(1,2,2)
% imagesc(fermi_levels)

sim_points = 50;
results = struct('NEGF_results',cell(1,sim_points),'complete',zeros(1,sim_points));

filename = 'testing.mat';
if isfile(filename)
    load(filename);
else
    results = struct('complete',zeros(1,sim_points));
    results.NEGF_results = cell(1,sim_points);
end
E = linspace(1,1.37,sim_points);
B = 0.7;
k = 0;
disp("0/" + sim_points)
calc_timewo = zeros(1,sim_points);
for x = 1:sim_points
    tic
    if ~results.complete(x)
        G0 = 0;
%         if x > 1
%             G0 = results.NEGF_results{x-1}.G;
%         end
        results.NEGF_results{x} = NEGF(sample,E(x),B,1e-6,0.8,1000,false,G0);
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
%plot(calc_time/max(calc_time))
plot(T)
hold off
%% ANimations of NEGF Sims

if true
    for x = 1:sim_points
        electrons = NEGF_result_remap(results.NEGF_results{x},'electrons');
        fermi_levels = NEGF_result_remap(results.NEGF_results{x},'fermi');
        %disp(max(fermi_levels,[],'all'));
        subplot(1,2,1)
        imagesc(electrons,[0,100])
        subplot(1,2,2)
        imagesc(fermi_levels,[0,1])
        title("E = " + results.NEGF_results{x}.E);
        pause(0.1);
    end
end
figure(2)
plot(T,E)

%% Looking for patterns in G for potential compression methods
x = 20
m = sample.length;
for j = 1:sim_points
    subplot(1,2,1)
    image(real(results.NEGF_results{j}.G))
    subplot(1,2,2)
    image(imag(results.NEGF_results{j}.G))
    pause(0.3)
end
testmat = zeros(sample.M,15);
for j = 1:sample.length
diff = results.NEGF_results{x}.G((m*(j-1)+1):m*j,1:15)- ...
    results.NEGF_results{x}.G((m*j+1):m*(j+1),1:15);
testmat((m*(j-1)+1):m*j,1:15) = diff;
subplot(1,3,1)
image(real(diff))
subplot(1,3,2)
image(imag(diff))
subplot(1,3,3)
image(abs(diff))
pause(1)
end


%% Discretization of matrix test
x = 20;
res = results.NEGF_results{x}.G; %Original
temp = whos('G');
fileSize = temp.bytes;
time = 0;

tic
G_dis = lin_discretize(res);  %Discrete values
range = G_dis.range;
temp = whos('G_dis');
fileSize = [fileSize temp.bytes];
time = [time toc]; tic


Y = fft2(G_dis.matrix);     %Fourier compression on discrete values
Y = Y.*(abs(Y)>(max(abs(Y),[],'all'))/600);
K = sparse(Y);
temp = whos('K');
fileSize = [fileSize temp.bytes];
time = [time toc]; tic
%G_dis.matrix = ifft2(Y);

Z = fft2(res);                %Fourier compression on continous matrix
Z = Z.*(abs(Z)>(max(abs(Z),[],'all'))/20);
K = sparse(Z);
temp = whos('K');
fileSize = [fileSize temp.bytes];
time = [time toc]; tic

disp(fileSize)
%% Evaluation of compression methods
G_res = ifft2(Z);
G_res = contin_mat(G_dis);
figure(4)
subplot(3,2,1)
imagesc(real(res),range(1,:))
subplot(3,2,2)
imagesc(imag(res),range(2,:))
subplot(3,2,3)
imagesc(real(G_res),range(1,:))
subplot(3,2,4)
imagesc(imag(G_res),range(2,:))
subplot(3,2,5)
imagesc(real(G_res)-real(res))
subplot(3,2,6)
imagesc(imag(G_res)-imag(res))

max_error = max(abs(G_res-res),[],'all');
error_ratio = max_error / max(abs(res),[],'all')

%%
clear
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

sample.D = eye(sample.M)*1e-4;
sample.applyNoise(0.01,3);

sim_points = 50;
results = struct('complete',zeros(1,sim_points));
results.NEGF_results = cell(1,sim_points);

E = linspace(1,1.5,sim_points);
B = 0;

disp("0/" + sim_points)
for x = 1:sim_points
    tic
    if ~results.complete(x)
        G0 = 0;
%         if x > 1
%             G0 = results.NEGF_results{x-1}.G;
%         end
        results.NEGF_results{x} = NEGF(sample,E(x),B,1e-6,0.8,100,false,G0);
        results.complete(x) = 1;
        disp(x + " / " + sim_points)
    end
end
z = zeros(1,sim_points);
comparison = struct('original', z,'linear_discretiziation',z,'QOI',z);


%% Fourier compression
x = 10;
G = results.NEGF_results{x}.G(1:15,1:15);
fG = fft2(G);
rfG = fG .* (abs(fG) > 0.05*max(abs(fG),[],"all"));
drfG = lin_discretize(rfG);
subplot(2,2,1)
imagesc(abs(G))
subplot(2,2,2)
imagesc(abs(fG))
subplot(2,2,3)
imagesc(abs(rfG))
subplot(2,2,4)
imagesc(abs(double(drfG.matrix)))
A = compress(rfG,'QOI');
whos('drfG','A')

%%
t = zeros(1,50);
for k = 40:50
G = results.NEGF_results{k}.G;
diffMat = blockDifference(results.NEGF_results{k},'diagonal');
subplot(1,3,1)
imagesc(abs(G))
subplot(1,3,2)
imagesc(abs(diffMat))
subplot(1,3,3)
imagesc(abs(contin_mat(lin_discretize(G-G.'))))
tic
discMat = lin_discretize(G);
p(k) = toc; tic;
%compMat = compress(diffMat,'QOI');
compMatOrig = compress(G,'QOI');
Q(k) = toc;
whos('G','diffMat','discMat','compMat','compMatOrig')
t(k) = sum(abs(G-G.'),'all');
pause(10);
end
subplot(1,3,1)
plot(t)
subplot(1,3,2)
plot(Q)
subplot(1,3,3)
plot(p)


%%
for x = 1:50
G = results.NEGF_results{x}.G;
GC = compress(G,'qoi');
GL = lin_discretize(G);
Gr = decompress(GC,'qoi',G);
subplot(1,3,1)
imagesc(abs(G))
grid on
subplot(1,3,2)
imagesc(abs(Gr))
grid on
subplot(1,3,3)
imagesc(abs(G-Gr))
grid on
disp(sum(abs(G-Gr),"all")/length(G)^2)
whos('G',"GL","GC")
pause(0.5)
end

%%
getSize(results,2)
tic
for x = 1:50
    results.NEGF_results{x}.reduce(1);
    disp(x + " / " + 50)
end
disp("elapsed compression time : " +toc);
getSize(results,2)
tic
for x = 1:50
    results.NEGF_results{x}.reduce(0);
    disp(x + " / " + 50)
end
disp("elapsed decompression time : " +toc);
getSize(results,2)


%% Evaluating if 8bit is enough acc.
clear
load('physical_constants.mat');

Ec = 0;
a = 25E-10;
effective_mass = 0.1;
t0 =  h_bar^2/(2*electron_mass*effective_mass* a^2 * e_charge);
eps = 4*t0;
t =-t0;

sample = Sample(25,3,eps,t,a);
%sample.append(ones(sajmple.width,2)*eps*1.1);
%sample.append(ones(sample.width,5)*eps);
%sample.append(ones(sample.width,2)*eps*1.1);
%sample.append(ones(sample.width,5)*eps);
sample.addContact(ones(sample.width,1)*eps,t,[1,1]);
sample.addContact(ones(sample.width,1)*eps,t,[1,sample.length]);
sample.contacts{end}.fermi = 0;
sample.contacts{end}.face = -1;

sample.D = 0;%ones(sample.M)*1e-6;
%sample.applyNoise(0.1,3);


E = t0;
B = 0:0.1:10;%linspace(0,10,100);
NEGF_param = NEGF_param(sample,E,B,true);
NEGF_param.errorMarg = 1e-9;
NEGF_param.print = true;
NEGF_param.error_halt = true;

tic
result = NEGF_map(NEGF_param);
toc
disp(getSize(result,2))
 %%
X_vals = result.B;
I = zeros(1,length(result.E));
I_tot = zeros(1,length(X_vals));
V = zeros(1,length(X_vals));

midpoint = ceil(result.NEGF_result{1,1}.sample.length/2);
landua = zeros(length(X_vals),result.NEGF_result{1,1}.sample.width);
for y = 1:length(result.B)
    disp(result.B(y));
    for x = 1:length(result.E)
        res = result.NEGF_result{y,x};
        fermi = NEGF_result_remap(res,"fermi");
        
        imagesc(fermi,[0 1]);
        pause(0.5)

        I(x) = NEGF_transmission(res);
    end
    ibums = diag(res.getA());
    landua(y,:) = ibums(1:result.NEGF_result{1,1}.sample.width);
    I_tot(y) = I(end);%trapz(result.E,I);
    V(y) = fermi(1,midpoint)-fermi(end,midpoint);
end
figure(2)
plot(X_vals,V./I_tot);







function  [diff, bestIndex] = blockDifference(result,angle)
%This method does not seem to save any siginificant amount of data.
G = result.G;
sample = result.sample;
blockSide = sample.width;
diff = zeros(size(G));
if angle == "diagonal"
    bestIndex = uint16(zeros(2*sample.length-1,1));
    diagLength = [1:sample.length, sample.length-1:-1:1];
    for i = 1:2*sample.length-1 %Length of bestIndex
        offsetY = max([length(G)-i*blockSide+1, 1]);
        offsetX = max([i*blockSide-length(G)+1, 1]);
        smallestComb = ones(size(G))*Inf;
        for j = 1:diagLength(i)
            temp = zeros(size(G));
            diffIndStart = [offsetY+(j-1)*blockSide, offsetX+(j-1)*blockSide];
            diffIndEnd = [offsetY+(j)*blockSide-1, offsetX+(j)*blockSide-1];

            temp(diffIndStart(1):diffIndEnd(1),diffIndStart(2):diffIndEnd(2)) = ...
                G(diffIndStart(1):diffIndEnd(1),diffIndStart(2):diffIndEnd(2));

            for k = 1:diagLength(i)
                if k ~= j
                    startIndex = [offsetY+(k-1)*blockSide, offsetX+(k-1)*blockSide];
                    endIndex = [offsetY+(k)*blockSide-1, offsetX+(k)*blockSide-1];

                    temp(startIndex(1):endIndex(1),startIndex(2):endIndex(2)) = ...
                        G(startIndex(1):endIndex(1),startIndex(2):endIndex(2)) - ...
                        G(diffIndStart(1):diffIndEnd(1),diffIndStart(2):diffIndEnd(2));
                end
            end

            if sum(abs(temp),'all') < sum(abs(smallestComb),'all')
                bestIndex(i) = j;
                smallestComb = temp;
            end
        end
    end
    for i = 1:2*sample.length-1 %Length of bestIndex
        offsetY = max([length(G)-i*blockSide+1, 1]);
        offsetX = max([i*blockSide-length(G)+1, 1]);
        sIC = [offsetY+(bestIndex(i)-1)*blockSide, offsetX+(bestIndex(i)-1)*blockSide];
        eIC = [offsetY+(bestIndex(i))*blockSide-1, offsetX+(bestIndex(i))*blockSide-1];

        diff(sIC(1):eIC(1),sIC(2):eIC(2)) = ...
            G(sIC(1):eIC(1),sIC(2):eIC(2));

        for k = 1:diagLength(i)
            if k ~= bestIndex(i)
                sI = [offsetY+(k-1)*blockSide, offsetX+(k-1)*blockSide];
                eI = [offsetY+(k)*blockSide-1, offsetX+(k)*blockSide-1];

                diff(sI(1):eI(1),sI(2):eI(2)) = ...
                    G(sI(1):eI(1),sI(2):eI(2)) - ...
                    G(sIC(1):eIC(1),sIC(2):eIC(2));
            end
        end
    end
elseif angle == "vertical"
    bestIndex = uint16(zeros(sample.length,1));
    for i = 1:sample.length %Length of bestIndex
        offsetY = 1;
        offsetX = 1+blockSide*(i-1);
        smallestComb = ones(size(G))*Inf;
        for j = 1:sample.length
            temp = zeros(size(G));
            diffIndStart = [offsetY+(j-1)*blockSide, offsetX];
            diffIndEnd = [offsetY+(j)*blockSide-1, offsetX+blockSide-1];

            temp(diffIndStart(1):diffIndEnd(1),diffIndStart(2):diffIndEnd(2)) = ...
                G(diffIndStart(1):diffIndEnd(1),diffIndStart(2):diffIndEnd(2));

            for k = 1:sample.length
                if k ~= j
                    startIndex = [offsetY+(k-1)*blockSide, offsetX];
                    endIndex = [offsetY+(k)*blockSide-1, offsetX+blockSide-1];

                    temp(startIndex(1):endIndex(1),startIndex(2):endIndex(2)) = ...
                        G(startIndex(1):endIndex(1),startIndex(2):endIndex(2)) - ...
                        G(diffIndStart(1):diffIndEnd(1),diffIndStart(2):diffIndEnd(2));
                end
            end

            if sum(abs(temp),'all') < sum(abs(smallestComb),'all')
                bestIndex(i) = j;
                smallestComb = temp;
                imagesc(abs(temp))
            end
        end
    end
    for i = 1:sample.length
        offsetY = 1;
        offsetX = 1+blockSide*(i-1);
        sIC = [offsetY+(bestIndex(i)-1)*blockSide, offsetX];
        eIC = [offsetY+(bestIndex(i))*blockSide-1, offsetX+blockSide-1];

        diff(sIC(1):eIC(1),sIC(2):eIC(2)) = ...
            G(sIC(1):eIC(1),sIC(2):eIC(2));

        for k = 1:sample.length
            if k ~= bestIndex(i)
                sI = [offsetY+(k-1)*blockSide, offsetX];
                eI = [offsetY+(k)*blockSide-1, offsetX+blockSide-1];

                diff(sI(1):eI(1),sI(2):eI(2)) = ...
                    G(sI(1):eI(1),sI(2):eI(2)) - ...
                    G(sIC(1):eIC(1),sIC(2):eIC(2));
            end
        end
    end
end
end
