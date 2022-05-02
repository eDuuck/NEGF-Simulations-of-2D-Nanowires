contacts = 2;
load('physical_constants.mat')
m = electron_mass;
Ec = 1.12 * eV;
a = 4*1e-10; %Just taken out of the air, 4 ?.

% t0 = -h_bar^2/(2*m*a^2);
% eps0 = Ec - 2*t0;

%t2 = t0^2;
Ec = 1;
a = 1;
t0 = 0.05/a^2;
eps0 = Ec + 4*t0;

%Create a simple structure of 5x10 units.
temp = eps0 * ones(50,200);
struct = barrier(temp, 0, 1);
struct(5) = struct(10)*10;
struct(15) = struct(10)*10;
[width, length] = size(struct);
% 

E0 = Ec;
E_end = Ec + 2*t0;
sim_points = 100;
E = linspace(E0,E_end,sim_points);
% E = Ec+1*t0;
t = -t0;


[H,V] = twoDimmRed(struct,t);
D = 0;
%% Method 1
T = zeros(1,sim_points);
parCon = zeros(width,length,sim_points); %Parallel concentration, before change of basis.
con = zeros(width,length,sim_points);
Sigma = zeros(length, length, contacts, width);
for x = 1:sim_points
    for j = 1:width
        eps = H(1,1,j);
        k = acos((E(x) - eps) / (2 * t))/a;
        Sigma(1,1,1,j) = t*exp(1i*k*a);
        Sigma(end,end,contacts,j) = t*exp(1i*k*a);
        %parCon(j,1,x) = Sigma(1,1,1,j);
        %parCon(j,end,x) = Sigma(end,end,contacts,j);
        
        [G,Gn,Sigma0,Sigma0In] = dephase_mat(E(x),H(:,:,j),Sigma(:,:,:,j),[1 0],D,1E-3*t0^2);
        A = real(1i * (G - G'));
        if (min(abs(diag(A))) > 0)
            parCon(j,:,x) = diag(Gn)./diag(A);
        end
        %disp(parCon(j,:,x))
        Gamma = gammCalc(Sigma(:,:,:,j));
        T(x) = T(x) + real(trace(Gamma(:,:,1) *  G * Gamma(:,:,2) * G'));
    end
    
    for j = 1:length
        con(:,j,x) = diag(V(:,:,j)^-1*diag(parCon(:,j,x))*V(:,:,j));
    end
    disp(x+"/"+sim_points)
end
%% Method 2



%% Plotting

for x = 1:sim_points
    imagesc(real(con(:,:,x)), [0 1])
    title("E = " + E(x));
    pause(1/30);
end

function t_struct = barrier(structure, width, U)
    t_struct = structure;
    [h,l] = size(structure);
    for j = 0:width
        e = floor(0.9*l);
        s = e-width+j;
        t_struct(1+j, s:e) = t_struct(1+j, s:e) + U;
        t_struct(h-j, s:e) = t_struct(h-j, s:e) + U;
    end
end