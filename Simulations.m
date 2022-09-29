M = 40;
contacts = 2;
animation = false;
%create_constants
load('physical_constants.mat')
m = electron_mass;
Ec = 1.12 * eV;
a = 2*1e-10; %Just taken out of the air, 4 ?.
t = -h_bar^2/(2*m*a^2);
eps = Ec - 2*t;

E0 = Ec + 0 * t;
E_end = Ec - 4 * t;
sim_points = 200;

E = linspace(E0,E_end,sim_points);



H = diag(ones(1,M)*eps) + diag(ones(1,M-1)*t,1) + diag(ones(1,M-1)*t,-1);
e = ones(M,1);
%H = spdiags([t*e, eps*e, t*e], -1:1,M,M);
B1 = floor(M/3);
B2 = floor(2*M/3);
H(B1,B1) = H(B1,B1) + 1 * eV;
H(B2,B2) = H(B2,B2) + 1 * eV;

Sigma = zeros(M,M,contacts);
Sigma0 = zeros(M,M);
Sigma0In = zeros(M,M);
Gamma = zeros(M,M,contacts);
G = zeros(M,M,sim_points);
Gn = zeros(M,M,sim_points);
A = zeros(M,M,sim_points);
T = zeros(1,sim_points);


D = 100e-5*t^2*diag(ones(1,M));

tic
disp("Starting simulations, 0/" + sim_points)
for R = 1:sim_points
    k = acos((E(R) - eps) / (2 * t))/a;
    Sigma(1,1,1) = t*exp(1i*k*a);
    Sigma(end,end,2) = t*exp(1i*k*a);
    Gamma = gammCalc(Sigma);
    
    %G(:,:,R) = (E(R)*eye(M) - H - sum(Sigma,3))^-1;
    %Gn(:,:,R) = G(:,:,R)*Gamma(:,:,1)*G(:,:,R)';
    if false
        [G(:,:,R),Gn(:,:,R),Sigma0,Sigma0In] = dephase_mat(E(R),H,Sigma,[1 0],D,G(:,:,R-1),1E-3*t^2);
    else
        [G(:,:,R),Gn(:,:,R),Sigma0,Sigma0In] = dephase_mat(E(R),H,Sigma,[1 0],D,0,1E-3*t^2);
    end
    A(:,:,R) = real(1i * (G(:,:,R) - G(:,:,R)'));
    T(R) = trace(Gamma(:,:,1) * G(:,:,R) * Gamma(:,:,2) * G(:,:,R)');
    disp(R+"/"+sim_points)
end
toc
%%
if animation
    scale1 = [0 1 E0 E_end];
    scale2 = 0;
    scale3 = [0 M 0 1];
    for R = 2:sim_points
        tic
        subplot(1,2,1);
        axis(scale1);
        plot(real(T),E);
        hold on
        plot(real(T(R)),E(R),'r*');
        hold off
%         subplot(1,2,2);
%         plot(E,real(2*pi*h_bar/(e_charge^2)*(2*(1-T)./T + 1))); 
%         axis([min(E) max(E) 0 1E6]);
        subplot(1,2,2);
        plot(real(diag(Gn(:,:,R))./diag(A(:,:,R))));
        hold on
        plot(1:M,0.5*ones(1,M),'x')
        plot([B1 B2],[0.5 0.5],'or')
        hold off
        axis(scale3);
        pause(15/sim_points - toc);
    end
else
    subplot(1,2,1);
    plot(real(T),E);
    subplot(1,2,2);
    plot(E,real(2*pi*h_bar/(e_charge^2)*(2*(1-T)./T + 1))); 
    axis([min(E) max(E) 0 1E6]);
end

