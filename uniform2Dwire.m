clear

Ec = 1;
a = 1;
t0 = 0.05/a^2;
eps0 = Ec + 4*t0;
t = -t0;


width = 15; length = 10; M = width * length;
sample = struct('width',width,'length',length,'units',[],'conn',t,'H',[]);
sample.units = ones(width,length) * eps0;

alpha = @(i) diag(ones(width,1)).*diag(sample.units(:,i)) + ...
    sample.conn*(diag(ones(width-1,1),1) + diag(ones(width-1,1),-1));
beta = diag(ones(width,1))*t;
H = zeros(M);
H(1:width,1:width) = alpha(1);
for j = 1:length-1
    H(((width*j)+1):(width*(j+1)),((width*j)+1):(width*(j+1))) = alpha(j+1);
    H(((width*j)+1):(width*(j+1)),((width*(j-1))+1):(width*j)) = beta;
    H(((width*(j-1))+1):(width*(j)),((width*j)+1):(width*(j+1))) = beta;
end
%imagesc(H)
sample.H = H;

%%
[V,D] = eig(alpha(1));
T = [];

sim_points = 200;
E = linspace(Ec,Ec+4*t0,sim_points);


sim_points = 1;
E = 1.1;

Gn = zeros(M,M,sim_points);
n_e = zeros(M,sim_points);
k = @(E,eps) acos((E-eps)/(-2*t0))/a;

v = 0;%1e-3*t^2*diag(ones(1,M));
disp("Starting calculations: 0/" + sim_points);
for x = 1:sim_points
    Sigma = zeros(M,M,2);
    mode_contact = diag(t*exp(1i*k(E(x),diag(D))*a));
    Sigma(1:width,1:width,1) = V*mode_contact*V';
    imagesc(abs(Sigma(1:15,1:15,1)))
    Sigma((M-width+1):end,(M-width+1):end,2) = V*mode_contact*V';
    Gamma = gammCalc(Sigma);
    [G,Gn(:,:,x),Sigma0,Sigma0In] = dephase_mat(E(x),H,Sigma,[1 0],v);
    A = 1;%(1i * (G - G'));
    n_e(:,x) = real(diag(Gn(:,:,x)./diag(A)));
    T = [T ,real(trace(Gamma(:,:,1) * G * Gamma(:,:,2) * G'))];
    disp("                       " + x + "/" + sim_points);
end
%%
img = zeros(width,length);
figure(1)
for x = 1:sim_points
    for j = 1:length
        img(:,j) = n_e(((j-1)*width+1):(width*j),x);
    end
    imagesc(img);
    pause(0.05);
end

figure(2)
plot(T,E)
