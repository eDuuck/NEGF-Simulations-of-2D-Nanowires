clear
tic

Ec = 1;
a = 1;
t0 = 0.05/a^2;
eps0 = Ec + 4*t0;
t = -t0;
sample = struct();
sample.width = 100;
sample.length = 100;
sample.dim = [sample.width ,sample.length];
sample.M = prod(sample.dim);
sample.units = ones(sample.width,sample.length) * eps0;
sample.compressed = false;
sample.conn = t;
sample.arch = 'rectangular';

H = hamiltonian(sample);

filename = 'non_uni2D.mat';
if isfile(filename)
    load(filename);
    if ~isequaln(data.sample.units,sample.units)
        disp("Loaded file does not correspond to sample");
    end
else
    data.sim_points = 200;
    data.sample = sample_compress(sample);
    data.E = linspace(Ec,Ec+4*t0,data.sim_points);
    data.H = H;
    data.result(1:data.sim_points) = struct('complete',0,'E', ...
        0,'G',sparse(data.sample.M,data.sample.M));
    save(filename,'data');
end
toc
%%
for j = 1:data.sim_points
    if ~data.result(j).complete
        data.result(j).E = data.E(j);
        %Set Complete to 1
        %Save file'
        
        %InAs
    end
end
