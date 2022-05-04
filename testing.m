clear
sa = Sample(20,20,1,linspace(1,0.001,20));
H = hamiltonian(sa);
imagesc(H);