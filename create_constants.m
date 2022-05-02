save("temp");
clear;

eV = 1.602176634e-19; %J
h = 6.62607015e-34; %J * Hz^-1
h_bar = h / (2*pi);
electron_mass = 9.1093837015E-31; %kg
proton_mass = 1.67262192369E-27; %kg
e_charge = 1.602176634E-19; %C

save("physical_constants.mat");
clear;
load("temp");
delete("temp.mat");
