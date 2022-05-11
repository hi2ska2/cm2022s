clear; close all; clc;
eps0=8.854187817e-12;
q=1.602192e-19; %Elementary charge,C
mobility_n = 518e-4;  % electron mobility
Nd=2e23; % /m^3

f=logspace(10,14, 10)';

w=2*pi*f;

Y_Re=q*mobility_n*Nd*(5e-9*1e-6)/(120e-9);
Y_Re=Y_Re*ones(length(f),1);
Y_Im=w*11.7*eps0*(5e-9*1e-6)/(120e-9);

loglog(f*1e-12,Y_Re,'-o', f*1e-12,Y_Im,'-o')
xlabel("Frequency (THz)")
ylabel("Admittance")
legend("Re", "Im")