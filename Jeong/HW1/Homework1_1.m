clear; close all; clc;

w1= 6; w2=11; w3=8; % width
Width=w1+w2+w3;

%단순히 mesh사이즈를 다르게 (delta_x1=0.5, delta_x2=0.2, delta_x3=0.4)
dx1=0.5; dx2=0.2; dx3=0.4;
x1=0:dx1:w1;  x2=w1+dx2:dx2:w1+w2; x3=w1+w2+dx3:dx3:w1+w2+w3;
x=[x1 x2 x3];
N=length(x);
interface1=length(x1);
interface2=N-length(x3);

dx(1,1)=0.5;
for i=2:length(x)-1
    dx(i,1)=0;
    dx(i,1)= x(1,i+1)-x(1,i);
end

A=zeros(N,N);
A(1,1)=1; A(N,N)=1; 
for i=2:N-1
    v=[1/dx(i-1,1) -1/dx(i-1,1)-1/dx(i,1) 1/dx(i,1)];
    for j=1:3
        A(i,i+j-2)=v(1,j);
    end
end

b=[zeros(N-1,1); ones(1,1)];
phi=A\b;

%Analytic solution
x_anl=0:0.001:Width;
phi_analytic=x_anl/(Width);

%%% Graph %%%
figure('Name','HW1')
plot(x_anl ,phi_analytic, '-',x,phi, 'ro');
xlabel("position(nm)");
ylabel("Electrostatic potential");
legend('Analytic Solution','Numerical Solution')
title("HomeWork1-1");
