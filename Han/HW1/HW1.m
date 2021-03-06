% HW1
% written by Seong-Min, Han  
% Due: AM 08:00, March 8, 2022
% 1. Solve the Laplace equation for a 1D structure.
% Length: Supplied by users
% Boundary conditions (Dirichlet boundary conditions at both ends): Supplied by users
% Mesh: Invent a clever way to build a mesh (non-equidistant grid points)
% ------------------------------------------------------------------------------------

% region : region1 | region2 | region3 | region3 | region2 | region1
region1 = 50; 
region2 = 30;
region3 = 20;
L = 2*(region1+region2+region3); 

% mesh
dx1 = 10; 
dx2 = dx1 / 5; 
dx3 = dx1 / 10; 

% node 
N = 2*(region1/dx1 + region2/dx2 + region3/dx3)+1;

% interface
it1 = region1 / dx1 + 1; 
it2 = it1 + region2 / dx2;
it4 = it2 + 2*(region3 / dx3);
it5 = it4 + region2 / dx2;

A = zeros(N,N);
b = zeros(N,1);
X = zeros(N,1);

for i=1:1:N
    if i==1
        A(i,i) = 1;
        b(i,1) = 0;
    elseif i==N
        A(i,i) = 1;
        b(i,1) = 1;
    elseif i == it1 
        A(i,i-1) = 1/dx1;
        A(i,i) = -1/dx1 -1/dx2;
        A(i,i+1) = 1/dx2;   
        b(i,1) = 0;
    elseif i == it5
        A(i,i-1) = 1/dx2;
        A(i,i) = -1/dx1 -1/dx2;
        A(i,i+1) = 1/dx1;   
        b(i,1) = 0;
    elseif  i < it1 || i > it5
        A(i,i-1) = 1/dx1;
        A(i,i) = -2/dx1;
        A(i,i+1) = 1/dx1;
        b(i,1) = 0;
    elseif (it1 < i && i < it2) || (it4 < i && i < it5)
        A(i,i-1) = 1/dx2;
        A(i,i) = -2/dx2;
        A(i,i+1) = 1/dx2;
        b(i,1) = 0;
    elseif i == it2 
        A(i,i-1) = 1/dx2;
        A(i,i) = -1/dx2 -1/dx3;
        A(i,i+1) = 1/dx3;
        b(i,1) = 0;
    elseif i == it4
        A(i,i-1) = 1/dx3;
        A(i,i) = -1/dx2 -1/dx3;
        A(i,i+1) = 1/dx2;
        b(i,1) = 0;
    else
        A(i,i-1) = 1/dx3;
        A(i,i) = -2/dx3;
        A(i,i+1) = 1/dx3;
        b(i,1) = 0;
    end
end

phi = A\ b;

for i=1:1:N
    if i == 1
        X(i,1) = 0;
    elseif i <= it1
        X(i,1) = X(i-1,1) + dx1;
    elseif i > it1 && i <= it2
        X(i,1) = X(i-1,1) + dx2;
    elseif i > it4 && i <= it5
        X(i,1) = X(i-1,1) + dx2;
    elseif i > it5
        X(i,1) = X(i-1,1) + dx1;
    else
        X(i,1) = X(i-1,1) + dx3;   
    end
end

plot(X, phi, 'o');
xlabel('Position[nm]')
ylabel('Electrostatic potential(V)')
title('1D Laplacian equation solver')
