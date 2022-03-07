% due to : 2022.03.08 (Thu) / HW1

% 1D Laplace equation (non-equidistant grid points)

N=31;       % initial node number
L=300;      % Length

interface1=150;      % Dense
interface2=240;

DFactor = 4;    % dense factor

dx1=3*L/(N-1);
dx2=dx1/DFactor;
dx3=dx2/DFactor;

N1=interface1/dx1;   % Total Node number of Region1
N2=(interface2-interface1)/dx2;     % Total Node number of Region1
N3=(L-interface2)/dx3;      % Total Node number of Region1

Nt=N1+N2+N3+1;        %Total Node Number

A=sparse(Nt,Nt);
b=zeros(Nt,1);

% Dirichlet Boundary condition
A(1,1)=1;
b(1,1)=1;

A(Nt,Nt)=1;

for ii=2:Nt-1

    if (ii-1) < N1
        A(ii,ii-1)=(1/dx1);
        A(ii,ii)=(-2/dx1);
        A(ii,ii+1)=(1/dx1);
        
    elseif (ii-1) == N1
        A(ii,ii-1)=(1/dx1);
        A(ii,ii)=-(1/dx1)-(1/dx2);
        A(ii,ii+1)=(1/dx2);
        
    elseif N1 < ii-1 && ii-1 < N1+N2
        A(ii,ii-1)=(1/dx2);
        A(ii,ii)=(-2/dx2);
        A(ii,ii+1)=(1/dx2);

    elseif ii-1 == N1+N2
        A(ii,ii-1)=(1/dx2);
        A(ii,ii)=-(1/dx2)-(1/dx3);
        A(ii,ii+1)=(1/dx3);

    else    % Region3
        A(ii,ii-1)=(1/dx3);
        A(ii,ii)=(-2/dx3);
        A(ii,ii+1)=(1/dx3);
    end
end

Phi = A\b;

% Indexing

X=zeros(Nt,1);
X(1,1)=0;
X(Nt,1)=L;

for jj=2:Nt-1
    
    if jj <= N1+1
        X(jj,1)= dx1*(jj-1);
    end
    
    if N1+1 < jj && jj <= N1+N2+1
        X(jj,1)= interface1+dx2*(jj-(N1+1));
    end
    
    if N1+N2+1 < jj
        X(jj,1)= interface2+dx3*(jj-(N1+N2+1));
    end
end


%%%%%%% Uniform mesh %%%%%%%%

dx=L/(N-1);
J=sparse(N,N);
r=zeros(N,1);

J(1,1)=1;
r(1,1)=1;

J(N,N)=1;

for ii=2:N-1        
        J(ii,ii-1)=1;
        J(ii,ii)=(-2);
        J(ii,ii+1)=(1);
end

Phi_Uni = J\r;

plot(X,Phi,'r--o', 0:10:300, Phi_Uni,'b--','LineWidth',1.5)
title('Laplace equation','FontSize',20)
xlabel('Position (nm)','FontSize',15)
ylabel('Electrostatic Potential (V)','FontSize',15)

legend('Non-equidistant grid points','Equidistant grid points','FontSize',15)