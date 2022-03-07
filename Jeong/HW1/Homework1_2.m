clear; close all; clc;

w1= 6; w2=11; w3=8; Width=w1+w2+w3; % width

%기본적인 mesh사이즈는 동일하게 0.5로 설정
%경계점에서 t이내에서는 점차 줄어들게(m배로 감소로) 설정
t=2; m=0.8;
dx(1,1)=0.5;
interface1=w1; interface2=w1+w2;

% x설정
a=1; b=1; i=1; x(1,1)=0;
while x(i,1)<Width
    i=i+1; x(i,1)=0;
    if x(i-1,1)>=interface1-t && x(i-1,1)<interface1
        x(i,1)=x(i-1,1)+m*(x(i-1,1)-x(i-2,1));
        if x(i,1)-x(i-1,1)<=0.0001 || x(i,1)>interface1
            x(i,1)=x(i-1,1)+0.5*(interface1-x(i-1,1));
            i=i+1;
            x(i,1)=interface1;
        end

    elseif x(i-1,1)>=interface2-t && x(i-1,1)<interface2
        x(i,1)=x(i-1,1)+m*(x(i-1,1)-x(i-2,1));
        if x(i,1)-x(i-1,1)<=0.0001 || x(i,1)>interface2
            x(i,1)=x(i-1,1)+0.5*(interface2-x(i-1,1));
            i=i+1;
            x(i,1)=interface2;
        end
        
    elseif x(i-1,1)>=interface1 && x(i-1,1)<interface1+t
        if a<=3
            x(i,1)=x(i-1,1)+(x(i-1,1)-x(i-2,1));
            if a==3
                x(i,1)=x(i-1,1)+(x(i-5,1)-x(i-6,1));
            end
            a=a+1;
        else
            x(i,1)=x(i-1,1)+1/m*(x(i-1,1)-x(i-2,1));
            
            if x(i,1)>interface1+t && x(i,1)<interface2
               x(i,1)=interface1+t;
            end
        end
        
    elseif x(i-1,1)>=interface2 && x(i-1,1)<interface2+t
        if b<=3
            x(i,1)=x(i-1,1)+(x(i-1,1)-x(i-2,1));
            if b==3
                x(i,1)=x(i-1,1)+(x(i-5,1)-x(i-6,1));
            end
            b=b+1;
        else
            x(i,1)=x(i-1,1)+1/m*(x(i-1,1)-x(i-2,1));
            
            if x(i,1)>interface2+t
               x(i,1)=interface2+t;
            end
        end
    else 
    x(i,1)=x(i-1,1)+dx(1,1);
    end
end

N=length(x);
for i=2:length(x)-1
    dx(i,1)= x(i+1,1)-x(i,1);
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
title("HomeWork1-2");