clear; close all; clc;

Case=4; % Case 선택
if Case==1
    p1=[0 0 0]; p2=[1 0 0]; p3=[1 1 0]; % case1
elseif Case==2
    p1=[0 0 0]; p2=[3 0 0]; p3=[2 3 0]; % case2
elseif Case==3
    p1=[0 0 0]; p2=[3 0 0]; p3=[3 4 0]; % case3
elseif Case==4
p1=[0 0 0]; p2=[2 1 0]; p3=[1 2 0]; % case3
end
p=[p1; p2; p3];

%Vector
v12= p2-p1; v13= p3-p1; v23= p3-p2; % 벡터, cross함수를 사용하기 위해 z=0으로 고정 후 3차원으로 바꿈

%삼각형 각 변 길이
L12=norm(v12); L23=norm(v23); L31=norm(v13); L=[L12; L23; L31];

% 삼각형 면적
Area = norm(0.5*cross(v12, v13));

%외접원의 반지름
R = (L12*L31*L23)/(4*Area); % Area = (L1*L2*L3)/(4*Area) 이용

%외심에서 최단거리까지의 거리
A=zeros(3,1);
for ii=1:3
    A(ii,1)= sqrt(R^2-(0.5*L(ii,1))^2);
end

jaco=zeros(3,3);
jaco(1,1)=1; jaco(3,3)=1;
for i=2
    v=[A(i-1,1)/L(i-1,1) -A(i-1,1)/L(i-1,1)-A(i,1)/L(i,1) A(i,1)/L(i,1)];
    for j=1:3
        jaco(i,i+j-2)=v(1,j);
    end
end

res=[0; 0; 1];
phi=jaco\res;

%%% Graph %%%
figure('Name','HW2')
plot(1:3, phi, '--b*');
xlabel("Vector");
ylabel("Electrostatic potential");
title("HomeWork2");