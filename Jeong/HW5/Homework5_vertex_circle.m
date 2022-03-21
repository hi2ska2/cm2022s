clear; close all; clc;

% 변수
r=2; % radius
N=16; % N은 4의 배수로 해주어야 오류가 안나도록 설정했다.

%%%%%%%%%%%%%%%%%%%%%%% Vertex %%%%%%%%%%%%%%%%%%%%%%%
center=[r r];
theta=linspace(0,2*pi,N+1)';
theta_tmp=theta;
theta=theta(1 : (length(theta)-1),1);
x_circle_tmp=r*cos(theta_tmp)+center(1,1); y_circle_tmp=(r*sin(theta_tmp)+center(1,2));
x_circle=sort(r*cos(theta)+center(1,1)); y_circle=sort(r*sin(theta)+center(1,2));
x_circle=unique((round(x_circle,8))); y_circle=unique((round(y_circle,8)));

% 원 쪼개준것을 기준으로 좌표점 설정
point_tmp=zeros(length(x_circle)*length(y_circle),2);
for ii=1:length(y_circle)
    for j=1:length(x_circle)
        for m=1:2
            v=[x_circle(j) y_circle(ii)];
            k=length(x_circle)*ii-(length(x_circle)-1)+(j-1);
            point_tmp(k,m)=v(:,m);
        end
    end
end

% 좌표가 원 내부에 있는 점인지 확인
n=1;
for ii=1:length(point_tmp)
    a=(point_tmp(ii,1)-center(1,1))^2+(point_tmp(ii,2)-center(1,2))^2; % x^2 + y^2 = r^2 이용
    a=round(a,6);
    if a<=r^2
        point(n,:)=[point_tmp(ii,:) 0];
        n=n+1;
    end
end

writematrix(point(:,[1 2]),'Vertex.txt'); % vertex 파일 출력
writematrix(Element,'Element.txt'); % element 파일 출력
writematrix(phi,'potential'); % element 파일 출력
% writematrix(Dirichlet_BC,'Dirichlet_BC'); % Dirichlet_BC 파일 출력

figure(1); plot(x_circle_tmp,y_circle_tmp, '-*');
figure(2); plot(point_tmp(:,1), point_tmp(:,2), '*');
figure(3); plot(x_circle_tmp,y_circle_tmp, '-*', point(:,1), point(:,2), '*');
%%%%%%%%%%%%%%%%%%%%%%% Vertex end %%%%%%%%%%%%%%%%%%%%%%%