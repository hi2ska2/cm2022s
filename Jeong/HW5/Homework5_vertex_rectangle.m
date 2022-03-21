clear; close all; clc;

% 변수
a=[5 5]; % size
N=0.5; % 점 사이의 간격
Regionx=[0 5];
Regiony=[2 3];

%%%%%%%%%%%%%%%%%%%%%%% Vertex %%%%%%%%%%%%%%%%%%%%%%%
x=(0:N:a(1,1))';
y=(0:N:a(1,2))';
% 원 쪼개준것을 기준으로 좌표점 설정
point_tmp=zeros(length(x)*length(y),2);
for ii=1:length(y)
    for j=1:length(x)
        for m=1:2
            v=[x(j) y(ii)];
            k=length(x)*ii-(length(x)-1)+(j-1);
            point_tmp(k,m)=v(:,m);
        end
    end
end

figure(1); plot(point_tmp(:,1),point_tmp(:,2), '*');

% 좌표가 Region 내부의 점인지 확인
n=1;
for ii=1:length(point_tmp)
    x_tmp=point_tmp(ii,1); y_tmp=point_tmp(ii,2);
    if x_tmp>=Regionx(1,1) && x_tmp<=Regionx(1,2) && y_tmp>=Regiony(1,1) && y_tmp<=Regiony(1,2)
        point(n,:)=[point_tmp(ii,:) 0];
        n=n+1;
    end
end

writematrix(point(:,[1 2]),'Vertex_2.txt'); % vertex 파일 출력
% writematrix(Element,'Element.txt'); % element 파일 출력
% writematrix(phi,'potential'); % element 파일 출력
% writematrix(Dirichlet_BC,'Dirichlet_BC'); % Dirichlet_BC 파일 출력

figure(2); plot(point(:,1), point(:,2), '*');
axis([0 a(1,1) 0 a(1,2)])
%%%%%%%%%%%%%%%%%%%%%%% Vertex end %%%%%%%%%%%%%%%%%%%%%%%