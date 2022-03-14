clear; close all; clc;

%변수
r=2; % radius
N=16; % N은 4의 배수로 해주어야 오류가 안나도록 설정했다.
Dirichlet_BC = [1 0; 41 1]; % Dirichlet_BC = [index potential]

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

figure(1); plot(x_circle_tmp,y_circle_tmp, '-*');
figure(2); plot(point_tmp(:,1), point_tmp(:,2), '*');
figure(3); plot(x_circle_tmp,y_circle_tmp, '-*', point(:,1), point(:,2), '*');
%%%%%%%%%%%%%%%%%%%%%%% Vertex end %%%%%%%%%%%%%%%%%%%%%%%

% 원의 경계점들의 인덱스 계산, 계차수열 이용.
ind_circle=zeros((N/2+1),2);
for n=1:(N/2+1) % 아래쪽 반원
   if n<= (N/2+1)/2+1.5
       ind_circle(n,1)=1+(n-1)^2;
   else
       a_n1 = ind_circle((N/2+1)/2+1.5,1);
       a_n2 = ind_circle((N/2+1)/2+0.5,1);
       m=n-((N/2+1)/2+1.5)+1;
       ind_circle(n,1) = a_n1 + (a_n1-a_n2-m)*(m-1);
   end
end
for n=1:(N/2+1) % 윗쪽 반원
   if n<= (N/2+1)/2+0.5
       ind_circle(n,2)=n^2;
   else
       a_n1 = ind_circle((N/2+1)/2+0.5,2);
       a_n2 = ind_circle((N/2+1)/2-0.5,2);
       m=n-((N/2+1)/2+0.5)+1;
       ind_circle(n,2) = a_n1 + (a_n1-a_n2-m)*(m-1);
   end
end

%%%%%%%%%%%%%%%%%%%%%%% Element %%%%%%%%%%%%%%%%%%%%%%%
m=1;
for n=1:N/2
    if n<=N/4 % 원 아랫부분
        p=ind_circle(n,1);
        Element(m,1)=ind_circle(n,1); Element(m,2)=ind_circle(n+1,1)+1; Element(m,3)=ind_circle(n+1,1);
        m=m+1;
        while p ~= ind_circle(n,2)
            q=p+(ind_circle(n+1,1)-ind_circle(n,1))+1;
            Element(m,1)=p; Element(m,2)=q+1; Element(m,3)=q;
            Element(m+1,1)=p; Element(m+1,2)=p+1; Element(m+1,3)=q+1;
            m=m+2;
            p=p+1;
        end
        Element(m,1)=ind_circle(n,2); Element(m,2)=ind_circle(n+1,2); Element(m,3)=ind_circle(n+1,2)-1;
        m=m+1;
    else % 원 윗부분
        p=ind_circle(n,1)+1;
        Element(m,1)=ind_circle(n,1); Element(m,2)=ind_circle(n,1)+1; Element(m,3)=ind_circle(n+1,1);
        m=m+1;
        while p ~= ind_circle(n,2)-1
            q=p+(ind_circle(n+1,1)-ind_circle(n,1))-1;
            Element(m,1)=p; Element(m,2)=q+1; Element(m,3)=q;
            Element(m+1,1)=p; Element(m+1,2)=p+1; Element(m+1,3)=q+1;
            m=m+2;
            p=p+1;
        end
        Element(m,1)=ind_circle(n,2)-1; Element(m,2)=ind_circle(n,2); Element(m,3)=ind_circle(n+1,2);
        m=m+1;
    end
end
%%%%%%%%%%%%%%%%%%%%%%% Element end %%%%%%%%%%%%%%%%%%%%%%%

%각 벡터 구하기
vec=cell(length(Element),3);
for ii=1:length(Element)
    vec{ii,1}=point(Element(ii,2),:)-point(Element(ii,1),:);
    vec{ii,2}=point(Element(ii,3),:)-point(Element(ii,2),:);
    vec{ii,3}=point(Element(ii,1),:)-point(Element(ii,3),:);
end

%삼각형의 길이 구하기
L=zeros(length(Element),3);
for ii=1:length(Element)
    for j=1:3
        L(ii,j) = norm(vec{ii,j});
    end
end

% 삼각형 면적
Area=zeros(length(Element),1);
for ii=1:length(Element)
    Area(ii,1)=norm(0.5*cross(vec{ii,1}, -vec{ii,3}));
end

%외접의 반지름
R=zeros(length(Element),1);
for ii=1:length(Element)
   R(ii,1) = L(ii,1)*L(ii,2)*L(ii,3)/(4*Area(ii,1));
end

%외심에서 최단거리까지의 거리
A=zeros(length(Element),3);
for ii=1:length(Element)
    v=[0.5*L(ii,1) 0.5*L(ii,2) 0.5*L(ii,3)];
    for j=1:3
        A(ii,j) = real(sqrt(R(ii,1)^2-v(:,j)^2));
        A(A<=1e-7)=0;
    end
end

% Jaco matrix 구하기
J=zeros(length(point),length(point));
for ii=1:length(Element) % 일단 모든 점에서의 자코비안 행렬 구하기
    v=[-A(ii,1)/L(ii,1)-A(ii,3)/L(ii,3) A(ii,1)/L(ii,1) A(ii,3)/L(ii,3);
        A(ii,1)/L(ii,1) -A(ii,1)/L(ii,1)-A(ii,2)/L(ii,2) A(ii,2)/L(ii,2);
        A(ii,3)/L(ii,3) A(ii,2)/L(ii,2) -A(ii,2)/L(ii,2)-A(ii,3)/L(ii,3)];
    for j=1:3
        for k=1:3
        J(Element(ii,j),Element(ii,k))=J(Element(ii,j),Element(ii,k))+v(j,k);
        end
    end
end
for ii=1:size(Dirichlet_BC,1)
    J(Dirichlet_BC(ii,1),:)=0;
end
for ii=1:size(Dirichlet_BC,1)
    J(Dirichlet_BC(ii,1),Dirichlet_BC(ii,1))=1;
end

%res matrix 구하기ㅏ
res=zeros(length(point),1);
for ii=1:size(Dirichlet_BC,1)
    res(Dirichlet_BC(ii,1),1)=Dirichlet_BC(ii,2);
end

phi=J\res;

writematrix(point(:,[1 2]),'Vertex.txt'); % vertex 파일 출력
writematrix(Element,'Element.txt'); % element 파일 출력
writematrix(phi,'potential'); % element 파일 출력


figure(4); plot3(point(:,1),point(:,2),phi,'*')
xlabel("x")
ylabel("y")
zlabel("potential")
axis equal;