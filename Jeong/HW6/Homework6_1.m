clear; close all; clc;

%유전율(변수)
eps1=11.7; eps2=3.9; eps3=11.7; eps4=3.9;

% element로 이루어진 경계조건 계산하기 위해 변환
Dirichlet_BC_tmp = importdata("Dirichlet_BC.txt"); % Dirichlet_BC = [index potential]
m=0; n=0;
for ii=1:size(Dirichlet_BC_tmp,1)
    for j=size(Dirichlet_BC_tmp,2):-1:1
        if ~isnan(Dirichlet_BC_tmp(ii,j))
            for k=1:j-1
                Dirichlet_BC(k+n,:)=[Dirichlet_BC_tmp(ii,k) Dirichlet_BC_tmp(ii,j)];
                m=m+1;
            end
            n=m;
            break;
        end
    end
end

% element 값 입력
F_1 = importdata("element_1_1.txt");
F_2 = importdata("element_1_2.txt");
F_3 = importdata("element_1_3.txt");
F_4 = importdata("element_1_4.txt");
eps1=ones(length(F_1),1).*eps1;
eps2=ones(length(F_2),1).*eps2;
eps3=ones(length(F_3),1).*eps3;
eps4=ones(length(F_4),1).*eps4;
Element = [F_1 eps1; F_2 eps2; F_3 eps3; F_4 eps4];

point = importdata("Vertex_circle.txt");
point_tmp=zeros(length(point),1);
point = [point point_tmp];

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
    v=Element(ii,4).*v;
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

%res matrix 구하기
res=zeros(length(point),1);
for ii=1:size(Dirichlet_BC,1)
    res(Dirichlet_BC(ii,1),1)=Dirichlet_BC(ii,2);
end

phi=J\res;

% 각 Region의 edge 정렬
n=1;
for ii=1:length(F_1) % Si region
   for j=1:3
      edge_1(n,:)= [F_1(ii,2) F_1(ii,1)];
      edge_1(n+1,:)= [F_1(ii,3) F_1(ii,2)];
      edge_1(n+2,:)= [F_1(ii,1) F_1(ii,3)];
      n=n+3;
   end
end
n=1;
for ii=1:length(F_2) % ox region
   for j=1:3
      edge_2(n,:)= [F_2(ii,1) F_2(ii,2)];
      edge_2(n+1,:)= [F_2(ii,2) F_2(ii,3)];
      edge_2(n+2,:)= [F_2(ii,3) F_2(ii,1)];
      n=n+3;
   end
end
n=1;
for ii=1:length(F_3) % Si region
   for j=1:3
      edge_3(n,:)= [F_3(ii,2) F_3(ii,1)];
      edge_3(n+1,:)= [F_3(ii,3) F_3(ii,2)];
      edge_3(n+2,:)= [F_3(ii,1) F_3(ii,3)];
      n=n+3;
   end
end
n=1;
for ii=1:length(F_4) % ox region
   for j=1:3
      edge_4(n,:)= [F_4(ii,1) F_4(ii,2)];
      edge_4(n+1,:)= [F_4(ii,2) F_4(ii,3)];
      edge_4(n+2,:)= [F_4(ii,3) F_4(ii,1)];
      n=n+3;
   end
end

% 중복 edge 정보 제거
edge_1=sort(edge_1, 2); % edge_si의 각 행을 오름차순으로 정렬(중복값 제거 위해)
edge_1=unique(edge_1,'rows');
edge_2=sort(edge_2, 2); % edge_ox의 각 행을 오름차순으로 정렬(중복값 제거 위해)
edge_2=unique(edge_2,'rows');
edge_3=sort(edge_3, 2); % edge_si의 각 행을 오름차순으로 정렬(중복값 제거 위해)
edge_3=unique(edge_3,'rows');
edge_4=sort(edge_4, 2); % edge_ox의 각 행을 오름차순으로 정렬(중복값 제거 위해)
edge_4=unique(edge_4,'rows');


% 각 region에서 동일한 edge가 있는지 확인, Interface
edge_interface1 = intersect(edge_1,edge_2,'rows');
edge_interface2 = intersect(edge_1,edge_3,'rows');
edge_interface3 = intersect(edge_1,edge_4,'rows');
edge_interface4 = intersect(edge_2,edge_3,'rows');
edge_interface5 = intersect(edge_2,edge_4,'rows');
edge_interface6 = intersect(edge_3,edge_4,'rows');
edge_interface=[edge_interface1; edge_interface2; edge_interface3; edge_interface4; edge_interface5; edge_interface6];

figure(4); plot3(point(:,1),point(:,2),phi,'*')
xlabel("x")
ylabel("y")
zlabel("potential")
axis equal;

% Visualize 
figure(5);
patch('Faces',Element(:,1:3),'Vertices',point(:,[1 2]), 'EdgeColor','black','FaceColor','none','LineWidth',0.5);
xlabel('X');
ylabel('Y');
title('Structure (N=20)')

figure(6);
patch('Faces',Element(:,1:3),'Vertices',point(:,[1 2]), 'FaceVertexCData',phi, 'EdgeColor','black','FaceColor','interp');
title('Structure (N=20)')
xlabel('X');
ylabel('Y');
colorbar
hold on
patch('Faces',F_1,'Vertices',point, 'EdgeColor','Blue','FaceColor','none','LineWidth',2)
hold on
patch('Faces',F_2,'Vertices',point, 'EdgeColor','Red','FaceColor','none','LineWidth',2)
hold on
patch('Faces',F_3,'Vertices',point, 'EdgeColor','Yellow','FaceColor','none','LineWidth',2)
hold on
patch('Faces',F_4,'Vertices',point, 'EdgeColor','magenta','FaceColor','none','LineWidth',2)
xlabel('X')
ylabel('Y')
title('Structure')
hold on
patch('Faces',edge_interface,'Vertices',point, 'EdgeColor','Green','FaceColor','none','LineWidth',5)
hold off
% caxis([0 2]);