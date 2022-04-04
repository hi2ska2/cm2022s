clear; close all; clc;

%유전율(변수)
eps1=3.9; eps2=11.7; eps3=3.9;

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
element1 = importdata("element_1.txt");
element2 = importdata("element_2.txt");
element3 = importdata("element_3.txt");
eps1=ones(length(element1),1).*eps1;
eps2=ones(length(element2),1).*eps2;
eps3=ones(length(element3),1).*eps3;
Element = [element1 eps1; element2 eps2; element3 eps3];

point = importdata("Vertex.txt");
point_tmp1=zeros(length(point),1);
point_tmp2 = [point point_tmp1];


%%%% 각 Region의 edge 정렬 %%%%
n=1;
for ii=1:length(element1) % Si region
   for j=1:3
      edge_1(n,:)= [element1(ii,2) element1(ii,1)];
      edge_1(n+1,:)= [element1(ii,3) element1(ii,2)];
      edge_1(n+2,:)= [element1(ii,1) element1(ii,3)];
      n=n+3;
   end
end
n=1;
for ii=1:length(element2) % ox region
   for j=1:3
      edge_2(n,:)= [element2(ii,1) element2(ii,2)];
      edge_2(n+1,:)= [element2(ii,2) element2(ii,3)];
      edge_2(n+2,:)= [element2(ii,3) element2(ii,1)];
      n=n+3;
   end
end
n=1;
for ii=1:length(element3) % Si region
   for j=1:3
      edge_3(n,:)= [element3(ii,2) element3(ii,1)];
      edge_3(n+1,:)= [element3(ii,3) element3(ii,2)];
      edge_3(n+2,:)= [element3(ii,1) element3(ii,3)];
      n=n+3;
   end
end
n=1;

% 중복 edge 정보 제거
edge_1=sort(edge_1, 2); % edge_si의 각 행을 오름차순으로 정렬(중복값 제거 위해)
edge_1=unique(edge_1,'rows');
edge_2=sort(edge_2, 2); % edge_ox의 각 행을 오름차순으로 정렬(중복값 제거 위해)
edge_2=unique(edge_2,'rows');
edge_3=sort(edge_3, 2); % edge_si의 각 행을 오름차순으로 정렬(중복값 제거 위해)
edge_3=unique(edge_3,'rows');


% 각 region에서 동일한 edge가 있는지 확인, Interface
edge_interface1 = intersect(edge_1,edge_2,'rows');
edge_interface2 = intersect(edge_1,edge_3,'rows');
edge_interface3 = intersect(edge_2,edge_3,'rows');
edge_interface=[edge_interface1; edge_interface2; edge_interface3];

% 각 region을 구성하는 vertex의 index를 확인하는 과정
Vertex_region1 = unique(sort(reshape(element1,[],1))); 
Vertex_region2 = unique(sort(reshape(element2,[],1)));
Vertex_region3 = unique(sort(reshape(element3,[],1)));

% 각 Reigon 내의 vertex 갯수 저장
Number_of_vertex=[length(Vertex_region1); length(Vertex_region2); length(Vertex_region3)];

% Re-point
Vertex_region=[Vertex_region1; Vertex_region2; Vertex_region3];
for ii=1:length(Vertex_region)
    for j=1:3
    point(ii,j)=point_tmp2(Vertex_region(ii,1),j);
    end
end
point = sortrows(point,2);

% re-element
vertex_interface=unique(sort(reshape(edge_interface,[],1)));
n=0;    % Region 1
for ii=1:length(element1)
    for j=1:3
        if element1(ii,j)>=1 && element1(ii,j)<=11
            element1(ii,j)=element1(ii,j)+(n);
        elseif element1(ii,j)>=12 && element1(ii,j)<=22
            element1(ii,j)=element1(ii,j)+(2+n);
        elseif element1(ii,j)>=23 && element1(ii,j)<=33
            element1(ii,j)=element1(ii,j)+(4+n);
        end
    end
end
n=1;    % Region 2
for ii=1:length(element2)
    for j=1:3
        if element2(ii,j)>=1 && element2(ii,j)<=11
            element2(ii,j)=element2(ii,j)+(n);
        elseif element2(ii,j)>=12 && element2(ii,j)<=22
            element2(ii,j)=element2(ii,j)+(2+n);
        elseif element2(ii,j)>=23 && element2(ii,j)<=33
            element2(ii,j)=element2(ii,j)+(4+n);
        end
    end
end
n=2;
for ii=1:length(element3)
    for j=1:3
        if element3(ii,j)>=1 && element3(ii,j)<=11
            element3(ii,j)=element3(ii,j)+(n);
        elseif element3(ii,j)>=12 && element3(ii,j)<=22
            element3(ii,j)=element3(ii,j)+(2+n);
        elseif element3(ii,j)>=23 && element3(ii,j)<=33
            element3(ii,j)=element3(ii,j)+(4+n);
        end
    end
end
Element = [element1 eps1; element2 eps2; element3 eps3];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
        %A(A<=1e-7)=0;
    end
end

% Jaco matrix 구하기
J=zeros(sum(Number_of_vertex),sum(Number_of_vertex));
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

% Interface_BC
Interface_BC=[3; 16; 29; 11; 24; 37];
for ii=1:3
    J(Interface_BC(ii,1)+1,:)=J(Interface_BC(ii,1)+1,:)+J(Interface_BC(ii,1),:);
end
for ii=4:6
    J(Interface_BC(ii,1)-1,:)=J(Interface_BC(ii,1)-1,:)+J(Interface_BC(ii,1),:);
end

for ii=1:size(Interface_BC,1)
    J(Interface_BC(ii,1),:)=0;
end
for ii=1:3
    J(Interface_BC(ii,1),Interface_BC(ii,1))=1;
    J(Interface_BC(ii,1),Interface_BC(ii,1)+1)=-1;
end
for ii=4:6
    J(Interface_BC(ii,1),Interface_BC(ii,1))=1;
    J(Interface_BC(ii,1),Interface_BC(ii,1)-1)=-1;
end

% Dirichlet_BC
for ii=1:size(Dirichlet_BC,1)
    J(Dirichlet_BC(ii,1),:)=0;
end
for ii=1:size(Dirichlet_BC,1)
    J(Dirichlet_BC(ii,1),Dirichlet_BC(ii,1))=1;
end

%res matrix 구하기
res=zeros(sum(Number_of_vertex),1);
for ii=1:size(Dirichlet_BC,1)
    res(Dirichlet_BC(ii,1),1)=Dirichlet_BC(ii,2);
end

phi=J\res;



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
title('Structure')

figure(6);
patch('Faces',Element(:,1:3),'Vertices',point(:,[1 2]), 'FaceVertexCData',phi, 'EdgeColor','black','FaceColor','interp');
title('Structure')
xlabel('X');
ylabel('Y');
colorbar
hold on
patch('Faces',element1,'Vertices',point, 'EdgeColor','Blue','FaceColor','none','LineWidth',2)
hold on
patch('Faces',element2,'Vertices',point, 'EdgeColor','Red','FaceColor','none','LineWidth',2)
hold on
patch('Faces',element3,'Vertices',point, 'EdgeColor','Yellow','FaceColor','none','LineWidth',2)
xlabel('X')
ylabel('Y')
title('Structure')
hold on
% patch('Faces',edge_interface,'Vertices',point, 'EdgeColor','Green','FaceColor','none','LineWidth',5)
hold off
% caxis([0 2]);