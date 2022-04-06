clear; close all; clc;

%Variable
eps1=3.9; eps2=11.7; eps3=3.9; % Epsilon
q=1.602192e-19; %Elementary charge,C
eps0=8.854187817e-12; %Vacuum permittivity
N_acc=1e24; % 1e18 / cm^3
k=1.3806488e-23; % Bolzmann constant, J/K
n_int=1.075e16; % 1.0e10 /cm^-3
T=300; % Temperture, K
V_T=k*T/q; % Thermal voltage

% Dirichlet_BC 불러오기
Dirichlet_BC_tmp = importdata("Dirichlet_BC.txt"); % Dirichlet_BC = [index potential]
m=0; n=0; % element로 이루어진 경계조건 계산하기 위해 변환
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

% Vertex 값 불러오기
Vertex_ox1 = importdata("Vertex_ox1.txt");
Vertex_si = importdata("Vertex_si.txt");
Vertex_ox2 = importdata("Vertex_ox2.txt");

% Indexing [x y z Region Index] 순서로 배치함
for ii=1:length(Vertex_ox1)
    Index_vertex_ox1(ii,:) = [Vertex_ox1(ii,:) 1 ii];
end
for ii=1:length(Vertex_si)
    Index_vertex_si(ii,:) = [Vertex_si(ii,:) 2 ii];
end
for ii=1:length(Vertex_ox2)
    Index_vertex_ox2(ii,:) = [Vertex_ox2(ii,:) 3 ii];
end
Index_vertex_tmp1=[Index_vertex_ox1; Index_vertex_si; Index_vertex_ox2];
Index_vertex_tmp2=sortrows(Index_vertex_tmp1, [2 1 3]);
Index_vertex=[Index_vertex_tmp2 transpose(1:length(Index_vertex_tmp2))];

% Vertex, index 확정
Vertex=Index_vertex(:,1:3); % [x y z]
Index=Index_vertex(:,4:6);  % [Region original-index re-index]

% Index_Region
n=1; % Region 1(ox1)
for ii=1:length(Index)
    if Index(ii,1)==1
        Index_Region_ox1(n,:)=Index(ii,:);
        n=n+1;
    end
end
n=1; % Region 2(si)
for ii=1:length(Index)
    if Index(ii,1)==2
        Index_Region_si(n,:)=Index(ii,:);
        n=n+1;
    end
end
n=1; % Region 3(ox2)
for ii=1:length(Index)
    if Index(ii,1)==3
        Index_Region_ox2(n,:)=Index(ii,:);
        n=n+1;
    end
end

% element 값 입력
Element_ox1_tmp = importdata("Element_ox1.txt"); Element_ox1=Element_ox1_tmp;
Element_si_tmp = importdata("Element_si.txt"); Element_si=Element_si_tmp;
Element_ox2_tmp = importdata("Element_ox2.txt"); Element_ox2=Element_ox2_tmp;
Element_tmp=[Element_ox1_tmp; Element_si_tmp; Element_ox2_tmp];

% re-element 작성
for ii=1:length(Element_ox1) % Region 1(ox1)
    for j=1:3
        tmp=Element_ox1_tmp(ii,j);
        Element_ox1(ii,j)=Index_Region_ox1(tmp,3);
    end
end
for ii=1:length(Element_si) % Region 2(si)
    for j=1:3
        tmp=Element_si_tmp(ii,j);
        Element_si(ii,j)=Index_Region_si(tmp,3);
    end
end
for ii=1:length(Element_ox2) % Region 3(ox2)
    for j=1:3
        tmp=Element_ox2_tmp(ii,j);
        Element_ox2(ii,j)=Index_Region_ox2(tmp,3);
    end
end

% Element 확정
eps1=ones(length(Element_ox1),1).*eps1;
eps2=ones(length(Element_si),1).*eps2;
eps3=ones(length(Element_ox2),1).*eps3;
Element = [Element_ox1 eps1; Element_si eps2; Element_ox2 eps3];

%%%% 각 Region의 edge 정렬 %%%%
n=1;
for ii=1:length(Element_ox1) % Si region
   for j=1:3
      edge_1(n,:)= [Element_ox1(ii,2) Element_ox1(ii,1)];
      edge_1(n+1,:)= [Element_ox1(ii,3) Element_ox1(ii,2)];
      edge_1(n+2,:)= [Element_ox1(ii,1) Element_ox1(ii,3)];
      n=n+3;
   end
end
n=1;
for ii=1:length(Element_si) % ox region
   for j=1:3
      edge_2(n,:)= [Element_si(ii,1) Element_si(ii,2)];
      edge_2(n+1,:)= [Element_si(ii,2) Element_si(ii,3)];
      edge_2(n+2,:)= [Element_si(ii,3) Element_si(ii,1)];
      n=n+3;
   end
end
n=1;
for ii=1:length(Element_ox2) % Si region
   for j=1:3
      edge_3(n,:)= [Element_ox2(ii,2) Element_ox2(ii,1)];
      edge_3(n+1,:)= [Element_ox2(ii,3) Element_ox2(ii,2)];
      edge_3(n+2,:)= [Element_ox2(ii,1) Element_ox2(ii,3)];
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

% 각 region에서 동일한 edge가 있는지 확인, Interface
edge_interface1 = intersect(edge_1,edge_2,'rows');
edge_interface2 = intersect(edge_1,edge_3,'rows');
edge_interface3 = intersect(edge_2,edge_3,'rows');
edge_interface=[edge_interface1; edge_interface2; edge_interface3];

% 각 region을 구성하는 vertex의 index를 확인하는 과정
Vertex_region1 = unique(sort(reshape(Element_ox1,[],1))); 
Vertex_region2 = unique(sort(reshape(Element_si,[],1)));
Vertex_region3 = unique(sort(reshape(Element_ox2,[],1)));

% 각 Reigon 내의 vertex 갯수 저장
Number_of_vertex=[length(Vertex_region1); length(Vertex_region2); length(Vertex_region3)];

%% 초기값 구하기
%각 벡터 구하기
vec=cell(length(Element),3);
for ii=1:length(Element)
    vec{ii,1}=Vertex(Element(ii,2),:)-Vertex(Element(ii,1),:);
    vec{ii,2}=Vertex(Element(ii,3),:)-Vertex(Element(ii,2),:);
    vec{ii,3}=Vertex(Element(ii,1),:)-Vertex(Element(ii,3),:);
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

% A matrix 구하기
AA=zeros(length(Vertex),length(Vertex));
for ii=1:length(Element) % 일단 모든 점에서의 자코비안 행렬 구하기
    v=[(-A(ii,1)/L(ii,1)-A(ii,3)/L(ii,3)) A(ii,1)/L(ii,1) A(ii,3)/L(ii,3);
        A(ii,1)/L(ii,1) -A(ii,1)/L(ii,1)-A(ii,2)/L(ii,2) A(ii,2)/L(ii,2);
        A(ii,3)/L(ii,3) A(ii,2)/L(ii,2) -A(ii,2)/L(ii,2)-A(ii,3)/L(ii,3)];
    v=Element(ii,4)*v;
    for j=1:3
        for k=1:3
        AA(Element(ii,j),Element(ii,k))=AA(Element(ii,j),Element(ii,k))+v(j,k);
        end
    end
end

% Interface_BC
Interface_BC=[3; 16; 29; 11; 24; 37];
for ii=1:3
    AA(Interface_BC(ii,1)+1,:)=AA(Interface_BC(ii,1)+1,:)+AA(Interface_BC(ii,1),:);
end
for ii=4:6
    AA(Interface_BC(ii,1)-1,:)=AA(Interface_BC(ii,1)-1,:)+AA(Interface_BC(ii,1),:);
end

for ii=1:size(Interface_BC,1)
    AA(Interface_BC(ii,1),:)=0;
end
for ii=1:3
    AA(Interface_BC(ii,1),Interface_BC(ii,1))=1;
    AA(Interface_BC(ii,1),Interface_BC(ii,1)+1)=-1;
end
for ii=4:6
    AA(Interface_BC(ii,1),Interface_BC(ii,1))=1;
    AA(Interface_BC(ii,1),Interface_BC(ii,1)-1)=-1;
end

% Dirichlet_BC
for ii=1:size(Dirichlet_BC,1)
    AA(Dirichlet_BC(ii,1),:)=0;
end
for ii=1:size(Dirichlet_BC,1)
    AA(Dirichlet_BC(ii,1),Dirichlet_BC(ii,1))=1;
end

%res matrix 구하기
b=zeros(length(Vertex),1);
for ii=1:size(Dirichlet_BC,1)
    b(Dirichlet_BC(ii,1),1)=Dirichlet_BC(ii,2);
end

phi=AA\b;

%Electron concentration
% n_elc=zeros(N,1);
for ii=1:length(Index_Region_si)
    n_elc0(ii,1)=n_int*exp(phi(Index_Region_si(ii,3),1)/V_T); % /m^3
end
%Hole concentration
% n_elc=zeros(N,1);
for ii=1:length(Index_Region_si)
    p_hole0(ii,1)=n_int*exp(-phi(Index_Region_si(ii,3),1)/V_T); % /m^3
end

%% Newton Method
% Variable, region 정보 저장
Variable_region=[1 1 1; % potential is in Region1,2,3(all) 
                 0 1 0; % electron is in Region2(si) 
                 0 1 0]; % hole is in Region2(si) 
% row: Variable. ex)aaa,bbb,ccc), columm: Region

% Solution vector x 생성
potential=NaN(sum(Number_of_vertex),1);
electron=NaN(sum(Number_of_vertex),1);
hole=NaN(sum(Number_of_vertex),1);
%potential
potential=phi; 
%electron
for ii=Number_of_vertex(1,1)+1:Number_of_vertex(1,1)+Number_of_vertex(2,1)
    j=ii-Number_of_vertex(1,1);
    electron(ii,1)=n_elc0(j,1);
end
%hole
for ii=Number_of_vertex(1,1)+1:Number_of_vertex(1,1)+Number_of_vertex(2,1)
    j=ii-Number_of_vertex(1,1);
    hole(ii,1)=p_hole0(j,1);
end
update_vector_tmp=[potential electron hole];

%x를 행백터로 변환
n=1;
for ii=1:sum(Number_of_vertex)
    for j=1:size(Variable_region,1)
        update_vector0(n,1)=update_vector_tmp(ii,j);
        n=n+1;
    end
end
%Nan값 제거
update_vector0=update_vector0(~isnan(update_vector0));

X_Vector_Newton=update_vector0;

% Jacobian 작성
% ox1
Jaco_ox1=zeros(length(update_vector0),length(update_vector0));
for ii=1:length(Element_ox1) % 일단 모든 점에서의 자코비안 행렬 구하기
    v=[(-A(ii,1)/(2*L(ii,1))-A(ii,3)/(2*L(ii,3))) A(ii,1)/(2*L(ii,1)) A(ii,3)/(2*L(ii,3));
        A(ii,1)/(2*L(ii,1)) -A(ii,1)/(2*L(ii,1))-A(ii,2)/(2*L(ii,2)) A(ii,2)/(2*L(ii,2));
        A(ii,3)/(2*L(ii,3)) A(ii,2)/(2*L(ii,2)) -A(ii,2)/(2*L(ii,2))-A(ii,3)/(2*L(ii,3))];
    v=Element(ii,4)*v;
    for j=1:3
        for k=1:3
        Jaco_ox1(Element(ii,j),Element(ii,k))=Jaco_ox1(Element(ii,j),Element(ii,k))+v(j,k);
        end
    end
end

% si
Jaco_si=zeros(length(update_vector0),length(update_vector0));
for ii=length(Element_ox1)+1:length(Element_ox1)+length(Element_si) % 일단 모든 점에서의 자코비안 행렬 구하기
    v=[(-A(ii,1)/(2*L(ii,1))-A(ii,3)/(2*L(ii,3))) A(ii,1)/(2*L(ii,1)) A(ii,3)/(2*L(ii,3));
        A(ii,1)/(2*L(ii,1)) -A(ii,1)/(2*L(ii,1))-A(ii,2)/(2*L(ii,2)) A(ii,2)/(2*L(ii,2));
        A(ii,3)/(2*L(ii,3)) A(ii,2)/(2*L(ii,2)) -A(ii,2)/(2*L(ii,2))-A(ii,3)/(2*L(ii,3))];
    v=Element(ii,4)*v;
    for j=1:3
        for k=1:3
        Jaco_si(Element(ii,j),Element(ii,k))=Jaco_si(Element(ii,j),Element(ii,k))+v(j,k);
        end
    end
end

% ox2
Jaco_ox2=zeros(length(update_vector0),length(update_vector0));
for ii=length(Element_ox1)+length(Element_si)+1:length(Element_ox1)+length(Element_si)+length(Element_ox2) % 일단 모든 점에서의 자코비안 행렬 구하기
    v=[(-A(ii,1)/(2*L(ii,1))-A(ii,3)/(2*L(ii,3))) A(ii,1)/(2*L(ii,1)) A(ii,3)/(2*L(ii,3));
        A(ii,1)/(2*L(ii,1)) -A(ii,1)/(2*L(ii,1))-A(ii,2)/(2*L(ii,2)) A(ii,2)/(2*L(ii,2));
        A(ii,3)/(2*L(ii,3)) A(ii,2)/(2*L(ii,2)) -A(ii,2)/(2*L(ii,2))-A(ii,3)/(2*L(ii,3))];
    v=Element(ii,4)*v;
    for j=1:3
        for k=1:3
        Jaco_ox2(Element(ii,j),Element(ii,k))=Jaco_ox2(Element(ii,j),Element(ii,k))+v(j,k);
        end
    end
end
% % Visualize 
% 
% figure(4); plot3(Vertex(:,1),Vertex(:,2),phi,'*')
% xlabel("x")
% ylabel("y")
% zlabel("potential")
% axis equal;
% 
% figure(5);
% patch('Faces',Element(:,1:3),'Vertices',Vertex(:,[1 2]), 'EdgeColor','black','FaceColor','none','LineWidth',0.5);
% xlabel('X');
% ylabel('Y');
% title('Structure')
% 
% figure(6);
% patch('Faces',Element(:,1:3),'Vertices',Vertex(:,[1 2]), 'FaceVertexCData',phi, 'EdgeColor','black','FaceColor','interp');
% title('Structure')
% xlabel('X');
% ylabel('Y');
% colorbar
% hold on
% patch('Faces',Element_ox1,'Vertices',Vertex, 'EdgeColor','Blue','FaceColor','none','LineWidth',2)
% hold on
% patch('Faces',Element_si,'Vertices',Vertex, 'EdgeColor','Red','FaceColor','none','LineWidth',2)
% hold on
% patch('Faces',Element_ox2,'Vertices',Vertex, 'EdgeColor','Yellow','FaceColor','none','LineWidth',2)
% xlabel('X')
% ylabel('Y')
% title('Structure')
% hold on
% % patch('Faces',edge_interface,'Vertices',point, 'EdgeColor','Green','FaceColor','none','LineWidth',5)
% hold off
% % caxis([0 2]);