clear; close all; clc;
if exist("data","dir")==0
    mkdir('./data')
end

%% Variable %%
eps=1; % 기존 코드 유지용

%% Build a Structure %%
% Vertex 값 작성
Vertex=vertex_doublegate(0,0, 1200,100, 50, 10);    % vertex = vertex_doublegate(x1, y1, x2, y2, dx, dy), [nm]

% anode_BC 불러오기
anode_BC_tmp=Contact(0,0,0,100,0);
m=0; n=0; % element로 이루어진 경계조건 계산하기 위해 변환
for ii=1:size(anode_BC_tmp,1)
    for j=size(anode_BC_tmp,2):-1:1
        if ~isnan(anode_BC_tmp(ii,j))
            for k=1:j-1
                anode_BC(k+n,:)=[anode_BC_tmp(ii,k) anode_BC_tmp(ii,j)];
                m=m+1;
            end
            n=m;
            break;
        end
    end
end

for ii=1:size(anode_BC_tmp,1)
    Element_anode_BC(ii,:)=anode_BC_tmp(ii,1:size(anode_BC_tmp,2)-1);
end


cathode_BC_tmp=Contact(1200,0, 1200,100, 0);
m=0; n=0; % element로 이루어진 경계조건 계산하기 위해 변환
for ii=1:size(cathode_BC_tmp,1)
    for j=size(cathode_BC_tmp,2):-1:1
        if ~isnan(cathode_BC_tmp(ii,j))
            for k=1:j-1
                cathode_BC(k+n,:)=[cathode_BC_tmp(ii,k) cathode_BC_tmp(ii,j)];
                m=m+1;
            end
            n=m;
            break;
        end
    end
end

for ii=1:size(cathode_BC_tmp,1)
    Element_cathode_BC(ii,:)=cathode_BC_tmp(ii,1:size(cathode_BC_tmp,2)-1);
end

% Dirac_function_BC=Contact(600,0, 600, 100, 1);
% m=0; n=0; % element로 이루어진 경계조건 계산하기 위해 변환
% for ii=1:size(Dirac_function_BC,1)
%     for j=size(Dirac_function_BC,2):-1:1
%         if ~isnan(Dirac_function_BC(ii,j))
%             for k=1:j-1
%                 Dirac_function_BC(k+n,:)=[Dirac_function_BC(ii,k) Dirac_function_BC(ii,j)];
%                 m=m+1;
%             end
%             n=m;
%             break;
%         end
%     end
% end
% 
% for ii=1:size(Dirac_function_BC,1)
%     Element_cathode_BC(ii,:)=Dirac_function_BC(ii,1:size(Dirac_function_BC,2)-1);
% end

% Element 값 입력
Element_si = element(0,0, 1200,100); % Element = element(x1, y1, x2, y2), [nm]
Element_cell={Element_si}; %% Region 변경시 수정필요!!
Element = vertcat(Element_cell{:}); % element 한 행렬로 합침.

for ii=1:size(Element_cell,1)
    nElmenet(ii,1) = length(Element_cell{ii,1});
end

% Element region table : 해당 element가 어떤 region에 속해있는지 알려주는 정보
count1=1;
for ii=1:size(nElmenet,1)
    count2=1;
    for j=1:nElmenet(ii,1)
        Table_element_region(count1,:)=[ii count2];
        count1=count1+1; count2=count2+1;
    end
end

% Vertex_region : 각 region을 내 vertex 확인, Ndop에 따라 초기 potential값 정해줄 건데 이때 필요.
for ii=1:size(nElmenet,1)
    Vertex_in_region{ii,1}=unique(sort(reshape(Element_cell{ii,1},[],1)));
end

% 각 Region 내 vertex 숫자.
for ii=1:size(nElmenet,1)
    Number_of_vertex(ii,1)=length(Vertex_in_region{ii,1});
end

Vertex_region_tmp=vertcat(Vertex_in_region{:});
Table_Vertex_region=zeros(max(Number_of_vertex),length(Number_of_vertex));
n=1;
for j=1:size(nElmenet,1)
    m=1;
    for ii=n:Number_of_vertex(j,1)+(n-1)
        Table_Vertex_region(m,j)=Vertex_region_tmp(ii,1);
        m=m+1;
    end
    n=n+Number_of_vertex(j,1);
end

%% Caculate the triangle mesh information %%
% 각 벡터 구하기
vec=cell(length(Element),3);
for ii=1:length(Element)
    vec{ii,1}=Vertex(Element(ii,2),:)-Vertex(Element(ii,1),:);
    vec{ii,2}=Vertex(Element(ii,3),:)-Vertex(Element(ii,2),:);
    vec{ii,3}=Vertex(Element(ii,1),:)-Vertex(Element(ii,3),:);
end

% 삼각형의 길이 구하기
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

% 외접의 반지름
R=zeros(length(Element),1);
for ii=1:length(Element)
    R(ii,1) = L(ii,1)*L(ii,2)*L(ii,3)/(4*Area(ii,1));
end

% 외심에서 각 변까지의 최단거리
edge=zeros(length(Element),3);
for ii=1:length(Element)
    edge_tmp2=[0.5*L(ii,1) 0.5*L(ii,2) 0.5*L(ii,3)];
    for j=1:3
        edge(ii,j) = real(sqrt(R(ii,1)^2-edge_tmp2(:,j)^2));
        if edge(ii,j)<=1e-15
            edge(ii,j)=0;
        end
    end
end

%% Jacobian Table 생성
n=1;
for nRegion=1:length(Number_of_vertex)
    for nVertex=1:Number_of_vertex(nRegion,1)
        for variable=1:3
            Table_Jaco(n,:)=[nRegion Table_Vertex_region(nVertex,nRegion) variable];
            n=n+1;
        end
    end
end

%% Green function
solution_vector_possion=zeros(length(Table_Jaco),length(Table_Jaco));
% Jacobian matrix / res vector

%     Jaco_possion=sparse(zeros(size(solution_vector_possion,1),size(solution_vector_possion,1))); res_possion=zeros(size(solution_vector_possion,1),1);
    Jaco_possion=zeros(size(solution_vector_possion,1),size(solution_vector_possion,1)); res_possion=zeros(size(solution_vector_possion,1),size(solution_vector_possion,1)); %% 희소행렬 아님. Jaco matrix 확인용.
    
% Potential %
    for ii=1:length(Element)
        % Control_Volume 계산
        Control_Volume=zeros(3,1);
        Control_Volume(1,1)=(0.25)*L(ii,3)*edge(ii,3)+(0.25)*L(ii,1)*edge(ii,1);
        Control_Volume(2,1)=(0.25)*L(ii,1)*edge(ii,1)+(0.25)*L(ii,2)*edge(ii,2);
        Control_Volume(3,1)=(0.25)*L(ii,2)*edge(ii,2)+(0.25)*L(ii,3)*edge(ii,3);

        % index : 현재 element에 속해있는 vertex가 resdue에 몇번째 있는가 저장
        index=zeros(9,1);
        n=1;
        for j=1:3
            index(n,1)=find(Table_Jaco(:,1)==Table_element_region(ii,1) & Table_Jaco(:,2)==Element(ii,j) & Table_Jaco(:,3)==1);
            index(n+1,1)=index(n,1)+1;
            index(n+2,1)=index(n,1)+2;
            n=n+3;
        end

        %% res vector
        res_possion=diag(ones(length(solution_vector_possion),1));

        %% Jaco matirx
        Jaco_tmp_possion=zeros(9,9);
        Jaco_tmp_possion(1,1)=eps(Table_element_region(ii,1),1)*(-edge(ii,1)/L(ii,1)-edge(ii,3)/L(ii,3));
%         Jaco_tmp_possion(1,2)=-coeff*Control_Volume(1,1);
%         Jaco_tmp_possion(1,3)=coeff*Control_Volume(1,1);
        Jaco_tmp_possion(1,4)=eps(Table_element_region(ii,1),1)*edge(ii,1)/L(ii,1);
        Jaco_tmp_possion(1,7)=eps(Table_element_region(ii,1),1)*edge(ii,3)/L(ii,3);

%         Jaco_tmp_possion(2,1)=-1/V_T*n_int*exp(solution_vector_possion(index(1,1),1)/V_T);
%         Jaco_tmp_possion(2,2)=1;

%         Jaco_tmp_possion(3,1)=1/V_T*n_int*exp(-solution_vector_possion(index(1,1),1)/V_T);
%         Jaco_tmp_possion(3,3)=1;

        Jaco_tmp_possion(4,1)=eps(Table_element_region(ii,1),1)*edge(ii,1)/L(ii,1);
        Jaco_tmp_possion(4,4)=eps(Table_element_region(ii,1),1)*(-edge(ii,1)/L(ii,1)-edge(ii,2)/L(ii,2));
%         Jaco_tmp_possion(4,5)=-coeff*Control_Volume(2,1);
%         Jaco_tmp_possion(4,6)=coeff*Control_Volume(2,1);
        Jaco_tmp_possion(4,7)=eps(Table_element_region(ii,1),1)*edge(ii,2)/L(ii,2);

%         Jaco_tmp_possion(5,4)=-1/V_T*n_int*exp(solution_vector_possion(index(4,1),1)/V_T);
%         Jaco_tmp_possion(5,5)=1;

%         Jaco_tmp_possion(6,4)=1/V_T*n_int*exp(-solution_vector_possion(index(4,1),1)/V_T);
%         Jaco_tmp_possion(6,6)=1;

        Jaco_tmp_possion(7,1)=eps(Table_element_region(ii,1),1)*edge(ii,3)/L(ii,3);
        Jaco_tmp_possion(7,4)=eps(Table_element_region(ii,1),1)*edge(ii,2)/L(ii,2);
        Jaco_tmp_possion(7,7)=eps(Table_element_region(ii,1),1)*(-edge(ii,2)/L(ii,2)-edge(ii,3)/L(ii,3));
%         Jaco_tmp_possion(7,8)=-coeff*Control_Volume(3,1);
%         Jaco_tmp_possion(7,9)=coeff*Control_Volume(3,1);

%         Jaco_tmp_possion(8,7)=-1/V_T*n_int*exp(solution_vector_possion(index(7,1),1)/V_T);
%         Jaco_tmp_possion(8,8)=1;

%         Jaco_tmp_possion(9,7)=1/V_T*n_int*exp(-solution_vector_possion(index(7,1),1)/V_T);
%         Jaco_tmp_possion(9,9)=1;

        for j=1:9
            for k=1:9
                Jaco_possion(index(j,1),index(k,1))=Jaco_possion(index(j,1),index(k,1))+Jaco_tmp_possion(j,k);
            end
        end
    end

    % cathode_BC
    for ii=1:size(anode_BC,1)
        BC_index=find(Table_Jaco(:,1)==1 & Table_Jaco(:,2)==anode_BC(ii,1) & Table_Jaco(:,3)==1);
        Jaco_possion(BC_index,:)=0; Jaco_possion(BC_index+1,:)=0; Jaco_possion(BC_index+2,:)=0;
        Jaco_possion(BC_index,BC_index)=1; Jaco_possion(BC_index+1,BC_index+1)=1; Jaco_possion(BC_index+2,BC_index+2)=1;
        res_possion(BC_index,:)=0; res_possion(BC_index+1,:)=0; res_possion(BC_index+2,:)=0;
        res_possion(BC_index,:)=solution_vector_possion(BC_index,1);
        res_possion(BC_index+1,:)=solution_vector_possion(BC_index+1,1);
        res_possion(BC_index+2,:)=solution_vector_possion(BC_index+2,1);
    end

    % cathode_BC
    for ii=1:size(cathode_BC,1)
        BC_index=find(Table_Jaco(:,1)==1 & Table_Jaco(:,2)==cathode_BC(ii,1) & Table_Jaco(:,3)==1);
        Jaco_possion(BC_index,:)=0; Jaco_possion(BC_index+1,:)=0; Jaco_possion(BC_index+2,:)=0;
        Jaco_possion(BC_index,BC_index)=1; Jaco_possion(BC_index+1,BC_index+1)=1; Jaco_possion(BC_index+2,BC_index+2)=1;
        res_possion(BC_index,:)=0; res_possion(BC_index+1,:)=0; res_possion(BC_index+2,:)=0;
        res_possion(BC_index,:)=solution_vector_possion(BC_index,1);
        res_possion(BC_index+1,:)=solution_vector_possion(BC_index+1,1);
        res_possion(BC_index+2,:)=solution_vector_possion(BC_index+2,1);
    end

%     % Scaling
%     Cvector_possion=zeros(size(solution_vector_possion,1),1); v=[V_T; max(abs(N_dop)); max(abs(N_dop))];
%     for ii=1:size(solution_vector_possion,1)
%         Cvector_possion(ii,1)=v(Table_Jaco(ii,3),1);
%     end
%     Jaco_possion=sparse(Jaco_possion);
%     Cmatrix_possion=spdiags(Cvector_possion,0,size(solution_vector_possion,1),size(solution_vector_possion,1));
%     Jaco_scaled_possion=Jaco_possion*Cmatrix_possion;
%     Rvector_possion=1./sum(abs(Jaco_scaled_possion),2);
%     Rmatrix_possion=spdiags(Rvector_possion,0,size(solution_vector_possion,1),size(solution_vector_possion,1));
%     Jaco_scaled_possion=Rmatrix_possion*Jaco_scaled_possion;
%     res_scaled_possion=Rmatrix_possion*res_possion;
%     update_scaled_possion=Jaco_scaled_possion\(-res_scaled_possion);
%     update_possion(:,1)=Cmatrix_possion*update_scaled_possion;
    
    res_possion=res_possion(1:3:length(res_possion),1:3:length(res_possion));
    Jaco_possion=Jaco_possion(1:3:length(Jaco_possion),1:3:length(Jaco_possion));

    solution_vector_possion=Jaco_possion\res_possion;

%% Save
save('./data/Homework19_green.mat')

%% Visualize
node=193;

Visual_solution_vector_possion=solution_vector_possion(:,node);

figure(1) % 점으로 3D
plot3(Vertex(:,1),Vertex(:,2),Visual_solution_vector_possion(:,1),'*')

figure(2); % mesh 모양
patch('Faces',Element,'Vertices',Vertex(:,[1 2]),'EdgeColor','black','FaceColor','none','LineWidth',1);
title('Mesh (1200nm*100nm)')
xlim([0 1.2*1e-6])
hold on
patch('Faces',Element_si,'Vertices',Vertex(:,[1 2]), 'EdgeColor','black','FaceColor','blue','LineWidth',0.5)
hold on
patch('Faces',Element_anode_BC,'Vertices',Vertex(:,[1 2]), 'EdgeColor','m','FaceColor','none','LineWidth',5)
hold on
patch('Faces',Element_cathode_BC,'Vertices',Vertex(:,[1 2]), 'EdgeColor','g','FaceColor','none','LineWidth',5)
hold off

figure(3); % potential
patch('Faces',Element,'Vertices',Vertex(:,[1 2]), 'FaceVertexCData',Visual_solution_vector_possion(:,1), 'EdgeColor','black','FaceColor','interp');
title('Initial potential')
hold on
patch('Faces',Element_anode_BC,'Vertices',Vertex(:,[1 2]), 'EdgeColor','m','FaceColor','none','LineWidth',5)
hold on
patch('Faces',Element_cathode_BC,'Vertices',Vertex(:,[1 2]), 'EdgeColor','g','FaceColor','none','LineWidth',5)
hold off
xlabel('X');
ylabel('Y');
xlim([0 1.2*1e-6])
colorbar