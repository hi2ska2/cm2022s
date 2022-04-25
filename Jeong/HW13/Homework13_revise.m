clear; close all; clc;

%% Variable %%
eps=[3.9; 11.7; 11.7; 11.7; 3.9]; % Epsilon=[eps1; eps2; eps3]
q=1.602192e-19; %Elementary charge,C
eps0=8.854187817e-12; %Vacuum permittivity
k=1.3806488e-23; % Bolzmann constant, J/K
n_int=1e16; % 1.0e10 /cm^-3
T=300; % Temperture, K
V_T=k*T/q; % Thermal voltage
V_gate=0.33374; % Voltage
Na=2e21; % /m^3
Nd=1e26; % /m^3
N_dop=[0; Nd; -Na; Nd; 0]; % -1e18 /cm^3, 행을 region별 도핑 농도
nm=10^-9; %길이단위
mobility_n = 1417e-4;  % electron mobility
mobility_p = 470.5e-4;     % Hole mobility

%% Build a Structure %%
% Vertex 값 작성
Vertex=vertex_doublegate(0,0, 30,10, 0.5, 0.5);    % vertex = vertex_doublegate(x1, y1, x2, y2, dx, dy), [nm]

% Dirichlet_BC 불러오기
% Dirichlet_BC_tmp = importdata("Dirichlet_BC.txt"); % Dirichlet_BC = [index potential]
Dirichlet_BC_tmp=[Contact(6,0,24,0,0); Contact(6,10,24,10,0)];
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

% Source_datin_BC_tmp = importdata("Source_datin_BC.txt"); % Dirichlet_BC = [index potential]
Source_darin_BC_tmp=[Contact(0,2, 0,8, V_T*log(Nd/n_int)); Contact(30,2, 30,8, V_T*log(Nd/n_int))];
m=0; n=0; % element로 이루어진 경계조건 계산하기 위해 변환
for ii=1:size(Source_darin_BC_tmp,1)
    for j=size(Source_darin_BC_tmp,2):-1:1
        if ~isnan(Source_darin_BC_tmp(ii,j))
            for k=1:j-1
                Source_drain_BC(k+n,:)=[Source_darin_BC_tmp(ii,k) Source_darin_BC_tmp(ii,j)];
                m=m+1;
            end
            n=m;
            break;
        end
    end
end

Drain_BC_tmp=[Contact(30,2, 30,8, V_T*log(Nd/n_int))];
m=0; n=0; % element로 이루어진 경계조건 계산하기 위해 변환
for ii=1:size(Drain_BC_tmp,1)
    for j=size(Drain_BC_tmp,2):-1:1
        if ~isnan(Drain_BC_tmp(ii,j))
            for k=1:j-1
                Drain_BC(k+n,:)=[Drain_BC_tmp(ii,k) Drain_BC_tmp(ii,j)];
                m=m+1;
            end
            n=m;
            break;
        end
    end
end

% Element 값 입력
Element_ox1 = element(0,0, 30,2); % Element = element(x1, y1, x2, y2), [nm]
Element_si_source = element(0,2, 6,8);
Element_si_channel = element(6,2, 24,8);
Element_si_drain = element(24,2, 30,8);
Element_ox2 = element(0,8, 30,10);
Element_cell={Element_ox1; Element_si_source; Element_si_channel; Element_si_drain; Element_ox2}; %% Region 변경시 수정필요!!
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

%% Find the Interface %%
% step 1. 각 Region의 edge 정렬 : element를 돌면서 edge정보를 얻은 후 중복되는 edge를 제거한다. 각 region별로 정리된 edge를 cell 내부에 저장한다.
for ii=1:size(nElmenet,1)
    n=1; clearvars Element_tmp edge_tmp1;  % 변수 초기화
    Element_tmp=Element_cell{ii,1};
    for j=1:nElmenet(ii,1)
        edge_tmp1(n,:)= [Element_tmp(j,2) Element_tmp(j,1)];
        edge_tmp1(n+1,:)= [Element_tmp(j,3) Element_tmp(j,2)];
        edge_tmp1(n+2,:)= [Element_tmp(j,1) Element_tmp(j,3)];
        n=n+3;
    end
    edge_cell{ii,1}=unique(sort(edge_tmp1, 2),'rows');
end

% Interface 찾기 : 각 egde별로 돌아가면서 겹치는 edge가 있는지 확인한다. 있다면 그것이 interface이다.
for ii=1:size(nElmenet,1)-1    % 1-1 1-2 1-3 1-4 ... 3-4 3-5 4-5 이렇게 edge를 확인하기 위한 for문
    for j=ii+1:size(nElmenet,1)
        edge_interface_tmp{ii,j}=intersect(edge_cell{ii,1},edge_cell{j,1},'rows');
    end
end

% 각 Interface 내부에 있는 vertex 정보 저장
for ii=1:size(nElmenet,1)-1
    for j=ii+1:size(nElmenet,1)
        Vertex_interface{ii,j}=unique(sort(reshape(edge_interface_tmp{ii,j},[],1)));
        Vertex_interface{j,ii}=Vertex_interface{ii,j};
    end
end

% 따로 떨어져 있는 edgd_interface를 한 행렬로 합쳐주는 과정 -> visulaize 할때 interface edge로 활용 가능
edge_interface=vertcat(edge_interface_tmp{:});

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
    end
end

%% Set a initial value %%
% Initial value table 작성(potential 값만 포함)
n=1;
for nRegion=1:length(Number_of_vertex)
    for nVertex=1:Number_of_vertex(nRegion,1)
        Table_initial(n,:)=[nRegion Table_Vertex_region(nVertex,nRegion)];
        n=n+1;
    end
end

% n-type: V_T*ln(Nd/n_int), p-type: -V_T*ln(Nd/n_int) 사용
phi=zeros(size(Table_initial,1),1);
for ii=1:size(Table_initial,1)
    if Table_initial(ii,1)==1 || Table_initial(ii,1)==5
        phi(ii,1)=0;
    elseif Table_initial(ii,1)==2 || Table_initial(ii,1)==4 % n-type
        phi(ii,1)=V_T*log(Nd/n_int);
    else
        phi(ii,1)=V_T*log(Na/n_int);
    end
end

% Chanel 의 값과 Source/Drain의 값을 동일하게 맞춰주기.
for ii=3
    for j=[2 4]
        clearvars Vertex_interface_tmp
        Vertex_interface_tmp=Vertex_interface{ii,j};
        for k=1:size(Vertex_interface_tmp,1)
            index_Channel=find(Table_initial(:,1)==ii & Table_initial(:,2)==Vertex_interface_tmp(k,1));
            index_Source_drain=find(Table_initial(:,1)==j & Table_initial(:,2)==Vertex_interface_tmp(k,1));

            mean_phi=(phi(index_Source_drain,1)+phi(index_Channel))/2;
            phi(index_Source_drain,1)=mean_phi;
            phi(index_Channel)=mean_phi;
        end
    end
end

% si 의 값과 ox의 값을 동일하게 맞춰주기.
for ii=2:4
    for j=[1 5]
        clearvars Vertex_interface_tmp
        Vertex_interface_tmp=Vertex_interface{ii,j};
        for k=1:size(Vertex_interface_tmp,1)
            index_si=find(Table_initial(:,1)==ii & Table_initial(:,2)==Vertex_interface_tmp(k,1));
            index_ox=find(Table_initial(:,1)==j & Table_initial(:,2)==Vertex_interface_tmp(k,1));

            phi(index_ox,1)=phi(index_si);
        end
    end
end

% Dirichlet_BC
for ii=1:size(Dirichlet_BC,1)
    for j=1:size(Table_initial,1)
        if Table_initial(j,2)==Dirichlet_BC(ii,1)
            phi(j,1)=Dirichlet_BC(ii,2)+V_gate;
        end
    end
end

% Source_darin_BC
for ii=1:size(Source_drain_BC,1)
    for j=1:size(Table_initial,1)
        if (Table_initial(j,1)==2 || Table_initial(j,1)==4) && Table_initial(j,2)==Source_drain_BC(ii,1)
            phi(j,1)=Source_drain_BC(ii,2);
        end
        if (Table_initial(j,1)==1 || Table_initial(j,1)==5) && Table_initial(j,2)==Source_drain_BC(ii,1)
            phi(j,1)=Source_drain_BC(ii,2);
        end
    end
end

% source, channel, drain 영역에서 initial carrrier density를 구해줌.
n_0=zeros(size(Table_initial,1),1); p_0=zeros(size(Table_initial,1),1);
for ii=1:size(Table_initial,1)
    if Table_initial(ii,1)==2 || Table_initial(ii,1)==3 || Table_initial(ii,1)==4
        n_0(ii,1)=n_int*exp(phi(ii,1)/V_T); % /m^3
        p_0(ii,1)=n_int*exp(-phi(ii,1)/V_T); % /m^3
    end
end

% Initial value
initial_value=[phi n_0 p_0];

n=1;
for ii=1:size(Table_initial,1)
    if Table_initial(ii,1)==1 || Table_initial(ii,1)==5
        solution_vector_possion(n,1)=initial_value(ii,1);
        n=n+1;
    else
        for j=1:3
            solution_vector_possion(n,1)=initial_value(ii,j);
            n=n+1;
        end
    end
end


% % Visualize
%
% % Jacobian Table 생성
% n=1;
% for nRegion=1:length(Number_of_vertex)
%     for nVertex=1:Number_of_vertex(nRegion,1)
%         if nRegion==2 || nRegion==3 || nRegion==4 %% Region 변하면 수정
%             for variable=1:3
%                 Table_Jaco(n,:)=[nRegion Table_Vertex_region(nVertex,nRegion) variable];
%                 n=n+1;
%             end
%         elseif nRegion==1 || nRegion==5
%             Table_Jaco(n,:)=[nRegion Table_Vertex_region(nVertex,nRegion) 1];
%             n=n+1;
%         end
%     end
% end
%
% Visual_solution_vector=zeros(size(Vertex,1),3);
% Visual_initial_vector=zeros(size(Vertex,1),3);
% for ii=1:size(Table_Jaco,1)
%     Visual_initial_vector(Table_Jaco(ii,2),Table_Jaco(ii,3))=solution_vector(ii,1);
% end
%
% figure(1);
% patch('Faces',Element,'Vertices',Vertex(:,[1 2]), 'FaceVertexCData',Visual_solution_vector(:,1), 'EdgeColor','black','FaceColor','interp');
% title('Initial potential')
% xlabel('X');
% ylabel('Y');
% colorbar
%
% figure(2);
% patch('Faces',Element,'Vertices',Vertex(:,[1 2]), 'FaceVertexCData',Visual_solution_vector(:,2), 'EdgeColor','black','FaceColor','interp');
% title('elctron')
% xlabel('X');
% ylabel('Y');
% colorbar
%
% figure(3);
% patch('Faces',Element,'Vertices',Vertex(:,[1 2]), 'FaceVertexCData',Visual_solution_vector(:,3), 'EdgeColor','black','FaceColor','interp');
% title('hole')
% xlabel('X');
% ylabel('Y');
% colorbar
%% Newton Method
% Jacobian Table 생성
n=1;
for nRegion=1:length(Number_of_vertex)
    for nVertex=1:Number_of_vertex(nRegion,1)
        if nRegion==2 || nRegion==3 || nRegion==4 %% Region 변하면 수정
            for variable=1:3
                Table_Jaco(n,:)=[nRegion Table_Vertex_region(nVertex,nRegion) variable];
                n=n+1;
            end
        elseif nRegion==1 || nRegion==5
            Table_Jaco(n,:)=[nRegion Table_Vertex_region(nVertex,nRegion) 1];
            n=n+1;
        end
    end
end

% Jacobian matrix / res vector
for Newton_possion=1:20
    %     Jaco=sparse(zeros(size(solution_vector_possion,1),size(solution_vector_possion,1))); res=zeros(size(solution_vector_possion,1),1);
    Jaco=zeros(size(solution_vector_possion,1),size(solution_vector_possion,1)); res=zeros(size(solution_vector_possion,1),1); %% 희소행렬 아님. Jaco matrix 확인용.
    % Potential
    for ii=1:length(Element)
        % Control_Volume 계산
        Control_Volume=zeros(3,1);
        Control_Volume=[(0.25)*L(ii,3)*edge(ii,3)+(0.25)*L(ii,1)*edge(ii,1);
            (0.25)*L(ii,1)*edge(ii,1)+(0.25)*L(ii,2)*edge(ii,2);
            (0.25)*L(ii,2)*edge(ii,2)+(0.25)*L(ii,3)*edge(ii,3)];
        coeff=q/eps0;     % m^2 -> nm^2으로 단위 변경, q/eps0 해줌

        if Table_element_region(ii,1)==2 || Table_element_region(ii,1)==3 || Table_element_region(ii,1)==4 %%% Region 변경되면 꼭 수정되어야함!!
            % index : 현재 element에 속해있는 vertex가 resdue에 몇번째 있는가 저장
            index=zeros(9,1);
            n=1;
            for j=1:3
                index(n,1)=find(Table_Jaco(:,1)==Table_element_region(ii,1) & Table_Jaco(:,2)==Element(ii,j) & Table_Jaco(:,3)==1);
                index(n+1,1)=index(n,1)+1;
                index(n+2,1)=index(n,1)+2;
                n=n+3;
            end
            % charge 계산
            for j=1:3
                %                 charge(j,1)=(q/eps0)*(N_dop(Table_element_region(ii,1),1)-solution_vector(index(3*j-1,1),Newton)+solution_vector(index(3*j,1),Newton));
            end
            Jaco_tmp=zeros(9,9);
            % Jaco_potential matrix
            Jaco_tmp=[eps(Table_element_region(ii,1),1)*(-edge(ii,1)/L(ii,1)-edge(ii,3)/L(ii,3)) -coeff*Control_Volume(1,1) coeff*Control_Volume(1,1) eps(Table_element_region(ii,1),1)*edge(ii,1)/L(ii,1) 0 0 eps(Table_element_region(ii,1),1)*edge(ii,3)/L(ii,3) 0 0;
                -1/V_T*n_int*exp(solution_vector_possion(index(1,1),Newton_possion)/V_T) 1 0 0 0 0 0 0 0;
                1/V_T*n_int*exp(-solution_vector_possion(index(1,1),Newton_possion)/V_T) 0 1 0 0 0 0 0 0;
                eps(Table_element_region(ii,1),1)*edge(ii,1)/L(ii,1) 0 0 eps(Table_element_region(ii,1),1)*(-edge(ii,1)/L(ii,1)-edge(ii,2)/L(ii,2)) -coeff*Control_Volume(2,1) coeff*Control_Volume(2,1) eps(Table_element_region(ii,1),1)*edge(ii,2)/L(ii,2) 0 0;
                0 0 0 -1/V_T*n_int*exp(solution_vector_possion(index(4,1),Newton_possion)/V_T) 1 0 0 0 0;
                0 0 0 1/V_T*n_int*exp(-solution_vector_possion(index(4,1),Newton_possion)/V_T) 0 1 0 0 0;
                eps(Table_element_region(ii,1),1)*edge(ii,3)/L(ii,3) 0 0 eps(Table_element_region(ii,1),1)*edge(ii,2)/L(ii,2) 0 0 eps(Table_element_region(ii,1),1)*(-edge(ii,2)/L(ii,2)-edge(ii,3)/L(ii,3)) -coeff*Control_Volume(3,1) coeff*Control_Volume(3,1);
                0 0 0 0 0 0 -1/V_T*n_int*exp(solution_vector_possion(index(7,1),Newton_possion)/V_T) 1 0;
                0 0 0 0 0 0 1/V_T*n_int*exp(-solution_vector_possion(index(7,1),Newton_possion)/V_T) 0 1];
            for j=1:9
                for k=1:9
                    Jaco(index(j,1),index(k,1))=Jaco(index(j,1),index(k,1))+Jaco_tmp(j,k);
                end
            end
            % n,p 값은 변하지 않는다.
            for j=1:3
                Jaco(index(3*j-1),:)=0;
                Jaco(index(3*j),:)=0;
                Jaco(index(3*j-1),index(3*j-2))=-1/V_T*n_int*exp(solution_vector_possion(index(3*j-2,1),Newton_possion)/V_T);
                Jaco(index(3*j-1),index(3*j-1))=1;
                Jaco(index(3*j),index(3*j-2))=1/V_T*n_int*exp(-solution_vector_possion(index(3*j-2,1),Newton_possion)/V_T);
                Jaco(index(3*j),index(3*j))=1;
            end

            % res vector
            res_tmp=zeros(9,1);
            %             res_potential_tmp1=eps(Table_element_region(ii,1),1).*[(-edge(ii,1)/L(ii,1)-edge(ii,3)/L(ii,3))*solution_vector(index(1,1),Newton)+edge(ii,1)/L(ii,1)*solution_vector(index(4,1),Newton)+edge(ii,3)/L(ii,3)*solution_vector(index(7,1),Newton);
            %                                                                      edge(ii,1)/L(ii,1)*solution_vector(index(1,1),Newton)+(-edge(ii,1)/L(ii,1)-edge(ii,2)/L(ii,2))*solution_vector(index(4,1),Newton)+edge(ii,2)/L(ii,2)*solution_vector(index(7,1),Newton);
            %                                                                      edge(ii,3)/L(ii,3)*solution_vector(index(1,1),Newton)+edge(ii,2)/L(ii,2)*solution_vector(index(4,1),Newton)+(-edge(ii,2)/L(ii,2)-edge(ii,3)/L(ii,3))*solution_vector(index(7,1),Newton)];

            res_potential_tmp1=eps(Table_element_region(ii,1),1).*[edge(ii,1)/L(ii,1)*(solution_vector_possion(index(4,1),Newton_possion)-solution_vector_possion(index(1,1),Newton_possion))+edge(ii,3)/L(ii,3)*(solution_vector_possion(index(7,1),Newton_possion)-solution_vector_possion(index(1,1),Newton_possion));
                edge(ii,2)/L(ii,2)*(solution_vector_possion(index(7,1),Newton_possion)-solution_vector_possion(index(4,1),Newton_possion))+edge(ii,1)/L(ii,1)*(solution_vector_possion(index(1,1),Newton_possion)-solution_vector_possion(index(4,1),Newton_possion));
                edge(ii,3)/L(ii,3)*(solution_vector_possion(index(1,1),Newton_possion)-solution_vector_possion(index(7,1),Newton_possion))+edge(ii,2)/L(ii,2)*(solution_vector_possion(index(4,1),Newton_possion)-solution_vector_possion(index(7,1),Newton_possion))];
            for j=1:3
                res_potential_tmp2(j,1)=(N_dop(Table_element_region(ii,1),1)-solution_vector_possion(index(3*j-1,1),Newton_possion)+solution_vector_possion(index(3*j,1),Newton_possion))*coeff*Control_Volume(j,1);
            end
            res_potential=res_potential_tmp1+res_potential_tmp2;
            for j=1:3
                res_electron(j,1)=solution_vector_possion(index(3*j-1,1),Newton_possion)-n_int*exp(solution_vector_possion(index(3*j-2,1),Newton_possion)/V_T);
            end
            for j=1:3
                res_hole(j,1)=solution_vector_possion(index(3*j,1),Newton_possion)-n_int*exp(-solution_vector_possion(index(3*j-2,1),Newton_possion)/V_T);
            end
            n=1;
            for j=1:3
                res_tmp(n,1)=res_potential(j,1);
                res_tmp(n+1,1)=res_electron(j,1);
                res_tmp(n+2,1)=res_hole(j,1);
                n=n+3;
            end

            for j=1:9
                res(index(j,1),1)=res(index(j,1),1)+res_tmp(j,1);
            end
            %n,p값은 변하지 않는다.
            for j=1:3
                res(index(3*j-1,1),1)=res_electron(j,1);
                res(index(3*j,1),1)=res_hole(j,1);
            end

        else
            % index : 현재 element에 속해있는 vertex가 resdue에 몇번째 있는가 저장
            index=zeros(3,1);
            n=1;
            for j=1:3
                index(j,1)=find(Table_Jaco(:,1)==Table_element_region(ii,1) & Table_Jaco(:,2)==Element(ii,j) & Table_Jaco(:,3)==1);
            end

            % Jaco_potential matrix
            Jaco_tmp=zeros(3,3);
            Jaco_tmp=[(-edge(ii,1)/L(ii,1)-edge(ii,3)/L(ii,3)) edge(ii,1)/L(ii,1) edge(ii,3)/L(ii,3);
                edge(ii,1)/L(ii,1) -edge(ii,1)/L(ii,1)-edge(ii,2)/L(ii,2) edge(ii,2)/L(ii,2);
                edge(ii,3)/L(ii,3) edge(ii,2)/L(ii,2) -edge(ii,2)/L(ii,2)-edge(ii,3)/L(ii,3)];
            Jaco_tmp=eps(Table_element_region(ii,1),1)*Jaco_tmp;
            for j=1:3
                for k=1:3
                    Jaco(index(j,1),index(k,1))=Jaco(index(j,1),index(k,1))+Jaco_tmp(j,k);
                end
            end

            % res vector
            res_tmp=zeros(3,1);
            res_tmp=eps(Table_element_region(ii,1),1).*[(-edge(ii,1)/L(ii,1)-edge(ii,3)/L(ii,3))*solution_vector_possion(index(1,1),Newton_possion)+edge(ii,1)/L(ii,1)*solution_vector_possion(index(2,1),Newton_possion)+edge(ii,3)/L(ii,3)*solution_vector_possion(index(3,1),Newton_possion);
                edge(ii,1)/L(ii,1)*solution_vector_possion(index(1,1),Newton_possion)+(-edge(ii,1)/L(ii,1)-edge(ii,2)/L(ii,2))*solution_vector_possion(index(2,1),Newton_possion)+edge(ii,2)/L(ii,2)*solution_vector_possion(index(3,1),Newton_possion);
                edge(ii,3)/L(ii,3)*solution_vector_possion(index(1,1),Newton_possion)+edge(ii,2)/L(ii,2)*solution_vector_possion(index(2,1),Newton_possion)+(-edge(ii,2)/L(ii,2)-edge(ii,3)/L(ii,3))*solution_vector_possion(index(3,1),Newton_possion)];
            for j=1:3
                res(index(j,1),1)=res(index(j,1),1)+res_tmp(j,1);
            end
        end
    end

    % ox 의 값을 si에 더해준 후 ox값은 1,-1로 변경하기
    for ii=[1 5]
        for j=[3 2 4]
            clearvars Vertex_interface_tmp
            Vertex_interface_tmp=Vertex_interface{ii,j};
            for k=1:size(Vertex_interface_tmp,1)
                index_ox=find(Table_Jaco(:,1)==ii & Table_Jaco(:,2)==Vertex_interface_tmp(k,1) & Table_Jaco(:,3)==1);
                index_si=find(Table_Jaco(:,1)==j & Table_Jaco(:,2)==Vertex_interface_tmp(k,1) & Table_Jaco(:,3)==1);
                if length(sort(nonzeros(Jaco(index_ox,:))))>2
                    Jaco(index_si,:)=Jaco(index_si,:)+Jaco(index_ox,:);
                    Jaco(index_ox,:)=0; Jaco(index_ox,index_ox)=1; Jaco(index_ox,index_si)=-1;
                    res(index_si,1)=res(index_si,1)+res(index_ox,1);
                    res(index_ox,1)=0;
                elseif ~isequal(sort(nonzeros(Jaco(index_ox,:))),[-1; 1]) % 점 3개가 접해있는 지역에서 중복적으로 합쳐지는것을 막기위해 추가함. jaco의 ox-channel행이 -1과1이라면 이미 다른행과 합쳐진 것이므로 pass.
                    Jaco(index_si,:)=Jaco(index_si,:)+Jaco(index_ox,:);
                    Jaco(index_ox,:)=0; Jaco(index_ox,index_ox)=1; Jaco(index_ox,index_si)=-1;
                    res(index_si,1)=res(index_si,1)+res(index_ox,1);
                    res(index_ox,1)=0;
                end
            end
        end
    end
    % source/drain 값을 channel에 더해주기. Si부분이니까 electron/hole density도 모두 바꾸어주어야함.
    for ii=[2 4]
        for j=3
            clearvars Vertex_interface_tmp
            Vertex_interface_tmp=Vertex_interface{ii,j};
            for k=1:size(Vertex_interface_tmp,1)
                % index 설정
                index_source_drain_potential=find(Table_Jaco(:,1)==ii & Table_Jaco(:,2)==Vertex_interface_tmp(k,1) & Table_Jaco(:,3)==1);
                index_source_drain_electron=index_source_drain_potential+1;
                index_source_drain_hole=index_source_drain_potential+2;
                index_channel_potential=find(Table_Jaco(:,1)==j & Table_Jaco(:,2)==Vertex_interface_tmp(k,1) & Table_Jaco(:,3)==1);
                index_channel_electron=index_channel_potential+1;
                index_channel_hole=index_channel_potential+2;
                % Jaco 변경, potential/electron/hole 순서
                Jaco(index_channel_potential,:)=Jaco(index_channel_potential,:)+Jaco(index_source_drain_potential,:);
                Jaco(index_channel_electron,:)=Jaco(index_channel_electron,:)+Jaco(index_source_drain_electron,:);
                Jaco(index_channel_hole,:)=Jaco(index_channel_hole,:)+Jaco(index_source_drain_hole,:);
                Jaco(index_source_drain_potential,:)=0; Jaco(index_source_drain_potential,index_source_drain_potential)=1; Jaco(index_source_drain_potential,index_channel_potential)=-1;
                Jaco(index_source_drain_electron,:)=0; Jaco(index_source_drain_electron,index_source_drain_electron)=1; Jaco(index_source_drain_electron,index_channel_electron)=-1;
                Jaco(index_source_drain_hole,:)=0; Jaco(index_source_drain_hole,index_source_drain_hole)=1; Jaco(index_source_drain_hole,index_channel_hole)=-1;

                % res 변경, potential/electron/hole 순서
                res(index_channel_potential,1)=res(index_channel_potential,1)+res(index_source_drain_potential,1);
                res(index_channel_electron,1)=res(index_channel_electron,1)+res(index_source_drain_electron,1);
                res(index_channel_hole,1)=res(index_channel_hole,1)+res(index_source_drain_hole,1);
                res(index_source_drain_potential,1)=0;
                res(index_source_drain_electron,1)=0;
                res(index_source_drain_hole,1)=0;
            end
        end
    end

    % Dirichlet_BC
    for ii=1:size(Dirichlet_BC,1)
        BC_index=find(Table_Jaco(:,2)==Dirichlet_BC(ii,1) & Table_Jaco(:,3)==1);
        Jaco(BC_index,:)=0;
        Jaco(BC_index,BC_index)=1;
        res(BC_index,1)=solution_vector_possion(BC_index,Newton_possion)-V_gate;
    end

    % Source_drain_BC
    for ii=1:size(Source_drain_BC,1)
        BC_index=find((Table_Jaco(:,1)==2 | Table_Jaco(:,1)==4) & Table_Jaco(:,2)==Source_drain_BC(ii,1) & Table_Jaco(:,3)==1);
        Jaco(BC_index,:)=0; Jaco(BC_index+1,:)=0; Jaco(BC_index+2,:)=0;
        Jaco(BC_index,BC_index)=1; Jaco(BC_index+1,BC_index+1)=1; Jaco(BC_index+2,BC_index+2)=1;
        res(BC_index,1)=solution_vector_possion(BC_index,Newton_possion)-V_T*log(Nd/n_int);
        res(BC_index+1,1)=solution_vector_possion(BC_index+1,Newton_possion)-Nd;
        res(BC_index+2,1)=solution_vector_possion(BC_index+2,Newton_possion)-n_int^2/Nd;
    end

    res_save(:,Newton_possion)=res(:,1);

    % Scaling
    Cvector=zeros(size(solution_vector_possion,1),1); v=[V_T; max(abs(N_dop)); max(abs(N_dop))];
    for ii=1:size(solution_vector_possion,1)
        Cvector(ii,1)=v(Table_Jaco(ii,3),1);
    end
    Cmatrix=spdiags(Cvector,0,size(solution_vector_possion,1),size(solution_vector_possion,1));
    Jaco_scaled=Jaco*Cmatrix;
    Rvector=1./sum(abs(Jaco_scaled),2);
    Rmatrix=spdiags(Rvector,0,size(solution_vector_possion,1),size(solution_vector_possion,1));
    Jaco_scaled=Rmatrix*Jaco_scaled;
    res_scaled=Rmatrix*res;
    update_scaled=Jaco_scaled\(-res_scaled);
    update(:,Newton_possion)=Cmatrix*update_scaled;
    solution_vector_possion(:,Newton_possion+1)=solution_vector_possion(:,Newton_possion)+update(:,Newton_possion);

    %     % non-Scaling
    %     update(:,Newton)=Jaco\-res;
    %     solution_vector(:,Newton+1)=solution_vector(:,Newton)+update(:,Newton);

end

% % Visualize
% Visual_solution_vector=zeros(size(Vertex,1),3);
% Visual_initial_vector=zeros(size(Vertex,1),3);
% for ii=1:size(Table_Jaco,1)
%     Visual_solution_vector(Table_Jaco(ii,2),Table_Jaco(ii,3))=solution_vector_possion(ii,Newton_possion);
%     Visual_initial_vector(Table_Jaco(ii,2),Table_Jaco(ii,3))=solution_vector_possion(ii,1);
% end
%
% dif=abs(Visual_initial_vector-Visual_solution_vector);
% plot3(Vertex(:,1),Vertex(:,2),Visual_solution_vector(:,1),'*')
% xlim([0 30*1e-9])
%
% figure(4);
% patch('Faces',Element,'Vertices',Vertex(:,[1 2]), 'FaceVertexCData',Visual_solution_vector(:,1), 'EdgeColor','black','FaceColor','interp');
% title('Initial potential')
% xlim([0 30*1e-9])
% xlabel('X');
% ylabel('Y');
% colorbar
%
% figure(5);
% patch('Faces',Element,'Vertices',Vertex(:,[1 2]), 'FaceVertexCData',Visual_solution_vector(:,2), 'EdgeColor','black','FaceColor','interp');
% title('elctron')
% xlim([0 30*1e-9])
% xlabel('X');
% ylabel('Y');
% colorbar
%
% figure(6);
% patch('Faces',Element,'Vertices',Vertex(:,[1 2]), 'FaceVertexCData',Visual_solution_vector(:,3), 'EdgeColor','black','FaceColor','interp');
% title('hole')
% xlim([0 30*1e-9])
% xlabel('X');
% ylabel('Y');
% colorbar
%
% % figure(7);
% % patch('Faces',Element,'Vertices',Vertex(:,[1 2]), 'FaceVertexCData',dif(:,1), 'EdgeColor','black','FaceColor','interp');
% % title('Initial potential')
% % xlabel('X');
% % ylabel('Y');
% % colorbar
% %
% % figure(8);
% % patch('Faces',Element,'Vertices',Vertex(:,[1 2]), 'FaceVertexCData',dif(:,2), 'EdgeColor','black','FaceColor','interp');
% % title('elctron')
% % xlabel('X');
% % ylabel('Y');
% % colorbar
% %
% % figure(9);
% % patch('Faces',Element,'Vertices',Vertex(:,[1 2]), 'FaceVertexCData',dif(:,3), 'EdgeColor','black','FaceColor','interp');
% % title('hole')
% % xlabel('X');
% % ylabel('Y');
% % colorbar

%% Drift-Diffusion Model
solution_vector_DD(:,1)=solution_vector_possion(:,Newton_possion);

% Jacobian matrix / res vector
for Newton_DD=1:10
    %     Jaco_DD=sparse(zeros(size(solution_vector_DD,1),size(solution_vector_DD,1))); res_DD=zeros(size(solution_vector_DD,1),1);
    Jaco_DD=zeros(size(solution_vector_DD,1),size(solution_vector_DD,1)); res_DD=zeros(size(solution_vector_DD,1),1); %% 희소행렬 아님. Jaco matrix 확인용.
    % Potential
    for ii=1:length(Element)
        % Control_Volume 계산
        Control_Volume=zeros(3,1);
        Control_Volume=[(0.25)*L(ii,3)*edge(ii,3)+(0.25)*L(ii,1)*edge(ii,1);
            (0.25)*L(ii,1)*edge(ii,1)+(0.25)*L(ii,2)*edge(ii,2);
            (0.25)*L(ii,2)*edge(ii,2)+(0.25)*L(ii,3)*edge(ii,3)];
        coeff=q/eps0;     % m^2 -> nm^2으로 단위 변경, q/eps0 해줌
        coeff_Jn=q*mobility_n*V_T;
        coeff_Jp=-q*mobility_p*V_T;

        if Table_element_region(ii,1)==2 || Table_element_region(ii,1)==3 || Table_element_region(ii,1)==4 %%% Region 변경되면 꼭 수정되어야함!!
            % index : 현재 element에 속해있는 vertex가 resdue에 몇번째 있는가 저장
            index=zeros(9,1);
            n=1;
            for j=1:3
                index(n,1)=find(Table_Jaco(:,1)==Table_element_region(ii,1) & Table_Jaco(:,2)==Element(ii,j) & Table_Jaco(:,3)==1);
                index(n+1,1)=index(n,1)+1;
                index(n+2,1)=index(n,1)+2;
                n=n+3;
            end

            % 식 간단히 표기하기 위해 미리 계산.
            x12=solution_vector_DD(index(1,1),1)-solution_vector_DD(index(4,1),1); x21=-x12;
            x23=solution_vector_DD(index(4,1),1)-solution_vector_DD(index(7,1),1); x32=-x23;
            x13=solution_vector_DD(index(1,1),1)-solution_vector_DD(index(7,1),1); x31=-x13;
            n1=solution_vector_DD(index(1,1)+1); n2=solution_vector_DD(index(4,1)+1); n3=solution_vector_DD(index(7,1)+1);
            p1=solution_vector_DD(index(1,1)+2); p2=solution_vector_DD(index(4,1)+2); p3=solution_vector_DD(index(7,1)+2);

            Jaco_tmp_DD=zeros(9,9);
            % Jaco_potential matrix
            %             Jaco_tmp=[eps(Table_element_region(ii,1),1)*(-edge(ii,1)/L(ii,1)-edge(ii,3)/L(ii,3)) -coeff*Control_Volume(1,1) coeff*Control_Volume(1,1) eps(Table_element_region(ii,1),1)*edge(ii,1)/L(ii,1) 0 0 eps(Table_element_region(ii,1),1)*edge(ii,3)/L(ii,3) 0 0;
            %                     -1/V_T*n_int*exp(solution_vector_DD(index(1,1),Newton_DD)/V_T) 1 0 0 0 0 0 0 0;
            %                     1/V_T*n_int*exp(-solution_vector_DD(index(1,1),Newton_DD)/V_T) 0 1 0 0 0 0 0 0;
            %                     eps(Table_element_region(ii,1),1)*edge(ii,1)/L(ii,1) 0 0 eps(Table_element_region(ii,1),1)*(-edge(ii,1)/L(ii,1)-edge(ii,2)/L(ii,2)) -coeff*Control_Volume(2,1) coeff*Control_Volume(2,1) eps(Table_element_region(ii,1),1)*edge(ii,2)/L(ii,2) 0 0;
            %                     0 0 0 -1/V_T*n_int*exp(solution_vector_DD(index(4,1),Newton_DD)/V_T) 1 0 0 0 0;
            %                     0 0 0 1/V_T*n_int*exp(-solution_vector_DD(index(4,1),Newton_DD)/V_T) 0 1 0 0 0;
            %                     eps(Table_element_region(ii,1),1)*edge(ii,3)/L(ii,3) 0 0 eps(Table_element_region(ii,1),1)*edge(ii,2)/L(ii,2) 0 0 eps(Table_element_region(ii,1),1)*(-edge(ii,2)/L(ii,2)-edge(ii,3)/L(ii,3)) -coeff*Control_Volume(3,1) coeff*Control_Volume(3,1);
            %                     0 0 0 0 0 0 -1/V_T*n_int*exp(solution_vector_DD(index(7,1),Newton_DD)/V_T) 1 0;
            %                     0 0 0 0 0 0 1/V_T*n_int*exp(-solution_vector_DD(index(7,1),Newton_DD)/V_T) 0 1];

            Jaco_tmp_DD(1,1)=eps(Table_element_region(ii,1),1)*(-edge(ii,1)/L(ii,1)-edge(ii,3)/L(ii,3));
            Jaco_tmp_DD(1,2)=-coeff*Control_Volume(1,1);
            Jaco_tmp_DD(1,3)=coeff*Control_Volume(1,1);
            Jaco_tmp_DD(1,4)=eps(Table_element_region(ii,1),1)*edge(ii,1)/L(ii,1);
            Jaco_tmp_DD(1,7)=eps(Table_element_region(ii,1),1)*edge(ii,3)/L(ii,3);

            Jaco_tmp_DD(2,1)=coeff_Jn*((edge(ii,1)/L(ii,1))*(n2*(-V_T)*Ber_d(x21/V_T)-n1*(V_T)*Ber_d(-x21/V_T))+((edge(ii,3)/L(ii,3))*(n3*(-V_T)*Ber_d(x31/V_T)-n1*(V_T)*Ber_d(-x31/V_T))));
            Jaco_tmp_DD(2,2)=coeff_Jn*((edge(ii,1)/L(ii,1))*(-Ber(-x21/V_T))+((edge(ii,3)/L(ii,3))*(-Ber(-x31/V_T))));
            Jaco_tmp_DD(2,4)=coeff_Jn*((edge(ii,1)/L(ii,1))*(n2*(V_T)*Ber_d(x21/V_T)-n1*(-V_T)*Ber_d(-x21/V_T)));
            Jaco_tmp_DD(2,5)=coeff_Jn*((edge(ii,1)/L(ii,1))*(Ber(x21/V_T)));
            Jaco_tmp_DD(2,7)=coeff_Jn*((edge(ii,3)/L(ii,3))*(n3*(V_T)*Ber_d(x31/V_T)-n1*(-V_T)*Ber_d(-x31/V_T)));
            Jaco_tmp_DD(2,8)=coeff_Jn*((edge(ii,3)/L(ii,3))*(Ber(x31/V_T)));

            Jaco_tmp_DD(3,1)=coeff_Jp*((edge(ii,1)/L(ii,1))*(p2*(V_T)*Ber_d(-x21/V_T)-p1*(-V_T)*Ber_d(x21/V_T))+((edge(ii,3)/L(ii,3))*(p3*(V_T)*Ber_d(-x31/V_T)-p1*(-V_T)*Ber_d(x31/V_T))));
            Jaco_tmp_DD(3,3)=coeff_Jp*((edge(ii,1)/L(ii,1))*(-Ber(x21/V_T))+(edge(ii,3)/L(ii,3))*(-Ber(x31/V_T)));
            Jaco_tmp_DD(3,4)=coeff_Jp*((edge(ii,1)/L(ii,1))*(p2*(V_T)*Ber_d(-x21/V_T)-p1*(V_T)*Ber_d(x21/V_T)));
            Jaco_tmp_DD(3,6)=coeff_Jp*((edge(ii,1)/L(ii,1))*(Ber(-x21/V_T)));
            Jaco_tmp_DD(3,7)=coeff_Jp*((edge(ii,3)/L(ii,3))*(p3*(-V_T)*Ber_d(-x31/V_T)-p1*(V_T)*Ber_d(x31/V_T)));
            Jaco_tmp_DD(3,9)=coeff_Jp*((edge(ii,3)/L(ii,3))*(Ber(-x31/V_T)));

            Jaco_tmp_DD(4,1)=eps(Table_element_region(ii,1),1)*edge(ii,1)/L(ii,1);
            Jaco_tmp_DD(4,4)=eps(Table_element_region(ii,1),1)*(-edge(ii,1)/L(ii,1)-edge(ii,2)/L(ii,2));
            Jaco_tmp_DD(4,5)=-coeff*Control_Volume(2,1);
            Jaco_tmp_DD(4,6)=coeff*Control_Volume(2,1);
            Jaco_tmp_DD(4,7)=eps(Table_element_region(ii,1),1)*edge(ii,2)/L(ii,2);

            Jaco_tmp_DD(5,1)=coeff_Jn*((edge(ii,1)/L(ii,1))*(n1*(V_T)*Ber_d(x12/V_T)-n2*(-V_T)*Ber_d(-x12/V_T)));
            Jaco_tmp_DD(5,2)=coeff_Jn*((edge(ii,1)/L(ii,1))*(Ber(x12/V_T)));
            Jaco_tmp_DD(5,4)=coeff_Jn*((edge(ii,2)/L(ii,2))*(n3*(-V_T)*Ber_d(x32/V_T)-n2*(V_T)*Ber_d(-x32/V_T))+((edge(ii,1)/L(ii,1))*(n1*(-V_T)*Ber_d(x12/V_T)-n2*(V_T)*Ber_d(-x12/V_T))));
            Jaco_tmp_DD(5,5)=coeff_Jn*((edge(ii,2)/L(ii,2))*(-Ber(-x32/V_T))+((edge(ii,1)/L(ii,1))*(-Ber(-x12/V_T))));
            Jaco_tmp_DD(5,7)=coeff_Jn*((edge(ii,2)/L(ii,2))*(n3*(V_T)*Ber_d(x32/V_T)-n2*(-V_T)*Ber_d(-x32/V_T)));
            Jaco_tmp_DD(5,8)=coeff_Jn*((edge(ii,2)/L(ii,2))*(Ber(x32/V_T)));

            Jaco_tmp_DD(6,1)=coeff_Jp*((edge(ii,1)/L(ii,1))*(p1*(-V_T)*Ber_d(-x12/V_T)-p2*(V_T)*Ber_d(x12/V_T)));
            Jaco_tmp_DD(6,3)=coeff_Jp*((edge(ii,1)/L(ii,1))*(Ber(-x12/V_T)));
            Jaco_tmp_DD(6,4)=coeff_Jp*((edge(ii,2)/L(ii,2))*(p3*(V_T)*Ber_d(-x32/V_T)-p2*(-V_T)*Ber_d(x32/V_T))+((edge(ii,1)/L(ii,1))*(p1*(V_T)*Ber_d(-x12/V_T)-p2*(-V_T)*Ber_d(x12/V_T))));
            Jaco_tmp_DD(6,6)=coeff_Jp*((edge(ii,2)/L(ii,2))*(-Ber(x32/V_T))+(edge(ii,1)/L(ii,1))*(-Ber(x12/V_T)));
            Jaco_tmp_DD(6,7)=coeff_Jp*((edge(ii,2)/L(ii,2))*(p3*(V_T)*Ber_d(-x32/V_T)-p2*(V_T)*Ber_d(x32/V_T)));
            Jaco_tmp_DD(6,9)=coeff_Jp*((edge(ii,2)/L(ii,2))*(Ber(-x32/V_T)));

            Jaco_tmp_DD(7,1)=eps(Table_element_region(ii,1),1)*edge(ii,3)/L(ii,3);
            Jaco_tmp_DD(7,4)=eps(Table_element_region(ii,1),1)*edge(ii,2)/L(ii,2);
            Jaco_tmp_DD(7,7)=eps(Table_element_region(ii,1),1)*(-edge(ii,2)/L(ii,2)-edge(ii,3)/L(ii,3));
            Jaco_tmp_DD(7,8)=-coeff*Control_Volume(3,1);
            Jaco_tmp_DD(7,9)=coeff*Control_Volume(3,1);

            Jaco_tmp_DD(8,1)=coeff_Jn*((edge(ii,3)/L(ii,3))*(n1*(V_T)*Ber_d(x13/V_T)-n3*(-V_T)*Ber_d(-x13/V_T)));
            Jaco_tmp_DD(8,2)=coeff_Jn*((edge(ii,3)/L(ii,3))*(Ber(x13/V_T)));
            Jaco_tmp_DD(8,4)=coeff_Jn*((edge(ii,2)/L(ii,2))*(n2*(V_T)*Ber_d(x23/V_T)-n3*(-V_T)*Ber_d(-x23/V_T)));
            Jaco_tmp_DD(8,5)=coeff_Jn*((edge(ii,2)/L(ii,2))*(Ber(x23/V_T)));
            Jaco_tmp_DD(8,7)=coeff_Jn*((edge(ii,3)/L(ii,3))*(n1*(-V_T)*Ber_d(x13/V_T)-n3*(V_T)*Ber_d(-x13/V_T))+((edge(ii,2)/L(ii,2))*(n2*(-V_T)*Ber_d(x23/V_T)-n3*(V_T)*Ber_d(-x23/V_T))));
            Jaco_tmp_DD(8,8)=coeff_Jn*((edge(ii,3)/L(ii,3))*(-Ber(-x13/V_T))+((edge(ii,2)/L(ii,2))*(-Ber(-x23/V_T))));

            Jaco_tmp_DD(9,1)=coeff_Jp*((edge(ii,3)/L(ii,3))*(p1*(V_T)*Ber_d(-x13/V_T)-p3*(V_T)*Ber_d(x13/V_T)));
            Jaco_tmp_DD(9,3)=coeff_Jp*((edge(ii,3)/L(ii,3))*(Ber(-x13/V_T)));
            Jaco_tmp_DD(9,4)=coeff_Jp*((edge(ii,2)/L(ii,2))*(p2*(-V_T)*Ber_d(-x23/V_T)-p3*(V_T)*Ber_d(x23/V_T)));
            Jaco_tmp_DD(9,6)=coeff_Jp*((edge(ii,2)/L(ii,2))*(Ber(-x23/V_T)));
            Jaco_tmp_DD(9,7)=coeff_Jp*((edge(ii,3)/L(ii,3))*(p1*(V_T)*Ber_d(-x13/V_T)-p3*(-V_T)*Ber_d(x13/V_T))+((edge(ii,2)/L(ii,2))*(p2*(V_T)*Ber_d(-x23/V_T)-p3*(-V_T)*Ber_d(x23/V_T))));
            Jaco_tmp_DD(9,9)=coeff_Jp*((edge(ii,3)/L(ii,3))*(-Ber(x13/V_T))+(edge(ii,2)/L(ii,2))*(-Ber(x23/V_T)));


            for j=1:9
                for k=1:9
                    Jaco_DD(index(j,1),index(k,1))=Jaco_DD(index(j,1),index(k,1))+Jaco_tmp_DD(j,k);
                end
            end

            % res vector %
            % potential
            res_tmp_DD=zeros(9,1);
            res_potential_tmp1_DD=eps(Table_element_region(ii,1),1).*[edge(ii,1)/L(ii,1)*(solution_vector_DD(index(4,1),Newton_DD)-solution_vector_DD(index(1,1),Newton_DD))+edge(ii,3)/L(ii,3)*(solution_vector_DD(index(7,1),Newton_DD)-solution_vector_DD(index(1,1),Newton_DD));
                edge(ii,2)/L(ii,2)*(solution_vector_DD(index(7,1),Newton_DD)-solution_vector_DD(index(4,1),Newton_DD))+edge(ii,1)/L(ii,1)*(solution_vector_DD(index(1,1),Newton_DD)-solution_vector_DD(index(4,1),Newton_DD));
                edge(ii,3)/L(ii,3)*(solution_vector_DD(index(1,1),Newton_DD)-solution_vector_DD(index(7,1),Newton_DD))+edge(ii,2)/L(ii,2)*(solution_vector_DD(index(4,1),Newton_DD)-solution_vector_DD(index(7,1),Newton_DD))];
            for j=1:3
                res_potential_tmp2_DD(j,1)=(N_dop(Table_element_region(ii,1),1)-solution_vector_DD(index(3*j-1,1),Newton_DD)+solution_vector_DD(index(3*j,1),Newton_DD))*coeff*Control_Volume(j,1);
            end
            res_potential=res_potential_tmp1_DD+res_potential_tmp2_DD;

            % Jn
            res_Jn=[coeff_Jn*((edge(ii,1)/L(ii,1))*(n2*Ber(x21/V_T)-n1*Ber(-x21/V_T))+(edge(ii,3)/L(ii,3))*(n3*Ber(x31/V_T)-n1*Ber(-x31/V_T)));
                coeff_Jn*((edge(ii,2)/L(ii,2))*(n3*Ber(x32/V_T)-n2*Ber(-x32/V_T))+(edge(ii,1)/L(ii,1))*(n1*Ber(x12/V_T)-n2*Ber(-x12/V_T)));
                coeff_Jn*((edge(ii,3)/L(ii,3))*(n1*Ber(x13/V_T)-n3*Ber(-x13/V_T))+(edge(ii,2)/L(ii,2))*(n2*Ber(x23/V_T)-n3*Ber(-x23/V_T)))];


            % Jp
            res_Jp=[coeff_Jp*((edge(ii,1)/L(ii,1))*(p2*Ber(-x21/V_T)-p1*Ber(x21/V_T))+(edge(ii,3)/L(ii,3))*(p3*Ber(-x31/V_T)-p1*Ber(x31/V_T)));
                coeff_Jp*((edge(ii,2)/L(ii,2))*(p3*Ber(-x32/V_T)-p2*Ber(x32/V_T))+(edge(ii,1)/L(ii,1))*(p1*Ber(-x12/V_T)-p2*Ber(x12/V_T)));
                coeff_Jp*((edge(ii,3)/L(ii,3))*(p1*Ber(-x13/V_T)-p3*Ber(x13/V_T))+(edge(ii,2)/L(ii,2))*(p2*Ber(-x23/V_T)-p3*Ber(x23/V_T)))];

            n=1;
            for j=1:3
                res_tmp_DD(n,1)=res_potential(j,1);
                res_tmp_DD(n+1,1)=res_Jn(j,1);
                res_tmp_DD(n+2,1)=res_Jp(j,1);
                n=n+3;
            end

            for j=1:9
                res_DD(index(j,1),1)=res_DD(index(j,1),1)+res_tmp_DD(j,1);
            end

        else
            % index : 현재 element에 속해있는 vertex가 resdue에 몇번째 있는가 저장
            index=zeros(3,1);
            n=1;
            for j=1:3
                index(j,1)=find(Table_Jaco(:,1)==Table_element_region(ii,1) & Table_Jaco(:,2)==Element(ii,j) & Table_Jaco(:,3)==1);
            end

            % Jaco_potential matrix
            Jaco_tmp_DD=zeros(3,3);
            Jaco_tmp_DD=[(-edge(ii,1)/L(ii,1)-edge(ii,3)/L(ii,3)) edge(ii,1)/L(ii,1) edge(ii,3)/L(ii,3);
                edge(ii,1)/L(ii,1) -edge(ii,1)/L(ii,1)-edge(ii,2)/L(ii,2) edge(ii,2)/L(ii,2);
                edge(ii,3)/L(ii,3) edge(ii,2)/L(ii,2) -edge(ii,2)/L(ii,2)-edge(ii,3)/L(ii,3)];
            Jaco_tmp_DD=eps(Table_element_region(ii,1),1)*Jaco_tmp_DD;
            for j=1:3
                for k=1:3
                    Jaco_DD(index(j,1),index(k,1))=Jaco_DD(index(j,1),index(k,1))+Jaco_tmp_DD(j,k);
                end
            end

            % res vector
            res_tmp_DD=zeros(3,1);
            res_tmp_DD=eps(Table_element_region(ii,1),1).*[(-edge(ii,1)/L(ii,1)-edge(ii,3)/L(ii,3))*solution_vector_DD(index(1,1),Newton_DD)+edge(ii,1)/L(ii,1)*solution_vector_DD(index(2,1),Newton_DD)+edge(ii,3)/L(ii,3)*solution_vector_DD(index(3,1),Newton_DD);
                edge(ii,1)/L(ii,1)*solution_vector_DD(index(1,1),Newton_DD)+(-edge(ii,1)/L(ii,1)-edge(ii,2)/L(ii,2))*solution_vector_DD(index(2,1),Newton_DD)+edge(ii,2)/L(ii,2)*solution_vector_DD(index(3,1),Newton_DD);
                edge(ii,3)/L(ii,3)*solution_vector_DD(index(1,1),Newton_DD)+edge(ii,2)/L(ii,2)*solution_vector_DD(index(2,1),Newton_DD)+(-edge(ii,2)/L(ii,2)-edge(ii,3)/L(ii,3))*solution_vector_DD(index(3,1),Newton_DD)];
            for j=1:3
                res_DD(index(j,1),1)=res_DD(index(j,1),1)+res_tmp_DD(j,1);
            end
        end
    end

    % ox 의 값을 si에 더해준 후 ox값은 1,-1로 변경하기
    for ii=[1 5]
        for j=[3 2 4]
            clearvars Vertex_interface_tmp
            Vertex_interface_tmp=Vertex_interface{ii,j};
            for k=1:size(Vertex_interface_tmp,1)
                index_ox=find(Table_Jaco(:,1)==ii & Table_Jaco(:,2)==Vertex_interface_tmp(k,1) & Table_Jaco(:,3)==1);
                index_si=find(Table_Jaco(:,1)==j & Table_Jaco(:,2)==Vertex_interface_tmp(k,1) & Table_Jaco(:,3)==1);
                if length(sort(nonzeros(Jaco_DD(index_ox,:))))>2
                    Jaco_DD(index_si,:)=Jaco_DD(index_si,:)+Jaco_DD(index_ox,:);
                    Jaco_DD(index_ox,:)=0; Jaco_DD(index_ox,index_ox)=1; Jaco_DD(index_ox,index_si)=-1;
                    res_DD(index_si,1)=res_DD(index_si,1)+res_DD(index_ox,1);
                    res_DD(index_ox,1)=0;
                elseif ~isequal(sort(nonzeros(Jaco_DD(index_ox,:))),[-1; 1]) % 점 3개가 접해있는 지역에서 중복적으로 합쳐지는것을 막기위해 추가함. jaco의 ox-channel행이 -1과1이라면 이미 다른행과 합쳐진 것이므로 pass.
                    Jaco_DD(index_si,:)=Jaco_DD(index_si,:)+Jaco_DD(index_ox,:);
                    Jaco_DD(index_ox,:)=0; Jaco_DD(index_ox,index_ox)=1; Jaco_DD(index_ox,index_si)=-1;
                    res_DD(index_si,1)=res_DD(index_si,1)+res_DD(index_ox,1);
                    res_DD(index_ox,1)=0;
                end
            end
        end
    end
    % source/drain 값을 channel에 더해주기. Si부분이니까 electron/hole density도 모두 바꾸어주어야함.
    for ii=[2 4]
        for j=3
            clearvars Vertex_interface_tmp
            Vertex_interface_tmp=Vertex_interface{ii,j};
            for k=1:size(Vertex_interface_tmp,1)
                % index 설정
                index_source_drain_potential=find(Table_Jaco(:,1)==ii & Table_Jaco(:,2)==Vertex_interface_tmp(k,1) & Table_Jaco(:,3)==1);
                index_source_drain_electron=index_source_drain_potential+1;
                index_source_drain_hole=index_source_drain_potential+2;
                index_channel_potential=find(Table_Jaco(:,1)==j & Table_Jaco(:,2)==Vertex_interface_tmp(k,1) & Table_Jaco(:,3)==1);
                index_channel_electron=index_channel_potential+1;
                index_channel_hole=index_channel_potential+2;
                % Jaco 변경, potential/electron/hole 순서
                Jaco_DD(index_channel_potential,:)=Jaco_DD(index_channel_potential,:)+Jaco_DD(index_source_drain_potential,:);
                Jaco_DD(index_channel_electron,:)=Jaco_DD(index_channel_electron,:)+Jaco_DD(index_source_drain_electron,:);
                Jaco_DD(index_channel_hole,:)=Jaco_DD(index_channel_hole,:)+Jaco_DD(index_source_drain_hole,:);
                Jaco_DD(index_source_drain_potential,:)=0; Jaco_DD(index_source_drain_potential,index_source_drain_potential)=1; Jaco_DD(index_source_drain_potential,index_channel_potential)=-1;
                Jaco_DD(index_source_drain_electron,:)=0; Jaco_DD(index_source_drain_electron,index_source_drain_electron)=1; Jaco_DD(index_source_drain_electron,index_channel_electron)=-1;
                Jaco_DD(index_source_drain_hole,:)=0; Jaco_DD(index_source_drain_hole,index_source_drain_hole)=1; Jaco_DD(index_source_drain_hole,index_channel_hole)=-1;

                % res 변경, potential/electron/hole 순서
                res_DD(index_channel_potential,1)=res_DD(index_channel_potential,1)+res_DD(index_source_drain_potential,1);
                res_DD(index_channel_electron,1)=res_DD(index_channel_electron,1)+res_DD(index_source_drain_electron,1);
                res_DD(index_channel_hole,1)=res_DD(index_channel_hole,1)+res_DD(index_source_drain_hole,1);
                res_DD(index_source_drain_potential,1)=0;
                res_DD(index_source_drain_electron,1)=0;
                res_DD(index_source_drain_hole,1)=0;
            end
        end
    end

    % Dirichlet_BC
    for ii=1:size(Dirichlet_BC,1)
        BC_index=find(Table_Jaco(:,2)==Dirichlet_BC(ii,1) & Table_Jaco(:,3)==1);
        Jaco_DD(BC_index,:)=0;
        Jaco_DD(BC_index,BC_index)=1;
        res_DD(BC_index,1)=solution_vector_DD(BC_index,Newton_DD)-V_gate;
    end

    % Source_drain_BC
    for ii=1:size(Source_drain_BC,1)
        BC_index=find((Table_Jaco(:,1)==2 | Table_Jaco(:,1)==4) & Table_Jaco(:,2)==Source_drain_BC(ii,1) & Table_Jaco(:,3)==1);
        Jaco_DD(BC_index,:)=0; Jaco_DD(BC_index+1,:)=0; Jaco_DD(BC_index+2,:)=0;
        Jaco_DD(BC_index,BC_index)=1; Jaco_DD(BC_index+1,BC_index+1)=1; Jaco_DD(BC_index+2,BC_index+2)=1;
        res_DD(BC_index,1)=solution_vector_DD(BC_index,Newton_DD)-V_T*log(Nd/n_int);
        res_DD(BC_index+1,1)=solution_vector_DD(BC_index+1,Newton_DD)-Nd;
        res_DD(BC_index+2,1)=solution_vector_DD(BC_index+2,Newton_DD)-n_int^2/Nd;
    end

    res_save_DD(:,Newton_DD)=res_DD(:,1);

    % Scaling
    Cvector=zeros(size(solution_vector_DD,1),1); v=[V_T; max(abs(N_dop)); max(abs(N_dop))];
    for ii=1:size(solution_vector_DD,1)
        Cvector(ii,1)=v(Table_Jaco(ii,3),1);
    end
    Cmatrix_DD=spdiags(Cvector,0,size(solution_vector_DD,1),size(solution_vector_DD,1));
    Jaco_scaled_DD=Jaco_DD*Cmatrix_DD;
    Rvector_DD=1./sum(abs(Jaco_scaled_DD),2);
    Rmatrix_DD=spdiags(Rvector_DD,0,size(solution_vector_DD,1),size(solution_vector_DD,1));
    Jaco_scaled_DD=Rmatrix_DD*Jaco_scaled_DD;
    res_scaled_DD=Rmatrix_DD*res_DD;
    update_scaled_DD=Jaco_scaled_DD\(-res_scaled_DD);
    update_DD(:,Newton_DD)=Cmatrix_DD*update_scaled_DD;
    solution_vector_DD(:,Newton_DD+1)=solution_vector_DD(:,Newton_DD)+update_DD(:,Newton_DD);

    %     % non-Scaling
    %     update(:,Newton)=Jaco\-res;
    %     solution_vector(:,Newton+1)=solution_vector(:,Newton)+update(:,Newton);

end

% % % Visualize
% Visual_solution_vector_DD=zeros(size(Vertex,1),3);
% Visual_initial_vector_DD=zeros(size(Vertex,1),3);
% for ii=1:size(Table_Jaco,1)
%     Visual_solution_vector_DD(Table_Jaco(ii,2),Table_Jaco(ii,3))=solution_vector_DD(ii,Newton_DD);
%     Visual_initial_vector_DD(Table_Jaco(ii,2),Table_Jaco(ii,3))=solution_vector_DD(ii,1);
% end
%
% dif_DD=abs(Visual_initial_vector_DD-Visual_solution_vector_DD);
% plot3(Vertex(:,1),Vertex(:,2),Visual_solution_vector_DD(:,1),'*')
% xlim([0 30*1e-9])
%
% figure(4);
% patch('Faces',Element,'Vertices',Vertex(:,[1 2]), 'FaceVertexCData',Visual_solution_vector_DD(:,1), 'EdgeColor','black','FaceColor','interp');
% title('Initial potential')
% xlim([0 30*1e-9])
% xlabel('X');
% ylabel('Y');
% colorbar
%
% figure(5);
% patch('Faces',Element,'Vertices',Vertex(:,[1 2]), 'FaceVertexCData',Visual_solution_vector_DD(:,2), 'EdgeColor','black','FaceColor','interp');
% title('elctron')
% xlim([0 30*1e-9])
% xlabel('X');
% ylabel('Y');
% colorbar
%
% figure(6);
% patch('Faces',Element,'Vertices',Vertex(:,[1 2]), 'FaceVertexCData',Visual_solution_vector_DD(:,3), 'EdgeColor','black','FaceColor','interp');
% title('hole')
% xlim([0 30*1e-9])
% xlabel('X');
% ylabel('Y');
% colorbar
%
% % figure(7);
% % patch('Faces',Element,'Vertices',Vertex(:,[1 2]), 'FaceVertexCData',dif_DD(:,1), 'EdgeColor','black','FaceColor','interp');
% % title('Initial potential')
% % xlabel('X');
% % ylabel('Y');
% % colorbar
% %
% % figure(8);
% % patch('Faces',Element,'Vertices',Vertex(:,[1 2]), 'FaceVertexCData',dif_DD(:,2), 'EdgeColor','black','FaceColor','interp');
% % title('elctron')
% % xlabel('X');
% % ylabel('Y');
% % colorbar
% %
% % figure(9);
% % patch('Faces',Element,'Vertices',Vertex(:,[1 2]), 'FaceVertexCData',dif_DD(:,3), 'EdgeColor','black','FaceColor','interp');
% % title('hole')
% % xlabel('X');
% % ylabel('Y');
% % colorbar

%% Ramping Drain
solution_vector_DD(:,1)=solution_vector_possion(:,Newton_possion);

for bias=1:11
    Vd=bias*0.1-0.1;

    % Jacobian matrix / res vector
    iteration=10;
    for Newton_DD=1:iteration
        %     Jaco_DD=sparse(zeros(size(solution_vector_DD,1),size(solution_vector_DD,1))); res_DD=zeros(size(solution_vector_DD,1),1);
        Jaco_DD=zeros(size(solution_vector_DD,1),size(solution_vector_DD,1)); res_DD=zeros(size(solution_vector_DD,1),1); %% 희소행렬 아님. Jaco matrix 확인용.
        % Potential
        for ii=1:length(Element)
            % Control_Volume 계산
            Control_Volume=zeros(3,1);
            Control_Volume=[(0.25)*L(ii,3)*edge(ii,3)+(0.25)*L(ii,1)*edge(ii,1);
                (0.25)*L(ii,1)*edge(ii,1)+(0.25)*L(ii,2)*edge(ii,2);
                (0.25)*L(ii,2)*edge(ii,2)+(0.25)*L(ii,3)*edge(ii,3)];
            coeff=q/eps0;     % m^2 -> nm^2으로 단위 변경, q/eps0 해줌
            coeff_Jn=q*mobility_n*V_T;
            coeff_Jp=-q*mobility_p*V_T;

            if Table_element_region(ii,1)==2 || Table_element_region(ii,1)==3 || Table_element_region(ii,1)==4 %%% Region 변경되면 꼭 수정되어야함!!
                % index : 현재 element에 속해있는 vertex가 resdue에 몇번째 있는가 저장
                index=zeros(9,1);
                n=1;
                for j=1:3
                    index(n,1)=find(Table_Jaco(:,1)==Table_element_region(ii,1) & Table_Jaco(:,2)==Element(ii,j) & Table_Jaco(:,3)==1);
                    index(n+1,1)=index(n,1)+1;
                    index(n+2,1)=index(n,1)+2;
                    n=n+3;
                end

                % 식 간단히 표기하기 위해 미리 계산.
                x12=solution_vector_DD(index(1,1),1)-solution_vector_DD(index(4,1),1); x21=-x12;
                x23=solution_vector_DD(index(4,1),1)-solution_vector_DD(index(7,1),1); x32=-x23;
                x13=solution_vector_DD(index(1,1),1)-solution_vector_DD(index(7,1),1); x31=-x13;
                n1=solution_vector_DD(index(1,1)+1); n2=solution_vector_DD(index(4,1)+1); n3=solution_vector_DD(index(7,1)+1);
                p1=solution_vector_DD(index(1,1)+2); p2=solution_vector_DD(index(4,1)+2); p3=solution_vector_DD(index(7,1)+2);

                Jaco_tmp_DD=zeros(9,9); % 9*9 Jaco matrix 만듬.
                Jaco_tmp_DD(1,1)=eps(Table_element_region(ii,1),1)*(-edge(ii,1)/L(ii,1)-edge(ii,3)/L(ii,3));
                Jaco_tmp_DD(1,2)=-coeff*Control_Volume(1,1);
                Jaco_tmp_DD(1,3)=coeff*Control_Volume(1,1);
                Jaco_tmp_DD(1,4)=eps(Table_element_region(ii,1),1)*edge(ii,1)/L(ii,1);
                Jaco_tmp_DD(1,7)=eps(Table_element_region(ii,1),1)*edge(ii,3)/L(ii,3);

                Jaco_tmp_DD(2,1)=coeff_Jn*((edge(ii,1)/L(ii,1))*(n2*(-V_T)*Ber_d(x21/V_T)-n1*(V_T)*Ber_d(-x21/V_T))+((edge(ii,3)/L(ii,3))*(n3*(-V_T)*Ber_d(x31/V_T)-n1*(V_T)*Ber_d(-x31/V_T))));
                Jaco_tmp_DD(2,2)=coeff_Jn*((edge(ii,1)/L(ii,1))*(-Ber(-x21/V_T))+((edge(ii,3)/L(ii,3))*(-Ber(-x31/V_T))));
                Jaco_tmp_DD(2,4)=coeff_Jn*((edge(ii,1)/L(ii,1))*(n2*(V_T)*Ber_d(x21/V_T)-n1*(-V_T)*Ber_d(-x21/V_T)));
                Jaco_tmp_DD(2,5)=coeff_Jn*((edge(ii,1)/L(ii,1))*(Ber(x21/V_T)));
                Jaco_tmp_DD(2,7)=coeff_Jn*((edge(ii,3)/L(ii,3))*(n3*(V_T)*Ber_d(x31/V_T)-n1*(-V_T)*Ber_d(-x31/V_T)));
                Jaco_tmp_DD(2,8)=coeff_Jn*((edge(ii,3)/L(ii,3))*(Ber(x31/V_T)));

                Jaco_tmp_DD(3,1)=coeff_Jp*((edge(ii,1)/L(ii,1))*(p2*(V_T)*Ber_d(-x21/V_T)-p1*(-V_T)*Ber_d(x21/V_T))+((edge(ii,3)/L(ii,3))*(p3*(V_T)*Ber_d(-x31/V_T)-p1*(-V_T)*Ber_d(x31/V_T))));
                Jaco_tmp_DD(3,3)=coeff_Jp*((edge(ii,1)/L(ii,1))*(-Ber(x21/V_T))+(edge(ii,3)/L(ii,3))*(-Ber(x31/V_T)));
                Jaco_tmp_DD(3,4)=coeff_Jp*((edge(ii,1)/L(ii,1))*(p2*(V_T)*Ber_d(-x21/V_T)-p1*(V_T)*Ber_d(x21/V_T)));
                Jaco_tmp_DD(3,6)=coeff_Jp*((edge(ii,1)/L(ii,1))*(Ber(-x21/V_T)));
                Jaco_tmp_DD(3,7)=coeff_Jp*((edge(ii,3)/L(ii,3))*(p3*(-V_T)*Ber_d(-x31/V_T)-p1*(V_T)*Ber_d(x31/V_T)));
                Jaco_tmp_DD(3,9)=coeff_Jp*((edge(ii,3)/L(ii,3))*(Ber(-x31/V_T)));

                Jaco_tmp_DD(4,1)=eps(Table_element_region(ii,1),1)*edge(ii,1)/L(ii,1);
                Jaco_tmp_DD(4,4)=eps(Table_element_region(ii,1),1)*(-edge(ii,1)/L(ii,1)-edge(ii,2)/L(ii,2));
                Jaco_tmp_DD(4,5)=-coeff*Control_Volume(2,1);
                Jaco_tmp_DD(4,6)=coeff*Control_Volume(2,1);
                Jaco_tmp_DD(4,7)=eps(Table_element_region(ii,1),1)*edge(ii,2)/L(ii,2);

                Jaco_tmp_DD(5,1)=coeff_Jn*((edge(ii,1)/L(ii,1))*(n1*(V_T)*Ber_d(x12/V_T)-n2*(-V_T)*Ber_d(-x12/V_T)));
                Jaco_tmp_DD(5,2)=coeff_Jn*((edge(ii,1)/L(ii,1))*(Ber(x12/V_T)));
                Jaco_tmp_DD(5,4)=coeff_Jn*((edge(ii,2)/L(ii,2))*(n3*(-V_T)*Ber_d(x32/V_T)-n2*(V_T)*Ber_d(-x32/V_T))+((edge(ii,1)/L(ii,1))*(n1*(-V_T)*Ber_d(x12/V_T)-n2*(V_T)*Ber_d(-x12/V_T))));
                Jaco_tmp_DD(5,5)=coeff_Jn*((edge(ii,2)/L(ii,2))*(-Ber(-x32/V_T))+((edge(ii,1)/L(ii,1))*(-Ber(-x12/V_T))));
                Jaco_tmp_DD(5,7)=coeff_Jn*((edge(ii,2)/L(ii,2))*(n3*(V_T)*Ber_d(x32/V_T)-n2*(-V_T)*Ber_d(-x32/V_T)));
                Jaco_tmp_DD(5,8)=coeff_Jn*((edge(ii,2)/L(ii,2))*(Ber(x32/V_T)));

                Jaco_tmp_DD(6,1)=coeff_Jp*((edge(ii,1)/L(ii,1))*(p1*(-V_T)*Ber_d(-x12/V_T)-p2*(V_T)*Ber_d(x12/V_T)));
                Jaco_tmp_DD(6,3)=coeff_Jp*((edge(ii,1)/L(ii,1))*(Ber(-x12/V_T)));
                Jaco_tmp_DD(6,4)=coeff_Jp*((edge(ii,2)/L(ii,2))*(p3*(V_T)*Ber_d(-x32/V_T)-p2*(-V_T)*Ber_d(x32/V_T))+((edge(ii,1)/L(ii,1))*(p1*(V_T)*Ber_d(-x12/V_T)-p2*(-V_T)*Ber_d(x12/V_T))));
                Jaco_tmp_DD(6,6)=coeff_Jp*((edge(ii,2)/L(ii,2))*(-Ber(x32/V_T))+(edge(ii,1)/L(ii,1))*(-Ber(x12/V_T)));
                Jaco_tmp_DD(6,7)=coeff_Jp*((edge(ii,2)/L(ii,2))*(p3*(V_T)*Ber_d(-x32/V_T)-p2*(V_T)*Ber_d(x32/V_T)));
                Jaco_tmp_DD(6,9)=coeff_Jp*((edge(ii,2)/L(ii,2))*(Ber(-x32/V_T)));

                Jaco_tmp_DD(7,1)=eps(Table_element_region(ii,1),1)*edge(ii,3)/L(ii,3);
                Jaco_tmp_DD(7,4)=eps(Table_element_region(ii,1),1)*edge(ii,2)/L(ii,2);
                Jaco_tmp_DD(7,7)=eps(Table_element_region(ii,1),1)*(-edge(ii,2)/L(ii,2)-edge(ii,3)/L(ii,3));
                Jaco_tmp_DD(7,8)=-coeff*Control_Volume(3,1);
                Jaco_tmp_DD(7,9)=coeff*Control_Volume(3,1);

                Jaco_tmp_DD(8,1)=coeff_Jn*((edge(ii,3)/L(ii,3))*(n1*(V_T)*Ber_d(x13/V_T)-n3*(-V_T)*Ber_d(-x13/V_T)));
                Jaco_tmp_DD(8,2)=coeff_Jn*((edge(ii,3)/L(ii,3))*(Ber(x13/V_T)));
                Jaco_tmp_DD(8,4)=coeff_Jn*((edge(ii,2)/L(ii,2))*(n2*(V_T)*Ber_d(x23/V_T)-n3*(-V_T)*Ber_d(-x23/V_T)));
                Jaco_tmp_DD(8,5)=coeff_Jn*((edge(ii,2)/L(ii,2))*(Ber(x23/V_T)));
                Jaco_tmp_DD(8,7)=coeff_Jn*((edge(ii,3)/L(ii,3))*(n1*(-V_T)*Ber_d(x13/V_T)-n3*(V_T)*Ber_d(-x13/V_T))+((edge(ii,2)/L(ii,2))*(n2*(-V_T)*Ber_d(x23/V_T)-n3*(V_T)*Ber_d(-x23/V_T))));
                Jaco_tmp_DD(8,8)=coeff_Jn*((edge(ii,3)/L(ii,3))*(-Ber(-x13/V_T))+((edge(ii,2)/L(ii,2))*(-Ber(-x23/V_T))));

                Jaco_tmp_DD(9,1)=coeff_Jp*((edge(ii,3)/L(ii,3))*(p1*(V_T)*Ber_d(-x13/V_T)-p3*(V_T)*Ber_d(x13/V_T)));
                Jaco_tmp_DD(9,3)=coeff_Jp*((edge(ii,3)/L(ii,3))*(Ber(-x13/V_T)));
                Jaco_tmp_DD(9,4)=coeff_Jp*((edge(ii,2)/L(ii,2))*(p2*(-V_T)*Ber_d(-x23/V_T)-p3*(V_T)*Ber_d(x23/V_T)));
                Jaco_tmp_DD(9,6)=coeff_Jp*((edge(ii,2)/L(ii,2))*(Ber(-x23/V_T)));
                Jaco_tmp_DD(9,7)=coeff_Jp*((edge(ii,3)/L(ii,3))*(p1*(V_T)*Ber_d(-x13/V_T)-p3*(-V_T)*Ber_d(x13/V_T))+((edge(ii,2)/L(ii,2))*(p2*(V_T)*Ber_d(-x23/V_T)-p3*(-V_T)*Ber_d(x23/V_T))));
                Jaco_tmp_DD(9,9)=coeff_Jp*((edge(ii,3)/L(ii,3))*(-Ber(x13/V_T))+(edge(ii,2)/L(ii,2))*(-Ber(x23/V_T)));


                for j=1:9
                    for k=1:9
                        Jaco_DD(index(j,1),index(k,1))=Jaco_DD(index(j,1),index(k,1))+Jaco_tmp_DD(j,k);
                    end
                end

                % res vector %
                % potential
                res_tmp_DD=zeros(9,1);
                res_potential_tmp1_DD=eps(Table_element_region(ii,1),1).*[edge(ii,1)/L(ii,1)*(solution_vector_DD(index(4,1),Newton_DD)-solution_vector_DD(index(1,1),Newton_DD))+edge(ii,3)/L(ii,3)*(solution_vector_DD(index(7,1),Newton_DD)-solution_vector_DD(index(1,1),Newton_DD));
                    edge(ii,2)/L(ii,2)*(solution_vector_DD(index(7,1),Newton_DD)-solution_vector_DD(index(4,1),Newton_DD))+edge(ii,1)/L(ii,1)*(solution_vector_DD(index(1,1),Newton_DD)-solution_vector_DD(index(4,1),Newton_DD));
                    edge(ii,3)/L(ii,3)*(solution_vector_DD(index(1,1),Newton_DD)-solution_vector_DD(index(7,1),Newton_DD))+edge(ii,2)/L(ii,2)*(solution_vector_DD(index(4,1),Newton_DD)-solution_vector_DD(index(7,1),Newton_DD))];
                for j=1:3
                    res_potential_tmp2_DD(j,1)=(N_dop(Table_element_region(ii,1),1)-solution_vector_DD(index(3*j-1,1),Newton_DD)+solution_vector_DD(index(3*j,1),Newton_DD))*coeff*Control_Volume(j,1);
                end
                res_potential=res_potential_tmp1_DD+res_potential_tmp2_DD;

                % Jn
                res_Jn=[coeff_Jn*((edge(ii,1)/L(ii,1))*(n2*Ber(x21/V_T)-n1*Ber(-x21/V_T))+(edge(ii,3)/L(ii,3))*(n3*Ber(x31/V_T)-n1*Ber(-x31/V_T)));
                    coeff_Jn*((edge(ii,2)/L(ii,2))*(n3*Ber(x32/V_T)-n2*Ber(-x32/V_T))+(edge(ii,1)/L(ii,1))*(n1*Ber(x12/V_T)-n2*Ber(-x12/V_T)));
                    coeff_Jn*((edge(ii,3)/L(ii,3))*(n1*Ber(x13/V_T)-n3*Ber(-x13/V_T))+(edge(ii,2)/L(ii,2))*(n2*Ber(x23/V_T)-n3*Ber(-x23/V_T)))];


                % Jp
                res_Jp=[coeff_Jp*((edge(ii,1)/L(ii,1))*(p2*Ber(-x21/V_T)-p1*Ber(x21/V_T))+(edge(ii,3)/L(ii,3))*(p3*Ber(-x31/V_T)-p1*Ber(x31/V_T)));
                    coeff_Jp*((edge(ii,2)/L(ii,2))*(p3*Ber(-x32/V_T)-p2*Ber(x32/V_T))+(edge(ii,1)/L(ii,1))*(p1*Ber(-x12/V_T)-p2*Ber(x12/V_T)));
                    coeff_Jp*((edge(ii,3)/L(ii,3))*(p1*Ber(-x13/V_T)-p3*Ber(x13/V_T))+(edge(ii,2)/L(ii,2))*(p2*Ber(-x23/V_T)-p3*Ber(x23/V_T)))];

                n=1;
                for j=1:3
                    res_tmp_DD(n,1)=res_potential(j,1);
                    res_tmp_DD(n+1,1)=res_Jn(j,1);
                    res_tmp_DD(n+2,1)=res_Jp(j,1);
                    n=n+3;
                end

                for j=1:9
                    res_DD(index(j,1),1)=res_DD(index(j,1),1)+res_tmp_DD(j,1);
                end

            else
                % index : 현재 element에 속해있는 vertex가 resdue에 몇번째 있는가 저장
                index=zeros(3,1);
                n=1;
                for j=1:3
                    index(j,1)=find(Table_Jaco(:,1)==Table_element_region(ii,1) & Table_Jaco(:,2)==Element(ii,j) & Table_Jaco(:,3)==1);
                end

                % Jaco_potential matrix
                Jaco_tmp_DD=zeros(3,3);
                Jaco_tmp_DD=[(-edge(ii,1)/L(ii,1)-edge(ii,3)/L(ii,3)) edge(ii,1)/L(ii,1) edge(ii,3)/L(ii,3);
                    edge(ii,1)/L(ii,1) -edge(ii,1)/L(ii,1)-edge(ii,2)/L(ii,2) edge(ii,2)/L(ii,2);
                    edge(ii,3)/L(ii,3) edge(ii,2)/L(ii,2) -edge(ii,2)/L(ii,2)-edge(ii,3)/L(ii,3)];
                Jaco_tmp_DD=eps(Table_element_region(ii,1),1)*Jaco_tmp_DD;
                for j=1:3
                    for k=1:3
                        Jaco_DD(index(j,1),index(k,1))=Jaco_DD(index(j,1),index(k,1))+Jaco_tmp_DD(j,k);
                    end
                end

                % res vector
                res_tmp_DD=zeros(3,1);
                res_tmp_DD=eps(Table_element_region(ii,1),1).*[(-edge(ii,1)/L(ii,1)-edge(ii,3)/L(ii,3))*solution_vector_DD(index(1,1),Newton_DD)+edge(ii,1)/L(ii,1)*solution_vector_DD(index(2,1),Newton_DD)+edge(ii,3)/L(ii,3)*solution_vector_DD(index(3,1),Newton_DD);
                    edge(ii,1)/L(ii,1)*solution_vector_DD(index(1,1),Newton_DD)+(-edge(ii,1)/L(ii,1)-edge(ii,2)/L(ii,2))*solution_vector_DD(index(2,1),Newton_DD)+edge(ii,2)/L(ii,2)*solution_vector_DD(index(3,1),Newton_DD);
                    edge(ii,3)/L(ii,3)*solution_vector_DD(index(1,1),Newton_DD)+edge(ii,2)/L(ii,2)*solution_vector_DD(index(2,1),Newton_DD)+(-edge(ii,2)/L(ii,2)-edge(ii,3)/L(ii,3))*solution_vector_DD(index(3,1),Newton_DD)];
                for j=1:3
                    res_DD(index(j,1),1)=res_DD(index(j,1),1)+res_tmp_DD(j,1);
                end
            end
        end

        % ox 의 값을 si에 더해준 후 ox값은 1,-1로 변경하기
        for ii=[1 5]
            for j=[3 2 4]
                clearvars Vertex_interface_tmp
                Vertex_interface_tmp=Vertex_interface{ii,j};
                for k=1:size(Vertex_interface_tmp,1)
                    index_ox=find(Table_Jaco(:,1)==ii & Table_Jaco(:,2)==Vertex_interface_tmp(k,1) & Table_Jaco(:,3)==1);
                    index_si=find(Table_Jaco(:,1)==j & Table_Jaco(:,2)==Vertex_interface_tmp(k,1) & Table_Jaco(:,3)==1);
                    if length(sort(nonzeros(Jaco_DD(index_ox,:))))>2
                        Jaco_DD(index_si,:)=Jaco_DD(index_si,:)+Jaco_DD(index_ox,:);
                        Jaco_DD(index_ox,:)=0; Jaco_DD(index_ox,index_ox)=1; Jaco_DD(index_ox,index_si)=-1;
                        res_DD(index_si,1)=res_DD(index_si,1)+res_DD(index_ox,1);
                        res_DD(index_ox,1)=0;
                    elseif ~isequal(sort(nonzeros(Jaco_DD(index_ox,:))),[-1; 1]) % 점 3개가 접해있는 지역에서 중복적으로 합쳐지는것을 막기위해 추가함. jaco의 ox-channel행이 -1과1이라면 이미 다른행과 합쳐진 것이므로 pass.
                        Jaco_DD(index_si,:)=Jaco_DD(index_si,:)+Jaco_DD(index_ox,:);
                        Jaco_DD(index_ox,:)=0; Jaco_DD(index_ox,index_ox)=1; Jaco_DD(index_ox,index_si)=-1;
                        res_DD(index_si,1)=res_DD(index_si,1)+res_DD(index_ox,1);
                        res_DD(index_ox,1)=0;
                    end
                end
            end
        end
        % source/drain 값을 channel에 더해주기. Si부분이니까 electron/hole density도 모두 바꾸어주어야함.
        for ii=[2 4]
            for j=3
                clearvars Vertex_interface_tmp
                Vertex_interface_tmp=Vertex_interface{ii,j};
                for k=1:size(Vertex_interface_tmp,1)
                    % index 설정
                    index_source_drain_potential=find(Table_Jaco(:,1)==ii & Table_Jaco(:,2)==Vertex_interface_tmp(k,1) & Table_Jaco(:,3)==1);
                    index_source_drain_electron=index_source_drain_potential+1;
                    index_source_drain_hole=index_source_drain_potential+2;
                    index_channel_potential=find(Table_Jaco(:,1)==j & Table_Jaco(:,2)==Vertex_interface_tmp(k,1) & Table_Jaco(:,3)==1);
                    index_channel_electron=index_channel_potential+1;
                    index_channel_hole=index_channel_potential+2;
                    % Jaco 변경, potential/electron/hole 순서
                    Jaco_DD(index_channel_potential,:)=Jaco_DD(index_channel_potential,:)+Jaco_DD(index_source_drain_potential,:);
                    Jaco_DD(index_channel_electron,:)=Jaco_DD(index_channel_electron,:)+Jaco_DD(index_source_drain_electron,:);
                    Jaco_DD(index_channel_hole,:)=Jaco_DD(index_channel_hole,:)+Jaco_DD(index_source_drain_hole,:);
                    Jaco_DD(index_source_drain_potential,:)=0; Jaco_DD(index_source_drain_potential,index_source_drain_potential)=1; Jaco_DD(index_source_drain_potential,index_channel_potential)=-1;
                    Jaco_DD(index_source_drain_electron,:)=0; Jaco_DD(index_source_drain_electron,index_source_drain_electron)=1; Jaco_DD(index_source_drain_electron,index_channel_electron)=-1;
                    Jaco_DD(index_source_drain_hole,:)=0; Jaco_DD(index_source_drain_hole,index_source_drain_hole)=1; Jaco_DD(index_source_drain_hole,index_channel_hole)=-1;

                    % res 변경, potential/electron/hole 순서
                    res_DD(index_channel_potential,1)=res_DD(index_channel_potential,1)+res_DD(index_source_drain_potential,1);
                    res_DD(index_channel_electron,1)=res_DD(index_channel_electron,1)+res_DD(index_source_drain_electron,1);
                    res_DD(index_channel_hole,1)=res_DD(index_channel_hole,1)+res_DD(index_source_drain_hole,1);
                    res_DD(index_source_drain_potential,1)=0;
                    res_DD(index_source_drain_electron,1)=0;
                    res_DD(index_source_drain_hole,1)=0;
                end
            end
        end

        % Dirichlet_BC
        for ii=1:size(Dirichlet_BC,1)
            BC_index=find(Table_Jaco(:,2)==Dirichlet_BC(ii,1) & Table_Jaco(:,3)==1);
            Jaco_DD(BC_index,:)=0;
            Jaco_DD(BC_index,BC_index)=1;
            res_DD(BC_index,1)=solution_vector_DD(BC_index,Newton_DD)-V_gate;
        end

        % Source_BC
        for ii=1:size(Source_drain_BC,1)/2
            BC_index=find((Table_Jaco(:,1)==2 | Table_Jaco(:,1)==4) & Table_Jaco(:,2)==Source_drain_BC(ii,1) & Table_Jaco(:,3)==1);
            Jaco_DD(BC_index,:)=0; Jaco_DD(BC_index+1,:)=0; Jaco_DD(BC_index+2,:)=0;
            Jaco_DD(BC_index,BC_index)=1; Jaco_DD(BC_index+1,BC_index+1)=1; Jaco_DD(BC_index+2,BC_index+2)=1;
            res_DD(BC_index,1)=solution_vector_DD(BC_index,Newton_DD)-V_T*log(Nd/n_int);
            res_DD(BC_index+1,1)=solution_vector_DD(BC_index+1,Newton_DD)-Nd;
            res_DD(BC_index+2,1)=solution_vector_DD(BC_index+2,Newton_DD)-n_int^2/Nd;
        end

        % Darin_BC
        for ii=size(Source_drain_BC,1)/2+1:size(Source_drain_BC,1)
            BC_index=find((Table_Jaco(:,1)==2 | Table_Jaco(:,1)==4) & Table_Jaco(:,2)==Source_drain_BC(ii,1) & Table_Jaco(:,3)==1);
            Jaco_DD(BC_index,:)=0; Jaco_DD(BC_index+1,:)=0; Jaco_DD(BC_index+2,:)=0;
            Jaco_DD(BC_index,BC_index)=1; Jaco_DD(BC_index+1,BC_index+1)=1; Jaco_DD(BC_index+2,BC_index+2)=1;
            res_DD(BC_index,1)=solution_vector_DD(BC_index,Newton_DD)-V_T*log(Nd/n_int)-Vd;
            res_DD(BC_index+1,1)=solution_vector_DD(BC_index+1,Newton_DD)-Nd;
            res_DD(BC_index+2,1)=solution_vector_DD(BC_index+2,Newton_DD)-n_int^2/Nd;
        end

        res_save_DD(:,Newton_DD)=res_DD(:,1);

        % Scaling
        Cvector=zeros(size(solution_vector_DD,1),1); v=[V_T; max(abs(N_dop)); max(abs(N_dop))];
        for ii=1:size(solution_vector_DD,1)
            Cvector(ii,1)=v(Table_Jaco(ii,3),1);
        end
        Cmatrix_DD=spdiags(Cvector,0,size(solution_vector_DD,1),size(solution_vector_DD,1));
        Jaco_scaled_DD=Jaco_DD*Cmatrix_DD;
        Rvector_DD=1./sum(abs(Jaco_scaled_DD),2);
        Rmatrix_DD=spdiags(Rvector_DD,0,size(solution_vector_DD,1),size(solution_vector_DD,1));
        Jaco_scaled_DD=Rmatrix_DD*Jaco_scaled_DD;
        res_scaled_DD=Rmatrix_DD*res_DD;
        update_scaled_DD=Jaco_scaled_DD\(-res_scaled_DD);
        update_DD(:,Newton_DD)=Cmatrix_DD*update_scaled_DD;
        solution_vector_DD(:,Newton_DD+1)=solution_vector_DD(:,Newton_DD)+update_DD(:,Newton_DD);
        
        if Newton_DD==iteration
            solution_vector_Ramping_save(:,bias)=solution_vector_DD(:,Newton_DD);
        end
        %     % non-Scaling
        %     update(:,Newton)=Jaco\-res;
        %     solution_vector(:,Newton+1)=solution_vector(:,Newton)+update(:,Newton);

    end
    %current
end

%% Visualize
Visual_solution_vector_DD=zeros(size(Vertex,1),3);
Visual_initial_vector_DD=zeros(size(Vertex,1),3);
for ii=1:size(Table_Jaco,1)
    Visual_solution_vector_DD(Table_Jaco(ii,2),Table_Jaco(ii,3))=solution_vector_Ramping_save(ii,11);
    Visual_initial_vector_DD(Table_Jaco(ii,2),Table_Jaco(ii,3))=solution_vector_DD(ii,1);
end

dif_DD=abs(Visual_initial_vector_DD-Visual_solution_vector_DD);
plot3(Vertex(:,1),Vertex(:,2),Visual_solution_vector_DD(:,1),'*')
xlim([0 30*1e-9])


figure(2);
patch('Faces',Element,'Vertices',Vertex(:,[1 2]),'EdgeColor','black','FaceColor','none','LineWidth',1);
title('Mesh')
xlim([0 30*1e-9])
xlabel('X');
ylabel('Y');
hold on
patch('Faces',edge_interface,'Vertices',Vertex(:,[1 2]),'EdgeColor','green','FaceColor','none','LineWidth',3);
hold on
patch('Faces',Element_si_source,'Vertices',Vertex(:,[1 2]), 'EdgeColor','black','FaceColor','yellow','LineWidth',0.5)
hold on
patch('Faces',Element_si_drain,'Vertices',Vertex(:,[1 2]), 'EdgeColor','black','FaceColor','yellow','LineWidth',0.5)
hold on
patch('Faces',Element_si_channel,'Vertices',Vertex(:,[1 2]), 'EdgeColor','black','FaceColor','red','LineWidth',0.5)
hold on
patch('Faces',Element_ox1,'Vertices',Vertex(:,[1 2]), 'EdgeColor','black','FaceColor','blue','LineWidth',0.5)
hold on
patch('Faces',Element_ox2,'Vertices',Vertex(:,[1 2]), 'EdgeColor','black','FaceColor','blue','LineWidth',0.5)
hold on
patch('Faces',Element_ox2,'Vertices',Vertex(:,[1 2]), 'EdgeColor','black','FaceColor','blue','LineWidth',0.5)
hold off

figure(4);
patch('Faces',Element,'Vertices',Vertex(:,[1 2]), 'FaceVertexCData',Visual_solution_vector_DD(:,1), 'EdgeColor','black','FaceColor','interp');
title('Initial potential')
hold on
patch('Faces',edge_interface,'Vertices',Vertex(:,[1 2]),'EdgeColor','c','FaceColor','none','LineWidth',5);
hold off
xlim([0 30*1e-9])
xlabel('X');
ylabel('Y');
colorbar

figure(5);
patch('Faces',Element,'Vertices',Vertex(:,[1 2]), 'FaceVertexCData',Visual_solution_vector_DD(:,2), 'EdgeColor','black','FaceColor','interp');
title('elctron')
xlim([0 30*1e-9])
xlabel('X');
ylabel('Y');
colorbar

figure(6);
patch('Faces',Element,'Vertices',Vertex(:,[1 2]), 'FaceVertexCData',Visual_solution_vector_DD(:,3), 'EdgeColor','black','FaceColor','interp');
title('hole')
xlim([0 30*1e-9])
xlabel('X');
ylabel('Y');
colorbar

figure(7)
semilogy(abs(update(175,:)))
xlabel('Iteration');
ylabel('Maximum potential update (V)');
title('Nonlinear-possion')
xlim([0 15])
%
% % figure(7);
% % patch('Faces',Element,'Vertices',Vertex(:,[1 2]), 'FaceVertexCData',dif_DD(:,1), 'EdgeColor','black','FaceColor','interp');
% % title('Initial potential')
% % xlabel('X');
% % ylabel('Y');
% % colorbar
% %
% % figure(8);
% % patch('Faces',Element,'Vertices',Vertex(:,[1 2]), 'FaceVertexCData',dif_DD(:,2), 'EdgeColor','black','FaceColor','interp');
% % title('elctron')
% % xlabel('X');
% % ylabel('Y');
% % colorbar
% %
% % figure(9);
% % patch('Faces',Element,'Vertices',Vertex(:,[1 2]), 'FaceVertexCData',dif_DD(:,3), 'EdgeColor','black','FaceColor','interp');
% % title('hole')
% % xlabel('X');
% % ylabel('Y');
% % colorbar