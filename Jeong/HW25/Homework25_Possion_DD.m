clear; close all; clc;
mkdir('./data')

%% Variable %%
eps=[3.9; 11.7; 11.7; 11.7; 3.9]; % Epsilon=[eps1; eps2; eps3]
q=1.602192e-19; %Elementary charge,C
eps0=8.854187817e-12; %Vacuum permittivity
k=1.3806488e-23; % Bolzmann constant, J/K
n_int=1e16; % 1.0e10 /cm^-3
T=300; % Temperture, K
V_T=k*T/q; % Thermal voltage
V_gate_workfunction=0.33374; % Voltage, si_ox workfunction.
Na=1e23; % /m^3
Nd=5e26; % /m^3
N_dop=[0; Nd; -Na; Nd; 0]; % -1e18 /cm^3, 행을 region별 도핑 농도
nm=10^-9; %길이단위
mobility_n = 1417e-4;  % electron mobility
mobility_p = 470.5e-4;     % Hole mobility
width=1e-6;

coeff=q/eps0;
coeff_Jn=q*mobility_n*V_T;
coeff_Jp=-q*mobility_p*V_T;

Vdd=0;
R1=100;
V_gate=0;

%% Build a Structure %%
% Vertex 값 작성
length_x=50;
length_y=14;
dx=1;
dy=0.5;
Vertex=vertex_doublegate(0,0, length_x,length_y, dx, dy);    % vertex = vertex_doublegate(x1, y1, x2, y2, dx, dy), [nm]

n=1;
for j=1:length_y/dy+1
    for ii=1:length_x/dx+1
        vertex_position(j,ii)=n;
        n=n+1;
    end
end
vertex_position=flip(vertex_position,1);

% Dirichlet_BC 설정
Dirichlet_BC_tmp=[Contact(15,0,35,0,0); Contact(15,14,35,14,0)];
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

for ii=1:size(Dirichlet_BC_tmp,1)
    Element_Dirichlet_BC(ii,:)=Dirichlet_BC_tmp(ii,1:size(Dirichlet_BC_tmp,2)-1);
end

% Source_BC 설정
Source_BC_tmp=[Contact(0,2, 0,12, V_T*log(Nd/n_int))];
m=0; n=0; % element로 이루어진 경계조건 계산하기 위해 변환
for ii=1:size(Source_BC_tmp,1)
    for j=size(Source_BC_tmp,2):-1:1
        if ~isnan(Source_BC_tmp(ii,j))
            for k=1:j-1
                Source_BC(k+n,:)=[Source_BC_tmp(ii,k) Source_BC_tmp(ii,j)];
                m=m+1;
            end
            n=m;
            break;
        end
    end
end

for ii=1:size(Source_BC_tmp,1)
    Element_source(ii,:)=Source_BC_tmp(ii,1:size(Source_BC_tmp,2)-1);
end

% Drain_BC 설정
Drain_BC_tmp=[Contact(50,2, 50,12, V_T*log(Nd/n_int))];
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

for ii=1:size(Drain_BC_tmp,1)
    Element_drain(ii,:)=Drain_BC_tmp(ii,1:size(Drain_BC_tmp,2)-1);
end

% Element 값 입력
Element_ox1 = element(0,0, 50,2); % Element = element(x1, y1, x2, y2), [nm]
Element_si_source = element(0,2, 14,12);
Element_si_channel = element(14,2, 36,12);
Element_si_drain = element(36,2, 50,12);
Element_ox2 = element(0,12, 50,14);
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
        if edge(ii,j)<=1e-15
            edge(ii,j)=0;
        end
    end
end

%% Jacobian Table 생성
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

            phi(index_Channel,1)=phi(index_Source_drain,1);
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
            phi(j,1)=Dirichlet_BC(ii,2)+V_gate_workfunction;
        end
    end
end

% Source_BC
for ii=1:size(Source_BC,1)
    for j=1:size(Table_initial,1)
        if (Table_initial(j,1)==2 || Table_initial(j,1)==4) && Table_initial(j,2)==Source_BC(ii,1)
            phi(j,1)=Source_BC(ii,2);
        end
        if (Table_initial(j,1)==1 || Table_initial(j,1)==5) && Table_initial(j,2)==Source_BC(ii,1)
            phi(j,1)=Source_BC(ii,2);
        end
    end
end

% Drain_BC
for ii=1:size(Drain_BC,1)
    for j=1:size(Table_initial,1)
        if (Table_initial(j,1)==2 || Table_initial(j,1)==4) && Table_initial(j,2)==Drain_BC(ii,1)
            phi(j,1)=Drain_BC(ii,2);
        end
        if (Table_initial(j,1)==1 || Table_initial(j,1)==5) && Table_initial(j,2)==Drain_BC(ii,1)
            phi(j,1)=Drain_BC(ii,2);
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
        solution_vector_poisson(n,1)=initial_value(ii,1);
        n=n+1;
    else
        for j=1:3
            solution_vector_poisson(n,1)=initial_value(ii,j);
            n=n+1;
        end
    end
end

%% Visualize
% Visual_solution_vector_poisson=zeros(size(Vertex,1),3);
% for ii=1:size(Table_Jaco,1)
%     Visual_solution_vector_poisson(Table_Jaco(ii,2),Table_Jaco(ii,3))=solution_vector_poisson(ii,1);
% end
% 
% figure(1) % 점으로 3D
% plot3(Vertex(:,1),Vertex(:,2),Visual_solution_vector_poisson(:,1),'*')
% xlim([0 50*1e-9])
% 
% figure(2); % mesh 모양
% patch('Faces',Element,'Vertices',Vertex(:,[1 2]),'EdgeColor','black','FaceColor','none','LineWidth',1);
% title('Mesh')
% xlim([0 50*1e-9])
% xlabel('X');
% ylabel('Y');
% hold on
% patch('Faces',edge_interface,'Vertices',Vertex(:,[1 2]),'EdgeColor','green','FaceColor','none','LineWidth',3);
% hold on
% patch('Faces',Element_si_source,'Vertices',Vertex(:,[1 2]), 'EdgeColor','black','FaceColor','yellow','LineWidth',0.5)
% hold on
% patch('Faces',Element_si_drain,'Vertices',Vertex(:,[1 2]), 'EdgeColor','black','FaceColor','yellow','LineWidth',0.5)
% hold on
% patch('Faces',Element_si_channel,'Vertices',Vertex(:,[1 2]), 'EdgeColor','black','FaceColor','red','LineWidth',0.5)
% hold on
% patch('Faces',Element_ox1,'Vertices',Vertex(:,[1 2]), 'EdgeColor','black','FaceColor','blue','LineWidth',0.5)
% hold on
% patch('Faces',Element_ox2,'Vertices',Vertex(:,[1 2]), 'EdgeColor','black','FaceColor','blue','LineWidth',0.5)
% hold on
% patch('Faces',Element_ox2,'Vertices',Vertex(:,[1 2]), 'EdgeColor','black','FaceColor','blue','LineWidth',0.5)
% hold on
% patch('Faces',Element_source,'Vertices',Vertex(:,[1 2]), 'EdgeColor','m','FaceColor','none','LineWidth',5)
% hold on
% patch('Faces',Element_drain,'Vertices',Vertex(:,[1 2]), 'EdgeColor','m','FaceColor','none','LineWidth',5)
% hold on
% patch('Faces',Element_Dirichlet_BC,'Vertices',Vertex(:,[1 2]), 'EdgeColor','g','FaceColor','none','LineWidth',5)
% hold off
% 
% figure(3); % potential
% patch('Faces',Element,'Vertices',Vertex(:,[1 2]), 'FaceVertexCData',Visual_solution_vector_poisson(:,1), 'EdgeColor','black','FaceColor','interp');
% title('Initial potential')
% hold on
% patch('Faces',edge_interface,'Vertices',Vertex(:,[1 2]),'EdgeColor','black','FaceColor','none','LineWidth',5);
% hold on
% patch('Faces',Element_source,'Vertices',Vertex(:,[1 2]), 'EdgeColor','m','FaceColor','none','LineWidth',5)
% hold on
% patch('Faces',Element_drain,'Vertices',Vertex(:,[1 2]), 'EdgeColor','m','FaceColor','none','LineWidth',5)
% hold on
% patch('Faces',Element_Dirichlet_BC,'Vertices',Vertex(:,[1 2]), 'EdgeColor','g','FaceColor','none','LineWidth',5)
% hold off
% xlim([0 50*1e-9])
% xlabel('X');
% ylabel('Y');
% colorbar
% 
% figure(4); % electron
% patch('Faces',Element,'Vertices',Vertex(:,[1 2]), 'FaceVertexCData',Visual_solution_vector_poisson(:,2), 'EdgeColor','black','FaceColor','interp');
% title('elctron')
% xlim([0 50*1e-9])
% xlabel('X');
% ylabel('Y');
% colorbar
% 
% figure(5); % hole
% patch('Faces',Element,'Vertices',Vertex(:,[1 2]), 'FaceVertexCData',Visual_solution_vector_poisson(:,3), 'EdgeColor','black','FaceColor','interp');
% title('hole')
% xlim([0 50*1e-9])
% xlabel('X');
% ylabel('Y');
% colorbar

%% Newton Method
solution_vector_poisson_saved(:,1)=solution_vector_poisson(:,1);
coeff=q/eps0;     % m^2 -> nm^2으로 단위 변경, q/eps0 해줌

% Jacobian matrix / res vector
Iteration_poisson=40;
for Newton_poisson=1:Iteration_poisson
    fprintf("Newton Method, Newton_poisson=%d\n" , Newton_poisson)
%     Jaco_poisson=sparse(zeros(size(solution_vector_poisson,1),size(solution_vector_poisson,1))); res_poisson=zeros(size(solution_vector_poisson,1),1);
    Jaco_poisson=zeros(size(solution_vector_poisson,1),size(solution_vector_poisson,1)); res_poisson=zeros(size(solution_vector_poisson,1),1); %% 희소행렬 아님. Jaco matrix 확인용.
    
% Potential %
    for ii=1:length(Element)
        %% 각종 변수 사전계산
        % Control_Volume 계산
        Control_Volume=zeros(3,1);
        Control_Volume(1,1)=(0.25)*L(ii,3)*edge(ii,3)+(0.25)*L(ii,1)*edge(ii,1);
        Control_Volume(2,1)=(0.25)*L(ii,1)*edge(ii,1)+(0.25)*L(ii,2)*edge(ii,2);
        Control_Volume(3,1)=(0.25)*L(ii,2)*edge(ii,2)+(0.25)*L(ii,3)*edge(ii,3);

        % 현재 region의 eps_r 불러오기
        eps_now=eps(Table_element_region(ii,1),1);

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

            % potential
            x1_poisson=solution_vector_poisson(index(1,1),1);
            x2_poisson=solution_vector_poisson(index(4,1),1);
            x3_poisson=solution_vector_poisson(index(7,1),1);
            x12_poisson=x1_poisson-x2_poisson; x21_poisson=-x12_poisson;
            x23_poisson=x2_poisson-x3_poisson; x32_poisson=-x23_poisson;
            x31_poisson=x3_poisson-x1_poisson; x13_poisson=-x31_poisson;
            x12_poisson_avr=(x1_poisson+x2_poisson)/2; x23_poisson_avr=(x2_poisson+x3_poisson)/2; x31_poisson_avr=(x3_poisson+x1_poisson)/2;
            x21_poisson_avr=x12_poisson_avr; x32_poisson_avr=x23_poisson_avr; x13_poisson_avr=x31_poisson_avr;


            Jaco_tmp_poisson=zeros(9,9);
            Jaco_tmp_poisson(1,1)=eps_now*(-edge(ii,1)/L(ii,1)-edge(ii,3)/L(ii,3));
            Jaco_tmp_poisson(1,2)=-coeff*Control_Volume(1,1);
            Jaco_tmp_poisson(1,3)=coeff*Control_Volume(1,1);
            Jaco_tmp_poisson(1,4)=eps_now*edge(ii,1)/L(ii,1);
            Jaco_tmp_poisson(1,7)=eps_now*edge(ii,3)/L(ii,3);

            Jaco_tmp_poisson(2,1)=-1/V_T*n_int*exp(x1_poisson/V_T);
            Jaco_tmp_poisson(2,2)=1;

            Jaco_tmp_poisson(3,1)=1/V_T*n_int*exp(-x1_poisson/V_T);
            Jaco_tmp_poisson(3,3)=1;

            Jaco_tmp_poisson(4,1)=eps_now*edge(ii,1)/L(ii,1);
            Jaco_tmp_poisson(4,4)=eps_now*(-edge(ii,1)/L(ii,1)-edge(ii,2)/L(ii,2));
            Jaco_tmp_poisson(4,5)=-coeff*Control_Volume(2,1);
            Jaco_tmp_poisson(4,6)=coeff*Control_Volume(2,1);
            Jaco_tmp_poisson(4,7)=eps_now*edge(ii,2)/L(ii,2);

            Jaco_tmp_poisson(5,4)=-1/V_T*n_int*exp(x2_poisson/V_T);
            Jaco_tmp_poisson(5,5)=1;

            Jaco_tmp_poisson(6,4)=1/V_T*n_int*exp(-x2_poisson/V_T);
            Jaco_tmp_poisson(6,6)=1;

            Jaco_tmp_poisson(7,1)=eps_now*edge(ii,3)/L(ii,3);
            Jaco_tmp_poisson(7,4)=eps_now*edge(ii,2)/L(ii,2);
            Jaco_tmp_poisson(7,7)=eps_now*(-edge(ii,2)/L(ii,2)-edge(ii,3)/L(ii,3));
            Jaco_tmp_poisson(7,8)=-coeff*Control_Volume(3,1);
            Jaco_tmp_poisson(7,9)=coeff*Control_Volume(3,1);

            Jaco_tmp_poisson(8,7)=-1/V_T*n_int*exp(x3_poisson/V_T);
            Jaco_tmp_poisson(8,8)=1;

            Jaco_tmp_poisson(9,7)=1/V_T*n_int*exp(-x3_poisson/V_T);
            Jaco_tmp_poisson(9,9)=1;

            for j=1:9
                for k=1:9
                    Jaco_poisson(index(j,1),index(k,1))=Jaco_poisson(index(j,1),index(k,1))+Jaco_tmp_poisson(j,k);
                end
            end
            % n,p 값은 변하지 않는다.
            for j=1:3
                Jaco_poisson(index(3*j-1),:)=0;
                Jaco_poisson(index(3*j),:)=0;
                Jaco_poisson(index(3*j-1),index(3*j-2))=-1/V_T*n_int*exp(solution_vector_poisson(index(3*j-2,1),1)/V_T);
                Jaco_poisson(index(3*j-1),index(3*j-1))=1;
                Jaco_poisson(index(3*j),index(3*j-2))=1/V_T*n_int*exp(-solution_vector_poisson(index(3*j-2,1),1)/V_T);
                Jaco_poisson(index(3*j),index(3*j))=1;
            end

            % res vector
            res_tmp_poisson=zeros(9,1);
            res_potential_tmp1(1,1)=eps_now*(edge(ii,1)/L(ii,1)*(solution_vector_poisson(index(4,1),1)-solution_vector_poisson(index(1,1),1))+edge(ii,3)/L(ii,3)*(solution_vector_poisson(index(7,1),1)-solution_vector_poisson(index(1,1),1)));
            res_potential_tmp1(2,1)=eps_now*(edge(ii,2)/L(ii,2)*(solution_vector_poisson(index(7,1),1)-solution_vector_poisson(index(4,1),1))+edge(ii,1)/L(ii,1)*(solution_vector_poisson(index(1,1),1)-solution_vector_poisson(index(4,1),1)));
            res_potential_tmp1(3,1)=eps_now*(edge(ii,3)/L(ii,3)*(solution_vector_poisson(index(1,1),1)-solution_vector_poisson(index(7,1),1))+edge(ii,2)/L(ii,2)*(solution_vector_poisson(index(4,1),1)-solution_vector_poisson(index(7,1),1)));
            for j=1:3
                res_potential_tmp2(j,1)=(N_dop(Table_element_region(ii,1),1)-solution_vector_poisson(index(3*j-1,1),1)+solution_vector_poisson(index(3*j,1),1))*coeff*Control_Volume(j,1);
            end
            res_potential=res_potential_tmp1+res_potential_tmp2;
            for j=1:3
                res_electron(j,1)=solution_vector_poisson(index(3*j-1,1),1)-n_int*exp(solution_vector_poisson(index(3*j-2,1),1)/V_T);
            end
            for j=1:3
                res_hole(j,1)=solution_vector_poisson(index(3*j,1),1)-n_int*exp(-solution_vector_poisson(index(3*j-2,1),1)/V_T);
            end
            n=1;
            for j=1:3

                res_tmp_poisson(n,1)=res_potential(j,1);
                res_tmp_poisson(n+1,1)=res_electron(j,1);
                res_tmp_poisson(n+2,1)=res_hole(j,1);
                n=n+3;
            end

            for j=1:9
                res_poisson(index(j,1),1)=res_poisson(index(j,1),1)+res_tmp_poisson(j,1);
            end
            %n,p값은 변하지 않는다.
            for j=1:3
                res_poisson(index(3*j-1,1),1)=res_electron(j,1);
                res_poisson(index(3*j,1),1)=res_hole(j,1);
            end

        else
            % index : 현재 element에 속해있는 vertex가 resdue에 몇번째 있는가 저장
            index=zeros(3,1);
            n=1;
            for j=1:3
                index(j,1)=find(Table_Jaco(:,1)==Table_element_region(ii,1) & Table_Jaco(:,2)==Element(ii,j) & Table_Jaco(:,3)==1);
            end

            % Jaco_potential matrix
            Jaco_tmp_poisson=zeros(3,3);
            Jaco_tmp_poisson=[(-edge(ii,1)/L(ii,1)-edge(ii,3)/L(ii,3)) edge(ii,1)/L(ii,1) edge(ii,3)/L(ii,3);
                edge(ii,1)/L(ii,1) -edge(ii,1)/L(ii,1)-edge(ii,2)/L(ii,2) edge(ii,2)/L(ii,2);
                edge(ii,3)/L(ii,3) edge(ii,2)/L(ii,2) -edge(ii,2)/L(ii,2)-edge(ii,3)/L(ii,3)];
            Jaco_tmp_poisson=eps_now*Jaco_tmp_poisson;
            for j=1:3
                for k=1:3
                    Jaco_poisson(index(j,1),index(k,1))=Jaco_poisson(index(j,1),index(k,1))+Jaco_tmp_poisson(j,k);
                end
            end

            % res vector
            res_tmp_poisson=zeros(3,1);
            res_tmp_poisson=eps_now.*[(-edge(ii,1)/L(ii,1)-edge(ii,3)/L(ii,3))*solution_vector_poisson(index(1,1),1)+edge(ii,1)/L(ii,1)*solution_vector_poisson(index(2,1),1)+edge(ii,3)/L(ii,3)*solution_vector_poisson(index(3,1),1);
                edge(ii,1)/L(ii,1)*solution_vector_poisson(index(1,1),1)+(-edge(ii,1)/L(ii,1)-edge(ii,2)/L(ii,2))*solution_vector_poisson(index(2,1),1)+edge(ii,2)/L(ii,2)*solution_vector_poisson(index(3,1),1);
                edge(ii,3)/L(ii,3)*solution_vector_poisson(index(1,1),1)+edge(ii,2)/L(ii,2)*solution_vector_poisson(index(2,1),1)+(-edge(ii,2)/L(ii,2)-edge(ii,3)/L(ii,3))*solution_vector_poisson(index(3,1),1)];
            for j=1:3
                res_poisson(index(j,1),1)=res_poisson(index(j,1),1)+res_tmp_poisson(j,1);
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
                if length(sort(nonzeros(Jaco_poisson(index_ox,:))))>2
                    Jaco_poisson(index_si,:)=Jaco_poisson(index_si,:)+Jaco_poisson(index_ox,:);
                    Jaco_poisson(index_ox,:)=0; Jaco_poisson(index_ox,index_ox)=1; Jaco_poisson(index_ox,index_si)=-1;
                    res_poisson(index_si,1)=res_poisson(index_si,1)+res_poisson(index_ox,1);
                    res_poisson(index_ox,1)=0;
                elseif ~isequal(sort(nonzeros(Jaco_poisson(index_ox,:))),[-1; 1]) % 점 3개가 접해있는 지역에서 중복적으로 합쳐지는것을 막기위해 추가함. jaco의 ox-channel행이 -1과1이라면 이미 다른행과 합쳐진 것이므로 pass.
                    Jaco_poisson(index_si,:)=Jaco_poisson(index_si,:)+Jaco_poisson(index_ox,:);
                    Jaco_poisson(index_ox,:)=0; Jaco_poisson(index_ox,index_ox)=1; Jaco_poisson(index_ox,index_si)=-1;
                    res_poisson(index_si,1)=res_poisson(index_si,1)+res_poisson(index_ox,1);
                    res_poisson(index_ox,1)=0;
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
                Jaco_poisson(index_channel_potential,:)=Jaco_poisson(index_channel_potential,:)+Jaco_poisson(index_source_drain_potential,:);
                Jaco_poisson(index_channel_electron,:)=Jaco_poisson(index_channel_electron,:)+Jaco_poisson(index_source_drain_electron,:);
                Jaco_poisson(index_channel_hole,:)=Jaco_poisson(index_channel_hole,:)+Jaco_poisson(index_source_drain_hole,:);
                Jaco_poisson(index_source_drain_potential,:)=0; Jaco_poisson(index_source_drain_potential,index_source_drain_potential)=1; Jaco_poisson(index_source_drain_potential,index_channel_potential)=-1;
                Jaco_poisson(index_source_drain_electron,:)=0; Jaco_poisson(index_source_drain_electron,index_source_drain_electron)=1; Jaco_poisson(index_source_drain_electron,index_channel_electron)=-1;
                Jaco_poisson(index_source_drain_hole,:)=0; Jaco_poisson(index_source_drain_hole,index_source_drain_hole)=1; Jaco_poisson(index_source_drain_hole,index_channel_hole)=-1;

                % res 변경, potential/electron/hole 순서
                res_poisson(index_channel_potential,1)=res_poisson(index_channel_potential,1)+res_poisson(index_source_drain_potential,1);
                res_poisson(index_channel_electron,1)=res_poisson(index_channel_electron,1)+res_poisson(index_source_drain_electron,1);
                res_poisson(index_channel_hole,1)=res_poisson(index_channel_hole,1)+res_poisson(index_source_drain_hole,1);
                res_poisson(index_source_drain_potential,1)=0;
                res_poisson(index_source_drain_electron,1)=0;
                res_poisson(index_source_drain_hole,1)=0;
            end
        end
    end

    % Dirichlet_BC
    for ii=1:size(Dirichlet_BC,1)
        BC_index=find(Table_Jaco(:,2)==Dirichlet_BC(ii,1) & Table_Jaco(:,3)==1);
        Jaco_poisson(BC_index,:)=0;
        Jaco_poisson(BC_index,BC_index)=1;
        res_poisson(BC_index,1)=solution_vector_poisson(BC_index,1)-V_gate_workfunction;
    end

    % Source_BC
    for ii=1:size(Source_BC,1)
        BC_index=find((Table_Jaco(:,1)==2 | Table_Jaco(:,1)==4) & Table_Jaco(:,2)==Source_BC(ii,1) & Table_Jaco(:,3)==1);
        Jaco_poisson(BC_index,:)=0; Jaco_poisson(BC_index+1,:)=0; Jaco_poisson(BC_index+2,:)=0;
        Jaco_poisson(BC_index,BC_index)=1; Jaco_poisson(BC_index+1,BC_index+1)=1; Jaco_poisson(BC_index+2,BC_index+2)=1;
        res_poisson(BC_index,1)=solution_vector_poisson(BC_index,1)-V_T*log(Nd/n_int);
        res_poisson(BC_index+1,1)=solution_vector_poisson(BC_index+1,1)-Nd;
        res_poisson(BC_index+2,1)=solution_vector_poisson(BC_index+2,1)-n_int^2/Nd;
    end

    % Drain_BC
    for ii=1:size(Drain_BC,1)
        BC_index=find((Table_Jaco(:,1)==2 | Table_Jaco(:,1)==4) & Table_Jaco(:,2)==Drain_BC(ii,1) & Table_Jaco(:,3)==1);
        Jaco_poisson(BC_index,:)=0; Jaco_poisson(BC_index+1,:)=0; Jaco_poisson(BC_index+2,:)=0;
        Jaco_poisson(BC_index,BC_index)=1; Jaco_poisson(BC_index+1,BC_index+1)=1; Jaco_poisson(BC_index+2,BC_index+2)=1;
        res_poisson(BC_index,1)=solution_vector_poisson(BC_index,1)-V_T*log(Nd/n_int);
        res_poisson(BC_index+1,1)=solution_vector_poisson(BC_index+1,1)-Nd;
        res_poisson(BC_index+2,1)=solution_vector_poisson(BC_index+2,1)-n_int^2/Nd;
    end

    % Scaling
    Cvector_poisson=zeros(size(solution_vector_poisson,1),1); v=[V_T; max(abs(N_dop)); max(abs(N_dop))];
    for ii=1:size(solution_vector_poisson,1)
        Cvector_poisson(ii,1)=v(Table_Jaco(ii,3),1);
    end
    Jaco_poisson=sparse(Jaco_poisson);
    Cmatrix_poisson=spdiags(Cvector_poisson,0,size(solution_vector_poisson,1),size(solution_vector_poisson,1));
    Jaco_scaled_poisson=Jaco_poisson*Cmatrix_poisson;
    Rvector_poisson=1./sum(abs(Jaco_scaled_poisson),2);
    Rmatrix_poisson=spdiags(Rvector_poisson,0,size(solution_vector_poisson,1),size(solution_vector_poisson,1));
    Jaco_scaled_poisson=Rmatrix_poisson*Jaco_scaled_poisson;
    res_scaled_poisson=Rmatrix_poisson*res_poisson;
    update_scaled_poisson=Jaco_scaled_poisson\(-res_scaled_poisson);
    update_poisson(:,Newton_poisson)=Cmatrix_poisson*update_scaled_poisson;
    
    solution_vector_poisson(:,1)=solution_vector_poisson(:,1)+update_poisson(:,Newton_poisson);
    solution_vector_poisson_saved(:,Newton_poisson+1)=solution_vector_poisson(:,1);

    % break
    update_poisson_phi_for_break=update_poisson(find(Table_Jaco(:,1)==2, 1 , 'first'):3:find(Table_Jaco(:,1)==4, 1 , 'last'),Newton_poisson);
    fprintf("Max Updatevector=%d\n\n" , max(abs(update_poisson_phi_for_break)))
    if max(abs(update_poisson_phi_for_break))<1e-15
        break;
    end
end

    update_poisson_phi=update_poisson(find(Table_Jaco(:,1)==2, 1 , 'first'):3:find(Table_Jaco(:,1)==4, 1 , 'last'),:);

%% Save

save('./data/homework25_poisson.mat')

%% Visualize
% Visual_solution_vector_poisson=zeros(size(Vertex,1),3);
% for ii=1:size(Table_Jaco,1)
%     Visual_solution_vector_poisson(Table_Jaco(ii,2),Table_Jaco(ii,3))=solution_vector_poisson_saved(ii,Newton_poisson);
% end
% 
% figure(1) % 점으로 3D
% plot3(Vertex(:,1),Vertex(:,2),Visual_solution_vector_poisson(:,1),'*')
% xlim([0 50*1e-9])
% 
% figure(2); % mesh 모양
% patch('Faces',Element,'Vertices',Vertex(:,[1 2]),'EdgeColor','black','FaceColor','none','LineWidth',1);
% title('Mesh')
% xlim([0 50*1e-9])
% xlabel('X');
% ylabel('Y');
% hold on
% patch('Faces',edge_interface,'Vertices',Vertex(:,[1 2]),'EdgeColor','green','FaceColor','none','LineWidth',3);
% hold on
% patch('Faces',Element_si_source,'Vertices',Vertex(:,[1 2]), 'EdgeColor','black','FaceColor','yellow','LineWidth',0.5)
% hold on
% patch('Faces',Element_si_drain,'Vertices',Vertex(:,[1 2]), 'EdgeColor','black','FaceColor','yellow','LineWidth',0.5)
% hold on
% patch('Faces',Element_si_channel,'Vertices',Vertex(:,[1 2]), 'EdgeColor','black','FaceColor','red','LineWidth',0.5)
% hold on
% patch('Faces',Element_ox1,'Vertices',Vertex(:,[1 2]), 'EdgeColor','black','FaceColor','blue','LineWidth',0.5)
% hold on
% patch('Faces',Element_ox2,'Vertices',Vertex(:,[1 2]), 'EdgeColor','black','FaceColor','blue','LineWidth',0.5)
% hold on
% patch('Faces',Element_ox2,'Vertices',Vertex(:,[1 2]), 'EdgeColor','black','FaceColor','blue','LineWidth',0.5)
% hold on
% patch('Faces',Element_source,'Vertices',Vertex(:,[1 2]), 'EdgeColor','m','FaceColor','none','LineWidth',5)
% hold on
% patch('Faces',Element_drain,'Vertices',Vertex(:,[1 2]), 'EdgeColor','m','FaceColor','none','LineWidth',5)
% hold on
% patch('Faces',Element_Dirichlet_BC,'Vertices',Vertex(:,[1 2]), 'EdgeColor','g','FaceColor','none','LineWidth',5)
% hold off
% axis equal
% 
% figure(3); % potential
% patch('Faces',Element,'Vertices',Vertex(:,[1 2]), 'FaceVertexCData',Visual_solution_vector_poisson(:,1), 'EdgeColor','black','FaceColor','interp');
% title('Initial potential')
% hold on
% patch('Faces',edge_interface,'Vertices',Vertex(:,[1 2]),'EdgeColor','black','FaceColor','none','LineWidth',5);
% hold on
% patch('Faces',Element_source,'Vertices',Vertex(:,[1 2]), 'EdgeColor','m','FaceColor','none','LineWidth',5)
% hold on
% patch('Faces',Element_drain,'Vertices',Vertex(:,[1 2]), 'EdgeColor','m','FaceColor','none','LineWidth',5)
% hold on
% patch('Faces',Element_Dirichlet_BC,'Vertices',Vertex(:,[1 2]), 'EdgeColor','g','FaceColor','none','LineWidth',5)
% hold off
% xlim([0 50*1e-9])
% xlabel('X');
% ylabel('Y');
% colorbar
% 
% figure(4); % electron
% patch('Faces',Element,'Vertices',Vertex(:,[1 2]), 'FaceVertexCData',Visual_solution_vector_poisson(:,2), 'EdgeColor','black','FaceColor','interp');
% title('elctron')
% xlim([0 50*1e-9])
% xlabel('X');
% ylabel('Y');
% colorbar
% 
% figure(5); % hole
% patch('Faces',Element,'Vertices',Vertex(:,[1 2]), 'FaceVertexCData',Visual_solution_vector_poisson(:,3), 'EdgeColor','black','FaceColor','interp');
% title('hole')
% xlim([0 50*1e-9])
% xlabel('X');
% ylabel('Y');
% colorbar

%% Drift-Diffusion Model + mixed mode simulation
% mixed mode simulation
solution_vector_circit=zeros(11,1);
Table_mixed_simulation=(size(solution_vector_poisson)+1:size(solution_vector_poisson)+size(solution_vector_circit))';
ind_Id=Table_mixed_simulation(1,1); ind_Is=Table_mixed_simulation(2,1); ind_Ig=Table_mixed_simulation(3,1);
ind_I1=Table_mixed_simulation(4,1); ind_I2=Table_mixed_simulation(5,1); ind_Vd=Table_mixed_simulation(6,1); ind_Vs=Table_mixed_simulation(7,1);
ind_Vg=Table_mixed_simulation(8,1); ind_V1=Table_mixed_simulation(9,1); ind_V2=Table_mixed_simulation(10,1); ind_Vout=Table_mixed_simulation(11,1);




Iteration_DD=10000;

% Jacobian matrix / res vector
solution_vector_DD(:,1)=[solution_vector_poisson(:,1); solution_vector_circit];
solution_vector_DD_saved(:,1)=solution_vector_DD(:,1);
for Newton_DD=1:Iteration_DD
    fprintf("Drift-Diffusion Model, Newton_DD=%d\n\n" , Newton_DD)
    clearvars Jaco_DD; clearvars res_DD; 
%     Jaco_DD=sparse(zeros(size(solution_vector_DD,1),size(solution_vector_DD,1))); res_DD=zeros(size(solution_vector_DD,1),1);
    Jaco_DD=zeros(size(solution_vector_DD,1),size(solution_vector_DD,1)); res_DD=zeros(size(solution_vector_DD,1),1); %% 희소행렬 아님. Jaco matrix 확인용.
    
% Jaco/res %
    for ii=1:length(Element)
        % Control_Volume 계산
        Control_Volume=zeros(3,1);
        Control_Volume(1,1)=(0.25)*L(ii,3)*edge(ii,3)+(0.25)*L(ii,1)*edge(ii,1);
        Control_Volume(2,1)=(0.25)*L(ii,1)*edge(ii,1)+(0.25)*L(ii,2)*edge(ii,2);
        Control_Volume(3,1)=(0.25)*L(ii,2)*edge(ii,2)+(0.25)*L(ii,3)*edge(ii,3);
        
        % 현재 region의 eps_r 불러오기
        eps_now=eps(Table_element_region(ii,1),1);

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
            x12_DD=solution_vector_DD(index(1,1),1)-solution_vector_DD(index(4,1),1); x21_DD=-x12_DD;
            x23_DD=solution_vector_DD(index(4,1),1)-solution_vector_DD(index(7,1),1); x32_DD=-x23_DD;
            x13_DD=solution_vector_DD(index(1,1),1)-solution_vector_DD(index(7,1),1); x31_DD=-x13_DD;
            n1_DD=solution_vector_DD(index(1,1)+1,1); n2_DD=solution_vector_DD(index(4,1)+1,1); n3_DD=solution_vector_DD(index(7,1)+1,1);
            p1_DD=solution_vector_DD(index(1,1)+2,1); p2_DD=solution_vector_DD(index(4,1)+2,1); p3_DD=solution_vector_DD(index(7,1)+2,1);

            % Jaco %
            Jaco_tmp_DD=zeros(9,9);
            Jaco_tmp_DD(1,1)=eps_now*(-edge(ii,1)/L(ii,1)-edge(ii,3)/L(ii,3));
            Jaco_tmp_DD(1,2)=-coeff*Control_Volume(1,1);
            Jaco_tmp_DD(1,3)=coeff*Control_Volume(1,1);
            Jaco_tmp_DD(1,4)=eps_now*edge(ii,1)/L(ii,1);
            Jaco_tmp_DD(1,7)=eps_now*edge(ii,3)/L(ii,3);

            Jaco_tmp_DD(2,1)=coeff_Jn*((edge(ii,1)/L(ii,1))*(n2_DD/(-V_T)*Ber_d(x21_DD/V_T)-n1_DD/(V_T)*Ber_d(-x21_DD/V_T))+((edge(ii,3)/L(ii,3))*(n3_DD/(-V_T)*Ber_d(x31_DD/V_T)-n1_DD/(V_T)*Ber_d(-x31_DD/V_T))));
            Jaco_tmp_DD(2,2)=coeff_Jn*((edge(ii,1)/L(ii,1))*(-Ber(-x21_DD/V_T))+((edge(ii,3)/L(ii,3))*(-Ber(-x31_DD/V_T))));
            Jaco_tmp_DD(2,4)=coeff_Jn*((edge(ii,1)/L(ii,1))*(n2_DD/(V_T)*Ber_d(x21_DD/V_T)-n1_DD/(-V_T)*Ber_d(-x21_DD/V_T)));
            Jaco_tmp_DD(2,5)=coeff_Jn*((edge(ii,1)/L(ii,1))*(Ber(x21_DD/V_T)));
            Jaco_tmp_DD(2,7)=coeff_Jn*((edge(ii,3)/L(ii,3))*(n3_DD/(V_T)*Ber_d(x31_DD/V_T)-n1_DD/(-V_T)*Ber_d(-x31_DD/V_T)));
            Jaco_tmp_DD(2,8)=coeff_Jn*((edge(ii,3)/L(ii,3))*(Ber(x31_DD/V_T)));

            Jaco_tmp_DD(3,1)=coeff_Jp*((edge(ii,1)/L(ii,1))*(p2_DD/(V_T)*Ber_d(-x21_DD/V_T)-p1_DD/(-V_T)*Ber_d(x21_DD/V_T))+((edge(ii,3)/L(ii,3))*(p3_DD/(V_T)*Ber_d(-x31_DD/V_T)-p1_DD/(-V_T)*Ber_d(x31_DD/V_T))));
            Jaco_tmp_DD(3,3)=coeff_Jp*((edge(ii,1)/L(ii,1))*(-Ber(x21_DD/V_T))+(edge(ii,3)/L(ii,3))*(-Ber(x31_DD/V_T)));
            Jaco_tmp_DD(3,4)=coeff_Jp*((edge(ii,1)/L(ii,1))*(p2_DD/(-V_T)*Ber_d(-x21_DD/V_T)-p1_DD/(V_T)*Ber_d(x21_DD/V_T)));
            Jaco_tmp_DD(3,6)=coeff_Jp*((edge(ii,1)/L(ii,1))*(Ber(-x21_DD/V_T)));
            Jaco_tmp_DD(3,7)=coeff_Jp*((edge(ii,3)/L(ii,3))*(p3_DD/(-V_T)*Ber_d(-x31_DD/V_T)-p1_DD/(V_T)*Ber_d(x31_DD/V_T)));
            Jaco_tmp_DD(3,9)=coeff_Jp*((edge(ii,3)/L(ii,3))*(Ber(-x31_DD/V_T)));

            Jaco_tmp_DD(4,1)=eps_now*edge(ii,1)/L(ii,1);
            Jaco_tmp_DD(4,4)=eps_now*(-edge(ii,1)/L(ii,1)-edge(ii,2)/L(ii,2));
            Jaco_tmp_DD(4,5)=-coeff*Control_Volume(2,1);
            Jaco_tmp_DD(4,6)=coeff*Control_Volume(2,1);
            Jaco_tmp_DD(4,7)=eps_now*edge(ii,2)/L(ii,2);

            Jaco_tmp_DD(5,1)=coeff_Jn*((edge(ii,1)/L(ii,1))*(n1_DD/(V_T)*Ber_d(x12_DD/V_T)-n2_DD/(-V_T)*Ber_d(-x12_DD/V_T)));
            Jaco_tmp_DD(5,2)=coeff_Jn*((edge(ii,1)/L(ii,1))*(Ber(x12_DD/V_T)));
            Jaco_tmp_DD(5,4)=coeff_Jn*((edge(ii,2)/L(ii,2))*(n3_DD/(-V_T)*Ber_d(x32_DD/V_T)-n2_DD/(V_T)*Ber_d(-x32_DD/V_T))+((edge(ii,1)/L(ii,1))*(n1_DD/(-V_T)*Ber_d(x12_DD/V_T)-n2_DD/(V_T)*Ber_d(-x12_DD/V_T))));
            Jaco_tmp_DD(5,5)=coeff_Jn*((edge(ii,2)/L(ii,2))*(-Ber(-x32_DD/V_T))+((edge(ii,1)/L(ii,1))*(-Ber(-x12_DD/V_T))));
            Jaco_tmp_DD(5,7)=coeff_Jn*((edge(ii,2)/L(ii,2))*(n3_DD/(V_T)*Ber_d(x32_DD/V_T)-n2_DD/(-V_T)*Ber_d(-x32_DD/V_T)));
            Jaco_tmp_DD(5,8)=coeff_Jn*((edge(ii,2)/L(ii,2))*(Ber(x32_DD/V_T)));

            Jaco_tmp_DD(6,1)=coeff_Jp*((edge(ii,1)/L(ii,1))*(p1_DD/(-V_T)*Ber_d(-x12_DD/V_T)-p2_DD/(V_T)*Ber_d(x12_DD/V_T)));
            Jaco_tmp_DD(6,3)=coeff_Jp*((edge(ii,1)/L(ii,1))*(Ber(-x12_DD/V_T)));
            Jaco_tmp_DD(6,4)=coeff_Jp*((edge(ii,2)/L(ii,2))*(p3_DD/(V_T)*Ber_d(-x32_DD/V_T)-p2_DD/(-V_T)*Ber_d(x32_DD/V_T))+((edge(ii,1)/L(ii,1))*(p1_DD/(V_T)*Ber_d(-x12_DD/V_T)-p2_DD/(-V_T)*Ber_d(x12_DD/V_T))));
            Jaco_tmp_DD(6,6)=coeff_Jp*((edge(ii,2)/L(ii,2))*(-Ber(x32_DD/V_T))+(edge(ii,1)/L(ii,1))*(-Ber(x12_DD/V_T)));
            Jaco_tmp_DD(6,7)=coeff_Jp*((edge(ii,2)/L(ii,2))*(p3_DD/(-V_T)*Ber_d(-x32_DD/V_T)-p2_DD/(V_T)*Ber_d(x32_DD/V_T)));
            Jaco_tmp_DD(6,9)=coeff_Jp*((edge(ii,2)/L(ii,2))*(Ber(-x32_DD/V_T)));

            Jaco_tmp_DD(7,1)=eps_now*edge(ii,3)/L(ii,3);
            Jaco_tmp_DD(7,4)=eps_now*edge(ii,2)/L(ii,2);
            Jaco_tmp_DD(7,7)=eps_now*(-edge(ii,2)/L(ii,2)-edge(ii,3)/L(ii,3));
            Jaco_tmp_DD(7,8)=-coeff*Control_Volume(3,1);
            Jaco_tmp_DD(7,9)=coeff*Control_Volume(3,1);

            Jaco_tmp_DD(8,1)=coeff_Jn*((edge(ii,3)/L(ii,3))*(n1_DD/(V_T)*Ber_d(x13_DD/V_T)-n3_DD/(-V_T)*Ber_d(-x13_DD/V_T)));
            Jaco_tmp_DD(8,2)=coeff_Jn*((edge(ii,3)/L(ii,3))*(Ber(x13_DD/V_T)));
            Jaco_tmp_DD(8,4)=coeff_Jn*((edge(ii,2)/L(ii,2))*(n2_DD/(V_T)*Ber_d(x23_DD/V_T)-n3_DD/(-V_T)*Ber_d(-x23_DD/V_T)));
            Jaco_tmp_DD(8,5)=coeff_Jn*((edge(ii,2)/L(ii,2))*(Ber(x23_DD/V_T)));
            Jaco_tmp_DD(8,7)=coeff_Jn*((edge(ii,3)/L(ii,3))*(n1_DD/(-V_T)*Ber_d(x13_DD/V_T)-n3_DD/(V_T)*Ber_d(-x13_DD/V_T))+((edge(ii,2)/L(ii,2))*(n2_DD/(-V_T)*Ber_d(x23_DD/V_T)-n3_DD/(V_T)*Ber_d(-x23_DD/V_T))));
            Jaco_tmp_DD(8,8)=coeff_Jn*((edge(ii,3)/L(ii,3))*(-Ber(-x13_DD/V_T))+((edge(ii,2)/L(ii,2))*(-Ber(-x23_DD/V_T))));

            Jaco_tmp_DD(9,1)=coeff_Jp*((edge(ii,3)/L(ii,3))*(p1_DD/(-V_T)*Ber_d(-x13_DD/V_T)-p3_DD/(V_T)*Ber_d(x13_DD/V_T)));
            Jaco_tmp_DD(9,3)=coeff_Jp*((edge(ii,3)/L(ii,3))*(Ber(-x13_DD/V_T)));
            Jaco_tmp_DD(9,4)=coeff_Jp*((edge(ii,2)/L(ii,2))*(p2_DD/(-V_T)*Ber_d(-x23_DD/V_T)-p3_DD/(V_T)*Ber_d(x23_DD/V_T)));
            Jaco_tmp_DD(9,6)=coeff_Jp*((edge(ii,2)/L(ii,2))*(Ber(-x23_DD/V_T)));
            Jaco_tmp_DD(9,7)=coeff_Jp*((edge(ii,3)/L(ii,3))*(p1_DD/(V_T)*Ber_d(-x13_DD/V_T)-p3_DD/(-V_T)*Ber_d(x13_DD/V_T))+((edge(ii,2)/L(ii,2))*(p2_DD/(V_T)*Ber_d(-x23_DD/V_T)-p3_DD/(-V_T)*Ber_d(x23_DD/V_T))));
            Jaco_tmp_DD(9,9)=coeff_Jp*((edge(ii,3)/L(ii,3))*(-Ber(x13_DD/V_T))+(edge(ii,2)/L(ii,2))*(-Ber(x23_DD/V_T)));


            for j=1:9
                for k=1:9
                    Jaco_DD(index(j,1),index(k,1))=Jaco_DD(index(j,1),index(k,1))+Jaco_tmp_DD(j,k);
                end
            end

            % res vector
            res_tmp_DD=zeros(9,1);
            res_potential_tmp1_DD(1,1)=eps_now*(edge(ii,1)/L(ii,1)*(x21_DD)+edge(ii,3)/L(ii,3)*(x31_DD));
            res_potential_tmp1_DD(2,1)=eps_now*(edge(ii,2)/L(ii,2)*(x32_DD)+edge(ii,1)/L(ii,1)*(x12_DD));
            res_potential_tmp1_DD(3,1)=eps_now*(edge(ii,3)/L(ii,3)*(x13_DD)+edge(ii,2)/L(ii,2)*(x23_DD));
            for j=1:3
                res_potential_tmp2_DD(j,1)=(N_dop(Table_element_region(ii,1),1)-solution_vector_DD(index(3*j-1,1),1)+solution_vector_DD(index(3*j,1),1))*coeff*Control_Volume(j,1);
            end
            
            res_potential_DD=res_potential_tmp1_DD+res_potential_tmp2_DD;
           
            % Jn
            res_Jn=zeros(3,1);
            res_Jn(1,1)=coeff_Jn*((edge(ii,1)/L(ii,1))*(n2_DD*Ber(x21_DD/V_T)-n1_DD*Ber(-x21_DD/V_T))+(edge(ii,3)/L(ii,3))*(n3_DD*Ber(x31_DD/V_T)-n1_DD*Ber(-x31_DD/V_T)));
            res_Jn(2,1)=coeff_Jn*((edge(ii,2)/L(ii,2))*(n3_DD*Ber(x32_DD/V_T)-n2_DD*Ber(-x32_DD/V_T))+(edge(ii,1)/L(ii,1))*(n1_DD*Ber(x12_DD/V_T)-n2_DD*Ber(-x12_DD/V_T)));
            res_Jn(3,1)=coeff_Jn*((edge(ii,3)/L(ii,3))*(n1_DD*Ber(x13_DD/V_T)-n3_DD*Ber(-x13_DD/V_T))+(edge(ii,2)/L(ii,2))*(n2_DD*Ber(x23_DD/V_T)-n3_DD*Ber(-x23_DD/V_T)));

            % Jp
            res_Jp=zeros(3,1);
            res_Jp(1,1)=coeff_Jp*((edge(ii,1)/L(ii,1))*(p2_DD*Ber(-x21_DD/V_T)-p1_DD*Ber(x21_DD/V_T))+(edge(ii,3)/L(ii,3))*(p3_DD*Ber(-x31_DD/V_T)-p1_DD*Ber(x31_DD/V_T)));
            res_Jp(2,1)=coeff_Jp*((edge(ii,2)/L(ii,2))*(p3_DD*Ber(-x32_DD/V_T)-p2_DD*Ber(x32_DD/V_T))+(edge(ii,1)/L(ii,1))*(p1_DD*Ber(-x12_DD/V_T)-p2_DD*Ber(x12_DD/V_T)));
            res_Jp(3,1)=coeff_Jp*((edge(ii,3)/L(ii,3))*(p1_DD*Ber(-x13_DD/V_T)-p3_DD*Ber(x13_DD/V_T))+(edge(ii,2)/L(ii,2))*(p2_DD*Ber(-x23_DD/V_T)-p3_DD*Ber(x23_DD/V_T)));

            n=1;
            for j=1:3
                res_tmp_DD(n,1)=res_potential_DD(j,1);
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
            Jaco_tmp_DD=eps_now*Jaco_tmp_DD;
            for j=1:3
                for k=1:3
                    Jaco_DD(index(j,1),index(k,1))=Jaco_DD(index(j,1),index(k,1))+Jaco_tmp_DD(j,k);
                end
            end

            % res vector
            res_tmp_DD=zeros(3,1);
            res_tmp_DD=eps_now.*[(-edge(ii,1)/L(ii,1)-edge(ii,3)/L(ii,3))*solution_vector_DD(index(1,1),1)+edge(ii,1)/L(ii,1)*solution_vector_DD(index(2,1),1)+edge(ii,3)/L(ii,3)*solution_vector_DD(index(3,1),1);
                edge(ii,1)/L(ii,1)*solution_vector_DD(index(1,1),1)+(-edge(ii,1)/L(ii,1)-edge(ii,2)/L(ii,2))*solution_vector_DD(index(2,1),1)+edge(ii,2)/L(ii,2)*solution_vector_DD(index(3,1),1);
                edge(ii,3)/L(ii,3)*solution_vector_DD(index(1,1),1)+edge(ii,2)/L(ii,2)*solution_vector_DD(index(2,1),1)+(-edge(ii,2)/L(ii,2)-edge(ii,3)/L(ii,3))*solution_vector_DD(index(3,1),1)];
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


% Boundary condition %
    % Dirichlet_BC
    for ii=1:size(Dirichlet_BC,1)
        BC_index=find(Table_Jaco(:,2)==Dirichlet_BC(ii,1) & Table_Jaco(:,3)==1);
        Jaco_DD(BC_index,:)=0;
        Jaco_DD(BC_index,BC_index)=1;
        Jaco_DD(BC_index,ind_Vg)=-1;
        res_DD(BC_index,1)=solution_vector_DD(BC_index,1)-V_gate_workfunction-V_gate;
    end

    % Source_BC
    for ii=1:size(Source_BC,1)
        BC_index=find(Table_Jaco(:,1)==2 & Table_Jaco(:,2)==Source_BC(ii,1) & Table_Jaco(:,3)==1);
        Jaco_DD(BC_index,:)=0; Jaco_DD(BC_index+1,:)=0; Jaco_DD(BC_index+2,:)=0;
        Jaco_DD(BC_index,BC_index)=1; Jaco_DD(BC_index+1,BC_index+1)=1; Jaco_DD(BC_index+2,BC_index+2)=1;
        Jaco_DD(BC_index,ind_Vs)=-1;
        res_DD(BC_index,1)=solution_vector_DD(BC_index,1)-V_T*log(Nd/n_int)-solution_vector_DD(ind_Vs,1);
        res_DD(BC_index+1,1)=solution_vector_DD(BC_index+1,1)-Nd;
        res_DD(BC_index+2,1)=solution_vector_DD(BC_index+2,1)-n_int^2/Nd;
    end

    % Darin_BC
    for ii=1:size(Drain_BC,1)
        BC_index=find(Table_Jaco(:,1)==4 & Table_Jaco(:,2)==Drain_BC(ii,1) & Table_Jaco(:,3)==1);
        
        % for mixed simulation
        Jaco_DD(ind_Id,:)=Jaco_DD(ind_Id,:)-Jaco_DD(BC_index+1,:)*width;
        Jaco_DD(ind_Id,:)=Jaco_DD(ind_Id,:)-Jaco_DD(BC_index+2,:)*width;
        res_DD(ind_Id,:)=res_DD(ind_Id,:)-res_DD(BC_index+1,:)*width;
        res_DD(ind_Id,:)=res_DD(ind_Id,:)-res_DD(BC_index+2,:)*width;
        
        % Boundary condition
        Jaco_DD(BC_index,:)=0; Jaco_DD(BC_index+1,:)=0; Jaco_DD(BC_index+2,:)=0;
        Jaco_DD(BC_index,BC_index)=1; Jaco_DD(BC_index+1,BC_index+1)=1; Jaco_DD(BC_index+2,BC_index+2)=1;
        Jaco_DD(BC_index,ind_Vd)=-1;
        res_DD(BC_index,1)=solution_vector_DD(BC_index,1)-V_T*log(Nd/n_int)-solution_vector_DD(ind_Vd,1);
        res_DD(BC_index+1,1)=solution_vector_DD(BC_index+1,1)-Nd;
        res_DD(BC_index+2,1)=solution_vector_DD(BC_index+2,1)-n_int^2/Nd;
    end

    %% mixed mode simulation
    % Id
    Jaco_DD(ind_Id,ind_Id)=1;
    res_DD(ind_Id,1)=solution_vector_DD(ind_Id,1)+res_DD(ind_Id,1);

    % Is
    Jaco_DD(ind_Is,ind_Is)=1; Jaco_DD(ind_Is,ind_Id)=1; Jaco_DD(ind_Is,ind_Ig)=1;
    res_DD(ind_Is,1)=0;

    % Ig
    Jaco_DD(ind_Ig,ind_Ig)=1;
    res_DD(ind_Ig,1)=solution_vector_DD(ind_Ig,1)-0;

    % I1
    Jaco_DD(ind_I1,ind_I1)=1; Jaco_DD(ind_I1,ind_I2)=1;
    res_DD(ind_I1,1)=0;

    % I2
    Jaco_DD(ind_I2,ind_I2)=1; Jaco_DD(ind_I2,ind_Vout)=1/R1; Jaco_DD(ind_I2,ind_V2)=-1/R1;
    res_DD(ind_I2,1)=solution_vector_DD(ind_I2,1)-(solution_vector_DD(ind_V2,1)-solution_vector_DD(ind_Vout,1))/R1;

    % Vd
    Jaco_DD(ind_Vd,ind_Vd)=1; Jaco_DD(ind_Vd,ind_Vout)=-1;
    res_DD(ind_Vd,1)=0;
    
    % Vs
    Jaco_DD(ind_Vs,ind_Vs)=1;
    res_DD(ind_Vs,1)=solution_vector_DD(ind_Vs,1)-0;

    % Vg
    Jaco_DD(ind_Vg,ind_Vg)=1;
    res_DD(ind_Vg,1)=solution_vector_DD(ind_Vg,1)-V_gate;

    % V1
    Jaco_DD(ind_V1,ind_V1)=1; Jaco_DD(ind_V1,ind_Vout)=-1;
    res_DD(ind_V1,1)=0;
    
    % V2
    Jaco_DD(ind_V2,ind_V2)=1;
    res_DD(ind_V2,1)=solution_vector_DD(ind_V2,1)-Vdd;

    % Vout
    Jaco_DD(ind_Vout,ind_Id)=1; Jaco_DD(ind_Vout,ind_I1)=1;
    res_DD(ind_Vout,1)=0;
 

% Scaling %
    Table_Jaco_mixed_tmp=zeros(11,3);
    Table_Jaco_mixed_tmp(1:5,3)=4; Table_Jaco_mixed_tmp(6:11,3)=5;
    Table_Jaco_mixed=[Table_Jaco; Table_Jaco_mixed_tmp];
    Cvector_DD=zeros(size(solution_vector_DD,1),1); v=[V_T; max(abs(N_dop)); max(abs(N_dop)); 1; 1];
    for ii=1:size(solution_vector_DD,1)
        Cvector_DD(ii,1)=v(Table_Jaco_mixed(ii,3),1);
    end
    Jaco_DD=sparse(Jaco_DD);
    Cmatrix_DD=spdiags(Cvector_DD,0,size(solution_vector_DD,1),size(solution_vector_DD,1));
    Jaco_scaled_DD=Jaco_DD*Cmatrix_DD;
    Rvector_DD=1./sum(abs(Jaco_scaled_DD),2);
    Rmatrix_DD=spdiags(Rvector_DD,0,size(solution_vector_DD,1),size(solution_vector_DD,1));
    Jaco_scaled_DD=Rmatrix_DD*Jaco_scaled_DD;
    res_scaled_DD=Rmatrix_DD*res_DD;
    update_scaled_DD=Jaco_scaled_DD\(-res_scaled_DD);
    update_DD(:,Newton_DD)=Cmatrix_DD*update_scaled_DD;
    
    solution_vector_DD(:,1)=solution_vector_DD(:,1)+update_DD(:,Newton_DD);
    solution_vector_DD_saved(:,Newton_DD+1)=solution_vector_DD(:,1);
    
    % break
    update_DD_phi_for_break=update_DD(find(Table_Jaco(:,1)==2, 1 , 'first'):3:find(Table_Jaco(:,1)==4, 1 , 'last'),Newton_DD);
    fprintf("Max Updatevector=%d\n\n" , max(abs(update_DD_phi_for_break)))
    if max(abs(update_DD_phi_for_break))<1e-15
        break;
    end

%     % break
%     update_DD_phi_for_break=update_DD(size(solution_vector_poisson,1)+1,Newton_DD);
%     fprintf("Max Updatevector=%d\n\n" , max(abs(update_DD_phi_for_break)))
%     if max(abs(update_DD_phi_for_break))<1e-10
%         break;
%     end

end
Result_name(:,1)=["Id"; "Is"; "Ig"; "I1"; "I2"; "Vd"; "Vs"; "Vg"; "V1"; "V2"; "Vout"];
Result_DD(:,1)=solution_vector_DD(size(solution_vector_poisson,1)+1:size(solution_vector_DD,1));

%% Save
save('./data/homework25_dd.mat')

%% Visualize
% Visual_solution_vector_DD=zeros(size(Vertex,1),3);
% for ii=1:size(Table_Jaco,1)
%     Visual_solution_vector_DD(Table_Jaco(ii,2),Table_Jaco(ii,3))=solution_vector_DD_saved(ii,Newton_DD);
% end
% 
% figure(1) % 점으로 3D
% plot3(Vertex(:,1),Vertex(:,2),Visual_solution_vector_DD(:,1),'*')
% xlim([0 50*1e-9])
% 
% figure(2); % mesh 모양
% patch('Faces',Element,'Vertices',Vertex(:,[1 2]),'EdgeColor','black','FaceColor','none','LineWidth',1);
% title('Mesh')
% xlim([0 50*1e-9])
% xlabel('X');
% ylabel('Y');
% hold on
% patch('Faces',edge_interface,'Vertices',Vertex(:,[1 2]),'EdgeColor','green','FaceColor','none','LineWidth',3);
% hold on
% patch('Faces',Element_si_source,'Vertices',Vertex(:,[1 2]), 'EdgeColor','black','FaceColor','yellow','LineWidth',0.5)
% hold on
% patch('Faces',Element_si_drain,'Vertices',Vertex(:,[1 2]), 'EdgeColor','black','FaceColor','yellow','LineWidth',0.5)
% hold on
% patch('Faces',Element_si_channel,'Vertices',Vertex(:,[1 2]), 'EdgeColor','black','FaceColor','red','LineWidth',0.5)
% hold on
% patch('Faces',Element_ox1,'Vertices',Vertex(:,[1 2]), 'EdgeColor','black','FaceColor','blue','LineWidth',0.5)
% hold on
% patch('Faces',Element_ox2,'Vertices',Vertex(:,[1 2]), 'EdgeColor','black','FaceColor','blue','LineWidth',0.5)
% hold on
% patch('Faces',Element_ox2,'Vertices',Vertex(:,[1 2]), 'EdgeColor','black','FaceColor','blue','LineWidth',0.5)
% hold on
% patch('Faces',Element_source,'Vertices',Vertex(:,[1 2]), 'EdgeColor','m','FaceColor','none','LineWidth',5)
% hold on
% patch('Faces',Element_drain,'Vertices',Vertex(:,[1 2]), 'EdgeColor','m','FaceColor','none','LineWidth',5)
% hold on
% patch('Faces',Element_Dirichlet_BC,'Vertices',Vertex(:,[1 2]), 'EdgeColor','g','FaceColor','none','LineWidth',5)
% hold off
% 
% figure(3); % potential
% patch('Faces',Element,'Vertices',Vertex(:,[1 2]), 'FaceVertexCData',Visual_solution_vector_DD(:,1), 'EdgeColor','black','FaceColor','interp');
% title('Initial potential')
% hold on
% patch('Faces',edge_interface,'Vertices',Vertex(:,[1 2]),'EdgeColor','black','FaceColor','none','LineWidth',5);
% hold on
% patch('Faces',Element_source,'Vertices',Vertex(:,[1 2]), 'EdgeColor','m','FaceColor','none','LineWidth',5)
% hold on
% patch('Faces',Element_drain,'Vertices',Vertex(:,[1 2]), 'EdgeColor','m','FaceColor','none','LineWidth',5)
% hold on
% patch('Faces',Element_Dirichlet_BC,'Vertices',Vertex(:,[1 2]), 'EdgeColor','g','FaceColor','none','LineWidth',5)
% hold off
% xlim([0 50*1e-9])
% xlabel('X');
% ylabel('Y');
% colorbar
% 
% figure(4); % electron
% patch('Faces',Element,'Vertices',Vertex(:,[1 2]), 'FaceVertexCData',Visual_solution_vector_DD(:,2), 'EdgeColor','black','FaceColor','interp');
% title('elctron')
% xlim([0 50*1e-9])
% xlabel('X');
% ylabel('Y');
% colorbar
% 
% figure(5); % hole
% patch('Faces',Element,'Vertices',Vertex(:,[1 2]), 'FaceVertexCData',Visual_solution_vector_DD(:,3), 'EdgeColor','black','FaceColor','interp');
% title('hole')
% xlim([0 50*1e-9])
% xlabel('X');
% ylabel('Y');
% colorbar