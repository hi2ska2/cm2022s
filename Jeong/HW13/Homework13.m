clear; close all; clc;

%%%%%%%%%%%%%%%% Variable %%%%%%%%%%%%%%%%%
eps=[3.9; 11.7; 11.7; 11.7; 3.9]; % Epsilon=[eps1; eps2; eps3]
q=1.602192e-19; %Elementary charge,C
eps0=8.854187817e-12; %Vacuum permittivity
k=1.3806488e-23; % Bolzmann constant, J/K
n_int=1e16; % 1.0e10 /cm^-3
T=300; % Temperture, K
V_T=k*T/q; % Thermal voltage
V_gate=0.33374; % Voltage
Na=2e21; % /m^3
Nd=5e23; % /m^3
N_dop=[0; 1e24; -1e24; 1e24; 0]; % -1e18 /cm^3, 행을 region별 도핑 농도
nm=10^-9; %길이단위
mobility_n = 1417e-4;  % electron mobility
mobility_p = 470.5e-4;     % Hole mobility
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 사전작업
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

Source_datin_BC_tmp = importdata("Source_datin_BC.txt"); % Dirichlet_BC = [index potential]
m=0; n=0; % element로 이루어진 경계조건 계산하기 위해 변환
for ii=1:size(Source_datin_BC_tmp,1)
    for j=size(Source_datin_BC_tmp,2):-1:1
        if ~isnan(Source_datin_BC_tmp(ii,j))
            for k=1:j-1
                Source_datin_BC(k+n,:)=[Source_datin_BC_tmp(ii,k) Source_datin_BC_tmp(ii,j)];
                m=m+1;
            end
            n=m;
            break;
        end
    end
end

% Vertex 값 작성
Vertex=vertex_doublegate(0,0, 10,7, 1, 1);    % vertex = vertex_doublegate(x1, y1, x2, y2, dx, dy), [nm]
% vertex=1e-9.*Vertex;

% element 값 입력
Element_ox1 = element(0,0, 10,2);                 % Element = element(x1, y1, x2, y2), [nm]
Element_si_source = element(0,2, 2,5);
Element_si_channel = element(2,2, 8,5);
Element_si_drain = element(8,2, 10,5);
Element_ox2 = element(0,5, 10,7);
Element = [Element_ox1; Element_si_source; Element_si_channel; Element_si_drain; Element_ox2];

% element region table
Table_element_region1=zeros(length(Element_ox1),1)+1;
Table_element_region2=zeros(length(Element_si_source),1)+2;
Table_element_region3=zeros(length(Element_si_channel),1)+3;
Table_element_region4=zeros(length(Element_si_drain),1)+4;
Table_element_region5=zeros(length(Element_ox2),1)+5;
Table_element_region=[Table_element_region1; Table_element_region2; Table_element_region3; Table_element_region4; Table_element_region5];


%%%% 각 Region의 edge 정렬 %%%%
n=1;
for ii=1:length(Element_ox1) % ox region
    for j=1:3
        edge_1(n,:)= [Element_ox1(ii,2) Element_ox1(ii,1)];
        edge_1(n+1,:)= [Element_ox1(ii,3) Element_ox1(ii,2)];
        edge_1(n+2,:)= [Element_ox1(ii,1) Element_ox1(ii,3)];
        n=n+3;
    end
end
n=1;
for ii=1:length(Element_si_source) % si region
    for j=1:3
        edge_2(n,:)= [Element_si_source(ii,1) Element_si_source(ii,2)];
        edge_2(n+1,:)= [Element_si_source(ii,2) Element_si_source(ii,3)];
        edge_2(n+2,:)= [Element_si_source(ii,3) Element_si_source(ii,1)];
        n=n+3;
    end
end
n=1;
for ii=1:length(Element_si_channel) % si region
    for j=1:3
        edge_3(n,:)= [Element_si_channel(ii,1) Element_si_channel(ii,2)];
        edge_3(n+1,:)= [Element_si_channel(ii,2) Element_si_channel(ii,3)];
        edge_3(n+2,:)= [Element_si_channel(ii,3) Element_si_channel(ii,1)];
        n=n+3;
    end
end
n=1;
for ii=1:length(Element_si_drain) % si region
    for j=1:3
        edge_4(n,:)= [Element_si_drain(ii,1) Element_si_drain(ii,2)];
        edge_4(n+1,:)= [Element_si_drain(ii,2) Element_si_drain(ii,3)];
        edge_4(n+2,:)= [Element_si_drain(ii,3) Element_si_drain(ii,1)];
        n=n+3;
    end
end
n=1;
for ii=1:length(Element_ox2) % ox region
    for j=1:3
        edge_5(n,:)= [Element_ox2(ii,2) Element_ox2(ii,1)];
        edge_5(n+1,:)= [Element_ox2(ii,3) Element_ox2(ii,2)];
        edge_5(n+2,:)= [Element_ox2(ii,1) Element_ox2(ii,3)];
        n=n+3;
    end
end

% 중복 edge 정보 제거
edge_1=sort(edge_1, 2); edge_1=unique(edge_1,'rows');
edge_2=sort(edge_2, 2); edge_2=unique(edge_2,'rows');
edge_3=sort(edge_3, 2); edge_3=unique(edge_3,'rows');
edge_4=sort(edge_4, 2); edge_4=unique(edge_4,'rows');
edge_5=sort(edge_5, 2); edge_5=unique(edge_5,'rows');

% 각 region에서 동일한 edge가 있는지 확인, Interface
edge_interface1_2 = intersect(edge_1,edge_2,'rows');
edge_interface1_3 = intersect(edge_1,edge_3,'rows');
edge_interface1_4 = intersect(edge_1,edge_4,'rows');
edge_interface1_5 = intersect(edge_1,edge_5,'rows');
edge_interface2_3 = intersect(edge_2,edge_3,'rows');
edge_interface2_4 = intersect(edge_2,edge_4,'rows');
edge_interface2_5 = intersect(edge_2,edge_5,'rows');
edge_interface3_4 = intersect(edge_3,edge_4,'rows');
edge_interface3_5 = intersect(edge_3,edge_5,'rows');
edge_interface4_5 = intersect(edge_4,edge_5,'rows');
edge_interface=[edge_interface1_2; edge_interface1_3; edge_interface1_4; edge_interface1_5; edge_interface2_3; edge_interface2_4; edge_interface2_5; edge_interface3_4; edge_interface3_5; edge_interface4_5;];
Vertex_interface1_2 = unique(sort(reshape(edge_interface1_2,[],1)));
Vertex_interface1_3 = unique(sort(reshape(edge_interface1_3,[],1)));
Vertex_interface1_4 = unique(sort(reshape(edge_interface1_4,[],1)));
Vertex_interface1_5 = unique(sort(reshape(edge_interface1_5,[],1)));
Vertex_interface2_3 = unique(sort(reshape(edge_interface2_3,[],1)));
Vertex_interface2_4 = unique(sort(reshape(edge_interface2_4,[],1)));
Vertex_interface2_5 = unique(sort(reshape(edge_interface2_5,[],1)));
Vertex_interface3_4 = unique(sort(reshape(edge_interface3_4,[],1)));
Vertex_interface3_5 = unique(sort(reshape(edge_interface3_5,[],1)));
Vertex_interface4_5 = unique(sort(reshape(edge_interface4_5,[],1)));

% 각 region을 내 vertex 갯수 확인, Number of vertex
Vertex_region1 = unique(sort(reshape(Element_ox1,[],1)));
Vertex_region2 = unique(sort(reshape(Element_si_source,[],1)));
Vertex_region3 = unique(sort(reshape(Element_si_channel,[],1)));
Vertex_region4 = unique(sort(reshape(Element_si_drain,[],1)));
Vertex_region5 = unique(sort(reshape(Element_ox2,[],1)));
Number_of_vertex=[length(Vertex_region1); length(Vertex_region2); length(Vertex_region3); length(Vertex_region4); length(Vertex_region5)];

% Region 별 vertex
Vertex_region_tmp=[Vertex_region1; Vertex_region2; Vertex_region3; Vertex_region4; Vertex_region5];
Vertex_region=zeros(max(Number_of_vertex),length(Number_of_vertex));
n=1;
for j=1:size(Number_of_vertex,1)
    m=1;
    for ii=n:Number_of_vertex(j,1)+(n-1)
        Vertex_region(m,j)=Vertex_region_tmp(ii,1);
        m=m+1;
    end
    n=n+Number_of_vertex(j,1);
end

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
edge=zeros(length(Element),3);
for ii=1:length(Element)
    edge_tmp=[0.5*L(ii,1) 0.5*L(ii,2) 0.5*L(ii,3)];
    for j=1:3
        edge(ii,j) = real(sqrt(R(ii,1)^2-edge_tmp(:,j)^2));
%         if edge(ii,j)<=1e-17
%             edge(ii,j)=0;
%         end
    end
end

% AA table 작성
n=1;
for nRegion=1:length(Number_of_vertex)
    for nVertex=1:Number_of_vertex(nRegion,1)
        Table_AA(n,:)=[nRegion Vertex_region(nVertex,nRegion)];
        n=n+1;
    end
end

% A matrix / b vector 구하기
AA=zeros(sum(Number_of_vertex),sum(Number_of_vertex)); b=zeros(sum(Number_of_vertex),1); AA_tmp=zeros(sum(Number_of_vertex),sum(Number_of_vertex));
for ii=1:length(Element)
    % AA matrix
    AA_tmp=[(-edge(ii,1)/L(ii,1)-edge(ii,3)/L(ii,3)) edge(ii,1)/L(ii,1) edge(ii,3)/L(ii,3);
        edge(ii,1)/L(ii,1) -edge(ii,1)/L(ii,1)-edge(ii,2)/L(ii,2) edge(ii,2)/L(ii,2);
        edge(ii,3)/L(ii,3) edge(ii,2)/L(ii,2) -edge(ii,2)/L(ii,2)-edge(ii,3)/L(ii,3)];
    AA_tmp=eps(Table_element_region(ii,1),1)*AA_tmp;
    for j=1:3
        for k=1:3
            p=find(Table_AA(:,1)==Table_element_region(ii,1) & Table_AA(:,2)==Element(ii,j)); % table 활용 index 찾기
            r=find(Table_AA(:,1)==Table_element_region(ii,1) & Table_AA(:,2)==Element(ii,k)); % table 활용 index 찾기
            AA(p,r)=AA(p,r)+AA_tmp(j,k);
        end
    end

    % b vector
    Control_Volume=[(0.25)*L(ii,3)*edge(ii,3)+(0.25)*L(ii,1)*edge(ii,1);
        (0.25)*L(ii,1)*edge(ii,1)+(0.25)*L(ii,2)*edge(ii,2);
        (0.25)*L(ii,2)*edge(ii,2)+(0.25)*L(ii,3)*edge(ii,3)];
    Control_Volume=1e-18.*Control_Volume;
    for j=1:3
        p=find(Table_AA(:,1)==Table_element_region(ii,1) & Table_AA(:,2)==Element(ii,j)); % table 활용 index 찾기
        b(p,1)=b(p,1)+(-q/eps0)*N_dop(Table_element_region(ii,1),1)*Control_Volume(j,1);
    end
end

%%% Interface_BC %%%
% ox 1의 값을 si_source에 더해준 후 ox값은 1,-1로 변경하기
for ii=1:length(Vertex_interface1_2)
    ox_index=find(Table_AA(:,1)==1 & Table_AA(:,2)==Vertex_interface1_2(ii,1));
    si_index=find(Table_AA(:,1)==2 & Table_AA(:,2)==Vertex_interface1_2(ii,1));
    AA(si_index,:)=AA(si_index,:)+AA(ox_index,:);
    AA(ox_index,:)=0; AA(ox_index,ox_index)=1; AA(ox_index,si_index)=-1;
end
% ox 1의 값을 si_draion에 더해주기.
for ii=1:length(Vertex_interface1_4)
    ox_index=find(Table_AA(:,1)==1 & Table_AA(:,2)==Vertex_interface1_4(ii,1));
    si_index=find(Table_AA(:,1)==4 & Table_AA(:,2)==Vertex_interface1_4(ii,1));
    AA(si_index,:)=AA(si_index,:)+AA(ox_index,:);
    AA(ox_index,:)=0; AA(ox_index,ox_index)=1; AA(ox_index,si_index)=-1;
end
% ox 2의 값을 si_source에 더해주기.
for ii=1:length(Vertex_interface2_5)
    ox_index=find(Table_AA(:,1)==5 & Table_AA(:,2)==Vertex_interface2_5(ii,1));
    si_index=find(Table_AA(:,1)==2 & Table_AA(:,2)==Vertex_interface2_5(ii,1));
    AA(si_index,:)=AA(si_index,:)+AA(ox_index,:);
    AA(ox_index,:)=0; AA(ox_index,ox_index)=1; AA(ox_index,si_index)=-1;
end
% ox 2의 값을 si_draion에 더해주기.
for ii=1:length(Vertex_interface4_5)
    ox_index=find(Table_AA(:,1)==5 & Table_AA(:,2)==Vertex_interface4_5(ii,1));
    si_index=find(Table_AA(:,1)==4 & Table_AA(:,2)==Vertex_interface4_5(ii,1));
    AA(si_index,:)=AA(si_index,:)+AA(ox_index,:);
    AA(ox_index,:)=0; AA(ox_index,ox_index)=1; AA(ox_index,si_index)=-1;
end
% si_source의 값을 si_channel에 더해주기.
for ii=1:length(Vertex_interface2_3)
    source_index=find(Table_AA(:,1)==2 & Table_AA(:,2)==Vertex_interface2_3(ii,1));
    channel_index=find(Table_AA(:,1)==3 & Table_AA(:,2)==Vertex_interface2_3(ii,1));
    AA(channel_index,:)=AA(channel_index,:)+AA(source_index,:);
    AA(source_index,:)=0; AA(source_index,source_index)=1; AA(source_index,channel_index)=-1;
end
% si_drain의 값을 si_channel에 더해주기.
for ii=1:length(Vertex_interface3_4)
    drain_index=find(Table_AA(:,1)==4 & Table_AA(:,2)==Vertex_interface3_4(ii,1));
    channel_index=find(Table_AA(:,1)==3 & Table_AA(:,2)==Vertex_interface3_4(ii,1));
    AA(channel_index,:)=AA(channel_index,:)+AA(drain_index,:);
    AA(drain_index,:)=0; AA(drain_index,drain_index)=1; AA(drain_index,channel_index)=-1;
end

% Dirichlet_BC
for ii=1:size(Dirichlet_BC,1)
    BC_index=find(Table_AA(:,2)==Dirichlet_BC(ii,1));
    AA(BC_index,:)=0;
    AA(BC_index,BC_index)=1;
end

% Source_datin_BC
for ii=1:size(Source_datin_BC,1)
    BC_index=find(Table_AA(:,2)==Source_datin_BC(ii,1));
    AA(BC_index,:)=0;
    AA(BC_index,BC_index)=1;
end

%b matrix_BC
for ii=1:size(Dirichlet_BC,1)
    BC_index=find(Table_AA(:,2)==Dirichlet_BC(ii,1));
    b(BC_index,1)=Dirichlet_BC(ii,2)+V_gate;
end

%b matrix_BC(Source_datin_BC)
for ii=1:size(Source_datin_BC,1)
    BC_index=find(Table_AA(:,2)==Source_datin_BC(ii,1));
    b(BC_index,1)=Source_datin_BC(ii,2);
end

phi=AA\b;

%Electron/Hole concentration
n_0=zeros(length(phi),1); p_0=zeros(length(phi),1);
for ii=1:length(phi)
    if Table_AA(ii,1)==2 || Table_AA(ii,1)==3 || Table_AA(ii,1)==4
        n_0(ii,1)=n_int*exp(phi(ii,1)/V_T); % /m^3
        p_0(ii,1)=n_int*exp(-phi(ii,1)/V_T); % /m^3
    end
end

% Initial value
initial_value_0=[phi n_0 p_0];
initial_value=[phi n_0 p_0];
for ii=1:length(Vertex_interface1_2)
    n=find(Table_AA(:,1)==1 & Table_AA(:,2)==Vertex_interface1_2(ii,1));
    initial_value(n,:)=nan;
end
for ii=1:length(Vertex_interface1_3)
    n=find(Table_AA(:,1)==1 & Table_AA(:,2)==Vertex_interface1_3(ii,1));
    initial_value(n,:)=nan;
end
for ii=1:length(Vertex_interface1_4)
    n=find(Table_AA(:,1)==1 & Table_AA(:,2)==Vertex_interface1_4(ii,1));
    initial_value(n,:)=nan;
end
for ii=1:length(Vertex_interface2_5)
    n=find(Table_AA(:,1)==5 & Table_AA(:,2)==Vertex_interface2_5(ii,1));
    initial_value(n,:)=nan;
end
for ii=1:length(Vertex_interface3_5)
    n=find(Table_AA(:,1)==5 & Table_AA(:,2)==Vertex_interface3_5(ii,1));
    initial_value(n,:)=nan;
end
for ii=1:length(Vertex_interface4_5)
    n=find(Table_AA(:,1)==5 & Table_AA(:,2)==Vertex_interface4_5(ii,1));
    initial_value(n,:)=nan;
end
for ii=1:length(Vertex_interface2_3)
    n=find(Table_AA(:,1)==2 & Table_AA(:,2)==Vertex_interface2_3(ii,1));
    initial_value(n,:)=nan;
end
for ii=1:length(Vertex_interface3_4)
    n=find(Table_AA(:,1)==4 & Table_AA(:,2)==Vertex_interface3_4(ii,1));
    initial_value(n,:)=nan;
end
initial_value = rmmissing(initial_value);

% % Viusalize
figure(1);
patch('Faces',Element,'Vertices',Vertex(:,[1 2]), 'FaceVertexCData',initial_value(:,1), 'EdgeColor','black','FaceColor','interp');
title('Initial potential')
xlabel('X');
ylabel('Y');
colorbar

figure(2);
patch('Faces',Element,'Vertices',Vertex(:,[1 2]), 'FaceVertexCData',initial_value(:,2), 'EdgeColor','black','FaceColor','interp');
title('elctron')
xlabel('X');
ylabel('Y');
colorbar

figure(3);
patch('Faces',Element,'Vertices',Vertex(:,[1 2]), 'FaceVertexCData',initial_value(:,3), 'EdgeColor','black','FaceColor','interp');
title('hole')
xlabel('X');
ylabel('Y');
colorbar

%% initial value 구하기/ n형:V_T*ln(Nd/n_int), p형:-V_T*ln(Nd/n_int) 사용
phi=zeros(size(Vertex,1),1);
for ii=[23:25 34:36 45:47 56:58 31:33 42:44 53:55 64:66]
    phi(ii,1)=V_T*log(Nd/n_int);
end
for ii=[26:30 37:41 46:52 59:63]
    phi(ii,1)=-V_T*log(Na/n_int);
end

% Dirichlet_BC
for ii=1:size(Dirichlet_BC,1)
    phi(Dirichlet_BC(ii,1),1)=Dirichlet_BC(ii,2)+V_gate;
end
% Source_datin_BC
for ii=1:size(Source_datin_BC,1)
    phi(Source_datin_BC(ii,1),1)=Source_datin_BC(ii,2);
end

n_0=zeros(size(Vertex,1),1); p_0=zeros(size(Vertex,1),1);
for ii=23:66
    n_0(ii,1)=n_int*exp(phi(ii,1)/V_T); % /m^3
    p_0(ii,1)=n_int*exp(-phi(ii,1)/V_T); % /m^3
end

% Initial value
initial_value=[phi n_0 p_0];


%% Newton Method
% Jacobian Table 생성
n=1;
for nRegion=1:length(Number_of_vertex)
    for nVertex=1:Number_of_vertex(nRegion,1)
        if nRegion==2
            for variable=1:3
                Table_Jaco(n,:)=[nRegion Vertex_region(nVertex,nRegion) variable];
                n=n+1;
            end
        else
            Table_Jaco(n,:)=[nRegion Vertex_region(nVertex,nRegion) 1];
            n=n+1;
        end
    end
end

% Solution vector x 생성
n=1;
for ii=1:length(Table_Jaco)
    solution_vector(n,1)=initial_value(Table_Jaco(ii,2),Table_Jaco(ii,3));
    n=n+1;
end

% Jacobian matrix / res vector

for Newton=1:10
    Jaco=zeros(size(solution_vector,1),size(solution_vector,1)); res=zeros(size(solution_vector,1),1);
    % Potential
    for ii=1:length(Element)
        % Control_Volume 계산
        Control_Volume=[(0.25)*L(ii,3)*edge(ii,3)+(0.25)*L(ii,1)*edge(ii,1);
            (0.25)*L(ii,1)*edge(ii,1)+(0.25)*L(ii,2)*edge(ii,2);
            (0.25)*L(ii,2)*edge(ii,2)+(0.25)*L(ii,3)*edge(ii,3)];
        Control_Volume=1e-18.*Control_Volume; % nm^2 -> m^2으로 단위 변경

        if Table_element_region(ii,1)==2
            % index : 현재 element에 속해있는 vertex가 resdue에 몇번째 있는가 저장
            index=zeros(9,1);
            n=1;
            for j=1:3
                index(n,1)=find(Table_Jaco(:,1)==Table_element_region(ii,1) & Table_Jaco(:,2)==Element(ii,j) & Table_Jaco(:,3)==1);
                index(n+1,1)=index(n,1)+1;
                index(n+2,1)=index(n,1)+2;
                n=n+3;
            end
            %schafetter용 상수
            for j=1:3
            coeff(j,1)=(q*mobility_n*V_T)*(Area(ii,index(j,1))/L(ii,index(j,1)));
            end
            % charge 계산
            for j=1:3
                charge(j,1)=(q/eps0)*(N_dop(Table_element_region(ii,1),1)-solution_vector(index(3*j-1,1),Newton)+solution_vector(index(3*j,1),Newton));
            end
            Jaco_tmp=zeros(9,9);
            % Jaco_potential matrix
            Jaco_tmp=[eps(Table_element_region(ii,1),1)*(-edge(ii,1)/L(ii,1)-edge(ii,3)/L(ii,3)) -q*Control_Volume(1,1) q*Control_Volume(1,1) eps(Table_element_region(ii,1),1)*edge(ii,1)/L(ii,1) 0 0 eps(Table_element_region(ii,1),1)*edge(ii,3)/L(ii,3) 0 0;
                -1/V_T*n_int*exp(solution_vector(index(1,1),Newton)/V_T) 1 0 0 0 0 0 0 0;
                1/V_T*n_int*exp(-solution_vector(index(1,1),Newton)/V_T) 0 1 0 0 0 0 0 0;
                eps(Table_element_region(ii,1),1)*edge(ii,1)/L(ii,1) 0 0 eps(Table_element_region(ii,1),1)*(-edge(ii,1)/L(ii,1)-edge(ii,2)/L(ii,2)) -q*Control_Volume(2,1) q*Control_Volume(2,1) eps(Table_element_region(ii,1),1)*edge(ii,2)/L(ii,2) 0 0;
                0 0 0 -1/V_T*n_int*exp(solution_vector(index(4,1),Newton)/V_T) 1 0 0 0 0;
                0 0 0 1/V_T*n_int*exp(-solution_vector(index(4,1),Newton)/V_T) 0 1 0 0 0;
                eps(Table_element_region(ii,1),1)*edge(ii,3)/L(ii,3) 0 0 eps(Table_element_region(ii,1),1)*edge(ii,2)/L(ii,2) 0 0 eps(Table_element_region(ii,1),1)*(-edge(ii,2)/L(ii,2)-edge(ii,3)/L(ii,3)) -q*Control_Volume(3,1) q*Control_Volume(3,1);
                0 0 0 0 0 0 -1/V_T*n_int*exp(solution_vector(index(7,1),Newton)/V_T) 1 0;
                0 0 0 0 0 0 1/V_T*n_int*exp(-solution_vector(index(7,1),Newton)/V_T) 0 1];
            for j=1:9
                for k=1:9
                    Jaco(index(j,1),index(k,1))=Jaco(index(j,1),index(k,1))+Jaco_tmp(j,k);
                end
            end
            % n,p 값은 변하지 않는다.
            for j=1:3
                Jaco(index(j*3-1),:)=0;
                Jaco(index(j*3),:)=0;
                Jaco(index(j*3-1),index(j*3-2))=-1/V_T*n_int*exp(solution_vector(index(3*j-2,1),Newton)/V_T);
                Jaco(index(j*3-1),index(j*3-1))=1;
                Jaco(index(j*3),index(j*3-2))=1/V_T*n_int*exp(-solution_vector(index(3*j-2,1),Newton)/V_T);
                Jaco(index(j*3),index(j*3))=1;
            end

            % res vector
            res_tmp=zeros(9,1);
            res_potential_tmp1=eps(Table_element_region(ii,1),1).*[(-edge(ii,1)/L(ii,1)-edge(ii,3)/L(ii,3))*solution_vector(index(1,1),Newton)+edge(ii,1)/L(ii,1)*solution_vector(index(4,1),Newton)+edge(ii,3)/L(ii,3)*solution_vector(index(7,1),Newton);
                                                                     edge(ii,1)/L(ii,1)*solution_vector(index(1,1),Newton)+(-edge(ii,1)/L(ii,1)-edge(ii,2)/L(ii,2))*solution_vector(index(4,1),Newton)+edge(ii,2)/L(ii,2)*solution_vector(index(7,1),Newton);
                                                                     edge(ii,3)/L(ii,3)*solution_vector(index(1,1),Newton)+edge(ii,2)/L(ii,2)*solution_vector(index(4,1),Newton)+(-edge(ii,2)/L(ii,2)-edge(ii,3)/L(ii,3))*solution_vector(index(7,1),Newton)];
            for j=1:3
                res_potential_tmp2(j,1)=charge(j,1)*Control_Volume(j,1);
            end
            res_potential=res_potential_tmp1+res_potential_tmp2;
            for j=1:3
                res_electron(j,1)=solution_vector(index(3*j-1,1),Newton)-n_int*exp(solution_vector(index(3*j-2,1),Newton)/V_T);
            end
            for j=1:3
                res_hole(j,1)=solution_vector(index(3*j,1),Newton)-n_int*exp(-solution_vector(index(3*j-2,1),Newton)/V_T);
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

        elseif Table_element_region(ii,1)==1 || Table_element_region(ii,1)==3
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
            res_tmp=eps(Table_element_region(ii,1),1).*[(-edge(ii,1)/L(ii,1)-edge(ii,3)/L(ii,3))*solution_vector(index(1,1),Newton)+edge(ii,1)/L(ii,1)*solution_vector(index(2,1),Newton)+edge(ii,3)/L(ii,3)*solution_vector(index(3,1),Newton);
                edge(ii,1)/L(ii,1)*solution_vector(index(1,1),Newton)+(-edge(ii,1)/L(ii,1)-edge(ii,2)/L(ii,2))*solution_vector(index(2,1),Newton)+edge(ii,2)/L(ii,2)*solution_vector(index(3,1),Newton);
                edge(ii,3)/L(ii,3)*solution_vector(index(1,1),Newton)+edge(ii,2)/L(ii,2)*solution_vector(index(2,1),Newton)+(-edge(ii,2)/L(ii,2)-edge(ii,3)/L(ii,3))*solution_vector(index(3,1),Newton)];
            for j=1:3
                res(index(j,1),1)=res(index(j,1),1)+res_tmp(j,1);
            end
        end
    end

    %%% Interface_BC %%%
    % ox 1의 값을 si에 더해준 후 ox값은 1,-1로 변경하기
    for ii=1:length(Vertex_interface1_2)
        ox_index=find(Table_Jaco(:,1)==1 & Table_Jaco(:,2)==Vertex_interface1_2(ii,1) & Table_Jaco(:,3)==1);
        si_index=find(Table_Jaco(:,1)==2 & Table_Jaco(:,2)==Vertex_interface1_2(ii,1) & Table_Jaco(:,3)==1);
        Jaco(si_index,:)=Jaco(si_index,:)+Jaco(ox_index,:);
        Jaco(ox_index,:)=0; Jaco(ox_index,ox_index)=1; Jaco(ox_index,si_index)=-1;
        res(si_index,1)=res(si_index,1)+res(ox_index,1);
        res(ox_index,1)=0;
    end
    % ox 2의 값을 si에 더해주기.
    for ii=1:length(Vertex_interface3)
        ox_index=find(Table_Jaco(:,1)==3 & Table_Jaco(:,2)==Vertex_interface3(ii,1) & Table_Jaco(:,3)==1);
        si_index=find(Table_Jaco(:,1)==2 & Table_Jaco(:,2)==Vertex_interface3(ii,1) & Table_Jaco(:,3)==1);
        Jaco(si_index,:)=Jaco(si_index,:)+Jaco(ox_index,:);
        Jaco(ox_index,:)=0; Jaco(ox_index,ox_index)=1; Jaco(ox_index,si_index)=-1;
        res(si_index,1)=res(si_index,1)+res(ox_index,1);
        res(ox_index,1)=0;
    end

    % Dirichlet_BC
    for ii=1:size(Dirichlet_BC,1)
        BC_index=find(Table_Jaco(:,2)==Dirichlet_BC(ii,1) & Table_Jaco(:,3)==1);
        Jaco(BC_index,:)=0;
        Jaco(BC_index,BC_index)=1;
    end

    %res matrix_BC
    for ii=1:size(Dirichlet_BC,1)
        BC_index=find(Table_Jaco(:,2)==Dirichlet_BC(ii,1) & Table_Jaco(:,3)==1);
        res(BC_index,1)=solution_vector(BC_index,Newton)-0.33374;
    end

    % Scaling
    Cvector=zeros(size(solution_vector,1),1); v=[V_T; max(abs(N_dop)); max(abs(N_dop))];
    for ii=1:size(solution_vector,1)
        Cvector(ii,1)=v(Table_Jaco(ii,3),1);
    end
    Cmatrix=spdiags(Cvector,0,size(solution_vector,1),size(solution_vector,1));
    Jaco_scaled=Jaco*Cmatrix;
    Rvector=1./sum(abs(Jaco_scaled),2);
    Rmatrix=spdiags(Rvector,0,size(solution_vector,1),size(solution_vector,1));
    Jaco_scaled=Rmatrix*Jaco_scaled;
    res_scaled=Rmatrix*res;
    update_scaled=Jaco_scaled\(-res_scaled);
    update(:,Newton)=Cmatrix*update_scaled;
    solution_vector(:,Newton+1)=solution_vector(:,Newton)+update(:,Newton);

    % % non-Scaling
    % update(:,Newton)=Jaco\-res;
    % solution_vector(:,Newton+1)=solution_vector(:,Newton)+update(:,Newton);

end

% % Visualize
% visualize_value=zeros(size(Vertex,1),3);
% for ii=1:size(Table_Jaco,1)
%     visualize_value(Table_Jaco(ii,2),Table_Jaco(ii,3))=solution_vector(ii,Newton);
% end
% 
% dif=abs(initial_value-visualize_value);
% 
% % figure(4);
% patch('Faces',Element,'Vertices',Vertex(:,[1 2]), 'FaceVertexCData',visualize_value(:,1), 'EdgeColor','black','FaceColor','interp');
% title('Initial potential')
% xlabel('X');
% ylabel('Y');
% colorbar
% 
% figure(5);
% patch('Faces',Element,'Vertices',Vertex(:,[1 2]), 'FaceVertexCData',visualize_value(:,2), 'EdgeColor','black','FaceColor','interp');
% title('elctron')
% xlabel('X');
% ylabel('Y');
% colorbar
% 
% figure(6);
% patch('Faces',Element,'Vertices',Vertex(:,[1 2]), 'FaceVertexCData',visualize_value(:,3), 'EdgeColor','black','FaceColor','interp');
% title('hole')
% xlabel('X');
% ylabel('Y');
% colorbar
% 
% figure(7);
% patch('Faces',Element,'Vertices',Vertex(:,[1 2]), 'FaceVertexCData',dif(:,1), 'EdgeColor','black','FaceColor','interp');
% title('Initial potential')
% xlabel('X');
% ylabel('Y');
% colorbar
% 
% figure(8);
% patch('Faces',Element,'Vertices',Vertex(:,[1 2]), 'FaceVertexCData',dif(:,2), 'EdgeColor','black','FaceColor','interp');
% title('elctron')
% xlabel('X');
% ylabel('Y');
% colorbar
% 
% figure(9);
% patch('Faces',Element,'Vertices',Vertex(:,[1 2]), 'FaceVertexCData',dif(:,3), 'EdgeColor','black','FaceColor','interp');
% title('hole')
% xlabel('X');
% ylabel('Y');
% colorbar
