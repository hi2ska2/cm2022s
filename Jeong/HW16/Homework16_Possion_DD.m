clear; close all; clc;
mkdir('./data')

%% Variable %%
eps=11.7; % Epsilon=[eps1; eps2; eps3]
q=1.602192e-19; %Elementary charge,C
eps0=8.854187817e-12; %Vacuum permittivity
k=1.3806488e-23; % Bolzmann constant, J/K
n_int=1.075e16; % 1.0e10 /cm^-3
T=300; % Temperture, K
V_T=k*T/q; % Thermal voltage
V_gate_workfunction=0; % Voltage, si_ox workfunction.
% Na=1e25; % /m^3
Nd=2e23; % /m^3
N_dop=Nd; % -1e18 /cm^3, 행을 region별 도핑 농도
nm=10^-9; %길이단위
mobility_n = 518e-4;  % electron mobility
mobility_p = 470.5e-4;     % Hole mobility

coeff=q/eps0;
coeff_Jn=q*mobility_n*V_T;
coeff_Jp=-q*mobility_p*V_T;


%% Build a Structure %%
% Vertex 값 작성
Vertex=vertex_doublegate(0,0, 1200,100, 50, 50);    % vertex = vertex_doublegate(x1, y1, x2, y2, dx, dy), [nm]

% anode_BC 불러오기
anode_BC_tmp=Contact(0,0,0,100,V_T*log(Nd/n_int));
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


cathode_BC_tmp=Contact(1200,0, 1200,100, V_T*log(Nd/n_int));
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
    phi(ii,1)=V_T*log(Nd/n_int);
end

% anode_BC
for ii=1:size(anode_BC,1)
    for j=1:size(Table_initial,1)
        if Table_initial(j,1)==1 && Table_initial(j,2)==anode_BC(ii,1)
            phi(j,1)=anode_BC(ii,2);
        end
    end
end

% cathode_BC
for ii=1:size(cathode_BC,1)
    for j=1:size(Table_initial,1)
        if Table_initial(j,1)==1 && Table_initial(j,2)==cathode_BC(ii,1)
            phi(j,1)=cathode_BC(ii,2);
        end
    end
end

% source, channel, drain 영역에서 initial carrrier density를 구해줌.
n_0=zeros(size(Table_initial,1),1); p_0=zeros(size(Table_initial,1),1);
for ii=1:size(Table_initial,1)
    n_0(ii,1)=n_int*exp(phi(ii,1)/V_T); % /m^3
    p_0(ii,1)=n_int*exp(-phi(ii,1)/V_T); % /m^3
end

% Initial value
initial_value=[phi n_0 p_0];

n=1;
for ii=1:size(Table_initial,1)
    for j=1:3
        solution_vector_possion(n,1)=initial_value(ii,j);
        n=n+1;
    end
end

%% Visualize
Visual_solution_vector_possion=zeros(size(Vertex,1),3);
for ii=1:size(Table_Jaco,1)
    Visual_solution_vector_possion(Table_Jaco(ii,2),Table_Jaco(ii,3))=solution_vector_possion(ii,1);
end

figure(1) % 점으로 3D
plot3(Vertex(:,1),Vertex(:,2),Visual_solution_vector_possion(:,1),'*')

%% Newton Method
solution_vector_possion_saved(:,1)=solution_vector_possion(:,1);
coeff=q/eps0;     % m^2 -> nm^2으로 단위 변경, q/eps0 해줌

% Jacobian matrix / res vector
Iteration_possion=40;
for Newton_possion=1:Iteration_possion
    fprintf("Newton Method, Newton_possion=%d\n" , Newton_possion)
%     Jaco_possion=sparse(zeros(size(solution_vector_possion,1),size(solution_vector_possion,1))); res_possion=zeros(size(solution_vector_possion,1),1);
    Jaco_possion=zeros(size(solution_vector_possion,1),size(solution_vector_possion,1)); res_possion=zeros(size(solution_vector_possion,1),1); %% 희소행렬 아님. Jaco matrix 확인용.
    
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

        Jaco_tmp_possion=zeros(9,9);
        Jaco_tmp_possion(1,1)=eps(Table_element_region(ii,1),1)*(-edge(ii,1)/L(ii,1)-edge(ii,3)/L(ii,3));
        Jaco_tmp_possion(1,2)=-coeff*Control_Volume(1,1);
        Jaco_tmp_possion(1,3)=coeff*Control_Volume(1,1);
        Jaco_tmp_possion(1,4)=eps(Table_element_region(ii,1),1)*edge(ii,1)/L(ii,1);
        Jaco_tmp_possion(1,7)=eps(Table_element_region(ii,1),1)*edge(ii,3)/L(ii,3);

        Jaco_tmp_possion(2,1)=-1/V_T*n_int*exp(solution_vector_possion(index(1,1),1)/V_T);
        Jaco_tmp_possion(2,2)=1;

        Jaco_tmp_possion(3,1)=1/V_T*n_int*exp(-solution_vector_possion(index(1,1),1)/V_T);
        Jaco_tmp_possion(3,3)=1;

        Jaco_tmp_possion(4,1)=eps(Table_element_region(ii,1),1)*edge(ii,1)/L(ii,1);
        Jaco_tmp_possion(4,4)=eps(Table_element_region(ii,1),1)*(-edge(ii,1)/L(ii,1)-edge(ii,2)/L(ii,2));
        Jaco_tmp_possion(4,5)=-coeff*Control_Volume(2,1);
        Jaco_tmp_possion(4,6)=coeff*Control_Volume(2,1);
        Jaco_tmp_possion(4,7)=eps(Table_element_region(ii,1),1)*edge(ii,2)/L(ii,2);

        Jaco_tmp_possion(5,4)=-1/V_T*n_int*exp(solution_vector_possion(index(4,1),1)/V_T);
        Jaco_tmp_possion(5,5)=1;

        Jaco_tmp_possion(6,4)=1/V_T*n_int*exp(-solution_vector_possion(index(4,1),1)/V_T);
        Jaco_tmp_possion(6,6)=1;

        Jaco_tmp_possion(7,1)=eps(Table_element_region(ii,1),1)*edge(ii,3)/L(ii,3);
        Jaco_tmp_possion(7,4)=eps(Table_element_region(ii,1),1)*edge(ii,2)/L(ii,2);
        Jaco_tmp_possion(7,7)=eps(Table_element_region(ii,1),1)*(-edge(ii,2)/L(ii,2)-edge(ii,3)/L(ii,3));
        Jaco_tmp_possion(7,8)=-coeff*Control_Volume(3,1);
        Jaco_tmp_possion(7,9)=coeff*Control_Volume(3,1);

        Jaco_tmp_possion(8,7)=-1/V_T*n_int*exp(solution_vector_possion(index(7,1),1)/V_T);
        Jaco_tmp_possion(8,8)=1;

        Jaco_tmp_possion(9,7)=1/V_T*n_int*exp(-solution_vector_possion(index(7,1),1)/V_T);
        Jaco_tmp_possion(9,9)=1;

        for j=1:9
            for k=1:9
                Jaco_possion(index(j,1),index(k,1))=Jaco_possion(index(j,1),index(k,1))+Jaco_tmp_possion(j,k);
            end
        end
        % n,p 값은 변하지 않는다.
        for j=1:3
            Jaco_possion(index(3*j-1),:)=0;
            Jaco_possion(index(3*j),:)=0;
            Jaco_possion(index(3*j-1),index(3*j-2))=-1/V_T*n_int*exp(solution_vector_possion(index(3*j-2,1),1)/V_T);
            Jaco_possion(index(3*j-1),index(3*j-1))=1;
            Jaco_possion(index(3*j),index(3*j-2))=1/V_T*n_int*exp(-solution_vector_possion(index(3*j-2,1),1)/V_T);
            Jaco_possion(index(3*j),index(3*j))=1;
        end

        % res vector
        res_tmp_possion=zeros(9,1);
        res_potential_tmp1(1,1)=eps(Table_element_region(ii,1),1)*(edge(ii,1)/L(ii,1)*(solution_vector_possion(index(4,1),1)-solution_vector_possion(index(1,1),1))+edge(ii,3)/L(ii,3)*(solution_vector_possion(index(7,1),1)-solution_vector_possion(index(1,1),1)));
        res_potential_tmp1(2,1)=eps(Table_element_region(ii,1),1)*(edge(ii,2)/L(ii,2)*(solution_vector_possion(index(7,1),1)-solution_vector_possion(index(4,1),1))+edge(ii,1)/L(ii,1)*(solution_vector_possion(index(1,1),1)-solution_vector_possion(index(4,1),1)));
        res_potential_tmp1(3,1)=eps(Table_element_region(ii,1),1)*(edge(ii,3)/L(ii,3)*(solution_vector_possion(index(1,1),1)-solution_vector_possion(index(7,1),1))+edge(ii,2)/L(ii,2)*(solution_vector_possion(index(4,1),1)-solution_vector_possion(index(7,1),1)));
        for j=1:3
            res_potential_tmp2(j,1)=(N_dop(Table_element_region(ii,1),1)-solution_vector_possion(index(3*j-1,1),1)+solution_vector_possion(index(3*j,1),1))*coeff*Control_Volume(j,1);
        end
        res_potential=res_potential_tmp1+res_potential_tmp2;
        for j=1:3
            res_electron(j,1)=solution_vector_possion(index(3*j-1,1),1)-n_int*exp(solution_vector_possion(index(3*j-2,1),1)/V_T);
        end
        for j=1:3
            res_hole(j,1)=solution_vector_possion(index(3*j,1),1)-n_int*exp(-solution_vector_possion(index(3*j-2,1),1)/V_T);
        end
        n=1;
        for j=1:3
            res_tmp_possion(n,1)=res_potential(j,1);
            res_tmp_possion(n+1,1)=res_electron(j,1);
            res_tmp_possion(n+2,1)=res_hole(j,1);
            n=n+3;
        end

        for j=1:9
            res_possion(index(j,1),1)=res_possion(index(j,1),1)+res_tmp_possion(j,1);
        end
        %n,p값은 변하지 않는다.
        for j=1:3
            res_possion(index(3*j-1,1),1)=res_electron(j,1);
            res_possion(index(3*j,1),1)=res_hole(j,1);
        end

    end

    % cathode_BC
    for ii=1:size(anode_BC,1)
        BC_index=find(Table_Jaco(:,1)==1 & Table_Jaco(:,2)==anode_BC(ii,1) & Table_Jaco(:,3)==1);
        Jaco_possion(BC_index,:)=0; Jaco_possion(BC_index+1,:)=0; Jaco_possion(BC_index+2,:)=0;
        Jaco_possion(BC_index,BC_index)=1; Jaco_possion(BC_index+1,BC_index+1)=1; Jaco_possion(BC_index+2,BC_index+2)=1;
        res_possion(BC_index,1)=solution_vector_possion(BC_index,1)-V_T*log(Nd/n_int);
        res_possion(BC_index+1,1)=solution_vector_possion(BC_index+1,1)-Nd;
        res_possion(BC_index+2,1)=solution_vector_possion(BC_index+2,1)-n_int^2/Nd;
    end

    % cathode_BC
    for ii=1:size(cathode_BC,1)
        BC_index=find(Table_Jaco(:,1)==1 & Table_Jaco(:,2)==cathode_BC(ii,1) & Table_Jaco(:,3)==1);
        Jaco_possion(BC_index,:)=0; Jaco_possion(BC_index+1,:)=0; Jaco_possion(BC_index+2,:)=0;
        Jaco_possion(BC_index,BC_index)=1; Jaco_possion(BC_index+1,BC_index+1)=1; Jaco_possion(BC_index+2,BC_index+2)=1;
        res_possion(BC_index,1)=solution_vector_possion(BC_index,1)-V_T*log(Nd/n_int);
        res_possion(BC_index+1,1)=solution_vector_possion(BC_index+1,1)-Nd;
        res_possion(BC_index+2,1)=solution_vector_possion(BC_index+2,1)-n_int^2/Nd;
    end

    % Scaling
    Cvector_possion=zeros(size(solution_vector_possion,1),1); v=[V_T; max(abs(N_dop)); max(abs(N_dop))];
    for ii=1:size(solution_vector_possion,1)
        Cvector_possion(ii,1)=v(Table_Jaco(ii,3),1);
    end
    Jaco_possion=sparse(Jaco_possion);
    Cmatrix_possion=spdiags(Cvector_possion,0,size(solution_vector_possion,1),size(solution_vector_possion,1));
    Jaco_scaled_possion=Jaco_possion*Cmatrix_possion;
    Rvector_possion=1./sum(abs(Jaco_scaled_possion),2);
    Rmatrix_possion=spdiags(Rvector_possion,0,size(solution_vector_possion,1),size(solution_vector_possion,1));
    Jaco_scaled_possion=Rmatrix_possion*Jaco_scaled_possion;
    res_scaled_possion=Rmatrix_possion*res_possion;
    update_scaled_possion=Jaco_scaled_possion\(-res_scaled_possion);
    update_possion(:,Newton_possion)=Cmatrix_possion*update_scaled_possion;
    
    solution_vector_possion(:,1)=solution_vector_possion(:,1)+update_possion(:,Newton_possion);
    solution_vector_possion_saved(:,Newton_possion+1)=solution_vector_possion(:,1);

    % break
    update_possion_phi_for_break=update_possion(find(Table_Jaco(:,1)==1, 1 , 'first'):3:find(Table_Jaco(:,1)==1, 1 , 'last'),Newton_possion);
    fprintf("Max Updatevector=%d\n\n" , max(abs(update_possion_phi_for_break)))
    if max(abs(update_possion_phi_for_break))<1e-15
        break;
    end
end

    update_possion_phi=update_possion(find(Table_Jaco(:,1)==1, 1 , 'first'):3:find(Table_Jaco(:,1)==1, 1 , 'last'),:);

%% Save

save('./data/Homework16_possion.mat')

%% Visualize
Visual_solution_vector_possion=zeros(size(Vertex,1),3);
for ii=1:size(Table_Jaco,1)
    Visual_solution_vector_possion(Table_Jaco(ii,2),Table_Jaco(ii,3))=solution_vector_possion_saved(ii,Newton_possion);
end

figure(1) % 점으로 3D
plot3(Vertex(:,1),Vertex(:,2),Visual_solution_vector_possion(:,1),'*')
%% Drift-Diffusion Model

Iteration_DD=40;

% Jacobian matrix / res vector
solution_vector_DD(:,1)=solution_vector_possion(:,1);
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
        n1=solution_vector_DD(index(1,1)+1,1); n2=solution_vector_DD(index(4,1)+1,1); n3=solution_vector_DD(index(7,1)+1,1);
        p1=solution_vector_DD(index(1,1)+2,1); p2=solution_vector_DD(index(4,1)+2,1); p3=solution_vector_DD(index(7,1)+2,1);

        % Jaco %
        Jaco_tmp_DD=zeros(9,9);
        Jaco_tmp_DD(1,1)=eps(Table_element_region(ii,1),1)*(-edge(ii,1)/L(ii,1)-edge(ii,3)/L(ii,3));
        Jaco_tmp_DD(1,2)=-coeff*Control_Volume(1,1);
        Jaco_tmp_DD(1,3)=coeff*Control_Volume(1,1);
        Jaco_tmp_DD(1,4)=eps(Table_element_region(ii,1),1)*edge(ii,1)/L(ii,1);
        Jaco_tmp_DD(1,7)=eps(Table_element_region(ii,1),1)*edge(ii,3)/L(ii,3);

        Jaco_tmp_DD(2,1)=coeff_Jn*((edge(ii,1)/L(ii,1))*(n2/(-V_T)*Ber_d(x21/V_T)-n1/(V_T)*Ber_d(-x21/V_T))+((edge(ii,3)/L(ii,3))*(n3/(-V_T)*Ber_d(x31/V_T)-n1/(V_T)*Ber_d(-x31/V_T))));
        Jaco_tmp_DD(2,2)=coeff_Jn*((edge(ii,1)/L(ii,1))*(-Ber(-x21/V_T))+((edge(ii,3)/L(ii,3))*(-Ber(-x31/V_T))));
        Jaco_tmp_DD(2,4)=coeff_Jn*((edge(ii,1)/L(ii,1))*(n2/(V_T)*Ber_d(x21/V_T)-n1/(-V_T)*Ber_d(-x21/V_T)));
        Jaco_tmp_DD(2,5)=coeff_Jn*((edge(ii,1)/L(ii,1))*(Ber(x21/V_T)));
        Jaco_tmp_DD(2,7)=coeff_Jn*((edge(ii,3)/L(ii,3))*(n3/(V_T)*Ber_d(x31/V_T)-n1/(-V_T)*Ber_d(-x31/V_T)));
        Jaco_tmp_DD(2,8)=coeff_Jn*((edge(ii,3)/L(ii,3))*(Ber(x31/V_T)));

        Jaco_tmp_DD(3,1)=coeff_Jp*((edge(ii,1)/L(ii,1))*(p2/(V_T)*Ber_d(-x21/V_T)-p1/(-V_T)*Ber_d(x21/V_T))+((edge(ii,3)/L(ii,3))*(p3/(V_T)*Ber_d(-x31/V_T)-p1/(-V_T)*Ber_d(x31/V_T))));
        Jaco_tmp_DD(3,3)=coeff_Jp*((edge(ii,1)/L(ii,1))*(-Ber(x21/V_T))+(edge(ii,3)/L(ii,3))*(-Ber(x31/V_T)));
        Jaco_tmp_DD(3,4)=coeff_Jp*((edge(ii,1)/L(ii,1))*(p2/(-V_T)*Ber_d(-x21/V_T)-p1/(V_T)*Ber_d(x21/V_T)));
        Jaco_tmp_DD(3,6)=coeff_Jp*((edge(ii,1)/L(ii,1))*(Ber(-x21/V_T)));
        Jaco_tmp_DD(3,7)=coeff_Jp*((edge(ii,3)/L(ii,3))*(p3/(-V_T)*Ber_d(-x31/V_T)-p1/(V_T)*Ber_d(x31/V_T)));
        Jaco_tmp_DD(3,9)=coeff_Jp*((edge(ii,3)/L(ii,3))*(Ber(-x31/V_T)));

        Jaco_tmp_DD(4,1)=eps(Table_element_region(ii,1),1)*edge(ii,1)/L(ii,1);
        Jaco_tmp_DD(4,4)=eps(Table_element_region(ii,1),1)*(-edge(ii,1)/L(ii,1)-edge(ii,2)/L(ii,2));
        Jaco_tmp_DD(4,5)=-coeff*Control_Volume(2,1);
        Jaco_tmp_DD(4,6)=coeff*Control_Volume(2,1);
        Jaco_tmp_DD(4,7)=eps(Table_element_region(ii,1),1)*edge(ii,2)/L(ii,2);

        Jaco_tmp_DD(5,1)=coeff_Jn*((edge(ii,1)/L(ii,1))*(n1/(V_T)*Ber_d(x12/V_T)-n2/(-V_T)*Ber_d(-x12/V_T)));
        Jaco_tmp_DD(5,2)=coeff_Jn*((edge(ii,1)/L(ii,1))*(Ber(x12/V_T)));
        Jaco_tmp_DD(5,4)=coeff_Jn*((edge(ii,2)/L(ii,2))*(n3/(-V_T)*Ber_d(x32/V_T)-n2/(V_T)*Ber_d(-x32/V_T))+((edge(ii,1)/L(ii,1))*(n1/(-V_T)*Ber_d(x12/V_T)-n2/(V_T)*Ber_d(-x12/V_T))));
        Jaco_tmp_DD(5,5)=coeff_Jn*((edge(ii,2)/L(ii,2))*(-Ber(-x32/V_T))+((edge(ii,1)/L(ii,1))*(-Ber(-x12/V_T))));
        Jaco_tmp_DD(5,7)=coeff_Jn*((edge(ii,2)/L(ii,2))*(n3/(V_T)*Ber_d(x32/V_T)-n2/(-V_T)*Ber_d(-x32/V_T)));
        Jaco_tmp_DD(5,8)=coeff_Jn*((edge(ii,2)/L(ii,2))*(Ber(x32/V_T)));

        Jaco_tmp_DD(6,1)=coeff_Jp*((edge(ii,1)/L(ii,1))*(p1/(-V_T)*Ber_d(-x12/V_T)-p2/(V_T)*Ber_d(x12/V_T)));
        Jaco_tmp_DD(6,3)=coeff_Jp*((edge(ii,1)/L(ii,1))*(Ber(-x12/V_T)));
        Jaco_tmp_DD(6,4)=coeff_Jp*((edge(ii,2)/L(ii,2))*(p3/(V_T)*Ber_d(-x32/V_T)-p2/(-V_T)*Ber_d(x32/V_T))+((edge(ii,1)/L(ii,1))*(p1/(V_T)*Ber_d(-x12/V_T)-p2/(-V_T)*Ber_d(x12/V_T))));
        Jaco_tmp_DD(6,6)=coeff_Jp*((edge(ii,2)/L(ii,2))*(-Ber(x32/V_T))+(edge(ii,1)/L(ii,1))*(-Ber(x12/V_T)));
        Jaco_tmp_DD(6,7)=coeff_Jp*((edge(ii,2)/L(ii,2))*(p3/(-V_T)*Ber_d(-x32/V_T)-p2/(V_T)*Ber_d(x32/V_T)));
        Jaco_tmp_DD(6,9)=coeff_Jp*((edge(ii,2)/L(ii,2))*(Ber(-x32/V_T)));

        Jaco_tmp_DD(7,1)=eps(Table_element_region(ii,1),1)*edge(ii,3)/L(ii,3);
        Jaco_tmp_DD(7,4)=eps(Table_element_region(ii,1),1)*edge(ii,2)/L(ii,2);
        Jaco_tmp_DD(7,7)=eps(Table_element_region(ii,1),1)*(-edge(ii,2)/L(ii,2)-edge(ii,3)/L(ii,3));
        Jaco_tmp_DD(7,8)=-coeff*Control_Volume(3,1);
        Jaco_tmp_DD(7,9)=coeff*Control_Volume(3,1);

        Jaco_tmp_DD(8,1)=coeff_Jn*((edge(ii,3)/L(ii,3))*(n1/(V_T)*Ber_d(x13/V_T)-n3/(-V_T)*Ber_d(-x13/V_T)));
        Jaco_tmp_DD(8,2)=coeff_Jn*((edge(ii,3)/L(ii,3))*(Ber(x13/V_T)));
        Jaco_tmp_DD(8,4)=coeff_Jn*((edge(ii,2)/L(ii,2))*(n2/(V_T)*Ber_d(x23/V_T)-n3/(-V_T)*Ber_d(-x23/V_T)));
        Jaco_tmp_DD(8,5)=coeff_Jn*((edge(ii,2)/L(ii,2))*(Ber(x23/V_T)));
        Jaco_tmp_DD(8,7)=coeff_Jn*((edge(ii,3)/L(ii,3))*(n1/(-V_T)*Ber_d(x13/V_T)-n3/(V_T)*Ber_d(-x13/V_T))+((edge(ii,2)/L(ii,2))*(n2/(-V_T)*Ber_d(x23/V_T)-n3/(V_T)*Ber_d(-x23/V_T))));
        Jaco_tmp_DD(8,8)=coeff_Jn*((edge(ii,3)/L(ii,3))*(-Ber(-x13/V_T))+((edge(ii,2)/L(ii,2))*(-Ber(-x23/V_T))));

        Jaco_tmp_DD(9,1)=coeff_Jp*((edge(ii,3)/L(ii,3))*(p1/(-V_T)*Ber_d(-x13/V_T)-p3/(V_T)*Ber_d(x13/V_T)));
        Jaco_tmp_DD(9,3)=coeff_Jp*((edge(ii,3)/L(ii,3))*(Ber(-x13/V_T)));
        Jaco_tmp_DD(9,4)=coeff_Jp*((edge(ii,2)/L(ii,2))*(p2/(-V_T)*Ber_d(-x23/V_T)-p3/(V_T)*Ber_d(x23/V_T)));
        Jaco_tmp_DD(9,6)=coeff_Jp*((edge(ii,2)/L(ii,2))*(Ber(-x23/V_T)));
        Jaco_tmp_DD(9,7)=coeff_Jp*((edge(ii,3)/L(ii,3))*(p1/(V_T)*Ber_d(-x13/V_T)-p3/(-V_T)*Ber_d(x13/V_T))+((edge(ii,2)/L(ii,2))*(p2/(V_T)*Ber_d(-x23/V_T)-p3/(-V_T)*Ber_d(x23/V_T))));
        Jaco_tmp_DD(9,9)=coeff_Jp*((edge(ii,3)/L(ii,3))*(-Ber(x13/V_T))+(edge(ii,2)/L(ii,2))*(-Ber(x23/V_T)));


        for j=1:9
            for k=1:9
                Jaco_DD(index(j,1),index(k,1))=Jaco_DD(index(j,1),index(k,1))+Jaco_tmp_DD(j,k);
            end
        end

        % res vector
        res_tmp_DD=zeros(9,1);
        res_potential_tmp1_DD(1,1)=eps(Table_element_region(ii,1),1)*(edge(ii,1)/L(ii,1)*(solution_vector_DD(index(4,1),1)-solution_vector_DD(index(1,1),1))+edge(ii,3)/L(ii,3)*(solution_vector_DD(index(7,1),1)-solution_vector_DD(index(1,1),1)));
        res_potential_tmp1_DD(2,1)=eps(Table_element_region(ii,1),1)*(edge(ii,2)/L(ii,2)*(solution_vector_DD(index(7,1),1)-solution_vector_DD(index(4,1),1))+edge(ii,1)/L(ii,1)*(solution_vector_DD(index(1,1),1)-solution_vector_DD(index(4,1),1)));
        res_potential_tmp1_DD(3,1)=eps(Table_element_region(ii,1),1)*(edge(ii,3)/L(ii,3)*(solution_vector_DD(index(1,1),1)-solution_vector_DD(index(7,1),1))+edge(ii,2)/L(ii,2)*(solution_vector_DD(index(4,1),1)-solution_vector_DD(index(7,1),1)));
        for j=1:3
            res_potential_tmp2_DD(j,1)=(N_dop(Table_element_region(ii,1),1)-solution_vector_DD(index(3*j-1,1),1)+solution_vector_DD(index(3*j,1),1))*coeff*Control_Volume(j,1);
        end

        res_potential_DD=res_potential_tmp1_DD+res_potential_tmp2_DD;

        % Jn
        res_Jn=zeros(3,1);
        res_Jn(1,1)=coeff_Jn*((edge(ii,1)/L(ii,1))*(n2*Ber(x21/V_T)-n1*Ber(-x21/V_T))+(edge(ii,3)/L(ii,3))*(n3*Ber(x31/V_T)-n1*Ber(-x31/V_T)));
        res_Jn(2,1)=coeff_Jn*((edge(ii,2)/L(ii,2))*(n3*Ber(x32/V_T)-n2*Ber(-x32/V_T))+(edge(ii,1)/L(ii,1))*(n1*Ber(x12/V_T)-n2*Ber(-x12/V_T)));
        res_Jn(3,1)=coeff_Jn*((edge(ii,3)/L(ii,3))*(n1*Ber(x13/V_T)-n3*Ber(-x13/V_T))+(edge(ii,2)/L(ii,2))*(n2*Ber(x23/V_T)-n3*Ber(-x23/V_T)));

        % Jp
        res_Jp=zeros(3,1);
        res_Jp(1,1)=coeff_Jp*((edge(ii,1)/L(ii,1))*(p2*Ber(-x21/V_T)-p1*Ber(x21/V_T))+(edge(ii,3)/L(ii,3))*(p3*Ber(-x31/V_T)-p1*Ber(x31/V_T)));
        res_Jp(2,1)=coeff_Jp*((edge(ii,2)/L(ii,2))*(p3*Ber(-x32/V_T)-p2*Ber(x32/V_T))+(edge(ii,1)/L(ii,1))*(p1*Ber(-x12/V_T)-p2*Ber(x12/V_T)));
        res_Jp(3,1)=coeff_Jp*((edge(ii,3)/L(ii,3))*(p1*Ber(-x13/V_T)-p3*Ber(x13/V_T))+(edge(ii,2)/L(ii,2))*(p2*Ber(-x23/V_T)-p3*Ber(x23/V_T)));

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
    end

% Boundary condition %
    % cathode_BC
    for ii=1:size(anode_BC,1)
        BC_index=find(Table_Jaco(:,1)==1 & Table_Jaco(:,2)==anode_BC(ii,1) & Table_Jaco(:,3)==1);
        Jaco_DD(BC_index,:)=0; Jaco_DD(BC_index+1,:)=0; Jaco_DD(BC_index+2,:)=0;
        Jaco_DD(BC_index,BC_index)=1; Jaco_DD(BC_index+1,BC_index+1)=1; Jaco_DD(BC_index+2,BC_index+2)=1;
        res_DD(BC_index,1)=solution_vector_DD(BC_index,1)-V_T*log(Nd/n_int);
        res_DD(BC_index+1,1)=solution_vector_DD(BC_index+1,1)-Nd;
        res_DD(BC_index+2,1)=solution_vector_DD(BC_index+2,1)-n_int^2/Nd;
    end

    % cathode_BC
    for ii=1:size(cathode_BC,1)
        BC_index=find(Table_Jaco(:,1)==1 & Table_Jaco(:,2)==cathode_BC(ii,1) & Table_Jaco(:,3)==1);
        Jaco_DD(BC_index,:)=0; Jaco_DD(BC_index+1,:)=0; Jaco_DD(BC_index+2,:)=0;
        Jaco_DD(BC_index,BC_index)=1; Jaco_DD(BC_index+1,BC_index+1)=1; Jaco_DD(BC_index+2,BC_index+2)=1;
        res_DD(BC_index,1)=solution_vector_DD(BC_index,1)-V_T*log(Nd/n_int);
        res_DD(BC_index+1,1)=solution_vector_DD(BC_index+1,1)-Nd;
        res_DD(BC_index+2,1)=solution_vector_DD(BC_index+2,1)-n_int^2/Nd;
    end

% Scaling %
    Cvector_DD=zeros(size(solution_vector_DD,1),1); v=[V_T; max(abs(N_dop)); max(abs(N_dop))];
    for ii=1:size(solution_vector_DD,1)
        Cvector_DD(ii,1)=v(Table_Jaco(ii,3),1);
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
    update_DD_phi_for_break=update_DD(find(Table_Jaco(:,1)==1, 1 , 'first'):3:find(Table_Jaco(:,1)==1, 1 , 'last'),Newton_DD);
    fprintf("Max Updatevector=%d\n\n" , max(abs(update_DD_phi_for_break)))
    if max(abs(update_DD_phi_for_break))<1e-15
        break;
    end

end

update_DD_phi=update_DD(find(Table_Jaco(:,1)==1, 1 , 'first'):3:find(Table_Jaco(:,1)==1, 1 , 'last'),:);

%% Save
save('./data/Homework16_dd.mat')

%% Visualize
Visual_solution_vector_DD=zeros(size(Vertex,1),3);
for ii=1:size(Table_Jaco,1)
    Visual_solution_vector_DD(Table_Jaco(ii,2),Table_Jaco(ii,3))=solution_vector_DD_saved(ii,Newton_DD);
end

figure(1) % 점으로 3D
plot3(Vertex(:,1),Vertex(:,2),Visual_solution_vector_DD(:,1),'*')