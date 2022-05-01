clear; close all; clc;

% due to : 2022.05.3. (Tue) / HW15

%%%%% HW15 %%%%%%%
% Transient simulation

%%%% parameter %%%%

q=1.602e-19;
nint=1.075e+16;     %1/m^3
Na= 1e+25;          %1/m^3   P-type
Nd = 5e+26;         %1/m^3   N-type
k_B=1.380662e-23;
T=300;      %K
esi = 11.7; eox = 3.9; e0=8.854187817e-12;  % Permittivity, F/m
m_n = 1417e-4;      % electron mobility
m_p = 470.5e-4;     % Hole mobility
coeff=(0.25*q)/e0;

Thermal_V = k_B*T/q;

coeff_Jn = q*m_n*Thermal_V;
coeff_Jp = -q*m_p*Thermal_V;

%%%%%%%% Read Files %%%%%%%%%%

Vertex = importdata("Vertex.txt");
[V_row,V_col]=size(Vertex);

Element= importdata('Element.txt');
[E_row,col]=size(Element);

Contact = importdata('Contact.txt');

Contact_Si = importdata('Contact_Si.txt');

[C_row,C_col]=size(Contact);

index = V_row;

%%%%%%%% Specify Region 1,2 %%%%%%%%%%

E_R1 = importdata('Element_Region1.txt');

E_R2L = importdata('Element_Region2_L.txt');
E_R2M = importdata('Element_Region2_M.txt');
E_R2R = importdata('Element_Region2_R.txt');

E_R3 = importdata('Element_Region3.txt');

R_row2_L = size(E_R2L,1);
R_row2_M = size(E_R2M,1);
R_row2_R = size(E_R2R,1);

T_2 = R_row2_R+R_row2_M+R_row2_L;

for ii =1:T_2
    if ii <=R_row2_L
        E_R2(ii,:) = E_R2L(ii,:);
    elseif ii <= R_row2_M+R_row2_L
        E_R2(ii,:) = E_R2M(ii-R_row2_L,:);
    else
        E_R2(ii,:) = E_R2R(ii-R_row2_L-R_row2_M,:);

    end
end
R_row1 = size(E_R1,1);
R_row2 = size(E_R2,1);
R_row3 = size(E_R3,1);

%%%%%% Geometry information %%%%%%

L=zeros(E_row,3);

for ii = 1:R_row1
    for jj = 1:3
        if jj == 1 % between 1 and 2
            L(ii,jj) = sqrt(abs((Vertex(E_R1(ii,jj+1),1)-Vertex(E_R1(ii,jj),1))^2+(Vertex(E_R1(ii,jj+1),2)-Vertex(E_R1(ii,jj),2))^2));
        elseif jj == 2 % between 2 and 3
            L(ii,jj) = sqrt(abs((Vertex(E_R1(ii,jj+1),1)-Vertex(E_R1(ii,jj),1))^2+(Vertex(E_R1(ii,jj+1),2)-Vertex(E_R1(ii,jj),2))^2));
        elseif jj ==3 % between 3 and 1
            L(ii,jj) = sqrt(abs((Vertex(E_R1(ii,jj),1)-Vertex(E_R1(ii,jj-2),1))^2+(Vertex(E_R1(ii,jj),2)-Vertex(E_R1(ii,jj-2),2))^2));
        end
    end
end

for ii = 1:R_row2
    for jj = 1:3
        if jj == 1 % between 1 and 2
            L(ii+R_row1,jj) = sqrt(abs((Vertex(E_R2(ii,jj+1),1)-Vertex(E_R2(ii,jj),1))^2+(Vertex(E_R2(ii,jj+1),2)-Vertex(E_R2(ii,jj),2))^2));
        elseif jj == 2 % between 2 and 3
            L(ii+R_row1,jj) = sqrt(abs((Vertex(E_R2(ii,jj+1),1)-Vertex(E_R2(ii,jj),1))^2+(Vertex(E_R2(ii,jj+1),2)-Vertex(E_R2(ii,jj),2))^2));
        elseif jj ==3 % between 3 and 1
            L(ii+R_row1,jj) = sqrt(abs((Vertex(E_R2(ii,jj),1)-Vertex(E_R2(ii,jj-2),1))^2+(Vertex(E_R2(ii,jj),2)-Vertex(E_R2(ii,jj-2),2))^2));
        end
    end
end

for ii = 1:R_row3
    for jj = 1:3
        if jj == 1 % between 1 and 2
            L(ii+R_row2+R_row1,jj) = sqrt(abs((Vertex(E_R3(ii,jj+1),1)-Vertex(E_R3(ii,jj),1))^2+(Vertex(E_R3(ii,jj+1),2)-Vertex(E_R3(ii,jj),2))^2));
        elseif jj == 2 % between 2 and 3
            L(ii+R_row2+R_row1,jj) = sqrt(abs((Vertex(E_R3(ii,jj+1),1)-Vertex(E_R3(ii,jj),1))^2+(Vertex(E_R3(ii,jj+1),2)-Vertex(E_R3(ii,jj),2))^2));
        elseif jj ==3 % between 3 and 1
            L(ii+R_row2+R_row1,jj) = sqrt(abs((Vertex(E_R3(ii,jj),1)-Vertex(E_R3(ii,jj-2),1))^2+(Vertex(E_R3(ii,jj),2)-Vertex(E_R3(ii,jj-2),2))^2));
        end
    end
end

% Vector
Volume=zeros(E_row,1);
R=zeros(E_row,1);

for ii = 1:R_row1
    a12 = [(Vertex(E_R1(ii,2),1)-Vertex(E_R1(ii,1),1)),(Vertex(E_R1(ii,2),2)-Vertex(E_R1(ii,1),2))];
    a13 = [(Vertex(E_R1(ii,3),1)-Vertex(E_R1(ii,1),1)),(Vertex(E_R1(ii,3),2)-Vertex(E_R1(ii,1),2))];
    a23 = [(Vertex(E_R1(ii,3),1)-Vertex(E_R1(ii,2),1)),(Vertex(E_R1(ii,3),2)-Vertex(E_R1(ii,2),2))];
    a21=-a12;
    a31=-a13;
    a32=-a23;

    Volume(ii,1) = (1/2)*abs(a12(1,1)*a13(1,2)-a12(1,2)*a13(1,1)); % m^2
    R(ii,1) = (L(ii,1)*L(ii,2)*L(ii,3))/(4*Volume(ii,1));
end

for ii = 1:R_row2
    a12 = [(Vertex(E_R2(ii,2),1)-Vertex(E_R2(ii,1),1)),(Vertex(E_R2(ii,2),2)-Vertex(E_R2(ii,1),2))];
    a13 = [(Vertex(E_R2(ii,3),1)-Vertex(E_R2(ii,1),1)),(Vertex(E_R2(ii,3),2)-Vertex(E_R2(ii,1),2))];
    a23 = [(Vertex(E_R2(ii,3),1)-Vertex(E_R2(ii,2),1)),(Vertex(E_R2(ii,3),2)-Vertex(E_R2(ii,2),2))];
    a21=-a12;
    a31=-a13;
    a32=-a23;

    Volume(ii+R_row1,1) = (1/2)*abs(a12(1,1)*a13(1,2)-a12(1,2)*a13(1,1)); % m^2
    R(ii+R_row1,1) = (L(ii+R_row1,1)*L(ii+R_row1,2)*L(ii+R_row1,3))/(4*Volume(ii+R_row1,1));
end

for ii = 1:R_row3
    a12 = [(Vertex(E_R3(ii,2),1)-Vertex(E_R3(ii,1),1)),(Vertex(E_R3(ii,2),2)-Vertex(E_R3(ii,1),2))];
    a13 = [(Vertex(E_R3(ii,3),1)-Vertex(E_R3(ii,1),1)),(Vertex(E_R3(ii,3),2)-Vertex(E_R3(ii,1),2))];
    a23 = [(Vertex(E_R3(ii,3),1)-Vertex(E_R3(ii,2),1)),(Vertex(E_R3(ii,3),2)-Vertex(E_R3(ii,2),2))];
    a21=-a12;
    a31=-a13;
    a32=-a23;

    Volume(ii+R_row2+R_row1,1) = (1/2)*abs(a12(1,1)*a13(1,2)-a12(1,2)*a13(1,1)); % m^2
    R(ii+R_row2+R_row1,1) = (L(ii+R_row2+R_row1,1)*L(ii+R_row2+R_row1,2)*L(ii+R_row2+R_row1,3))/(4*Volume(ii+R_row2+R_row1,1));
end

%%%%%%% Area %%%%%%%%
Area=zeros(E_row,3);

for ii = 1:R_row1
    for jj = 1:3
        if jj == 1 % between 1 and 2
            Area(ii,jj) = sqrt(abs(R(ii,1)^2-(L(ii,jj)/2)^2));
        elseif jj == 2 % between 2 and 3
            Area(ii,jj) = sqrt(abs(R(ii,1)^2-(L(ii,jj)/2)^2));
        elseif jj ==3 % between 3 and 1
            Area(ii,jj) = sqrt(abs(R(ii,1)^2-(L(ii,jj)/2)^2));
        end
    end
end

for ii = 1:R_row2
    for jj = 1:3
        if jj == 1 % between 1 and 2
            Area(ii+R_row1,jj) = sqrt(abs(R(ii+R_row1,1)^2-(L(ii+R_row1,jj)/2)^2));
        elseif jj == 2 % between 2 and 3
            Area(ii+R_row1,jj) = sqrt(abs(R(ii+R_row1,1)^2-(L(ii+R_row1,jj)/2)^2));
        elseif jj ==3 % between 3 and 1
            Area(ii+R_row1,jj) = sqrt(abs(R(ii+R_row1,1)^2-(L(ii+R_row1,jj)/2)^2));
        end
    end
end

for ii = 1:R_row3
    for jj = 1:3
        if jj == 1 % between 1 and 2
            Area(ii+R_row2+R_row1,jj) = sqrt(abs(R(ii+R_row2+R_row1,1)^2-(L(ii+R_row2+R_row1,jj)/2)^2));
        elseif jj == 2 % between 2 and 3
            Area(ii+R_row2+R_row1,jj) = sqrt(abs(R(ii+R_row2+R_row1,1)^2-(L(ii+R_row2+R_row1,jj)/2)^2));
        elseif jj ==3 % between 3 and 1
            Area(ii+R_row2+R_row1,jj) = sqrt(abs(R(ii+R_row2+R_row1,1)^2-(L(ii+R_row2+R_row1,jj)/2)^2));
        end
    end
end

V_R1 = unique(E_R1);
V_R2 = unique(E_R2);
V_R3 = unique(E_R3);

VR1 = size(V_R1,1);
VR2 = size(V_R2,1);
VR3 = size(V_R3,1);

T = VR1+VR2+VR3;

%%%%%%%%% Find the interface Edge %%%%%%%%%%%

S1=intersect(E_R1, E_R2);
j =1;

for ii = 1:R_row2
    if size(intersect(E_R2(ii,:),S1'),2) == 2
        Edge1(j,:) = intersect(E_R2(ii,:),S1');
        j = j+1;
    end
end

S2=intersect(E_R3, E_R2);
j =1;

for ii = 1:R_row2
    if size(intersect(E_R2(ii,:),S2'),2) == 2
        Edge2(j,:) = intersect(E_R2(ii,:),S2);
        j = j+1;
    end
end

interface1 = unique(Edge1);
interface2 = unique(Edge2);

i_row1 = size(interface1,1);
i_row2 = size(interface2,1);

%%%%% Make the table %%%%%%%

Potential = 1;
eDensity = 2;
hDensity = 3;

Table = zeros(T,3);
T_row = size(Table,1);

for iRegion = 1:3
    if iRegion ==1  % Oxide , only potential

        Table(1:VR1,1) = 1;

        for iVertex =  1:VR1
            Table(iVertex,2) = V_R1(iVertex,1);
            Table(iVertex,3) = Potential;   % only Potential
        end

    elseif iRegion == 2     % Silicon, Potential, eDensity, hDensity

        Table(1+VR1:VR1+VR2,1) = 2;

        for iVertex =  1:VR2
            Table(VR1+iVertex,2) = V_R2(iVertex,1);
            Table(VR1+iVertex,3) = Potential;
        end

    elseif iRegion ==3  % Oxide , only potential

        Table(1+VR1+VR2:T,1) = 3;
        
        for iVertex =  1:VR3
            Table(VR1+VR2+iVertex,2) = V_R3(iVertex,1);
            Table(VR1+VR2+iVertex,3) = Potential;
        end
    end
end

%%%%%%%% Find the Silicon interface Edge %%%%%%%%%%%

S1=intersect(E_R2L, E_R2M);
j =1;

for ii = 1:R_row2_M
    if size(intersect(E_R2M(ii,:),S1'),2) == 2
        intery1(j,:) = intersect(E_R2M(ii,:),S1');
        j = j+1;
    end
end

S2=intersect(E_R2R, E_R2M);
j =1;

for ii = 1:R_row2_M
    if size(intersect(E_R2M(ii,:),S2'),2) == 2
        intery2(j,:) = intersect(E_R2M(ii,:),S2);
        j = j+1;
    end
end

iy11 = unique(intery1);
iy22 = unique(intery2);

ii1_row = size(iy11,1);

for ii = 1:ii1_row
    iy1 = find(Table(1:T,2)==iy11(ii,1));
    iy2 = find(Table(1:T,2)==iy22(ii,1));

    iy1_1(1:length(iy1),ii)=iy1;
    iy2_2(1:length(iy2),ii)=iy2;
end

% % region 2 interface vertex index number
iy1 = nonzeros(unique(iy1_1));
iy2 = nonzeros(unique(iy2_2));

%%%%%%%% Initial Poisson %%%%%%%%%%

phi = zeros(T_row,1);
elec = zeros(T_row,1);
hole = zeros(T_row,1);

phi(:,1) = Thermal_V*log(Nd/nint);

for ii = 1:size(iy1,1)

    phi(iy1(ii,1)+1:iy2(ii,1)-1,1) = Thermal_V*log((Na)/1e+4/nint);
    phi(iy2(ii,1)-3:iy2(ii,1)-1,1) = Thermal_V*log((Na)/10/nint);
    phi(iy1(ii,1)+1:iy1(ii,1)+3,1) = Thermal_V*log((Na)/10/nint);
    phi(iy2(ii,1)-5:iy2(ii,1)-3,1) = Thermal_V*log((Na)/5e+2/nint);
    phi(iy1(ii,1)+3:iy1(ii,1)+5,1) = Thermal_V*log((Na)/5e+2/nint);
end

ry1 = unique(Contact_Si(1:10,:)); % Source
ry2 = unique(Contact_Si(11:20,:)); % Drain

rr_row = size(ry1,1);
ii_row = size(iy1,1);

elec(1+VR1:1:VR2+VR1,1) = nint*exp(phi(1+VR1:1:VR2+VR1,1)/Thermal_V);
hole(1+VR1:1:VR2+VR1,1) = nint*exp(-phi(1+VR1:1:VR2+VR1,1)/Thermal_V);

%%%%%% Fully-Coupled %%%%%%%%

T1 = VR1+3*VR2+VR3;
Table = zeros(T1,3);

T_row = size(Table,1);

for iRegion = 1:3
    if iRegion ==1  % Oxide , only potential

        Table(1:VR1,1) = 1;

        for iVertex =  1:VR1
            Table(iVertex,2) = V_R1(iVertex,1);
            Table(iVertex,3) = Potential;   % only Potential
        end

    elseif iRegion == 2     % Silicon, Potential, eDensity, hDensity

        Table(1+VR1:VR1+3*VR2,1) = 2;

        for iVertex =  1:VR2

            Table(3*iVertex+(VR1-2),2) = V_R2(iVertex,1);
            Table(3*iVertex+(VR1-1),2) = V_R2(iVertex,1);
            Table(3*iVertex+(VR1),2) = V_R2(iVertex,1);

            Table(3*iVertex+(VR1-2),3) = Potential;
            Table(3*iVertex+(VR1-1),3) = eDensity;
            Table(3*iVertex+(VR1),3) = hDensity;
        end

    elseif iRegion ==3  % Oxide , only potential

        Table(1+VR1+3*VR2:T1,1) = 3;
        
        for iVertex =  1:VR3
            Table(VR1+3*VR2+iVertex,2) = V_R3(iVertex,1);
            Table(VR1+3*VR2+iVertex,3) = Potential;
        end
    end
end

%%%%%%%  Poisson, eDensity, hDensity %%%%%%%%%%

for Newton = 1:50

    Jaco = sparse(T_row,T_row);
    res = zeros(T_row,1);
    
    for ii=1:T_row
    
        if Table(ii,1) == 1     % Region1
            K = find(E_R1 == Table(ii,2));
            K_row = size(K,1);
    
            for rr = 1:K_row
    
                if Table(ii,3) == Potential
    
                    if K(rr,1) <= R_row1
    
                        jj = 1;
                        V3 = find(Table(1:VR1,2)==E_R1(K(rr,1),jj+2)) ;
                        V2 = find(Table(1:VR1,2)==E_R1(K(rr,1),jj+1));
                        V1 = find(Table(1:VR1,2)==E_R1(K(rr,1),jj));
                       
                        res(ii,1) = res(ii,1) + eox*(phi(V2,1)-phi(V1,1))*(Area(K(rr,1),jj)/L(K(rr,1),jj)) + eox*(phi(V3,1)-phi(V1,1))*(Area(K(rr,1),jj+2)/L(K(rr,1),jj+2));
                        Jaco(ii, V3) = Jaco(ii, V3)+ eox*(Area(K(rr,1),jj+2)/L(K(rr,1),jj+2));
                        Jaco(ii, V2) = Jaco(ii, V2)+eox*(Area(K(rr,1),jj)/L(K(rr,1),jj));
                        Jaco(ii, V1) = Jaco(ii, V1)-eox*(Area(K(rr,1),jj+2)/L(K(rr,1),jj+2))-eox*(Area(K(rr,1),jj)/L(K(rr,1),jj));
    
                    elseif K(rr,1) > R_row1 && K(rr,1) <= 2*R_row1
                        jj = 2;
                        K(rr,1) = K(rr,1)-R_row1;
    
                        V3 = find(Table(1:VR1,2)==E_R1(K(rr,1),jj+1)) ;
                        V2 = find(Table(1:VR1,2)==E_R1(K(rr,1),jj));
                        V1 = find(Table(1:VR1,2)==E_R1(K(rr,1),jj-1));
    
                        res(ii,1) = res(ii,1) + eox*(phi(V1,1)-phi(V2,1))*(Area(K(rr,1),jj-1)/L(K(rr,1),jj-1)) + eox*(phi(V3,1)-phi(V2,1))*(Area(K(rr,1),jj)/L(K(rr,1),jj));
                        Jaco(ii, V3) = Jaco(ii, V3)+eox*(Area(K(rr,1),jj)/L(K(rr,1),jj));
                        Jaco(ii, V2) =  Jaco(ii, V2)-eox*(Area(K(rr,1),jj)/L(K(rr,1),jj))-eox*(Area(K(rr,1),jj-1)/L(K(rr,1),jj-1));
                        Jaco(ii, V1) = Jaco(ii, V1)+eox*(Area(K(rr,1),jj-1)/L(K(rr,1),jj-1));
    
                    elseif K(rr,1) > 2*R_row1 && K(rr,1) <= 3*R_row1
                        jj = 3;
                        K(rr,1) = K(rr,1)-2*R_row1;
    
                        V3 = find(Table(1:VR1,2)==E_R1(K(rr,1),jj-2)) ;
                        V2 = find(Table(1:VR1,2)==E_R1(K(rr,1),jj-1));
                        V1 = find(Table(1:VR1,2)==E_R1(K(rr,1),jj));
    
                        res(ii,1) = res(ii,1) + eox*(phi(V3,1)-phi(V1,1))*(Area(K(rr,1),jj)/L(K(rr,1),jj)) + eox*(phi(V2,1)-phi(V1,1))*(Area(K(rr,1),jj-1)/L(K(rr,1),jj-1));
                        Jaco(ii, V3) = Jaco(ii, V3)+eox*(Area(K(rr,1),jj)/L(K(rr,1),jj));
                        Jaco(ii, V2) = Jaco(ii, V2)+eox*(Area(K(rr,1),jj-1)/L(K(rr,1),jj-1));
                        Jaco(ii, V1) = Jaco(ii, V1)-eox*(Area(K(rr,1),jj)/L(K(rr,1),jj))-eox*(Area(K(rr,1),jj-1)/L(K(rr,1),jj-1));
    
                    end
                end
            end

        elseif Table(ii,1) == 2     %region2

            for kk = 1:rr_row

                if Table(ii,2) >= iy1(kk,1) && Table(ii,2) <= iy2(kk,1)

                    K = find(E_R2 == Table(ii,2));
                    K_row = size(K,1);


                    if Table(ii,3) == Potential

                        for rr = 1:K_row

                            if K(rr,1) <= R_row2
                                jj = 1;

                                V3 = find(Table(1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj+2)) ;
                                V2 = find(Table(1+1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj+1));
                                V1 = find(Table(1+1+1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj));

                                res(ii,1) = res(ii,1) + esi*(phi(VR1+V2,1)-phi(VR1+V1,1))*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj)) + esi*(phi(VR1+V3,1)-phi(VR1+V1,1))*(Area(K(rr,1)+R_row1,jj+2)/L(K(rr,1)+R_row1,jj+2));
                                res(ii,1) = res(ii,1) + coeff*(Area(K(rr,1)+R_row1,jj)*L(K(rr,1)+R_row1,jj)+Area(K(rr,1)+R_row1,jj+2)*L(K(rr,1)+R_row1,jj+2))*(-Na-elec(VR1+V1,1)+hole(VR1+V1,1));

                                Jaco(ii, 3*(V3)+(VR1-2)) = Jaco(ii, 3*(V3)+(VR1-2))+ esi*(Area(K(rr,1)+R_row1,jj+2)/L(K(rr,1)+R_row1,jj+2));
                                Jaco(ii, 3*(V2)+(VR1-2)) = Jaco(ii, 3*(V2)+(VR1-2))+esi*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj));
                                Jaco(ii, 3*(V1)+(VR1-2)) = Jaco(ii, 3*(V1)+(VR1-2))-esi*(Area(K(rr,1)+R_row1,jj+2)/L(K(rr,1)+R_row1,jj+2))-esi*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj));

                                Jaco(ii, 3*(V1)+(VR1-1)) = Jaco(ii, 3*(V1)+(VR1-1)) - coeff*(Area(K(rr,1)+R_row1,jj)*L(K(rr,1)+R_row1,jj)+Area(K(rr,1)+R_row1,jj+2)*L(K(rr,1)+R_row1,jj+2));
                                Jaco(ii, 3*(V1)+(VR1)) = Jaco(ii, 3*(V1)+(VR1)) + coeff*(Area(K(rr,1)+R_row1,jj)*L(K(rr,1)+R_row1,jj)+Area(K(rr,1)+R_row1,jj+2)*L(K(rr,1)+R_row1,jj+2));

                            elseif K(rr,1) > R_row2 && K(rr,1) <= 2*R_row2
                                jj = 2;
                                K(rr,1) = K(rr,1)-R_row2;

                                V3 = find(Table(1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj+1)) ;
                                V2 = find(Table(1+1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj));
                                V1 = find(Table(1+1+1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj-1));

                                res(ii,1) = res(ii,1) + esi*(phi(VR1+V1,1)-phi(VR1+V2,1))*(Area(K(rr,1)+R_row1,jj-1)/L(K(rr,1)+R_row1,jj-1)) + esi*(phi(VR1+V3,1)-phi(VR1+V2,1))*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj));
                                res(ii,1) = res(ii,1) + coeff*(Area(K(rr,1)+R_row1,jj)*L(K(rr,1)+R_row1,jj)+Area(K(rr,1)+R_row1,jj-1)*L(K(rr,1)+R_row1,jj-1))*(-Na-elec(VR1+V2,1)+hole(VR1+V2,1));

                                Jaco(ii, 3*(V3)+(VR1-2)) = Jaco(ii, 3*(V3)+(VR1-2))+esi*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj));
                                Jaco(ii, 3*(V2)+(VR1-2)) =  Jaco(ii, 3*(V2)+(VR1-2))-esi*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj))-esi*(Area(K(rr,1)+R_row1,jj-1)/L(K(rr,1)+R_row1,jj-1));
                                Jaco(ii, 3*(V1)+(VR1-2)) = Jaco(ii, 3*(V1)+(VR1-2))+esi*(Area(K(rr,1)+R_row1,jj-1)/L(K(rr,1)+R_row1,jj-1));

                                Jaco(ii, 3*(V2)+(VR1-1)) = Jaco(ii, 3*(V2)+(VR1-1)) - coeff*(Area(K(rr,1)+R_row1,jj)*L(K(rr,1)+R_row1,jj)+Area(K(rr,1)+R_row1,jj-1)*L(K(rr,1)+R_row1,jj-1));
                                Jaco(ii, 3*(V2)+(VR1)) = Jaco(ii, 3*(V2)+(VR1)) + coeff*(Area(K(rr,1)+R_row1,jj)*L(K(rr,1)+R_row1,jj)+Area(K(rr,1)+R_row1,jj-1)*L(K(rr,1)+R_row1,jj-1));

                            elseif K(rr,1) > 2*R_row2 && K(rr,1) <= 3*R_row2
                                jj = 3;
                                K(rr,1) = K(rr,1)-2*R_row2;

                                V3 = find(Table(1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj-2)) ;
                                V2 = find(Table(1+1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj-1));
                                V1 = find(Table(1+1+1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj));

                                res(ii,1) = res(ii,1) + esi*(phi(VR1+V3,1)-phi(VR1+V1,1))*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj)) + esi*(phi(VR1+V2,1)-phi(VR1+V1,1))*(Area(K(rr,1)+R_row1,jj-1)/L(K(rr,1)+R_row1,jj-1));
                                res(ii,1) = res(ii,1) + coeff*(Area(K(rr,1)+R_row1,jj)*L(K(rr,1)+R_row1,jj)+Area(K(rr,1)+R_row1,jj-1)*L(K(rr,1)+R_row1,jj-1))*(-Na-elec(VR1+V1,1)+hole(VR1+V1,1));

                                Jaco(ii, 3*(V3)+(VR1-2)) = Jaco(ii, 3*(V3)+(VR1-2))+esi*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj));
                                Jaco(ii, 3*(V2)+(VR1-2)) = Jaco(ii, 3*(V2)+(VR1-2))+esi*(Area(K(rr,1)+R_row1,jj-1)/L(K(rr,1)+R_row1,jj-1));
                                Jaco(ii, 3*(V1)+(VR1-2)) = Jaco(ii, 3*(V1)+(VR1-2))-esi*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj))-esi*(Area(K(rr,1)+R_row1,jj-1)/L(K(rr,1)+R_row1,jj-1));

                                Jaco(ii, 3*(V1)+(VR1-1)) = Jaco(ii, 3*(V1)+(VR1-1)) - coeff*(Area(K(rr,1)+R_row1,jj)*L(K(rr,1)+R_row1,jj)+Area(K(rr,1)+R_row1,jj-1)*L(K(rr,1)+R_row1,jj-1));
                                Jaco(ii, 3*(V1)+(VR1)) = Jaco(ii, 3*(V1)+(VR1)) + coeff*(Area(K(rr,1)+R_row1,jj)*L(K(rr,1)+R_row1,jj)+Area(K(rr,1)+R_row1,jj-1)*L(K(rr,1)+R_row1,jj-1));
                            end
                        end

                    elseif Table(ii,3) == eDensity

                        for rr = 1:K_row
                            if K(rr,1) <= R_row2
                                jj = 1;

                                V3 = find(Table(1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj+2)) ;  % node 3 
                                V2 = find(Table(1+1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj+1)); % node 2
                                V1 = find(Table(2+1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj));   % node 1

                                res(ii,1) =  elec(VR1+V1,1)-nint*exp(phi(VR1+V1,1)/Thermal_V);

                                Jaco(ii, 3*(V1)+(VR1-1)) = 1 ;
                                Jaco(ii, 3*(V1)+(VR1-2)) = -(nint/Thermal_V)*exp(phi(VR1+V1,1)/Thermal_V);


                            elseif K(rr,1) > R_row2 && K(rr,1) <= 2*R_row2
                                jj = 2;
                                K(rr,1) = K(rr,1)-R_row2;

                                V3 = find(Table(1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj+1)) ; % node 3
                                V2 = find(Table(1+1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj));  % node 2
                                V1 = find(Table(1+1+1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj-1)); % node 1

                                res(ii,1) =  elec(VR1+V2,1) - nint*exp(phi(VR1+V2,1)/Thermal_V);

                                Jaco(ii, 3*(V2)+(VR1-1)) = 1;
                                Jaco(ii, 3*(V2)+(VR1-2)) = -(nint/Thermal_V)*exp(phi(VR1+V2,1)/Thermal_V);


                            elseif K(rr,1) > 2*R_row2 && K(rr,1) <= 3*R_row2
                                jj = 3; 
                                K(rr,1) = K(rr,1)-2*R_row2;

                                V3 = find(Table(1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj-2)) ;  % node 1
                                V2 = find(Table(1+1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj-1)); % node 2
                                V1 = find(Table(1+1+1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj)); % node 3

                                res(ii,1) =  elec(VR1+V1,1) - nint*exp(phi(VR1+V1,1)/Thermal_V);

                                Jaco(ii, 3*(V1)+(VR1-1)) = 1;
                                Jaco(ii, 3*(V1)+(VR1-2)) = -(nint/Thermal_V)*exp(phi(VR1+V1,1)/Thermal_V);

                            end
                        end


                    elseif Table(ii,3) == hDensity

                        for rr = 1:K_row

                            if K(rr,1) <= R_row2
                                jj = 1;

                                V3 = find(Table(1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj+2)) ;
                                V2 = find(Table(1+1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj+1));
                                V1 = find(Table(2+1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj));

                                res(ii,1) = hole(VR1+V1,1)-nint*exp(-phi(VR1+V1,1)/Thermal_V);

                                Jaco(ii, 3*(V1)+(VR1)) = 1;
                                Jaco(ii, 3*(V1)+(VR1-2)) = (nint/Thermal_V)*exp(-phi(VR1+V1,1)/Thermal_V);


                            elseif K(rr,1) > R_row2 && K(rr,1) <= 2*R_row2
                                jj = 2;
                                K(rr,1) = K(rr,1)-R_row2;

                                V3 = find(Table(1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj+1)) ;
                                V2 = find(Table(1+1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj));
                                V1 = find(Table(1+1+1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj-1));

                                res(ii,1) =  hole(VR1+V2,1)-nint*exp(-phi(VR1+V2,1)/Thermal_V);

                                Jaco(ii, 3*(V2)+(VR1)) =  1;
                                Jaco(ii, 3*(V2)+(VR1-2)) = (nint/Thermal_V)*exp(-phi(VR1+V2,1)/Thermal_V);


                            elseif K(rr,1) > 2*R_row2 && K(rr,1) <= 3*R_row2
                                jj = 3;
                                K(rr,1) = K(rr,1)-2*R_row2;

                                V3 = find(Table(1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj-2)) ;
                                V2 = find(Table(1+1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj-1));
                                V1 = find(Table(1+1+1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj));

                                res(ii,1) =  hole(VR1+V1,1)-nint*exp(-phi(VR1+V1,1)/Thermal_V);

                                Jaco(ii, 3*(V1)+(VR1)) =  1;
                                Jaco(ii, 3*(V1)+(VR1-2)) =  (nint/Thermal_V)*exp(-phi(VR1+V1,1)/Thermal_V);
                            end
                        end
                    
                    end
                end

                if (ry1(kk,1) <= Table(ii,2) && Table(ii,2) <= iy1(kk,1)) || (ry2(kk,1) >= Table(ii,2) && Table(ii,2) >= iy2(kk,1))

                    K = find(E_R2 == Table(ii,2));
                    K_row = size(K,1);
                    if Table(ii,3) == Potential

                        for rr = 1:K_row

                            if K(rr,1) <= R_row2
                                jj = 1;

                                V3 = find(Table(1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj+2)) ;
                                V2 = find(Table(1+1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj+1));
                                V1 = find(Table(1+1+1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj));

                                res(ii,1) = res(ii,1) + esi*(phi(VR1+V2,1)-phi(VR1+V1,1))*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj)) + esi*(phi(VR1+V3,1)-phi(VR1+V1,1))*(Area(K(rr,1)+R_row1,jj+2)/L(K(rr,1)+R_row1,jj+2));
                                res(ii,1) = res(ii,1) + coeff*(Area(K(rr,1)+R_row1,jj)*L(K(rr,1)+R_row1,jj)+Area(K(rr,1)+R_row1,jj+2)*L(K(rr,1)+R_row1,jj+2))*(Nd-elec(VR1+V1,1)+hole(VR1+V1,1));

                                Jaco(ii, 3*(V3)+(VR1-2)) = Jaco(ii, 3*(V3)+(VR1-2))+ esi*(Area(K(rr,1)+R_row1,jj+2)/L(K(rr,1)+R_row1,jj+2));
                                Jaco(ii, 3*(V2)+(VR1-2)) = Jaco(ii, 3*(V2)+(VR1-2))+esi*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj));
                                Jaco(ii, 3*(V1)+(VR1-2)) = Jaco(ii, 3*(V1)+(VR1-2))-esi*(Area(K(rr,1)+R_row1,jj+2)/L(K(rr,1)+R_row1,jj+2))-esi*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj));

                                Jaco(ii, 3*(V1)+(VR1-1)) = Jaco(ii, 3*(V1)+(VR1-1)) - coeff*(Area(K(rr,1)+R_row1,jj)*L(K(rr,1)+R_row1,jj)+Area(K(rr,1)+R_row1,jj+2)*L(K(rr,1)+R_row1,jj+2));
                                Jaco(ii, 3*(V1)+(VR1)) = Jaco(ii, 3*(V1)+(VR1)) + coeff*(Area(K(rr,1)+R_row1,jj)*L(K(rr,1)+R_row1,jj)+Area(K(rr,1)+R_row1,jj+2)*L(K(rr,1)+R_row1,jj+2));

                            elseif K(rr,1) > R_row2 && K(rr,1) <= 2*R_row2
                                jj = 2;
                                K(rr,1) = K(rr,1)-R_row2;

                                V3 = find(Table(1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj+1)) ;
                                V2 = find(Table(1+1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj));
                                V1 = find(Table(1+1+1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj-1));

                                res(ii,1) = res(ii,1) + esi*(phi(VR1+V1,1)-phi(VR1+V2,1))*(Area(K(rr,1)+R_row1,jj-1)/L(K(rr,1)+R_row1,jj-1)) + esi*(phi(VR1+V3,1)-phi(VR1+V2,1))*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj));
                                res(ii,1) = res(ii,1) + coeff*(Area(K(rr,1)+R_row1,jj)*L(K(rr,1)+R_row1,jj)+Area(K(rr,1)+R_row1,jj-1)*L(K(rr,1)+R_row1,jj-1))*(Nd-elec(VR1+V2,1)+hole(VR1+V2,1));

                                Jaco(ii, 3*(V3)+(VR1-2)) = Jaco(ii, 3*(V3)+(VR1-2))+esi*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj));
                                Jaco(ii, 3*(V2)+(VR1-2)) =  Jaco(ii, 3*(V2)+(VR1-2))-esi*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj))-esi*(Area(K(rr,1)+R_row1,jj-1)/L(K(rr,1)+R_row1,jj-1));
                                Jaco(ii, 3*(V1)+(VR1-2)) = Jaco(ii, 3*(V1)+(VR1-2))+esi*(Area(K(rr,1)+R_row1,jj-1)/L(K(rr,1)+R_row1,jj-1));

                                Jaco(ii, 3*(V2)+(VR1-1)) = Jaco(ii, 3*(V2)+(VR1-1)) - coeff*(Area(K(rr,1)+R_row1,jj)*L(K(rr,1)+R_row1,jj)+Area(K(rr,1)+R_row1,jj-1)*L(K(rr,1)+R_row1,jj-1));
                                Jaco(ii, 3*(V2)+(VR1)) = Jaco(ii, 3*(V2)+(VR1)) + coeff*(Area(K(rr,1)+R_row1,jj)*L(K(rr,1)+R_row1,jj)+Area(K(rr,1)+R_row1,jj-1)*L(K(rr,1)+R_row1,jj-1));

                            elseif K(rr,1) > 2*R_row2 && K(rr,1) <= 3*R_row2
                                jj = 3;
                                K(rr,1) = K(rr,1)-2*R_row2;

                                V3 = find(Table(1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj-2)) ;
                                V2 = find(Table(1+1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj-1));
                                V1 = find(Table(1+1+1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj));

                                res(ii,1) = res(ii,1) + esi*(phi(VR1+V3,1)-phi(VR1+V1,1))*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj)) + esi*(phi(VR1+V2,1)-phi(VR1+V1,1))*(Area(K(rr,1)+R_row1,jj-1)/L(K(rr,1)+R_row1,jj-1));
                                res(ii,1) = res(ii,1) + coeff*(Area(K(rr,1)+R_row1,jj)*L(K(rr,1)+R_row1,jj)+Area(K(rr,1)+R_row1,jj-1)*L(K(rr,1)+R_row1,jj-1))*(Nd-elec(VR1+V1,1)+hole(VR1+V1,1));

                                Jaco(ii, 3*(V3)+(VR1-2)) = Jaco(ii, 3*(V3)+(VR1-2))+esi*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj));
                                Jaco(ii, 3*(V2)+(VR1-2)) = Jaco(ii, 3*(V2)+(VR1-2))+esi*(Area(K(rr,1)+R_row1,jj-1)/L(K(rr,1)+R_row1,jj-1));
                                Jaco(ii, 3*(V1)+(VR1-2)) = Jaco(ii, 3*(V1)+(VR1-2))-esi*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj))-esi*(Area(K(rr,1)+R_row1,jj-1)/L(K(rr,1)+R_row1,jj-1));

                                Jaco(ii, 3*(V1)+(VR1-1)) = Jaco(ii, 3*(V1)+(VR1-1)) - coeff*(Area(K(rr,1)+R_row1,jj)*L(K(rr,1)+R_row1,jj)+Area(K(rr,1)+R_row1,jj-1)*L(K(rr,1)+R_row1,jj-1));
                                Jaco(ii, 3*(V1)+(VR1)) = Jaco(ii, 3*(V1)+(VR1)) + coeff*(Area(K(rr,1)+R_row1,jj)*L(K(rr,1)+R_row1,jj)+Area(K(rr,1)+R_row1,jj-1)*L(K(rr,1)+R_row1,jj-1));
                            end
                        end

                    elseif Table(ii,3) == eDensity

                        for rr = 1:K_row
                            if K(rr,1) <= R_row2
                                jj = 1;

                                V3 = find(Table(1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj+2)) ;  % node 3 
                                V2 = find(Table(1+1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj+1)); % node 2
                                V1 = find(Table(2+1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj));   % node 1

                                res(ii,1) =  elec(VR1+V1,1)-nint*exp(phi(VR1+V1,1)/Thermal_V);

                                Jaco(ii, 3*(V1)+(VR1-1)) = 1 ;
                                Jaco(ii, 3*(V1)+(VR1-2)) = -(nint/Thermal_V)*exp(phi(VR1+V1,1)/Thermal_V);


                            elseif K(rr,1) > R_row2 && K(rr,1) <= 2*R_row2
                                jj = 2;
                                K(rr,1) = K(rr,1)-R_row2;

                                V3 = find(Table(1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj+1)) ; % node 3
                                V2 = find(Table(1+1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj));  % node 2
                                V1 = find(Table(1+1+1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj-1)); % node 1

                                res(ii,1) =  elec(VR1+V2,1) - nint*exp(phi(VR1+V2,1)/Thermal_V);

                                Jaco(ii, 3*(V2)+(VR1-1)) = 1;
                                Jaco(ii, 3*(V2)+(VR1-2)) = -(nint/Thermal_V)*exp(phi(VR1+V2,1)/Thermal_V);


                            elseif K(rr,1) > 2*R_row2 && K(rr,1) <= 3*R_row2
                                jj = 3; 
                                K(rr,1) = K(rr,1)-2*R_row2;

                                V3 = find(Table(1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj-2)) ;  % node 1
                                V2 = find(Table(1+1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj-1)); % node 2
                                V1 = find(Table(1+1+1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj)); % node 3

                                res(ii,1) =  elec(VR1+V1,1) - nint*exp(phi(VR1+V1,1)/Thermal_V);

                                Jaco(ii, 3*(V1)+(VR1-1)) = 1;
                                Jaco(ii, 3*(V1)+(VR1-2)) = -(nint/Thermal_V)*exp(phi(VR1+V1,1)/Thermal_V);

                            end
                        end


                    elseif Table(ii,3) == hDensity

                        for rr = 1:K_row

                            if K(rr,1) <= R_row2
                                jj = 1;

                                V3 = find(Table(1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj+2)) ;
                                V2 = find(Table(1+1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj+1));
                                V1 = find(Table(2+1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj));

                                res(ii,1) = hole(VR1+V1,1)-nint*exp(-phi(VR1+V1,1)/Thermal_V);

                                Jaco(ii, 3*(V1)+(VR1)) = 1;
                                Jaco(ii, 3*(V1)+(VR1-2)) = (nint/Thermal_V)*exp(-phi(VR1+V1,1)/Thermal_V);


                            elseif K(rr,1) > R_row2 && K(rr,1) <= 2*R_row2
                                jj = 2;
                                K(rr,1) = K(rr,1)-R_row2;

                                V3 = find(Table(1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj+1)) ;
                                V2 = find(Table(1+1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj));
                                V1 = find(Table(1+1+1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj-1));

                                res(ii,1) =  hole(VR1+V2,1)-nint*exp(-phi(VR1+V2,1)/Thermal_V);

                                Jaco(ii, 3*(V2)+(VR1)) =  1;
                                Jaco(ii, 3*(V2)+(VR1-2)) = (nint/Thermal_V)*exp(-phi(VR1+V2,1)/Thermal_V);


                            elseif K(rr,1) > 2*R_row2 && K(rr,1) <= 3*R_row2
                                jj = 3;
                                K(rr,1) = K(rr,1)-2*R_row2;

                                V3 = find(Table(1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj-2)) ;
                                V2 = find(Table(1+1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj-1));
                                V1 = find(Table(1+1+1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj));

                                res(ii,1) =  hole(VR1+V1,1)-nint*exp(-phi(VR1+V1,1)/Thermal_V);

                                Jaco(ii, 3*(V1)+(VR1)) =  1;
                                Jaco(ii, 3*(V1)+(VR1-2)) =  (nint/Thermal_V)*exp(-phi(VR1+V1,1)/Thermal_V);
                            end
                        end
                    end
                end
            end
        
        elseif Table(ii,1) == 3  % Region3
    
            K = find(E_R3 == Table(ii,2));
            K_row = size(K,1);
    
            for rr = 1:K_row

                if Table(ii,3) == Potential
    
                    if K(rr,1) <= R_row3
    
                        jj = 1;
    
                        V3 = find(Table(1+VR1+3*VR2:T1,2)==E_R3(K(rr,1),jj+2)) ;
                        V2 = find(Table(1+VR1+3*VR2:T1,2)==E_R3(K(rr,1),jj+1));
                        V1 = find(Table(1+VR1+3*VR2:T1,2)==E_R3(K(rr,1),jj));
    
                        res(ii,1) = res(ii,1) + eox*(phi(VR1+VR2+V2,1)-phi(VR1+VR2+V1,1))*(Area(K(rr,1)+R_row2+R_row1,jj)/L(K(rr,1)+R_row2+R_row1,jj)) + eox*(phi(VR1+VR2+V3,1)-phi(VR1+VR2+V1,1))*(Area(K(rr,1)+R_row2+R_row1,jj+2)/L(K(rr,1)+R_row2+R_row1,jj+2));
    
                        Jaco(ii, VR1+3*VR2+V3) = Jaco(ii, VR1+3*VR2+V3)+ eox*(Area(K(rr,1)+R_row2+R_row1,jj+2)/L(K(rr,1)+R_row2+R_row1,jj+2));
                        Jaco(ii, VR1+3*VR2+V2) = Jaco(ii, VR1+3*VR2+V2)+eox*(Area(K(rr,1)+R_row2+R_row1,jj)/L(K(rr,1)+R_row2+R_row1,jj));
                        Jaco(ii, VR1+3*VR2+V1) = Jaco(ii, VR1+3*VR2+V1)-eox*(Area(K(rr,1)+R_row2+R_row1,jj+2)/L(K(rr,1)+R_row2+R_row1,jj+2))-eox*(Area(K(rr,1)+R_row2+R_row1,jj)/L(K(rr,1)+R_row2+R_row1,jj));
    
                    elseif K(rr,1) > R_row3 && K(rr,1) <= 2*R_row3
                        
                        jj = 2;
                        K(rr,1) = K(rr,1)-R_row3;
    
                        V3 = find(Table(1+VR1+3*VR2:T1,2)==E_R3(K(rr,1),jj+1)) ;
                        V2 = find(Table(1+VR1+3*VR2:T1,2)==E_R3(K(rr,1),jj));
                        V1 = find(Table(1+VR1+3*VR2:T1,2)==E_R3(K(rr,1),jj-1));
    
                        res(ii,1) = res(ii,1) + eox*(phi(VR1+VR2+V1,1)-phi(VR1+VR2+V2,1))*(Area(K(rr,1)+R_row2+R_row1,jj-1)/L(K(rr,1)+R_row2+R_row1,jj-1)) + eox*(phi(VR1+VR2+V3,1)-phi(VR1+VR2+V2,1))*(Area(K(rr,1)+R_row2+R_row1,jj)/L(K(rr,1)+R_row2+R_row1,jj));
    
                        Jaco(ii, VR1+3*VR2+V3) = Jaco(ii, VR1+3*VR2+V3)+eox*(Area(K(rr,1)+R_row2+R_row1,jj)/L(K(rr,1)+R_row2+R_row1,jj));
                        Jaco(ii, VR1+3*VR2+V2) =  Jaco(ii, VR1+3*VR2+V2)-eox*(Area(K(rr,1)+R_row2+R_row1,jj)/L(K(rr,1)+R_row2+R_row1,jj))-eox*(Area(K(rr,1)+R_row2+R_row1,jj-1)/L(K(rr,1)+R_row2+R_row1,jj-1));
                        Jaco(ii, VR1+3*VR2+V1) = Jaco(ii, VR1+3*VR2+V1)+eox*(Area(K(rr,1)+R_row2+R_row1,jj-1)/L(K(rr,1)+R_row2+R_row1,jj-1));
    
    
                    elseif K(rr,1) > 2*R_row3 && K(rr,1) <= 3*R_row3
                        
                        jj = 3;
                        K(rr,1) = K(rr,1)-2*R_row3;
    
                        V3 = find(Table(1+VR1+3*VR2:T1,2)==E_R3(K(rr,1),jj-2)) ;
                        V2 = find(Table(1+VR1+3*VR2:T1,2)==E_R3(K(rr,1),jj-1));
                        V1 = find(Table(1+VR1+3*VR2:T1,2)==E_R3(K(rr,1),jj));

                        res(ii,1) = res(ii,1) + eox*(phi(VR1+VR2+V3,1)-phi(VR1+VR2+V1,1))*(Area(K(rr,1)+R_row2+R_row1,jj)/L(K(rr,1)+R_row2+R_row1,jj)) + eox*(phi(VR1+VR2+V2,1)-phi(VR1+VR2+V1,1))*(Area(K(rr,1)+R_row2+R_row1,jj-1)/L(K(rr,1)+R_row2+R_row1,jj-1));

                        Jaco(ii, VR1+3*VR2+V3) = Jaco(ii, VR1+3*VR2+V3) + eox*(Area(K(rr,1)+R_row2+R_row1,jj)/L(K(rr,1)+R_row2+R_row1,jj));
                        Jaco(ii, VR1+3*VR2+V2) = Jaco(ii, VR1+3*VR2+V2) + eox*(Area(K(rr,1)+R_row2+R_row1,jj-1)/L(K(rr,1)+R_row2+R_row1,jj-1));
                        Jaco(ii, VR1+3*VR2+V1) = Jaco(ii, VR1+3*VR2+V1) - eox*(Area(K(rr,1)+R_row2+R_row1,jj)/L(K(rr,1)+R_row2+R_row1,jj))-eox*(Area(K(rr,1)+R_row2+R_row1,jj-1)/L(K(rr,1)+R_row2+R_row1,jj-1));

                    end
                end
            end
        end
    end

%    Boundary Condition

    for rr = 1:i_row1

        Vi11 = find(Table(1:VR1,2)==interface1(rr,1));
        Vi21 = find(Table(1+VR1:3:VR1+3*VR2,2)==interface1(rr,1));

        Jaco(3*Vi21+(VR1-2),:)=Jaco(3*Vi21+(VR1-2),:)+Jaco(Vi11,:);
        Jaco(Vi11,:)=0; Jaco(Vi11,Vi11)=1; Jaco(Vi11,3*Vi21+(VR1-2))=-1;

        res(3*Vi21+(VR1-2),1)=res(3*Vi21+(VR1-2),1)+res(Vi11,1);
        res(Vi11,1)=0;
    end

    for rr = 1:i_row2

        Vi12 = find(Table(1+VR1+3*VR2:T1,2)==interface2(rr,1));
        Vi22 = find(Table(1+VR1:3:VR1+3*VR2,2)==interface2(rr,1));

        Jaco(3*Vi22+(VR1-2),:)=Jaco(3*Vi22+(VR1-2),:)+Jaco(VR1+3*VR2+Vi12,:);
        Jaco(VR1+3*VR2+Vi12,:)=0; Jaco(VR1+3*VR2+Vi12,VR1+3*VR2+Vi12)=1; Jaco(VR1+3*VR2+Vi12,3*Vi22+(VR1-2))=-1;

        res(3*Vi22+(VR1-2),1)=res(3*Vi22+(VR1-2),1)+res(VR1+3*VR2+Vi12,1);
        res(VR1+3*VR2+Vi12,1)=0;

    end

    %%%% Dirichlet Boundary Condition %%%%%

    

    for n=1:C_col
        for i =1: C_row
            if i <= (C_row/2) % Bottom Gate
                V = find(Table(1:VR1,2)==Contact(i,n));

                Jaco(V,:) =0;
                Jaco(V,V) = 1;
                res(V,1)=phi(V,1)-0.33374;

            else % Top Gate

                V = find(Table(1+VR1+3*VR2:T_row,2)==Contact(i,n));

                Jaco(VR1+3*VR2+V,:) =0;
                Jaco(VR1+3*VR2+V,VR1+3*VR2+V) = 1;
                res(VR1+3*VR2+V,1)= phi(VR1+VR2+V,1)-0.33374;
            end
        end
    end

    %%%% Source, Drain %%%%
    q_row = size(ry1,1);
    w_row = size(ry2,1);

    for rr = 1:q_row  % Source

        V = find(Table(1+VR1:3:VR1+3*VR2,2)==ry1(rr,1));
        Jaco(3*V+(VR1-2),:) = 0;
        Jaco(3*V+(VR1-2),3*V+(VR1-2)) = 1;
        Jaco(3*V+(VR1-1),:) = 0;
        Jaco(3*V+(VR1-1),3*V+(VR1-1)) = 1;
        Jaco(3*V+(VR1),:) = 0;
        Jaco(3*V+(VR1),3*V+(VR1)) = 1;

        res(3*V+(VR1-2),1) = phi(VR1+V,1)-Thermal_V*log(Nd/nint);
        res(3*V+(VR1-1),1) = elec(VR1+V,1)-Nd;
        res(3*V+(VR1),1) = hole(VR1+V,1)-(nint^2/Nd);

    end

    for rr = 1:w_row  % Drain

        V = find(Table(1+VR1:3:VR1+3*VR2,2)==ry2(rr,1));
 
        Jaco(3*V+(VR1-2),:) = 0;
        Jaco(3*V+(VR1-2),3*V+(VR1-2)) = 1;
        Jaco(3*V+(VR1-1),:) = 0;
        Jaco(3*V+(VR1-1),3*V+(VR1-1)) = 1;
        Jaco(3*V+(VR1),:) = 0;
        Jaco(3*V+(VR1),3*V+(VR1)) = 1;

        res(3*V+(VR1-2),1) = phi(VR1+V,1)-Thermal_V*log(Nd/nint);
        res(3*V+(VR1-1),1) = elec(VR1+V,1)-Nd;
        res(3*V+(VR1),1) = hole(VR1+V,1)-(nint^2/Nd);

    end

        %     scaling
        Cvector = zeros(T_row,1);
        Cvector(1:VR1,1) = Thermal_V;
        Cvector(VR1+3*VR2:T_row,1) = Thermal_V;
        Cvector(1+VR1:3:VR1+3*VR2,1) = Thermal_V;
        Cvector(2+VR1:3:VR1+3*VR2,1) = Nd;
        Cvector(3+VR1:3:VR1+3*VR2,1) = Nd;

        Cmatrix = spdiags(Cvector,0,T_row,T_row);
        Jaco_scaled = Jaco * Cmatrix;
        Rvector = 1./sum(abs(Jaco_scaled),2);
        Rmatrix = spdiags(Rvector,0,T_row,T_row);
        Jaco_scaled = Rmatrix* Jaco_scaled;
        res_scaled = Rmatrix *res;
        update_scaled = Jaco_scaled \ (-res_scaled);
        update_vector(:,Newton) = Cmatrix* update_scaled;

        phi(1:VR1,1) = phi(1:VR1,1) + update_vector(1:VR1,Newton);
        phi(1+VR1:VR2+VR1,1) = phi(1+VR1:VR2+VR1,1) + update_vector(1+VR1:3:VR1+3*VR2,Newton);
        phi(1+VR1+VR2:T,1) = phi(1+VR1+VR2:T,1) + update_vector(1+VR1+3*VR2:T_row,Newton);

        elec(1+VR1:1:VR2+VR1,1) = elec(1+VR1:1:VR2+VR1,1) + update_vector(2+VR1:3:VR1+3*VR2,Newton);
        hole(1+VR1:1:VR2+VR1,1) = hole(1+VR1:1:VR2+VR1,1) + update_vector(3+VR1:3:VR1+3*VR2,Newton);

        update_poisson(:,Newton) = abs(update_vector(1+VR1:3:VR1+3*VR2,Newton));

        if norm(abs(update_vector(1+VR1:3:VR1+3*VR2,Newton)),inf) < 1e-15
            break;
        end
end

set(gcf,'Color','w')
patch('Faces',E_R1,'Vertices',Vertex,'EdgeColor','black','FaceColor','Green','LineWidth',2);
patch('Faces',E_R2,'Vertices',Vertex,'EdgeColor','black','FaceColor','blue','LineWidth',2);
patch('Faces',E_R3,'Vertices',Vertex,'EdgeColor','black','FaceColor','Green','LineWidth',2);
patch('Faces',Edge1,'Vertices',Vertex,'EdgeColor','yellow','FaceColor','none','LineWidth',3);
patch('Faces',Edge2,'Vertices',Vertex,'EdgeColor','yellow','FaceColor','none','LineWidth',3);
patch('Faces',Contact,'Vertices',Vertex,'EdgeColor','Cyan','FaceColor','none','LineWidth',3);
patch('Faces',Contact_Si,'Vertices',Vertex,'EdgeColor','magenta','FaceColor','none','LineWidth',3);

%%%%%% Drift-Diffusion %%%%%%%%

T1 = VR1+3*VR2+VR3;
Table = zeros(T1,3);

T_row = size(Table,1);

for iRegion = 1:3
    if iRegion ==1  % Oxide , only potential

        Table(1:VR1,1) = 1;

        for iVertex =  1:VR1
            Table(iVertex,2) = V_R1(iVertex,1);
            Table(iVertex,3) = Potential;   % only Potential
        end

    elseif iRegion == 2     % Silicon, Potential, eDensity, hDensity

        Table(1+VR1:VR1+3*VR2,1) = 2;

        for iVertex =  1:VR2

            Table(3*iVertex+(VR1-2),2) = V_R2(iVertex,1);
            Table(3*iVertex+(VR1-1),2) = V_R2(iVertex,1);
            Table(3*iVertex+(VR1),2) = V_R2(iVertex,1);

            Table(3*iVertex+(VR1-2),3) = Potential;
            Table(3*iVertex+(VR1-1),3) = eDensity;
            Table(3*iVertex+(VR1),3) = hDensity;
        end

    elseif iRegion ==3  % Oxide , only potential

        Table(1+VR1+3*VR2:T1,1) = 3;
        
        for iVertex = 1:VR3
            Table(VR1+3*VR2+iVertex,2) = V_R3(iVertex,1);
            Table(VR1+3*VR2+iVertex,3) = Potential;
        end
         
    end
end