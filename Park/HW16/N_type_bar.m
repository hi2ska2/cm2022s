clear; close all; clc;

% due to : 2022.05.5. (Thur) / HW16

%%%%% HW16 %%%%%%%
% Transient simulation 2

%%%% parameter %%%%

q=1.602e-19;
nint=1.075e+16; 
Na = -2e+23;          %1/m^3   N+-type
Na1 = 2e+23;          %1/m^3   N-type
k_B=1.380662e-23;
T=300;      %K
esi = 11.7; eox = 3.9; e0=8.854187817e-12;  % Permittivity, F/m
%m_n = 1417e-4;      % electron mobility
m_n = 518e-4;      % electron mobility
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

Contact_Si = importdata('Contact_Si.txt');

index = V_row;

E_R2 = Element;

R_row1 = 0;
R_row2 = size(E_R2,1);
R_row3 = 0;

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

V_R1 = 0;
V_R2 = unique(E_R2);
V_R3 = 0;

VR1 = 0;
VR2 = size(V_R2,1);
VR3 = 0;

T = VR1+VR2+VR3;

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

%%%%%%%% Initial Poisson %%%%%%%%%%

phi = zeros(T_row,1);
elec = zeros(T_row,1);
hole = zeros(T_row,1);

phi(:,1) = Thermal_V*log(-Na/nint);

ry1 = unique(Contact_Si(1:10,:)); % Source
ry2 = unique(Contact_Si(11:20,:)); % Drain

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
                        res(ii,1) = res(ii,1) + coeff*(Area(K(rr,1)+R_row1,jj)*L(K(rr,1)+R_row1,jj)+Area(K(rr,1)+R_row1,jj+2)*L(K(rr,1)+R_row1,jj+2))*(Na1-elec(VR1+V1,1)+hole(VR1+V1,1));

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
                        res(ii,1) = res(ii,1) + coeff*(Area(K(rr,1)+R_row1,jj)*L(K(rr,1)+R_row1,jj)+Area(K(rr,1)+R_row1,jj-1)*L(K(rr,1)+R_row1,jj-1))*(Na1-elec(VR1+V2,1)+hole(VR1+V2,1));

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
                        res(ii,1) = res(ii,1) + coeff*(Area(K(rr,1)+R_row1,jj)*L(K(rr,1)+R_row1,jj)+Area(K(rr,1)+R_row1,jj-1)*L(K(rr,1)+R_row1,jj-1))*(Na1-elec(VR1+V1,1)+hole(VR1+V1,1));

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

        res(3*V+(VR1-2),1) = phi(VR1+V,1)-Thermal_V*log(-Na/nint);
        res(3*V+(VR1-1),1) = elec(VR1+V,1)+Na;
        res(3*V+(VR1),1) = hole(VR1+V,1)-(nint^2/-Na);

    end

    for rr = 1:w_row  % Drain

        V = find(Table(1+VR1:3:VR1+3*VR2,2)==ry2(rr,1));
 
        Jaco(3*V+(VR1-2),:) = 0;
        Jaco(3*V+(VR1-2),3*V+(VR1-2)) = 1;
        Jaco(3*V+(VR1-1),:) = 0;
        Jaco(3*V+(VR1-1),3*V+(VR1-1)) = 1;
        Jaco(3*V+(VR1),:) = 0;
        Jaco(3*V+(VR1),3*V+(VR1)) = 1;

        res(3*V+(VR1-2),1) = phi(VR1+V,1)-Thermal_V*log(-Na/nint);
        res(3*V+(VR1-1),1) = elec(VR1+V,1)+Na;
        res(3*V+(VR1),1) = hole(VR1+V,1)-(nint^2/-Na);
    end

        %     Scaling
        Cvector = zeros(T_row,1);
        Cvector(1:VR1,1) = Thermal_V;
        Cvector(VR1+3*VR2:T_row,1) = Thermal_V;
        Cvector(1+VR1:3:VR1+3*VR2,1) = Thermal_V;
        Cvector(2+VR1:3:VR1+3*VR2,1) = -Na;
        Cvector(3+VR1:3:VR1+3*VR2,1) = -Na;

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
patch('Faces',E_R2,'Vertices',Vertex,'EdgeColor','black','FaceColor','blue','LineWidth',2);
patch('Faces',Contact_Si,'Vertices',Vertex,'EdgeColor','magenta','FaceColor','none','LineWidth',3);

%%%%%% Drift-Diffusion %%%%%%%%

T1 = VR1+3*VR2+VR3;
Table = zeros(T1,3);

T_row = size(Table,1);

for iRegion = 1:3
    if iRegion ==1  % Oxide , only Potential

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

    elseif iRegion ==3  % Oxide , only Potential

        Table(1+VR1+3*VR2:T1,1) = 3;
        
        for iVertex = 1:VR3
            Table(VR1+3*VR2+iVertex,2) = V_R3(iVertex,1);
            Table(VR1+3*VR2+iVertex,3) = Potential;
        end
         
    end
end