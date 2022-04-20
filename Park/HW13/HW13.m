clear; close all; clc;

% due to : 2022.04.16 (Tue) / HW13

%%%%% HW13 %%%%%%%

%%%% parameter %%%%

q=1.602e-19;
nint=1e+16;     %1/m^3
Na= 2e21;   %1/m^3   P-type
Nd = 5e23 ;  %1/m^3   N-type
k_B=1.38065e-23; % 
T=300; %K
esi = 11.7; eox = 3.9; e0=8.854e-12;  % Permittivity, F/m
m_n = 1417e-4;  % electron mobility
m_p = 470.5e-4;     % Hole mobility
coeff=(0.25*q)/e0;

Thermal_V = k_B*T/q;

coeff_Jn = 1;% q*m_n*Thermal_V; %q?
coeff_Jp = -1;%-q*m_p*Thermal_V; %q?

%%%%%%%% Read Files %%%%%%%%%%
fileID = fopen('Vertex.txt');

vertex = textscan(fileID,'%f%f','Delimiter','\t');

Vertex = horzcat(vertex{1},vertex{2});
fclose(fileID);
[V_row,V_col]=size(Vertex);

fileID1 = fopen('Element.txt');
element = textscan(fileID1,'%f%f%f','Delimiter','\t');

Element = horzcat(element{1},element{2},element{3});
fclose(fileID1);
[E_row,col]=size(Element);

fileID2 = fopen(['Contact.txt']);
contact = textscan(fileID2,'%f%f','Delimiter','\t');

Contact = horzcat(contact{1},contact{2});
fclose(fileID2);
[C_row,C_col]=size(Contact);

index = V_row;

%%%%%%%% Specify Region 1,2 %%%%%%%%%%

fileID2 = fopen('Element_Region1.txt');
element_Region1 = textscan(fileID2,'%f%f%f','Delimiter','\t');

E_R1 = horzcat(element_Region1{1},element_Region1{2},element_Region1{3});
fclose(fileID2);

fileID3 = fopen('Element_Region2.txt');
element_Region2 = textscan(fileID3,'%f%f%f','Delimiter','\t');

E_R2 = horzcat(element_Region2{1},element_Region2{2},element_Region2{3});
fclose(fileID3);

fileID4 = fopen('Element_Region3.txt');
element_Region3 = textscan(fileID4,'%f%f%f','Delimiter','\t');

E_R3 = horzcat(element_Region3{1},element_Region3{2},element_Region3{3});
fclose(fileID4);

fileID3 = fopen('Element_Region4.txt');
element_Region4 = textscan(fileID3,'%f%f%f','Delimiter','\t');

E_R4 = horzcat(element_Region4{1},element_Region4{2},element_Region4{3});
fclose(fileID3);

fileID3 = fopen('Element_Region5.txt');
element_Region5 = textscan(fileID3,'%f%f%f','Delimiter','\t');

E_R5 = horzcat(element_Region5{1},element_Region5{2},element_Region5{3});
fclose(fileID3);

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

%%%%% Region vertex information %%%%%

V_R1 = unique(E_R1);
V_R2 = unique(E_R2);
V_R3 = unique(E_R3);
V_R4 = unique(E_R4);
V_R5 = unique(E_R5);

VR1 = size(V_R1,1);
VR2 = size(V_R2,1);
VR3 = size(V_R3,1);

T = VR1+VR2+VR3;

iy1 = [15; 26; 37; 48];
iy2 = [19; 30; 41; 52];

ry1 = [12; 23; 34; 45];
ry2 = [22; 33; 44; 55];

Area = round(Area*1e14)/1e14;

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

%%%%%%%% Initial Poisson %%%%%%%%%%

%%%%% Need to revise %%%%%

phi = zeros(T_row,1);
elec = zeros(T_row,1);
hole = zeros(T_row,1);

phi(1:11,1) = Thermal_V*log((Nd+Na)/2/nint);
phi(12:34,1) = Thermal_V*log((Nd+Na)/2/nint);

phi(56:77,1) = Thermal_V*log((Nd+Na)/2/nint);
phi(78:T_row,1) = Thermal_V*log((Nd+Na)/2/nint);

phi(35:55,1) = Thermal_V*log((Nd+Na)/2/nint);

elec(23:1:66,1) = nint*exp(phi(23:1:66,1)/Thermal_V);
hole(23:1:66,1) = nint*exp(-phi(23:1:66,1)/Thermal_V);

%%%%%% Fully-Coupled %%%%%%%%
%%% paste before HW

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

            Table(3*iVertex+20,2) = V_R2(iVertex,1);
            Table(3*iVertex+21,2) = V_R2(iVertex,1);
            Table(3*iVertex+22,2) = V_R2(iVertex,1);

            Table(3*iVertex+20,3) = Potential;
            Table(3*iVertex+21,3) = eDensity;
            Table(3*iVertex+22,3) = hDensity;
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

for Newton = 1:10

    Jaco = zeros(T_row,T_row);
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

            for kk = 1:4
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

                                Jaco(ii, 3*(V3)+20) = Jaco(ii, 3*(V3)+20)+ esi*(Area(K(rr,1)+R_row1,jj+2)/L(K(rr,1)+R_row1,jj+2));
                                Jaco(ii, 3*(V2)+20) = Jaco(ii, 3*(V2)+20)+esi*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj));
                                Jaco(ii, 3*(V1)+20) = Jaco(ii, 3*(V1)+20)-esi*(Area(K(rr,1)+R_row1,jj+2)/L(K(rr,1)+R_row1,jj+2))-esi*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj));

                                Jaco(ii, 3*(V1)+21) = Jaco(ii, 3*(V1)+21) - coeff*(Area(K(rr,1)+R_row1,jj)*L(K(rr,1)+R_row1,jj)+Area(K(rr,1)+R_row1,jj+2)*L(K(rr,1)+R_row1,jj+2));
                                Jaco(ii, 3*(V1)+22) = Jaco(ii, 3*(V1)+22) + coeff*(Area(K(rr,1)+R_row1,jj)*L(K(rr,1)+R_row1,jj)+Area(K(rr,1)+R_row1,jj+2)*L(K(rr,1)+R_row1,jj+2));

                            elseif K(rr,1) > R_row2 && K(rr,1) <= 2*R_row2
                                jj = 2;
                                K(rr,1) = K(rr,1)-R_row2;

                                V3 = find(Table(1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj+1)) ;
                                V2 = find(Table(1+1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj));
                                V1 = find(Table(1+1+1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj-1));

                                res(ii,1) = res(ii,1) + esi*(phi(VR1+V1,1)-phi(VR1+V2,1))*(Area(K(rr,1)+R_row1,jj-1)/L(K(rr,1)+R_row1,jj-1)) + esi*(phi(VR1+V3,1)-phi(VR1+V2,1))*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj));
                                res(ii,1) = res(ii,1) + coeff*(Area(K(rr,1)+R_row1,jj)*L(K(rr,1)+R_row1,jj)+Area(K(rr,1)+R_row1,jj-1)*L(K(rr,1)+R_row1,jj-1))*(-Na-elec(VR1+V2,1)+hole(VR1+V2,1));

                                Jaco(ii, 3*(V3)+20) = Jaco(ii, 3*(V3)+20)+esi*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj));
                                Jaco(ii, 3*(V2)+20) =  Jaco(ii, 3*(V2)+20)-esi*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj))-esi*(Area(K(rr,1)+R_row1,jj-1)/L(K(rr,1)+R_row1,jj-1));
                                Jaco(ii, 3*(V1)+20) = Jaco(ii, 3*(V1)+20)+esi*(Area(K(rr,1)+R_row1,jj-1)/L(K(rr,1)+R_row1,jj-1));

                                Jaco(ii, 3*(V2)+21) = Jaco(ii, 3*(V2)+21) - coeff*(Area(K(rr,1)+R_row1,jj)*L(K(rr,1)+R_row1,jj)+Area(K(rr,1)+R_row1,jj-1)*L(K(rr,1)+R_row1,jj-1));
                                Jaco(ii, 3*(V2)+22) = Jaco(ii, 3*(V2)+22) + coeff*(Area(K(rr,1)+R_row1,jj)*L(K(rr,1)+R_row1,jj)+Area(K(rr,1)+R_row1,jj-1)*L(K(rr,1)+R_row1,jj-1));

                            elseif K(rr,1) > 2*R_row2 && K(rr,1) <= 3*R_row2
                                jj = 3;
                                K(rr,1) = K(rr,1)-2*R_row2;

                                V3 = find(Table(1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj-2)) ;
                                V2 = find(Table(1+1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj-1));
                                V1 = find(Table(1+1+1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj));

                                res(ii,1) = res(ii,1) + esi*(phi(VR1+V3,1)-phi(VR1+V1,1))*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj)) + esi*(phi(VR1+V2,1)-phi(VR1+V1,1))*(Area(K(rr,1)+R_row1,jj-1)/L(K(rr,1)+R_row1,jj-1));
                                res(ii,1) = res(ii,1) + coeff*(Area(K(rr,1)+R_row1,jj)*L(K(rr,1)+R_row1,jj)+Area(K(rr,1)+R_row1,jj-1)*L(K(rr,1)+R_row1,jj-1))*(-Na-elec(VR1+V1,1)+hole(VR1+V1,1));

                                Jaco(ii, 3*(V3)+20) = Jaco(ii, 3*(V3)+20)+esi*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj));
                                Jaco(ii, 3*(V2)+20) = Jaco(ii, 3*(V2)+20)+esi*(Area(K(rr,1)+R_row1,jj-1)/L(K(rr,1)+R_row1,jj-1));
                                Jaco(ii, 3*(V1)+20) = Jaco(ii, 3*(V1)+20)-esi*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj))-esi*(Area(K(rr,1)+R_row1,jj-1)/L(K(rr,1)+R_row1,jj-1));

                                Jaco(ii, 3*(V1)+21) = Jaco(ii, 3*(V1)+21) - coeff*(Area(K(rr,1)+R_row1,jj)*L(K(rr,1)+R_row1,jj)+Area(K(rr,1)+R_row1,jj-1)*L(K(rr,1)+R_row1,jj-1));
                                Jaco(ii, 3*(V1)+22) = Jaco(ii, 3*(V1)+22) + coeff*(Area(K(rr,1)+R_row1,jj)*L(K(rr,1)+R_row1,jj)+Area(K(rr,1)+R_row1,jj-1)*L(K(rr,1)+R_row1,jj-1));
                            end
                        end

                    elseif Table(ii,3) == eDensity

                        for rr = 1:K_row
                            if K(rr,1) <= R_row2
                                jj = 1;

                                V3 = find(Table(1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj+2)) ;  % node 3
                                V2 = find(Table(1+1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj+1)); % node 2
                                V1 = find(Table(2+1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj));   % node 1

                                x21 = (phi(VR1+V2,1)-phi(VR1+V1,1))/Thermal_V;
                                x12 = -x21;
                                x31 = (phi(VR1+V3,1)-phi(VR1+V1,1))/Thermal_V;
                                x13 = -x31;

                                res(ii,1) = res(ii,1) + coeff_Jn*((Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj))*(elec(VR1+V2,1)*Ber(x21) - elec(VR1+V1,1)*Ber(x12)));
                                res(ii,1) = res(ii,1) + coeff_Jn*((Area(K(rr,1)+R_row1,jj+2)/L(K(rr,1)+R_row1,jj+2))*(elec(VR1+V3,1)*Ber(x31) - elec(VR1+V1,1)*Ber(x13)));

                                % eDensity
                                Jaco(ii, 3*(V1)+21) = Jaco(ii, 3*(V1)+21) + coeff_Jn*((Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj))*(-Ber(x12)));
                                Jaco(ii, 3*(V1)+21) = Jaco(ii, 3*(V1)+21) + coeff_Jn*((Area(K(rr,1)+R_row1,jj+2)/L(K(rr,1)+R_row1,jj+2))*(-Ber(x13)));

                                Jaco(ii, 3*(V2)+21) = Jaco(ii, 3*(V2)+21) + coeff_Jn*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj))*Ber(x21);
                                Jaco(ii, 3*(V3)+21) = Jaco(ii, 3*(V3)+21) + coeff_Jn*(Area(K(rr,1)+R_row1,jj+2)/L(K(rr,1)+R_row1,jj+2))*Ber(x31);

                                % Potential
                                Jaco(ii, 3*(V1)+20) = Jaco(ii, 3*(V1)+20) + coeff_Jn/Thermal_V*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj))*(-elec(VR1+V2,1)*Ber_d(x21)-elec(VR1+V1,1)*Ber_d(x12));
                                Jaco(ii, 3*(V1)+20) = Jaco(ii, 3*(V1)+20) + coeff_Jn/Thermal_V*(Area(K(rr,1)+R_row1,jj+2)/L(K(rr,1)+R_row1,jj+2))*(-elec(VR1+V3,1)*Ber_d(x31)-elec(VR1+V1,1)*Ber_d(x13));

                                Jaco(ii, 3*(V2)+20) = Jaco(ii, 3*(V2)+20) + coeff_Jn/Thermal_V*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj))*(elec(VR1+V2,1)*Ber_d(x21)+elec(VR1+V1,1)*Ber_d(x12));
                                Jaco(ii, 3*(V3)+20) = Jaco(ii, 3*(V3)+20) + coeff_Jn/Thermal_V*(Area(K(rr,1)+R_row1,jj+2)/L(K(rr,1)+R_row1,jj+2))*(elec(VR1+V3,1)*Ber_d(x31)+elec(VR1+V1,1)*Ber_d(x13));


                            elseif K(rr,1) > R_row2 && K(rr,1) <= 2*R_row2
                                jj = 2;
                                K(rr,1) = K(rr,1)-R_row2;

                                V3 = find(Table(1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj+1)) ; % node 3
                                V2 = find(Table(1+1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj));  % node 2
                                V1 = find(Table(2+1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj-1)); % node 1

                                x12 = (phi(VR1+V1,1)-phi(VR1+V2,1))/Thermal_V;
                                x21 = -x12;
                                x32 = (phi(VR1+V3,1)-phi(VR1+V2,1))/Thermal_V;
                                x23 = -x32;

                                res(ii,1) =  res(ii,1) + coeff_Jn*((Area(K(rr,1)+R_row1,jj-1)/L(K(rr,1)+R_row1,jj-1))*(elec(VR1+V1,1)*Ber(x12) - elec(VR1+V2,1)*Ber(x21)));
                                res(ii,1) =  res(ii,1) + coeff_Jn*((Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj))*(elec(VR1+V3,1)*Ber(x32) - elec(VR1+V2,1)*Ber(x23)));

                                % eDensity
                                Jaco(ii, 3*(V2)+21) = Jaco(ii, 3*(V2)+21) + coeff_Jn*((Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj))*(-Ber(x23)));
                                Jaco(ii, 3*(V2)+21) = Jaco(ii, 3*(V2)+21) + coeff_Jn*((Area(K(rr,1)+R_row1,jj-1)/L(K(rr,1)+R_row1,jj-1))*(-Ber(x21)));

                                Jaco(ii, 3*(V1)+21) = Jaco(ii, 3*(V1)+21) + coeff_Jn*(Area(K(rr,1)+R_row1,jj-1)/L(K(rr,1)+R_row1,jj-1))*Ber(x12);
                                Jaco(ii, 3*(V3)+21) = Jaco(ii, 3*(V3)+21) + coeff_Jn*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj))*Ber(x32);

                                % Potential
                                Jaco(ii, 3*(V2)+20) = Jaco(ii, 3*(V2)+20) + coeff_Jn/Thermal_V*(Area(K(rr,1)+R_row1,jj-1)/L(K(rr,1)+R_row1,jj-1))*(-elec(VR1+V1,1)*Ber_d(x12)-elec(VR1+V2,1)*Ber_d(x21));
                                Jaco(ii, 3*(V2)+20) = Jaco(ii, 3*(V2)+20) + coeff_Jn/Thermal_V*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj))*(-elec(VR1+V3,1)*Ber_d(x32)-elec(VR1+V2,1)*Ber_d(x23));

                                Jaco(ii, 3*(V1)+20) = Jaco(ii, 3*(V1)+20) + coeff_Jn/Thermal_V*(Area(K(rr,1)+R_row1,jj-1)/L(K(rr,1)+R_row1,jj-1))*(elec(VR1+V1,1)*Ber_d(x12)+elec(VR1+V2,1)*Ber_d(x21));
                                Jaco(ii, 3*(V3)+20) = Jaco(ii, 3*(V3)+20) + coeff_Jn/Thermal_V*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj))*(elec(VR1+V3,1)*Ber_d(x32)+elec(VR1+V2,1)*Ber_d(x23));

                            elseif K(rr,1) > 2*R_row2 && K(rr,1) <= 3*R_row2
                                jj = 3;
                                K(rr,1) = K(rr,1)-2*R_row2;

                                V3 = find(Table(2+1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj)); % node 3
                                V2 = find(Table(1+1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj-1)); % node 2
                                V1 = find(Table(1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj-2)) ;  % node 1

                                x13 = (phi(VR1+V1,1)-phi(VR1+V3,1))/Thermal_V;
                                x31 = -x13;
                                x23 = (phi(VR1+V2,1)-phi(VR1+V3,1))/Thermal_V;
                                x32 = -x23;

                                res(ii,1) = res(ii,1) + coeff_Jn*((Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj))*(elec(VR1+V1,1)*Ber(x13) - elec(VR1+V3,1)*Ber(x31)));
                                res(ii,1) = res(ii,1) + coeff_Jn*((Area(K(rr,1)+R_row1,jj-1)/L(K(rr,1)+R_row1,jj-1))*(elec(VR1+V2,1)*Ber(x23) - elec(VR1+V3,1)*Ber(x32)));

                                % eDensity
                                Jaco(ii, 3*(V1)+21) = Jaco(ii, 3*(V1)+21) + coeff_Jn*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj))*Ber(x13);
                                Jaco(ii, 3*(V2)+21) = Jaco(ii, 3*(V2)+21) + coeff_Jn*(Area(K(rr,1)+R_row1,jj-1)/L(K(rr,1)+R_row1,jj-1))*Ber(x23);

                                Jaco(ii, 3*(V3)+21) = Jaco(ii, 3*(V3)+21) + coeff_Jn*((Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj))*(-Ber(x31)));
                                Jaco(ii, 3*(V3)+21) = Jaco(ii, 3*(V3)+21) + coeff_Jn*((Area(K(rr,1)+R_row1,jj-1)/L(K(rr,1)+R_row1,jj-1))*(-Ber(x32)));

                                % Potential
                                Jaco(ii, 3*(V3)+20) = Jaco(ii, 3*(V3)+20) + coeff_Jn*(Area(K(rr,1)+R_row1,jj-1)/L(K(rr,1)+R_row1,jj-1))*(-elec(VR1+V2,1)*Ber_d(x23)-elec(VR1+V3,1)*Ber_d(x32))/Thermal_V;
                                Jaco(ii, 3*(V3)+20) = Jaco(ii, 3*(V3)+20) + coeff_Jn*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj))*(-elec(VR1+V1,1)*Ber_d(x13)-elec(VR1+V3,1)*Ber_d(x31))/Thermal_V;

                                Jaco(ii, 3*(V2)+20) = Jaco(ii, 3*(V2)+20) + coeff_Jn*(Area(K(rr,1)+R_row1,jj-1)/L(K(rr,1)+R_row1,jj-1))*(elec(VR1+V2,1)*Ber_d(x23)+elec(VR1+V3,1)*Ber_d(x32))/Thermal_V;
                                Jaco(ii, 3*(V1)+20) = Jaco(ii, 3*(V1)+20) + coeff_Jn*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj))*(elec(VR1+V1,1)*Ber_d(x13)+elec(VR1+V3,1)*Ber_d(x31))/Thermal_V;
                            end
                        end


                    elseif Table(ii,3) == hDensity

                        for rr = 1:K_row

                            if K(rr,1) <= R_row2
                                jj = 1;

                                V3 = find(Table(1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj+2)) ;  % node 3
                                V2 = find(Table(1+1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj+1));  % node 2
                                V1 = find(Table(2+1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj));  % node 1

                                x21 = (phi(VR1+V2,1)-phi(VR1+V1,1))/Thermal_V;
                                x12 = -x21;
                                x31 = (phi(VR1+V3,1)-phi(VR1+V1,1))/Thermal_V;
                                x13 = -x31;

                                res(ii,1) = res(ii,1) + coeff_Jp*((Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj))*(hole(VR1+V2,1)*Ber(x12) - hole(VR1+V1,1)*Ber(x21)));
                                res(ii,1) = res(ii,1) + coeff_Jp*((Area(K(rr,1)+R_row1,jj+2)/L(K(rr,1)+R_row1,jj+2))*(hole(VR1+V3,1)*Ber(x13) - hole(VR1+V1,1)*Ber(x31)));

                                %hDensity
                                Jaco(ii, 3*(V1)+22) = Jaco(ii, 3*(V1)+22) + coeff_Jp*((Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj))*(-Ber(x21))+(Area(K(rr,1)+R_row1,jj+2)/L(K(rr,1)+R_row1,jj+2))*(-Ber(x31)));
                                Jaco(ii, 3*(V2)+22) = Jaco(ii, 3*(V2)+22) + coeff_Jp*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj))*Ber(x12);
                                Jaco(ii, 3*(V3)+22) = Jaco(ii, 3*(V3)+22) + coeff_Jp*(Area(K(rr,1)+R_row1,jj+2)/L(K(rr,1)+R_row1,jj+2))*Ber(x13);

                                %Potential
                                Jaco(ii, 3*(V1)+20) = Jaco(ii, 3*(V1)+20) + coeff_Jp*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj))*(hole(VR1+V2,1)*Ber_d(x12)+hole(VR1+V1,1)*Ber_d(x21))/Thermal_V;
                                Jaco(ii, 3*(V1)+20) = Jaco(ii, 3*(V1)+20) + coeff_Jp*(Area(K(rr,1)+R_row1,jj+2)/L(K(rr,1)+R_row1,jj+2))*(hole(VR1+V3,1)*Ber_d(x13)+hole(VR1+V1,1)*Ber_d(x31))/Thermal_V;
                                Jaco(ii, 3*(V2)+20) = Jaco(ii, 3*(V2)+20) + coeff_Jp*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj))*(-hole(VR1+V2,1)*Ber_d(x12)-hole(VR1+V1,1)*Ber_d(x21))/Thermal_V;
                                Jaco(ii, 3*(V3)+20) = Jaco(ii, 3*(V3)+20) + coeff_Jp*(Area(K(rr,1)+R_row1,jj+2)/L(K(rr,1)+R_row1,jj+2))*(-hole(VR1+V3,1)*Ber_d(x13)-hole(VR1+V1,1)*Ber_d(x31))/Thermal_V;


                            elseif K(rr,1) > R_row2 && K(rr,1) <= 2*R_row2
                                jj = 2;
                                K(rr,1) = K(rr,1)-R_row2;

                                V3 = find(Table(1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj+1)) ;  % node 3
                                V2 = find(Table(1+1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj));  % node 2
                                V1 = find(Table(2+1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj-1));  % node 1

                                x12 = (phi(VR1+V1,1)-phi(VR1+V2,1))/Thermal_V;
                                x21 = -x12;
                                x32 = (phi(VR1+V3,1)-phi(VR1+V2,1))/Thermal_V;
                                x23 = -x32;

                                res(ii,1) =  res(ii,1) + coeff_Jp*((Area(K(rr,1)+R_row1,jj-1)/L(K(rr,1)+R_row1,jj-1))*(hole(VR1+V1,1)*Ber(x21) - hole(VR1+V2,1)*Ber(x12)));
                                res(ii,1) =  res(ii,1) + coeff_Jp*((Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj))*(hole(VR1+V3,1)*Ber(x23) - hole(VR1+V2,1)*Ber(x32)));

                                % hDensity
                                Jaco(ii, 3*(V1)+22) = Jaco(ii, 3*(V1)+22) + coeff_Jp*(Area(K(rr,1)+R_row1,jj-1)/L(K(rr,1)+R_row1,jj-1))*Ber(x21);
                                Jaco(ii, 3*(V3)+22) = Jaco(ii, 3*(V3)+22) + coeff_Jp*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj))*Ber(x23);

                                Jaco(ii, 3*(V2)+22) = Jaco(ii, 3*(V2)+22) + coeff_Jp*((Area(K(rr,1)+R_row1,jj-1)/L(K(rr,1)+R_row1,jj-1))*(-Ber(x12)));
                                Jaco(ii, 3*(V2)+22) = Jaco(ii, 3*(V2)+22) + coeff_Jp*((Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj))*(-Ber(x32)));


                                %Potential
                                Jaco(ii, 3*(V2)+20) = Jaco(ii, 3*(V2)+20) + coeff_Jp*(Area(K(rr,1)+R_row1,jj-1)/L(K(rr,1)+R_row1,jj-1))*(hole(VR1+V1,1)*Ber_d(x21)+hole(VR1+V2,1)*Ber_d(x12))/Thermal_V;
                                Jaco(ii, 3*(V2)+20) = Jaco(ii, 3*(V2)+20) + coeff_Jp*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj))*(hole(VR1+V3,1)*Ber_d(x23)+hole(VR1+V2,1)*Ber_d(x32))/Thermal_V;

                                Jaco(ii, 3*(V1)+20) = Jaco(ii, 3*(V1)+20) + coeff_Jp*(Area(K(rr,1)+R_row1,jj-1)/L(K(rr,1)+R_row1,jj-1))*(-hole(VR1+V1,1)*Ber_d(x21)-hole(VR1+V2,1)*Ber_d(x12))/Thermal_V;
                                Jaco(ii, 3*(V3)+20) = Jaco(ii, 3*(V3)+20) + coeff_Jp*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj))*(-hole(VR1+V3,1)*Ber_d(x23)-hole(VR1+V2,1)*Ber_d(x32))/Thermal_V;


                            elseif K(rr,1) > 2*R_row2 && K(rr,1) <= 3*R_row2
                                jj = 3;
                                K(rr,1) = K(rr,1)-2*R_row2;

                                V3 = find(Table(2+1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj));  % node 3
                                V2 = find(Table(1+1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj-1));  % node 2
                                V1 = find(Table(1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj-2)) ;  % node 1

                                x13 = (phi(VR1+V1,1)-phi(VR1+V3,1))/Thermal_V;
                                x31 = -x13;
                                x23 = (phi(VR1+V2,1)-phi(VR1+V3,1))/Thermal_V;
                                x32 = -x23;

                                res(ii,1) = res(ii,1) + coeff_Jp*((Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj))*(hole(VR1+V1,1)*Ber(x31) - hole(VR1+V3,1)*Ber(x13)));
                                res(ii,1) = res(ii,1) + coeff_Jp*((Area(K(rr,1)+R_row1,jj-1)/L(K(rr,1)+R_row1,jj-1))*(hole(VR1+V2,1)*Ber(x32) - hole(VR1+V3,1)*Ber(x23)));

                                %hDensity
                                Jaco(ii, 3*(V1)+22) = Jaco(ii, 3*(V1)+22) + coeff_Jp*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj))*Ber(x31);
                                Jaco(ii, 3*(V2)+22) = Jaco(ii, 3*(V2)+22) + coeff_Jp*(Area(K(rr,1)+R_row1,jj-1)/L(K(rr,1)+R_row1,jj-1))*Ber(x32);

                                Jaco(ii, 3*(V3)+22) = Jaco(ii, 3*(V3)+22) + coeff_Jp*((Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj))*(-Ber(x13)));
                                Jaco(ii, 3*(V3)+22) = Jaco(ii, 3*(V3)+22) + coeff_Jp*((Area(K(rr,1)+R_row1,jj-1)/L(K(rr,1)+R_row1,jj-1))*(-Ber(x23)));

                                %Potential
                                Jaco(ii, 3*(V3)+20) = Jaco(ii, 3*(V3)+20) + coeff_Jp*(Area(K(rr,1)+R_row1,jj-1)/L(K(rr,1)+R_row1,jj-1))*(hole(VR1+V2,1)*Ber_d(x32)+hole(VR1+V3,1)*Ber_d(x23))/Thermal_V;
                                Jaco(ii, 3*(V3)+20) = Jaco(ii, 3*(V3)+20) + coeff_Jp*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj))*(hole(VR1+V1,1)*Ber_d(x31)+hole(VR1+V3,1)*Ber_d(x13))/Thermal_V;
                                Jaco(ii, 3*(V2)+20) = Jaco(ii, 3*(V2)+20) + coeff_Jp*(Area(K(rr,1)+R_row1,jj-1)/L(K(rr,1)+R_row1,jj-1))*(-hole(VR1+V2,1)*Ber_d(x32)-hole(VR1+V3,1)*Ber_d(x23))/Thermal_V;
                                Jaco(ii, 3*(V1)+20) = Jaco(ii, 3*(V1)+20) + coeff_Jp*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj))*(-hole(VR1+V1,1)*Ber_d(x31)-hole(VR1+V3,1)*Ber_d(x13))/Thermal_V;
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

                                Jaco(ii, 3*(V3)+20) = Jaco(ii, 3*(V3)+20)+ esi*(Area(K(rr,1)+R_row1,jj+2)/L(K(rr,1)+R_row1,jj+2));
                                Jaco(ii, 3*(V2)+20) = Jaco(ii, 3*(V2)+20)+esi*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj));
                                Jaco(ii, 3*(V1)+20) = Jaco(ii, 3*(V1)+20)-esi*(Area(K(rr,1)+R_row1,jj+2)/L(K(rr,1)+R_row1,jj+2))-esi*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj));

                                Jaco(ii, 3*(V1)+21) = Jaco(ii, 3*(V1)+21) - coeff*(Area(K(rr,1)+R_row1,jj)*L(K(rr,1)+R_row1,jj)+Area(K(rr,1)+R_row1,jj+2)*L(K(rr,1)+R_row1,jj+2));
                                Jaco(ii, 3*(V1)+22) = Jaco(ii, 3*(V1)+22) + coeff*(Area(K(rr,1)+R_row1,jj)*L(K(rr,1)+R_row1,jj)+Area(K(rr,1)+R_row1,jj+2)*L(K(rr,1)+R_row1,jj+2));

                            elseif K(rr,1) > R_row2 && K(rr,1) <= 2*R_row2
                                jj = 2;
                                K(rr,1) = K(rr,1)-R_row2;

                                V3 = find(Table(1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj+1)) ;
                                V2 = find(Table(1+1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj));
                                V1 = find(Table(1+1+1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj-1));

                                res(ii,1) = res(ii,1) + esi*(phi(VR1+V1,1)-phi(VR1+V2,1))*(Area(K(rr,1)+R_row1,jj-1)/L(K(rr,1)+R_row1,jj-1)) + esi*(phi(VR1+V3,1)-phi(VR1+V2,1))*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj));
                                res(ii,1) = res(ii,1) + coeff*(Area(K(rr,1)+R_row1,jj)*L(K(rr,1)+R_row1,jj)+Area(K(rr,1)+R_row1,jj-1)*L(K(rr,1)+R_row1,jj-1))*(Nd-elec(VR1+V2,1)+hole(VR1+V2,1));

                                Jaco(ii, 3*(V3)+20) = Jaco(ii, 3*(V3)+20)+esi*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj));
                                Jaco(ii, 3*(V2)+20) =  Jaco(ii, 3*(V2)+20)-esi*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj))-esi*(Area(K(rr,1)+R_row1,jj-1)/L(K(rr,1)+R_row1,jj-1));
                                Jaco(ii, 3*(V1)+20) = Jaco(ii, 3*(V1)+20)+esi*(Area(K(rr,1)+R_row1,jj-1)/L(K(rr,1)+R_row1,jj-1));

                                Jaco(ii, 3*(V2)+21) = Jaco(ii, 3*(V2)+21) - coeff*(Area(K(rr,1)+R_row1,jj)*L(K(rr,1)+R_row1,jj)+Area(K(rr,1)+R_row1,jj-1)*L(K(rr,1)+R_row1,jj-1));
                                Jaco(ii, 3*(V2)+22) = Jaco(ii, 3*(V2)+22) + coeff*(Area(K(rr,1)+R_row1,jj)*L(K(rr,1)+R_row1,jj)+Area(K(rr,1)+R_row1,jj-1)*L(K(rr,1)+R_row1,jj-1));

                            elseif K(rr,1) > 2*R_row2 && K(rr,1) <= 3*R_row2
                                jj = 3;
                                K(rr,1) = K(rr,1)-2*R_row2;

                                V3 = find(Table(1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj-2)) ;  % node 1
                                V2 = find(Table(1+1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj-1));  % node 1
                                V1 = find(Table(1+1+1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj));  % node 1

                                res(ii,1) = res(ii,1) + esi*(phi(VR1+V3,1)-phi(VR1+V1,1))*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj)) + esi*(phi(VR1+V2,1)-phi(VR1+V1,1))*(Area(K(rr,1)+R_row1,jj-1)/L(K(rr,1)+R_row1,jj-1));
                                res(ii,1) = res(ii,1) + coeff*(Area(K(rr,1)+R_row1,jj)*L(K(rr,1)+R_row1,jj)+Area(K(rr,1)+R_row1,jj-1)*L(K(rr,1)+R_row1,jj-1))*(Nd-elec(VR1+V1,1)+hole(VR1+V1,1));

                                Jaco(ii, 3*(V3)+20) = Jaco(ii, 3*(V3)+20)+esi*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj));
                                Jaco(ii, 3*(V2)+20) = Jaco(ii, 3*(V2)+20)+esi*(Area(K(rr,1)+R_row1,jj-1)/L(K(rr,1)+R_row1,jj-1));
                                Jaco(ii, 3*(V1)+20) = Jaco(ii, 3*(V1)+20)-esi*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj))-esi*(Area(K(rr,1)+R_row1,jj-1)/L(K(rr,1)+R_row1,jj-1));

                                Jaco(ii, 3*(V1)+21) = Jaco(ii, 3*(V1)+21) - coeff*(Area(K(rr,1)+R_row1,jj)*L(K(rr,1)+R_row1,jj)+Area(K(rr,1)+R_row1,jj-1)*L(K(rr,1)+R_row1,jj-1));
                                Jaco(ii, 3*(V1)+22) = Jaco(ii, 3*(V1)+22) + coeff*(Area(K(rr,1)+R_row1,jj)*L(K(rr,1)+R_row1,jj)+Area(K(rr,1)+R_row1,jj-1)*L(K(rr,1)+R_row1,jj-1));
                            end
                        end

                    elseif Table(ii,3) == eDensity

                        for rr = 1:K_row
                            if K(rr,1) <= R_row2
                                jj = 1;

                                V3 = find(Table(1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj+2)) ; % node 3
                                V2 = find(Table(1+1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj+1)); % node 2
                                V1 = find(Table(2+1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj)); %node 1

                                x21 = (phi(VR1+V2,1)-phi(VR1+V1,1))/Thermal_V;
                                x12 = -x21;
                                x31 = (phi(VR1+V3,1)-phi(VR1+V1,1))/Thermal_V;
                                x13 = -x31;

                                res(ii,1) = res(ii,1) + coeff_Jn*((Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj))*(elec(VR1+V2,1)*Ber(x21) - elec(VR1+V1,1)*Ber(x12)));
                                res(ii,1) = res(ii,1) + coeff_Jn*((Area(K(rr,1)+R_row1,jj+2)/L(K(rr,1)+R_row1,jj+2))*(elec(VR1+V3,1)*Ber(x31) - elec(VR1+V1,1)*Ber(x13)));

                                % eDensity
                                Jaco(ii, 3*(V1)+21) = Jaco(ii, 3*(V1)+21) + coeff_Jn*((Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj))*(-Ber(x12))+(Area(K(rr,1)+R_row1,jj+2)/L(K(rr,1)+R_row1,jj+2))*(-Ber(x13)));
                                Jaco(ii, 3*(V2)+21) = Jaco(ii, 3*(V2)+21) + coeff_Jn*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj))*Ber(x21);
                                Jaco(ii, 3*(V3)+21) = Jaco(ii, 3*(V3)+21) + coeff_Jn*(Area(K(rr,1)+R_row1,jj+2)/L(K(rr,1)+R_row1,jj+2))*Ber(x31);

                                % Potential
                                Jaco(ii, 3*(V1)+20) = Jaco(ii, 3*(V1)+20) + coeff_Jn*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj))*(-elec(VR1+V2,1)*Ber_d(x21)-elec(VR1+V1,1)*Ber_d(x12))/Thermal_V;
                                Jaco(ii, 3*(V1)+20) = Jaco(ii, 3*(V1)+20) + coeff_Jn*(Area(K(rr,1)+R_row1,jj+2)/L(K(rr,1)+R_row1,jj+2))*(-elec(VR1+V3,1)*Ber_d(x31)-elec(VR1+V1,1)*Ber_d(x13))/Thermal_V;
                                Jaco(ii, 3*(V2)+20) = Jaco(ii, 3*(V2)+20) + coeff_Jn*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj))*(elec(VR1+V2,1)*Ber_d(x21)+elec(VR1+V1,1)*Ber_d(x12))/Thermal_V;
                                Jaco(ii, 3*(V3)+20) = Jaco(ii, 3*(V3)+20) + coeff_Jn*(Area(K(rr,1)+R_row1,jj+2)/L(K(rr,1)+R_row1,jj+2))*(elec(VR1+V3,1)*Ber_d(x31)+elec(VR1+V1,1)*Ber_d(x13))/Thermal_V;

                            elseif K(rr,1) > R_row2 && K(rr,1) <= 2*R_row2
                                jj = 2;
                                K(rr,1) = K(rr,1)-R_row2;

                                V3 = find(Table(1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj+1)) ; % node 3
                                V2 = find(Table(1+1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj)); % node 2
                                V1 = find(Table(1+1+1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj-1)); % node 1

                                x12 = (phi(VR1+V1,1)-phi(VR1+V2,1))/Thermal_V;
                                x21 = -x12;
                                x32 = (phi(VR1+V3,1)-phi(VR1+V2,1))/Thermal_V;
                                x23 = -x32;

                                res(ii,1) =  res(ii,1) + coeff_Jn*((Area(K(rr,1)+R_row1,jj-1)/L(K(rr,1)+R_row1,jj-1))*(elec(VR1+V1,1)*Ber(x12) - elec(VR1+V2,1)*Ber(x21)));
                                res(ii,1) =  res(ii,1) + coeff_Jn*((Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj))*(elec(VR1+V3,1)*Ber(x32) - elec(VR1+V2,1)*Ber(x23)));

                                % eDensity
                                Jaco(ii, 3*(V1)+21) = Jaco(ii, 3*(V1)+21) + coeff_Jn*(Area(K(rr,1)+R_row1,jj-1)/L(K(rr,1)+R_row1,jj-1))*Ber(x12);
                                Jaco(ii, 3*(V2)+21) = Jaco(ii, 3*(V2)+21) + coeff_Jn*((Area(K(rr,1)+R_row1,jj-1)/L(K(rr,1)+R_row1,jj-1))*(-Ber(x21))+(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj))*(-Ber(x23)));
                                Jaco(ii, 3*(V3)+21) = Jaco(ii, 3*(V3)+21) + coeff_Jn*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj))*Ber(x32);

                                % Potential
                                Jaco(ii, 3*(V2)+20) = Jaco(ii, 3*(V2)+20) + coeff_Jn*(Area(K(rr,1)+R_row1,jj-1)/L(K(rr,1)+R_row1,jj-1))*(-elec(VR1+V1,1)*Ber_d(x12)-elec(VR1+V2,1)*Ber_d(x21))/Thermal_V;
                                Jaco(ii, 3*(V2)+20) = Jaco(ii, 3*(V2)+20) + coeff_Jn*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj))*(-elec(VR1+V3,1)*Ber_d(x32)-elec(VR1+V2,1)*Ber_d(x23))/Thermal_V;
                                Jaco(ii, 3*(V1)+20) = Jaco(ii, 3*(V1)+20) + coeff_Jn*(Area(K(rr,1)+R_row1,jj-1)/L(K(rr,1)+R_row1,jj-1))*(elec(VR1+V1,1)*Ber_d(x12)+elec(VR1+V2,1)*Ber_d(x21))/Thermal_V;
                                Jaco(ii, 3*(V3)+20) = Jaco(ii, 3*(V3)+20) + coeff_Jn*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj))*(elec(VR1+V3,1)*Ber_d(x32)+elec(VR1+V2,1)*Ber_d(x23))/Thermal_V;

                            elseif K(rr,1) > 2*R_row2 && K(rr,1) <= 3*R_row2
                                jj = 3;
                                K(rr,1) = K(rr,1)-2*R_row2;

                                V1 = find(Table(1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj-2)) ; % node 1
                                V2 = find(Table(1+1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj-1)); % node 2
                                V3 = find(Table(1+1+1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj)); % node 3

                                x13 = (phi(VR1+V1,1)-phi(VR1+V3,1))/Thermal_V;
                                x31 = -x13;
                                x23 = (phi(VR1+V2,1)-phi(VR1+V3,1))/Thermal_V;
                                x32 = -x23;

                                res(ii,1) = res(ii,1) + coeff_Jn*((Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj))*(elec(VR1+V1,1)*Ber(x13) - elec(VR1+V3,1)*Ber(x31)));
                                res(ii,1) = res(ii,1) + coeff_Jn*((Area(K(rr,1)+R_row1,jj-1)/L(K(rr,1)+R_row1,jj-1))*(elec(VR1+V2,1)*Ber(x23) - elec(VR1+V3,1)*Ber(x32)));

                                % eDensity
                                Jaco(ii, 3*(V1)+21) = Jaco(ii, 3*(V1)+21) + coeff_Jn*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj))*Ber(x13);
                                Jaco(ii, 3*(V2)+21) = Jaco(ii, 3*(V2)+21) + coeff_Jn*(Area(K(rr,1)+R_row1,jj-1)/L(K(rr,1)+R_row1,jj-1))*Ber(x23);
                                Jaco(ii, 3*(V3)+21) = Jaco(ii, 3*(V3)+21) + coeff_Jn*((Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj))*(-Ber(x31))+(Area(K(rr,1)+R_row1,jj-1)/L(K(rr,1)+R_row1,jj-1))*(-Ber(x32)));

                                % Potential
                                Jaco(ii, 3*(V3)+20) = Jaco(ii, 3*(V3)+20) + coeff_Jn*(Area(K(rr,1)+R_row1,jj-1)/L(K(rr,1)+R_row1,jj-1))*(-elec(VR1+V2,1)*Ber_d(x23)-elec(VR1+V3,1)*Ber_d(x32))/Thermal_V;
                                Jaco(ii, 3*(V3)+20) = Jaco(ii, 3*(V3)+20) + coeff_Jn*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj))*(-elec(VR1+V1,1)*Ber_d(x13)-elec(VR1+V3,1)*Ber_d(x31))/Thermal_V;
                                Jaco(ii, 3*(V2)+20) = Jaco(ii, 3*(V2)+20) + coeff_Jn*(Area(K(rr,1)+R_row1,jj-1)/L(K(rr,1)+R_row1,jj-1))*(elec(VR1+V2,1)*Ber_d(x23)+elec(VR1+V3,1)*Ber_d(x32))/Thermal_V;
                                Jaco(ii, 3*(V1)+20) = Jaco(ii, 3*(V1)+20) + coeff_Jn*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj))*(elec(VR1+V1,1)*Ber_d(x13)+elec(VR1+V3,1)*Ber_d(x31))/Thermal_V;

                            end
                        end


                    elseif Table(ii,3) == hDensity

                        for rr = 1:K_row

                            if K(rr,1) <= R_row2
                                jj = 1;

                                V3 = find(Table(1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj+2)) ;  % node 3
                                V2 = find(Table(1+1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj+1));  % node 2
                                V1 = find(Table(2+1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj));  % node 1

                                x21 = (phi(VR1+V2,1)-phi(VR1+V1,1))/Thermal_V;
                                x12 = -x21;
                                x31 = (phi(VR1+V3,1)-phi(VR1+V1,1))/Thermal_V;
                                x13 = -x31;

                                res(ii,1) = res(ii,1) + coeff_Jp*((Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj))*(hole(VR1+V2,1)*Ber(x12) - hole(VR1+V1,1)*Ber(x21)));
                                res(ii,1) = res(ii,1) + coeff_Jp*((Area(K(rr,1)+R_row1,jj+2)/L(K(rr,1)+R_row1,jj+2))*(hole(VR1+V3,1)*Ber(x13) - hole(VR1+V1,1)*Ber(x31)));

                                %hDensity
                                Jaco(ii, 3*(V1)+22) = Jaco(ii, 3*(V1)+22) + coeff_Jp*((Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj))*(-Ber(x21))+(Area(K(rr,1)+R_row1,jj+2)/L(K(rr,1)+R_row1,jj+2))*(-Ber(x31)));
                                Jaco(ii, 3*(V2)+22) = Jaco(ii, 3*(V2)+22) + coeff_Jp*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj))*Ber(x12);
                                Jaco(ii, 3*(V3)+22) = Jaco(ii, 3*(V3)+22) + coeff_Jp*(Area(K(rr,1)+R_row1,jj+2)/L(K(rr,1)+R_row1,jj+2))*Ber(x13);

                                %Potential
                                Jaco(ii, 3*(V1)+20) = Jaco(ii, 3*(V1)+20) + coeff_Jp*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj))*(hole(VR1+V2,1)*Ber_d(x12)+hole(VR1+V1,1)*Ber_d(x21))/Thermal_V;
                                Jaco(ii, 3*(V1)+20) = Jaco(ii, 3*(V1)+20) + coeff_Jp*(Area(K(rr,1)+R_row1,jj+2)/L(K(rr,1)+R_row1,jj+2))*(hole(VR1+V3,1)*Ber_d(x13)+hole(VR1+V1,1)*Ber_d(x31))/Thermal_V;
                                Jaco(ii, 3*(V2)+20) = Jaco(ii, 3*(V2)+20) + coeff_Jp*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj))*(-hole(VR1+V2,1)*Ber_d(x12)-hole(VR1+V1,1)*Ber_d(x21))/Thermal_V;
                                Jaco(ii, 3*(V3)+20) = Jaco(ii, 3*(V3)+20) + coeff_Jp*(Area(K(rr,1)+R_row1,jj+2)/L(K(rr,1)+R_row1,jj+2))*(-hole(VR1+V3,1)*Ber_d(x13)-hole(VR1+V1,1)*Ber_d(x31))/Thermal_V;


                            elseif K(rr,1) > R_row2 && K(rr,1) <= 2*R_row2
                                jj = 2;
                                K(rr,1) = K(rr,1)-R_row2;

                                V3 = find(Table(1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj+1)) ;  % node 3
                                V2 = find(Table(1+1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj));  % node 2
                                V1 = find(Table(1+1+1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj-1));  % node 1

                                x12 = (phi(VR1+V1,1)-phi(VR1+V2,1))/Thermal_V;
                                x21 = -x12;
                                x32 = (phi(VR1+V3,1)-phi(VR1+V2,1))/Thermal_V;
                                x23 = -x32;

                                res(ii,1) =  res(ii,1) + coeff_Jp*((Area(K(rr,1)+R_row1,jj-1)/L(K(rr,1)+R_row1,jj-1))*(hole(VR1+V1,1)*Ber(x21) - hole(VR1+V2,1)*Ber(x12)));
                                res(ii,1) =  res(ii,1) + coeff_Jp*((Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj))*(hole(VR1+V3,1)*Ber(x23) - hole(VR1+V2,1)*Ber(x32)));

                                % hDensity
                                Jaco(ii, 3*(V1)+22) = Jaco(ii, 3*(V1)+22) + coeff_Jp*(Area(K(rr,1)+R_row1,jj-1)/L(K(rr,1)+R_row1,jj-1))*Ber(x21);
                                Jaco(ii, 3*(V2)+22) = Jaco(ii, 3*(V2)+22) + coeff_Jp*((Area(K(rr,1)+R_row1,jj-1)/L(K(rr,1)+R_row1,jj-1))*(-Ber(x12))+(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj))*(-Ber(x32)));
                                Jaco(ii, 3*(V3)+22) = Jaco(ii, 3*(V3)+22) + coeff_Jp*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj))*Ber(x23);

                                %Potential
                                Jaco(ii, 3*(V2)+20) = Jaco(ii, 3*(V2)+20) + coeff_Jp*(Area(K(rr,1)+R_row1,jj-1)/L(K(rr,1)+R_row1,jj-1))*(hole(VR1+V1,1)*Ber_d(x21)+hole(VR1+V2,1)*Ber_d(x12))/Thermal_V;
                                Jaco(ii, 3*(V2)+20) = Jaco(ii, 3*(V2)+20) + coeff_Jp*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj))*(hole(VR1+V3,1)*Ber_d(x23)+hole(VR1+V2,1)*Ber_d(x32))/Thermal_V;
                                Jaco(ii, 3*(V1)+20) = Jaco(ii, 3*(V1)+20) + coeff_Jp*(Area(K(rr,1)+R_row1,jj-1)/L(K(rr,1)+R_row1,jj-1))*(-hole(VR1+V1,1)*Ber_d(x21)-hole(VR1+V2,1)*Ber_d(x12))/Thermal_V;
                                Jaco(ii, 3*(V3)+20) = Jaco(ii, 3*(V3)+20) + coeff_Jp*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj))*(-hole(VR1+V3,1)*Ber_d(x23)-hole(VR1+V2,1)*Ber_d(x32))/Thermal_V;



                            elseif K(rr,1) > 2*R_row2 && K(rr,1) <= 3*R_row2
                                jj = 3;
                                K(rr,1) = K(rr,1)-2*R_row2;

                                V1 = find(Table(1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj-2)) ;  % node 1
                                V2 = find(Table(1+1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj-1));  % node 2
                                V3 = find(Table(1+1+1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj));  % node 3

                                x13 = (phi(VR1+V1,1)-phi(VR1+V3,1))/Thermal_V;
                                x31 = -x13;
                                x23 = (phi(VR1+V2,1)-phi(VR1+V3,1))/Thermal_V;
                                x32 = -x23;

                                res(ii,1) = res(ii,1) + coeff_Jp*((Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj))*(hole(VR1+V1,1)*Ber(x31) - hole(VR1+V3,1)*Ber(x13)));
                                res(ii,1) = res(ii,1) + coeff_Jp*((Area(K(rr,1)+R_row1,jj-1)/L(K(rr,1)+R_row1,jj-1))*(hole(VR1+V2,1)*Ber(x32) - hole(VR1+V3,1)*Ber(x23)));

                                %hDensity
                                Jaco(ii, 3*(V1)+22) = Jaco(ii, 3*(V1)+22) + coeff_Jp*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj))*Ber(x31);
                                Jaco(ii, 3*(V2)+22) = Jaco(ii, 3*(V2)+22) + coeff_Jp*(Area(K(rr,1)+R_row1,jj-1)/L(K(rr,1)+R_row1,jj-1))*Ber(x32);
                                Jaco(ii, 3*(V3)+22) = Jaco(ii, 3*(V3)+22) + coeff_Jp*((Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj))*(-Ber(x13))+(Area(K(rr,1)+R_row1,jj-1)/L(K(rr,1)+R_row1,jj-1))*(-Ber(x23)));

                                %Potential
                                Jaco(ii, 3*(V3)+20) = Jaco(ii, 3*(V3)+20) +coeff_Jp*(Area(K(rr,1)+R_row1,jj-1)/L(K(rr,1)+R_row1,jj-1))*(hole(VR1+V2,1)*Ber_d(x32)+hole(VR1+V3,1)*Ber_d(x23))/Thermal_V;
                                Jaco(ii, 3*(V3)+20) = Jaco(ii, 3*(V3)+20) +coeff_Jp*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj))*(hole(VR1+V1,1)*Ber_d(x31)+hole(VR1+V3,1)*Ber_d(x13))/Thermal_V;
                                Jaco(ii, 3*(V2)+20) = Jaco(ii, 3*(V2)+20) +coeff_Jp*(Area(K(rr,1)+R_row1,jj-1)/L(K(rr,1)+R_row1,jj-1))*(-hole(VR1+V2,1)*Ber_d(x32)-hole(VR1+V3,1)*Ber_d(x23))/Thermal_V;
                                Jaco(ii, 3*(V1)+20) = Jaco(ii, 3*(V1)+20) +coeff_Jp*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj))*(-hole(VR1+V1,1)*Ber_d(x31)-hole(VR1+V3,1)*Ber_d(x13))/Thermal_V;

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

        Jaco(3*Vi21+20,:)=Jaco(3*Vi21+20,:)+Jaco(Vi11,:);
        Jaco(Vi11,:)=0; Jaco(Vi11,Vi11)=1; Jaco(Vi11,3*Vi21+20)=-1;

        res(3*Vi21+20,1)=res(3*Vi21+20,1)+res(Vi11,1);
        res(Vi11,1)=0;
    end

    for rr = 1:i_row2

        Vi12 = find(Table(1+VR1+3*VR2:T1,2)==interface2(rr,1));
        Vi22 = find(Table(1+VR1:3:VR1+3*VR2,2)==interface2(rr,1));

        Jaco(3*Vi22+20,:)=Jaco(3*Vi22+20,:)+Jaco(VR1+3*VR2+Vi12,:);
        Jaco(VR1+3*VR2+Vi12,:)=0; Jaco(VR1+3*VR2+Vi12,VR1+3*VR2+Vi12)=1; Jaco(VR1+3*VR2+Vi12,3*Vi22+20)=-1;

        res(3*Vi22+20,1)=res(3*Vi22+20,1)+res(VR1+3*VR2+Vi12,1);
        res(VR1+3*VR2+Vi12,1)=0;

    end

    %%%% Dirichlet Boundary Condition %%%%%

    for n=1:C_col
        for i =1: C_row
            if i <= 4 % Bottom Gate
                Jaco(Contact(i,n),:) =0;
                Jaco(Contact(i,n),Contact(i,n)) = 1;
                res(Contact(i,n),1)=phi(Contact(i,n),1)-0.33374;

            elseif i > 4 && i <= 8 % Top Gate
                Jaco(110+Contact(i,n),:) =0;
                Jaco(110+Contact(i,n),110+Contact(i,n)) = 1;
                res(110+Contact(i,n),1)= phi(22+Contact(i,n),1)-0.33374;

            end
        end
    end

    %%%% Source, Drain %%%%
    q_row = size(ry1,1);
    w_row = size(ry2,1);

    for rr = 1:q_row  % Source

        V = find(Table(1+VR1:3:VR1+3*VR2,2)==ry1(rr,1));
        Jaco(3*V+20,:) = 0;
        Jaco(3*V+20,3*V+20) = 1;
        res(3*V+20,1) = phi(VR1+V,1)-0;

    end

    for rr = 1:w_row  % Drain

        V = find(Table(1+VR1:3:VR1+3*VR2,2)==ry2(rr,1));
        Jaco(3*V+20,:) = 0;
        Jaco(3*V+20,3*V+20) = 1;
        res(3*V+20,1) = phi(VR1+V,1)-0;

    end

    %     scaling
    Cvector = zeros(T_row,1);
    Cvector(1:22,1) = Thermal_V;
    Cvector(155:176,1) = Thermal_V;
    Cvector(23:3:154,1) = Thermal_V;
    Cvector(24:3:154,1) = Nd;
    Cvector(25:3:154,1) = Na;

    Cmatrix = spdiags(Cvector,0,T_row,T_row);
    Jaco_scaled = Jaco * Cmatrix;
    Rvector = 1./sum(abs(Jaco_scaled),2);
    Rmatrix = spdiags(Rvector,0,T_row,T_row);
    Jaco_scaled = Rmatrix* Jaco_scaled;
    res_scaled = Rmatrix *res;
    update_scaled = Jaco_scaled \ (-res_scaled);
    update_vector_DD(:,Newton) = Cmatrix* update_scaled;

    phi(1:22,1) = phi(1:22,1) + update_vector_DD(1:22,Newton);
    phi(23:66,1) = phi(23:66,1) + update_vector_DD(23:3:154,Newton);
    phi(67:88,1) = phi(67:88,1) + update_vector_DD(155:176,Newton);

    elec(23:66,1) = elec(23:66,1) + update_vector_DD(24:3:154,Newton);
    hole(23:66,1) = hole(23:66,1) + update_vector_DD(25:3:154,Newton);

    update_poisson_DD(:,Newton) = abs(update_vector_DD(23:3:154,Newton));
end
