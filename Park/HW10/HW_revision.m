clear; close all; clc;

% due to : 2022.04.07 (Thur) / HW10

%%%%% HW edit %%%%%%%

%%%% parameter %%%%

q=1.602e-19;
nint=1e+16;     %1/m^3
Nd=-1e24;   %1/m^3
k_B=1.38065e-23; % 
T=300; %K
esi = 11.7; eox = 3.9; e0=8.854e-12;  % Permittivity, F/m
coeff=(0.25*q)/e0;

Thermal_V = k_B*T/q;

%%%%%%%% Read Files %%%%%%%%%%
fileID = fopen('Vertex1.txt');

vertex = textscan(fileID,'%f%f','Delimiter','\t');

Vertex = horzcat(vertex{1},vertex{2});
fclose(fileID);
[V_row,V_col]=size(Vertex);

fileID1 = fopen('Element1.txt');
element = textscan(fileID1,'%f%f%f','Delimiter','\t');

Element = horzcat(element{1},element{2},element{3});
fclose(fileID1);
[E_row,col]=size(Element);

fileID2 = fopen('contact.txt');
contact = textscan(fileID2,'%f%f','Delimiter','\t');

Contact = horzcat(contact{1},contact{2});
fclose(fileID2);
[C_row,C_col]=size(Contact);

index = V_row;

%%%%%%%% Specify Region 1,2 %%%%%%%%%%

fileID2 = fopen('Element_region1.txt');
element_Region1 = textscan(fileID2,'%f%f%f','Delimiter','\t');

E_R1 = horzcat(element_Region1{1},element_Region1{2},element_Region1{3});
fclose(fileID2);

fileID3 = fopen('Element_region2.txt');
element_Region2 = textscan(fileID3,'%f%f%f','Delimiter','\t');

E_R2 = horzcat(element_Region2{1},element_Region2{2},element_Region2{3});
fclose(fileID3);

fileID4 = fopen('Element_region3.txt');
element_Region3 = textscan(fileID4,'%f%f%f','Delimiter','\t');

E_R3 = horzcat(element_Region3{1},element_Region3{2},element_Region3{3});
fclose(fileID4);

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
Area=zeros(E_row,1);
R=zeros(E_row,1);

for ii = 1:R_row1
    a12 = [(Vertex(E_R1(ii,2),1)-Vertex(E_R1(ii,1),1)),(Vertex(E_R1(ii,2),2)-Vertex(E_R1(ii,1),2))];
    a13 = [(Vertex(E_R1(ii,3),1)-Vertex(E_R1(ii,1),1)),(Vertex(E_R1(ii,3),2)-Vertex(E_R1(ii,1),2))];
    a23 = [(Vertex(E_R1(ii,3),1)-Vertex(E_R1(ii,2),1)),(Vertex(E_R1(ii,3),2)-Vertex(E_R1(ii,2),2))];
    a21=-a12;
    a31=-a13;
    a32=-a23;

    Area(ii,1) = (1/2)*abs(a12(1,1)*a13(1,2)-a12(1,2)*a13(1,1)); % m^2
    R(ii,1) = (L(ii,1)*L(ii,2)*L(ii,3))/(4*Area(ii,1));
end

for ii = 1:R_row2
    a12 = [(Vertex(E_R2(ii,2),1)-Vertex(E_R2(ii,1),1)),(Vertex(E_R2(ii,2),2)-Vertex(E_R2(ii,1),2))];
    a13 = [(Vertex(E_R2(ii,3),1)-Vertex(E_R2(ii,1),1)),(Vertex(E_R2(ii,3),2)-Vertex(E_R2(ii,1),2))];
    a23 = [(Vertex(E_R2(ii,3),1)-Vertex(E_R2(ii,2),1)),(Vertex(E_R2(ii,3),2)-Vertex(E_R2(ii,2),2))];
    a21=-a12;
    a31=-a13;
    a32=-a23;

    Area(ii+R_row1,1) = (1/2)*abs(a12(1,1)*a13(1,2)-a12(1,2)*a13(1,1)); % m^2
    R(ii+R_row1,1) = (L(ii+R_row1,1)*L(ii+R_row1,2)*L(ii+R_row1,3))/(4*Area(ii+R_row1,1));
end

for ii = 1:R_row3
    a12 = [(Vertex(E_R3(ii,2),1)-Vertex(E_R3(ii,1),1)),(Vertex(E_R3(ii,2),2)-Vertex(E_R3(ii,1),2))];
    a13 = [(Vertex(E_R3(ii,3),1)-Vertex(E_R3(ii,1),1)),(Vertex(E_R3(ii,3),2)-Vertex(E_R3(ii,1),2))];
    a23 = [(Vertex(E_R3(ii,3),1)-Vertex(E_R3(ii,2),1)),(Vertex(E_R3(ii,3),2)-Vertex(E_R3(ii,2),2))];
    a21=-a12;
    a31=-a13;
    a32=-a23;

    Area(ii+R_row2+R_row1,1) = (1/2)*abs(a12(1,1)*a13(1,2)-a12(1,2)*a13(1,1)); % m^2
    R(ii+R_row2+R_row1,1) = (L(ii+R_row2+R_row1,1)*L(ii+R_row2+R_row1,2)*L(ii+R_row2+R_row1,3))/(4*Area(ii+R_row2+R_row1,1));
end

%%%%%%% Area %%%%%%%%
edge=zeros(E_row,3);

for ii = 1:R_row1
    for jj = 1:3
        if jj == 1 % between 1 and 2
            edge(ii,jj) = sqrt(abs(R(ii,1)^2-(L(ii,jj)/2)^2));
        elseif jj == 2 % between 2 and 3
            edge(ii,jj) = sqrt(abs(R(ii,1)^2-(L(ii,jj)/2)^2));
        elseif jj ==3 % between 3 and 1
            edge(ii,jj) = sqrt(abs(R(ii,1)^2-(L(ii,jj)/2)^2));
        end
    end
end

for ii = 1:R_row2
    for jj = 1:3
        if jj == 1 % between 1 and 2
            edge(ii+R_row1,jj) = sqrt(abs(R(ii+R_row1,1)^2-(L(ii+R_row1,jj)/2)^2));
        elseif jj == 2 % between 2 and 3
            edge(ii+R_row1,jj) = sqrt(abs(R(ii+R_row1,1)^2-(L(ii+R_row1,jj)/2)^2));
        elseif jj ==3 % between 3 and 1
            edge(ii+R_row1,jj) = sqrt(abs(R(ii+R_row1,1)^2-(L(ii+R_row1,jj)/2)^2));
        end
    end
end

for ii = 1:R_row3
    for jj = 1:3
        if jj == 1 % between 1 and 2
            edge(ii+R_row2+R_row1,jj) = sqrt(abs(R(ii+R_row2+R_row1,1)^2-(L(ii+R_row2+R_row1,jj)/2)^2));
        elseif jj == 2 % between 2 and 3
            edge(ii+R_row2+R_row1,jj) = sqrt(abs(R(ii+R_row2+R_row1,1)^2-(L(ii+R_row2+R_row1,jj)/2)^2));
        elseif jj ==3 % between 3 and 1
            edge(ii+R_row2+R_row1,jj) = sqrt(abs(R(ii+R_row2+R_row1,1)^2-(L(ii+R_row2+R_row1,jj)/2)^2));
        end
    end
end


L = L* 1e-9; %m
edge = edge*1e-9; %m
Area=Area*1e-18; %m^2
R=R*1e-9; %m 

%%%%% Region vertex information %%%%%

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

%%%%%%%% Initial Poisson %%%%%%%%%%

Jaco = zeros(T_row,T_row);
res = zeros(T_row,1);

K = find(E_R1 == Table(1,2));

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

                    Jaco(ii, V3) = Jaco(ii, V3)+ eox*(edge(K(rr,1),jj+2)/L(K(rr,1),jj+2));
                    Jaco(ii, V2) = Jaco(ii, V2)+eox*(edge(K(rr,1),jj)/L(K(rr,1),jj));
                    Jaco(ii, V1) = Jaco(ii, V1)-eox*(edge(K(rr,1),jj+2)/L(K(rr,1),jj+2))-eox*(edge(K(rr,1),jj)/L(K(rr,1),jj));

                elseif K(rr,1) > R_row1 && K(rr,1) <= 2*R_row1
                    jj = 2;
                    K(rr,1) = K(rr,1)-R_row1;

                    V3 = find(Table(1:VR1,2)==E_R1(K(rr,1),jj+1)) ;
                    V2 = find(Table(1:VR1,2)==E_R1(K(rr,1),jj));
                    V1 = find(Table(1:VR1,2)==E_R1(K(rr,1),jj-1));

                    Jaco(ii, V3) = Jaco(ii, V3)+eox*(edge(K(rr,1),jj)/L(K(rr,1),jj));
                    Jaco(ii, V2) =  Jaco(ii, V2)-eox*(edge(K(rr,1),jj)/L(K(rr,1),jj))-eox*(edge(K(rr,1),jj-1)/L(K(rr,1),jj-1));
                    Jaco(ii, V1) = Jaco(ii, V1)+eox*(edge(K(rr,1),jj-1)/L(K(rr,1),jj-1));



                elseif K(rr,1) > 2*R_row1 && K(rr,1) <= 3*R_row1
                    jj = 3;
                    K(rr,1) = K(rr,1)-2*R_row1;

                    V3 = find(Table(1:VR1,2)==E_R1(K(rr,1),jj-2)) ;
                    V2 = find(Table(1:VR1,2)==E_R1(K(rr,1),jj-1));
                    V1 = find(Table(1:VR1,2)==E_R1(K(rr,1),jj));

                    Jaco(ii, V3) = Jaco(ii, V3)+eox*(edge(K(rr,1),jj)/L(K(rr,1),jj));
                    Jaco(ii, V2) = Jaco(ii, V2)+eox*(edge(K(rr,1),jj-1)/L(K(rr,1),jj-1));
                    Jaco(ii, V1) = Jaco(ii, V1)-eox*(edge(K(rr,1),jj)/L(K(rr,1),jj))-eox*(edge(K(rr,1),jj-1)/L(K(rr,1),jj-1));

                end
            end
        end

    elseif Table(ii,1) == 2     %region2

        K = find(E_R2 == Table(ii,2));
        K_row = size(K,1);

        for rr = 1:K_row
            if Table(ii,3) == Potential

                if K(rr,1) <= R_row2
                    jj = 1;

                    V3 = find(Table(1+VR1:VR1+VR2,2)==E_R2(K(rr,1),jj+2)) ;
                    V2 = find(Table(1+VR1:VR1+VR2,2)==E_R2(K(rr,1),jj+1));
                    V1 = find(Table(1+VR1:VR1+VR2,2)==E_R2(K(rr,1),jj));

                    res(ii,1)=res(ii,1)+coeff*(edge(K(rr,1)+R_row1,jj)*L(K(rr,1)+R_row1,jj)+edge(K(rr,1)+R_row1,jj+2)*L(K(rr,1)+R_row1,jj+2))*Nd;
                    Jaco(ii, VR1+V3) = Jaco(ii, VR1+V3)+ esi*(edge(K(rr,1)+R_row1,jj+2)/L(K(rr,1)+R_row1,jj+2));
                    Jaco(ii, VR1+V2) = Jaco(ii, VR1+V2)+esi*(edge(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj));
                    Jaco(ii, VR1+V1) = Jaco(ii, VR1+V1)-esi*(edge(K(rr,1)+R_row1,jj+2)/L(K(rr,1)+R_row1,jj+2))-esi*(edge(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj));

                elseif K(rr,1) > R_row2 && K(rr,1) <= 2*R_row2
                    jj = 2;
                    K(rr,1) = K(rr,1)-R_row2;

                    V3 = find(Table(1+VR1:VR1+VR2,2)==E_R2(K(rr,1),jj+1)) ;
                    V2 = find(Table(1+VR1:VR1+VR2,2)==E_R2(K(rr,1),jj));
                    V1 = find(Table(1+VR1:VR1+VR2,2)==E_R2(K(rr,1),jj-1));

                    res(ii,1)=res(ii,1)+coeff*(edge(K(rr,1)+R_row1,jj)*L(K(rr,1)+R_row1,jj)+edge(K(rr,1)+R_row1,jj-1)*L(K(rr,1)+R_row1,jj-1))*Nd;
                    Jaco(ii, VR1+V3) = Jaco(ii,VR1+V3)+esi*(edge(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj));
                    Jaco(ii, VR1+V2) =  Jaco(ii, VR1+V2)-esi*(edge(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj))-esi*(edge(K(rr,1)+R_row1,jj-1)/L(K(rr,1)+R_row1,jj-1));
                    Jaco(ii, VR1+V1) = Jaco(ii, VR1+V1)+esi*(edge(K(rr,1)+R_row1,jj-1)/L(K(rr,1)+R_row1,jj-1));

                elseif K(rr,1) > 2*R_row2 && K(rr,1) <= 3*R_row2
                    jj = 3;
                    K(rr,1) = K(rr,1)-2*R_row2;

                    V3 = find(Table(1+VR1:VR1+VR2,2)==E_R2(K(rr,1),jj-2)) ;
                    V2 = find(Table(1+VR1:VR1+VR2,2)==E_R2(K(rr,1),jj-1));
                    V1 = find(Table(1+VR1:VR1+VR2,2)==E_R2(K(rr,1),jj));

                    res(ii,1)=res(ii,1)+coeff*(edge(K(rr,1)+R_row1,jj)*L(K(rr,1)+R_row1,jj)+edge(K(rr,1)+R_row1,jj-1)*L(K(rr,1)+R_row1,jj-1))*Nd;
                    Jaco(ii, VR1+V3) = Jaco(ii, VR1+V3)+esi*(edge(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj));
                    Jaco(ii, VR1+V2) = Jaco(ii, VR1+V2)+esi*(edge(K(rr,1)+R_row1,jj-1)/L(K(rr,1)+R_row1,jj-1));
                    Jaco(ii, VR1+V1) = Jaco(ii, VR1+V1)-esi*(edge(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj))-esi*(edge(K(rr,1)+R_row1,jj-1)/L(K(rr,1)+R_row1,jj-1));
                end
            end
        end


    elseif Table(ii,1) == 3     % Region3

        K = find(E_R3 == Table(ii,2));
        K_row = size(K,1);

        for rr = 1:K_row
            if Table(ii,3) == Potential

                if K(rr,1) <= R_row3
                    jj = 1;

                    V3 = find(Table(1+VR1+VR2:T,2)==E_R3(K(rr,1),jj+2)) ;
                    V2 = find(Table(1+VR1+VR2:T,2)==E_R3(K(rr,1),jj+1));
                    V1 = find(Table(1+VR1+VR2:T,2)==E_R3(K(rr,1),jj));

                    Jaco(ii, VR1+VR2+V3) = Jaco(ii, VR1+VR2+V3)+ eox*(edge(K(rr,1)+R_row2+R_row1,jj+2)/L(K(rr,1)+R_row2+R_row1,jj+2));
                    Jaco(ii, VR1+VR2+V2) = Jaco(ii, VR1+VR2+V2)+eox*(edge(K(rr,1)+R_row2+R_row1,jj)/L(K(rr,1)+R_row2+R_row1,jj));
                    Jaco(ii, VR1+VR2+V1) = Jaco(ii, VR1+VR2+V1)-eox*(edge(K(rr,1)+R_row2+R_row1,jj+2)/L(K(rr,1)+R_row2+R_row1,jj+2))-eox*(edge(K(rr,1)+R_row2+R_row1,jj)/L(K(rr,1)+R_row2+R_row1,jj));

                elseif K(rr,1) > R_row3 && K(rr,1) <= 2*R_row3
                    jj = 2;
                    K(rr,1) = K(rr,1)-R_row3;

                    V3 = find(Table(1+VR1+VR2:T,2)==E_R3(K(rr,1),jj+1)) ;
                    V2 = find(Table(1+VR1+VR2:T,2)==E_R3(K(rr,1),jj));
                    V1 = find(Table(1+VR1+VR2:T,2)==E_R3(K(rr,1),jj-1));

                    Jaco(ii, VR1+VR2+V3) = Jaco(ii, VR1+VR2+V3)+eox*(edge(K(rr,1)+R_row2+R_row1,jj)/L(K(rr,1)+R_row2+R_row1,jj));
                    Jaco(ii, VR1+VR2+V2) =  Jaco(ii, VR1+VR2+V2)-eox*(edge(K(rr,1)+R_row2+R_row1,jj)/L(K(rr,1)+R_row2+R_row1,jj))-eox*(edge(K(rr,1)+R_row2+R_row1,jj-1)/L(K(rr,1)+R_row2+R_row1,jj-1));
                    Jaco(ii, VR1+VR2+V1) = Jaco(ii, VR1+VR2+V1)+eox*(edge(K(rr,1)+R_row2+R_row1,jj-1)/L(K(rr,1)+R_row2+R_row1,jj-1));

                elseif K(rr,1) > 2*R_row3 && K(rr,1) <= 3*R_row3
                    jj = 3;
                    K(rr,1) = K(rr,1)-2*R_row3;

                    V3 = find(Table(1+VR1+VR2:T,2)==E_R3(K(rr,1),jj-2)) ;
                    V2 = find(Table(1+VR1+VR2:T,2)==E_R3(K(rr,1),jj-1));
                    V1 = find(Table(1+VR1+VR2:T,2)==E_R3(K(rr,1),jj));

                    Jaco(ii, VR1+VR2+V3) = Jaco(ii, VR1+VR2+V3)+eox*(edge(K(rr,1)+R_row2+R_row1,jj)/L(K(rr,1)+R_row2+R_row1,jj));
                    Jaco(ii, VR1+VR2+V2) = Jaco(ii, VR1+VR2+V2)+eox*(edge(K(rr,1)+R_row2+R_row1,jj-1)/L(K(rr,1)+R_row2+R_row1,jj-1));
                    Jaco(ii, VR1+VR2+V1) = Jaco(ii, VR1+VR2+V1)-eox*(edge(K(rr,1)+R_row2+R_row1,jj)/L(K(rr,1)+R_row2+R_row1,jj))-eox*(edge(K(rr,1)+R_row2+R_row1,jj-1)/L(K(rr,1)+R_row2+R_row1,jj-1));
                end
            end
        end
    end

end

%Boundary Condition

    Jaco(11,:)=Jaco(11,:)+Jaco(7,:);
    Jaco(7,:)=0; Jaco(7,7)=1; Jaco(7,11)=-1;
    Jaco(12,:)=Jaco(12,:)+Jaco(8,:);
    Jaco(8,:)=0; Jaco(8,8)=1; Jaco(8,12)=-1;
    Jaco(13,:)=Jaco(13,:)+Jaco(9,:);
    Jaco(9,:)=0; Jaco(9,9)=1; Jaco(9,13)=-1;

    Jaco(21,:)=Jaco(21,:)+Jaco(25,:);
    Jaco(25,:)=0; Jaco(25,25)=1; Jaco(25,21)=-1;
    Jaco(22,:)=Jaco(22,:)+Jaco(26,:);
    Jaco(26,:)=0; Jaco(26,26)=1; Jaco(26,22)=-1;
    Jaco(23,:)=Jaco(23,:)+Jaco(27,:);
    Jaco(27,:)=0; Jaco(27,27)=1; Jaco(27,23)=-1;

%%%%% Dirichlet Boundary Condition %%%%%

for n=1:C_col
    for i =1: C_row
        if i <= 2
            Jaco(Contact(i,n),:) =0;
            Jaco(Contact(i,n),Contact(i,n)) = 1;
            res(Contact(i,n),1)=0.33374;

        elseif i > 2 && i <= 4
            Jaco(6+Contact(i,n),:) =0;
            Jaco(6+Contact(i,n),6+Contact(i,n)) = 1;
            res(6+Contact(i,n),1)= 0.33374;

        end
    end
end

phi = Jaco \ (res);

elec = zeros(33,1);
hole = zeros(33,1);

elec(7:1:27,1) = nint*exp(phi(7:1:27,1)/Thermal_V);
hole(7:1:27,1) = nint*exp(-phi(7:1:27,1)/Thermal_V);

%%%%%%%% Fully-Coupled %%%%%%%%

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

            Table(3*iVertex+7,2) = V_R2(iVertex,1);
            Table(3*iVertex+8,2) = V_R2(iVertex,1);
            Table(3*iVertex+9,2) = V_R2(iVertex,1);

            Table(3*iVertex+7,3) = Potential;
            Table(3*iVertex+8,3) = eDensity;
            Table(3*iVertex+9,3) = hDensity;
        end

    elseif iRegion ==3  % Oxide , only potential

        Table(1+VR1+3*VR2:T1,1) = 3;
        
        for iVertex =  1:VR3
            Table(VR1+3*VR2+iVertex,2) = V_R3(iVertex,1);
            Table(VR1+3*VR2+iVertex,3) = Potential;
        end
         
    end
end

%%%%%%%%  Poisson, eDensity, hDensity %%%%%%%%%%

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
    
                        res(ii,1) = res(ii,1) + eox*(phi(V2,1)-phi(V1,1))*(edge(K(rr,1),jj)/L(K(rr,1),jj)) + eox*(phi(V3,1)-phi(V1,1))*(edge(K(rr,1),jj+2)/L(K(rr,1),jj+2));
                        Jaco(ii, V3) = Jaco(ii, V3)+ eox*(edge(K(rr,1),jj+2)/L(K(rr,1),jj+2));
                        Jaco(ii, V2) = Jaco(ii, V2)+eox*(edge(K(rr,1),jj)/L(K(rr,1),jj));
                        Jaco(ii, V1) = Jaco(ii, V1)-eox*(edge(K(rr,1),jj+2)/L(K(rr,1),jj+2))-eox*(edge(K(rr,1),jj)/L(K(rr,1),jj));
    
                    elseif K(rr,1) > R_row1 && K(rr,1) <= 2*R_row1
                        jj = 2;
                        K(rr,1) = K(rr,1)-R_row1;
    
                        V3 = find(Table(1:VR1,2)==E_R1(K(rr,1),jj+1)) ;
                        V2 = find(Table(1:VR1,2)==E_R1(K(rr,1),jj));
                        V1 = find(Table(1:VR1,2)==E_R1(K(rr,1),jj-1));
    
                        res(ii,1) = res(ii,1) + eox*(phi(V1,1)-phi(V2,1))*(edge(K(rr,1),jj-1)/L(K(rr,1),jj-1)) + eox*(phi(V3,1)-phi(V2,1))*(edge(K(rr,1),jj)/L(K(rr,1),jj));
                        Jaco(ii, V3) = Jaco(ii, V3)+eox*(edge(K(rr,1),jj)/L(K(rr,1),jj));
                        Jaco(ii, V2) =  Jaco(ii, V2)-eox*(edge(K(rr,1),jj)/L(K(rr,1),jj))-eox*(edge(K(rr,1),jj-1)/L(K(rr,1),jj-1));
                        Jaco(ii, V1) = Jaco(ii, V1)+eox*(edge(K(rr,1),jj-1)/L(K(rr,1),jj-1));
    
                    elseif K(rr,1) > 2*R_row1 && K(rr,1) <= 3*R_row1
                        jj = 3;
                        K(rr,1) = K(rr,1)-2*R_row1;
    
                        V3 = find(Table(1:VR1,2)==E_R1(K(rr,1),jj-2)) ;
                        V2 = find(Table(1:VR1,2)==E_R1(K(rr,1),jj-1));
                        V1 = find(Table(1:VR1,2)==E_R1(K(rr,1),jj));
    
                        res(ii,1) = res(ii,1) + eox*(phi(V3,1)-phi(V1,1))*(edge(K(rr,1),jj)/L(K(rr,1),jj)) + eox*(phi(V2,1)-phi(V1,1))*(edge(K(rr,1),jj-1)/L(K(rr,1),jj-1));
                        Jaco(ii, V3) = Jaco(ii, V3)+eox*(edge(K(rr,1),jj)/L(K(rr,1),jj));
                        Jaco(ii, V2) = Jaco(ii, V2)+eox*(edge(K(rr,1),jj-1)/L(K(rr,1),jj-1));
                        Jaco(ii, V1) = Jaco(ii, V1)-eox*(edge(K(rr,1),jj)/L(K(rr,1),jj))-eox*(edge(K(rr,1),jj-1)/L(K(rr,1),jj-1));
    
                    end
                end
            end
    
        elseif Table(ii,1) == 2     %region2
    
            K = find(E_R2 == Table(ii,2));
            K_row = size(K,1);
            i2 = rem(ii,3);
    
            if Table(ii,3) == Potential && i2 == 1
    
                for rr = 1:K_row
    
                    if K(rr,1) <= R_row2
                        jj = 1;
    
                        V3 = find(Table(1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj+2)) ;
                        V2 = find(Table(1+1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj+1));
                        V1 = find(Table(1+1+1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj));
    
                        res(ii,1) = res(ii,1) + esi*(phi(VR1+V2,1)-phi(VR1+V1,1))*(edge(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj)) + esi*(phi(VR1+V3,1)-phi(VR1+V1,1))*(edge(K(rr,1)+R_row1,jj+2)/L(K(rr,1)+R_row1,jj+2));
                        res(ii,1) = res(ii,1) + coeff*(edge(K(rr,1)+R_row1,jj)*L(K(rr,1)+R_row1,jj)+edge(K(rr,1)+R_row1,jj+2)*L(K(rr,1)+R_row1,jj+2))*(Nd-elec(VR1+V1,1)+hole(VR1+V1,1));
    
                        Jaco(ii, 3*(V3)+7) = Jaco(ii, 3*(V3)+7)+ esi*(edge(K(rr,1)+R_row1,jj+2)/L(K(rr,1)+R_row1,jj+2));
                        Jaco(ii, 3*(V2)+7) = Jaco(ii, 3*(V2)+7)+esi*(edge(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj));
                        Jaco(ii, 3*(V1)+7) = Jaco(ii, 3*(V1)+7)-esi*(edge(K(rr,1)+R_row1,jj+2)/L(K(rr,1)+R_row1,jj+2))-esi*(edge(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj));
    
                        Jaco(ii, 3*(V1)+8) = Jaco(ii, 3*(V1)+8) - coeff*(edge(K(rr,1)+R_row1,jj)*L(K(rr,1)+R_row1,jj)+edge(K(rr,1)+R_row1,jj+2)*L(K(rr,1)+R_row1,jj+2));
                        Jaco(ii, 3*(V1)+9) = Jaco(ii, 3*(V1)+9) + coeff*(edge(K(rr,1)+R_row1,jj)*L(K(rr,1)+R_row1,jj)+edge(K(rr,1)+R_row1,jj+2)*L(K(rr,1)+R_row1,jj+2));
    
                    elseif K(rr,1) > R_row2 && K(rr,1) <= 2*R_row2
                        jj = 2;
                        K(rr,1) = K(rr,1)-R_row2;
    
                        V3 = find(Table(1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj+1)) ;
                        V2 = find(Table(1+1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj));
                        V1 = find(Table(1+1+1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj-1));
    
                        res(ii,1) = res(ii,1) + esi*(phi(VR1+V1,1)-phi(VR1+V2,1))*(edge(K(rr,1)+R_row1,jj-1)/L(K(rr,1)+R_row1,jj-1)) + esi*(phi(VR1+V3,1)-phi(VR1+V2,1))*(edge(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj));
                        res(ii,1) = res(ii,1) + coeff*(edge(K(rr,1)+R_row1,jj)*L(K(rr,1)+R_row1,jj)+edge(K(rr,1)+R_row1,jj-1)*L(K(rr,1)+R_row1,jj-1))*(Nd-elec(VR1+V2,1)+hole(VR1+V2,1));
    
                        Jaco(ii, 3*(V3)+7) = Jaco(ii, 3*(V3)+7)+esi*(edge(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj));
                        Jaco(ii, 3*(V2)+7) =  Jaco(ii, 3*(V2)+7)-esi*(edge(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj))-esi*(edge(K(rr,1)+R_row1,jj-1)/L(K(rr,1)+R_row1,jj-1));
                        Jaco(ii, 3*(V1)+7) = Jaco(ii, 3*(V1)+7)+esi*(edge(K(rr,1)+R_row1,jj-1)/L(K(rr,1)+R_row1,jj-1));
    
                        Jaco(ii, 3*(V2)+8) = Jaco(ii, 3*(V2)+8) - coeff*(edge(K(rr,1)+R_row1,jj)*L(K(rr,1)+R_row1,jj)+edge(K(rr,1)+R_row1,jj-1)*L(K(rr,1)+R_row1,jj-1));
                        Jaco(ii, 3*(V2)+9) = Jaco(ii, 3*(V2)+9) + coeff*(edge(K(rr,1)+R_row1,jj)*L(K(rr,1)+R_row1,jj)+edge(K(rr,1)+R_row1,jj-1)*L(K(rr,1)+R_row1,jj-1));
    
                    elseif K(rr,1) > 2*R_row2 && K(rr,1) <= 3*R_row2
                        jj = 3;
                        K(rr,1) = K(rr,1)-2*R_row2;
    
                        V3 = find(Table(1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj-2)) ;
                        V2 = find(Table(1+1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj-1));
                        V1 = find(Table(1+1+1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj));
    
                        res(ii,1) = res(ii,1) + esi*(phi(VR1+V3,1)-phi(VR1+V1,1))*(edge(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj)) + esi*(phi(VR1+V2,1)-phi(VR1+V1,1))*(edge(K(rr,1)+R_row1,jj-1)/L(K(rr,1)+R_row1,jj-1));
                        res(ii,1) = res(ii,1) + coeff*(edge(K(rr,1)+R_row1,jj)*L(K(rr,1)+R_row1,jj)+edge(K(rr,1)+R_row1,jj-1)*L(K(rr,1)+R_row1,jj-1))*(Nd-elec(VR1+V1,1)+hole(VR1+V1,1));
    
                        Jaco(ii, 3*(V3)+7) = Jaco(ii, 3*(V3)+7)+esi*(edge(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj));
                        Jaco(ii, 3*(V2)+7) = Jaco(ii, 3*(V2)+7)+esi*(edge(K(rr,1)+R_row1,jj-1)/L(K(rr,1)+R_row1,jj-1));
                        Jaco(ii, 3*(V1)+7) = Jaco(ii, 3*(V1)+7)-esi*(edge(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj))-esi*(edge(K(rr,1)+R_row1,jj-1)/L(K(rr,1)+R_row1,jj-1));
    
                        Jaco(ii, 3*(V1)+8) = Jaco(ii, 3*(V1)+8) - coeff*(edge(K(rr,1)+R_row1,jj)*L(K(rr,1)+R_row1,jj)+edge(K(rr,1)+R_row1,jj-1)*L(K(rr,1)+R_row1,jj-1));
                        Jaco(ii, 3*(V1)+9) = Jaco(ii, 3*(V1)+9) + coeff*(edge(K(rr,1)+R_row1,jj)*L(K(rr,1)+R_row1,jj)+edge(K(rr,1)+R_row1,jj-1)*L(K(rr,1)+R_row1,jj-1));
                    end
                end
    
            elseif Table(ii,3) == eDensity && i2 == 2
    
                for rr = 1:K_row
                    if K(rr,1) <= R_row2
                        jj = 1;
    
                        V3 = find(Table(1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj+2)) ;
                        V2 = find(Table(1+1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj+1));
                        V1 = find(Table(2+1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj));
    
                        res(ii,1) =  elec(VR1+V1,1) - nint*exp(phi(VR1+V1,1)/Thermal_V);
    
                        Jaco(ii, 3*(V1)+8) = 1 ;
                        Jaco(ii, 3*(V1)+7) = -(nint/Thermal_V)*exp(phi(VR1+V1,1)/Thermal_V);
    
    
                    elseif K(rr,1) > R_row2 && K(rr,1) <= 2*R_row2
                        jj = 2;
                        K(rr,1) = K(rr,1)-R_row2;
    
                        V3 = find(Table(1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj+1)) ;
                        V2 = find(Table(1+1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj));
                        V1 = find(Table(1+1+1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj-1));
    
                        res(ii,1) =  elec(VR1+V2,1) - nint*exp(phi(VR1+V2,1)/Thermal_V);
    
                        Jaco(ii, 3*(V2)+8) = 1;
                        Jaco(ii, 3*(V2)+7) = -(nint/Thermal_V)*exp(phi(VR1+V2,1)/Thermal_V);
    
                    elseif K(rr,1) > 2*R_row2 && K(rr,1) <= 3*R_row2
                        jj = 3;
                        K(rr,1) = K(rr,1)-2*R_row2;
    
                        V3 = find(Table(1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj-2)) ;
                        V2 = find(Table(1+1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj-1));
                        V1 = find(Table(1+1+1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj));
    
                        res(ii,1) =  elec(VR1+V1,1) - nint*exp(phi(VR1+V1,1)/Thermal_V);
    
                        Jaco(ii, 3*(V1)+8) = 1;
                        Jaco(ii, 3*(V1)+7) = -(nint/Thermal_V)*exp(phi(VR1+V1,1)/Thermal_V);
                    end
                end
    
    
            elseif Table(ii,3) == hDensity && i2 == 0
    
                for rr = 1:K_row
    
                    if K(rr,1) <= R_row2
                        jj = 1;
    
                        V3 = find(Table(1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj+2)) ;
                        V2 = find(Table(1+1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj+1));
                        V1 = find(Table(2+1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj));
    
                        res(ii,1) = hole(VR1+V1,1)-nint*exp(-phi(VR1+V1,1)/Thermal_V);
    
                        Jaco(ii, 3*(V1)+9) = 1; 
                        Jaco(ii, 3*(V1)+7) = (nint/Thermal_V)*exp(-phi(VR1+V1,1)/Thermal_V);
    
    
                    elseif K(rr,1) > R_row2 && K(rr,1) <= 2*R_row2
                        jj = 2;
                        K(rr,1) = K(rr,1)-R_row2;
    
                        V3 = find(Table(1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj+1)) ;
                        V2 = find(Table(1+1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj));
                        V1 = find(Table(1+1+1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj-1));
    
                        res(ii,1) =  hole(VR1+V2,1)-nint*exp(-phi(VR1+V2,1)/Thermal_V);
    
                        Jaco(ii, 3*(V2)+9) =  1; 
                        Jaco(ii, 3*(V2)+7) = (nint/Thermal_V)*exp(-phi(VR1+V2,1)/Thermal_V);
    
    
                    elseif K(rr,1) > 2*R_row2 && K(rr,1) <= 3*R_row2
                        jj = 3;
                        K(rr,1) = K(rr,1)-2*R_row2;
    
                        V3 = find(Table(1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj-2)) ;
                        V2 = find(Table(1+1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj-1));
                        V1 = find(Table(1+1+1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj));
    
                        res(ii,1) =  hole(VR1+V1,1)-nint*exp(-phi(VR1+V1,1)/Thermal_V);
    
                        Jaco(ii, 3*(V1)+9) =  1; 
                        Jaco(ii, 3*(V1)+7) =  (nint/Thermal_V)*exp(-phi(VR1+V1,1)/Thermal_V);
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
    
                        res(ii,1) = res(ii,1) + eox*(phi(VR1+VR2+V2,1)-phi(VR1+VR2+V1,1))*(edge(K(rr,1)+R_row2+R_row1,jj)/L(K(rr,1)+R_row2+R_row1,jj)) + eox*(phi(VR1+VR2+V3,1)-phi(VR1+VR2+V1,1))*(edge(K(rr,1)+R_row2+R_row1,jj+2)/L(K(rr,1)+R_row2+R_row1,jj+2));
    
                        Jaco(ii, VR1+3*VR2+V3) = Jaco(ii, VR1+3*VR2+V3)+ eox*(edge(K(rr,1)+R_row2+R_row1,jj+2)/L(K(rr,1)+R_row2+R_row1,jj+2));
                        Jaco(ii, VR1+3*VR2+V2) = Jaco(ii, VR1+3*VR2+V2)+eox*(edge(K(rr,1)+R_row2+R_row1,jj)/L(K(rr,1)+R_row2+R_row1,jj));
                        Jaco(ii, VR1+3*VR2+V1) = Jaco(ii, VR1+3*VR2+V1)-eox*(edge(K(rr,1)+R_row2+R_row1,jj+2)/L(K(rr,1)+R_row2+R_row1,jj+2))-eox*(edge(K(rr,1)+R_row2+R_row1,jj)/L(K(rr,1)+R_row2+R_row1,jj));
    
                    elseif K(rr,1) > R_row3 && K(rr,1) <= 2*R_row3
                        
                        jj = 2;
                        K(rr,1) = K(rr,1)-R_row3;
    
                        V3 = find(Table(1+VR1+3*VR2:T1,2)==E_R3(K(rr,1),jj+1)) ;
                        V2 = find(Table(1+VR1+3*VR2:T1,2)==E_R3(K(rr,1),jj));
                        V1 = find(Table(1+VR1+3*VR2:T1,2)==E_R3(K(rr,1),jj-1));
    
                        res(ii,1) = res(ii,1) + eox*(phi(VR1+VR2+V1,1)-phi(VR1+VR2+V2,1))*(edge(K(rr,1)+R_row2+R_row1,jj-1)/L(K(rr,1)+R_row2+R_row1,jj-1)) + eox*(phi(VR1+VR2+V3,1)-phi(VR1+VR2+V2,1))*(edge(K(rr,1)+R_row2+R_row1,jj)/L(K(rr,1)+R_row2+R_row1,jj));
    
                        Jaco(ii, VR1+3*VR2+V3) = Jaco(ii, VR1+3*VR2+V3)+eox*(edge(K(rr,1)+R_row2+R_row1,jj)/L(K(rr,1)+R_row2+R_row1,jj));
                        Jaco(ii, VR1+3*VR2+V2) =  Jaco(ii, VR1+3*VR2+V2)-eox*(edge(K(rr,1)+R_row2+R_row1,jj)/L(K(rr,1)+R_row2+R_row1,jj))-eox*(edge(K(rr,1)+R_row2+R_row1,jj-1)/L(K(rr,1)+R_row2+R_row1,jj-1));
                        Jaco(ii, VR1+3*VR2+V1) = Jaco(ii, VR1+3*VR2+V1)+eox*(edge(K(rr,1)+R_row2+R_row1,jj-1)/L(K(rr,1)+R_row2+R_row1,jj-1));
    
    
                    elseif K(rr,1) > 2*R_row3 && K(rr,1) <= 3*R_row3
                        
                        jj = 3;
                        K(rr,1) = K(rr,1)-2*R_row3;
    
                        V3 = find(Table(1+VR1+3*VR2:T1,2)==E_R3(K(rr,1),jj-2)) ;
                        V2 = find(Table(1+VR1+3*VR2:T1,2)==E_R3(K(rr,1),jj-1));
                        V1 = find(Table(1+VR1+3*VR2:T1,2)==E_R3(K(rr,1),jj));

                        res(ii,1) = res(ii,1) + eox*(phi(VR1+VR2+V3,1)-phi(VR1+VR2+V1,1))*(edge(K(rr,1)+R_row2+R_row1,jj)/L(K(rr,1)+R_row2+R_row1,jj)) + eox*(phi(VR1+VR2+V2,1)-phi(VR1+VR2+V1,1))*(edge(K(rr,1)+R_row2+R_row1,jj-1)/L(K(rr,1)+R_row2+R_row1,jj-1));

                        Jaco(ii, VR1+3*VR2+V3) = Jaco(ii, VR1+3*VR2+V3) + eox*(edge(K(rr,1)+R_row2+R_row1,jj)/L(K(rr,1)+R_row2+R_row1,jj));
                        Jaco(ii, VR1+3*VR2+V2) = Jaco(ii, VR1+3*VR2+V2) + eox*(edge(K(rr,1)+R_row2+R_row1,jj-1)/L(K(rr,1)+R_row2+R_row1,jj-1));
                        Jaco(ii, VR1+3*VR2+V1) = Jaco(ii, VR1+3*VR2+V1) - eox*(edge(K(rr,1)+R_row2+R_row1,jj)/L(K(rr,1)+R_row2+R_row1,jj))-eox*(edge(K(rr,1)+R_row2+R_row1,jj-1)/L(K(rr,1)+R_row2+R_row1,jj-1));

                    end
                end
            end
        end
    end


    %Boundary Condition

    for rr = 1:i_row1

        Vi11 = find(Table(1:VR1,2)==interface1(rr,1));
        Vi21 = find(Table(1+VR1:3:VR1+3*VR2,2)==interface1(rr,1));

        Jaco(3*Vi21+7,:)=Jaco(3*Vi21+7,:)+Jaco(Vi11,:);
        Jaco(Vi11,:)=0; Jaco(Vi11,Vi11)=1; Jaco(Vi11,3*Vi21+7)=-1;

        res(3*Vi21+7,1)=res(3*Vi21+7,1)+res(Vi11,1);
        res(Vi11,1)=0;
    end

    for rr = 1:i_row2

        Vi12 = find(Table(1+VR1+3*VR2:T1,2)==interface2(rr,1));
        Vi22 = find(Table(1+VR1:3:VR1+3*VR2,2)==interface2(rr,1));

        Jaco(3*Vi22+7,:)=Jaco(3*Vi22+7,:)+Jaco(VR1+3*VR2+Vi12,:);
        Jaco(VR1+3*VR2+Vi12,:)=0; Jaco(VR1+3*VR2+Vi12,VR1+3*VR2+Vi12)=1; Jaco(VR1+3*VR2+Vi12,3*Vi22+7)=-1;

        res(3*Vi22+7,1)=res(3*Vi22+7,1)+res(VR1+3*VR2+Vi12,1);
        res(VR1+3*VR2+Vi12,1)=0;

    end

    %%%%% Dirichlet Boundary Condition %%%%%

    for n=1:C_col
        for i =1: C_row
            if i <= 2
                Jaco(Contact(i,n),:) =0;
                Jaco(Contact(i,n),Contact(i,n)) = 1;
                res(Contact(i,n),1)=phi(Contact(i,n),1)-0.33374;

            elseif i > 2 && i <= 4
                Jaco(36+Contact(i,n),:) =0;
                Jaco(36+Contact(i,n),36+Contact(i,n)) = 1;
                res(36+Contact(i,n),1)= phi(6+Contact(i,n),1)-0.33374;

            end
        end
    end

        %scaling
        Cvector = zeros(T_row,1);
        Cvector(1:9,1) = Thermal_V;
        Cvector(55:63,1) = Thermal_V;
        Cvector(10:3:54,1) = Thermal_V;
        Cvector(11:3:54,1) = Nd;
        Cvector(12:3:54,1) = Nd;
    
        Cmatrix = spdiags(Cvector,0,T_row,T_row);
        Jaco_scaled = Jaco * Cmatrix;
        Rvector = 1./sum(abs(Jaco_scaled),2);
        Rmatrix = spdiags(Rvector,0,T_row,T_row);
        Jaco_scaled = Rmatrix* Jaco_scaled;
        res_scaled = Rmatrix *res;
        update_scaled = Jaco_scaled \ (-res_scaled);
        update_vector(:,Newton) = Cmatrix* update_scaled;
    
        phi(1:9,1) = phi(1:9,1) + update_vector(1:9,Newton);
        phi(10:24,1) = phi(10:24,1) + update_vector(10:3:54,Newton);
        phi(25:33,1) = phi(25:33,1) + update_vector(55:63,Newton);
    
        elec(10:24,1) = elec(10:24,1) + update_vector(11:3:54,Newton);
        hole(10:24,1) = hole(10:24,1) + update_vector(12:3:54,Newton);

        update_poisson(:,Newton) = abs(update_vector(10:3:54,Newton));
end