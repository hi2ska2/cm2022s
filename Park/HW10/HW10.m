clear; close all; clc;

% due to : 2022.04.07 (Thur) / HW10

%%%%% HW10, repeat 2.11.1 %%%%%%%

%%%% parameter %%%%

q=1.602e-19;
nint=1e+16;     %1/m^3
Na=1e24;   %1/m^3
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

%%%% Potential %%%%
Jaco_i = zeros(index,index);
res_i = zeros(index,1);

for ii=1:R_row1
    for jj=1:3

        if jj==1
            Jaco_i(E_R1(ii,jj), E_R1(ii,jj+2)) = Jaco_i(E_R1(ii,jj), E_R1(ii,jj+2))+ eox*(edge(ii,jj+2)/L(ii,jj+2));
            Jaco_i(E_R1(ii,jj), E_R1(ii,jj+1)) = Jaco_i(E_R1(ii,jj), E_R1(ii,jj+1))+eox*(edge(ii,jj)/L(ii,jj));
            Jaco_i(E_R1(ii,jj), E_R1(ii,jj)) = Jaco_i(E_R1(ii,jj), E_R1(ii,jj))-eox*(edge(ii,jj+2)/L(ii,jj+2))-eox*(edge(ii,jj)/L(ii,jj));

        elseif jj ==2
            Jaco_i(E_R1(ii,jj), E_R1(ii,jj+1)) = Jaco_i(E_R1(ii,jj), E_R1(ii,jj+1))+eox*(edge(ii,jj)/L(ii,jj));
            Jaco_i(E_R1(ii,jj), E_R1(ii,jj)) =  Jaco_i(E_R1(ii,jj), E_R1(ii,jj))-eox*(edge(ii,jj)/L(ii,jj))-eox*(edge(ii,jj-1)/L(ii,jj-1));
            Jaco_i(E_R1(ii,jj), E_R1(ii,jj-1)) = Jaco_i(E_R1(ii,jj), E_R1(ii,jj-1))+eox*(edge(ii,jj-1)/L(ii,jj-1));

        elseif jj == 3
            Jaco_i(E_R1(ii,jj), E_R1(ii,jj-2)) = Jaco_i(E_R1(ii,jj), E_R1(ii,jj-2))+eox*(edge(ii,jj)/L(ii,jj));
            Jaco_i(E_R1(ii,jj), E_R1(ii,jj-1)) = Jaco_i(E_R1(ii,jj), E_R1(ii,jj-1))+eox*(edge(ii,jj-1)/L(ii,jj-1));
            Jaco_i(E_R1(ii,jj), E_R1(ii,jj)) = Jaco_i(E_R1(ii,jj), E_R1(ii,jj))-eox*(edge(ii,jj)/L(ii,jj))-eox*(edge(ii,jj-1)/L(ii,jj-1));
        end
    end
end

for ii=1:R_row2
    for jj=1:3

        if jj==1
            res_i(E_R2(ii,jj),1)=coeff*(edge(ii+R_row1,jj)*L(ii+R_row1,jj)+edge(ii+R_row1,jj+2)*L(ii+R_row1,jj+2))*Na;
            Jaco_i(E_R2(ii,jj), E_R2(ii,jj+2)) = Jaco_i(E_R2(ii,jj), E_R2(ii,jj+2))+ esi*(edge(ii+R_row1,jj+2)/L(ii+R_row1,jj+2));
            Jaco_i(E_R2(ii,jj), E_R2(ii,jj+1)) = Jaco_i(E_R2(ii,jj), E_R2(ii,jj+1))+esi*(edge(ii+R_row1,jj)/L(ii+R_row1,jj));
            Jaco_i(E_R2(ii,jj), E_R2(ii,jj)) = Jaco_i(E_R2(ii,jj), E_R2(ii,jj))-esi*(edge(ii+R_row1,jj+2)/L(ii+R_row1,jj+2))-esi*(edge(ii+R_row1,jj)/L(ii+R_row1,jj));

        elseif jj ==2
            res_i(E_R2(ii,jj),1)=coeff*(edge(ii+R_row1,jj)*L(ii+R_row1,jj)+edge(ii+R_row1,jj-1)*L(ii+R_row1,jj-1))*Na;
            Jaco_i(E_R2(ii,jj), E_R2(ii,jj+1)) = Jaco_i(E_R2(ii,jj), E_R2(ii,jj+1))+esi*(edge(ii+R_row1,jj)/L(ii+R_row1,jj));
            Jaco_i(E_R2(ii,jj), E_R2(ii,jj)) =  Jaco_i(E_R2(ii,jj), E_R2(ii,jj))-esi*(edge(ii+R_row1,jj)/L(ii+R_row1,jj))-esi*(edge(ii+R_row1,jj-1)/L(ii+R_row1,jj-1));
            Jaco_i(E_R2(ii,jj), E_R2(ii,jj-1)) = Jaco_i(E_R2(ii,jj), E_R2(ii,jj-1))+esi*(edge(ii+R_row1,jj-1)/L(ii+R_row1,jj-1));
        elseif jj == 3
            res_i(E_R2(ii,jj),1)=coeff*(edge(ii+R_row1,jj)*L(ii+R_row1,jj)+edge(ii+R_row1,jj-1)*L(ii+R_row1,jj-1))*Na;
            Jaco_i(E_R2(ii,jj), E_R2(ii,jj-2)) = Jaco_i(E_R2(ii,jj), E_R2(ii,jj-2))+esi*(edge(ii+R_row1,jj)/L(ii+R_row1,jj));
            Jaco_i(E_R2(ii,jj), E_R2(ii,jj-1)) = Jaco_i(E_R2(ii,jj), E_R2(ii,jj-1))+esi*(edge(ii+R_row1,jj-1)/L(ii+R_row1,jj-1));
            Jaco_i(E_R2(ii,jj), E_R2(ii,jj)) = Jaco_i(E_R2(ii,jj), E_R2(ii,jj))-esi*(edge(ii+R_row1,jj)/L(ii+R_row1,jj))-esi*(edge(ii+R_row1,jj-1)/L(ii+R_row1,jj-1));
        end
    end
end

for ii=1:R_row3
    for jj=1:3

        if jj==1
            Jaco_i(E_R3(ii,jj), E_R3(ii,jj+2)) = Jaco_i(E_R3(ii,jj), E_R3(ii,jj+2))+ eox*(edge(ii+R_row2+R_row1,jj+2)/L(ii+R_row2+R_row1,jj+2));
            Jaco_i(E_R3(ii,jj), E_R3(ii,jj+1)) = Jaco_i(E_R3(ii,jj), E_R3(ii,jj+1))+eox*(edge(ii+R_row2+R_row1,jj)/L(ii+R_row2+R_row1,jj));
            Jaco_i(E_R3(ii,jj), E_R3(ii,jj)) = Jaco_i(E_R3(ii,jj), E_R3(ii,jj))-eox*(edge(ii+R_row2+R_row1,jj+2)/L(ii+R_row2+R_row1,jj+2))-eox*(edge(ii+R_row2+R_row1,jj)/L(ii+R_row2+R_row1,jj));

        elseif jj ==2
            Jaco_i(E_R3(ii,jj), E_R3(ii,jj+1)) = Jaco_i(E_R3(ii,jj), E_R3(ii,jj+1))+eox*(edge(ii+R_row2+R_row1,jj)/L(ii+R_row2+R_row1,jj));
            Jaco_i(E_R3(ii,jj), E_R3(ii,jj)) =  Jaco_i(E_R3(ii,jj), E_R3(ii,jj))-eox*(edge(ii+R_row2+R_row1,jj)/L(ii+R_row2+R_row1,jj))-eox*(edge(ii+R_row2+R_row1,jj-1)/L(ii+R_row2+R_row1,jj-1));
            Jaco_i(E_R3(ii,jj), E_R3(ii,jj-1)) = Jaco_i(E_R3(ii,jj), E_R3(ii,jj-1))+eox*(edge(ii+R_row2+R_row1,jj-1)/L(ii+R_row2+R_row1,jj-1));

        elseif jj == 3
            Jaco_i(E_R3(ii,jj), E_R3(ii,jj-2)) = Jaco_i(E_R3(ii,jj), E_R3(ii,jj-2))+eox*(edge(ii+R_row2+R_row1,jj)/L(ii+R_row2+R_row1,jj));
            Jaco_i(E_R3(ii,jj), E_R3(ii,jj-1)) = Jaco_i(E_R3(ii,jj), E_R3(ii,jj-1))+eox*(edge(ii+R_row2+R_row1,jj-1)/L(ii+R_row2+R_row1,jj-1));
            Jaco_i(E_R3(ii,jj), E_R3(ii,jj)) = Jaco_i(E_R3(ii,jj), E_R3(ii,jj))-eox*(edge(ii+R_row2+R_row1,jj)/L(ii+R_row2+R_row1,jj))-eox*(edge(ii+R_row2+R_row1,jj-1)/L(ii+R_row2+R_row1,jj-1));
        end
    end
end

%%%%% Dirichlet Boundary Condition %%%%%
for n=1:C_col
    for ii =1: C_row
        if ii <= 2
            Jaco_i(Contact(ii,n),:) =0;
            Jaco_i(Contact(ii,n),Contact(ii,n)) = 1;
            res_i(Contact(ii,n),1)=0;
        
         elseif ii > 2 && ii <= 4
            Jaco_i(Contact(ii,n),:) =0;
            Jaco_i(Contact(ii,n),Contact(ii,n)) = 1;
            res_i(Contact(ii,n),1)=1;
 
        end
    end
end

%%%% Result %%%%
phi = Jaco_i\res_i;

elec(7:1:21,1) = nint*exp(phi(7:1:21,1)/Thermal_V);
hole(7:1:21,1) = nint*exp(-phi(7:1:21,1)/Thermal_V);

%%%%% Nonlinear Poisson (Fully-Coupled) %%%%%%

for Newton=1:1

    Jaco1=zeros(index, index);
    res1=zeros(index,1);

    Jaco2=zeros(3*index,3*index);
    res2=zeros(3*index,1);

    Jaco3=zeros(index, index);
    res3=zeros(index,1);

    for ii=1:R_row1
        for jj=1:3

            if jj==1
                res1(E_R1(ii,jj),1) = res1(E_R1(ii,jj),1) + eox*(phi(E_R1(ii,jj+1),1)-phi(E_R1(ii,jj),1))*(edge(ii,jj)/L(ii,jj)) + eox*(phi(E_R1(ii,jj+2),1)-phi(E_R1(ii,jj),1))*(edge(ii,jj+2)/L(ii,jj+2));
                Jaco1(E_R1(ii,jj), E_R1(ii,jj+2)) = Jaco1(E_R1(ii,jj), E_R1(ii,jj+2))+ eox*(edge(ii,jj+2)/L(ii,jj+2));
                Jaco1(E_R1(ii,jj), E_R1(ii,jj+1)) = Jaco1(E_R1(ii,jj), E_R1(ii,jj+1))+eox*(edge(ii,jj)/L(ii,jj));
                Jaco1(E_R1(ii,jj), E_R1(ii,jj)) = Jaco1(E_R1(ii,jj), E_R1(ii,jj))-eox*(edge(ii,jj+2)/L(ii,jj+2))-eox*(edge(ii,jj)/L(ii,jj));

            elseif jj ==2
                res1(E_R1(ii,jj),1) = res1(E_R1(ii,jj),1) + eox*(phi(E_R1(ii,jj-1),1)-phi(E_R1(ii,jj),1))*(edge(ii,jj-1)/L(ii,jj-1)) + eox*(phi(E_R1(ii,jj+1),1)-phi(E_R1(ii,jj),1))*(edge(ii,jj)/L(ii,jj));
                Jaco1(E_R1(ii,jj), E_R1(ii,jj+1)) = Jaco1(E_R1(ii,jj), E_R1(ii,jj+1))+eox*(edge(ii,jj)/L(ii,jj));
                Jaco1(E_R1(ii,jj), E_R1(ii,jj)) =  Jaco1(E_R1(ii,jj), E_R1(ii,jj))-eox*(edge(ii,jj)/L(ii,jj))-eox*(edge(ii,jj-1)/L(ii,jj-1));
                Jaco1(E_R1(ii,jj), E_R1(ii,jj-1)) = Jaco1(E_R1(ii,jj), E_R1(ii,jj-1))+eox*(edge(ii,jj-1)/L(ii,jj-1));

            elseif jj == 3
                res1(E_R1(ii,jj),1) = res1(E_R1(ii,jj),1) + eox*(phi(E_R1(ii,jj-2),1)-phi(E_R1(ii,jj),1))*(edge(ii,jj)/L(ii,jj)) + eox*(phi(E_R1(ii,jj-1),1)-phi(E_R1(ii,jj),1))*(edge(ii,jj-1)/L(ii,jj-1));
                Jaco1(E_R1(ii,jj), E_R1(ii,jj-2)) = Jaco1(E_R1(ii,jj), E_R1(ii,jj-2))+eox*(edge(ii,jj)/L(ii,jj));
                Jaco1(E_R1(ii,jj), E_R1(ii,jj-1)) = Jaco1(E_R1(ii,jj), E_R1(ii,jj-1))+eox*(edge(ii,jj-1)/L(ii,jj-1));
                Jaco1(E_R1(ii,jj), E_R1(ii,jj)) = Jaco1(E_R1(ii,jj), E_R1(ii,jj))-eox*(edge(ii,jj)/L(ii,jj))-eox*(edge(ii,jj-1)/L(ii,jj-1));
            end
        end
    end

    for ii=1:R_row2
        for jj=1:3

            if jj==1

                %potential
                res2(3*E_R2(ii,jj)-2,1) = res2(3*E_R2(ii,jj)-2,1) + esi*(phi(E_R2(ii,jj+1),1)-phi(E_R2(ii,jj),1))*(edge(ii+R_row1,jj)/L(ii+R_row1,jj)) + esi*(phi(E_R2(ii,jj+2),1)-phi(E_R2(ii,jj),1))*(edge(ii+R_row1,jj+2)/L(ii+R_row1,jj+2));
                res2(3*E_R2(ii,jj)-2,1) = res2(3*E_R2(ii,jj)-2,1) + coeff*(edge(ii+R_row1,jj)*L(ii+R_row1,jj)+edge(ii+R_row1,jj+2)*L(ii+R_row1,jj+2))*(-Na-elec(E_R2(ii,jj),1)+hole(E_R2(ii,jj),1));

                Jaco2(3*E_R2(ii,jj)-2, 3*E_R2(ii,jj+2)-2) = Jaco2(3*E_R2(ii,jj)-2, 3*(E_R2(ii,jj+2))-2)+ esi*(edge(ii+R_row1,jj+2)/L(ii+R_row1,jj+2));
                Jaco2(3*E_R2(ii,jj)-2, 3*E_R2(ii,jj+1)-2) = Jaco2(3*E_R2(ii,jj)-2, 3*(E_R2(ii,jj+1))-2) + esi*(edge(ii+R_row1,jj)/L(ii+R_row1,jj));
                Jaco2(3*E_R2(ii,jj)-2, 3*E_R2(ii,jj)-2) = Jaco2(3*E_R2(ii,jj)-2, 3*E_R2(ii,jj)-2)-esi*(edge(ii+R_row1,jj+2)/L(ii+R_row1,jj+2))-esi*(edge(ii+R_row1,jj)/L(ii+R_row1,jj));

                Jaco2(3*E_R2(ii,jj)-2, 3*E_R2(ii,jj)-1) = Jaco2(3*E_R2(ii,jj)-2, 3*E_R2(ii,jj)-1) - coeff*(edge(ii+R_row1,jj)*L(ii+R_row1,jj)+edge(ii+R_row1,jj+2)*L(ii+R_row1,jj+2));
                Jaco2(3*E_R2(ii,jj)-2, 3*E_R2(ii,jj)) = Jaco2(3*E_R2(ii,jj)-2, 3*E_R2(ii,jj)) + coeff*(edge(ii+R_row1,jj)*L(ii+R_row1,jj)+edge(ii+R_row1,jj+2)*L(ii+R_row1,jj+2));

                %elec
                res2(3*E_R2(ii,jj)-1,1) = res2(3*E_R2(ii,jj)-1,1) + elec(E_R2(ii,jj),1) - nint*exp(phi(E_R2(ii,jj),1)/Thermal_V);

                Jaco2(3*E_R2(ii,jj)-1, 3*E_R2(ii,jj)-1) = Jaco2(3*E_R2(ii,jj)-1, 3*E_R2(ii,jj)-1) + 1;
                Jaco2(3*E_R2(ii,jj)-1, 3*E_R2(ii,jj)-1-1) = Jaco2(3*E_R2(ii,jj)-1, 3*E_R2(ii,jj)-1-1) - (nint/Thermal_V)*exp(phi(E_R2(ii,jj),1)/Thermal_V);

                %hole
                res2(3*E_R2(ii,jj),1) = res2(3*E_R2(ii,jj),1) + hole(E_R2(ii,jj),1)-nint*exp(-phi(E_R2(ii,jj),1)/Thermal_V);

                Jaco2(3*E_R2(ii,jj), 3*E_R2(ii,jj)) = Jaco2(3*E_R2(ii,jj), 3*E_R2(ii,jj)) + 1;
                Jaco2(3*E_R2(ii,jj), 3*E_R2(ii,jj)-1-1) = Jaco2(3*E_R2(ii,jj), 3*E_R2(ii,jj)-1-1) + (nint/Thermal_V)*exp(-phi(E_R2(ii,jj),1)/Thermal_V);


            elseif jj ==2

                %potential
                res2(3*E_R2(ii,jj)-2,1) = res2(3*E_R2(ii,jj)-2,1) + esi*(phi(E_R2(ii,jj-1),1)-phi(E_R2(ii,jj),1))*(edge(ii+R_row1,jj-1)/L(ii+R_row1,jj-1)) + esi*(phi(E_R2(ii,jj+1),1)-phi(E_R2(ii,jj),1))*(edge(ii+R_row1,jj)/L(ii+R_row1,jj));
                res2(3*E_R2(ii,jj)-2,1) = res2(3*E_R2(ii,jj)-2,1) + coeff*(edge(ii+R_row1,jj)*L(ii+R_row1,jj)+edge(ii+R_row1,jj-1)*L(ii+R_row1,jj-1))*(-Na-elec(E_R2(ii,jj),1)+hole(E_R2(ii,jj),1));

                Jaco2(3*E_R2(ii,jj)-2, 3*E_R2(ii,jj+1)-2) = Jaco2(3*E_R2(ii,jj)-2, 3*E_R2(ii,jj+1)-2)+esi*(edge(ii+R_row1,jj)/L(ii+R_row1,jj));
                Jaco2(3*E_R2(ii,jj)-2, 3*E_R2(ii,jj-1)-2) = Jaco2(3*E_R2(ii,jj)-2, 3*E_R2(ii,jj-1)-2)+esi*(edge(ii+R_row1,jj-1)/L(ii+R_row1,jj-1));
                Jaco2(3*E_R2(ii,jj)-2, 3*E_R2(ii,jj)-2) =  Jaco2(3*E_R2(ii,jj)-2, 3*E_R2(ii,jj)-2)-esi*(edge(ii+R_row1,jj)/L(ii+R_row1,jj))-esi*(edge(ii+R_row1,jj-1)/L(ii+R_row1,jj-1));

                Jaco2(3*E_R2(ii,jj)-2, 3*E_R2(ii,jj)-1) =  Jaco2(3*E_R2(ii,jj)-2, 3*E_R2(ii,jj)-1) - coeff*(edge(ii+R_row1,jj)*L(ii+R_row1,jj)+edge(ii+R_row1,jj-1)*L(ii+R_row1,jj-1));
                Jaco2(3*E_R2(ii,jj)-2, 3*E_R2(ii,jj)) =  Jaco2(3*E_R2(ii,jj)-2, 3*E_R2(ii,jj)) + coeff*(edge(ii+R_row1,jj)*L(ii+R_row1,jj)+edge(ii+R_row1,jj-1)*L(ii+R_row1,jj-1));

                %elec
                res2(3*E_R2(ii,jj)-1,1) = res2(3*E_R2(ii,jj)-1,1) + elec(E_R2(ii,jj),1) - nint*exp(phi(E_R2(ii,jj),1)/Thermal_V);

                Jaco2(3*E_R2(ii,jj)-1, 3*E_R2(ii,jj)-1) = Jaco2(3*E_R2(ii,jj)-1, 3*E_R2(ii,jj)-1) + 1;
                Jaco2(3*E_R2(ii,jj)-1, 3*E_R2(ii,jj)-1-1) = Jaco2(3*E_R2(ii,jj)-1, 3*E_R2(ii,jj)-1-1) - (nint/Thermal_V)*exp(phi(E_R2(ii,jj),1)/Thermal_V);

                %hole
                res2(3*E_R2(ii,jj),1) = res2(3*E_R2(ii,jj),1) + hole(E_R2(ii,jj),1)-nint*exp(-phi(E_R2(ii,jj),1)/Thermal_V);

                Jaco2(3*E_R2(ii,jj), 3*E_R2(ii,jj)) = Jaco2(3*E_R2(ii,jj), 3*E_R2(ii,jj)) + 1;
                Jaco2(3*E_R2(ii,jj), 3*E_R2(ii,jj)-1-1) = Jaco2(3*E_R2(ii,jj), 3*E_R2(ii,jj)-1-1) + (nint/Thermal_V)*exp(-phi(E_R2(ii,jj),1)/Thermal_V);


            elseif jj == 3
 
                %potential
                res2(3*E_R2(ii,jj)-2,1) = res2(3*E_R2(ii,jj)-2,1) + esi*(phi(E_R2(ii,jj-2),1)-phi(E_R2(ii,jj),1))*(edge(ii+R_row1,jj)/L(ii+R_row1,jj)) + esi*(phi(E_R2(ii,jj-1),1)-phi(E_R2(ii,jj),1))*(edge(ii+R_row1,jj-1)/L(ii+R_row1,jj-1));
                res2(3*E_R2(ii,jj)-2,1) = res2(3*E_R2(ii,jj)-2,1) + coeff*(edge(ii+R_row1,jj)*L(ii+R_row1,jj)+edge(ii+R_row1,jj-1)*L(ii+R_row1,jj-1))*(-Na-elec(E_R2(ii,jj),1)+hole(E_R2(ii,jj),1));

                Jaco2(3*E_R2(ii,jj)-2, 3*E_R2(ii,jj-2)-2) = Jaco2(3*E_R2(ii,jj)-2, 3*E_R2(ii,jj-2)-2)+esi*(edge(ii+R_row1,jj)/L(ii+R_row1,jj));
                Jaco2(3*E_R2(ii,jj)-2, 3*E_R2(ii,jj-1)-2) = Jaco2(3*E_R2(ii,jj)-2, 3*E_R2(ii,jj-1)-2)+esi*(edge(ii+R_row1,jj-1)/L(ii+R_row1,jj-1));
                Jaco2(3*E_R2(ii,jj)-2, 3*E_R2(ii,jj)-2) = Jaco2(3*E_R2(ii,jj)-2, 3*E_R2(ii,jj)-2)-esi*(edge(ii+R_row1,jj)/L(ii+R_row1,jj))-esi*(edge(ii+R_row1,jj-1)/L(ii+R_row1,jj-1));

                Jaco2(3*E_R2(ii,jj)-2, 3*E_R2(ii,jj)-1) =  Jaco2(3*E_R2(ii,jj)-2, 3*E_R2(ii,jj)-1) - coeff*(edge(ii+R_row1,jj)*L(ii+R_row1,jj)+edge(ii+R_row1,jj-1)*L(ii+R_row1,jj-1));
                Jaco2(3*E_R2(ii,jj)-2, 3*E_R2(ii,jj)) =  Jaco2(3*E_R2(ii,jj)-2, 3*E_R2(ii,jj)) + coeff*(edge(ii+R_row1,jj)*L(ii+R_row1,jj)+edge(ii+R_row1,jj-1)*L(ii+R_row1,jj-1));
 
                %elec
                res2(3*E_R2(ii,jj)-1,1) = res2(3*E_R2(ii,jj)-1,1) + elec(E_R2(ii,jj),1) - nint*exp(phi(E_R2(ii,jj),1)/Thermal_V);

                Jaco2(3*E_R2(ii,jj)-1, 3*E_R2(ii,jj)-1) = Jaco2(3*E_R2(ii,jj)-1, 3*E_R2(ii,jj)-1) + 1;
                Jaco2(3*E_R2(ii,jj)-1, 3*E_R2(ii,jj)-1-1) = Jaco2(3*E_R2(ii,jj)-1, 3*E_R2(ii,jj)-1-1) - (nint/Thermal_V)*exp(phi(E_R2(ii,jj),1)/Thermal_V);
            
                %hole
                res2(3*E_R2(ii,jj),1) = res2(3*E_R2(ii,jj),1) + hole(E_R2(ii,jj),1)-nint*exp(-phi(E_R2(ii,jj),1)/Thermal_V);

                Jaco2(3*E_R2(ii,jj), 3*E_R2(ii,jj)) = Jaco2(3*E_R2(ii,jj), 3*E_R2(ii,jj)) + 1;
                Jaco2(3*E_R2(ii,jj), 3*E_R2(ii,jj)-1-1) = Jaco2(3*E_R2(ii,jj), 3*E_R2(ii,jj)-1-1) + (nint/Thermal_V)*exp(-phi(E_R2(ii,jj),1)/Thermal_V);
            end
        end
    end

    for ii=1:R_row3
        for jj=1:3

            if jj==1
                res3(E_R3(ii,jj),1) = res3(E_R3(ii,jj),1) + eox*(phi(E_R3(ii,jj+1),1)-phi(E_R3(ii,jj),1))*(edge(ii+R_row2+R_row1,jj)/L(ii+R_row2+R_row1,jj)) + eox*(phi(E_R3(ii,jj+2),1)-phi(E_R3(ii,jj),1))*(edge(ii+R_row2+R_row1,jj+2)/L(ii+R_row2+R_row1,jj+2));
                Jaco3(E_R3(ii,jj), E_R3(ii,jj+2)) = Jaco3(E_R3(ii,jj), E_R3(ii,jj+2))+ eox*(edge(ii+R_row2+R_row1,jj+2)/L(ii+R_row2+R_row1,jj+2));
                Jaco3(E_R3(ii,jj), E_R3(ii,jj+1)) = Jaco3(E_R3(ii,jj), E_R3(ii,jj+1))+eox*(edge(ii+R_row2+R_row1,jj)/L(ii+R_row2+R_row1,jj));
                Jaco3(E_R3(ii,jj), E_R3(ii,jj)) = Jaco3(E_R3(ii,jj), E_R3(ii,jj))-eox*(edge(ii+R_row2+R_row1,jj+2)/L(ii+R_row2+R_row1,jj+2))-eox*(edge(ii+R_row2+R_row1,jj)/L(ii+R_row2+R_row1,jj));

            elseif jj ==2
                res3(E_R3(ii,jj),1) = res3(E_R3(ii,jj),1) + eox*(phi(E_R3(ii,jj-1),1)-phi(E_R3(ii,jj),1))*(edge(ii+R_row2+R_row1,jj-1)/L(ii+R_row2+R_row1,jj-1)) + eox*(phi(E_R3(ii,jj+1),1)-phi(E_R3(ii,jj),1))*(edge(ii+R_row2+R_row1,jj)/L(ii+R_row2+R_row1,jj));
                Jaco3(E_R3(ii,jj), E_R3(ii,jj+1)) = Jaco3(E_R3(ii,jj), E_R3(ii,jj+1))+eox*(edge(ii+R_row2+R_row1,jj)/L(ii+R_row2+R_row1,jj));
                Jaco3(E_R3(ii,jj), E_R3(ii,jj)) =  Jaco3(E_R3(ii,jj), E_R3(ii,jj))-eox*(edge(ii+R_row2+R_row1,jj)/L(ii+R_row2+R_row1,jj))-eox*(edge(ii+R_row2+R_row1,jj-1)/L(ii+R_row2+R_row1,jj-1));
                Jaco3(E_R3(ii,jj), E_R3(ii,jj-1)) = Jaco3(E_R3(ii,jj), E_R3(ii,jj-1))+eox*(edge(ii+R_row2+R_row1,jj-1)/L(ii+R_row2+R_row1,jj-1));

            elseif jj == 3
                res3(E_R3(ii,jj),1) = res3(E_R3(ii,jj),1) + eox*(phi(E_R3(ii,jj-2),1)-phi(E_R3(ii,jj),1))*(edge(ii+R_row2+R_row1,jj)/L(ii+R_row2+R_row1,jj)) + eox*(phi(E_R3(ii,jj-1),1)-phi(E_R3(ii,jj),1))*(edge(ii+R_row2+R_row1,jj-1)/L(ii+R_row2+R_row1,jj-1));
                Jaco3(E_R3(ii,jj), E_R3(ii,jj-2)) = Jaco3(E_R3(ii,jj), E_R3(ii,jj-2))+eox*(edge(ii+R_row2+R_row1,jj)/L(ii+R_row2+R_row1,jj));
                Jaco3(E_R3(ii,jj), E_R3(ii,jj-1)) = Jaco3(E_R3(ii,jj), E_R3(ii,jj-1))+eox*(edge(ii+R_row2+R_row1,jj-1)/L(ii+R_row2+R_row1,jj-1));
                Jaco3(E_R3(ii,jj), E_R3(ii,jj)) = Jaco3(E_R3(ii,jj), E_R3(ii,jj))-eox*(edge(ii+R_row2+R_row1,jj)/L(ii+R_row2+R_row1,jj))-eox*(edge(ii+R_row2+R_row1,jj-1)/L(ii+R_row2+R_row1,jj-1));
            end
        end
    end

    %%%%% Dirichlet Boundary Condition %%%%%
    for n=1:C_col
        for i =1: C_row
            if i <= 2
                Jaco1(Contact(i,n),:) =0;
                Jaco1(Contact(i,n),Contact(i,n)) = 1;
                res1(Contact(i,n),1)=phi(Contact(i,n),1)-0;

            elseif i > 2 && i <= 4
                Jaco3(Contact(i,n),:) =0;
                Jaco3(Contact(i,n),Contact(i,n)) = 1;
                res3(Contact(i,n),1)= phi(Contact(i,n),1)-1;

            end
        end
    end

    %%%%%%%%% Jacobian Reindexing %%%%%%%%%

    Jaco=zeros(63, 63);
    res=zeros(63,1);

    %%%% trimming the null row %%%%

    %%%% Jacobian %%%%

    Jaco1(7,:)=[];
    Jaco1(:,7)=[];

    Jaco3(21,:)=[];
    Jaco3(:,21)=[];

    Jaco(1:9,1:9) = Jaco(1:9,1:9)+Jaco1(1:9,1:9);
    Jaco(10:54,10:54) = Jaco(10:54,10:54)+Jaco2(19:63,19:63);
    Jaco(55:63,55:63) = Jaco(55:63,55:63)+Jaco3(18:26,18:26);

    %Boundary Condition

    Jaco(13,:)=Jaco(13,:)+Jaco(7,:);
    Jaco(7,:)=0; Jaco(7,7)=1; Jaco(7,13)=-1;
    Jaco(16,:)=Jaco(16,:)+Jaco(8,:);
    Jaco(8,:)=0; Jaco(8,8)=1; Jaco(8,16)=-1;
    Jaco(19,:)=Jaco(19,:)+Jaco(9,:);
    Jaco(9,:)=0; Jaco(9,9)=1; Jaco(9,19)=-1;

    Jaco(43,:)=Jaco(43,:)+Jaco(55,:);
    Jaco(55,:)=0; Jaco(55,55)=1; Jaco(55,43)=-1;
    Jaco(46,:)=Jaco(46,:)+Jaco(56,:);
    Jaco(56,:)=0; Jaco(56,56)=1; Jaco(56,46)=-1;
    Jaco(49,:)=Jaco(49,:)+Jaco(57,:);
    Jaco(57,:)=0; Jaco(57,57)=1; Jaco(57,49)=-1;


    %%%% Residue %%%%%

    res1(7,:)=[];
    res3(21,:)=[];

    res(1:9,1) = res(1:9,1)+res1(1:9,1);
    res(10:54,1) = res(10:54,1)+res2(19:63,1);
    res(55:63,1) = res(55:63,1)+res3(18:26,1);

    %Boundary Condition

    res(13,1)=res(13,1)+res(7,1);
    res(7,1)=0; res(7,1)=1;
    res(16,1)=res(16,1)+res(8,1);
    res(8,1)=0; res(8,1)=1;
    res(19,1)=res(19,1)+res(9,1);
    res(9,1)=0; res(9,1)=1;

    res(43,1)=res(43,1)+res(55,1);
    res(55,1)=0; res(55,1)=1;
    res(46,1)=res(46,1)+res(56,1);
    res(56,1)=0; res(56,1)=1;
    res(49,1)=res(49,1)+res(57,1);
    res(57,1)=0; res(57,1)=1;



    delphi_coupled(:,Newton) = Jaco \ (-res);
%     phi = phi+delphi_coupled(:,Newton);
end
