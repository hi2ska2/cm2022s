clear; close all; clc;

% due to : 2022.03.29 (Tue) / HW7
% Contstruct a solution vertor, x 

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
[R3_row,C_col]=size(Contact);

index = V_row;

%%%%%%%% Specify Region 1,2 %%%%%%%%%%

fileID2 = fopen('Element_region1.txt');
element_Region1 = textscan(fileID2,'%f%f%f','Delimiter','\t');

Element_Region1 = horzcat(element_Region1{1},element_Region1{2},element_Region1{3});
fclose(fileID2);

fileID3 = fopen('Element_region2.txt');
element_Region2 = textscan(fileID3,'%f%f%f','Delimiter','\t');

Element_Region2 = horzcat(element_Region2{1},element_Region2{2},element_Region2{3});
fclose(fileID3);

fileID4 = fopen('Element_region3.txt');
element_Region3 = textscan(fileID4,'%f%f%f','Delimiter','\t');

Element_Region3 = horzcat(element_Region3{1},element_Region3{2},element_Region3{3});
fclose(fileID4);

R_row1 = size(Element_Region1,1);
R_row2 = size(Element_Region2,1);
R_row3 = size(Element_Region3,1);

esi = 11.7; eox = 3.9;

%%%%%% Geometry information %%%%%%

L=zeros(E_row,3);

for ii = 1:R_row1
    for jj = 1:3
        if jj == 1 % between 1 and 2
            L(ii,jj) = sqrt(abs((Vertex(Element_Region1(ii,jj+1),1)-Vertex(Element_Region1(ii,jj),1))^2+(Vertex(Element_Region1(ii,jj+1),2)-Vertex(Element_Region1(ii,jj),2))^2));
        elseif jj == 2 % between 2 and 3
            L(ii,jj) = sqrt(abs((Vertex(Element_Region1(ii,jj+1),1)-Vertex(Element_Region1(ii,jj),1))^2+(Vertex(Element_Region1(ii,jj+1),2)-Vertex(Element_Region1(ii,jj),2))^2));
        elseif jj ==3 % between 3 and 1
            L(ii,jj) = sqrt(abs((Vertex(Element_Region1(ii,jj),1)-Vertex(Element_Region1(ii,jj-2),1))^2+(Vertex(Element_Region1(ii,jj),2)-Vertex(Element_Region1(ii,jj-2),2))^2));
        end
    end
end

for ii = 1:R_row2
    for jj = 1:3
        if jj == 1 % between 1 and 2
            L(ii+R_row1,jj) = sqrt(abs((Vertex(Element_Region2(ii,jj+1),1)-Vertex(Element_Region2(ii,jj),1))^2+(Vertex(Element_Region2(ii,jj+1),2)-Vertex(Element_Region2(ii,jj),2))^2));
        elseif jj == 2 % between 2 and 3
            L(ii+R_row1,jj) = sqrt(abs((Vertex(Element_Region2(ii,jj+1),1)-Vertex(Element_Region2(ii,jj),1))^2+(Vertex(Element_Region2(ii,jj+1),2)-Vertex(Element_Region2(ii,jj),2))^2));
        elseif jj ==3 % between 3 and 1
            L(ii+R_row1,jj) = sqrt(abs((Vertex(Element_Region2(ii,jj),1)-Vertex(Element_Region2(ii,jj-2),1))^2+(Vertex(Element_Region2(ii,jj),2)-Vertex(Element_Region2(ii,jj-2),2))^2));
        end
    end
end

for ii = 1:R_row3
    for jj = 1:3
        if jj == 1 % between 1 and 2
            L(ii+R_row2+R_row1,jj) = sqrt(abs((Vertex(Element_Region3(ii,jj+1),1)-Vertex(Element_Region3(ii,jj),1))^2+(Vertex(Element_Region3(ii,jj+1),2)-Vertex(Element_Region3(ii,jj),2))^2));
        elseif jj == 2 % between 2 and 3
            L(ii+R_row2+R_row1,jj) = sqrt(abs((Vertex(Element_Region3(ii,jj+1),1)-Vertex(Element_Region3(ii,jj),1))^2+(Vertex(Element_Region3(ii,jj+1),2)-Vertex(Element_Region3(ii,jj),2))^2));
        elseif jj ==3 % between 3 and 1
            L(ii+R_row2+R_row1,jj) = sqrt(abs((Vertex(Element_Region3(ii,jj),1)-Vertex(Element_Region3(ii,jj-2),1))^2+(Vertex(Element_Region3(ii,jj),2)-Vertex(Element_Region3(ii,jj-2),2))^2));
        end
    end
end

% Vector
Area=zeros(E_row,1);
R=zeros(E_row,1);

for ii = 1:R_row1
    a12 = [(Vertex(Element_Region1(ii,2),1)-Vertex(Element_Region1(ii,1),1)),(Vertex(Element_Region1(ii,2),2)-Vertex(Element_Region1(ii,1),2))];
    a13 = [(Vertex(Element_Region1(ii,3),1)-Vertex(Element_Region1(ii,1),1)),(Vertex(Element_Region1(ii,3),2)-Vertex(Element_Region1(ii,1),2))];
    a23 = [(Vertex(Element_Region1(ii,3),1)-Vertex(Element_Region1(ii,2),1)),(Vertex(Element_Region1(ii,3),2)-Vertex(Element_Region1(ii,2),2))];
    a21=-a12;
    a31=-a13;
    a32=-a23;

    Area(ii,1) = (1/2)*abs(a12(1,1)*a13(1,2)-a12(1,2)*a13(1,1));
    R(ii,1) = (L(ii,1)*L(ii,2)*L(ii,3))/(4*Area(ii,1));
end

for ii = 1:R_row2
    a12 = [(Vertex(Element_Region2(ii,2),1)-Vertex(Element_Region2(ii,1),1)),(Vertex(Element_Region2(ii,2),2)-Vertex(Element_Region2(ii,1),2))];
    a13 = [(Vertex(Element_Region2(ii,3),1)-Vertex(Element_Region2(ii,1),1)),(Vertex(Element_Region2(ii,3),2)-Vertex(Element_Region2(ii,1),2))];
    a23 = [(Vertex(Element_Region2(ii,3),1)-Vertex(Element_Region2(ii,2),1)),(Vertex(Element_Region2(ii,3),2)-Vertex(Element_Region2(ii,2),2))];
    a21=-a12;
    a31=-a13;
    a32=-a23;

    Area(ii+R_row1,1) = (1/2)*abs(a12(1,1)*a13(1,2)-a12(1,2)*a13(1,1));
    R(ii+R_row1,1) = (L(ii+R_row1,1)*L(ii+R_row1,2)*L(ii+R_row1,3))/(4*Area(ii+R_row1,1));
end

for ii = 1:R_row3
    a12 = [(Vertex(Element_Region3(ii,2),1)-Vertex(Element_Region3(ii,1),1)),(Vertex(Element_Region3(ii,2),2)-Vertex(Element_Region3(ii,1),2))];
    a13 = [(Vertex(Element_Region3(ii,3),1)-Vertex(Element_Region3(ii,1),1)),(Vertex(Element_Region3(ii,3),2)-Vertex(Element_Region3(ii,1),2))];
    a23 = [(Vertex(Element_Region3(ii,3),1)-Vertex(Element_Region3(ii,2),1)),(Vertex(Element_Region3(ii,3),2)-Vertex(Element_Region3(ii,2),2))];
    a21=-a12;
    a31=-a13;
    a32=-a23;

    Area(ii+R_row2+R_row1,1) = (1/2)*abs(a12(1,1)*a13(1,2)-a12(1,2)*a13(1,1));
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

%%%% Potential %%%%
Jaco = zeros(index,index);
res = zeros(index,1);

for ii=1:R_row1
    for jj=1:3

        if jj==1
            Jaco(Element_Region1(ii,jj), Element_Region1(ii,jj+2)) = Jaco(Element_Region1(ii,jj), Element_Region1(ii,jj+2))+ eox*(edge(ii,jj+2)/L(ii,jj+2));
            Jaco(Element_Region1(ii,jj), Element_Region1(ii,jj+1)) = Jaco(Element_Region1(ii,jj), Element_Region1(ii,jj+1))+eox*(edge(ii,jj)/L(ii,jj));
            Jaco(Element_Region1(ii,jj), Element_Region1(ii,jj)) = Jaco(Element_Region1(ii,jj), Element_Region1(ii,jj))-eox*(edge(ii,jj+2)/L(ii,jj+2))-eox*(edge(ii,jj)/L(ii,jj));

        elseif jj ==2
            Jaco(Element_Region1(ii,jj), Element_Region1(ii,jj+1)) = Jaco(Element_Region1(ii,jj), Element_Region1(ii,jj+1))+eox*(edge(ii,jj)/L(ii,jj));
            Jaco(Element_Region1(ii,jj), Element_Region1(ii,jj)) =  Jaco(Element_Region1(ii,jj), Element_Region1(ii,jj))-eox*(edge(ii,jj)/L(ii,jj))-eox*(edge(ii,jj-1)/L(ii,jj-1));
            Jaco(Element_Region1(ii,jj), Element_Region1(ii,jj-1)) = Jaco(Element_Region1(ii,jj), Element_Region1(ii,jj-1))+eox*(edge(ii,jj-1)/L(ii,jj-1));

        elseif jj == 3
            Jaco(Element_Region1(ii,jj), Element_Region1(ii,jj-2)) = Jaco(Element_Region1(ii,jj), Element_Region1(ii,jj-2))+eox*(edge(ii,jj)/L(ii,jj));
            Jaco(Element_Region1(ii,jj), Element_Region1(ii,jj-1)) = Jaco(Element_Region1(ii,jj), Element_Region1(ii,jj-1))+eox*(edge(ii,jj-1)/L(ii,jj-1));
            Jaco(Element_Region1(ii,jj), Element_Region1(ii,jj)) = Jaco(Element_Region1(ii,jj), Element_Region1(ii,jj))-eox*(edge(ii,jj)/L(ii,jj))-eox*(edge(ii,jj-1)/L(ii,jj-1));
        end
    end
end

for ii=1:R_row2
    for jj=1:3

        if jj==1
            Jaco(Element_Region2(ii,jj), Element_Region2(ii,jj+2)) = Jaco(Element_Region2(ii,jj), Element_Region2(ii,jj+2))+ esi*(edge(ii+R_row1,jj+2)/L(ii+R_row1,jj+2));
            Jaco(Element_Region2(ii,jj), Element_Region2(ii,jj+1)) = Jaco(Element_Region2(ii,jj), Element_Region2(ii,jj+1))+esi*(edge(ii+R_row1,jj)/L(ii+R_row1,jj));
            Jaco(Element_Region2(ii,jj), Element_Region2(ii,jj)) = Jaco(Element_Region2(ii,jj), Element_Region2(ii,jj))-esi*(edge(ii+R_row1,jj+2)/L(ii+R_row1,jj+2))-esi*(edge(ii+R_row1,jj)/L(ii+R_row1,jj));

        elseif jj ==2
            Jaco(Element_Region2(ii,jj), Element_Region2(ii,jj+1)) = Jaco(Element_Region2(ii,jj), Element_Region2(ii,jj+1))+esi*(edge(ii+R_row1,jj)/L(ii+R_row1,jj));
            Jaco(Element_Region2(ii,jj), Element_Region2(ii,jj)) =  Jaco(Element_Region2(ii,jj), Element_Region2(ii,jj))-esi*(edge(ii+R_row1,jj)/L(ii+R_row1,jj))-esi*(edge(ii+R_row1,jj-1)/L(ii+R_row1,jj-1));
            Jaco(Element_Region2(ii,jj), Element_Region2(ii,jj-1)) = Jaco(Element_Region2(ii,jj), Element_Region2(ii,jj-1))+esi*(edge(ii+R_row1,jj-1)/L(ii+R_row1,jj-1));
        elseif jj == 3
            Jaco(Element_Region2(ii,jj), Element_Region2(ii,jj-2)) = Jaco(Element_Region2(ii,jj), Element_Region2(ii,jj-2))+esi*(edge(ii+R_row1,jj)/L(ii+R_row1,jj));
            Jaco(Element_Region2(ii,jj), Element_Region2(ii,jj-1)) = Jaco(Element_Region2(ii,jj), Element_Region2(ii,jj-1))+esi*(edge(ii+R_row1,jj-1)/L(ii+R_row1,jj-1));
            Jaco(Element_Region2(ii,jj), Element_Region2(ii,jj)) = Jaco(Element_Region2(ii,jj), Element_Region2(ii,jj))-esi*(edge(ii+R_row1,jj)/L(ii+R_row1,jj))-esi*(edge(ii+R_row1,jj-1)/L(ii+R_row1,jj-1));
        end
    end
end

for ii=1:R_row3
    for jj=1:3

        if jj==1
            Jaco(Element_Region3(ii,jj), Element_Region3(ii,jj+2)) = Jaco(Element_Region3(ii,jj), Element_Region3(ii,jj+2))+ eox*(edge(ii+R_row2+R_row1,jj+2)/L(ii+R_row2+R_row1,jj+2));
            Jaco(Element_Region3(ii,jj), Element_Region3(ii,jj+1)) = Jaco(Element_Region3(ii,jj), Element_Region3(ii,jj+1))+eox*(edge(ii+R_row2+R_row1,jj)/L(ii+R_row2+R_row1,jj));
            Jaco(Element_Region3(ii,jj), Element_Region3(ii,jj)) = Jaco(Element_Region3(ii,jj), Element_Region3(ii,jj))-eox*(edge(ii+R_row2+R_row1,jj+2)/L(ii+R_row2+R_row1,jj+2))-eox*(edge(ii+R_row2+R_row1,jj)/L(ii+R_row2+R_row1,jj));

        elseif jj ==2
            Jaco(Element_Region3(ii,jj), Element_Region3(ii,jj+1)) = Jaco(Element_Region3(ii,jj), Element_Region3(ii,jj+1))+eox*(edge(ii+R_row2+R_row1,jj)/L(ii+R_row2+R_row1,jj));
            Jaco(Element_Region3(ii,jj), Element_Region3(ii,jj)) =  Jaco(Element_Region3(ii,jj), Element_Region3(ii,jj))-eox*(edge(ii+R_row2+R_row1,jj)/L(ii+R_row2+R_row1,jj))-eox*(edge(ii+R_row2+R_row1,jj-1)/L(ii+R_row2+R_row1,jj-1));
            Jaco(Element_Region3(ii,jj), Element_Region3(ii,jj-1)) = Jaco(Element_Region3(ii,jj), Element_Region3(ii,jj-1))+eox*(edge(ii+R_row2+R_row1,jj-1)/L(ii+R_row2+R_row1,jj-1));

        elseif jj == 3
            Jaco(Element_Region3(ii,jj), Element_Region3(ii,jj-2)) = Jaco(Element_Region3(ii,jj), Element_Region3(ii,jj-2))+eox*(edge(ii+R_row2+R_row1,jj)/L(ii+R_row2+R_row1,jj));
            Jaco(Element_Region3(ii,jj), Element_Region3(ii,jj-1)) = Jaco(Element_Region3(ii,jj), Element_Region3(ii,jj-1))+eox*(edge(ii+R_row2+R_row1,jj-1)/L(ii+R_row2+R_row1,jj-1));
            Jaco(Element_Region3(ii,jj), Element_Region3(ii,jj)) = Jaco(Element_Region3(ii,jj), Element_Region3(ii,jj))-eox*(edge(ii+R_row2+R_row1,jj)/L(ii+R_row2+R_row1,jj))-eox*(edge(ii+R_row2+R_row1,jj-1)/L(ii+R_row2+R_row1,jj-1));
        end
    end
end

%%%%% Dirichlet Boundary Condition %%%%%

for n=1:C_col
    for ii =1: R3_row
        if ii <= 2
            Jaco(Contact(ii,n),:) =0;
            Jaco(Contact(ii,n),Contact(ii,n)) = 1;
            res(Contact(ii,n),1)=0;
        
         elseif ii > 2 && ii <= 4
            Jaco(Contact(ii,n),:) =0;
            Jaco(Contact(ii,n),Contact(ii,n)) = 1;
            res(Contact(ii,n),1)=0;

         elseif ii > 4 && ii <= 6
            Jaco(Contact(ii,n),:) =0;
            Jaco(Contact(ii,n),Contact(ii,n)) = 1;
            res(Contact(ii,n),1)=0;
         elseif ii > 6
            Jaco(Contact(ii,n),:) =0;
            Jaco(Contact(ii,n),Contact(ii,n)) = 1;
            res(Contact(ii,n),1)=1;

        end
    end
end

%%%% Result %%%%
phi = Jaco\res;

%%%%%%%%% Find the interface Edge %%%%%%%%%%%

S1=intersect(Element_Region1, Element_Region2);
j =1;

for ii = 1:R_row2
    if size(intersect(Element_Region2(ii,:),S1'),2) == 2
        Edge1(j,:) = intersect(Element_Region2(ii,:),S1');
        j = j+1;
    end
end

S2=intersect(Element_Region3, Element_Region2);
j =1;

for ii = 1:R_row2
    if size(intersect(Element_Region2(ii,:),S2'),2) == 2
        Edge2(j,:) = intersect(Element_Region2(ii,:),S2);
        j = j+1;
    end
end

patch('Faces',Element_Region1,'Vertices',Vertex,'EdgeColor','Red','FaceColor','none','LineWidth',3);
patch('Faces',Element_Region2,'Vertices',Vertex,'EdgeColor','blue','FaceColor','none','LineWidth',3);
patch('Faces',Element_Region3,'Vertices',Vertex,'EdgeColor','green','FaceColor','none','LineWidth',3);
patch('Faces',Edge1,'Vertices',Vertex,'EdgeColor','yellow','FaceColor','none','LineWidth',3);
patch('Faces',Edge2,'Vertices',Vertex,'EdgeColor','yellow','FaceColor','none','LineWidth',3);
patch('Faces',Contact,'Vertices',Vertex,'EdgeColor','Cyan','FaceColor','none','LineWidth',3);

% %%%% Vertex index for each region %%%

A = unique(Element_Region1);
B = unique(Element_Region2);
C = unique(Element_Region3);

R1_row = size(A,1); % Region 1 (Oxide)
R2_row = size(B,1); % Region 2 (Si_Sub)
R3_row = size(C,1); % Region 3 (Oxide)

t=R1_row+R2_row+R3_row;
tsi=R2_row;
tox=R1_row+R3_row;

AAA = ones(t,1); % Total variables
BBB = rand(tsi,1); % variables only Si
CCC = rand(tox,1);  % variables only Ox

%%%%%%% Check the number of variables applied to the region %%%%%%%

Nv1=2; % aaa and ccc
Nv2=2; % aaa and bbb
Nv3=2; % aaa and ccc

T=Nv1*R1_row+Nv2*R2_row+Nv3*R3_row;

%%%%%%% Contstruct a solution vertor, x %%%%%

X = zeros(T,1);

for ii = 1:T
    if ii <= R1_row
        X(Nv1*ii-1) = AAA(ii);
        X(Nv1*ii) = CCC(ii);
    elseif  ii <= t-tox
        X(Nv2*ii-1) = AAA(ii-R1_row);
        X(Nv2*ii) = BBB(ii-R1_row);
    elseif ii <= t
        X(Nv3*ii-1) = AAA(ii-R2_row);
        X(Nv3*ii) = CCC(ii-R2_row);
    end
end

%%%%%%% Input Variable, Region, Vertex, Come out X index %%%%%%%

str = input('Enter Variable : ' ,'s');
r = input('Enter Region number (we have only 3 region) : ' );
v = input('Enter Vertex number : ' );

if true(str == 'bbb') == 1

    if (r == 1 || r == 3)
        P = sprintf('Error, Please Check your Region number');
        disp(P)

    elseif v <= R2_row && r==2
        V = v-1+R1_row;
        R = r-1;
        a=2*V+R+1;
        P = sprintf('X vector index : %d, X(index): %f ',a, X(a)');
        disp(P)
    else
        P = sprintf('Error, Please Check your input number');
        disp(P)

    end

elseif true(str == 'ccc') == 1

    if (r == 2)
        P = sprintf('Error, Please Check your Region number');
        disp(P)

    elseif v <= R1_row && r == 1
        V = v;
        R = r-1;
        a=2*V+R;
        P = sprintf('X vector index : %d, X(index): %f ',a, X(a)');
        disp(P)

    elseif v <= R3_row && (r == 3)
        V = v-1+R1_row+R2_row;
        R = r-1;
        a=2*V+R;
        P = sprintf('X vector index : %d, X(index): %f ',a, X(a)');
        disp(P)

    else
        P = sprintf('Error, Please Check your input number');
        disp(P)

    end

elseif true(str == 'aaa') == 1
    if v <= R1_row && r == 1
        V = v-1;
        R = r-1;
        a=2*V+R+1;
        P = sprintf('X vector index : %d, X(index): %f ',a, X(a)');
        disp(P)

    elseif v <= R2_row && r==2
        V = v-1+R1_row;
        R = r-1;
        a=2*V+R;
        P = sprintf('X vector index : %d, X(index): %f ',a, X(a)');
        disp(P)
    
    elseif v <= R3_row && r==3
        V = v-1+R1_row+R2_row;
        R = r-1;
        a=2*V+R-1;
        P = sprintf('X vector index : %d, X(index): %f ',a, X(a)');
        disp(P)

    else
        P = sprintf('Error, Please Check your input number');
        disp(P)

    end
end