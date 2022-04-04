clear; close all; clc;

% due to : 2022.04.05 (Tue) / HW9
% Repeat the Source-free poisson equation solver (advanced way)

clear; close all; clc;

% due to : 2022.03.24 (Thur) / HW6
% Solve the Source-free poisson equation for several structures.

%%%%%%%% Read Files %%%%%%%%%%
fileID = fopen('Vertex2.txt');

vertex = textscan(fileID,'%f%f','Delimiter','\t');

Vertex = horzcat(vertex{1},vertex{2});
fclose(fileID);
[V_row,V_col]=size(Vertex);

fileID1 = fopen('Element2.txt');
element = textscan(fileID1,'%f%f%f','Delimiter','\t');

Element = horzcat(element{1},element{2},element{3});
fclose(fileID1);
[E_row,col]=size(Element);

fileID2 = fopen('contact2.txt');
contact = textscan(fileID2,'%f%f','Delimiter','\t');

Contact = horzcat(contact{1},contact{2});
fclose(fileID2);
[C_row,C_col]=size(Contact);

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

fileID4 = fopen('Element_region_substrate.txt');
element_Region3 = textscan(fileID4,'%f%f%f','Delimiter','\t');

Element_Region3 = horzcat(element_Region3{1},element_Region3{2},element_Region3{3});
fclose(fileID4);

% fileID5 = fopen('Element_region_substrate.txt');
% element_Region4 = textscan(fileID5,'%f%f%f','Delimiter','\t');
% 
% Element_Region4 = horzcat(element_Region4{1},element_Region4{2},element_Region4{3});
% fclose(fileID5);

R_row1 = size(Element_Region1,1);
R_row2 = size(Element_Region2,1);
R_row3 = size(Element_Region3,1);
esi = 11.7; eox = 3.9;

T=R_row1+R_row2+R_row3;


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

%%%%%%%%% Find the interface Edge %%%%%%%%%%%

S1=intersect(Element_Region1, Element_Region2);
j =1;

for ii = 1:R_row2
    if size(intersect(Element_Region2(ii,:),S1'),2) == 2
        Edge1(j,:) = intersect(Element_Region2(ii,:),S1');
        j = j+1;
    end
end

S2=intersect(Element_Region1, Element_Region3);
j =1;

for ii = 1:R_row2
    if size(intersect(Element_Region1(ii,:),S2'),2) == 2
        Edge2(j,:) = intersect(Element_Region1(ii,:),S2);
        j = j+1;
    end
end

S3=intersect(Element_Region2, Element_Region3);
j =1;

for ii = 1:R_row3
    if size(intersect(Element_Region3(ii,:),S3'),2) == 2
        Edge3(j,:) = intersect(Element_Region3(ii,:),S3');
        j = j+1;
    end
end

% %%%% Vertex index for each region %%%

A = size(unique(Element_Region1),1);
B = size(unique(Element_Region2),1);
C = size(unique(Element_Region3),1);

%%%% Potential %%%%
Jaco = zeros(index+12,index+12);
res = zeros(index+12,1);

Jaco1 = zeros(index+12,index+12);
Jaco2 = zeros(index+12,index+12);
Jaco3 = zeros(index+12,index+12);

for ii=1:R_row1
    for jj=1:3

        if jj==1
            Jaco1(Element_Region1(ii,jj), Element_Region1(ii,jj+2)) = Jaco1(Element_Region1(ii,jj), Element_Region1(ii,jj+2))+ eox*(edge(ii,jj+2)/L(ii,jj+2));
            Jaco1(Element_Region1(ii,jj), Element_Region1(ii,jj+1)) = Jaco1(Element_Region1(ii,jj), Element_Region1(ii,jj+1))+eox*(edge(ii,jj)/L(ii,jj));
            Jaco1(Element_Region1(ii,jj), Element_Region1(ii,jj)) = Jaco1(Element_Region1(ii,jj), Element_Region1(ii,jj))-eox*(edge(ii,jj+2)/L(ii,jj+2))-eox*(edge(ii,jj)/L(ii,jj));

        elseif jj ==2
            Jaco1(Element_Region1(ii,jj), Element_Region1(ii,jj+1)) = Jaco1(Element_Region1(ii,jj), Element_Region1(ii,jj+1))+eox*(edge(ii,jj)/L(ii,jj));
            Jaco1(Element_Region1(ii,jj), Element_Region1(ii,jj)) =  Jaco1(Element_Region1(ii,jj), Element_Region1(ii,jj))-eox*(edge(ii,jj)/L(ii,jj))-eox*(edge(ii,jj-1)/L(ii,jj-1));
            Jaco1(Element_Region1(ii,jj), Element_Region1(ii,jj-1)) = Jaco1(Element_Region1(ii,jj), Element_Region1(ii,jj-1))+eox*(edge(ii,jj-1)/L(ii,jj-1));

        elseif jj == 3
            Jaco1(Element_Region1(ii,jj), Element_Region1(ii,jj-2)) = Jaco1(Element_Region1(ii,jj), Element_Region1(ii,jj-2))+eox*(edge(ii,jj)/L(ii,jj));
            Jaco1(Element_Region1(ii,jj), Element_Region1(ii,jj-1)) = Jaco1(Element_Region1(ii,jj), Element_Region1(ii,jj-1))+eox*(edge(ii,jj-1)/L(ii,jj-1));
            Jaco1(Element_Region1(ii,jj), Element_Region1(ii,jj)) = Jaco1(Element_Region1(ii,jj), Element_Region1(ii,jj))-eox*(edge(ii,jj)/L(ii,jj))-eox*(edge(ii,jj-1)/L(ii,jj-1));
        
        end
    end
end

for ii=1:R_row2
    for jj=1:3
        if jj==1
            Jaco2(Element_Region2(ii,jj), Element_Region2(ii,jj+2)) = Jaco2(Element_Region2(ii,jj), Element_Region2(ii,jj+2))+ eox*(edge(ii+R_row1,jj+2)/L(ii+R_row1,jj+2));
            Jaco2(Element_Region2(ii,jj), Element_Region2(ii,jj+1)) = Jaco2(Element_Region2(ii,jj), Element_Region2(ii,jj+1))+eox*(edge(ii+R_row1,jj)/L(ii+R_row1,jj));
            Jaco2(Element_Region2(ii,jj), Element_Region2(ii,jj)) = Jaco2(Element_Region2(ii,jj), Element_Region2(ii,jj))-eox*(edge(ii+R_row1,jj+2)/L(ii+R_row1,jj+2))-eox*(edge(ii+R_row1,jj)/L(ii+R_row1,jj));

        elseif jj ==2

            Jaco2(Element_Region2(ii,jj), Element_Region2(ii,jj+1)) = Jaco2(Element_Region2(ii,jj), Element_Region2(ii,jj+1))+esi*(edge(ii+R_row1,jj)/L(ii+R_row1,jj));
            Jaco2(Element_Region2(ii,jj), Element_Region2(ii,jj)) =  Jaco2(Element_Region2(ii,jj), Element_Region2(ii,jj))-esi*(edge(ii+R_row1,jj)/L(ii+R_row1,jj))-esi*(edge(ii+R_row1,jj-1)/L(ii+R_row1,jj-1));
            Jaco2(Element_Region2(ii,jj), Element_Region2(ii,jj-1)) = Jaco2(Element_Region2(ii,jj), Element_Region2(ii,jj-1))+esi*(edge(ii+R_row1,jj-1)/L(ii+R_row1,jj-1));

        elseif jj == 3

            Jaco2(Element_Region2(ii,jj), Element_Region2(ii,jj-2)) = Jaco2(Element_Region2(ii,jj), Element_Region2(ii,jj-2))+esi*(edge(ii+R_row1,jj)/L(ii+R_row1,jj));
            Jaco2(Element_Region2(ii,jj), Element_Region2(ii,jj-1)) = Jaco2(Element_Region2(ii,jj), Element_Region2(ii,jj-1))+esi*(edge(ii+R_row1,jj-1)/L(ii+R_row1,jj-1));
            Jaco2(Element_Region2(ii,jj), Element_Region2(ii,jj)) = Jaco2(Element_Region2(ii,jj), Element_Region2(ii,jj))-esi*(edge(ii+R_row1,jj)/L(ii+R_row1,jj))-esi*(edge(ii+R_row1,jj-1)/L(ii+R_row1,jj-1));

        end
    end
end

for ii=1:R_row3
    for jj=1:3

        if jj==1
            Jaco3(Element_Region3(ii,jj), Element_Region3(ii,jj+2)) = Jaco3(Element_Region3(ii,jj), Element_Region3(ii,jj+2))+ esi*(edge(ii+R_row2+R_row1,jj+2)/L(ii+R_row2+R_row1,jj+2));
            Jaco3(Element_Region3(ii,jj), Element_Region3(ii,jj+1)) = Jaco3(Element_Region3(ii,jj), Element_Region3(ii,jj+1))+esi*(edge(ii+R_row2+R_row1,jj)/L(ii+R_row2+R_row1,jj));
            Jaco3(Element_Region3(ii,jj), Element_Region3(ii,jj)) = Jaco3(Element_Region3(ii,jj), Element_Region3(ii,jj))-esi*(edge(ii+R_row2+R_row1,jj+2)/L(ii+R_row2+R_row1,jj+2))-esi*(edge(ii+R_row2+R_row1,jj)/L(ii+R_row2+R_row1,jj));

        elseif jj ==2
            Jaco3(Element_Region3(ii,jj), Element_Region3(ii,jj+1)) = Jaco3(Element_Region3(ii,jj), Element_Region3(ii,jj+1))+esi*(edge(ii+R_row2+R_row1,jj)/L(ii+R_row2+R_row1,jj));
            Jaco3(Element_Region3(ii,jj), Element_Region3(ii,jj)) =  Jaco3(Element_Region3(ii,jj), Element_Region3(ii,jj))-esi*(edge(ii+R_row2+R_row1,jj)/L(ii+R_row2+R_row1,jj))-esi*(edge(ii+R_row2+R_row1,jj-1)/L(ii+R_row2+R_row1,jj-1));
            Jaco3(Element_Region3(ii,jj), Element_Region3(ii,jj-1)) = Jaco3(Element_Region3(ii,jj), Element_Region3(ii,jj-1))+esi*(edge(ii+R_row2+R_row1,jj-1)/L(ii+R_row2+R_row1,jj-1));

        elseif jj == 3
            Jaco3(Element_Region3(ii,jj), Element_Region3(ii,jj-2)) = Jaco3(Element_Region3(ii,jj), Element_Region3(ii,jj-2))+esi*(edge(ii+R_row2+R_row1,jj)/L(ii+R_row2+R_row1,jj));
            Jaco3(Element_Region3(ii,jj), Element_Region3(ii,jj-1)) = Jaco3(Element_Region3(ii,jj), Element_Region3(ii,jj-1))+esi*(edge(ii+R_row2+R_row1,jj-1)/L(ii+R_row2+R_row1,jj-1));
            Jaco3(Element_Region3(ii,jj), Element_Region3(ii,jj)) = Jaco3(Element_Region3(ii,jj), Element_Region3(ii,jj))-esi*(edge(ii+R_row2+R_row1,jj)/L(ii+R_row2+R_row1,jj))-esi*(edge(ii+R_row2+R_row1,jj-1)/L(ii+R_row2+R_row1,jj-1));
        end
    end
end

%%%%% Dirichlet Boundary Condition %%%%%
for n=1:C_col
    for ii =1: C_row
        if ii <= 2
            Jaco(Contact(ii,n),:) = 0;
            Jaco(Contact(ii,n),Contact(ii,n)) = 1;
            res(Contact(ii,n),1) = 0;
        
         elseif ii > 2 && ii <= 4
            Jaco(Contact(ii,n),:) = 0;
            Jaco(Contact(ii,n),Contact(ii,n)) = 1;
            res(Contact(ii,n),1) = 0;
 
        elseif ii > 4 && ii <= 6
            Jaco(Contact(ii,n),:) =0;
            Jaco(Contact(ii,n),Contact(ii,n)) = 1;
            res(Contact(ii,n),1)=0;

        elseif ii > 6
            Jaco3(Contact(ii,n),:) = 0;
            Jaco3(Contact(ii,n),Contact(ii,n)) = 1;
            res(Contact(ii,n),1) = 1;

        
        end
    end
end

J_row=size(Jaco,1);

%%% Reindexing the Jacobian %%%%

for ii = 1:J_row


    if ii <= 45
        Jaco(ii,:)=Jaco(ii,:)+Jaco3(ii,:);

    elseif ii == 46
        Jaco(ii+4,:)=Jaco(ii+4,:)+circshift(Jaco3(ii,:),4,2);

    elseif ii > 46 && ii <= 50
        Jaco(ii+4,:)=Jaco(ii+4,:)+circshift(Jaco3(ii,:),4,2);

    elseif ii > 50 && ii <= 52
        Jaco(ii+8,:)=Jaco(ii+8,:)+circshift(Jaco3(ii,:),8,2);

    elseif ii >= 56 && ii <= 58 % 65, 67 
        Jaco(ii+9,:)=Jaco(ii+9,:)+circshift(Jaco3(ii,:),9,2);

    elseif ii >= 62 && ii <= 63
       Jaco(ii+10,:)=Jaco(ii+10,:)+circshift(Jaco1(ii,:),10,2); 

    elseif ii >= 66 && ii <= 67
        Jaco(ii+10,:)=Jaco(ii+10,:)+circshift(Jaco1(ii,:),10,2); 

    elseif ii >= 69 && ii <= 70
        Jaco(ii+12,:)=Jaco(ii+12,:)+circshift(Jaco2(ii,:),12,2);  

    end



    if ii >= 42 && ii <= 45
        Jaco(ii+4,:)=Jaco(ii+4,:)+circshift(Jaco1(ii,:),4,2);

    elseif ii > 46 && ii <= 50
        Jaco(ii+8,:)=Jaco(ii+8,:)+circshift(Jaco2(ii,:),8,2);

    elseif ii > 52 && ii <= 56
        Jaco(ii+8,:)=Jaco(ii+8,:)+circshift(Jaco1(ii,:),8,2);

    elseif ii >= 58 && ii <= 61
        Jaco(ii+10,:)=Jaco(ii+10,:)+circshift(Jaco2(ii,:),10,2);

    elseif ii >= 64 && ii <= 65
        Jaco(ii+10,:)=Jaco(ii+10,:)+circshift(Jaco2(ii,:),10,2);
    
    elseif ii >= 67 && ii <= 69
        Jaco(ii+11,:)=Jaco(ii+11,:)+circshift(Jaco3(ii,:),11,2);
    
  
    end
    for jj =1:6

        if ii < 46 && ii==S2(jj,1)
            Jaco(ii,:)=Jaco1(ii,:)+circshift(Jaco3(ii,:),4,2);
            Jaco(ii+4,:)=0;
            Jaco(ii+4,ii)= -1;
            Jaco(ii+4,ii+4)=1;

        elseif ii < 51 && ii == S3(jj,1) 
            Jaco(ii+4,:)=circshift(Jaco2(ii,:),4,2)+circshift(Jaco3(ii,:),8,2);
            Jaco(ii+8,:)=0;
            Jaco(ii+8,ii+4)= -1;
            Jaco(ii+8,ii+8)=1;

        elseif ii == 56
            Jaco(ii+8,:)=circshift(Jaco1(ii,:),8,2)+circshift(Jaco3(ii,:),9,2);
            Jaco(ii+9,:)=0;
            Jaco(ii+9,ii+8)= -1;
            Jaco(ii+9,ii+9)=1;

        elseif ii == 58
            Jaco(ii+9,:)=circshift(Jaco2(ii,:),9,2)+circshift(Jaco3(ii,:),10,2);
            Jaco(ii+10,:)=0;
            Jaco(ii+10,ii+9)= -1;
            Jaco(ii+10,ii+10)=1;  

        elseif ii == 67
            Jaco(ii+10,:)=circshift(Jaco1(ii,:),10,2)+circshift(Jaco3(ii,:),11,2);
            Jaco(ii+11,:)=0;
            Jaco(ii+11,ii+10)= -1;
            Jaco(ii+11,ii+11)=1;

        elseif ii == 69
            Jaco(ii+11,:)=circshift(Jaco2(ii,:),11,2)+circshift(Jaco3(ii,:),12,2);
            Jaco(ii+12,:)=0;
            Jaco(ii+12,ii+11)= -1;
            Jaco(ii+12,ii+12)=1;    
        end
    end
    
end

%%%% Result %%%%
phi = Jaco\res;
