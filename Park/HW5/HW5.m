clear; close all; clc;

% due to : 2022.03.22 (Thu) / HW5
% Demonstrate that your code can consider a structure with multiple (>=2) regions and contacts.


%%%%%%%% Read Files %%%%%%%%%%
fileID = fopen('Vertex_2.txt');
% buffer = fgetl(fileID);

vertex = textscan(fileID,'%f%f','Delimiter','\t');

Vertex = horzcat(vertex{1},vertex{2});
fclose(fileID);
[V_row,V_col]=size(Vertex);

fileID1 = fopen('Element_2.txt');
element = textscan(fileID1,'%f%f%f','Delimiter','\t');

Element = horzcat(element{1},element{2},element{3});
fclose(fileID1);
[E_row,col]=size(Element);

fileID2 = fopen('Contact2.txt');
contact = textscan(fileID2,'%f%f','Delimiter','\t');

Contact = horzcat(contact{1},contact{2});
fclose(fileID2);
[C_row,C_col]=size(Contact);

index = V_row;

% Geometry information

L=zeros(E_row,3);

for ii = 1:E_row
    for jj = 1:3
        if jj == 1 % between 1 and 2
            L(ii,jj) = sqrt(abs((Vertex(Element(ii,jj+1),1)-Vertex(Element(ii,jj),1))^2+(Vertex(Element(ii,jj+1),2)-Vertex(Element(ii,jj),2))^2));
        elseif jj == 2 % between 2 and 3
            L(ii,jj) = sqrt(abs((Vertex(Element(ii,jj+1),1)-Vertex(Element(ii,jj),1))^2+(Vertex(Element(ii,jj+1),2)-Vertex(Element(ii,jj),2))^2));
        elseif jj ==3 % between 3 and 1
            L(ii,jj) = sqrt(abs((Vertex(Element(ii,jj),1)-Vertex(Element(ii,jj-2),1))^2+(Vertex(Element(ii,jj),2)-Vertex(Element(ii,jj-2),2))^2));
        end
    end
end

% Vector
Area=zeros(E_row,1);
R=zeros(E_row,1);

for ii = 1:E_row
    a12 = [(Vertex(Element(ii,2),1)-Vertex(Element(ii,1),1)),(Vertex(Element(ii,2),2)-Vertex(Element(ii,1),2))];
    a13 = [(Vertex(Element(ii,3),1)-Vertex(Element(ii,1),1)),(Vertex(Element(ii,3),2)-Vertex(Element(ii,1),2))];
    a23 = [(Vertex(Element(ii,3),1)-Vertex(Element(ii,2),1)),(Vertex(Element(ii,3),2)-Vertex(Element(ii,2),2))];
    a21=-a12;
    a31=-a13;
    a32=-a23;

    Area(ii,1) = (1/2)*abs(a12(1,1)*a13(1,2)-a12(1,2)*a13(1,1));
    R(ii,1) = (L(ii,1)*L(ii,2)*L(ii,3))/(4*Area(ii,1));
end

%%%%%%% Area %%%%%%%%
edge=zeros(E_row,3);

for ii = 1:E_row
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

%%%% Potential %%%%
Jaco = zeros(index,index);
res = zeros(index,1);

for ii=1:E_row
    for jj=1:3

        if jj==1
            Jaco(Element(ii,jj), Element(ii,jj+2)) = Jaco(Element(ii,jj), Element(ii,jj+2))+ (edge(ii,jj+2)/L(ii,jj+2));
            Jaco(Element(ii,jj), Element(ii,jj+1)) = Jaco(Element(ii,jj), Element(ii,jj+1))+(edge(ii,jj)/L(ii,jj));
            Jaco(Element(ii,jj), Element(ii,jj)) = Jaco(Element(ii,jj), Element(ii,jj))-(edge(ii,jj+2)/L(ii,jj+2))-(edge(ii,jj)/L(ii,jj));

        elseif jj ==2
            Jaco(Element(ii,jj), Element(ii,jj+1)) = Jaco(Element(ii,jj), Element(ii,jj+1))+(edge(ii,jj)/L(ii,jj));
            Jaco(Element(ii,jj), Element(ii,jj)) =  Jaco(Element(ii,jj), Element(ii,jj))-(edge(ii,jj)/L(ii,jj))-(edge(ii,jj-1)/L(ii,jj-1));
            Jaco(Element(ii,jj), Element(ii,jj-1)) = Jaco(Element(ii,jj), Element(ii,jj-1))+(edge(ii,jj-1)/L(ii,jj-1));

        elseif jj == 3
            Jaco(Element(ii,jj), Element(ii,jj-2)) = Jaco(Element(ii,jj), Element(ii,jj-2))+(edge(ii,jj)/L(ii,jj));
            Jaco(Element(ii,jj), Element(ii,jj-1)) = Jaco(Element(ii,jj), Element(ii,jj-1))+(edge(ii,jj-1)/L(ii,jj-1));
            Jaco(Element(ii,jj), Element(ii,jj)) = Jaco(Element(ii,jj), Element(ii,jj))-(edge(ii,jj)/L(ii,jj))-(edge(ii,jj-1)/L(ii,jj-1));
        end
    end
end

%%%%% Dirichlet Boundary Condition %%%%%
for n=1:C_col
    for ii =1: C_row
        if ii <= 3
            Jaco(Contact(ii,n),:) =0;
            Jaco(Contact(ii,n),Contact(ii,n)) = 1;
            res(Contact(ii,n),1)=1;
        
         elseif ii > 3
            Jaco(Contact(ii,n),:) =0;
            Jaco(Contact(ii,n),Contact(ii,n)) = 1;
            res(Contact(ii,n),1)=0;
% 
%         elseif ii > 6 && ii <= 19
%             Jaco(Contact(ii,n),:) =0;
%             Jaco(Contact(ii,n),Contact(ii,n)) = 1;
%             res(Contact(ii,n),1)=1;
% 
%         elseif ii > 19
%             Jaco(Contact(ii,n),:) =0;
%             Jaco(Contact(ii,n),Contact(ii,n)) = 1;
%             res(Contact(ii,n),1)=1;

        
        end
    end
end

%%%% Result %%%%
phi = Jaco\res;

set(gcf,'Color','w')
patch('Faces',Element,'Vertices',Vertex,'FaceVertexCData',phi,'FaceColor','interp','Marker','o','MarkerFaceColor','Black','Markersize',8);
colorbar

%%%%%%%% Specify Region 1,2 %%%%%%%%%%

fileID2 = fopen('Element_2_region1.txt');
element_Region1 = textscan(fileID2,'%f%f%f','Delimiter','\t');

Element_Region1 = horzcat(element_Region1{1},element_Region1{2},element_Region1{3});
fclose(fileID2);

fileID3 = fopen('Element_2_region2.txt');
element_Region2 = textscan(fileID3,'%f%f%f','Delimiter','\t');

Element_Region2 = horzcat(element_Region2{1},element_Region2{2},element_Region2{3});
fclose(fileID3);

fileID4 = fopen('Element_2_region3.txt');
element_Region3 = textscan(fileID4,'%f%f%f','Delimiter','\t');

Element_Region3 = horzcat(element_Region3{1},element_Region3{2},element_Region3{3});
fclose(fileID4);


% fileID5 = fopen('Element_region_substrate.txt');
% element_Region4 = textscan(fileID5,'%f%f%f','Delimiter','\t');
% 
% Element_Region4 = horzcat(element_Region4{1},element_Region4{2},element_Region4{3});
% fclose(fileID5);

[R_row,R_Col] = size(Element_Region3);

%%%%%%%%% Find the interface Edge %%%%%%%%%%%

S1=intersect(Element_Region1, Element_Region3);
j =1;

for ii = 1:R_row
    if size(intersect(Element_Region3(ii,:),S1'),2) == 2
        Edge1(j,:) = intersect(Element_Region3(ii,:),S1');
        j = j+1;
    end
end

S2=intersect(Element_Region3, Element_Region2);
j =1;

for ii = 1:R_row
    if size(intersect(Element_Region3(ii,:),S2'),2) == 2
        Edge2(j,:) = intersect(Element_Region3(ii,:),S2);
        j = j+1;
    end
end

S3=intersect(Element_Region1, Element_Region2);
j =1;

for ii = 1:11
    if size(intersect(Element_Region1(ii,:),S3'),2) == 2
        Edge3(j,:) = intersect(Element_Region1(ii,:),S3');
        j = j+1;
    end
end

%%%%%%% Visualize %%%%%%%%%

set(gcf,'Color','w')
patch('Faces',Element_Region1,'Vertices',Vertex,'EdgeColor','Red','FaceColor','none','LineWidth',3);
patch('Faces',Element_Region2,'Vertices',Vertex,'EdgeColor','yellow','FaceColor','none','LineWidth',3);
patch('Faces',Element_Region3,'Vertices',Vertex,'EdgeColor','Blue','FaceColor','none','LineWidth',3);
% % patch('Faces',Element_Region4,'Vertices',Vertex,'EdgeColor','Black','FaceColor','none','LineWidth',3);
patch('Faces',Edge1,'Vertices',Vertex,'EdgeColor','magenta','FaceColor','none','LineWidth',3);
patch('Faces',Edge2,'Vertices',Vertex,'EdgeColor','Green','FaceColor','none','LineWidth',3);
patch('Faces',Edge3,'Vertices',Vertex,'EdgeColor','White','FaceColor','none','LineWidth',3);
patch('Faces',Contact,'Vertices',Vertex,'EdgeColor','Cyan','FaceColor','none','LineWidth',3);
