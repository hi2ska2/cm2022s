% due to : 2022.03.17 (Thur) / HW4-1

% Solve the Laplace equation for a circle (using Contact file).

%%%%%%%% Genearte Vertex file %%%%%%%%%%
clear; close all; clc;

index = 33;

Vertex = zeros(index,2);
R1 = 6;
R2 = 4;
R3 = 2;

for ii = 1:index
    for R = R3:2:R1

        if ii <= 16 && R == R1
            Vertex(ii,1) = R*cos(((ii-1)*pi)/8);
            Vertex(ii,2) = R*sin(((ii-1)*pi)/8);    
        
        elseif (ii > 16 && ii <= 24 ) && R == R2
            Vertex(ii,1) = R*cos((ii-1)*pi/4);
            Vertex(ii,2) = R*sin(((ii-1)*pi)/4);

        elseif (ii > 24 && ii <= 32) && R == R3
            Vertex(ii,1) = R*cos((ii-1)*pi/4-(pi/8));
            Vertex(ii,2) = R*sin(((ii-1)*pi)/4-(pi/8));

        end   
      
    end
end

result = [Vertex(:,1)'; Vertex(:,2)'];

fileID = fopen('Vertex4.txt','W');
fprintf(fileID,'%6s %12s\n','x','y');
fprintf(fileID,'%6f %12f\n', result);
fclose(fileID);

%%%%%%%% Read Files %%%%%%%%%%

fileID1 = fopen('element4.txt');
element = textscan(fileID1,'%f%f%f','Delimiter','\t');

Element = horzcat(element{1},element{2},element{3});
fclose(fileID1);
[E_row,col]=size(Element);

fileID2 = fopen('contact.txt');
contact = textscan(fileID2,'%f%f%f','Delimiter','\t');

Contact = horzcat(contact{1},contact{2},contact{3});
fclose(fileID2);
[row,C_col]=size(Contact);

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

%% Area 
edge=zeros(E_row,3);

for ii = 1:E_row
    for jj = 1:3
        if jj == 1 % between 1 and 2
            edge(ii,jj) = sqrt(R(ii,1)^2-(L(ii,jj)/2)^2);
        elseif jj == 2 % between 2 and 3
            edge(ii,jj) = sqrt(R(ii,1)^2-(L(ii,jj)/2)^2);
        elseif jj ==3 % between 3 and 1
            edge(ii,jj) = sqrt(R(ii,1)^2-(L(ii,jj)/2)^2);
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

%%%%% Dirichlet Boundary Condition
for n=1:C_col

    Jaco(Contact(1,n),:) =0;
    Jaco(Contact(1,n),Contact(1,n)) = 1;
    res(Contact(1,n),1)=1;

    Jaco(Contact(2,n),:) =0;
    Jaco(Contact(2,n),Contact(2,n)) = 1;
    res(Contact(2,n),1)=0;

    Jaco(Contact(3,n),:) =0;
    Jaco(Contact(3,n),Contact(3,n)) = 1;
    res(Contact(3,n),1)=-1;
end

%%%% Result %%%%
phi = Jaco\res;
set(gcf,'Color','w')
patch('Faces',Element,'Vertices',Vertex,'FaceVertexCData',phi,'FaceColor','interp');
colorbar