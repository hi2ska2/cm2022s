% due to : 2022.03.15 (Thu) / HW3-1

% Solve the Laplace equation for a circle.

%%%%%%%% Genearte Vertex file %%%%%%%%%%
index = 17;

Vertex = zeros(index,2);
R1 = 6;
R2 = 4;
R3 = 2;

for ii = 1:17
    for R = R3:2:R1

        if ii <= 8 && R == R1
            Vertex(ii,1) = R*cos(((ii-1)*pi)/4);
            Vertex(ii,2) = R*sin(((ii-1)*pi)/4);    
        
        elseif (ii > 8 && ii <= 12 ) && R == R2
            Vertex(ii,1) = R*cos((ii-1)*pi/2);
            Vertex(ii,2) = R*sin(((ii-1)*pi)/2);

        elseif (ii > 12 && ii <= 16) && R == R3
            Vertex(ii,1) = R*cos((ii-1)*pi/2-(pi/4));
            Vertex(ii,2) = R*sin(((ii-1)*pi)/2-(pi/4));

        end   
      
    end
end

result = [Vertex(:,1)'; Vertex(:,2)'];

fileID = fopen('Vertex.txt','W');
fprintf(fileID,'%6s %12s\n','x','y');
fprintf(fileID,'%6f %12f\n', result);
fclose(fileID);



fileID1 = fopen('element.txt');
element = textscan(fileID1,'%f%f%f','Delimiter','\t');

Element = horzcat(element{1},element{2},element{3});
fclose(fileID1);

% Geometry information

L=zeros(24,3);

for ii = 1:24
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
Area=zeros(24,1);
R=zeros(24,1);

for ii = 1:24
    a12 = [(Vertex(Element(ii,2),1)-Vertex(Element(ii,1),1)),(Vertex(Element(ii,2),2)-Vertex(Element(ii,1),2))];
    a13 = [(Vertex(Element(ii,3),1)-Vertex(Element(ii,1),1)),(Vertex(Element(ii,3),2)-Vertex(Element(ii,1),2))];
    a23 = [(Vertex(Element(ii,3),1)-Vertex(Element(ii,2),1)),(Vertex(Element(ii,3),2)-Vertex(Element(ii,2),2))];
    a21=-a12;
    a31=-a13;
    a32=-a23;

    Area(ii,1) = (1/2)*abs(a12(1,1)*a13(1,2)-a12(1,2)*a13(1,1));
    R(ii,1) = (L(ii,1)*L(ii,2)*L(ii,3))/(4*Area(ii,1));
end

edge=zeros(24,3);

for ii = 1:24
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

Jaco = zeros(index,index);
res = zeros(index,1);

for ii=1:24
    for jj=1:3
        if Element(ii,jj) == 3
           Jaco(Element(ii,jj),Element(ii,jj)) = 1;
           res(Element(ii,jj),1)=1;

        elseif Element(ii,jj) == 7
            Jaco(Element(ii,jj),Element(ii,jj)) = 1;
            
        elseif jj==1
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

phi = Jaco\res;

% patch('Faces',Element,'Vertices',Vertex,'FaceVertexCData',phi,'FaceColor','interp');
% colorbar