% due to : 2022.03.15 (Thu) / HW3-2

% Visulize the read structure

%%%%%%%% Read Vertex information %%%%%%%%%%

fileID = fopen('Rec_Vertex.txt');
buffer = fgetl(fileID);

vertex = textscan(fileID,'%f%f','Delimiter','\t');

Vertex = horzcat(vertex{1},vertex{2});
fclose(fileID);

%%%%%%%% Read Element Information %%%%%%%%%%

fileID1 = fopen('Rec_element.txt');
element = textscan(fileID1,'%f%f%f','Delimiter','\t');

Element = horzcat(element{1},element{2},element{3});
fclose(fileID1);

patch('Faces',Element,'Vertices',Vertex,'EdgeColor','Black','FaceColor','none','LineWidth',3);


