% due to : 2022.03.17 (Thur) / HW4-1

% Write *.element file to specifiy the yin and yang pattern.

%%%%%%%% Read Vertex information %%%%%%%%%%
clear; close all; clc;

fileID = fopen('Vertex4.txt');
buffer = fgetl(fileID);

vertex = textscan(fileID,'%f%f','Delimiter','\t');

Vertex = horzcat(vertex{1},vertex{2});
fclose(fileID);

%%%%%%%% Read Element Information %%%%%%%%%%

fileID1 = fopen('element4.txt');
element = textscan(fileID1,'%f%f%f','Delimiter','\t');

Element = horzcat(element{1},element{2},element{3});
fclose(fileID1);

%%%%%%%% Specify Region 1,2 %%%%%%%%%%

fileID2 = fopen('element_Region_1.txt');
element_Region1 = textscan(fileID2,'%f%f%f','Delimiter','\t');

Element_Region1 = horzcat(element_Region1{1},element_Region1{2},element_Region1{3});
fclose(fileID2);

fileID3 = fopen('element_Region_2.txt');
element_Region2 = textscan(fileID3,'%f%f%f','Delimiter','\t');

Element_Region2 = horzcat(element_Region2{1},element_Region2{2},element_Region2{3});
fclose(fileID3);

%%%%%%%%% Find the interface Edge %%%%%%%%%%%

S=intersect(Element_Region1, Element_Region2);
j =1;

for ii = 1:24
    if size(intersect(Element_Region2(ii,:),S'),2) == 2
        Edge(j,:) = intersect(Element_Region2(ii,:),S');
        j = j+1;
    end
end

%%%%%%% Visualize %%%%%%%%%

set(gcf,'Color','w')
patch('Faces',Element_Region1,'Vertices',Vertex,'EdgeColor','Red','FaceColor','none','LineWidth',3);
patch('Faces',Element_Region2,'Vertices',Vertex,'EdgeColor','Blue','FaceColor','none','LineWidth',3);
patch('Faces',Edge,'Vertices',Vertex,'EdgeColor','Green','FaceColor','none','LineWidth',3);
