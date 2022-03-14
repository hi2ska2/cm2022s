clear; clc; close all;

% Read a vertex, element, potential data.
F = importdata("Element.txt");
V = importdata("Vertex.txt");
P = importdata("potential.txt");

% Visualize 
figure(1);
patch('Faces',F,'Vertices',V, 'FaceVertexCData',P, 'EdgeColor','black','FaceColor','none','LineWidth',2);
xlabel('X');
ylabel('Y');
title('Structure')

figure(2);
patch('Faces',F,'Vertices',V, 'FaceVertexCData',P, 'EdgeColor','black','FaceColor','interp');
title('Structure (N=16)')
xlabel('X');
ylabel('Y');
colorbar