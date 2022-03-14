clear; clc; close all;

% Read a vertex, element, potential data.
F = importdata("Element.txt");
V = importdata("Vertex.txt");
P = importdata("potential.txt");

% Visualize 
figure(1);
patch('Faces', F , 'Vertices', V , 'EdgeColor','black','FaceColor','none','LineWidth',2);
xlabel('X');https://github.com/hi2ska2/cm2022s/blob/main/Jeong/HW3/Homework3_visualize.m
ylabel('Y');
title('Structure')

figure(2);
patch('Faces', F ,'Vertices', V , 'FaceVertexCData', P , 'EdgeColor','black','FaceColor','interp');
title('Structure (N=16)')
xlabel('X');
ylabel('Y');
colorbar
