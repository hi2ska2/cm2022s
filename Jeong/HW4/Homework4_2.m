clear; clc; close all;

% Read a vertex, element, potential data.
F_si = importdata("Element_Si.txt");
F_ox = importdata("Element_ox.txt");
F = [F_si; F_ox];
V = importdata("Vertex.txt");

% 각 Region의 edge 정렬
n=1;
for ii=1:length(F_si) % Si region
   for j=1:3
      edge_si(n,:)= [F_si(ii,2) F_si(ii,1)];
      edge_si(n+1,:)= [F_si(ii,3) F_si(ii,2)];
      edge_si(n+2,:)= [F_si(ii,1) F_si(ii,3)];
      n=n+3;
   end
end
n=1;
for ii=1:length(F_si) % ox region
   for j=1:3
      edge_ox(n,:)= [F_ox(ii,1) F_ox(ii,2)];
      edge_ox(n+1,:)= [F_ox(ii,2) F_ox(ii,3)];
      edge_ox(n+2,:)= [F_ox(ii,3) F_ox(ii,1)];
      n=n+3;
   end
end

% 중복 edge 정보 제거
edge_si=sort(edge_si, 2); % edge_si의 각 행을 오름차순으로 정렬(중복값 제거 위해)
edge_si=unique(edge_si,'rows');
edge_ox=sort(edge_ox, 2); % edge_ox의 각 행을 오름차순으로 정렬(중복값 제거 위해)
edge_ox=unique(edge_ox,'rows');

% 각 region에서 동일한 edge가 있는지 확인, Interface
edge_interface = intersect(edge_si,edge_ox,'rows');


% Visualize 
figure(1)
patch('Faces',F_si,'Vertices',V, 'EdgeColor','black','FaceColor','Blue','LineWidth',2)
xlabel('X')
ylabel('Y')
title('Region 1 (Si)')

figure(2)
patch('Faces',F_ox,'Vertices',V, 'EdgeColor','black','FaceColor','Red','LineWidth',2)
xlabel('X')
ylabel('Y')
title('Region 1 (ox)')

figure(3)
patch('Faces',F,'Vertices',V, 'EdgeColor','black','FaceColor','none','LineWidth',2)
xlabel('X')
ylabel('Y')
title('Structure')

figure(4)
patch('Faces',F_si,'Vertices',V, 'EdgeColor','black','FaceColor','Blue','LineWidth',2)
hold on
patch('Faces',F_ox,'Vertices',V, 'EdgeColor','black','FaceColor','Red','LineWidth',2)
xlabel('X')
ylabel('Y')
title('Structure')
hold on
patch('Faces',edge_interface,'Vertices',V, 'EdgeColor','Green','FaceColor','none','LineWidth',5)
hold off