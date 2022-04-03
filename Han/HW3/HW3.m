%%% HW3
%%% written by Seong-Min, Han
%%% Due: AM 08:00, March 15, 2022

N = 4; % 가중치 % 각 원의 점의 개수는 N * numc
numc = 4; % circle 개수 
numf = (2*numc) * (2*numc); % face 개수
numv = 2*numc*numc + 2*numc + 1; % Vertex 개수
V = importdata("vertex.txt");
F = importdata("element.txt");
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% vertex와 vertex 사이의 거리
len = zeros(numf,3); 
for i=1:numf
    for j =1:3
        if j==3
            len(i,j) = norm(V(F(i,j),:)-V(F(i,1),:));
        else
            len(i,j) = norm(V(F(i,j),:)-V(F(i,j+1),:));
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% vertex 1->2, 1->3을 연결하는 vector와 vec의 norm 값
v21 = zeros(numf,2); 
v31 = zeros(numf,2);
norm_v21 = zeros(numf,1); 
norm_v31 = zeros(numf,1); 
for i=1:numf
    v21(i,:) = V(F(i,2),:)-V(F(i,1),:);
    v31(i,:) = V(F(i,3),:)-V(F(i,1),:);
    norm_v21(i,1) = norm(v21(i,:));
    norm_v31(i,1) = norm(v31(i,:));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%element로 이루어진 face의 area와 face의 외접원 radius
area = zeros(numf,1); 
rad = zeros(numf,1);
for i=1:numf
    area(i,1) = sqrt(norm_v21(i,1)*norm_v21(i,1)*norm_v31(i,1)*norm_v31(i,1) - dot(v21(i,:),v31(i,:))*dot(v21(i,:),v31(i,:)))/2;
    rad(i,1) = (len(i,1)*len(i,2)*len(i,3)) / (4*area(i,1));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%face의 edge / edge1(1-2), edge2(2-3), edge3(3-1)
edge = zeros(numf,1);
for i=1:numf
    for j=1:3
        edge(i,j) = sqrt(rad(i,1)*rad(i,1) - (len(i,j)/2)*(len(i,j)/2));
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% A, b, phi matrix 
A = zeros(numv,numv);
b = zeros(numv,1);
for i=1:numf
    for j=1:3
        if j==1
            A(F(i,j),F(i,j)) = A(F(i,j),F(i,j)) - edge(i,1)/len(i,1) - edge(i,3)/len(i,3);
            A(F(i,j),F(i,2)) = A(F(i,j),F(i,2)) + edge(i,1)/len(i,1);
            A(F(i,j),F(i,3)) = A(F(i,j),F(i,3)) + edge(i,3)/len(i,3);
        elseif j==2
            A(F(i,j),F(i,j)) = A(F(i,j),F(i,j)) - edge(i,2)/len(i,2) - edge(i,1)/len(i,1);
            A(F(i,j),F(i,3)) = A(F(i,j),F(i,3)) + edge(i,2)/len(i,2);
            A(F(i,j),F(i,1)) = A(F(i,j),F(i,1)) + edge(i,1)/len(i,1);
        elseif j==3
            A(F(i,j),F(i,j)) = A(F(i,j),F(i,j)) - edge(i,3)/len(i,3) - edge(i,2)/len(i,2);
            A(F(i,j),F(i,1)) = A(F(i,j),F(i,1)) + edge(i,3)/len(i,3);
            A(F(i,j),F(i,2)) = A(F(i,j),F(i,2)) + edge(i,2)/len(i,2);        
        end
    end
end
A(2*numc*numc-numc+2,:)= 0; %potential이 1인 부분 - 가장 위 vertex
A(2*numc*numc+numc+2,:)= 0; %potential이 0인 부분 - 가장 아래 vertex
A(2*numc*numc-numc+2,2*numc*numc-numc+2) = 1;
A(2*numc*numc+numc+2,2*numc*numc+numc+2) = 1;
b(2*numc*numc-numc+2,1) = 1;
phi = A\ b;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% visualizing 
figure
patch('Faces',F,'Vertices',V, 'FaceVertexCData',phi, 'EdgeColor','black','FaceColor','interp','LineWidth',2, 'Marker','o');
title('structure visualizing')
colorbar
% % % plot(phi, 'b*')
% % % xlabel('vertex')
% % % ylabel('potential')
