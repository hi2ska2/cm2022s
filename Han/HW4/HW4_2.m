%%% HW4_2
%%% written by Seong-Min, Han

numf = 4*numc*numc+8*numc-4; % face 개수
numv = 2*numc*numc+6*numc+1; % Vertex 개수
V = importdata("vertex.txt");
F = importdata("element.txt");
contact = importdata("contact.txt");
Fox = importdata("element_ox.txt");
Fsi = importdata("element_si.txt");
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
        if j==3
            A(F(i,j),F(i,j)) = A(F(i,j),F(i,j)) - edge(i,j)/len(i,j);
            A(F(i,j),F(i,1)) = A(F(i,j),F(i,1)) + edge(i,j)/len(i,j);
        else
            A(F(i,j),F(i,j)) = A(F(i,j),F(i,j)) - edge(i,j)/len(i,j);
            A(F(i,j),F(i,j+1)) = A(F(i,j),F(i,j+1)) + edge(i,j)/len(i,j);
        end
    end
end

for ii=1:6
    if ii<=3 %potential이 1인 영역
        A(contact(ii,1),:)= 0; 
        A(contact(ii,1),contact(ii,1)) = 1;
        b(contact(1:3,1),1) = contact(1:3,2);
    elseif ii>3 && ii<=6  %potential이 0인 영역
        A(contact(ii,1),:)= 0; 
        A(contact(ii,1),contact(ii,1)) = 1;
        b(contact(4:6,1),1) = contact(4:6,2);
    end
end
phi = A\ b;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
inter_edge = zeros(numf/2,2);
for i=1:numf/2
    for ii=1:numf/2
        for j=1:3
            if Fox(i,j) == Fsi(ii,j)
                if j == 1
                    if Fox(i,2) == Fsi(ii,3)
                        inter_edge(i,1) = Fox(i,j);
                        inter_edge(i,2) = Fox(i,2);
                    elseif Fox(i,3) == Fsi(ii,2)
                        inter_edge(i,1) = Fox(i,j);
                        inter_edge(i,2) = Fox(i,3);
                    end

                elseif j == 2
                    if Fox(i,3) == Fsi(ii,1)
                        inter_edge(i,1) = Fox(i,j);
                        inter_edge(i,2) = Fox(i,3);
                    elseif Fox(i,1) == Fsi(ii,3)
                        inter_edge(i,1) = Fox(i,j);
                        inter_edge(i,2) = Fox(i,1);
                    end

                elseif j == 3
                    if Fox(i,1) == Fsi(ii,2)
                        inter_edge(i,1) = Fox(i,j);
                        inter_edge(i,2) = Fox(i,1);
                    elseif Fox(i,2) == Fsi(ii,1)
                        inter_edge(i,5) = Fox(i,j);
                        inter_edge(i,6) = Fox(i,2);
                    end
                end
            end
        end
    end
end
N = nnz(inter_edge); %% 0이 아닌 요소의 개수
save_edge = zeros(N/2,2); 
a=1;
for i=1:28
    if inter_edge(i,1)~=0
        save_edge(a,:) = inter_edge(i,:);
        a = a+1;
    end
end

t = linspace(0, 2);
r1 = 1; r2 = 2; r3 = 2.5;
r1x = r1*cospi(t); r2x = r2*cospi(t); r3x = r3*cospi(t);
r1y = r1*sinpi(t); r2y = r2*sinpi(t); r3y = r3*sinpi(t);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% visualizing
figure
% patch(r3x, r3y, 'k')
% patch(r2x, r2y, 'w')
% patch(r1x, r1y, 'k')
patch('Faces',F,'Vertices',V, 'FaceVertexCData',phi, 'EdgeColor','black','FaceColor','none','LineWidth',1, 'Marker','o');
patch('Faces',Fox,'Vertices',V, 'FaceVertexCData',phi, 'EdgeColor','black','FaceColor','red','LineWidth',1, 'Marker','o');
patch('Faces',Fsi,'Vertices',V, 'FaceVertexCData',phi, 'EdgeColor','black','FaceColor','blue','LineWidth',1, 'Marker','o');
patch('Faces',save_edge,'Vertices',V, 'EdgeColor','yellow','FaceColor','none','LineWidth', 3, 'Marker','o');
title('structure visualizing')
% % % plot(phi, 'b*')
% % % xlabel('vertex')
% % % ylabel('potential')
