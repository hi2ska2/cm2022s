%%% HW5_region_divide
%%% written by Seong-Min, Han
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
V = importdata("vertex.txt");
[Vrow, Vcol] = size(V); 
F = importdata("element.txt");
[Frow, Fcol] = size(F); 
contact = importdata("contact.txt");
F1 = importdata("element2_region1.txt");
[F1row, F1col] = size(F1); 
F2 = importdata("element2_region2.txt");
[F2row, F2col] = size(F2); 
F3 = importdata("element2_region3.txt");
[F3row, F3col] = size(F3); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% vertex와 vertex 사이의 length
len = zeros(Frow,3);
for i=1:Frow
    for j =1:3
        if j==3
            len(i,j) = norm(V(F(i,j),:)-V(F(i,1),:));
        else
            len(i,j) = norm(V(F(i,j),:)-V(F(i,j+1),:));
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% vertex 1->2, 1->3을 연결하는 vector와 vector의 norm 값
v21 = zeros(Frow,2);
v31 = zeros(Frow,2);
norm_v21 = zeros(Frow,1);
norm_v31 = zeros(Frow,1);
for i=1:Frow
    v21(i,:) = V(F(i,2),:)-V(F(i,1),:);
    v31(i,:) = V(F(i,3),:)-V(F(i,1),:);
    norm_v21(i,1) = norm(v21(i,:));
    norm_v31(i,1) = norm(v31(i,:));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%element로 이루어진 face의 area와 face의 외접원 radius
area = zeros(Frow,1);
rad = zeros(Frow,1);
for i=1:Frow
    area(i,1) = sqrt(norm_v21(i,1)*norm_v21(i,1)*norm_v31(i,1)*norm_v31(i,1) - dot(v21(i,:),v31(i,:))*dot(v21(i,:),v31(i,:)))/2;
    rad(i,1) = (len(i,1)*len(i,2)*len(i,3)) / (4*area(i,1));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%face의 edge / edge1(1-2), edge2(2-3), edge3(3-1)
edge = zeros(Frow,1);
for i=1:Frow
    for j=1:3
        edge(i,j) = real(sqrt(rad(i,1)*rad(i,1) - (len(i,j)/2)*(len(i,j)/2)));
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% check contact
checkcontact = zeros(Frow,3);
[Ncrow, Nccol] = size(contact); 
for i=1:Frow
    for ii=1:Ncrow
        for j=1:3
            for jj=1:Nccol
                if contact(ii,jj) == F(i,j)
                   if jj==1
                       if j==1
                           if contact(ii,2) == F(i,2)
                               checkcontact(i,j) = 1;
                           elseif contact(ii,2) == F(i,3)
                               checkcontact(i,j) = 1;
                           end
                       elseif j==2
                           if contact(ii,2) == F(i,3)
                               checkcontact(i,j) = 1;
                           elseif contact(ii,2) == F(i,1)
                               checkcontact(i,j) = 1;
                           end
                       elseif j==3
                           if contact(ii,2) == F(i,1)
                               checkcontact(i,j) = 1;
                           elseif contact(ii,2) == F(i,2)
                               checkcontact(i,j) = 1;
                           end
                       end
                   elseif jj==2
                       if j==1
                           if contact(ii,1) == F(i,2)
                               checkcontact(i,j) = 1;
                           elseif contact(ii,1) == F(i,3)
                               checkcontact(i,1) = 1;
                           end
                       elseif j==2
                           if contact(ii,1) == F(i,3)
                               checkcontact(i,j) = 1;
                           elseif contact(ii,1) == F(i,1)
                               checkcontact(i,j) = 1;
                           end
                       elseif j==3
                           if contact(ii,1) == F(i,1)
                               checkcontact(i,j) = 1;
                           elseif contact(ii,1) ==  F(i,2)
                               checkcontact(i,j) = 1;
                           end
                       end
                   end
                end
            end
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% A, b, phi matrix
A = zeros(Vrow,Vrow);
for i=1:Frow
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

for i=1:Frow
    for j=1:3
        if checkcontact(i,j)==1
            A(F(i,j),:)= 0;
            A(F(i,j),F(i,j))=1;
        end
    end
end

phi1 = 1;
phi2 = 2;
phi3 = 1;
phi4 = 10;

b = zeros(Vrow,1);
for i=1:Ncrow
    for j=1:Nccol
        if i>=1 && i<=phi1
            b(contact(i,j),1) = 1;
        elseif i>phi1 &&  i<=phi1+phi2    
            b(contact(i,j),1) = 0;
        elseif i>phi1+phi2 && i<=phi1+phi2+phi3
            b(contact(i,j),1) = 1.5;
        elseif i>phi1+phi2+phi3 && i<=phi1+phi2+phi3+phi4
            b(contact(i,j),1) = 0;
        end
    end
end

phi = A\b;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% interface edge 1 (between region1 &  region2)
inter_edge1 = zeros(F1row,2);
for i=1:F1row
    for ii=1:F2row
        for j=1:3
            if F1(i,j) == F2(ii,j)
                if j == 1
                    if F1(i,2) == F2(ii,3)
                        inter_edge1(i,1) = F1(i,j);
                        inter_edge1(i,2) = F1(i,2);
                    elseif F1(i,3) == F2(ii,2)
                        inter_edge1(i,1) = F1(i,j);
                        inter_edge1(i,2) = F1(i,3);
                    end
                elseif j == 2
                    if F1(i,3) == F2(ii,1)
                        inter_edge1(i,1) = F1(i,j);
                        inter_edge1(i,2) = F1(i,3);
                    elseif F1(i,1) == F2(ii,3)
                        inter_edge1(i,1) = F1(i,j);
                        inter_edge1(i,2) = F1(i,1);
                    end
                elseif j == 3
                    if F1(i,1) == F2(ii,2)
                        inter_edge1(i,1) = F1(i,j);
                        inter_edge1(i,2) = F1(i,1);
                    elseif F1(i,2) == F2(ii,1)
                        inter_edge1(i,1) = F1(i,j);
                        inter_edge1(i,2) = F1(i,2);
                    end
                end
            end
        end
    end
end

N1 = nnz(inter_edge1); %% 0이 아닌 요소의 개수
save_edge1 = zeros(N1/2,2); 
aa=1;
for i=1:F1row
    if inter_edge1(i,1)~=0
        save_edge1(aa,:) = inter_edge1(i,:);
        aa = aa+1;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% interface edge 2 (between region2 &  region3)
inter_edge2 = zeros(F2row,2);
for i=1:F2row
    for ii=1:F3row
        for j=1:3
            if F2(i,j) == F3(ii,j)
                if j == 1
                    if F2(i,2) == F3(ii,3)
                        inter_edge2(i,1) = F2(i,j);
                        inter_edge2(i,2) = F2(i,2);
                    elseif F2(i,3) == F3(ii,2)
                        inter_edge1(i,1) = F2(i,j);
                        inter_edge1(i,2) = F2(i,3);
                    end
                elseif j == 2
                    if F2(i,3) == F3(ii,1)
                        inter_edge2(i,1) = F2(i,j);
                        inter_edge2(i,2) = F2(i,3);
                    elseif F2(i,1) == F3(ii,3)
                        inter_edge2(i,1) = F2(i,j);
                        inter_edge2(i,2) = F2(i,1);
                    end
                elseif j == 3
                    if F2(i,1) == F3(ii,2)
                        inter_edge2(i,1) = F2(i,j);
                        inter_edge2(i,2) = F2(i,1);
                    elseif F2(i,2) == F3(ii,1)
                        inter_edge2(i,1) = F2(i,j);
                        inter_edge2(i,2) = F2(i,2);
                    end
                end
            end
        end
    end
end

N2 = nnz(inter_edge2); %% 0이 아닌 요소의 개수
save_edge2 = zeros(N2/2,2); 
aa=1;
for i=1:F2row
    if inter_edge2(i,1)~=0
        save_edge2(aa,:) = inter_edge2(i,:);
        aa = aa+1;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% interface edge 3 (between region3 & region1)
inter_edge3 = zeros(F3row,2);
for i=1:F3row
    for ii=1:F1row
        for j=1:3
            if F3(i,j) == F1(ii,j)
                if j == 1
                    if F3(i,2) == F1(ii,3)
                        inter_edge3(i,1) = F3(i,j);
                        inter_edge3(i,2) = F3(i,2);
                    elseif F3(i,3) == F1(ii,2)
                        inter_edge3(i,1) = F3(i,j);
                        inter_edge3(i,2) = F3(i,3);
                    end
                elseif j == 2
                    if F3(i,3) == F1(ii,1)
                        inter_edge3(i,1) = F3(i,j);
                        inter_edge3(i,2) = F3(i,3);
                    elseif F3(i,1) == F1(ii,3)
                        inter_edge3(i,1) = F3(i,j);
                        inter_edge3(i,2) = F3(i,1);
                    end
                elseif j == 3
                    if F3(i,1) == F1(ii,2)
                        inter_edge3(i,1) = F3(i,j);
                        inter_edge3(i,2) = F3(i,1);
                    elseif F3(i,2) == F1(ii,1)
                        inter_edge3(i,1) = F3(i,j);
                        inter_edge3(i,2) = F3(i,2);
                    end
                end
            end
        end
    end
end

N3 = nnz(inter_edge3); %% 0이 아닌 요소의 개수
save_edge3 = zeros(N3/2,2); 
aa=1;
for i=1:F3row
    if inter_edge3(i,1)~=0
        save_edge3(aa,:) = inter_edge3(i,:);
        aa = aa+1;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% visualizing
figure
%%% potential
patch('Faces',F,'Vertices',V, 'FaceVertexCData',phi, 'EdgeColor','black','FaceColor','none','LineWidth',1, 'Marker','o');
patch('Faces',F1,'Vertices',V, 'FaceVertexCData',phi, 'EdgeColor','black','FaceColor','interp','LineWidth',1, 'Marker','o');
patch('Faces',F2,'Vertices',V, 'FaceVertexCData',phi, 'EdgeColor','black','FaceColor','interp','LineWidth',1, 'Marker','o');
patch('Faces',F3,'Vertices',V, 'FaceVertexCData',phi, 'EdgeColor','black','FaceColor','interp','LineWidth',1, 'Marker','o');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% structure
% patch('Faces',F,'Vertices',V, 'FaceVertexCData',phi, 'EdgeColor','black','FaceColor','interp','LineWidth',1, 'Marker','o');
% patch('Faces',F1,'Vertices',V, 'FaceVertexCData',phi, 'EdgeColor','black','FaceColor','green','LineWidth',1, 'Marker','o');
% patch('Faces',F2,'Vertices',V, 'FaceVertexCData',phi, 'EdgeColor','black','FaceColor','cyan','LineWidth',1, 'Marker','o');
% patch('Faces',F3,'Vertices',V, 'FaceVertexCData',phi, 'EdgeColor','black','FaceColor','blue','LineWidth',1, 'Marker','o');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% interface
patch('Faces',save_edge1,'Vertices',V, 'EdgeColor','red','FaceColor','none','LineWidth', 3, 'Marker','o');
patch('Faces',save_edge2,'Vertices',V, 'EdgeColor','magenta','FaceColor','none','LineWidth', 3, 'Marker','o');
patch('Faces',save_edge3,'Vertices',V, 'EdgeColor','yellow','FaceColor','none','LineWidth', 3, 'Marker','o');

title('structure')
colorbar
