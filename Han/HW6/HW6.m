%%% HW6_calculate per region_edit
%%% Seong-Min,Han
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
V = importdata("vertex.txt");
[Vrow, Vcol] = size(V); 
F = importdata("element.txt");
[Frow, Fcol] = size(F); 
contact = importdata("contact.txt");
F1 = importdata("element_region1.txt");
[F1row, F1col] = size(F1); 
F2 = importdata("element_region2.txt");
[F2row, F2col] = size(F2); 
F3 = importdata("element_region3.txt");
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
    area(i,1) = sqrt(norm_v21(i,1)*norm_v21(i,1)*norm_v31(i,1)*norm_v31(i,1) ...
        - dot(v21(i,:),v31(i,:))*dot(v21(i,:),v31(i,:)))/2;
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
e1 = 11.7; e2 = 3.9; e3 = 11.7; 
% for i=1:Frow
%     for j=1:3
%         if j==3
%             A(F(i,j),F(i,j)) = A(F(i,j),F(i,j)) - edge(i,j)/len(i,j);
%             A(F(i,j),F(i,1)) = A(F(i,j),F(i,1)) + edge(i,j)/len(i,j);
%         else
%             A(F(i,j),F(i,j)) = A(F(i,j),F(i,j)) - edge(i,j)/len(i,j);
%             A(F(i,j),F(i,j+1)) = A(F(i,j),F(i,j+1)) + edge(i,j)/len(i,j);
%         end
%     end
% end

edgeF1=zeros(F1row,3); lenF1=zeros(F1row,3);
edgeF2=zeros(F2row,3); lenF2=zeros(F2row,3);
edgeF3=zeros(F3row,3); lenF3=zeros(F3row,3);
aa=1; bb=1; cc=1;

for i=1:Frow
    for j=1:F1row
        if F(i,:) == F1(j,:)
            edgeF1(aa,:) = edge(i,:);
            lenF1(aa,:) = len(i,:);
            aa=aa+1;
        end
    end
    for j=1:F2row
        if F(i,:) == F2(j,:)
            edgeF2(bb,:) = edge(i,:);
            lenF2(bb,:) = len(i,:);
            bb=bb+1;
        end
    end
    for j=1:F3row
        if F(i,:) == F3(j,:)
            edgeF3(cc,:) = edge(i,:);
            lenF3(cc,:) = len(i,:);
            cc=cc+1;
        end
    end
end
for i=1:F1row
    for j=1:3
        if j==3
            A(F1(i,j),F1(i,j)) = A(F1(i,j),F1(i,j)) - e1*edgeF1(i,j)/lenF1(i,j);
            A(F1(i,j),F1(i,1)) = A(F1(i,j),F1(i,1)) + e1*edgeF1(i,j)/lenF1(i,j);
        else
            A(F1(i,j),F1(i,j)) = A(F1(i,j),F1(i,j)) - e1*edgeF1(i,j)/lenF1(i,j);
            A(F1(i,j),F1(i,j+1)) = A(F1(i,j),F1(i,j+1)) + e1*edgeF1(i,j)/lenF1(i,j);
        end
    end
end

for i=1:F2row
    for j=1:3
        if j==3
            A(F2(i,j),F2(i,j)) = A(F2(i,j),F2(i,j)) - e2*edgeF2(i,j)/lenF2(i,j);
            A(F2(i,j),F2(i,1)) = A(F2(i,j),F2(i,1)) + e2*edgeF2(i,j)/lenF2(i,j);
        else
            A(F2(i,j),F2(i,j)) = A(F2(i,j),F2(i,j)) - e2*edgeF2(i,j)/lenF2(i,j);
            A(F2(i,j),F2(i,j+1)) = A(F2(i,j),F2(i,j+1)) + e2*edgeF2(i,j)/lenF2(i,j);
        end
    end
end

for i=1:F3row
    for j=1:3
        if j==3
            A(F3(i,j),F3(i,j)) = A(F3(i,j),F3(i,j)) - e3*edgeF3(i,j)/lenF3(i,j);
            A(F3(i,j),F3(i,1)) = A(F3(i,j),F3(i,1)) + e3*edgeF3(i,j)/lenF3(i,j);
        else
            A(F3(i,j),F3(i,j)) = A(F3(i,j),F3(i,j)) - e3*edgeF3(i,j)/lenF3(i,j);
            A(F3(i,j),F3(i,j+1)) = A(F3(i,j),F3(i,j+1)) + e3*edgeF3(i,j)/lenF3(i,j);
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

con1 = 1;
con2 = 1;
con3 = 1;

b1 = zeros(Vrow,1);
for i=1:Ncrow
    for j=1:Nccol
        if i>=1 && i<=con1
            b1(contact(i,j),1) = 1;
        elseif i>con1 &&  i<=con1+con2    
            b1(contact(i,j),1) = 0;
        elseif i>con1+con2 && i<=con1+con2+con3
            b1(contact(i,j),1) = 0;
        end
    end
end
phi1 = A\b1;

b2 = zeros(Vrow,1);
for i=1:Ncrow
    for j=1:Nccol
        if i>=1 && i<=con1
            b2(contact(i,j),1) = 0;
        elseif i>con1 &&  i<=con1+con2    
            b2(contact(i,j),1) = 1;
        elseif i>con1+con2 && i<=con1+con2+con3
            b2(contact(i,j),1) = 0;
        end
    end
end
phi2 = A\b2;

b3 = zeros(Vrow,1);
for i=1:Ncrow
    for j=1:Nccol
        if i>=1 && i<=con1
            b3(contact(i,j),1) = 0;
        elseif i>con1 &&  i<=con1+con2    
            b3(contact(i,j),1) = 0;
        elseif i>con1+con2 && i<=con1+con2+con3
            b3(contact(i,j),1) = 1;
        end
    end
end
phi3 = A\b3;

sumphi = phi1 + phi2 + phi3;

b4 = zeros(Vrow,1);
for i=1:Ncrow
    for j=1:Nccol
        if i>=1 && i<=con1
            b4(contact(i,j),1) = 1;
        elseif i>con1 &&  i<=con1+con2    
            b4(contact(i,j),1) = 1;
        elseif i>con1+con2 && i<=con1+con2+con3
            b4(contact(i,j),1) = 1;
        end
    end
end
phi4 = A\b4;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% interface edge
is1=intersect(F1, F2);
is1row = size(is1,1);
interedge1 = zeros(is1row-1,2);
j=1;
for i = 1:F1row
    if size(intersect(F1(i,:),is1'),2) == 2
        interedge1(j,:) = intersect(F1(i,:),is1');
        j=j+1;
    end
end

is2=intersect(F2, F3);
is2row = size(is2,1);
interedge2 = zeros(is2row-1,2);
j=1;
for i = 1:F2row
    if size(intersect(F2(i,:),is2'),2) == 2
        interedge2(j,:) = intersect(F2(i,:),is2');
        j=j+1;
    end
end

is3=intersect(F3, F1);
is3row = size(is3, 1);
interedge3 = zeros(is3row-1,2);
j=1;
for i = 1:F3row
    if size(intersect(F3(i,:),is3'),2) == 2
        interedge3(j,:) = intersect(F3(i,:),is3');
        j=j+1;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% visualizing
%%% potential & interfaces
figure
subplot(221)
patch('Faces',F,'Vertices',V, 'FaceVertexCData',phi1, 'EdgeColor','black','FaceColor','interp','LineWidth',1, 'Marker','o');
title('Φ1 visualizing')
colorbar
subplot(222)
patch('Faces',F,'Vertices',V, 'FaceVertexCData',phi2, 'EdgeColor','black','FaceColor','interp','LineWidth',1, 'Marker','o');
title('Φ2 visualizing')
colorbar
subplot(223)
patch('Faces',F,'Vertices',V, 'FaceVertexCData',phi3, 'EdgeColor','black','FaceColor','interp','LineWidth',1, 'Marker','o');
title('Φ3 visualizing')
colorbar
subplot(224)
patch('Faces',F,'Vertices',V, 'FaceVertexCData',phi4, 'EdgeColor','black','FaceColor','interp','LineWidth',1, 'Marker','o');
patch('Faces',interedge1,'Vertices',V, 'EdgeColor','red','FaceColor','none','LineWidth', 3, 'Marker','o');
patch('Faces',interedge2,'Vertices',V, 'EdgeColor','magenta','FaceColor','none','LineWidth', 3, 'Marker','o');
patch('Faces',interedge3,'Vertices',V, 'EdgeColor','yellow','FaceColor','none','LineWidth', 3, 'Marker','o');
colorbar
title('Φ4 visualizing')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% structure & interfaces
% figure
% patch('Faces',F,'Vertices',V, 'FaceVertexCData',phi, 'EdgeColor','black','FaceColor','interp','LineWidth',1, 'Marker','o');
% patch('Faces',F1,'Vertices',V, 'FaceVertexCData',phi, 'EdgeColor','black','FaceColor','green','LineWidth',1, 'Marker','o');
% patch('Faces',F2,'Vertices',V, 'FaceVertexCData',phi, 'EdgeColor','black','FaceColor','cyan','LineWidth',1, 'Marker','o');
% patch('Faces',F3,'Vertices',V, 'FaceVertexCData',phi, 'EdgeColor','black','FaceColor','blue','LineWidth',1, 'Marker','o');
% patch('Faces',save_edge1,'Vertices',V, 'EdgeColor','red','FaceColor','none','LineWidth', 3, 'Marker','o');
% patch('Faces',save_edge2,'Vertices',V, 'EdgeColor','magenta','FaceColor','none','LineWidth', 3, 'Marker','o');
% patch('Faces',save_edge3,'Vertices',V, 'EdgeColor','yellow','FaceColor','none','LineWidth', 3, 'Marker','o');
% title('structure & interfaces')
