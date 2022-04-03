%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% HW7_reupload
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
V = importdata("vertex.txt"); Vrow = size(V,1);
F = importdata("element.txt"); Frow = size(F,1);
F1 = importdata("element_region1.txt"); F1row = size(F1,1);
F2 = importdata("element_region2.txt"); F2row = size(F2,1);
F3 = importdata("element_region3.txt"); F3row = size(F3,1);
contact = importdata("contact.txt"); [Ncrow, Nccol] = size(contact);
VF = unique(F); VFrow= size(VF,1);
VF1 = unique(F1); VF1row = size(VF1, 1);
VF2 = unique(F2); VF2row = size(VF2, 1);
VF3 = unique(F3); VF3row = size(VF3, 1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% element 기반 vertex-vertex length
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
%%% vector % norm vector
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
%%% 외접원의 radius
area = zeros(Frow,1);
rad = zeros(Frow,1);
for i=1:Frow
    area(i,1) = sqrt(norm_v21(i,1)*norm_v21(i,1)*norm_v31(i,1)*norm_v31(i,1) ...
        - dot(v21(i,:),v31(i,:))*dot(v21(i,:),v31(i,:)))/2;
    rad(i,1) = (len(i,1)*len(i,2)*len(i,3)) / (4*area(i,1));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% element 기반 edge
edge = zeros(Frow,1);
for i=1:Frow
    for j=1:3
        edge(i,j) = real(sqrt(rad(i,1)*rad(i,1) - (len(i,j)/2)*(len(i,j)/2)));
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% vertex 기반 total edge & length
tot_edge = zeros(Vrow,Vrow);
tot_len = zeros(Vrow,Vrow);
for i=1:Frow
    for j=1:3
        if j==3
            tot_edge(F(i,j),F(i,1)) = edge(i,3);
            tot_len(F(i,j),F(i,1)) = len(i,3);
        else
            tot_edge(F(i,j),F(i,j+1)) = edge(i,j);
            tot_len(F(i,j),F(i,j+1)) = len(i,j);
        end
    end
end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% A, b, phi matrix
A = zeros(Vrow,Vrow);
e1 = 11.7;
e2 = 3.9;
e3 = 11.7;
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
        if j==1
            A(F1(i,j),F1(i,j)) = A(F1(i,j),F1(i,j)) - e1*edgeF1(i,1)/lenF1(i,1) - e1*edgeF1(i,3)/lenF1(i,3);
            A(F1(i,j),F1(i,2)) = A(F1(i,j),F1(i,2)) + e1*edgeF1(i,1)/lenF1(i,1);
            A(F1(i,j),F1(i,3)) = A(F1(i,j),F1(i,3)) + e1*edgeF1(i,3)/lenF1(i,3);
        elseif j==2
            A(F1(i,j),F1(i,j)) = A(F1(i,j),F1(i,j)) - e1*edgeF1(i,2)/lenF1(i,2) - e1*edgeF1(i,1)/lenF1(i,1);
            A(F1(i,j),F1(i,3)) = A(F1(i,j),F1(i,3)) + e1*edgeF1(i,2)/lenF1(i,2);
            A(F1(i,j),F1(i,1)) = A(F1(i,j),F1(i,1)) + e1*edgeF1(i,1)/lenF1(i,1);
        elseif j==3
            A(F1(i,j),F1(i,j)) = A(F1(i,j),F1(i,j)) - e1*edgeF1(i,3)/lenF1(i,3) - e1*edgeF1(i,2)/lenF1(i,2);
            A(F1(i,j),F1(i,1)) = A(F1(i,j),F1(i,1)) + e1*edgeF1(i,3)/lenF1(i,3);
            A(F1(i,j),F1(i,2)) = A(F1(i,j),F1(i,2)) + e1*edgeF1(i,2)/lenF1(i,2);        
        end
    end
end

for i=1:F2row
    for j=1:3
        if j==1
            A(F2(i,j),F2(i,j)) = A(F2(i,j),F2(i,j)) - e2*edgeF2(i,1)/lenF2(i,1) - e2*edgeF2(i,3)/lenF2(i,3);
            A(F2(i,j),F2(i,2)) = A(F2(i,j),F2(i,2)) + e2*edgeF2(i,1)/lenF2(i,1);
            A(F2(i,j),F2(i,3)) = A(F2(i,j),F2(i,3)) + e2*edgeF2(i,3)/lenF2(i,3);  
        elseif j==2
            A(F2(i,j),F2(i,j)) = A(F2(i,j),F2(i,j)) - e2*edgeF2(i,2)/lenF2(i,2) - e2*edgeF2(i,1)/lenF2(i,1);
            A(F2(i,j),F2(i,3)) = A(F2(i,j),F2(i,3)) + e2*edgeF2(i,2)/lenF2(i,2);
            A(F2(i,j),F2(i,1)) = A(F2(i,j),F2(i,1)) + e2*edgeF2(i,1)/lenF2(i,1);  
        elseif j==3
            A(F2(i,j),F2(i,j)) = A(F2(i,j),F2(i,j)) - e2*edgeF2(i,3)/lenF2(i,3) - e2*edgeF2(i,2)/lenF2(i,2);
            A(F2(i,j),F2(i,1)) = A(F2(i,j),F2(i,1)) + e2*edgeF2(i,3)/lenF2(i,3);
            A(F2(i,j),F2(i,2)) = A(F2(i,j),F2(i,2)) + e2*edgeF2(i,2)/lenF2(i,2);  
        end
    end
end

for i=1:F3row
    for j=1:3
        if j==1
            A(F3(i,j),F3(i,j)) = A(F3(i,j),F3(i,j)) - e3*edgeF3(i,1)/lenF3(i,1) - e3*edgeF3(i,3)/lenF3(i,3);
            A(F3(i,j),F3(i,2)) = A(F3(i,j),F3(i,2)) + e3*edgeF3(i,1)/lenF3(i,1);
            A(F3(i,j),F3(i,3)) = A(F3(i,j),F3(i,3)) + e3*edgeF3(i,3)/lenF3(i,3);
        elseif j==2
            A(F3(i,j),F3(i,j)) = A(F3(i,j),F3(i,j)) - e3*edgeF3(i,2)/lenF3(i,2) - e3*edgeF3(i,1)/lenF3(i,1);
            A(F3(i,j),F3(i,3)) = A(F3(i,j),F3(i,3)) + e3*edgeF3(i,2)/lenF3(i,2);
            A(F3(i,j),F3(i,1)) = A(F3(i,j),F3(i,1)) + e3*edgeF3(i,1)/lenF3(i,1);
        elseif j==3
            A(F3(i,j),F3(i,j)) = A(F3(i,j),F3(i,j)) - e3*edgeF3(i,3)/lenF3(i,3) - e3*edgeF3(i,2)/lenF3(i,2);
            A(F3(i,j),F3(i,1)) = A(F3(i,j),F3(i,1)) + e3*edgeF3(i,3)/lenF3(i,3);
            A(F3(i,j),F3(i,2)) = A(F3(i,j),F3(i,2)) + e3*edgeF3(i,2)/lenF3(i,2);
        end
    end
end

for i=1:Ncrow
    for j=1:Nccol
        A(contact(i,j),:)= 0;
        A(contact(i,j),contact(i,j))=1;
    end
end

con1 = 6; % top contact
con2 = 6; % bottom contact

b = zeros(Vrow,1);
for i=1:Ncrow
    for j=1:Nccol
        if i>=1 && i<=con1
            b(contact(i,j),1) = 1;
        elseif i>con1 &&  i<=con1+con2
            b(contact(i,j),1) = 0;
        end
    end
end

phi = A\b;

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

is = zeros(is1row+is2row+is3row,1);
isrow = size(is, 1);
for i=1:size(is,1)
    if i<=is1row
        is(i,1)=is1(i,1);
    elseif i>is1row && i<=(is1row+is2row)
        is(i,1)=is2(i-is1row,1);
    elseif is>(is1row+is2row) && i<=(is1row+is2row+is3row)
        is(i,1)=is2(i-(is1row+is2row),1);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
aaa = zeros(VF1row+VF2row+VF3row,1); %%% variable for total regions
VF_entire = zeros(size(aaa,1),1);
for i=1:size(aaa,1)
    if i<=VF1row
        aaa(i)=phi(VF1(i));
        VF_entire(i)=VF1(i);
    elseif i>VF1row && i<=VF1row+VF2row
        aaa(i)=phi(VF2(i-VF1row));
        VF_entire(i)=VF2(i-VF1row);
    else
        aaa(i)=phi(VF3(i-(VF1row+VF2row)));
        VF_entire(i)=VF3(i-(VF1row+VF2row));
    end
end

bbb = zeros(VF2row,1); %%% variable for semiconductor region
for i=1:size(bbb,1)
    bbb(i)=phi(VF2(i));
end

ccc = zeros(VF1row+VF3row,1); %%% variable for insulator region
VF_insulator=zeros(size(ccc,1),1);
for i=1:size(ccc,1)
    if i<=VF1row
        ccc(i)=phi(VF1(i));
        VF_insulator(i)=VF1(i);
    else
        ccc(i)=phi(VF3(i-VF1row));
        VF_insulator(i)=VF3(i-VF1row);
    end
end

%%%  solution X matrix 
X = zeros(2*size(aaa,1),1);
VF_entire = zeros(size(aaa,1),1);
for i=1:size(aaa,1)
    if i <= VF1row 
        X(2*i-1) = aaa(i);
        X(2*i) = ccc(i);
    elseif i > VF1row && i <= (VF1row + VF2row)
        X(2*i-1) = aaa(i);
        X(2*i) = bbb(i-VF1row);
    elseif i > (VF1row + VF2row) && i <= (VF1row + VF2row + VF3row)
        X(2*i-1) = aaa(i);
        X(2*i) = ccc(i-VF2row);  
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% scan part
while 1
    reg = input('region?\n (1)oxide (2)silicon\n : ');
    var = input('variable?\n aaa(entire)\n bbb(only semicinductor)\n ccc(only insulator)\n : ', 's');
    if (reg==1 && strcmp(var, 'bbb')) || (reg==2 && strcmp(var, 'ccc'))
        fprintf('Error!\n');
        continue
    end
    ver = input('vertex number? : ');
%%%%% print part
    if reg==1 
        if strcmp(var, 'aaa')
            if size(find(ver==VF1),1)==1
                xindex = 2*find(ver==VF1)-1;
                fprintf('X[%d] = %f\n',xindex, X(xindex));
            elseif size(find(ver==VF3),1)==1
                xindex = 2*(VF1row+VF2row)+2*find(ver==VF3)-1;
                fprintf('X[%d] = %f\n',xindex, X(xindex));  
            else 
                fprintf('Error !\n');
            end
        elseif strcmp(var, 'ccc')
            if size(find(ver==VF1),1)==1
                xindex = 2*find(ver==VF1);
                fprintf('X[%d] = %f\n',xindex, X(xindex));
            elseif size(find(ver==VF3),1)==1
                 xindex = 2*(VF1row+VF2row)+2*find(ver==VF3);
                fprintf('X[%d] = %f\n',xindex, X(xindex));              
            else 
                fprintf('Error !\n');
            end       
        end 
    elseif reg==2 && size(find(ver==VF2),1)==1
        if strcmp(var, 'aaa')
            xindex = 2*VF1row+2*find(ver==VF2)-1;
            fprintf('X[%d] = %f\n',xindex, X(xindex));
        elseif strcmp(var, 'bbb')
            xindex = 2*VF1row+2*find(ver==VF2);
            fprintf('X[%d] = %f\n',xindex, X(xindex));
        else 
                fprintf('Error !\n');
        end
    else
        fprintf('Error !\n');
    end
    break
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% visualizing
% figure
% patch('Faces',F1,'Vertices',V, 'FaceVertexCData',phi, 'EdgeColor','blue','FaceColor','interp','LineWidth',1, 'Marker','o');
% patch('Faces',F2,'Vertices',V, 'FaceVertexCData',phi, 'EdgeColor','green','FaceColor','interp','LineWidth',1, 'Marker','o');
% patch('Faces',F3,'Vertices',V, 'FaceVertexCData',phi, 'EdgeColor','blue','FaceColor','interp','LineWidth',1, 'Marker','o');
% title('Φ visualizing structure')
% patch('Faces',interedge1,'Vertices',V, 'EdgeColor','magenta','FaceColor','none','LineWidth', 3, 'Marker','o');
% patch('Faces',interedge2,'Vertices',V, 'EdgeColor','magenta','FaceColor','none','LineWidth', 3, 'Marker','o');
% patch('Faces',interedge3,'Vertices',V, 'EdgeColor','magenta','FaceColor','none','LineWidth', 3, 'Marker','o');
% colorbar
