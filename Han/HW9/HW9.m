%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% HW9
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
VFT = zeros(VF1row+VF2row+VF3row,1);
VFTrow = size(VFT,1);
for i=1:size(VFT,1)
    if i<=VF1row
        VFT(i)=VF1(i);
    elseif i>VF1row && i<=VF1row+VF2row
        VFT(i)=VF2(i-VF1row);
    else
        VFT(i)=VF3(i-(VF1row+VF2row));
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
VF1_aaa = zeros(VF1row,2);
aa=1;
for i=1:VF1row
    VF1_aaa(i,1)=VF1(i);
    VF1_aaa(i,2)=aa;
    aa=aa+1;
end
VF2_aaa = zeros(VF2row,2);
bb=VF1row+1;
for i=1:VF2row
    VF2_aaa(i,1)=VF2(i);
    VF2_aaa(i,2)=bb;
    bb=bb+1;
end
VF3_aaa = zeros(VF3row,2);
cc=VF1row+VF2row+1;
for i=1:VF3row
    VF3_aaa(i,1)=VF3(i);
    VF3_aaa(i,2)=cc;
    cc=cc+1;
end
VF2_bb1 = zeros(VF2row,1);
bb1 = VF1row+VF2row+VF3row+1;
for i=1:VF2row
    VF2_bb1(i,1)=VF2(i);
    VF2_bb1(i,2)=bb1;
    bb1=bb1+1;
end
VF2_bb2 = zeros(VF2row,1);
bb2 = VF1row+VF2row+VF3row+VF2row+1;
for i=1:VF2row
    VF2_bb2(i,1)=VF2(i);
    VF2_bb2(i,2)=bb2;
    bb2=bb2+1;
end
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
%%% A, b, phi matrix
jac = zeros(VFTrow,VFTrow);
e1 = 11.7;  e2 = 3.9;  e3 = 11.7;
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
        fa = find(F1(i,1)==VF1_aaa(:,1));
        fb = find(F1(i,2)==VF1_aaa(:,1));
        fc = find(F1(i,3)==VF1_aaa(:,1));
        if j==1
            jac(VF1_aaa(fa,2),VF1_aaa(fa,2)) = jac(VF1_aaa(fa,2),VF1_aaa(fa,2)) - e1*edgeF1(i,1)/lenF1(i,1) - e1*edgeF1(i,3)/lenF1(i,3);
            jac(VF1_aaa(fa,2),VF1_aaa(fb,2)) = jac(VF1_aaa(fa,2),VF1_aaa(fb,2)) + e1*edgeF1(i,1)/lenF1(i,1);
            jac(VF1_aaa(fa,2),VF1_aaa(fc,2)) = jac(VF1_aaa(fa,2),VF1_aaa(fc,2)) + e1*edgeF1(i,3)/lenF1(i,3);
        elseif j==2
            jac(VF1_aaa(fb,2),VF1_aaa(fb,2)) = jac(VF1_aaa(fb,2),VF1_aaa(fb,2)) - e1*edgeF1(i,2)/lenF1(i,2) - e1*edgeF1(i,1)/lenF1(i,1);
            jac(VF1_aaa(fb,2),VF1_aaa(fc,2)) = jac(VF1_aaa(fb,2),VF1_aaa(fc,2)) + e1*edgeF1(i,2)/lenF1(i,2);
            jac(VF1_aaa(fb,2),VF1_aaa(fa,2)) = jac(VF1_aaa(fb,2),VF1_aaa(fa,2)) + e1*edgeF1(i,1)/lenF1(i,1);
        elseif j==3
            jac(VF1_aaa(fc,2),VF1_aaa(fc,2)) = jac(VF1_aaa(fc,2),VF1_aaa(fc,2)) - e1*edgeF1(i,3)/lenF1(i,3) - e1*edgeF1(i,2)/lenF1(i,2);
            jac(VF1_aaa(fc,2),VF1_aaa(fa,2)) = jac(VF1_aaa(fc,2),VF1_aaa(fa,2)) + e1*edgeF1(i,3)/lenF1(i,3);
            jac(VF1_aaa(fc,2),VF1_aaa(fb,2)) = jac(VF1_aaa(fc,2),VF1_aaa(fb,2)) + e1*edgeF1(i,2)/lenF1(i,2);        
        end
    end
end

for i=1:F2row
    for j=1:3
        fa = find(F2(i,1)==VF2_aaa(:,1));
        fb = find(F2(i,2)==VF2_aaa(:,1));
        fc = find(F2(i,3)==VF2_aaa(:,1));
        if j==1
            jac(VF2_aaa(fa,2),VF2_aaa(fa,2)) = jac(VF2_aaa(fa,2),VF2_aaa(fa,2)) - e2*edgeF2(i,1)/lenF2(i,1) - e2*edgeF2(i,3)/lenF2(i,3);
            jac(VF2_aaa(fa,2),VF2_aaa(fb,2)) = jac(VF2_aaa(fa,2),VF2_aaa(fb,2)) + e2*edgeF2(i,1)/lenF2(i,1);
            jac(VF2_aaa(fa,2),VF2_aaa(fc,2)) = jac(VF2_aaa(fa,2),VF2_aaa(fc,2)) + e2*edgeF2(i,3)/lenF2(i,3);  
        elseif j==2
            jac(VF2_aaa(fb,2),VF2_aaa(fb,2)) = jac(VF2_aaa(fb,2),VF2_aaa(fb,2)) - e2*edgeF2(i,2)/lenF2(i,2) - e2*edgeF2(i,1)/lenF2(i,1);
            jac(VF2_aaa(fb,2),VF2_aaa(fc,2)) = jac(VF2_aaa(fb,2),VF2_aaa(fc,2)) + e2*edgeF2(i,2)/lenF2(i,2);
            jac(VF2_aaa(fb,2),VF2_aaa(fa,2)) = jac(VF2_aaa(fb,2),VF2_aaa(fa,2)) + e2*edgeF2(i,1)/lenF2(i,1);  
        elseif j==3
            jac(VF2_aaa(fc,2),VF2_aaa(fc,2)) = jac(VF2_aaa(fc,2),VF2_aaa(fc,2)) - e2*edgeF2(i,3)/lenF2(i,3) - e2*edgeF2(i,2)/lenF2(i,2);
            jac(VF2_aaa(fc,2),VF2_aaa(fa,2)) = jac(VF2_aaa(fc,2),VF2_aaa(fa,2)) + e2*edgeF2(i,3)/lenF2(i,3);
            jac(VF2_aaa(fc,2),VF2_aaa(fb,2)) = jac(VF2_aaa(fc,2),VF2_aaa(fb,2)) + e2*edgeF2(i,2)/lenF2(i,2);  
        end
    end
end

for i=1:F3row
    for j=1:3
        fa = find(F3(i,1)==VF3_aaa(:,1));
        fb = find(F3(i,2)==VF3_aaa(:,1));
        fc = find(F3(i,3)==VF3_aaa(:,1));
        if j==1
            jac(VF3_aaa(fa,2),VF3_aaa(fa,2)) = jac(VF3_aaa(fa,2),VF3_aaa(fa,2)) - e3*edgeF3(i,1)/lenF3(i,1) - e3*edgeF3(i,3)/lenF3(i,3);
            jac(VF3_aaa(fa,2),VF3_aaa(fb,2)) = jac(VF3_aaa(fa,2),VF3_aaa(fb,2)) + e3*edgeF3(i,1)/lenF3(i,1);
            jac(VF3_aaa(fa,2),VF3_aaa(fc,2)) = jac(VF3_aaa(fa,2),VF3_aaa(fc,2)) + e3*edgeF3(i,3)/lenF3(i,3);
        elseif j==2
            jac(VF3_aaa(fb,2),VF3_aaa(fb,2)) = jac(VF3_aaa(fb,2),VF3_aaa(fb,2)) - e3*edgeF3(i,2)/lenF3(i,2) - e3*edgeF3(i,1)/lenF3(i,1);
            jac(VF3_aaa(fb,2),VF3_aaa(fc,2)) = jac(VF3_aaa(fb,2),VF3_aaa(fc,2)) + e3*edgeF3(i,2)/lenF3(i,2);
            jac(VF3_aaa(fb,2),VF3_aaa(fa,2)) = jac(VF3_aaa(fb,2),VF3_aaa(fa,2)) + e3*edgeF3(i,1)/lenF3(i,1);
        elseif j==3
            jac(VF3_aaa(fc,2),VF3_aaa(fc,2)) = jac(VF3_aaa(fc,2),VF3_aaa(fc,2)) - e3*edgeF3(i,3)/lenF3(i,3) - e3*edgeF3(i,2)/lenF3(i,2);
            jac(VF3_aaa(fc,2),VF3_aaa(fa,2)) = jac(VF3_aaa(fc,2),VF3_aaa(fa,2)) + e3*edgeF3(i,3)/lenF3(i,3);
            jac(VF3_aaa(fc,2),VF3_aaa(fb,2)) = jac(VF3_aaa(fc,2),VF3_aaa(fb,2)) + e3*edgeF3(i,2)/lenF3(i,2);
        end
    end
end

for i=1:Ncrow
    for j=1:Nccol
        fcon=find(contact(i,j)==VFT);
        jac(fcon,:) = 0;
        jac(fcon,fcon) = 1;
    end
end

for ii=1:is1row
    minus1 = find(is1(ii)==VF1_aaa(:,1));
    plus1 = find(is1(ii)==VF2_aaa(:,1));
    jac(VF2_aaa(plus1,2),:) = jac(VF2_aaa(plus1,2),:) + jac(VF1_aaa(minus1,2),:);
    jac(VF1_aaa(minus1,2),:)=0;
    jac(VF1_aaa(minus1,2),VF1_aaa(minus1,2))=-1;
    jac(VF1_aaa(minus1,2),VF2_aaa(plus1,2))=1; 
end

for ii=1:is2row
    minus1 = find(is2(ii)==VF3_aaa(:,1));
    plus1 = find(is2(ii)==VF2_aaa(:,1));
    jac(VF2_aaa(plus1,2),:) = jac(VF2_aaa(plus1,2),:) + jac(VF3_aaa(minus1,2),:);
    jac(VF3_aaa(minus1,2),:)=0;
    jac(VF3_aaa(minus1,2),VF3_aaa(minus1,2))=-1;
    jac(VF3_aaa(minus1,2),VF2_aaa(plus1,2))=1; 
end

con1 = 5; % top contact
con2 = 5; % bottom contact

res = zeros(VFTrow,1);
for i=1:Ncrow
    for j=1:Nccol
        fcon=find(contact(i,j)==VFT);
        if i<=con1
            res(fcon) = 1;
        elseif i>con1 &&  i<=con1+con2
            res(fcon) = 0;
        end
    end
end

phi = jac\res;
phi_sorted = zeros(VFrow,1);
for i=1:size(VF1_aaa,1)
    phi_sorted(VF1_aaa(i,1))=phi(VF1_aaa(i,2));
end
for i=1:size(VF2_aaa,1)
    phi_sorted(VF2_aaa(i,1))=phi(VF2_aaa(i,2));
end
for i=1:size(VF3_aaa,1)
    phi_sorted(VF3_aaa(i,1))=phi(VF3_aaa(i,2));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
aaa = zeros(VF1row+VF2row+VF3row,1); %%% variable for total regions
for i=1:size(aaa,1)
    if i<=VF1row
        aaa(i)=phi_sorted(VF1(i));
    elseif i>VF1row && i<=VF1row+VF2row
        aaa(i)=phi_sorted(VF2(i-VF1row));
    else
        aaa(i)=phi_sorted(VF3(i-(VF1row+VF2row)));
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%%% visualize part
figure
patch('Faces',F1,'Vertices',V, 'FaceVertexCData',phi_sorted, 'EdgeColor','blue','FaceColor','interp','LineWidth',1, 'Marker','o');
patch('Faces',F2,'Vertices',V, 'FaceVertexCData',phi_sorted, 'EdgeColor','green','FaceColor','interp','LineWidth',1, 'Marker','o');
patch('Faces',F3,'Vertices',V, 'FaceVertexCData',phi_sorted, 'EdgeColor','blue','FaceColor','interp','LineWidth',1, 'Marker','o');
title('Φ visualizing structure')
patch('Faces',interedge1,'Vertices',V, 'EdgeColor','magenta','FaceColor','none','LineWidth', 3, 'Marker','o');
patch('Faces',interedge2,'Vertices',V, 'EdgeColor','magenta','FaceColor','none','LineWidth', 3, 'Marker','o');
patch('Faces',interedge3,'Vertices',V, 'EdgeColor','magenta','FaceColor','none','LineWidth', 3, 'Marker','o');
colorbar
