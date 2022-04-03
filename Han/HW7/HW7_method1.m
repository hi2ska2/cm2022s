clc
clear
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% HW7
%%% method 1: [aaa; bbb; ccc] 
%%% Seong-Min,Han
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
V = importdata("vertex.txt"); Vrow = size(V,1);
F = importdata("element.txt"); Frow = size(F,1);
F1 = importdata("element_region1.txt"); F1row = size(F1,1);
F2 = importdata("element_region2.txt"); F2row = size(F2,1);
F3 = importdata("element_region3.txt"); F3row = size(F3,1);
contact = importdata("contact.txt");
VF = unique(F); VFrow= size(VF,1);
VF1 = unique(F1); VF1row = size(VF1, 1);
VF2 = unique(F2); VF2row = size(VF2, 1);
VF3 = unique(F3); VF3row = size(VF3, 1);
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
%%% element로 이루어진 face의 area, face의 외접원 radius
area = zeros(Frow,1);
rad = zeros(Frow,1);
for i=1:Frow
    area(i,1) = sqrt(norm_v21(i,1)*norm_v21(i,1)*norm_v31(i,1)*norm_v31(i,1) ...
        - dot(v21(i,:),v31(i,:))*dot(v21(i,:),v31(i,:)))/2;
    rad(i,1) = (len(i,1)*len(i,2)*len(i,3)) / (4*area(i,1));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% face의 edge : edge1(1-2), edge2(2-3), edge3(3-1)
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

for i=1:Frow
    for j=1:3
        if checkcontact(i,j)==1
            A(F(i,j),:)= 0;
            A(F(i,j),F(i,j))=1;
        end
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% visualizing
figure
patch('Faces',F1,'Vertices',V, 'FaceVertexCData',phi, 'EdgeColor','blue','FaceColor','interp','LineWidth',1, 'Marker','o');
patch('Faces',F2,'Vertices',V, 'FaceVertexCData',phi, 'EdgeColor','green','FaceColor','interp','LineWidth',1, 'Marker','o');
patch('Faces',F3,'Vertices',V, 'FaceVertexCData',phi, 'EdgeColor','blue','FaceColor','interp','LineWidth',1, 'Marker','o');
title('Φ visualizing structure')
patch('Faces',interedge1,'Vertices',V, 'EdgeColor','magenta','FaceColor','none','LineWidth', 3, 'Marker','o');
patch('Faces',interedge2,'Vertices',V, 'EdgeColor','magenta','FaceColor','none','LineWidth', 3, 'Marker','o');
patch('Faces',interedge3,'Vertices',V, 'EdgeColor','magenta','FaceColor','none','LineWidth', 3, 'Marker','o');
colorbar
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% method 1
aaa = zeros(VFrow,1); %%% variable for total regionsaksdi
bbb = zeros(VF2row,1); %%% variable for semiconductor region
ccc = zeros(VF1row+VF3row,1); %%% variable for insulator regions

%%%  solution X matrix 
NX = VFrow+VF1row+VF2row+VF3row;
X = zeros(NX,1);
for i=1:NX
    if i <= VFrow
        X(i) = phi(i,1);
        aaa(i) = phi(i,1);
    elseif i > VFrow && i <= (VFrow+VF2row)
        X(i) = phi(VF2(i-VFrow));
        bbb(i-VFrow) = phi(VF2(i-VFrow));
    elseif i > (VFrow+VF2row) && i <= (VFrow+VF2row+VF1row)
        X(i) = phi(VF1(i-(VFrow+VF2row)));
        ccc(i-(VFrow+VF2row)) = phi(VF1(i-(VFrow+VF2row)));
    elseif i > (VFrow+VF2row+VF1row) && i <= (VFrow+VF2row+VF1row+VF3row)
        X(i) = phi(VF3(i-(VFrow+VF2row+VF1row)));
        ccc(i-(VFrow+VF2row+VF1row)) = phi(VF3(i-(VFrow+VF2row+VF1row)));
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% scan
while 1
    var = input('variable?\n (1)aaa(entire)\n (2)bbb(only semicinductor)\n (3)ccc(only insulator)\n : ');
    reg = input('region?\n (1)oxide\n (2)silicon\n : ');
    if var==2 && reg ~= 2 || var==3 && reg ~= 1
        fprintf('Error!\n');
        continue
    end
    ver = input('vertex number? : ');
%%% print 
    if var==1 &&    size(find(ver==VF),1)==1
        xindex = find(ver==VF);
        fprintf('X[%d] = %f\n',xindex, X(xindex));
    elseif var==2 && reg==2 && size(find(ver==VF2),1)==1
        xindex = VFrow+find(ver==VF2);
        fprintf('X[%d] = %f\n',xindex, X(xindex));
    elseif var==3
        if reg==1 && size(find(ver==VF1),1)==1
            xindex = VFrow+VF2row+find(ver==VF1);
            fprintf('X[%d] = %f\n',xindex, X(xindex));
        elseif reg==3 && size(find(ver==VF3),1)==1
            xindex = VFrow+VF2row+VF1row+find(ver==VF3);
            fprintf('X[%d] = %f\n',xindex, X(xindex));
        else
            fprintf('Error !\n');
        end
    else
        fprintf('Error !\n');
    end
    break
end
