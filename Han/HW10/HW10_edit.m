clc
clear
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% HW_edit
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Define
q = 1.602e-19;
e0 = 8.8542e-12;
nint = 1e16;
Ndop = -1e24;
Kb = 1.38e-23;
T = 300;
Vt = Kb*T/q;
con1 = 0.33374; % number of contact1
con2 = 0.33374; % number of contact2
Niter = 10; % number of iteration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
V = importdata("vertex.txt"); Vrow = size(V,1);
E = importdata("element.txt"); Erow = size(E,1);
E1 = importdata("element_region1.txt"); E1row = size(E1,1);
E2 = importdata("element_region2.txt"); E2row = size(E2,1);
E3 = importdata("element_region3.txt"); E3row = size(E3,1);

contact = importdata("contact.txt"); [conrow, concol] = size(contact);
unicon = unique(contact);
Ncon1 = 5; % number of top contact edge
Ncon2 = 5; % number bottom contact edge

VE = unique(E); VErow= size(VE,1);
VE1 = unique(E1); VE1row = size(VE1, 1);
VE2 = unique(E2); VE2row = size(VE2, 1);
VE3 = unique(E3); VE3row = size(VE3, 1);

Nregion = 3;

Nvertex = zeros(Nregion,1);
Nvertex(1,1) = VE1row;
Nvertex(2,1) = VE2row;
Nvertex(3,1) = VE3row;
sumvertex = zeros(Nregion+1,1);
sumvertex(1,1) = 0;

for i=1:Nregion
    sumvertex(i+1,1) = sumvertex(i,1) + Nvertex(i,1);
end

Nvariable = 3;

%%% initial potential indexing
X1 = zeros(sum(Nvertex(:,1)),Nregion);
X1(:,1) = vertcat(VE1,VE2,VE3);
count1 = 0;
index1 = 0;
for iregion = 1:Nregion
    for ivertex = 1:Nvertex(iregion,1)
        count1 = count1 + 1;
        index1 = index1 + 1;
        X1(count1,2) = index1;
    end
end

for i=1:Nvertex(1)
    X1(i,3) = X1(i,2);
end

for i=2:Nregion
    loc = sumvertex(i)+1;
    for j=sumvertex(i)+1:sumvertex(i+1)
        for k=sumvertex(i-1)+1:sumvertex(i)
            if X1(j,1)==X1(k,1)
                X1(k,i+2) = X1(j,2);
            else
                X1(loc,i+2) = X1(j,2);
            end
        end
        loc = loc+1;
    end
end

[C, ia, ic] = unique(X1(:,1));
X1_sorted = zeros(size(ia,1),Nregion);
for i=1:size(ia,1)
    X1_sorted(i,:) = X1(ia(i),3:Nregion+2);
end

%%% Total variable indexing
X2 = zeros(3*sum(Nvertex(:,1)),1);
count2 = 0;
index2 = 0;
for iregion = 1:Nregion
    for ivertex = 1:Nvertex(iregion,1)
        for ivariable = 1:3
            if iregion == 2
                count2 = count2 + 1;
                index2 = index2 + 1;
                X2(count2,1) = index2;
            else
                if ivariable == 1
                    count2 = count2 + 1;
                    index2 = index2 + 1;
                    X2(count2,1) = index2;
                else
                    count2 = count2 + 1;
                    X2(count2,1) = 0;
                end
            end
        end
    end
end

%%% vertex와 vertex 사이의 length
len = zeros(Erow,3);
for i=1:Erow
    for j =1:3
        if j==3
            len(i,j) = norm(V(E(i,j),:)-V(E(i,1),:));
        else
            len(i,j) = norm(V(E(i,j),:)-V(E(i,j+1),:));
        end
    end
end

%%% vertex 1->2, 1->3을 연결하는 vector와 vector의 norm 값
v21 = zeros(Erow,2);
v31 = zeros(Erow,2);
norm_v21 = zeros(Erow,1);
norm_v31 = zeros(Erow,1);
for i=1:Erow
    v21(i,:) = V(E(i,2),:)-V(E(i,1),:);
    v31(i,:) = V(E(i,3),:)-V(E(i,1),:);
    norm_v21(i,1) = norm(v21(i,:));
    norm_v31(i,1) = norm(v31(i,:));
end

%%% element로 이루어진 face의 area, face의 외접원 radius
area = zeros(Erow,1);
rad = zeros(Erow,1);
for i=1:Erow
    area(i,1) = sqrt(norm_v21(i,1)*norm_v21(i,1)*norm_v31(i,1)*norm_v31(i,1) ...
        - dot(v21(i,:),v31(i,:))*dot(v21(i,:),v31(i,:)))/2;
    rad(i,1) = (len(i,1)*len(i,2)*len(i,3)) / (4*area(i,1));
end

%%% element edge / edge1(1-2), edge2(2-3), edge3(3-1)
edge = zeros(Erow,1);
for i=1:Erow
    for j=1:3
        edge(i,j) = real(sqrt(rad(i,1)*rad(i,1) - (len(i,j)/2)*(len(i,j)/2)));
    end
end

%%% interface edge
is1=intersect(E1, E2);
is1row = size(is1,1);
interedge1 = zeros(is1row-1,2);
j=1;
for i = 1:E1row
    if size(intersect(E1(i,:),is1'),2) == 2
        interedge1(j,:) = intersect(E1(i,:),is1');
        j=j+1;
    end
end

is2=intersect(E2, E3);
is2row = size(is2,1);
interedge2 = zeros(is2row-1,2);
j=1;
for i = 1:E2row
    if size(intersect(E2(i,:),is2'),2) == 2
        interedge2(j,:) = intersect(E2(i,:),is2');
        j=j+1;
    end
end

is3=intersect(E3, E1);
is3row = size(is3, 1);
interedge3 = zeros(is3row-1,2);
j=1;
for i = 1:E3row
    if size(intersect(E3(i,:),is3'),2) == 2
        interedge3(j,:) = intersect(E3(i,:),is3');
        j=j+1;
    end
end

%%% initial value calculation
e1 = 11.7; e2 = 3.9;  e3 = 11.7;
edgeE1=zeros(E1row,3); lenE1=zeros(E1row,3);
edgeE2=zeros(E2row,3); lenE2=zeros(E2row,3);
edgeE3=zeros(E3row,3); lenE3=zeros(E3row,3);

aa=1; bb=1; cc=1;
for i=1:Erow
    for j=1:E1row
        if E(i,:) == E1(j,:)
            edgeE1(aa,:) = edge(i,:);
            lenE1(aa,:) = len(i,:);
            aa=aa+1;
        end
    end
    for j=1:E2row
        if E(i,:) == E2(j,:)
            edgeE2(bb,:) = edge(i,:);
            lenE2(bb,:) = len(i,:);
            bb=bb+1;
        end
    end
    for j=1:E3row
        if E(i,:) == E3(j,:)
            edgeE3(cc,:) = edge(i,:);
            lenE3(cc,:) = len(i,:);
            cc=cc+1;
        end
    end
end

jac1 = zeros(sum(Nvertex(:,1)),sum(Nvertex(:,1)));
jac1row = size(jac1,1);
res1 = zeros(sum(Nvertex(:,1)),1);
for i=1:E1row
    jac1(X1_sorted(E1(i,1),1),X1_sorted(E1(i,1),1)) = jac1(X1_sorted(E1(i,1),1),X1_sorted(E1(i,1),1))...
        - e1*(edgeE1(i,1)/lenE1(i,1) + edgeE1(i,3)/lenE1(i,3));
    jac1(X1_sorted(E1(i,1),1),X1_sorted(E1(i,2),1)) = jac1(X1_sorted(E1(i,1),1),X1_sorted(E1(i,2),1)) + e1*edgeE1(i,1)/lenE1(i,1);
    jac1(X1_sorted(E1(i,1),1),X1_sorted(E1(i,3),1)) = jac1(X1_sorted(E1(i,1),1),X1_sorted(E1(i,3),1)) + e1*edgeE1(i,3)/lenE1(i,3);

    jac1(X1_sorted(E1(i,2),1),X1_sorted(E1(i,2),1)) = jac1(X1_sorted(E1(i,2),1),X1_sorted(E1(i,2),1))...
        - e1*(edgeE1(i,1)/lenE1(i,1) + edgeE1(i,2)/lenE1(i,2));
    jac1(X1_sorted(E1(i,2),1),X1_sorted(E1(i,1),1)) = jac1(X1_sorted(E1(i,2),1),X1_sorted(E1(i,1),1)) + e1*edgeE1(i,1)/lenE1(i,1);
    jac1(X1_sorted(E1(i,2),1),X1_sorted(E1(i,3),1)) = jac1(X1_sorted(E1(i,2),1),X1_sorted(E1(i,3),1)) + e1*edgeE1(i,2)/lenE1(i,2);

    jac1(X1_sorted(E1(i,3),1),X1_sorted(E1(i,3),1)) = jac1(X1_sorted(E1(i,3),1),X1_sorted(E1(i,3),1))...
        - e1*(edgeE1(i,2)/lenE1(i,2) + edgeE1(i,3)/lenE1(i,3));
    jac1(X1_sorted(E1(i,3),1),X1_sorted(E1(i,2),1)) = jac1(X1_sorted(E1(i,3),1),X1_sorted(E1(i,2),1)) + e1*edgeE1(i,2)/lenE1(i,2);
    jac1(X1_sorted(E1(i,3),1),X1_sorted(E1(i,1),1)) = jac1(X1_sorted(E1(i,3),1),X1_sorted(E1(i,1),1)) + e1*edgeE1(i,3)/lenE1(i,3);
end

for i=1:E2row
    jac1(X1_sorted(E2(i,1),2),X1_sorted(E2(i,1),2)) = jac1(X1_sorted(E2(i,1),2),X1_sorted(E2(i,1),2))...
        - e2*(edgeE2(i,1)/lenE2(i,1) + edgeE2(i,3)/lenE2(i,3));
    jac1(X1_sorted(E2(i,1),2),X1_sorted(E2(i,2),2)) = jac1(X1_sorted(E2(i,1),2),X1_sorted(E2(i,2),2)) + e2*edgeE2(i,1)/lenE2(i,1);
    jac1(X1_sorted(E2(i,1),2),X1_sorted(E2(i,3),2)) = jac1(X1_sorted(E2(i,1),2),X1_sorted(E2(i,3),2)) + e2*edgeE2(i,3)/lenE2(i,3);
    res1(X1_sorted(E2(i,1),2),1) = res1(X1_sorted(E2(i,1),2),1) + (edgeE2(i,1)*lenE2(i,1) + edgeE2(i,3)*lenE2(i,3))/4*q/e0*(-Ndop);

    jac1(X1_sorted(E2(i,2),2),X1_sorted(E2(i,2),2)) = jac1(X1_sorted(E2(i,2),2),X1_sorted(E2(i,2),2))...
        - e2*(edgeE2(i,1)/lenE2(i,1) + edgeE2(i,2)/lenE2(i,2));
    jac1(X1_sorted(E2(i,2),2),X1_sorted(E2(i,1),2)) = jac1(X1_sorted(E2(i,2),2),X1_sorted(E2(i,1),2)) + e2*edgeE2(i,1)/lenE2(i,1);
    jac1(X1_sorted(E2(i,2),2),X1_sorted(E2(i,3),2)) = jac1(X1_sorted(E2(i,2),2),X1_sorted(E2(i,3),2)) + e2*edgeE2(i,2)/lenE2(i,2);
    res1(X1_sorted(E2(i,2),2),1) = res1(X1_sorted(E2(i,2),2),1) + (edgeE2(i,1)*lenE2(i,1) + edgeE2(i,2)*lenE2(i,2))/4*q/e0*(-Ndop);

    jac1(X1_sorted(E2(i,3),2),X1_sorted(E2(i,3),2)) = jac1(X1_sorted(E2(i,3),2),X1_sorted(E2(i,3),2))...
        - e2*(edgeE2(i,2)/lenE2(i,2) + edgeE2(i,3)/lenE2(i,3));
    jac1(X1_sorted(E2(i,3),2),X1_sorted(E2(i,2),2)) = jac1(X1_sorted(E2(i,3),2),X1_sorted(E2(i,2),2)) + e2*edgeE2(i,2)/lenE2(i,2);
    jac1(X1_sorted(E2(i,3),2),X1_sorted(E2(i,1),2)) = jac1(X1_sorted(E2(i,3),2),X1_sorted(E2(i,1),2)) + e2*edgeE2(i,3)/lenE2(i,3);
    res1(X1_sorted(E2(i,3),2),1) = res1(X1_sorted(E2(i,3),2),1) + (edgeE2(i,2)*lenE2(i,2) + edgeE2(i,3)*lenE2(i,3))/4*q/e0*(-Ndop);
end

for i=1:E3row
    jac1(X1_sorted(E3(i,1),3),X1_sorted(E3(i,1),3)) = jac1(X1_sorted(E3(i,1),3),X1_sorted(E3(i,1),3))...
        - e3*(edgeE3(i,1)/lenE3(i,1) + edgeE3(i,3)/lenE3(i,3));
    jac1(X1_sorted(E3(i,1),3),X1_sorted(E3(i,2),3)) = jac1(X1_sorted(E3(i,1),3),X1_sorted(E3(i,2),3)) + e3*edgeE3(i,1)/lenE3(i,1);
    jac1(X1_sorted(E3(i,1),3),X1_sorted(E3(i,3),3)) = jac1(X1_sorted(E3(i,1),3),X1_sorted(E3(i,3),3)) + e3*edgeE3(i,3)/lenE3(i,3);

    jac1(X1_sorted(E3(i,2),3),X1_sorted(E3(i,2),3)) = jac1(X1_sorted(E3(i,2),3),X1_sorted(E3(i,2),3))...
        - e3*(edgeE3(i,1)/lenE3(i,1) + edgeE3(i,2)/lenE3(i,2));
    jac1(X1_sorted(E3(i,2),3),X1_sorted(E3(i,1),3)) = jac1(X1_sorted(E3(i,2),3),X1_sorted(E3(i,1),3)) + e3*edgeE3(i,1)/lenE3(i,1);
    jac1(X1_sorted(E3(i,2),3),X1_sorted(E3(i,3),3)) = jac1(X1_sorted(E3(i,2),3),X1_sorted(E3(i,3),3)) + e3*edgeE3(i,2)/lenE3(i,2);

    jac1(X1_sorted(E3(i,3),3),X1_sorted(E3(i,3),3)) = jac1(X1_sorted(E3(i,3),3),X1_sorted(E3(i,3),3))...
        - e3*(edgeE3(i,2)/lenE3(i,2) + edgeE3(i,3)/lenE3(i,3));
    jac1(X1_sorted(E3(i,3),3),X1_sorted(E3(i,2),3)) = jac1(X1_sorted(E3(i,3),3),X1_sorted(E3(i,2),3)) + e3*edgeE3(i,2)/lenE3(i,2);
    jac1(X1_sorted(E3(i,3),3),X1_sorted(E3(i,1),3)) = jac1(X1_sorted(E3(i,3),3),X1_sorted(E3(i,1),3)) + e3*edgeE3(i,3)/lenE3(i,3);
end

for i=1:conrow
    for j=1:concol
        fcon=find(contact(i,j)==X1(:,1));
        jac1(X1(fcon,2),:) = 0;
        jac1(X1(fcon,2),X1(fcon,2)) = 1;
    end
end

for i=1:is1row
    minus1 = X1_sorted(is1(i),1);
    plus1 = X1_sorted(is1(i),2);
    jac1(X1(plus1,2),:) = jac1(X1(plus1,2),:) + jac1(X1(minus1,2),:);
    jac1(X1(minus1,2),:) = 0;
    jac1(X1(minus1,2),X1(minus1,2))=-1;
    jac1(X1(minus1,2),X1(plus1,2))=1;
end

for i=1:is2row
    minus1 = X1_sorted(is2(i),3);
    plus1 = X1_sorted(is2(i),2);
    jac1(X1(plus1,2),:) = jac1(X1(plus1,2),:) + jac1(X1(minus1,2),:);
    jac1(X1(minus1,2),:) = 0;
    jac1(X1(minus1,2), X1(minus1,2))= -1;
    jac1(X1(minus1,2), X1(plus1,2))= 1;
end

for i=1:conrow
    for j=1:concol
        fcon=find(contact(i,j)==X1(:,1));
        if i<=Ncon1
            res1(X1(fcon,2)) = con1;
        elseif i>Ncon1 &&  i<=Ncon1+Ncon2
            res1(X1(fcon,2)) = con2;
        end
    end
end

%%% initial potential
phi = jac1\res1;

%%% initial electron
n = zeros(VE2row,1);
for i=1:VE2row
    n(i) = nint*exp(phi(X1_sorted(VE2(i),2))/Vt);
end

%%% initial hole
p = zeros(VE2row,1);
for i=1:VE2row
    p(i) = nint*exp(-phi(X1_sorted(VE2(i),2))/Vt);
end

%%% jacobian, residue matrix
X2_sorted = find(X2);
for it=1:Niter
    jac2 = zeros(size(X2_sorted,1),size(X2_sorted,1));
    jac2row = size(size(X2_sorted,1),1);
    res2 = zeros(size(X2_sorted,1),1);
    save_phi_update = zeros(Niter,1);
    
    for i=1:E1row
        jac2(X2(3*(X1_sorted(E1(i,1),1))-2,1),X2(3*(X1_sorted(E1(i,1),1))-2,1)) = jac2(X2(3*(X1_sorted(E1(i,1),1))-2,1),X2(3*(X1_sorted(E1(i,1),1))-2,1))...
            - e1*(edgeE1(i,1)/lenE1(i,1) + edgeE1(i,3)/lenE1(i,3));
        jac2(X2(3*(X1_sorted(E1(i,1),1))-2,1),X2(3*(X1_sorted(E1(i,2),1))-2,1)) = jac2(X2(3*(X1_sorted(E1(i,1),1))-2,1),X2(3*(X1_sorted(E1(i,2),1))-2,1))...
            + e1*edgeE1(i,1)/lenE1(i,1);
        jac2(X2(3*(X1_sorted(E1(i,1),1))-2,1),X2(3*(X1_sorted(E1(i,3),1))-2,1)) = jac2(X2(3*(X1_sorted(E1(i,1),1))-2,1),X2(3*(X1_sorted(E1(i,3),1))-2,1))...
            + e1*edgeE1(i,3)/lenE1(i,3);
        res2(X2(3*(X1_sorted(E1(i,1),1))-2,1),1) = res2(X2(3*(X1_sorted(E1(i,1),1))-2,1),1)...
            - e1*(edgeE1(i,1)/lenE1(i,1) + edgeE1(i,3)/lenE1(i,3))*phi(X1_sorted(E1(i,1),1))...
            + e1*edgeE1(i,1)/lenE1(i,1)*phi(X1_sorted(E1(i,2),1))...
            + e1*edgeE1(i,3)/lenE1(i,3)*phi(X1_sorted(E1(i,3),1));

        jac2(X2(3*(X1_sorted(E1(i,2),1))-2,1),X2(3*(X1_sorted(E1(i,2),1))-2,1)) = jac2(X2(3*(X1_sorted(E1(i,2),1))-2,1),X2(3*(X1_sorted(E1(i,2),1))-2,1))...
            - e1*(edgeE1(i,1)/lenE1(i,1) + edgeE1(i,2)/lenE1(i,2));
        jac2(X2(3*(X1_sorted(E1(i,2),1))-2,1),X2(3*(X1_sorted(E1(i,1),1))-2,1)) = jac2(X2(3*(X1_sorted(E1(i,2),1))-2,1),X2(3*(X1_sorted(E1(i,1),1))-2,1))...
            + e1*edgeE1(i,1)/lenE1(i,1);
        jac2(X2(3*(X1_sorted(E1(i,2),1))-2,1),X2(3*(X1_sorted(E1(i,3),1))-2,1)) = jac2(X2(3*(X1_sorted(E1(i,2),1))-2,1),X2(3*(X1_sorted(E1(i,3),1))-2,1))...
            + e1*edgeE1(i,2)/lenE1(i,2);
        res2(X2(3*(X1_sorted(E1(i,2),1))-2,1),1) = res2(X2(3*(X1_sorted(E1(i,2),1))-2,1),1)...
            - e1*(edgeE1(i,1)/lenE1(i,1) + edgeE1(i,2)/lenE1(i,2))*phi(X1_sorted(E1(i,2),1))...
            + e1*edgeE1(i,1)/lenE1(i,1)*phi(X1_sorted(E1(i,1),1))...
            + e1*edgeE1(i,2)/lenE1(i,2)*phi(X1_sorted(E1(i,3),1));

        jac2(X2(3*(X1_sorted(E1(i,3),1))-2,1),X2(3*(X1_sorted(E1(i,3),1))-2,1)) = jac2(X2(3*(X1_sorted(E1(i,3),1))-2,1),X2(3*(X1_sorted(E1(i,3),1))-2,1))...
            - e1*(edgeE1(i,2)/lenE1(i,2) + edgeE1(i,3)/lenE1(i,3));
        jac2(X2(3*(X1_sorted(E1(i,3),1))-2,1),X2(3*(X1_sorted(E1(i,2),1))-2,1)) = jac2(X2(3*(X1_sorted(E1(i,3),1))-2,1),X2(3*(X1_sorted(E1(i,2),1))-2,1))...
            + e1*edgeE1(i,2)/lenE1(i,2);
        jac2(X2(3*(X1_sorted(E1(i,3),1))-2,1),X2(3*(X1_sorted(E1(i,1),1))-2,1)) = jac2(X2(3*(X1_sorted(E1(i,3),1))-2,1),X2(3*(X1_sorted(E1(i,1),1))-2,1))...
            + e1*edgeE1(i,3)/lenE1(i,3);
        res2(X2(3*(X1_sorted(E1(i,3),1))-2,1),1) = res2(X2(3*(X1_sorted(E1(i,3),1))-2,1),1)...
            - e1*(edgeE1(i,2)/lenE1(i,2) + edgeE1(i,3)/lenE1(i,3))*phi(X1_sorted(E1(i,3),1))...
            + e1*edgeE1(i,2)/lenE1(i,2)*phi(X1_sorted(E1(i,2),1))...
            + e1*edgeE1(i,3)/lenE1(i,3)*phi(X1_sorted(E1(i,1),1));
    end

    saveexl = zeros(size(X1,1),1);

    for i=1:E2row
        jac2(X2(3*(X1_sorted(E2(i,1),2))-2,1),X2(3*(X1_sorted(E2(i,1),2))-2,1)) = jac2(X2(3*(X1_sorted(E2(i,1),2))-2,1),X2(3*(X1_sorted(E2(i,1),2))-2,1))...
            - e2*(edgeE2(i,1)/lenE2(i,1) + edgeE2(i,3)/lenE2(i,3));
        jac2(X2(3*(X1_sorted(E2(i,1),2))-2,1),X2(3*(X1_sorted(E2(i,2),2))-2,1)) = jac2(X2(3*(X1_sorted(E2(i,1),2))-2,1),X2(3*(X1_sorted(E2(i,2),2))-2,1))...
            + e2*edgeE2(i,1)/lenE2(i,1);
        jac2(X2(3*(X1_sorted(E2(i,1),2))-2,1),X2(3*(X1_sorted(E2(i,3),2))-2,1)) = jac2(X2(3*(X1_sorted(E2(i,1),2))-2,1),X2(3*(X1_sorted(E2(i,3),2))-2,1))...
            + e2*edgeE2(i,3)/lenE2(i,3);
        res2(X2(3*(X1_sorted(E2(i,1),2))-2,1),1) =  res2(X2(3*(X1_sorted(E2(i,1),2))-2,1),1)...
            + (edgeE2(i,1)*lenE2(i,1) + edgeE2(i,3)*lenE2(i,3))/4*q/e0*Ndop...
            - e2*(edgeE2(i,1)/lenE2(i,1) + edgeE2(i,3)/lenE2(i,3))*phi(X1_sorted(E2(i,1),2))...
            + e2*edgeE2(i,1)/lenE2(i,1)*phi(X1_sorted(E2(i,2),2))...
            + e2*edgeE2(i,3)/lenE2(i,3)*phi(X1_sorted(E2(i,3),2));
        saveexl(X1_sorted(E2(i,1),2)) = saveexl(X1_sorted(E2(i,1),2)) + (edgeE2(i,1)*lenE2(i,1) + edgeE2(i,3)*lenE2(i,3))/4;

        jac2(X2(3*(X1_sorted(E2(i,2),2))-2,1),X2(3*(X1_sorted(E2(i,2),2))-2,1)) = jac2(X2(3*(X1_sorted(E2(i,2),2))-2,1),X2(3*(X1_sorted(E2(i,2),2))-2,1))...
            - e2*(edgeE2(i,1)/lenE2(i,1) + edgeE2(i,2)/lenE2(i,2));
        jac2(X2(3*(X1_sorted(E2(i,2),2))-2,1),X2(3*(X1_sorted(E2(i,1),2))-2,1)) = jac2(X2(3*(X1_sorted(E2(i,2),2))-2,1),X2(3*(X1_sorted(E2(i,1),2))-2,1))...
            + e2*edgeE2(i,1)/lenE2(i,1);
        jac2(X2(3*(X1_sorted(E2(i,2),2))-2,1),X2(3*(X1_sorted(E2(i,3),2))-2,1)) = jac2(X2(3*(X1_sorted(E2(i,2),2))-2,1),X2(3*(X1_sorted(E2(i,3),2))-2,1))...
            + e2*edgeE2(i,2)/lenE2(i,2);
        res2(X2(3*(X1_sorted(E2(i,2),2))-2,1),1) =  res2(X2(3*(X1_sorted(E2(i,2),2))-2,1),1)...
            + (edgeE2(i,1)*lenE2(i,1) + edgeE2(i,2)*lenE2(i,2))/4*q/e0*Ndop...
            - e2*(edgeE2(i,1)/lenE2(i,1) + edgeE2(i,2)/lenE2(i,2))*phi(X1_sorted(E2(i,2),2))...
            + e2*edgeE2(i,1)/lenE2(i,1)*phi(X1_sorted(E2(i,1),2))...
            + e2*edgeE2(i,2)/lenE2(i,2)*phi(X1_sorted(E2(i,3),2));
        saveexl(X1_sorted(E2(i,2),2)) = saveexl(X1_sorted(E2(i,2),2)) + (edgeE2(i,1)*lenE2(i,1) + edgeE2(i,2)*lenE2(i,2))/4;

        jac2(X2(3*(X1_sorted(E2(i,3),2))-2,1),X2(3*(X1_sorted(E2(i,3),2))-2,1)) = jac2(X2(3*(X1_sorted(E2(i,3),2))-2,1),X2(3*(X1_sorted(E2(i,3),2))-2,1))...
            - e2*(edgeE2(i,2)/lenE2(i,2) + edgeE2(i,3)/lenE2(i,3));
        jac2(X2(3*(X1_sorted(E2(i,3),2))-2,1),X2(3*(X1_sorted(E2(i,2),2))-2,1)) = jac2(X2(3*(X1_sorted(E2(i,3),2))-2,1),X2(3*(X1_sorted(E2(i,2),2))-2,1))...
            + e2*edgeE2(i,2)/lenE2(i,2);
        jac2(X2(3*(X1_sorted(E2(i,3),2))-2,1),X2(3*(X1_sorted(E2(i,1),2))-2,1)) = jac2(X2(3*(X1_sorted(E2(i,3),2))-2,1),X2(3*(X1_sorted(E2(i,1),2))-2,1))...
            + e2*edgeE2(i,3)/lenE2(i,3);
        res2(X2(3*(X1_sorted(E2(i,3),2))-2,1),1) =  res2(X2(3*(X1_sorted(E2(i,3),2))-2,1),1)...
            + (edgeE2(i,2)*lenE2(i,2) + edgeE2(i,3)*lenE2(i,3))/4*q/e0*Ndop...
            - e2*(edgeE2(i,2)/lenE2(i,2) + edgeE2(i,3)/lenE2(i,3))*phi(X1_sorted(E2(i,3),2))...
            + e2*edgeE2(i,2)/lenE2(i,2)*phi(X1_sorted(E2(i,2),2))...
            + e2*edgeE2(i,3)/lenE2(i,3)*phi(X1_sorted(E2(i,1),2));
        saveexl(X1_sorted(E2(i,3),2)) = saveexl(X1_sorted(E2(i,3),2)) + (edgeE2(i,2)*lenE2(i,2) + edgeE2(i,3)*lenE2(i,3))/4;
    end

    for i=1:E3row
        jac2(X2(3*(X1_sorted(E3(i,1),3))-2,1),X2(3*(X1_sorted(E3(i,1),3))-2,1)) = jac2(X2(3*(X1_sorted(E3(i,1),3))-2,1),X2(3*(X1_sorted(E3(i,1),3))-2,1))...
            - e3*(edgeE3(i,1)/lenE3(i,1) + edgeE3(i,3)/lenE3(i,3));
        jac2(X2(3*(X1_sorted(E3(i,1),3))-2,1),X2(3*(X1_sorted(E3(i,2),3))-2,1)) = jac2(X2(3*(X1_sorted(E3(i,1),3))-2,1),X2(3*(X1_sorted(E3(i,2),3))-2,1))...
            + e3*edgeE3(i,1)/lenE3(i,1);
        jac2(X2(3*(X1_sorted(E3(i,1),3))-2,1),X2(3*(X1_sorted(E3(i,3),3))-2,1)) = jac2(X2(3*(X1_sorted(E3(i,1),3))-2,1),X2(3*(X1_sorted(E3(i,3),3))-2,1))...
            + e3*edgeE3(i,3)/lenE3(i,3);
        res2(X2(3*(X1_sorted(E3(i,1),3))-2,1),1) = res2(X2(3*(X1_sorted(E3(i,1),3))-2,1),1)...
            - e1*(edgeE3(i,1)/lenE3(i,1) + edgeE3(i,3)/lenE3(i,3))*phi(X1_sorted(E3(i,1),3))...
            + e1*edgeE3(i,1)/lenE3(i,1)*phi(X1_sorted(E3(i,2),3))...
            + e1*edgeE3(i,3)/lenE3(i,3)*phi(X1_sorted(E3(i,3),3));

        jac2(X2(3*(X1_sorted(E3(i,2),3))-2,1),X2(3*(X1_sorted(E3(i,2),3))-2,1)) = jac2(X2(3*(X1_sorted(E3(i,2),3))-2,1),X2(3*(X1_sorted(E3(i,2),3))-2,1))...
            - e3*(edgeE3(i,1)/lenE3(i,1) + edgeE3(i,2)/lenE3(i,2));
        jac2(X2(3*(X1_sorted(E3(i,2),3))-2,1),X2(3*(X1_sorted(E3(i,1),3))-2,1)) = jac2(X2(3*(X1_sorted(E3(i,2),3))-2,1),X2(3*(X1_sorted(E3(i,1),3))-2,1))...
            + e3*edgeE3(i,1)/lenE3(i,1);
        jac2(X2(3*(X1_sorted(E3(i,2),3))-2,1),X2(3*(X1_sorted(E3(i,3),3))-2,1)) = jac2(X2(3*(X1_sorted(E3(i,2),3))-2,1),X2(3*(X1_sorted(E3(i,3),3))-2,1))...
            + e3*edgeE3(i,2)/lenE3(i,2);
        res2(X2(3*(X1_sorted(E3(i,2),3))-2,1),1) = res2(X2(3*(X1_sorted(E3(i,2),3))-2,1),1)...
            - e1*(edgeE3(i,1)/lenE3(i,1) + edgeE3(i,2)/lenE3(i,2))*phi(X1_sorted(E3(i,2),3))...
            + e1*edgeE3(i,1)/lenE3(i,1)*phi(X1_sorted(E3(i,1),3))...
            + e1*edgeE3(i,2)/lenE3(i,2)*phi(X1_sorted(E3(i,3),3));

        jac2(X2(3*(X1_sorted(E3(i,3),3))-2,1),X2(3*(X1_sorted(E3(i,3),3))-2,1)) = jac2(X2(3*(X1_sorted(E3(i,3),3))-2,1),X2(3*(X1_sorted(E3(i,3),3))-2,1))...
            - e3*(edgeE3(i,2)/lenE3(i,2) + edgeE3(i,3)/lenE3(i,3));
        jac2(X2(3*(X1_sorted(E3(i,3),3))-2,1),X2(3*(X1_sorted(E3(i,2),3))-2,1)) = jac2(X2(3*(X1_sorted(E3(i,3),3))-2,1),X2(3*(X1_sorted(E3(i,2),3))-2,1))...
            + e3*edgeE3(i,2)/lenE3(i,2);
        jac2(X2(3*(X1_sorted(E3(i,3),3))-2,1),X2(3*(X1_sorted(E3(i,1),3))-2,1)) = jac2(X2(3*(X1_sorted(E3(i,3),3))-2,1),X2(3*(X1_sorted(E3(i,1),3))-2,1))...
            + e3*edgeE3(i,3)/lenE3(i,3);
        res2(X2(3*(X1_sorted(E3(i,3),3))-2,1),1) = res2(X2(3*(X1_sorted(E3(i,3),3))-2,1),1)...
            - e1*(edgeE3(i,2)/lenE3(i,2) + edgeE3(i,3)/lenE3(i,3))*phi(X1_sorted(E3(i,3),3))...
            + e1*edgeE3(i,2)/lenE3(i,2)*phi(X1_sorted(E3(i,2),3))...
            + e1*edgeE3(i,3)/lenE3(i,3)*phi(X1_sorted(E3(i,1),3));
    end

    %%%  for electron
    for i=1:VE2row
        %%% n-nit*exp(phi/Vt)
        jac2(X2(3*X1_sorted(VE2(i),2)-1),X2(3*X1_sorted(VE2(i),2)-1)) = 1;
        jac2(X2(3*X1_sorted(VE2(i),2)-1),X2(3*X1_sorted(VE2(i),2)-2)) = - nint/Vt*exp(phi(X1_sorted(VE2(i),2))/Vt);
        res2(X2(3*X1_sorted(VE2(i),2)-1),1) = n(i) - nint*exp(phi(X1_sorted(VE2(i),2))/Vt);
        %%% - qn 
        jac2(X2(3*X1_sorted(VE2(i),2)-2),X2(3*X1_sorted(VE2(i),2)-1)) = - saveexl(X1_sorted(VE2(i),2))*q/e0;
        res2(X2(3*X1_sorted(VE2(i),2)-2),1) = res2(X2(3*X1_sorted(VE2(i),2)-2),1) - saveexl(X1_sorted(VE2(i),2))*q/e0*n(i);
    end

    %%%  for hole
    for i=1:VE2row
        %%% p-nit*exp(-phi/Vt)
        jac2(X2(3*X1_sorted(VE2(i),2)),X2(3*X1_sorted(VE2(i),2))) = 1;
        jac2(X2(3*X1_sorted(VE2(i),2)),X2(3*X1_sorted(VE2(i),2)-2)) = nint/Vt*exp(-phi(X1_sorted(VE2(i),2))/Vt);
        res2(X2(3*X1_sorted(VE2(i),2)),1) = p(i) - nint*exp(-phi(X1_sorted(VE2(i),2))/Vt);
        %%% + qp 
        jac2(X2(3*X1_sorted(VE2(i),2)-2),X2(3*X1_sorted(VE2(i),2))) = saveexl(X1_sorted(VE2(i),2))*q/e0;
        res2(X2(3*X1_sorted(VE2(i),2)-2),1) = res2(X2(3*X1_sorted(VE2(i),2)-2),1) + saveexl(X1_sorted(VE2(i),2))*q/e0*p(i);
    end

    for i=1:conrow
        for j=1:concol
            fcon=find(contact(i,j)==X1(:,1));
            jac2(X2(3*X1(fcon,2)-2),:) = 0;
            jac2(X2(3*X1(fcon,2)-2),X2(3*X1(fcon,2)-2)) = 1;
            if i<=Ncon1
                res2(X2(3*X1(fcon,2)-2),1) = phi(X1(fcon,2),1) - con1;
            elseif i>Ncon1 &&  i<=Ncon1+Ncon2
                res2(X2(3*X1(fcon,2)-2),1) = phi(X1(fcon,2),1) - con2;
            end
        end
    end

    for i=1:is1row
        minus1 = X1_sorted(is1(i),1);
        plus1 = X1_sorted(is1(i),2);
        jac2(X2(3*X1(plus1,2)-2),:) = jac2(X2(3*X1(plus1,2)-2),:) + jac2(X2(3*X1(minus1,2)-2),:);
        res2(X2(3*X1(plus1,2)-2)) = res2(X2(3*X1(plus1,2)-2)) + res2(X2(3*X1(minus1,2)-2));
        jac2(X2(3*X1(minus1,2)-2),:) = 0;
        res2(X2(3*X1(minus1,2)-2),1) = 0;
        jac2(X2(3*X1(minus1,2)-2), X2(3*X1(minus1,2)-2)) = -1;
        jac2(X2(3*X1(minus1,2)-2), X2(3*X1(plus1,2)-2)) = 1;
    end

    for i=1:is2row
        minus1 = X1_sorted(is2(i),3);
        plus1 = X1_sorted(is2(i),2);
        jac2(X2(3*X1(plus1,2)-2),:) = jac2(X2(3*X1(plus1,2)-2),:) + jac2(X2(3*X1(minus1,2)-2),:);
        res2(X2(3*X1(plus1,2)-2),1) = res2(X2(3*X1(plus1,2)-2),1) + res2(X2(3*X1(minus1,2)-2),1);
        jac2(X2(3*X1(minus1,2)-2),:) = 0;
        res2(X2(3*X1(minus1,2)-2),1) = 0;
        jac2(X2(3*X1(minus1,2)-2), X2(3*X1(minus1,2)-2)) = -1;
        jac2(X2(3*X1(minus1,2)-2), X2(3*X1(plus1,2)-2)) = 1;
    end

    %%% scaled
    Cvector = zeros(size(jac2,1),1);
    for i=1:size(jac2,1)
        if rem(X2_sorted(i),3) == 1
            Cvector(i) = Vt;
        elseif rem(X2_sorted(i),3) == 2
            Cvector(i) = Ndop;
        elseif rem(X2_sorted(i),3) == 0
            Cvector(i) = Ndop;
        end
    end

    Cmatrix = spdiags(Cvector,0,size(X2_sorted,1),size(X2_sorted,1));
    jac2_scaled = jac2 * Cmatrix;
    Rvector = 1./sum(abs(jac2_scaled),2);
    Rmatrix = spdiags(Rvector,0,size(X2_sorted,1),size(X2_sorted,1));
    jac2_scaled = Rmatrix * jac2_scaled;
    res2_scaled = Rmatrix * res2;
    update_scaled = jac2_scaled \ (-res2_scaled);
    update = Cmatrix * update_scaled;

    countphi=1;
    countn=1;
    countp=1;

    maxupdate = 0;
    for i=1:size(X2_sorted,1)
        if rem(X2_sorted(i),3) == 1
            phi(countphi) = phi(countphi) + update(i);
            if maxupdate < max(abs(update(i)))
                maxupdate = max(abs(update(i)));
            end
            countphi = countphi+1;
        elseif rem(X2_sorted(i),3) == 2
            n(countn) = n(countn) + update(i);
            countn = countn + 1;
        elseif rem(X2_sorted(i),3) == 0
            p(countp) = p(countp) + update(i);
            countp = countp + 1;
        end
    end

    save_phi_update(it,1) = maxupdate;

end

figure
semilogy(save_phi_update)
xlabel('iteration number')
ylabel('potential update [V]')
title('max potential update')

%%% variable
% pot = zeros(sum(Nvertex(:,1)),1); %%% potential total regions
% elec = zeros(Nvertex(2),1); %%% electron for semiconductor region
% hole = zeros(Nvertex(2),1); %%% hole for semiconductor region
