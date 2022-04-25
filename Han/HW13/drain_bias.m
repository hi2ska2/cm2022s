clc
clear
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% HW
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Define
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
width = 1e-6;
dx = 5e-9;
dy = 2e-9;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
q = 1.6021e-19;
e0 = 8.8542e-12;
eox = 3.9;
esi = 11.7;
nint = 1e16;
Kb = 1.38e-23;
T = 300;
Vt = Kb*T/q;
mun = 1417e-4; % electron mobility coefficient (@equilibrium)
mup = 470e-4; % hole mobility coefficient (@equilibrium)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Doping density
Nd = 5e26;
Na = -2e21;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Iteration number
Niter1 = 20;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
V = importdata("vertex.txt"); Vrow = size(V,1);
E = importdata("element.txt"); Erow = size(E,1);
E1 = importdata("element_region1.txt"); E1row = size(E1,1);
E2 = importdata("element_region2.txt"); E2row = size(E2,1);
E3 = importdata("element_region3.txt"); E3row = size(E3,1);
E4 = importdata("element_region4.txt"); E4row = size(E4,1);
E5 = importdata("element_region5.txt"); E5row = size(E5,1);
Esi = vertcat(E2, E3, E4); Esirow = size(Esi,1);

% contact information
contact = importdata("contact.txt"); [conrow, concol] = size(contact);
unicon = unique(contact);
Ncon1 = 11;
Ncon2 = 11;
Ncon3 = 11;
Ncon4 = 11;

cond = zeros(2*Ncon3,1);
ccc=1;
for i=Ncon1+Ncon2+Ncon3+1:Ncon1+Ncon2+Ncon3+Ncon4
    for j=1:concol
        cond(ccc,1) = contact(i,j);
        ccc = ccc+1;
    end
end
unicond = unique(cond);

%number of variable
Nvariable = 3;

%number of region
Nregion = 5;

% number of vertex
VE = unique(E); VErow= size(VE,1);
VE1 = unique(E1); VE1row = size(VE1, 1);
VE2 = unique(E2); VE2row = size(VE2, 1);
VE3 = unique(E3); VE3row = size(VE3, 1);
VE4 = unique(E4); VE4row = size(VE4, 1);
VE5 = unique(E5); VE5row = size(VE5, 1);

Nvertex = zeros(Nregion,1);
Nvertex(1,1) = VE1row;
Nvertex(2,1) = VE2row;
Nvertex(3,1) = VE3row;
Nvertex(4,1) = VE4row;
Nvertex(5,1) = VE5row;

V_silicon(:,1) = unique(Esi);
V_silicon_row = size(V_silicon, 1);
V_oxide = zeros(VE1row+VE5row, 2);
V_oxide_row = size(V_oxide,1);

materialvertex = zeros(3,1);
materialvertex(1,1) = VE1row;
materialvertex(2,1) = V_silicon_row;
materialvertex(3,1) = VE5row;

sumvertex = zeros(4,1);
sumvertex(1,1) = 0;
for i=1:3
    if i==1
        sumvertex(i+1,1) = sumvertex(i,1) + VE1row;
    elseif i==2
        sumvertex(i+1,1) = sumvertex(i,1) + V_silicon_row;
    elseif i==3
        sumvertex(i+1,1) = sumvertex(i,1) + VE5row;
    end
end

%%% initial potential indexing
X1 = zeros(VE1row + V_silicon_row + VE5row ,2);
X1(:,1) = vertcat(VE1,V_silicon,VE5);
count1 = 0;
index1 = 0;
for iregion = 1:3
    for ivertex = 1:materialvertex(iregion,1)
        count1 = count1 + 1;
        index1 = index1 + 1;
        X1(count1,2) = index1;
    end
end

X1_sorted = zeros(VErow,1);

for i=1:3
    for j=sumvertex(i)+1:sumvertex(i+1)
        X1_sorted(X1(j,1),i) = X1(j,2);
    end
end

%%% Initial potential re-indexing
loc_ox = 1;
loc_si = 1;
for i=1:3
    if i==1 || i==3
        for j=sumvertex(i)+1:sumvertex(i+1)
            V_oxide(loc_ox,:) = X1(j,:);
            loc_ox = loc_ox+1;
        end
    else
        for j=sumvertex(i)+1:sumvertex(i+1)
            V_silicon(loc_si,2) = X1(j,2);
            loc_si = loc_si+1;
        end
    end
end

%%% Total variable indexing
X2 = zeros(3*(size(X1,1)),1);
count2 = 0;
index2 = 0;
for iregion = 1:3
    for ivertex = 1:materialvertex(iregion,1)
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
X2_sorted = find(X2);

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
is1=intersect(E1, Esi);
is1row = size(is1,1);
interedge1 = zeros(is1row-1,2);
j=1;
for i = 1:E1row
    if size(intersect(E1(i,:),is1'),2) == 2
        interedge1(j,:) = intersect(E1(i,:),is1');
        j=j+1;
    end
end

is2=intersect(E5, Esi);
is2row = size(is2,1);
interedge2 = zeros(is2row-1,2);
j=1;
for i = 1:E5row
    if size(intersect(E5(i,:),is2'),2) == 2
        interedge2(j,:) = intersect(E5(i,:),is2');
        j=j+1;
    end
end

%%% edge per region
edgeE1=zeros(E1row,3); lenE1=zeros(E1row,3);
edgeE2=zeros(E2row,3); lenE2=zeros(E2row,3);
edgeE3=zeros(E3row,3); lenE3=zeros(E3row,3);
edgeE4=zeros(E4row,3); lenE4=zeros(E4row,3);
edgeE5=zeros(E5row,3); lenE5=zeros(E5row,3);
edgeEsi=zeros(Esirow,3); lenEsi=zeros(Esirow,3);
aa=1; bb=1; cc=1; dd=1; ee=1;
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
    for j=1:E4row
        if E(i,:) == E4(j,:)
            edgeE4(dd,:) = edge(i,:);
            lenE4(dd,:) = len(i,:);
            dd=dd+1;
        end
    end
    for j=1:E5row
        if E(i,:) == E5(j,:)
            edgeE5(ee,:) = edge(i,:);
            lenE5(ee,:) = len(i,:);
            ee=ee+1;
        end
    end
end

edgeEsi = vertcat(edgeE2, edgeE3, edgeE4);
lenEsi = vertcat(lenE2, lenE3, lenE4);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
phi = zeros(size(X1,1),1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Initial potential_edit by euqation : Vt*ln(Ndop/nint)
for i=1:VE3row
    phi(X1_sorted(VE3(i),2)) = Vt*log(-Na/nint);
    if X1_sorted(VE3(i),1) ~= 0
        phi(X1_sorted(VE3(i),1)) = Vt*log(-Na/nint);
    elseif X1_sorted(VE3(i),3) ~= 0
        phi(X1_sorted(VE3(i),3)) = Vt*log(-Na/nint);
    end
end

for i=1:VE2row
    phi(X1_sorted(VE2(i),2)) = Vt*log(Nd/nint);
    if X1_sorted(VE2(i),1) ~= 0
        phi(X1_sorted(VE2(i),1)) = Vt*log(Nd/nint);
    elseif X1_sorted(VE2(i),3) ~= 0
        phi(X1_sorted(VE2(i),3)) = Vt*log(Nd/nint);
    end
end

for i=1:VE4row
    phi(X1_sorted(VE4(i),2)) = Vt*log(Nd/nint);
    if X1_sorted(VE4(i),1) ~= 0
        phi(X1_sorted(VE4(i),1)) = Vt*log(Nd/nint);
    elseif X1_sorted(VE4(i),3) ~= 0
        phi(X1_sorted(VE4(i),3)) = Vt*log(Nd/nint);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Dirichlet boundary condition
for i=1:conrow
    for j=1:concol
        fcon = find(contact(i,j)==X1(:,1));
        if i<=Ncon1
            phi(X1(fcon,2)) = 0.33374;
        elseif i>Ncon1 && i<=(Ncon1+Ncon2)
            phi(X1(fcon,2)) = 0.33374;
        elseif i>(Ncon1+Ncon2) && i<=(Ncon1+Ncon2+Ncon3)
            phi(X1(fcon,2)) = Vt*log(Nd/nint);
        elseif i>(Ncon1+Ncon2+Ncon3) && i<=(Ncon1+Ncon2+Ncon3+Ncon4)
            phi(X1(fcon,2)) = Vt*log(Nd/nint);
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% initial potential
initial_phi = phi(:,1);
%%% initial electron
n = zeros(V_silicon_row,1);
for i=1:V_silicon_row
    n(i) = nint*exp(phi(V_silicon(i,2))/Vt);
end

%%% initial hole
p = zeros(V_silicon_row,1);
for i=1:V_silicon_row
    p(i) = nint*exp(-phi(V_silicon(i,2))/Vt);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% jacobian, residue matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
save_phi_upd = zeros(Niter1,1);
save_upd = zeros(size(X2_sorted,1),1);

for it=1:Niter1

    jac1 = sparse(size(X2_sorted,1),size(X2_sorted,1));
    res1 = zeros(size(X2_sorted,1),1);

    for i=1:E1row
        jac1(X2(3*(X1_sorted(E1(i,1),1))-2,1),X2(3*(X1_sorted(E1(i,1),1))-2,1)) = jac1(X2(3*(X1_sorted(E1(i,1),1))-2,1),X2(3*(X1_sorted(E1(i,1),1))-2,1))...
            - eox*(edgeE1(i,1)/lenE1(i,1) + edgeE1(i,3)/lenE1(i,3));
        jac1(X2(3*(X1_sorted(E1(i,1),1))-2,1),X2(3*(X1_sorted(E1(i,2),1))-2,1)) = jac1(X2(3*(X1_sorted(E1(i,1),1))-2,1),X2(3*(X1_sorted(E1(i,2),1))-2,1))...
            + eox*edgeE1(i,1)/lenE1(i,1);
        jac1(X2(3*(X1_sorted(E1(i,1),1))-2,1),X2(3*(X1_sorted(E1(i,3),1))-2,1)) = jac1(X2(3*(X1_sorted(E1(i,1),1))-2,1),X2(3*(X1_sorted(E1(i,3),1))-2,1))...
            + eox*edgeE1(i,3)/lenE1(i,3);
        res1(X2(3*(X1_sorted(E1(i,1),1))-2,1),1) = res1(X2(3*(X1_sorted(E1(i,1),1))-2,1),1)...
            - eox*(edgeE1(i,1)/lenE1(i,1) + edgeE1(i,3)/lenE1(i,3))*phi(X1_sorted(E1(i,1),1))...
            + eox*edgeE1(i,1)/lenE1(i,1)*phi(X1_sorted(E1(i,2),1))...
            + eox*edgeE1(i,3)/lenE1(i,3)*phi(X1_sorted(E1(i,3),1));

        jac1(X2(3*(X1_sorted(E1(i,2),1))-2,1),X2(3*(X1_sorted(E1(i,2),1))-2,1)) = jac1(X2(3*(X1_sorted(E1(i,2),1))-2,1),X2(3*(X1_sorted(E1(i,2),1))-2,1))...
            - eox*(edgeE1(i,1)/lenE1(i,1) + edgeE1(i,2)/lenE1(i,2));
        jac1(X2(3*(X1_sorted(E1(i,2),1))-2,1),X2(3*(X1_sorted(E1(i,1),1))-2,1)) = jac1(X2(3*(X1_sorted(E1(i,2),1))-2,1),X2(3*(X1_sorted(E1(i,1),1))-2,1))...
            + eox*edgeE1(i,1)/lenE1(i,1);
        jac1(X2(3*(X1_sorted(E1(i,2),1))-2,1),X2(3*(X1_sorted(E1(i,3),1))-2,1)) = jac1(X2(3*(X1_sorted(E1(i,2),1))-2,1),X2(3*(X1_sorted(E1(i,3),1))-2,1))...
            + eox*edgeE1(i,2)/lenE1(i,2);
        res1(X2(3*(X1_sorted(E1(i,2),1))-2,1),1) = res1(X2(3*(X1_sorted(E1(i,2),1))-2,1),1)...
            - eox*(edgeE1(i,1)/lenE1(i,1) + edgeE1(i,2)/lenE1(i,2))*phi(X1_sorted(E1(i,2),1))...
            + eox*edgeE1(i,1)/lenE1(i,1)*phi(X1_sorted(E1(i,1),1))...
            + eox*edgeE1(i,2)/lenE1(i,2)*phi(X1_sorted(E1(i,3),1));

        jac1(X2(3*(X1_sorted(E1(i,3),1))-2,1),X2(3*(X1_sorted(E1(i,3),1))-2,1)) = jac1(X2(3*(X1_sorted(E1(i,3),1))-2,1),X2(3*(X1_sorted(E1(i,3),1))-2,1))...
            - eox*(edgeE1(i,2)/lenE1(i,2) + edgeE1(i,3)/lenE1(i,3));
        jac1(X2(3*(X1_sorted(E1(i,3),1))-2,1),X2(3*(X1_sorted(E1(i,2),1))-2,1)) = jac1(X2(3*(X1_sorted(E1(i,3),1))-2,1),X2(3*(X1_sorted(E1(i,2),1))-2,1))...
            + eox*edgeE1(i,2)/lenE1(i,2);
        jac1(X2(3*(X1_sorted(E1(i,3),1))-2,1),X2(3*(X1_sorted(E1(i,1),1))-2,1)) = jac1(X2(3*(X1_sorted(E1(i,3),1))-2,1),X2(3*(X1_sorted(E1(i,1),1))-2,1))...
            + eox*edgeE1(i,3)/lenE1(i,3);
        res1(X2(3*(X1_sorted(E1(i,3),1))-2,1),1) = res1(X2(3*(X1_sorted(E1(i,3),1))-2,1),1)...
            - eox*(edgeE1(i,2)/lenE1(i,2) + edgeE1(i,3)/lenE1(i,3))*phi(X1_sorted(E1(i,3),1))...
            + eox*edgeE1(i,2)/lenE1(i,2)*phi(X1_sorted(E1(i,2),1))...
            + eox*edgeE1(i,3)/lenE1(i,3)*phi(X1_sorted(E1(i,1),1));
    end

    saveexl = zeros(size(X1,1),1);

    for i=1:Esirow
        if i <= E2row
            jac1(X2(3*(X1_sorted(Esi(i,1),2))-2,1),X2(3*(X1_sorted(Esi(i,1),2))-2,1)) = jac1(X2(3*(X1_sorted(Esi(i,1),2))-2,1),X2(3*(X1_sorted(Esi(i,1),2))-2,1))...
                - esi*(edgeEsi(i,1)/lenEsi(i,1) + edgeEsi(i,3)/lenEsi(i,3));
            jac1(X2(3*(X1_sorted(Esi(i,1),2))-2,1),X2(3*(X1_sorted(Esi(i,2),2))-2,1)) = jac1(X2(3*(X1_sorted(Esi(i,1),2))-2,1),X2(3*(X1_sorted(Esi(i,2),2))-2,1))...
                + esi*edgeEsi(i,1)/lenEsi(i,1);
            jac1(X2(3*(X1_sorted(Esi(i,1),2))-2,1),X2(3*(X1_sorted(Esi(i,3),2))-2,1)) = jac1(X2(3*(X1_sorted(Esi(i,1),2))-2,1),X2(3*(X1_sorted(Esi(i,3),2))-2,1))...
                + esi*edgeEsi(i,3)/lenEsi(i,3);
            res1(X2(3*(X1_sorted(Esi(i,1),2))-2,1),1) =  res1(X2(3*(X1_sorted(Esi(i,1),2))-2,1),1)...
                + (edgeEsi(i,1)*lenEsi(i,1) + edgeEsi(i,3)*lenEsi(i,3))/4*q/e0*Nd...
                - esi*(edgeEsi(i,1)/lenEsi(i,1) + edgeEsi(i,3)/lenEsi(i,3))*phi(X1_sorted(Esi(i,1),2))...
                + esi*edgeEsi(i,1)/lenEsi(i,1)*phi(X1_sorted(Esi(i,2),2))...
                + esi*edgeEsi(i,3)/lenEsi(i,3)*phi(X1_sorted(Esi(i,3),2));
            saveexl(X1_sorted(Esi(i,1),2)) = saveexl(X1_sorted(Esi(i,1),2)) + (edgeEsi(i,1)*lenEsi(i,1) + edgeEsi(i,3)*lenEsi(i,3))/4;

            jac1(X2(3*(X1_sorted(Esi(i,2),2))-2,1),X2(3*(X1_sorted(Esi(i,2),2))-2,1)) = jac1(X2(3*(X1_sorted(Esi(i,2),2))-2,1),X2(3*(X1_sorted(Esi(i,2),2))-2,1))...
                - esi*(edgeEsi(i,1)/lenEsi(i,1) + edgeEsi(i,2)/lenEsi(i,2));
            jac1(X2(3*(X1_sorted(Esi(i,2),2))-2,1),X2(3*(X1_sorted(Esi(i,1),2))-2,1)) = jac1(X2(3*(X1_sorted(Esi(i,2),2))-2,1),X2(3*(X1_sorted(Esi(i,1),2))-2,1))...
                + esi*edgeEsi(i,1)/lenEsi(i,1);
            jac1(X2(3*(X1_sorted(Esi(i,2),2))-2,1),X2(3*(X1_sorted(Esi(i,3),2))-2,1)) = jac1(X2(3*(X1_sorted(Esi(i,2),2))-2,1),X2(3*(X1_sorted(Esi(i,3),2))-2,1))...
                + esi*edgeEsi(i,2)/lenEsi(i,2);
            res1(X2(3*(X1_sorted(Esi(i,2),2))-2,1),1) =  res1(X2(3*(X1_sorted(Esi(i,2),2))-2,1),1)...
                + (edgeEsi(i,1)*lenEsi(i,1) + edgeEsi(i,2)*lenEsi(i,2))/4*q/e0*Nd...
                - esi*(edgeEsi(i,1)/lenEsi(i,1) + edgeEsi(i,2)/lenEsi(i,2))*phi(X1_sorted(Esi(i,2),2))...
                + esi*edgeEsi(i,1)/lenEsi(i,1)*phi(X1_sorted(Esi(i,1),2))...
                + esi*edgeEsi(i,2)/lenEsi(i,2)*phi(X1_sorted(Esi(i,3),2));
            saveexl(X1_sorted(Esi(i,2),2)) = saveexl(X1_sorted(Esi(i,2),2)) + (edgeEsi(i,1)*lenEsi(i,1) + edgeEsi(i,2)*lenEsi(i,2))/4;

            jac1(X2(3*(X1_sorted(Esi(i,3),2))-2,1),X2(3*(X1_sorted(Esi(i,3),2))-2,1)) = jac1(X2(3*(X1_sorted(Esi(i,3),2))-2,1),X2(3*(X1_sorted(Esi(i,3),2))-2,1))...
                - esi*(edgeEsi(i,2)/lenEsi(i,2) + edgeEsi(i,3)/lenEsi(i,3));
            jac1(X2(3*(X1_sorted(Esi(i,3),2))-2,1),X2(3*(X1_sorted(Esi(i,2),2))-2,1)) = jac1(X2(3*(X1_sorted(Esi(i,3),2))-2,1),X2(3*(X1_sorted(Esi(i,2),2))-2,1))...
                + esi*edgeEsi(i,2)/lenEsi(i,2);
            jac1(X2(3*(X1_sorted(Esi(i,3),2))-2,1),X2(3*(X1_sorted(Esi(i,1),2))-2,1)) = jac1(X2(3*(X1_sorted(Esi(i,3),2))-2,1),X2(3*(X1_sorted(Esi(i,1),2))-2,1))...
                + esi*edgeEsi(i,3)/lenEsi(i,3);
            res1(X2(3*(X1_sorted(Esi(i,3),2))-2,1),1) =  res1(X2(3*(X1_sorted(Esi(i,3),2))-2,1),1)...
                + (edgeEsi(i,2)*lenEsi(i,2) + edgeEsi(i,3)*lenEsi(i,3))/4*q/e0*Nd...
                - esi*(edgeEsi(i,2)/lenEsi(i,2) + edgeEsi(i,3)/lenEsi(i,3))*phi(X1_sorted(Esi(i,3),2))...
                + esi*edgeEsi(i,2)/lenEsi(i,2)*phi(X1_sorted(Esi(i,2),2))...
                + esi*edgeEsi(i,3)/lenEsi(i,3)*phi(X1_sorted(Esi(i,1),2));
            saveexl(X1_sorted(Esi(i,3),2)) = saveexl(X1_sorted(Esi(i,3),2)) + (edgeEsi(i,2)*lenEsi(i,2) + edgeEsi(i,3)*lenEsi(i,3))/4;

        elseif i > E2row && i <= (E2row + E3row)
            jac1(X2(3*(X1_sorted(Esi(i,1),2))-2,1),X2(3*(X1_sorted(Esi(i,1),2))-2,1)) = jac1(X2(3*(X1_sorted(Esi(i,1),2))-2,1),X2(3*(X1_sorted(Esi(i,1),2))-2,1))...
                - esi*(edgeEsi(i,1)/lenEsi(i,1) + edgeEsi(i,3)/lenEsi(i,3));
            jac1(X2(3*(X1_sorted(Esi(i,1),2))-2,1),X2(3*(X1_sorted(Esi(i,2),2))-2,1)) = jac1(X2(3*(X1_sorted(Esi(i,1),2))-2,1),X2(3*(X1_sorted(Esi(i,2),2))-2,1))...
                + esi*edgeEsi(i,1)/lenEsi(i,1);
            jac1(X2(3*(X1_sorted(Esi(i,1),2))-2,1),X2(3*(X1_sorted(Esi(i,3),2))-2,1)) = jac1(X2(3*(X1_sorted(Esi(i,1),2))-2,1),X2(3*(X1_sorted(Esi(i,3),2))-2,1))...
                + esi*edgeEsi(i,3)/lenEsi(i,3);
            res1(X2(3*(X1_sorted(Esi(i,1),2))-2,1),1) =  res1(X2(3*(X1_sorted(Esi(i,1),2))-2,1),1)...
                + (edgeEsi(i,1)*lenEsi(i,1) + edgeEsi(i,3)*lenEsi(i,3))/4*q/e0*Na...
                - esi*(edgeEsi(i,1)/lenEsi(i,1) + edgeEsi(i,3)/lenEsi(i,3))*phi(X1_sorted(Esi(i,1),2))...
                + esi*edgeEsi(i,1)/lenEsi(i,1)*phi(X1_sorted(Esi(i,2),2))...
                + esi*edgeEsi(i,3)/lenEsi(i,3)*phi(X1_sorted(Esi(i,3),2));
            saveexl(X1_sorted(Esi(i,1),2)) = saveexl(X1_sorted(Esi(i,1),2)) + (edgeEsi(i,1)*lenEsi(i,1) + edgeEsi(i,3)*lenEsi(i,3))/4;

            jac1(X2(3*(X1_sorted(Esi(i,2),2))-2,1),X2(3*(X1_sorted(Esi(i,2),2))-2,1)) = jac1(X2(3*(X1_sorted(Esi(i,2),2))-2,1),X2(3*(X1_sorted(Esi(i,2),2))-2,1))...
                - esi*(edgeEsi(i,1)/lenEsi(i,1) + edgeEsi(i,2)/lenEsi(i,2));
            jac1(X2(3*(X1_sorted(Esi(i,2),2))-2,1),X2(3*(X1_sorted(Esi(i,1),2))-2,1)) = jac1(X2(3*(X1_sorted(Esi(i,2),2))-2,1),X2(3*(X1_sorted(Esi(i,1),2))-2,1))...
                + esi*edgeEsi(i,1)/lenEsi(i,1);
            jac1(X2(3*(X1_sorted(Esi(i,2),2))-2,1),X2(3*(X1_sorted(Esi(i,3),2))-2,1)) = jac1(X2(3*(X1_sorted(Esi(i,2),2))-2,1),X2(3*(X1_sorted(Esi(i,3),2))-2,1))...
                + esi*edgeEsi(i,2)/lenEsi(i,2);
            res1(X2(3*(X1_sorted(Esi(i,2),2))-2,1),1) =  res1(X2(3*(X1_sorted(Esi(i,2),2))-2,1),1)...
                + (edgeEsi(i,1)*lenEsi(i,1) + edgeEsi(i,2)*lenEsi(i,2))/4*q/e0*Na...
                - esi*(edgeEsi(i,1)/lenEsi(i,1) + edgeEsi(i,2)/lenEsi(i,2))*phi(X1_sorted(Esi(i,2),2))...
                + esi*edgeEsi(i,1)/lenEsi(i,1)*phi(X1_sorted(Esi(i,1),2))...
                + esi*edgeEsi(i,2)/lenEsi(i,2)*phi(X1_sorted(Esi(i,3),2));
            saveexl(X1_sorted(Esi(i,2),2)) = saveexl(X1_sorted(Esi(i,2),2)) + (edgeEsi(i,1)*lenEsi(i,1) + edgeEsi(i,2)*lenEsi(i,2))/4;

            jac1(X2(3*(X1_sorted(Esi(i,3),2))-2,1),X2(3*(X1_sorted(Esi(i,3),2))-2,1)) = jac1(X2(3*(X1_sorted(Esi(i,3),2))-2,1),X2(3*(X1_sorted(Esi(i,3),2))-2,1))...
                - esi*(edgeEsi(i,2)/lenEsi(i,2) + edgeEsi(i,3)/lenEsi(i,3));
            jac1(X2(3*(X1_sorted(Esi(i,3),2))-2,1),X2(3*(X1_sorted(Esi(i,2),2))-2,1)) = jac1(X2(3*(X1_sorted(Esi(i,3),2))-2,1),X2(3*(X1_sorted(Esi(i,2),2))-2,1))...
                + esi*edgeEsi(i,2)/lenEsi(i,2);
            jac1(X2(3*(X1_sorted(Esi(i,3),2))-2,1),X2(3*(X1_sorted(Esi(i,1),2))-2,1)) = jac1(X2(3*(X1_sorted(Esi(i,3),2))-2,1),X2(3*(X1_sorted(Esi(i,1),2))-2,1))...
                + esi*edgeEsi(i,3)/lenEsi(i,3);
            res1(X2(3*(X1_sorted(Esi(i,3),2))-2,1),1) =  res1(X2(3*(X1_sorted(Esi(i,3),2))-2,1),1)...
                + (edgeEsi(i,2)*lenEsi(i,2) + edgeEsi(i,3)*lenEsi(i,3))/4*q/e0*Na...
                - esi*(edgeEsi(i,2)/lenEsi(i,2) + edgeEsi(i,3)/lenEsi(i,3))*phi(X1_sorted(Esi(i,3),2))...
                + esi*edgeEsi(i,2)/lenEsi(i,2)*phi(X1_sorted(Esi(i,2),2))...
                + esi*edgeEsi(i,3)/lenEsi(i,3)*phi(X1_sorted(Esi(i,1),2));
            saveexl(X1_sorted(Esi(i,3),2)) = saveexl(X1_sorted(Esi(i,3),2)) + (edgeEsi(i,2)*lenEsi(i,2) + edgeEsi(i,3)*lenEsi(i,3))/4;

        elseif i > (E2row+E3row) && i <= Esirow
            jac1(X2(3*(X1_sorted(Esi(i,1),2))-2,1),X2(3*(X1_sorted(Esi(i,1),2))-2,1)) = jac1(X2(3*(X1_sorted(Esi(i,1),2))-2,1),X2(3*(X1_sorted(Esi(i,1),2))-2,1))...
                - esi*(edgeEsi(i,1)/lenEsi(i,1) + edgeEsi(i,3)/lenEsi(i,3));
            jac1(X2(3*(X1_sorted(Esi(i,1),2))-2,1),X2(3*(X1_sorted(Esi(i,2),2))-2,1)) = jac1(X2(3*(X1_sorted(Esi(i,1),2))-2,1),X2(3*(X1_sorted(Esi(i,2),2))-2,1))...
                + esi*edgeEsi(i,1)/lenEsi(i,1);
            jac1(X2(3*(X1_sorted(Esi(i,1),2))-2,1),X2(3*(X1_sorted(Esi(i,3),2))-2,1)) = jac1(X2(3*(X1_sorted(Esi(i,1),2))-2,1),X2(3*(X1_sorted(Esi(i,3),2))-2,1))...
                + esi*edgeEsi(i,3)/lenEsi(i,3);
            res1(X2(3*(X1_sorted(Esi(i,1),2))-2,1),1) =  res1(X2(3*(X1_sorted(Esi(i,1),2))-2,1),1)...
                + (edgeEsi(i,1)*lenEsi(i,1) + edgeEsi(i,3)*lenEsi(i,3))/4*q/e0*Nd...
                - esi*(edgeEsi(i,1)/lenEsi(i,1) + edgeEsi(i,3)/lenEsi(i,3))*phi(X1_sorted(Esi(i,1),2))...
                + esi*edgeEsi(i,1)/lenEsi(i,1)*phi(X1_sorted(Esi(i,2),2))...
                + esi*edgeEsi(i,3)/lenEsi(i,3)*phi(X1_sorted(Esi(i,3),2));
            saveexl(X1_sorted(Esi(i,1),2)) = saveexl(X1_sorted(Esi(i,1),2)) + (edgeEsi(i,1)*lenEsi(i,1) + edgeEsi(i,3)*lenEsi(i,3))/4;

            jac1(X2(3*(X1_sorted(Esi(i,2),2))-2,1),X2(3*(X1_sorted(Esi(i,2),2))-2,1)) = jac1(X2(3*(X1_sorted(Esi(i,2),2))-2,1),X2(3*(X1_sorted(Esi(i,2),2))-2,1))...
                - esi*(edgeEsi(i,1)/lenEsi(i,1) + edgeEsi(i,2)/lenEsi(i,2));
            jac1(X2(3*(X1_sorted(Esi(i,2),2))-2,1),X2(3*(X1_sorted(Esi(i,1),2))-2,1)) = jac1(X2(3*(X1_sorted(Esi(i,2),2))-2,1),X2(3*(X1_sorted(Esi(i,1),2))-2,1))...
                + esi*edgeEsi(i,1)/lenEsi(i,1);
            jac1(X2(3*(X1_sorted(Esi(i,2),2))-2,1),X2(3*(X1_sorted(Esi(i,3),2))-2,1)) = jac1(X2(3*(X1_sorted(Esi(i,2),2))-2,1),X2(3*(X1_sorted(Esi(i,3),2))-2,1))...
                + esi*edgeEsi(i,2)/lenEsi(i,2);
            res1(X2(3*(X1_sorted(Esi(i,2),2))-2,1),1) =  res1(X2(3*(X1_sorted(Esi(i,2),2))-2,1),1)...
                + (edgeEsi(i,1)*lenEsi(i,1) + edgeEsi(i,2)*lenEsi(i,2))/4*q/e0*Nd...
                - esi*(edgeEsi(i,1)/lenEsi(i,1) + edgeEsi(i,2)/lenEsi(i,2))*phi(X1_sorted(Esi(i,2),2))...
                + esi*edgeEsi(i,1)/lenEsi(i,1)*phi(X1_sorted(Esi(i,1),2))...
                + esi*edgeEsi(i,2)/lenEsi(i,2)*phi(X1_sorted(Esi(i,3),2));
            saveexl(X1_sorted(Esi(i,2),2)) = saveexl(X1_sorted(Esi(i,2),2)) + (edgeEsi(i,1)*lenEsi(i,1) + edgeEsi(i,2)*lenEsi(i,2))/4;

            jac1(X2(3*X1_sorted(Esi(i,3),2)-2,1),X2(3*(X1_sorted(Esi(i,3),2))-2,1)) = jac1(X2(3*(X1_sorted(Esi(i,3),2))-2,1),X2(3*(X1_sorted(Esi(i,3),2))-2,1))...
                - esi*(edgeEsi(i,2)/lenEsi(i,2) + edgeEsi(i,3)/lenEsi(i,3));
            jac1(X2(3*X1_sorted(Esi(i,3),2)-2,1),X2(3*(X1_sorted(Esi(i,2),2))-2,1)) = jac1(X2(3*(X1_sorted(Esi(i,3),2))-2,1),X2(3*(X1_sorted(Esi(i,2),2))-2,1))...
                + esi*edgeEsi(i,2)/lenEsi(i,2);
            jac1(X2(3*X1_sorted(Esi(i,3),2)-2,1),X2(3*(X1_sorted(Esi(i,1),2))-2,1)) = jac1(X2(3*(X1_sorted(Esi(i,3),2))-2,1),X2(3*(X1_sorted(Esi(i,1),2))-2,1))...
                + esi*edgeEsi(i,3)/lenEsi(i,3);
            res1(X2(3*X1_sorted(Esi(i,3),2)-2,1),1) =  res1(X2(3*X1_sorted(Esi(i,3),2)-2,1),1)...
                + (edgeEsi(i,2)*lenEsi(i,2) + edgeEsi(i,3)*lenEsi(i,3))/4*q/e0*Nd...
                - esi*(edgeEsi(i,2)/lenEsi(i,2) + edgeEsi(i,3)/lenEsi(i,3))*phi(X1_sorted(Esi(i,3),2))...
                + esi*edgeEsi(i,2)/lenEsi(i,2)*phi(X1_sorted(Esi(i,2),2))...
                + esi*edgeEsi(i,3)/lenEsi(i,3)*phi(X1_sorted(Esi(i,1),2));
            saveexl(X1_sorted(Esi(i,3),2)) = saveexl(X1_sorted(Esi(i,3),2)) + (edgeEsi(i,2)*lenEsi(i,2) + edgeEsi(i,3)*lenEsi(i,3))/4;
        end
    end

    for i=1:E5row
        jac1(X2(3*(X1_sorted(E5(i,1),3))-2,1),X2(3*(X1_sorted(E5(i,1),3))-2,1)) = jac1(X2(3*(X1_sorted(E5(i,1),3))-2,1),X2(3*(X1_sorted(E5(i,1),3))-2,1))...
            - eox*(edgeE5(i,1)/lenE5(i,1) + edgeE5(i,3)/lenE5(i,3));
        jac1(X2(3*(X1_sorted(E5(i,1),3))-2,1),X2(3*(X1_sorted(E5(i,2),3))-2,1)) = jac1(X2(3*(X1_sorted(E5(i,1),3))-2,1),X2(3*(X1_sorted(E5(i,2),3))-2,1))...
            + eox*edgeE5(i,1)/lenE5(i,1);
        jac1(X2(3*(X1_sorted(E5(i,1),3))-2,1),X2(3*(X1_sorted(E5(i,3),3))-2,1)) = jac1(X2(3*(X1_sorted(E5(i,1),3))-2,1),X2(3*(X1_sorted(E5(i,3),3))-2,1))...
            + eox*edgeE5(i,3)/lenE5(i,3);
        res1(X2(3*(X1_sorted(E5(i,1),3))-2,1),1) = res1(X2(3*(X1_sorted(E5(i,1),3))-2,1),1)...
            - eox*(edgeE5(i,1)/lenE5(i,1) + edgeE5(i,3)/lenE5(i,3))*phi(X1_sorted(E5(i,1),3))...
            + eox*edgeE5(i,1)/lenE5(i,1)*phi(X1_sorted(E5(i,2),3))...
            + eox*edgeE5(i,3)/lenE5(i,3)*phi(X1_sorted(E5(i,3),3));

        jac1(X2(3*(X1_sorted(E5(i,2),3))-2,1),X2(3*(X1_sorted(E5(i,2),3))-2,1)) = jac1(X2(3*(X1_sorted(E5(i,2),3))-2,1),X2(3*(X1_sorted(E5(i,2),3))-2,1))...
            - eox*(edgeE5(i,1)/lenE5(i,1) + edgeE5(i,2)/lenE5(i,2));
        jac1(X2(3*(X1_sorted(E5(i,2),3))-2,1),X2(3*(X1_sorted(E5(i,1),3))-2,1)) = jac1(X2(3*(X1_sorted(E5(i,2),3))-2,1),X2(3*(X1_sorted(E5(i,1),3))-2,1))...
            + eox*edgeE5(i,1)/lenE5(i,1);
        jac1(X2(3*(X1_sorted(E5(i,2),3))-2,1),X2(3*(X1_sorted(E5(i,3),3))-2,1)) = jac1(X2(3*(X1_sorted(E5(i,2),3))-2,1),X2(3*(X1_sorted(E5(i,3),3))-2,1))...
            + eox*edgeE5(i,2)/lenE5(i,2);
        res1(X2(3*(X1_sorted(E5(i,2),3))-2,1),1) = res1(X2(3*(X1_sorted(E5(i,2),3))-2,1),1)...
            - eox*(edgeE5(i,1)/lenE5(i,1) + edgeE5(i,2)/lenE5(i,2))*phi(X1_sorted(E5(i,2),3))...
            + eox*edgeE5(i,1)/lenE5(i,1)*phi(X1_sorted(E5(i,1),3))...
            + eox*edgeE5(i,2)/lenE5(i,2)*phi(X1_sorted(E5(i,3),3));

        jac1(X2(3*(X1_sorted(E5(i,3),3))-2,1),X2(3*(X1_sorted(E5(i,3),3))-2,1)) = jac1(X2(3*(X1_sorted(E5(i,3),3))-2,1),X2(3*(X1_sorted(E5(i,3),3))-2,1))...
            - eox*(edgeE5(i,2)/lenE5(i,2) + edgeE5(i,3)/lenE5(i,3));
        jac1(X2(3*(X1_sorted(E5(i,3),3))-2,1),X2(3*(X1_sorted(E5(i,2),3))-2,1)) = jac1(X2(3*(X1_sorted(E5(i,3),3))-2,1),X2(3*(X1_sorted(E5(i,2),3))-2,1))...
            + eox*edgeE5(i,2)/lenE5(i,2);
        jac1(X2(3*(X1_sorted(E5(i,3),3))-2,1),X2(3*(X1_sorted(E5(i,1),3))-2,1)) = jac1(X2(3*(X1_sorted(E5(i,3),3))-2,1),X2(3*(X1_sorted(E5(i,1),3))-2,1))...
            + eox*edgeE5(i,3)/lenE5(i,3);
        res1(X2(3*(X1_sorted(E5(i,3),3))-2,1),1) = res1(X2(3*(X1_sorted(E5(i,3),3))-2,1),1)...
            - eox*(edgeE5(i,2)/lenE5(i,2) + edgeE5(i,3)/lenE5(i,3))*phi(X1_sorted(E5(i,3),3))...
            + eox*edgeE5(i,2)/lenE5(i,2)*phi(X1_sorted(E5(i,2),3))...
            + eox*edgeE5(i,3)/lenE5(i,3)*phi(X1_sorted(E5(i,1),3));
    end

    %%%  for electron
    for i=1:V_silicon_row
        %%% n-nit*exp(phi/Vt)
        jac1(X2(3*V_silicon(i,2)-1),X2(3*V_silicon(i,2)-1)) = 1;
        jac1(X2(3*V_silicon(i,2)-1),X2(3*V_silicon(i,2)-2)) = - nint/Vt*exp(phi(V_silicon(i,2))/Vt);
        res1(X2(3*V_silicon(i,2)-1),1) = n(i) - nint*exp(phi(V_silicon(i,2))/Vt);
        %%% - qn
        jac1(X2(3*V_silicon(i,2)-2),X2(3*V_silicon(i,2)-1)) = - saveexl(V_silicon(i,2))*q/e0;
        res1(X2(3*V_silicon(i,2)-2),1) = res1(X2(3*V_silicon(i,2)-2),1) - saveexl(V_silicon(i,2))*q/e0*n(i);
    end

    %%%  for hole
    for i=1:V_silicon_row
        %%% p-nit*exp(-phi/Vt)
        jac1(X2(3*V_silicon(i,2)),X2(3*V_silicon(i,2))) = 1;
        jac1(X2(3*V_silicon(i,2)),X2(3*V_silicon(i,2)-2)) = nint/Vt*exp(-phi(V_silicon(i,2))/Vt);
        res1(X2(3*V_silicon(i,2)),1) = p(i) - nint*exp(-phi(V_silicon(i,2))/Vt);
        %%% + qp
        jac1(X2(3*V_silicon(i,2)-2),X2(3*V_silicon(i,2))) = saveexl(V_silicon(i,2))*q/e0;
        res1(X2(3*V_silicon(i,2)-2),1) = res1(X2(3*V_silicon(i,2)-2),1) + saveexl(V_silicon(i,2))*q/e0*p(i);
    end

    for i=1:is1row
        minus1 = X1_sorted(is1(i),1);
        plus1 = X1_sorted(is1(i),2);
        jac1(X2(3*plus1-2),:) = jac1(X2(3*plus1-2),:) + jac1(X2(3*minus1-2),:);
        res1(X2(3*plus1-2)) = res1(X2(3*plus1-2)) + res1(X2(3*minus1-2));
        jac1(X2(3*minus1-2),:) = 0;
        res1(X2(3*minus1-2),1) = 0;
        jac1(X2(3*minus1-2), X2(3*minus1-2)) = -1;
        jac1(X2(3*minus1-2), X2(3*plus1-2)) = 1;
    end

    for i=1:is2row
        minus1 = X1_sorted(is2(i),3);
        plus1 = X1_sorted(is2(i),2);
        jac1(X2(3*plus1-2),:) = jac1(X2(3*plus1-2),:) + jac1(X2(3*minus1-2),:);
        res1(X2(3*plus1-2),1) = res1(X2(3*plus1-2),1) + res1(X2(3*minus1-2),1);
        jac1(X2(3*minus1-2),:) = 0;
        res1(X2(3*minus1-2),1) = 0;
        jac1(X2(3*minus1-2), X2(3*minus1-2)) = -1;
        jac1(X2(3*minus1-2), X2(3*plus1-2)) = 1;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Dirichlet boundary condition for 1st iteration
    % gate1
    for i=1:Ncon1
        for j=1:concol
            fcon = X1_sorted(contact(i,j),1);
            jac1(X2(3*fcon-2),:) = 0;
            jac1(X2(3*fcon-2), X2(3*fcon-2)) = 1;
            res1(X2(3*fcon-2),1) = phi(fcon,1) - 0.33374;
        end
    end

    % gate2
    for i=Ncon1+1:(Ncon1+Ncon2)
        for j=1:concol
            fcon = X1_sorted(contact(i,j),3);
            jac1(X2(3*fcon-2),:) = 0;
            jac1(X2(3*fcon-2), X2(3*fcon-2)) = 1;
            res1(X2(3*fcon-2),1) = phi(fcon) - 0.33374;
        end
    end

    % source & drain
    for i=(Ncon1+Ncon2)+1:(Ncon1+Ncon2+Ncon3+Ncon4)
        for j=1:concol
            fcon = X1_sorted(contact(i,j),2);
            fcon_np = find(contact(i,j)==V_silicon(:,1));

            jac1(X2(3*fcon-2),:) = 0;
            jac1(X2(3*fcon-1),:) = 0;
            jac1(X2(3*fcon),:) = 0;

            jac1(X2(3*fcon-2), X2(3*fcon-2)) = 1;
            jac1(X2(3*fcon-1), X2(3*fcon-1)) = 1;
            jac1(X2(3*fcon), X2(3*fcon)) = 1;

            % source
            if i<=(Ncon1+Ncon2+Ncon3)
                res1(X2(3*fcon-2),1) = phi(fcon) - Vt*log(Nd/nint);
                res1(X2(3*fcon-1),1) = n(fcon_np) - Nd;
                res1(X2(3*fcon),1) = p(fcon_np) - nint*nint/Nd;

                % drain
            elseif i>(Ncon1+Ncon2+Ncon3) &&  i<=(Ncon1+Ncon2+Ncon3+Ncon4)
                res1(X2(3*fcon-2),1) = phi(fcon) - Vt*log(Nd/nint);
                res1(X2(3*fcon-1),1) = n(fcon_np) - Nd;
                res1(X2(3*fcon),1) = p(fcon_np) - nint*nint/Nd;
            end
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% scaled part
    cvector = zeros(size(jac1,1),1);
    for i=1:size(jac1,1)
        if rem(X2_sorted(i),3) == 1
            cvector(i) = Vt;
        elseif rem(X2_sorted(i),3) == 2
            cvector(i) = Nd;
        elseif rem(X2_sorted(i),3) == 0
            cvector(i) = Nd;
        end
    end

    cmatrix = spdiags(cvector,0,size(X2_sorted,1),size(X2_sorted,1));
    jac1_scaled = jac1 * cmatrix;
    rvector = 1./sum(abs(jac1_scaled),2);
    rmatrix = spdiags(rvector,0,size(X2_sorted,1),size(X2_sorted,1));
    jac1_scaled = rmatrix * jac1_scaled;
    res1_scaled = rmatrix * res1;
    upd_scaled = jac1_scaled \ (-res1_scaled);
    upd = cmatrix * upd_scaled;

    countphi1=1;
    countn1=1;
    countp1=1;

    maxupdate1 = 0;
    for i=1:size(X2_sorted,1)
        if rem(X2_sorted(i),3) == 1
            phi(countphi1) = phi(countphi1) + upd(i);
            if maxupdate1 < max(abs(upd(i)))
                maxupdate1 = max(abs(upd(i)));
            end
            countphi1 = countphi1+1;
        elseif rem(X2_sorted(i),3) == 2
            n(countn1) = n(countn1) + upd(i);
            countn1 = countn1 + 1;
        elseif rem(X2_sorted(i),3) == 0
            p(countp1) = p(countp1) + upd(i);
            countp1 = countp1 + 1;
        end
    end

    save_upd(:,it) = upd(:,1);
    save_phi_upd(it,1) = maxupdate1;
end

initial_phi_1iter = phi(:,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% jacobian, residue matrix (DD)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Final iteration number
Niter2 = 50;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ndstep = 2;
save_drain_phi = zeros(size(phi,1),Ndstep);
save_drain_n = zeros(size(n,1), Ndstep);
save_drain_p = zeros(size(p,1), Ndstep);

save_drain_phi_update = zeros(Niter2,1);
save_drain_update = zeros(size(X2_sorted,1),1);

I = zeros(Ndstep,1);

for step=1:Ndstep

    Vd = 0.001*(step-1);

    for it=1:Niter2

        Jaco = sparse(size(X2_sorted,1),size(X2_sorted,1));
        Res = zeros(size(X2_sorted,1),1);

        for i=1:E1row
            Jaco(X2(3*(X1_sorted(E1(i,1),1))-2,1), X2(3*(X1_sorted(E1(i,1),1))-2,1)) = Jaco(X2(3*(X1_sorted(E1(i,1),1))-2,1), X2(3*(X1_sorted(E1(i,1),1))-2,1))...
                - eox*(edgeE1(i,1)/lenE1(i,1) + edgeE1(i,3)/lenE1(i,3));
            Jaco(X2(3*(X1_sorted(E1(i,1),1))-2,1), X2(3*(X1_sorted(E1(i,2),1))-2,1)) = Jaco(X2(3*(X1_sorted(E1(i,1),1))-2,1), X2(3*(X1_sorted(E1(i,2),1))-2,1))...
                + eox*edgeE1(i,1)/lenE1(i,1);
            Jaco(X2(3*(X1_sorted(E1(i,1),1))-2,1), X2(3*(X1_sorted(E1(i,3),1))-2,1)) = Jaco(X2(3*(X1_sorted(E1(i,1),1))-2,1), X2(3*(X1_sorted(E1(i,3),1))-2,1))...
                + eox*edgeE1(i,3)/lenE1(i,3);
            Res(X2(3*(X1_sorted(E1(i,1),1))-2,1),1) = Res(X2(3*(X1_sorted(E1(i,1),1))-2,1),1)...
                - eox*(edgeE1(i,1)/lenE1(i,1) + edgeE1(i,3)/lenE1(i,3))*phi(X1_sorted(E1(i,1),1))...
                + eox*edgeE1(i,1)/lenE1(i,1)*phi(X1_sorted(E1(i,2),1))...
                + eox*edgeE1(i,3)/lenE1(i,3)*phi(X1_sorted(E1(i,3),1));

            Jaco(X2(3*(X1_sorted(E1(i,2),1))-2,1), X2(3*(X1_sorted(E1(i,2),1))-2,1)) = Jaco(X2(3*(X1_sorted(E1(i,2),1))-2,1), X2(3*(X1_sorted(E1(i,2),1))-2,1))...
                - eox*(edgeE1(i,1)/lenE1(i,1) + edgeE1(i,2)/lenE1(i,2));
            Jaco(X2(3*(X1_sorted(E1(i,2),1))-2,1), X2(3*(X1_sorted(E1(i,1),1))-2,1)) = Jaco(X2(3*(X1_sorted(E1(i,2),1))-2,1), X2(3*(X1_sorted(E1(i,1),1))-2,1))...
                + eox*edgeE1(i,1)/lenE1(i,1);
            Jaco(X2(3*(X1_sorted(E1(i,2),1))-2,1), X2(3*(X1_sorted(E1(i,3),1))-2,1)) = Jaco(X2(3*(X1_sorted(E1(i,2),1))-2,1), X2(3*(X1_sorted(E1(i,3),1))-2,1))...
                + eox*edgeE1(i,2)/lenE1(i,2);
            Res(X2(3*(X1_sorted(E1(i,2),1))-2,1),1) = Res(X2(3*(X1_sorted(E1(i,2),1))-2,1),1)...
                - eox*(edgeE1(i,1)/lenE1(i,1) + edgeE1(i,2)/lenE1(i,2))*phi(X1_sorted(E1(i,2),1))...
                + eox*edgeE1(i,1)/lenE1(i,1)*phi(X1_sorted(E1(i,1),1))...
                + eox*edgeE1(i,2)/lenE1(i,2)*phi(X1_sorted(E1(i,3),1));

            Jaco(X2(3*(X1_sorted(E1(i,3),1))-2,1), X2(3*(X1_sorted(E1(i,3),1))-2,1)) = Jaco(X2(3*(X1_sorted(E1(i,3),1))-2,1), X2(3*(X1_sorted(E1(i,3),1))-2,1))...
                - eox*(edgeE1(i,2)/lenE1(i,2) + edgeE1(i,3)/lenE1(i,3));
            Jaco(X2(3*(X1_sorted(E1(i,3),1))-2,1), X2(3*(X1_sorted(E1(i,2),1))-2,1)) = Jaco(X2(3*(X1_sorted(E1(i,3),1))-2,1), X2(3*(X1_sorted(E1(i,2),1))-2,1))...
                + eox*edgeE1(i,2)/lenE1(i,2);
            Jaco(X2(3*(X1_sorted(E1(i,3),1))-2,1), X2(3*(X1_sorted(E1(i,1),1))-2,1)) = Jaco(X2(3*(X1_sorted(E1(i,3),1))-2,1), X2(3*(X1_sorted(E1(i,1),1))-2,1))...
                + eox*edgeE1(i,3)/lenE1(i,3);
            Res(X2(3*(X1_sorted(E1(i,3),1))-2,1),1) = Res(X2(3*(X1_sorted(E1(i,3),1))-2,1),1)...
                - eox*(edgeE1(i,2)/lenE1(i,2) + edgeE1(i,3)/lenE1(i,3))*phi(X1_sorted(E1(i,3),1))...
                + eox*edgeE1(i,2)/lenE1(i,2)*phi(X1_sorted(E1(i,2),1))...
                + eox*edgeE1(i,3)/lenE1(i,3)*phi(X1_sorted(E1(i,1),1));
        end

        saveexl = zeros(size(X1,1),1);

        for i=1:Esirow
            if i <= E2row
                Jaco(X2(3*(X1_sorted(Esi(i,1),2))-2,1), X2(3*(X1_sorted(Esi(i,1),2))-2,1)) = Jaco(X2(3*(X1_sorted(Esi(i,1),2))-2,1), X2(3*(X1_sorted(Esi(i,1),2))-2,1))...
                    - esi*(edgeEsi(i,1)/lenEsi(i,1) + edgeEsi(i,3)/lenEsi(i,3));
                Jaco(X2(3*(X1_sorted(Esi(i,1),2))-2,1), X2(3*(X1_sorted(Esi(i,2),2))-2,1)) = Jaco(X2(3*(X1_sorted(Esi(i,1),2))-2,1), X2(3*(X1_sorted(Esi(i,2),2))-2,1))...
                    + esi*edgeEsi(i,1)/lenEsi(i,1);
                Jaco(X2(3*(X1_sorted(Esi(i,1),2))-2,1), X2(3*(X1_sorted(Esi(i,3),2))-2,1)) = Jaco(X2(3*(X1_sorted(Esi(i,1),2))-2,1), X2(3*(X1_sorted(Esi(i,3),2))-2,1))...
                    + esi*edgeEsi(i,3)/lenEsi(i,3);
                Res(X2(3*(X1_sorted(Esi(i,1),2))-2,1),1) =  Res(X2(3*(X1_sorted(Esi(i,1),2))-2,1),1)...
                    + (edgeEsi(i,1)*lenEsi(i,1) + edgeEsi(i,3)*lenEsi(i,3))/4*q/e0*Nd...
                    - esi*(edgeEsi(i,1)/lenEsi(i,1) + edgeEsi(i,3)/lenEsi(i,3))*phi(X1_sorted(Esi(i,1),2))...
                    + esi*edgeEsi(i,1)/lenEsi(i,1)*phi(X1_sorted(Esi(i,2),2))...
                    + esi*edgeEsi(i,3)/lenEsi(i,3)*phi(X1_sorted(Esi(i,3),2));
                saveexl(X1_sorted(Esi(i,1),2)) = saveexl(X1_sorted(Esi(i,1),2)) + (edgeEsi(i,1)*lenEsi(i,1) + edgeEsi(i,3)*lenEsi(i,3))/4;

                Jaco(X2(3*(X1_sorted(Esi(i,2),2))-2,1), X2(3*(X1_sorted(Esi(i,2),2))-2,1)) = Jaco(X2(3*(X1_sorted(Esi(i,2),2))-2,1), X2(3*(X1_sorted(Esi(i,2),2))-2,1))...
                    - esi*(edgeEsi(i,1)/lenEsi(i,1) + edgeEsi(i,2)/lenEsi(i,2));
                Jaco(X2(3*(X1_sorted(Esi(i,2),2))-2,1), X2(3*(X1_sorted(Esi(i,1),2))-2,1)) = Jaco(X2(3*(X1_sorted(Esi(i,2),2))-2,1), X2(3*(X1_sorted(Esi(i,1),2))-2,1))...
                    + esi*edgeEsi(i,1)/lenEsi(i,1);
                Jaco(X2(3*(X1_sorted(Esi(i,2),2))-2,1), X2(3*(X1_sorted(Esi(i,3),2))-2,1)) = Jaco(X2(3*(X1_sorted(Esi(i,2),2))-2,1), X2(3*(X1_sorted(Esi(i,3),2))-2,1))...
                    + esi*edgeEsi(i,2)/lenEsi(i,2);
                Res(X2(3*(X1_sorted(Esi(i,2),2))-2,1),1) =  Res(X2(3*(X1_sorted(Esi(i,2),2))-2,1),1)...
                    + (edgeEsi(i,1)*lenEsi(i,1) + edgeEsi(i,2)*lenEsi(i,2))/4*q/e0*Nd...
                    - esi*(edgeEsi(i,1)/lenEsi(i,1) + edgeEsi(i,2)/lenEsi(i,2))*phi(X1_sorted(Esi(i,2),2))...
                    + esi*edgeEsi(i,1)/lenEsi(i,1)*phi(X1_sorted(Esi(i,1),2))...
                    + esi*edgeEsi(i,2)/lenEsi(i,2)*phi(X1_sorted(Esi(i,3),2));
                saveexl(X1_sorted(Esi(i,2),2)) = saveexl(X1_sorted(Esi(i,2),2)) + (edgeEsi(i,1)*lenEsi(i,1) + edgeEsi(i,2)*lenEsi(i,2))/4;

                Jaco(X2(3*(X1_sorted(Esi(i,3),2))-2,1), X2(3*(X1_sorted(Esi(i,3),2))-2,1)) = Jaco(X2(3*(X1_sorted(Esi(i,3),2))-2,1), X2(3*(X1_sorted(Esi(i,3),2))-2,1))...
                    - esi*(edgeEsi(i,2)/lenEsi(i,2) + edgeEsi(i,3)/lenEsi(i,3));
                Jaco(X2(3*(X1_sorted(Esi(i,3),2))-2,1), X2(3*(X1_sorted(Esi(i,2),2))-2,1)) = Jaco(X2(3*(X1_sorted(Esi(i,3),2))-2,1), X2(3*(X1_sorted(Esi(i,2),2))-2,1))...
                    + esi*edgeEsi(i,2)/lenEsi(i,2);
                Jaco(X2(3*(X1_sorted(Esi(i,3),2))-2,1), X2(3*(X1_sorted(Esi(i,1),2))-2,1)) = Jaco(X2(3*(X1_sorted(Esi(i,3),2))-2,1), X2(3*(X1_sorted(Esi(i,1),2))-2,1))...
                    + esi*edgeEsi(i,3)/lenEsi(i,3);
                Res(X2(3*(X1_sorted(Esi(i,3),2))-2,1),1) =  Res(X2(3*(X1_sorted(Esi(i,3),2))-2,1),1)...
                    + (edgeEsi(i,2)*lenEsi(i,2) + edgeEsi(i,3)*lenEsi(i,3))/4*q/e0*Nd...
                    - esi*(edgeEsi(i,2)/lenEsi(i,2) + edgeEsi(i,3)/lenEsi(i,3))*phi(X1_sorted(Esi(i,3),2))...
                    + esi*edgeEsi(i,2)/lenEsi(i,2)*phi(X1_sorted(Esi(i,2),2))...
                    + esi*edgeEsi(i,3)/lenEsi(i,3)*phi(X1_sorted(Esi(i,1),2));
                saveexl(X1_sorted(Esi(i,3),2)) = saveexl(X1_sorted(Esi(i,3),2)) + (edgeEsi(i,2)*lenEsi(i,2) + edgeEsi(i,3)*lenEsi(i,3))/4;

            elseif i > E2row && i <= (E2row + E3row)
                Jaco(X2(3*(X1_sorted(Esi(i,1),2))-2,1), X2(3*(X1_sorted(Esi(i,1),2))-2,1)) = Jaco(X2(3*(X1_sorted(Esi(i,1),2))-2,1), X2(3*(X1_sorted(Esi(i,1),2))-2,1))...
                    - esi*(edgeEsi(i,1)/lenEsi(i,1) + edgeEsi(i,3)/lenEsi(i,3));
                Jaco(X2(3*(X1_sorted(Esi(i,1),2))-2,1), X2(3*(X1_sorted(Esi(i,2),2))-2,1)) = Jaco(X2(3*(X1_sorted(Esi(i,1),2))-2,1), X2(3*(X1_sorted(Esi(i,2),2))-2,1))...
                    + esi*edgeEsi(i,1)/lenEsi(i,1);
                Jaco(X2(3*(X1_sorted(Esi(i,1),2))-2,1), X2(3*(X1_sorted(Esi(i,3),2))-2,1)) = Jaco(X2(3*(X1_sorted(Esi(i,1),2))-2,1), X2(3*(X1_sorted(Esi(i,3),2))-2,1))...
                    + esi*edgeEsi(i,3)/lenEsi(i,3);
                Res(X2(3*(X1_sorted(Esi(i,1),2))-2,1),1) =  Res(X2(3*(X1_sorted(Esi(i,1),2))-2,1),1)...
                    + (edgeEsi(i,1)*lenEsi(i,1) + edgeEsi(i,3)*lenEsi(i,3))/4*q/e0*Na...
                    - esi*(edgeEsi(i,1)/lenEsi(i,1) + edgeEsi(i,3)/lenEsi(i,3))*phi(X1_sorted(Esi(i,1),2))...
                    + esi*edgeEsi(i,1)/lenEsi(i,1)*phi(X1_sorted(Esi(i,2),2))...
                    + esi*edgeEsi(i,3)/lenEsi(i,3)*phi(X1_sorted(Esi(i,3),2));
                saveexl(X1_sorted(Esi(i,1),2)) = saveexl(X1_sorted(Esi(i,1),2)) + (edgeEsi(i,1)*lenEsi(i,1) + edgeEsi(i,3)*lenEsi(i,3))/4;

                Jaco(X2(3*(X1_sorted(Esi(i,2),2))-2,1), X2(3*(X1_sorted(Esi(i,2),2))-2,1)) = Jaco(X2(3*(X1_sorted(Esi(i,2),2))-2,1), X2(3*(X1_sorted(Esi(i,2),2))-2,1))...
                    - esi*(edgeEsi(i,1)/lenEsi(i,1) + edgeEsi(i,2)/lenEsi(i,2));
                Jaco(X2(3*(X1_sorted(Esi(i,2),2))-2,1), X2(3*(X1_sorted(Esi(i,1),2))-2,1)) = Jaco(X2(3*(X1_sorted(Esi(i,2),2))-2,1), X2(3*(X1_sorted(Esi(i,1),2))-2,1))...
                    + esi*edgeEsi(i,1)/lenEsi(i,1);
                Jaco(X2(3*(X1_sorted(Esi(i,2),2))-2,1), X2(3*(X1_sorted(Esi(i,3),2))-2,1)) = Jaco(X2(3*(X1_sorted(Esi(i,2),2))-2,1), X2(3*(X1_sorted(Esi(i,3),2))-2,1))...
                    + esi*edgeEsi(i,2)/lenEsi(i,2);
                Res(X2(3*(X1_sorted(Esi(i,2),2))-2,1),1) =  Res(X2(3*(X1_sorted(Esi(i,2),2))-2,1),1)...
                    + (edgeEsi(i,1)*lenEsi(i,1) + edgeEsi(i,2)*lenEsi(i,2))/4*q/e0*Na...
                    - esi*(edgeEsi(i,1)/lenEsi(i,1) + edgeEsi(i,2)/lenEsi(i,2))*phi(X1_sorted(Esi(i,2),2))...
                    + esi*edgeEsi(i,1)/lenEsi(i,1)*phi(X1_sorted(Esi(i,1),2))...
                    + esi*edgeEsi(i,2)/lenEsi(i,2)*phi(X1_sorted(Esi(i,3),2));
                saveexl(X1_sorted(Esi(i,2),2)) = saveexl(X1_sorted(Esi(i,2),2)) + (edgeEsi(i,1)*lenEsi(i,1) + edgeEsi(i,2)*lenEsi(i,2))/4;

                Jaco(X2(3*(X1_sorted(Esi(i,3),2))-2,1), X2(3*(X1_sorted(Esi(i,3),2))-2,1)) = Jaco(X2(3*(X1_sorted(Esi(i,3),2))-2,1), X2(3*(X1_sorted(Esi(i,3),2))-2,1))...
                    - esi*(edgeEsi(i,2)/lenEsi(i,2) + edgeEsi(i,3)/lenEsi(i,3));
                Jaco(X2(3*(X1_sorted(Esi(i,3),2))-2,1), X2(3*(X1_sorted(Esi(i,2),2))-2,1)) = Jaco(X2(3*(X1_sorted(Esi(i,3),2))-2,1), X2(3*(X1_sorted(Esi(i,2),2))-2,1))...
                    + esi*edgeEsi(i,2)/lenEsi(i,2);
                Jaco(X2(3*(X1_sorted(Esi(i,3),2))-2,1), X2(3*(X1_sorted(Esi(i,1),2))-2,1)) = Jaco(X2(3*(X1_sorted(Esi(i,3),2))-2,1), X2(3*(X1_sorted(Esi(i,1),2))-2,1))...
                    + esi*edgeEsi(i,3)/lenEsi(i,3);
                Res(X2(3*(X1_sorted(Esi(i,3),2))-2,1),1) =  Res(X2(3*(X1_sorted(Esi(i,3),2))-2,1),1)...
                    + (edgeEsi(i,2)*lenEsi(i,2) + edgeEsi(i,3)*lenEsi(i,3))/4*q/e0*Na...
                    - esi*(edgeEsi(i,2)/lenEsi(i,2) + edgeEsi(i,3)/lenEsi(i,3))*phi(X1_sorted(Esi(i,3),2))...
                    + esi*edgeEsi(i,2)/lenEsi(i,2)*phi(X1_sorted(Esi(i,2),2))...
                    + esi*edgeEsi(i,3)/lenEsi(i,3)*phi(X1_sorted(Esi(i,1),2));
                saveexl(X1_sorted(Esi(i,3),2)) = saveexl(X1_sorted(Esi(i,3),2)) + (edgeEsi(i,2)*lenEsi(i,2) + edgeEsi(i,3)*lenEsi(i,3))/4;

            elseif i > (E2row+E3row) && i <= Esirow
                Jaco(X2(3*(X1_sorted(Esi(i,1),2))-2,1), X2(3*(X1_sorted(Esi(i,1),2))-2,1)) = Jaco(X2(3*(X1_sorted(Esi(i,1),2))-2,1), X2(3*(X1_sorted(Esi(i,1),2))-2,1))...
                    - esi*(edgeEsi(i,1)/lenEsi(i,1) + edgeEsi(i,3)/lenEsi(i,3));
                Jaco(X2(3*(X1_sorted(Esi(i,1),2))-2,1), X2(3*(X1_sorted(Esi(i,2),2))-2,1)) = Jaco(X2(3*(X1_sorted(Esi(i,1),2))-2,1), X2(3*(X1_sorted(Esi(i,2),2))-2,1))...
                    + esi*edgeEsi(i,1)/lenEsi(i,1);
                Jaco(X2(3*(X1_sorted(Esi(i,1),2))-2,1), X2(3*(X1_sorted(Esi(i,3),2))-2,1)) = Jaco(X2(3*(X1_sorted(Esi(i,1),2))-2,1), X2(3*(X1_sorted(Esi(i,3),2))-2,1))...
                    + esi*edgeEsi(i,3)/lenEsi(i,3);
                Res(X2(3*(X1_sorted(Esi(i,1),2))-2,1),1) =  Res(X2(3*(X1_sorted(Esi(i,1),2))-2,1),1)...
                    + (edgeEsi(i,1)*lenEsi(i,1) + edgeEsi(i,3)*lenEsi(i,3))/4*q/e0*Nd...
                    - esi*(edgeEsi(i,1)/lenEsi(i,1) + edgeEsi(i,3)/lenEsi(i,3))*phi(X1_sorted(Esi(i,1),2))...
                    + esi*edgeEsi(i,1)/lenEsi(i,1)*phi(X1_sorted(Esi(i,2),2))...
                    + esi*edgeEsi(i,3)/lenEsi(i,3)*phi(X1_sorted(Esi(i,3),2));
                saveexl(X1_sorted(Esi(i,1),2)) = saveexl(X1_sorted(Esi(i,1),2)) + (edgeEsi(i,1)*lenEsi(i,1) + edgeEsi(i,3)*lenEsi(i,3))/4;

                Jaco(X2(3*(X1_sorted(Esi(i,2),2))-2,1), X2(3*(X1_sorted(Esi(i,2),2))-2,1)) = Jaco(X2(3*(X1_sorted(Esi(i,2),2))-2,1), X2(3*(X1_sorted(Esi(i,2),2))-2,1))...
                    - esi*(edgeEsi(i,1)/lenEsi(i,1) + edgeEsi(i,2)/lenEsi(i,2));
                Jaco(X2(3*(X1_sorted(Esi(i,2),2))-2,1), X2(3*(X1_sorted(Esi(i,1),2))-2,1)) = Jaco(X2(3*(X1_sorted(Esi(i,2),2))-2,1), X2(3*(X1_sorted(Esi(i,1),2))-2,1))...
                    + esi*edgeEsi(i,1)/lenEsi(i,1);
                Jaco(X2(3*(X1_sorted(Esi(i,2),2))-2,1), X2(3*(X1_sorted(Esi(i,3),2))-2,1)) = Jaco(X2(3*(X1_sorted(Esi(i,2),2))-2,1), X2(3*(X1_sorted(Esi(i,3),2))-2,1))...
                    + esi*edgeEsi(i,2)/lenEsi(i,2);
                Res(X2(3*(X1_sorted(Esi(i,2),2))-2,1),1) =  Res(X2(3*(X1_sorted(Esi(i,2),2))-2,1),1)...
                    + (edgeEsi(i,1)*lenEsi(i,1) + edgeEsi(i,2)*lenEsi(i,2))/4*q/e0*Nd...
                    - esi*(edgeEsi(i,1)/lenEsi(i,1) + edgeEsi(i,2)/lenEsi(i,2))*phi(X1_sorted(Esi(i,2),2))...
                    + esi*edgeEsi(i,1)/lenEsi(i,1)*phi(X1_sorted(Esi(i,1),2))...
                    + esi*edgeEsi(i,2)/lenEsi(i,2)*phi(X1_sorted(Esi(i,3),2));
                saveexl(X1_sorted(Esi(i,2),2)) = saveexl(X1_sorted(Esi(i,2),2)) + (edgeEsi(i,1)*lenEsi(i,1) + edgeEsi(i,2)*lenEsi(i,2))/4;

                Jaco(X2(3*(X1_sorted(Esi(i,3),2))-2,1), X2(3*(X1_sorted(Esi(i,3),2))-2,1)) = Jaco(X2(3*(X1_sorted(Esi(i,3),2))-2,1), X2(3*(X1_sorted(Esi(i,3),2))-2,1))...
                    - esi*(edgeEsi(i,2)/lenEsi(i,2) + edgeEsi(i,3)/lenEsi(i,3));
                Jaco(X2(3*(X1_sorted(Esi(i,3),2))-2,1), X2(3*(X1_sorted(Esi(i,2),2))-2,1)) = Jaco(X2(3*(X1_sorted(Esi(i,3),2))-2,1), X2(3*(X1_sorted(Esi(i,2),2))-2,1))...
                    + esi*edgeEsi(i,2)/lenEsi(i,2);
                Jaco(X2(3*(X1_sorted(Esi(i,3),2))-2,1), X2(3*(X1_sorted(Esi(i,1),2))-2,1)) = Jaco(X2(3*(X1_sorted(Esi(i,3),2))-2,1), X2(3*(X1_sorted(Esi(i,1),2))-2,1))...
                    + esi*edgeEsi(i,3)/lenEsi(i,3);
                Res(X2(3*(X1_sorted(Esi(i,3),2))-2,1),1) =  Res(X2(3*(X1_sorted(Esi(i,3),2))-2,1),1)...
                    + (edgeEsi(i,2)*lenEsi(i,2) + edgeEsi(i,3)*lenEsi(i,3))/4*q/e0*Nd...
                    - esi*(edgeEsi(i,2)/lenEsi(i,2) + edgeEsi(i,3)/lenEsi(i,3))*phi(X1_sorted(Esi(i,3),2))...
                    + esi*edgeEsi(i,2)/lenEsi(i,2)*phi(X1_sorted(Esi(i,2),2))...
                    + esi*edgeEsi(i,3)/lenEsi(i,3)*phi(X1_sorted(Esi(i,1),2));
                saveexl(X1_sorted(Esi(i,3),2)) = saveexl(X1_sorted(Esi(i,3),2)) + (edgeEsi(i,2)*lenEsi(i,2) + edgeEsi(i,3)*lenEsi(i,3))/4;
            end
        end

        for i=1:E5row
            Jaco(X2(3*(X1_sorted(E5(i,1),3))-2,1), X2(3*(X1_sorted(E5(i,1),3))-2,1)) = Jaco(X2(3*(X1_sorted(E5(i,1),3))-2,1), X2(3*(X1_sorted(E5(i,1),3))-2,1))...
                - eox*(edgeE5(i,1)/lenE5(i,1) + edgeE5(i,3)/lenE5(i,3));
            Jaco(X2(3*(X1_sorted(E5(i,1),3))-2,1), X2(3*(X1_sorted(E5(i,2),3))-2,1)) = Jaco(X2(3*(X1_sorted(E5(i,1),3))-2,1), X2(3*(X1_sorted(E5(i,2),3))-2,1))...
                + eox*edgeE5(i,1)/lenE5(i,1);
            Jaco(X2(3*(X1_sorted(E5(i,1),3))-2,1), X2(3*(X1_sorted(E5(i,3),3))-2,1)) = Jaco(X2(3*(X1_sorted(E5(i,1),3))-2,1), X2(3*(X1_sorted(E5(i,3),3))-2,1))...
                + eox*edgeE5(i,3)/lenE5(i,3);
            Res(X2(3*(X1_sorted(E5(i,1),3))-2,1),1) = Res(X2(3*(X1_sorted(E5(i,1),3))-2,1),1)...
                - eox*(edgeE5(i,1)/lenE5(i,1) + edgeE5(i,3)/lenE5(i,3))*phi(X1_sorted(E5(i,1),3))...
                + eox*edgeE5(i,1)/lenE5(i,1)*phi(X1_sorted(E5(i,2),3))...
                + eox*edgeE5(i,3)/lenE5(i,3)*phi(X1_sorted(E5(i,3),3));

            Jaco(X2(3*(X1_sorted(E5(i,2),3))-2,1), X2(3*(X1_sorted(E5(i,2),3))-2,1)) = Jaco(X2(3*(X1_sorted(E5(i,2),3))-2,1), X2(3*(X1_sorted(E5(i,2),3))-2,1))...
                - eox*(edgeE5(i,1)/lenE5(i,1) + edgeE5(i,2)/lenE5(i,2));
            Jaco(X2(3*(X1_sorted(E5(i,2),3))-2,1), X2(3*(X1_sorted(E5(i,1),3))-2,1)) = Jaco(X2(3*(X1_sorted(E5(i,2),3))-2,1), X2(3*(X1_sorted(E5(i,1),3))-2,1))...
                + eox*edgeE5(i,1)/lenE5(i,1);
            Jaco(X2(3*(X1_sorted(E5(i,2),3))-2,1), X2(3*(X1_sorted(E5(i,3),3))-2,1)) = Jaco(X2(3*(X1_sorted(E5(i,2),3))-2,1), X2(3*(X1_sorted(E5(i,3),3))-2,1))...
                + eox*edgeE5(i,2)/lenE5(i,2);
            Res(X2(3*(X1_sorted(E5(i,2),3))-2,1),1) = Res(X2(3*(X1_sorted(E5(i,2),3))-2,1),1)...
                - eox*(edgeE5(i,1)/lenE5(i,1) + edgeE5(i,2)/lenE5(i,2))*phi(X1_sorted(E5(i,2),3))...
                + eox*edgeE5(i,1)/lenE5(i,1)*phi(X1_sorted(E5(i,1),3))...
                + eox*edgeE5(i,2)/lenE5(i,2)*phi(X1_sorted(E5(i,3),3));

            Jaco(X2(3*(X1_sorted(E5(i,3),3))-2,1), X2(3*(X1_sorted(E5(i,3),3))-2,1)) = Jaco(X2(3*(X1_sorted(E5(i,3),3))-2,1), X2(3*(X1_sorted(E5(i,3),3))-2,1))...
                - eox*(edgeE5(i,2)/lenE5(i,2) + edgeE5(i,3)/lenE5(i,3));
            Jaco(X2(3*(X1_sorted(E5(i,3),3))-2,1), X2(3*(X1_sorted(E5(i,2),3))-2,1)) = Jaco(X2(3*(X1_sorted(E5(i,3),3))-2,1), X2(3*(X1_sorted(E5(i,2),3))-2,1))...
                + eox*edgeE5(i,2)/lenE5(i,2);
            Jaco(X2(3*(X1_sorted(E5(i,3),3))-2,1), X2(3*(X1_sorted(E5(i,1),3))-2,1)) = Jaco(X2(3*(X1_sorted(E5(i,3),3))-2,1), X2(3*(X1_sorted(E5(i,1),3))-2,1))...
                + eox*edgeE5(i,3)/lenE5(i,3);
            Res(X2(3*(X1_sorted(E5(i,3),3))-2,1),1) = Res(X2(3*(X1_sorted(E5(i,3),3))-2,1),1)...
                - eox*(edgeE5(i,2)/lenE5(i,2) + edgeE5(i,3)/lenE5(i,3))*phi(X1_sorted(E5(i,3),3))...
                + eox*edgeE5(i,2)/lenE5(i,2)*phi(X1_sorted(E5(i,2),3))...
                + eox*edgeE5(i,3)/lenE5(i,3)*phi(X1_sorted(E5(i,1),3));
        end

        %%%  for electron
        for i=1:V_silicon_row
            %%% - qn
            Jaco(X2(3*V_silicon(i,2)-2),X2(3*V_silicon(i,2)-1)) = Jaco(X2(3*V_silicon(i,2)-2),X2(3*V_silicon(i,2)-1)) - saveexl(V_silicon(i,2))*q/e0;
            Res(X2(3*V_silicon(i,2)-2),1) = Res(X2(3*V_silicon(i,2)-2),1) - saveexl(V_silicon(i,2))*q/e0*n(i);
            %%% + qp
            Jaco(X2(3*V_silicon(i,2)-2),X2(3*V_silicon(i,2))) =  Jaco(X2(3*V_silicon(i,2)-2),X2(3*V_silicon(i,2))) + saveexl(V_silicon(i,2))*q/e0;
            Res(X2(3*V_silicon(i,2)-2),1) = Res(X2(3*V_silicon(i,2)-2),1) + saveexl(V_silicon(i,2))*q/e0*p(i);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%  for electron DD
        for i=1:Esirow

            f1 = find(V_silicon(:,1)==Esi(i,1));
            f2 = find(V_silicon(:,1)==Esi(i,2));
            f3 = find(V_silicon(:,1)==Esi(i,3));

            ncoeff = q*mun*Vt;

            %%% 1st element
            % n part
            Jaco(X2(3*X1_sorted(Esi(i,1),2)-1,1), X2(3*X1_sorted(Esi(i,1),2)-1,1)) = Jaco(X2(3*X1_sorted(Esi(i,1),2)-1,1), X2(3*X1_sorted(Esi(i,1),2)-1,1))...
                - ncoeff*edgeEsi(i,1)/lenEsi(i,1)*bern(phi(X1_sorted(Esi(i,2),2)),phi(X1_sorted(Esi(i,1),2)))...
                - ncoeff*edgeEsi(i,3)/lenEsi(i,3)*bern(phi(X1_sorted(Esi(i,3),2)),phi(X1_sorted(Esi(i,1),2)));
            Jaco(X2(3*X1_sorted(Esi(i,1),2)-1,1), X2(3*X1_sorted(Esi(i,2),2)-1,1)) = Jaco(X2(3*X1_sorted(Esi(i,1),2)-1,1), X2(3*X1_sorted(Esi(i,2),2)-1,1))...
                + ncoeff*edgeEsi(i,1)/lenEsi(i,1)*bern(phi(X1_sorted(Esi(i,1),2)),phi(X1_sorted(Esi(i,2),2)));
            Jaco(X2(3*X1_sorted(Esi(i,1),2)-1,1), X2(3*X1_sorted(Esi(i,3),2)-1,1)) = Jaco(X2(3*X1_sorted(Esi(i,1),2)-1,1), X2(3*X1_sorted(Esi(i,3),2)-1,1))...
                + ncoeff*edgeEsi(i,3)/lenEsi(i,3)*bern(phi(X1_sorted(Esi(i,1),2)),phi(X1_sorted(Esi(i,3),2)));

            % phi part
            Jaco(X2(3*X1_sorted(Esi(i,1),2)-1,1), X2(3*X1_sorted(Esi(i,1),2)-2,1)) = Jaco(X2(3*X1_sorted(Esi(i,1),2)-1,1), X2(3*X1_sorted(Esi(i,1),2)-2,1))...
                - ncoeff/Vt*edgeEsi(i,1)/lenEsi(i,1)*dbern(phi(X1_sorted(Esi(i,1),2)), phi(X1_sorted(Esi(i,2),2)))*n(f2)...
                - ncoeff/Vt*edgeEsi(i,1)/lenEsi(i,1)*dbern(phi(X1_sorted(Esi(i,2),2)), phi(X1_sorted(Esi(i,1),2)))*n(f1)...
                - ncoeff/Vt*edgeEsi(i,3)/lenEsi(i,3)*dbern(phi(X1_sorted(Esi(i,1),2)), phi(X1_sorted(Esi(i,3),2)))*n(f3)...
                - ncoeff/Vt*edgeEsi(i,3)/lenEsi(i,3)*dbern(phi(X1_sorted(Esi(i,3),2)), phi(X1_sorted(Esi(i,1),2)))*n(f1);
            Jaco(X2(3*X1_sorted(Esi(i,1),2)-1,1), X2(3*X1_sorted(Esi(i,2),2)-2,1)) = Jaco(X2(3*X1_sorted(Esi(i,1),2)-1,1), X2(3*X1_sorted(Esi(i,2),2)-2,1))...
                + ncoeff/Vt*edgeEsi(i,1)/lenEsi(i,1)*dbern(phi(X1_sorted(Esi(i,1),2)), phi(X1_sorted(Esi(i,2),2)))*n(f2)...
                + ncoeff/Vt*edgeEsi(i,1)/lenEsi(i,1)*dbern(phi(X1_sorted(Esi(i,2),2)), phi(X1_sorted(Esi(i,1),2)))*n(f1);
            Jaco(X2(3*X1_sorted(Esi(i,1),2)-1,1), X2(3*X1_sorted(Esi(i,3),2)-2,1)) = Jaco(X2(3*X1_sorted(Esi(i,1),2)-1,1), X2(3*X1_sorted(Esi(i,3),2)-2,1))...
                + ncoeff/Vt*edgeEsi(i,3)/lenEsi(i,3)*dbern(phi(X1_sorted(Esi(i,1),2)), phi(X1_sorted(Esi(i,3),2)))*n(f3)...
                + ncoeff/Vt*edgeEsi(i,3)/lenEsi(i,3)*dbern(phi(X1_sorted(Esi(i,3),2)), phi(X1_sorted(Esi(i,1),2)))*n(f1);

            Res(X2(3*X1_sorted(Esi(i,1),2)-1,1),1) = Res(X2(3*X1_sorted(Esi(i,1),2)-1,1),1)...
                + ncoeff*edgeEsi(i,1)/lenEsi(i,1)*bern(phi(X1_sorted(Esi(i,1),2)),phi(X1_sorted(Esi(i,2),2)))*n(f2)...
                - ncoeff*edgeEsi(i,1)/lenEsi(i,1)*bern(phi(X1_sorted(Esi(i,2),2)),phi(X1_sorted(Esi(i,1),2)))*n(f1)...
                + ncoeff*edgeEsi(i,3)/lenEsi(i,3)*bern(phi(X1_sorted(Esi(i,1),2)),phi(X1_sorted(Esi(i,3),2)))*n(f3)...
                - ncoeff*edgeEsi(i,3)/lenEsi(i,3)*bern(phi(X1_sorted(Esi(i,3),2)),phi(X1_sorted(Esi(i,1),2)))*n(f1);

            %%% 2nd element
            % n part
            Jaco(X2(3*X1_sorted(Esi(i,2),2)-1,1), X2(3*X1_sorted(Esi(i,2),2)-1,1)) = Jaco(X2(3*X1_sorted(Esi(i,2),2)-1,1), X2(3*X1_sorted(Esi(i,2),2)-1,1))...
                - ncoeff*edgeEsi(i,2)/lenEsi(i,2)*bern(phi(X1_sorted(Esi(i,3),2)), phi(X1_sorted(Esi(i,2),2)))...
                - ncoeff*edgeEsi(i,1)/lenEsi(i,1)*bern(phi(X1_sorted(Esi(i,1),2)), phi(X1_sorted(Esi(i,2),2)));
            Jaco(X2(3*X1_sorted(Esi(i,2),2)-1,1), X2(3*X1_sorted(Esi(i,3),2)-1,1)) = Jaco(X2(3*X1_sorted(Esi(i,2),2)-1,1), X2(3*X1_sorted(Esi(i,3),2)-1,1))...
                + ncoeff*edgeEsi(i,2)/lenEsi(i,2)*bern(phi(X1_sorted(Esi(i,2),2)), phi(X1_sorted(Esi(i,3),2)));
            Jaco(X2(3*X1_sorted(Esi(i,2),2)-1,1), X2(3*X1_sorted(Esi(i,1),2)-1,1)) = Jaco(X2(3*X1_sorted(Esi(i,2),2)-1,1), X2(3*X1_sorted(Esi(i,1),2)-1,1))...
                + ncoeff*edgeEsi(i,1)/lenEsi(i,1)*bern(phi(X1_sorted(Esi(i,2),2)), phi(X1_sorted(Esi(i,1),2)));

            % phi part
            Jaco(X2(3*X1_sorted(Esi(i,2),2)-1,1), X2(3*X1_sorted(Esi(i,2),2)-2,1)) = Jaco(X2(3*X1_sorted(Esi(i,2),2)-1,1), X2(3*X1_sorted(Esi(i,2),2)-2,1))...
                - ncoeff/Vt*edgeEsi(i,2)/lenEsi(i,2)*dbern(phi(X1_sorted(Esi(i,2),2)), phi(X1_sorted(Esi(i,3),2)))*n(f3)...
                - ncoeff/Vt*edgeEsi(i,2)/lenEsi(i,2)*dbern(phi(X1_sorted(Esi(i,3),2)), phi(X1_sorted(Esi(i,2),2)))*n(f2)...
                - ncoeff/Vt*edgeEsi(i,1)/lenEsi(i,1)*dbern(phi(X1_sorted(Esi(i,2),2)), phi(X1_sorted(Esi(i,1),2)))*n(f1)...
                - ncoeff/Vt*edgeEsi(i,1)/lenEsi(i,1)*dbern(phi(X1_sorted(Esi(i,1),2)), phi(X1_sorted(Esi(i,2),2)))*n(f2);
            Jaco(X2(3*X1_sorted(Esi(i,2),2)-1,1), X2(3*X1_sorted(Esi(i,3),2)-2,1)) = Jaco(X2(3*X1_sorted(Esi(i,2),2)-1,1), X2(3*X1_sorted(Esi(i,3),2)-2,1))...
                + ncoeff/Vt*edgeEsi(i,2)/lenEsi(i,2)*dbern(phi(X1_sorted(Esi(i,2),2)), phi(X1_sorted(Esi(i,3),2)))*n(f3)...
                + ncoeff/Vt*edgeEsi(i,2)/lenEsi(i,2)*dbern(phi(X1_sorted(Esi(i,3),2)), phi(X1_sorted(Esi(i,2),2)))*n(f2);
            Jaco(X2(3*X1_sorted(Esi(i,2),2)-1,1), X2(3*X1_sorted(Esi(i,1),2)-2,1)) = Jaco(X2(3*X1_sorted(Esi(i,2),2)-1,1), X2(3*X1_sorted(Esi(i,1),2)-2,1))...
                + ncoeff/Vt*edgeEsi(i,1)/lenEsi(i,1)*dbern(phi(X1_sorted(Esi(i,2),2)), phi(X1_sorted(Esi(i,1),2)))*n(f1)...
                + ncoeff/Vt*edgeEsi(i,1)/lenEsi(i,1)*dbern(phi(X1_sorted(Esi(i,1),2)), phi(X1_sorted(Esi(i,2),2)))*n(f2);

            Res(X2(3*X1_sorted(Esi(i,2),2)-1,1),1) = Res(X2(3*X1_sorted(Esi(i,2),2)-1,1),1)...
                + ncoeff*edgeEsi(i,2)/lenEsi(i,2)*bern(phi(X1_sorted(Esi(i,2),2)),phi(X1_sorted(Esi(i,3),2)))*n(f3)...
                - ncoeff*edgeEsi(i,2)/lenEsi(i,2)*bern(phi(X1_sorted(Esi(i,3),2)),phi(X1_sorted(Esi(i,2),2)))*n(f2)...
                + ncoeff*edgeEsi(i,1)/lenEsi(i,1)*bern(phi(X1_sorted(Esi(i,2),2)),phi(X1_sorted(Esi(i,1),2)))*n(f1)...
                - ncoeff*edgeEsi(i,1)/lenEsi(i,1)*bern(phi(X1_sorted(Esi(i,1),2)),phi(X1_sorted(Esi(i,2),2)))*n(f2);

            %%% 3rd element
            % n part
            Jaco(X2(3*X1_sorted(Esi(i,3),2)-1,1), X2(3*X1_sorted(Esi(i,3),2)-1,1)) = Jaco(X2(3*X1_sorted(Esi(i,3),2)-1,1), X2(3*X1_sorted(Esi(i,3),2)-1,1))...
                - ncoeff*edgeEsi(i,3)/lenEsi(i,3)*bern(phi(X1_sorted(Esi(i,1),2)),phi(X1_sorted(Esi(i,3),2)))...
                - ncoeff*edgeEsi(i,2)/lenEsi(i,2)*bern(phi(X1_sorted(Esi(i,2),2)),phi(X1_sorted(Esi(i,3),2)));
            Jaco(X2(3*X1_sorted(Esi(i,3),2)-1,1), X2(3*X1_sorted(Esi(i,1),2)-1,1)) = Jaco(X2(3*X1_sorted(Esi(i,3),2)-1,1), X2(3*X1_sorted(Esi(i,1),2)-1,1))...
                + ncoeff*edgeEsi(i,3)/lenEsi(i,3)*bern(phi(X1_sorted(Esi(i,3),2)),phi(X1_sorted(Esi(i,1),2)));
            Jaco(X2(3*X1_sorted(Esi(i,3),2)-1,1), X2(3*X1_sorted(Esi(i,2),2)-1,1)) = Jaco(X2(3*X1_sorted(Esi(i,3),2)-1,1), X2(3*X1_sorted(Esi(i,2),2)-1,1))...
                + ncoeff*edgeEsi(i,2)/lenEsi(i,2)*bern(phi(X1_sorted(Esi(i,3),2)),phi(X1_sorted(Esi(i,2),2)));

            % phi part
            Jaco(X2(3*X1_sorted(Esi(i,3),2)-1,1), X2(3*X1_sorted(Esi(i,3),2)-2,1)) = Jaco(X2(3*X1_sorted(Esi(i,3),2)-1,1), X2(3*X1_sorted(Esi(i,3),2)-2,1))...
                - ncoeff/Vt*edgeEsi(i,3)/lenEsi(i,3)*dbern(phi(X1_sorted(Esi(i,3),2)), phi(X1_sorted(Esi(i,1),2)))*n(f1)...
                - ncoeff/Vt*edgeEsi(i,3)/lenEsi(i,3)*dbern(phi(X1_sorted(Esi(i,1),2)), phi(X1_sorted(Esi(i,3),2)))*n(f3)...
                - ncoeff/Vt*edgeEsi(i,2)/lenEsi(i,2)*dbern(phi(X1_sorted(Esi(i,3),2)), phi(X1_sorted(Esi(i,2),2)))*n(f2)...
                - ncoeff/Vt*edgeEsi(i,2)/lenEsi(i,2)*dbern(phi(X1_sorted(Esi(i,2),2)), phi(X1_sorted(Esi(i,3),2)))*n(f3);
            Jaco(X2(3*X1_sorted(Esi(i,3),2)-1,1), X2(3*X1_sorted(Esi(i,1),2)-2,1)) = Jaco(X2(3*X1_sorted(Esi(i,3),2)-1,1), X2(3*X1_sorted(Esi(i,1),2)-2,1))...
                + ncoeff/Vt*edgeEsi(i,3)/lenEsi(i,3)*dbern(phi(X1_sorted(Esi(i,3),2)), phi(X1_sorted(Esi(i,1),2)))*n(f1)...
                + ncoeff/Vt*edgeEsi(i,3)/lenEsi(i,3)*dbern(phi(X1_sorted(Esi(i,1),2)), phi(X1_sorted(Esi(i,3),2)))*n(f3);
            Jaco(X2(3*X1_sorted(Esi(i,3),2)-1,1), X2(3*X1_sorted(Esi(i,2),2)-2,1)) = Jaco(X2(3*X1_sorted(Esi(i,3),2)-1,1), X2(3*X1_sorted(Esi(i,2),2)-2,1))...
                + ncoeff/Vt*edgeEsi(i,2)/lenEsi(i,2)*dbern(phi(X1_sorted(Esi(i,3),2)), phi(X1_sorted(Esi(i,2),2)))*n(f2)...
                + ncoeff/Vt*edgeEsi(i,2)/lenEsi(i,2)*dbern(phi(X1_sorted(Esi(i,2),2)), phi(X1_sorted(Esi(i,3),2)))*n(f3);

            Res(X2(3*X1_sorted(Esi(i,3),2)-1,1),1) = Res(X2(3*X1_sorted(Esi(i,3),2)-1,1),1)...
                + ncoeff*edgeEsi(i,3)/lenEsi(i,3)*bern(phi(X1_sorted(Esi(i,3),2)),phi(X1_sorted(Esi(i,1),2)))*n(f1)...
                - ncoeff*edgeEsi(i,3)/lenEsi(i,3)*bern(phi(X1_sorted(Esi(i,1),2)),phi(X1_sorted(Esi(i,3),2)))*n(f3)...
                + ncoeff*edgeEsi(i,2)/lenEsi(i,2)*bern(phi(X1_sorted(Esi(i,3),2)),phi(X1_sorted(Esi(i,2),2)))*n(f2)...
                - ncoeff*edgeEsi(i,2)/lenEsi(i,2)*bern(phi(X1_sorted(Esi(i,2),2)),phi(X1_sorted(Esi(i,3),2)))*n(f3);
        end

        %%% for hole DD
        for i=1:Esirow

            f1 = find(V_silicon(:,1)==Esi(i,1));
            f2 = find(V_silicon(:,1)==Esi(i,2));
            f3 = find(V_silicon(:,1)==Esi(i,3));

            pcoeff = q*mup/Vt;

            %%% 1st element
            % p part
            Jaco(X2(3*X1_sorted(Esi(i,1),2),1), X2(3*X1_sorted(Esi(i,1),2),1)) = Jaco(X2(3*X1_sorted(Esi(i,1),2),1), X2(3*X1_sorted(Esi(i,1),2),1))...
                + pcoeff*edgeEsi(i,1)/lenEsi(i,1)*bern(phi(X1_sorted(Esi(i,1),2)),phi(X1_sorted(Esi(i,2),2)))...
                + pcoeff*edgeEsi(i,3)/lenEsi(i,3)*bern(phi(X1_sorted(Esi(i,1),2)),phi(X1_sorted(Esi(i,3),2)));
            Jaco(X2(3*X1_sorted(Esi(i,1),2),1), X2(3*X1_sorted(Esi(i,2),2),1)) = Jaco(X2(3*X1_sorted(Esi(i,1),2),1), X2(3*X1_sorted(Esi(i,2),2),1))...
                - pcoeff*edgeEsi(i,1)/lenEsi(i,1)*bern(phi(X1_sorted(Esi(i,2),2)),phi(X1_sorted(Esi(i,1),2)));
            Jaco(X2(3*X1_sorted(Esi(i,1),2),1), X2(3*X1_sorted(Esi(i,3),2),1)) = Jaco(X2(3*X1_sorted(Esi(i,1),2),1), X2(3*X1_sorted(Esi(i,3),2),1))...
                - pcoeff*edgeEsi(i,3)/lenEsi(i,3)*bern(phi(X1_sorted(Esi(i,3),2)),phi(X1_sorted(Esi(i,1),2)));

            % phi part
            Jaco(X2(3*X1_sorted(Esi(i,1),2),1), X2(3*X1_sorted(Esi(i,1),2)-2,1)) = Jaco(X2(3*X1_sorted(Esi(i,1),2),1), X2(3*X1_sorted(Esi(i,1),2)-2,1))...
                - pcoeff/Vt*edgeEsi(i,1)/lenEsi(i,1)*dbern(phi(X1_sorted(Esi(i,2),2)),phi(X1_sorted(Esi(i,1),2)))*p(f2)...
                - pcoeff/Vt*edgeEsi(i,1)/lenEsi(i,1)*dbern(phi(X1_sorted(Esi(i,1),2)),phi(X1_sorted(Esi(i,2),2)))*p(f1)...
                - pcoeff/Vt*edgeEsi(i,3)/lenEsi(i,3)*dbern(phi(X1_sorted(Esi(i,3),2)),phi(X1_sorted(Esi(i,1),2)))*p(f3)...
                - pcoeff/Vt*edgeEsi(i,3)/lenEsi(i,3)*dbern(phi(X1_sorted(Esi(i,1),2)),phi(X1_sorted(Esi(i,3),2)))*p(f1);
            Jaco(X2(3*X1_sorted(Esi(i,1),2),1), X2(3*X1_sorted(Esi(i,2),2)-2,1)) = Jaco(X2(3*X1_sorted(Esi(i,1),2),1), X2(3*X1_sorted(Esi(i,2),2)-2,1))...
                + pcoeff/Vt*edgeEsi(i,1)/lenEsi(i,1)*dbern(phi(X1_sorted(Esi(i,2),2)),phi(X1_sorted(Esi(i,1),2)))*p(f2)...
                + pcoeff/Vt*edgeEsi(i,1)/lenEsi(i,1)*dbern(phi(X1_sorted(Esi(i,1),2)),phi(X1_sorted(Esi(i,2),2)))*p(f1);
            Jaco(X2(3*X1_sorted(Esi(i,1),2),1), X2(3*X1_sorted(Esi(i,3),2)-2,1)) = Jaco(X2(3*X1_sorted(Esi(i,1),2),1), X2(3*X1_sorted(Esi(i,3),2)-2,1))...
                + pcoeff/Vt*edgeEsi(i,3)/lenEsi(i,3)*dbern(phi(X1_sorted(Esi(i,3),2)),phi(X1_sorted(Esi(i,1),2)))*p(f3)...
                + pcoeff/Vt*edgeEsi(i,3)/lenEsi(i,3)*dbern(phi(X1_sorted(Esi(i,1),2)),phi(X1_sorted(Esi(i,3),2)))*p(f1);

            Res(X2(3*X1_sorted(Esi(i,1),2),1),1) = Res(X2(3*X1_sorted(Esi(i,1),2),1),1)...
                - pcoeff*edgeEsi(i,1)/lenEsi(i,1)*bern(phi(X1_sorted(Esi(i,2),2)),phi(X1_sorted(Esi(i,1),2)))*p(f2)...
                + pcoeff*edgeEsi(i,1)/lenEsi(i,1)*bern(phi(X1_sorted(Esi(i,1),2)),phi(X1_sorted(Esi(i,2),2)))*p(f1)...
                - pcoeff*edgeEsi(i,3)/lenEsi(i,3)*bern(phi(X1_sorted(Esi(i,3),2)),phi(X1_sorted(Esi(i,1),2)))*p(f3)...
                + pcoeff*edgeEsi(i,3)/lenEsi(i,3)*bern(phi(X1_sorted(Esi(i,1),2)),phi(X1_sorted(Esi(i,3),2)))*p(f1);

            %%% 2nd element
            % p part
            Jaco(X2(3*X1_sorted(Esi(i,2),2),1), X2(3*X1_sorted(Esi(i,2),2),1)) = Jaco(X2(3*X1_sorted(Esi(i,2),2),1), X2(3*X1_sorted(Esi(i,2),2),1))...
                + pcoeff*edgeEsi(i,2)/lenEsi(i,2)*bern(phi(X1_sorted(Esi(i,2),2)),phi(X1_sorted(Esi(i,3),2)))...
                + pcoeff*edgeEsi(i,1)/lenEsi(i,1)*bern(phi(X1_sorted(Esi(i,2),2)),phi(X1_sorted(Esi(i,1),2)));
            Jaco(X2(3*X1_sorted(Esi(i,2),2),1), X2(3*X1_sorted(Esi(i,3),2),1)) = Jaco(X2(3*X1_sorted(Esi(i,2),2),1), X2(3*X1_sorted(Esi(i,3),2),1))...
                - pcoeff*edgeEsi(i,2)/lenEsi(i,2)*bern(phi(X1_sorted(Esi(i,3),2)),phi(X1_sorted(Esi(i,2),2)));
            Jaco(X2(3*X1_sorted(Esi(i,2),2),1), X2(3*X1_sorted(Esi(i,1),2),1)) = Jaco(X2(3*X1_sorted(Esi(i,2),2),1), X2(3*X1_sorted(Esi(i,1),2),1))...
                - pcoeff*edgeEsi(i,1)/lenEsi(i,1)*bern(phi(X1_sorted(Esi(i,1),2)),phi(X1_sorted(Esi(i,2),2)));

            % phi part
            Jaco(X2(3*X1_sorted(Esi(i,2),2),1), X2(3*X1_sorted(Esi(i,2),2)-2,1)) = Jaco(X2(3*X1_sorted(Esi(i,2),2),1), X2(3*X1_sorted(Esi(i,2),2)-2,1))...
                - pcoeff/Vt*edgeEsi(i,2)/lenEsi(i,2)*dbern(phi(X1_sorted(Esi(i,3),2)),phi(X1_sorted(Esi(i,2),2)))*p(f3)...
                - pcoeff/Vt*edgeEsi(i,2)/lenEsi(i,2)*dbern(phi(X1_sorted(Esi(i,2),2)),phi(X1_sorted(Esi(i,3),2)))*p(f2)...
                - pcoeff/Vt*edgeEsi(i,1)/lenEsi(i,1)*dbern(phi(X1_sorted(Esi(i,1),2)),phi(X1_sorted(Esi(i,2),2)))*p(f1)...
                - pcoeff/Vt*edgeEsi(i,1)/lenEsi(i,1)*dbern(phi(X1_sorted(Esi(i,2),2)),phi(X1_sorted(Esi(i,1),2)))*p(f2);
            Jaco(X2(3*X1_sorted(Esi(i,2),2),1), X2(3*X1_sorted(Esi(i,3),2)-2,1)) = Jaco(X2(3*X1_sorted(Esi(i,2),2),1), X2(3*X1_sorted(Esi(i,3),2)-2,1))...
                + pcoeff/Vt*edgeEsi(i,2)/lenEsi(i,2)*dbern(phi(X1_sorted(Esi(i,3),2)),phi(X1_sorted(Esi(i,2),2)))*p(f3)...
                + pcoeff/Vt*edgeEsi(i,2)/lenEsi(i,2)*dbern(phi(X1_sorted(Esi(i,2),2)),phi(X1_sorted(Esi(i,3),2)))*p(f2);
            Jaco(X2(3*X1_sorted(Esi(i,2),2),1), X2(3*X1_sorted(Esi(i,1),2)-2,1)) = Jaco(X2(3*X1_sorted(Esi(i,2),2),1), X2(3*X1_sorted(Esi(i,1),2)-2,1))...
                + pcoeff/Vt*edgeEsi(i,1)/lenEsi(i,1)*dbern(phi(X1_sorted(Esi(i,1),2)),phi(X1_sorted(Esi(i,2),2)))*p(f1)...
                + pcoeff/Vt*edgeEsi(i,1)/lenEsi(i,1)*dbern(phi(X1_sorted(Esi(i,2),2)),phi(X1_sorted(Esi(i,1),2)))*p(f2);

            Res(X2(3*X1_sorted(Esi(i,2),2),1),1) = Res(X2(3*X1_sorted(Esi(i,2),2),1),1)...
                - pcoeff*edgeEsi(i,2)/lenEsi(i,2)*bern(phi(X1_sorted(Esi(i,3),2)),phi(X1_sorted(Esi(i,2),2)))*p(f3)...
                + pcoeff*edgeEsi(i,2)/lenEsi(i,2)*bern(phi(X1_sorted(Esi(i,2),2)),phi(X1_sorted(Esi(i,3),2)))*p(f2)...
                - pcoeff*edgeEsi(i,1)/lenEsi(i,1)*bern(phi(X1_sorted(Esi(i,1),2)),phi(X1_sorted(Esi(i,2),2)))*p(f1)...
                + pcoeff*edgeEsi(i,1)/lenEsi(i,1)*bern(phi(X1_sorted(Esi(i,2),2)),phi(X1_sorted(Esi(i,1),2)))*p(f2);

            %%% 3rd element
            % p part
            Jaco(X2(3*X1_sorted(Esi(i,3),2),1), X2(3*X1_sorted(Esi(i,3),2),1)) = Jaco(X2(3*X1_sorted(Esi(i,3),2),1), X2(3*X1_sorted(Esi(i,3),2),1))...
                + pcoeff*edgeEsi(i,3)/lenEsi(i,3)*bern(phi(X1_sorted(Esi(i,3),2)),phi(X1_sorted(Esi(i,1),2)))...
                + pcoeff*edgeEsi(i,2)/lenEsi(i,2)*bern(phi(X1_sorted(Esi(i,3),2)),phi(X1_sorted(Esi(i,2),2)));
            Jaco(X2(3*X1_sorted(Esi(i,3),2),1), X2(3*X1_sorted(Esi(i,1),2),1)) = Jaco(X2(3*X1_sorted(Esi(i,3),2),1), X2(3*X1_sorted(Esi(i,1),2),1))...
                - pcoeff*edgeEsi(i,3)/lenEsi(i,3)*bern(phi(X1_sorted(Esi(i,1),2)),phi(X1_sorted(Esi(i,3),2)));
            Jaco(X2(3*X1_sorted(Esi(i,3),2),1), X2(3*X1_sorted(Esi(i,2),2),1)) = Jaco(X2(3*X1_sorted(Esi(i,3),2),1), X2(3*X1_sorted(Esi(i,2),2),1))...
                - pcoeff*edgeEsi(i,2)/lenEsi(i,2)*bern(phi(X1_sorted(Esi(i,2),2)),phi(X1_sorted(Esi(i,3),2)));

            % phi part
            Jaco(X2(3*X1_sorted(Esi(i,3),2),1), X2(3*X1_sorted(Esi(i,3),2)-2,1)) = Jaco(X2(3*X1_sorted(Esi(i,3),2),1), X2(3*X1_sorted(Esi(i,3),2)-2,1))...
                - pcoeff/Vt*edgeEsi(i,3)/lenEsi(i,3)*dbern(phi(X1_sorted(Esi(i,1),2)),phi(X1_sorted(Esi(i,3),2)))*p(f1)...
                - pcoeff/Vt*edgeEsi(i,3)/lenEsi(i,3)*dbern(phi(X1_sorted(Esi(i,3),2)),phi(X1_sorted(Esi(i,1),2)))*p(f3)...
                - pcoeff/Vt*edgeEsi(i,2)/lenEsi(i,2)*dbern(phi(X1_sorted(Esi(i,2),2)),phi(X1_sorted(Esi(i,3),2)))*p(f2)...
                - pcoeff/Vt*edgeEsi(i,2)/lenEsi(i,2)*dbern(phi(X1_sorted(Esi(i,3),2)),phi(X1_sorted(Esi(i,2),2)))*p(f3);
            Jaco(X2(3*X1_sorted(Esi(i,3),2),1), X2(3*X1_sorted(Esi(i,1),2)-2,1)) = Jaco(X2(3*X1_sorted(Esi(i,3),2),1), X2(3*X1_sorted(Esi(i,1),2)-2,1))...
                + pcoeff/Vt*edgeEsi(i,3)/lenEsi(i,3)*dbern(phi(X1_sorted(Esi(i,1),2)),phi(X1_sorted(Esi(i,3),2)))*p(f1)...
                + pcoeff/Vt*edgeEsi(i,3)/lenEsi(i,3)*dbern(phi(X1_sorted(Esi(i,3),2)),phi(X1_sorted(Esi(i,1),2)))*p(f3);
            Jaco(X2(3*X1_sorted(Esi(i,3),2),1), X2(3*X1_sorted(Esi(i,2),2)-2,1)) = Jaco(X2(3*X1_sorted(Esi(i,3),2),1), X2(3*X1_sorted(Esi(i,2),2)-2,1))...
                + pcoeff/Vt*edgeEsi(i,2)/lenEsi(i,2)*dbern(phi(X1_sorted(Esi(i,2),2)),phi(X1_sorted(Esi(i,3),2)))*p(f2)...
                + pcoeff/Vt*edgeEsi(i,2)/lenEsi(i,2)*dbern(phi(X1_sorted(Esi(i,3),2)),phi(X1_sorted(Esi(i,2),2)))*p(f3);

            Res(X2(3*X1_sorted(Esi(i,3),2),1),1) = Res(X2(3*X1_sorted(Esi(i,3),2),1),1)...
                - pcoeff*edgeEsi(i,3)/lenEsi(i,3)*bern(phi(X1_sorted(Esi(i,1),2)),phi(X1_sorted(Esi(i,3),2)))*p(f1)...
                + pcoeff*edgeEsi(i,3)/lenEsi(i,3)*bern(phi(X1_sorted(Esi(i,3),2)),phi(X1_sorted(Esi(i,1),2)))*p(f3)...
                - pcoeff*edgeEsi(i,2)/lenEsi(i,2)*bern(phi(X1_sorted(Esi(i,2),2)),phi(X1_sorted(Esi(i,3),2)))*p(f2)...
                + pcoeff*edgeEsi(i,2)/lenEsi(i,2)*bern(phi(X1_sorted(Esi(i,3),2)),phi(X1_sorted(Esi(i,2),2)))*p(f3);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% continuity
        for i=1:is1row
            minus1 = X1_sorted(is1(i),1);
            plus1 = X1_sorted(is1(i),2);
            Jaco(X2(3*plus1-2),:) = Jaco(X2(3*plus1-2),:) + Jaco(X2(3*minus1-2),:);
            Res(X2(3*plus1-2)) = Res(X2(3*plus1-2)) + Res(X2(3*minus1-2));
            Jaco(X2(3*minus1-2),:) = 0;
            Res(X2(3*minus1-2),1) = 0;
            Jaco(X2(3*minus1-2), X2(3*minus1-2)) = -1;
            Jaco(X2(3*minus1-2), X2(3*plus1-2)) = 1;
        end

        for i=1:is2row
            minus1 = X1_sorted(is2(i),3);
            plus1 = X1_sorted(is2(i),2);
            Jaco(X2(3*plus1-2),:) = Jaco(X2(3*plus1-2),:) + Jaco(X2(3*minus1-2),:);
            Res(X2(3*plus1-2),1) = Res(X2(3*plus1-2),1) + Res(X2(3*minus1-2),1);
            Jaco(X2(3*minus1-2),:) = 0;
            Res(X2(3*minus1-2),1) = 0;
            Jaco(X2(3*minus1-2), X2(3*minus1-2)) = -1;
            Jaco(X2(3*minus1-2), X2(3*plus1-2)) = 1;
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Dirichlet boundary condition for final iteration
        % gate1
        for i=1:Ncon1
            for j=1:concol
                fcon = X1_sorted(contact(i,j),1);
                Jaco(X2(3*fcon-2),:) = 0;
                Jaco(X2(3*fcon-2), X2(3*fcon-2)) = 1;
                Res(X2(3*fcon-2),1) = phi(fcon,1) - 0.33374;
            end
        end

        % gate2
        for i=Ncon1+1:(Ncon1+Ncon2)
            for j=1:concol
                fcon = X1_sorted(contact(i,j),3);
                Jaco(X2(3*fcon-2),:) = 0;
                Jaco(X2(3*fcon-2), X2(3*fcon-2)) = 1;
                Res(X2(3*fcon-2),1) = phi(fcon) - 0.33374;
            end
        end

        % source & drain
        for i=(Ncon1+Ncon2)+1:(Ncon1+Ncon2+Ncon3+Ncon4)
            for j=1:concol
                fcon = X1_sorted(contact(i,j),2);
                fcon_np = find(contact(i,j)==V_silicon(:,1));

                Jaco(X2(3*fcon-2),:) = 0;
                Jaco(X2(3*fcon-1),:) = 0;
                Jaco(X2(3*fcon),:) = 0;

                Jaco(X2(3*fcon-2), X2(3*fcon-2)) = 1;
                Jaco(X2(3*fcon-1), X2(3*fcon-1)) = 1;
                Jaco(X2(3*fcon), X2(3*fcon)) = 1;

                % source
                if i<=(Ncon1+Ncon2+Ncon3)
                    Res(X2(3*fcon-2),1) = phi(fcon) - Vt*log(Nd/nint);
                    Res(X2(3*fcon-1),1) = n(fcon_np) - Nd;
                    Res(X2(3*fcon),1) = p(fcon_np) - nint*nint/Nd;

                    % drain
                elseif i>(Ncon1+Ncon2+Ncon3) &&  i<=(Ncon1+Ncon2+Ncon3+Ncon4)
                    Res(X2(3*fcon-2),1) = phi(fcon) - Vt*log(Nd/nint) - Vd;
                    Res(X2(3*fcon-1),1) = n(fcon_np) - Nd;
                    Res(X2(3*fcon),1) = p(fcon_np) - nint*nint/Nd;
                end
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% scaled part
        Cvector = zeros(size(Jaco,1),1);
        for i=1:size(Jaco,1)
            if rem(X2_sorted(i),3) == 1
                Cvector(i) = Vt;
            elseif rem(X2_sorted(i),3) == 2
                Cvector(i) = Nd;
            elseif rem(X2_sorted(i),3) == 0
                Cvector(i) = Nd;
            end
        end

        Cmatrix = spdiags(Cvector,0,size(X2_sorted,1),size(X2_sorted,1));
        Jaco_scaled = Jaco * Cmatrix;
        Rvector = 1./sum(abs(Jaco_scaled),2);
        Rmatrix = spdiags(Rvector,0,size(X2_sorted,1),size(X2_sorted,1));
        Jaco_scaled = Rmatrix * Jaco_scaled;
        Res_scaled = Rmatrix * Res;
        update_scaled = Jaco_scaled \ (-Res_scaled);
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

        save_drain_update(:,it) = update(:,1);
        save_drain_phi_update(it,step) = maxupdate;

        if maxupdate < 1e-10
            break
        end
    end

    save_drain_phi(:,step) = phi(:,1);
    save_drain_n(:,step) = n(:,1);
    save_drain_p(:,step) = p(:,1);

    for i = 1:size(unicond,1)
        for j = 1:E4row
            fcon = find(unicond(i) == E4(j,:));
            f1 = find(V_silicon(:,1) == E4(j,1));
            f2 = find(V_silicon(:,1) == E4(j,2));
            f3 = find(V_silicon(:,1) == E4(j,3));

            if fcon == 1
                I(step,1) = I(step,1) ...
                + ncoeff*edgeE4(j,2)/lenE4(j,3)*(n(f2)*bern(phi(X1_sorted(E4(j,1),2)),phi(X1_sorted(E4(j,2),2))) - n(f1)*bern(phi(X1_sorted(E4(j,2),2)),phi(X1_sorted(E4(j,1),2))))*width...
                + ncoeff*edgeE4(j,3)/lenE4(j,3)*(n(f3)*bern(phi(X1_sorted(E4(j,1),2)),phi(X1_sorted(E4(j,3),2))) - n(f1)*bern(phi(X1_sorted(E4(j,3),2)),phi(X1_sorted(E4(j,1),2))))*width;
  
            elseif fcon == 2
                I(step,1) = I(step,1) ...
                + ncoeff*edgeE4(j,2)/lenE4(j,2)*(n(f3)*bern(phi(X1_sorted(E4(j,2),2)),phi(X1_sorted(E4(j,3),2))) - n(f2)*bern(phi(X1_sorted(E4(j,3),2)),phi(X1_sorted(E4(j,2),2))))*width...
                + ncoeff*edgeE4(j,1)/lenE4(j,1)*(n(f1)*bern(phi(X1_sorted(E4(j,2),2)),phi(X1_sorted(E4(j,1),2))) - n(f2)*bern(phi(X1_sorted(E4(j,1),2)),phi(X1_sorted(E4(j,2),2))))*width;
            
            elseif fcon == 3
                I(step,1) = I(step,1) ...
                + ncoeff*edgeE4(j,3)/lenE4(j,1)*(n(f3)*bern(phi(X1_sorted(E4(j,3),2)),phi(X1_sorted(E4(j,1),2))) - n(f3)*bern(phi(X1_sorted(E4(j,1),2)),phi(X1_sorted(E4(j,3),2))))*width...
                + ncoeff*edgeE4(j,2)/lenE4(j,2)*(n(f2)*bern(phi(X1_sorted(E4(j,3),2)),phi(X1_sorted(E4(j,2),2))) - n(f3)*bern(phi(X1_sorted(E4(j,2),2)),phi(X1_sorted(E4(j,3),2))))*width;
            end
        end
    end

end

phi_Vd_1mV = phi(:,1);
