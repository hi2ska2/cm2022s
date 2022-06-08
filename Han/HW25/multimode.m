clc
clear
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Double-gate MOSFET 1st gate ramping %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
q = 1.6021e-19;
e0 = 8.8542e-12;
eox = 3.9;
esi = 11.7;
nint = 1e16;
Kb = 1.3803e-23;
Temp = 300;
Vt = Kb*Temp/q;
mun = 1417e-4;
mup = 470e-4;
ncoeff = q*mun;
pcoeff = q*mup;
Vdd = 0;
resistor = 1000;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
width = 1e-6;
Nterminal = 4;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Doping density
Nd = 5e26;
Na = -2e21;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% variable
Vtg = 0;
Vbg = 0;
Vs = 0;
Vd = 0;
Itg = 0;
Ibg = 0;
Is = 0;
Id = 0;
I1 = 0;
I2 = 0;
V1 = 0;
V2 = 0;
Vout = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Iteration number
Niter = 30;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
V = importdata("dg_vertex.txt"); Vrow = size(V,1);
E = importdata("dg_element.txt"); Erow = size(E,1);
E1 = importdata("dg_element_region1.txt"); E1row = size(E1,1);
E2 = importdata("dg_element_region2.txt"); E2row = size(E2,1);
E3 = importdata("dg_element_region3.txt"); E3row = size(E3,1);
E4 = importdata("dg_element_region4.txt"); E4row = size(E4,1);
E5 = importdata("dg_element_region5.txt"); E5row = size(E5,1);
Esi = vertcat(E2, E3, E4); Esirow = size(Esi,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% contact information
contact = importdata("dg_contact.txt");
[conrow, concol] = size(contact);
unicon = unique(contact);
Ncon1 = 20; % Top gate contact
Ncon2 = 20; % Bottom gate contact
Ncon3 = 12; % source contact
Ncon4 = 12; % drain cintact

% Top gate terminal analysis
contg = zeros(2*Ncon1,1);
ttt=1;
for i=1:Ncon1
    for j=1:concol
        contg(ttt,1) = contact(i,j);
        ttt = ttt+1;
    end
end
unicontg = unique(contg);

% Bottom gate terminal analysis
conbg = zeros(2*Ncon2,1);
bbb=1;
for i=Ncon1+1:Ncon1+Ncon2
    for j=1:concol
        conbg(bbb,1) = contact(i,j);
        bbb = bbb+1;
    end
end
uniconbg = unique(conbg);

% Source terminal analysis
cons = zeros(2*Ncon3,1);
sss=1;
for i=Ncon1+Ncon2+1:Ncon1+Ncon2+Ncon3
    for j=1:concol
        cons(sss,1) = contact(i,j);
        sss = sss+1;
    end
end
unicons = unique(cons);

% Drain terminal analysis
cond = zeros(2*Ncon4,1);
ddd=1;
for i=Ncon1+Ncon2+Ncon3+1:Ncon1+Ncon2+Ncon3+Ncon4
    for j=1:concol
        cond(ddd,1) = contact(i,j);
        ddd = ddd+1;
    end
end
unicond = unique(cond);

% variable information
Nvariable = 3;
Nregion = 5;

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

V_silicon(:,1) = unique(Esi); V_silicon_row = size(V_silicon, 1);
V_oxide = zeros(VE1row+VE5row, 2); V_oxide_row = size(V_oxide,1);

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
X1(:,1) = vertcat(VE1, V_silicon, VE5);
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

%%% Total Variable indexing
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

%%% Length
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

%%% Vector % Vector norm
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

%%% Area, Radius
area = zeros(Erow,1);
rad = zeros(Erow,1);

for i=1:Erow
    area(i,1) = sqrt(norm_v21(i,1)*norm_v21(i,1)*norm_v31(i,1)*norm_v31(i,1) - dot(v21(i,:),v31(i,:))*dot(v21(i,:),v31(i,:)))/2;
    rad(i,1) = (len(i,1)*len(i,2)*len(i,3))/(4*area(i,1));
end

%%% element's edge / edge1(1-2), edge2(2-3), edge3(3-1)
edge = zeros(Erow,1);
for i=1:Erow
    for j=1:3
        if j==1
            edge(i,1) = real(sqrt(rad(i,1)*rad(i,1) - (len(i,1)/2)*(len(i,1)/2)));
            if edge(i,1) < 2e-17
                edge(i,1) = 0.0;
            end
        elseif j==2
            edge(i,2) = real(sqrt(rad(i,1)*rad(i,1) - (len(i,2)/2)*(len(i,2)/2)));
            if edge(i,2) < 2e-17
                edge(i,2) = 0.0;
            end
        elseif j==3
            edge(i,3) = real(sqrt(rad(i,1)*rad(i,1) - (len(i,3)/2)*(len(i,3)/2)));
            if edge(i,3) < 2e-17
                edge(i,3) = 0.0;
            end
        end
    end
end

%%% interface edge
is1=intersect(E1, Esi); is1row = size(is1,1);
interedge1 = zeros(is1row-1,2);
j=1;
for i = 1:E1row
    if size(intersect(E1(i,:),is1'),2) == 2
        interedge1(j,:) = intersect(E1(i,:),is1');
        j=j+1;
    end
end

is2=intersect(E5, Esi); is2row = size(is2,1);
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
%%% Initial potential
phi = zeros(size(X1,1),1);
for i=1:VE3row % substrate
    phi(X1_sorted(VE3(i),2)) = Vt*log(-Na/nint);
    if X1_sorted(VE3(i),1) ~= 0
        phi(X1_sorted(VE3(i),1)) = Vt*log(-Na/nint);
    elseif X1_sorted(VE3(i),3) ~= 0
        phi(X1_sorted(VE3(i),3)) = Vt*log(-Na/nint);
    end
end

for i=1:VE2row % source
    phi(X1_sorted(VE2(i),2)) = Vt*log(Nd/nint);
    if X1_sorted(VE2(i),1) ~= 0
        phi(X1_sorted(VE2(i),1)) = Vt*log(Nd/nint);
    elseif X1_sorted(VE2(i),3) ~= 0
        phi(X1_sorted(VE2(i),3)) = Vt*log(Nd/nint);
    end
end

for i=1:VE4row % drain
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
            phi(X1(fcon,2)) = 0.33374; % gate work function difference
        elseif i>Ncon1 && i<=(Ncon1+Ncon2)
            phi(X1(fcon,2)) = 0.33374; % gate work function difference
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
%%% Non-linear Possion equation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for it=1:Niter

    jac = sparse(size(X2_sorted,1),size(X2_sorted,1));
    res = zeros(size(X2_sorted,1),1);

    for i=1:E1row
        jac(X2(3*(X1_sorted(E1(i,1),1))-2,1),X2(3*(X1_sorted(E1(i,1),1))-2,1)) = jac(X2(3*(X1_sorted(E1(i,1),1))-2,1),X2(3*(X1_sorted(E1(i,1),1))-2,1))...
            - eox*(edgeE1(i,1)/lenE1(i,1) + edgeE1(i,3)/lenE1(i,3));
        jac(X2(3*(X1_sorted(E1(i,1),1))-2,1),X2(3*(X1_sorted(E1(i,2),1))-2,1)) = jac(X2(3*(X1_sorted(E1(i,1),1))-2,1),X2(3*(X1_sorted(E1(i,2),1))-2,1))...
            + eox*edgeE1(i,1)/lenE1(i,1);
        jac(X2(3*(X1_sorted(E1(i,1),1))-2,1),X2(3*(X1_sorted(E1(i,3),1))-2,1)) = jac(X2(3*(X1_sorted(E1(i,1),1))-2,1),X2(3*(X1_sorted(E1(i,3),1))-2,1))...
            + eox*edgeE1(i,3)/lenE1(i,3);
        res(X2(3*(X1_sorted(E1(i,1),1))-2,1),1) = res(X2(3*(X1_sorted(E1(i,1),1))-2,1),1)...
            - eox*(edgeE1(i,1)/lenE1(i,1) + edgeE1(i,3)/lenE1(i,3))*phi(X1_sorted(E1(i,1),1))...
            + eox*edgeE1(i,1)/lenE1(i,1)*phi(X1_sorted(E1(i,2),1))...
            + eox*edgeE1(i,3)/lenE1(i,3)*phi(X1_sorted(E1(i,3),1));

        jac(X2(3*(X1_sorted(E1(i,2),1))-2,1),X2(3*(X1_sorted(E1(i,2),1))-2,1)) = jac(X2(3*(X1_sorted(E1(i,2),1))-2,1),X2(3*(X1_sorted(E1(i,2),1))-2,1))...
            - eox*(edgeE1(i,1)/lenE1(i,1) + edgeE1(i,2)/lenE1(i,2));
        jac(X2(3*(X1_sorted(E1(i,2),1))-2,1),X2(3*(X1_sorted(E1(i,1),1))-2,1)) = jac(X2(3*(X1_sorted(E1(i,2),1))-2,1),X2(3*(X1_sorted(E1(i,1),1))-2,1))...
            + eox*edgeE1(i,1)/lenE1(i,1);
        jac(X2(3*(X1_sorted(E1(i,2),1))-2,1),X2(3*(X1_sorted(E1(i,3),1))-2,1)) = jac(X2(3*(X1_sorted(E1(i,2),1))-2,1),X2(3*(X1_sorted(E1(i,3),1))-2,1))...
            + eox*edgeE1(i,2)/lenE1(i,2);
        res(X2(3*(X1_sorted(E1(i,2),1))-2,1),1) = res(X2(3*(X1_sorted(E1(i,2),1))-2,1),1)...
            - eox*(edgeE1(i,1)/lenE1(i,1) + edgeE1(i,2)/lenE1(i,2))*phi(X1_sorted(E1(i,2),1))...
            + eox*edgeE1(i,1)/lenE1(i,1)*phi(X1_sorted(E1(i,1),1))...
            + eox*edgeE1(i,2)/lenE1(i,2)*phi(X1_sorted(E1(i,3),1));

        jac(X2(3*(X1_sorted(E1(i,3),1))-2,1),X2(3*(X1_sorted(E1(i,3),1))-2,1)) = jac(X2(3*(X1_sorted(E1(i,3),1))-2,1),X2(3*(X1_sorted(E1(i,3),1))-2,1))...
            - eox*(edgeE1(i,2)/lenE1(i,2) + edgeE1(i,3)/lenE1(i,3));
        jac(X2(3*(X1_sorted(E1(i,3),1))-2,1),X2(3*(X1_sorted(E1(i,2),1))-2,1)) = jac(X2(3*(X1_sorted(E1(i,3),1))-2,1),X2(3*(X1_sorted(E1(i,2),1))-2,1))...
            + eox*edgeE1(i,2)/lenE1(i,2);
        jac(X2(3*(X1_sorted(E1(i,3),1))-2,1),X2(3*(X1_sorted(E1(i,1),1))-2,1)) = jac(X2(3*(X1_sorted(E1(i,3),1))-2,1),X2(3*(X1_sorted(E1(i,1),1))-2,1))...
            + eox*edgeE1(i,3)/lenE1(i,3);
        res(X2(3*(X1_sorted(E1(i,3),1))-2,1),1) = res(X2(3*(X1_sorted(E1(i,3),1))-2,1),1)...
            - eox*(edgeE1(i,2)/lenE1(i,2) + edgeE1(i,3)/lenE1(i,3))*phi(X1_sorted(E1(i,3),1))...
            + eox*edgeE1(i,2)/lenE1(i,2)*phi(X1_sorted(E1(i,2),1))...
            + eox*edgeE1(i,3)/lenE1(i,3)*phi(X1_sorted(E1(i,1),1));
    end

    saveexl = zeros(V_silicon_row,1);

    for i=1:Esirow
        if i <= E2row
            jac(X2(3*(X1_sorted(Esi(i,1),2))-2,1),X2(3*(X1_sorted(Esi(i,1),2))-2,1)) = jac(X2(3*(X1_sorted(Esi(i,1),2))-2,1),X2(3*(X1_sorted(Esi(i,1),2))-2,1))...
                - esi*(edgeEsi(i,1)/lenEsi(i,1) + edgeEsi(i,3)/lenEsi(i,3));
            jac(X2(3*(X1_sorted(Esi(i,1),2))-2,1),X2(3*(X1_sorted(Esi(i,2),2))-2,1)) = jac(X2(3*(X1_sorted(Esi(i,1),2))-2,1),X2(3*(X1_sorted(Esi(i,2),2))-2,1))...
                + esi*edgeEsi(i,1)/lenEsi(i,1);
            jac(X2(3*(X1_sorted(Esi(i,1),2))-2,1),X2(3*(X1_sorted(Esi(i,3),2))-2,1)) = jac(X2(3*(X1_sorted(Esi(i,1),2))-2,1),X2(3*(X1_sorted(Esi(i,3),2))-2,1))...
                + esi*edgeEsi(i,3)/lenEsi(i,3);
            res(X2(3*(X1_sorted(Esi(i,1),2))-2,1),1) =  res(X2(3*(X1_sorted(Esi(i,1),2))-2,1),1)...
                + (edgeEsi(i,1)*lenEsi(i,1) + edgeEsi(i,3)*lenEsi(i,3))/4*q/e0*Nd...
                - esi*(edgeEsi(i,1)/lenEsi(i,1) + edgeEsi(i,3)/lenEsi(i,3))*phi(X1_sorted(Esi(i,1),2))...
                + esi*edgeEsi(i,1)/lenEsi(i,1)*phi(X1_sorted(Esi(i,2),2))...
                + esi*edgeEsi(i,3)/lenEsi(i,3)*phi(X1_sorted(Esi(i,3),2));
            exl = find(Esi(i,1)==V_silicon(:,1));
            saveexl(exl,1) =  saveexl(exl,1) + 0.25*(edgeEsi(i,1)*lenEsi(i,1) + edgeEsi(i,3)*lenEsi(i,3));

            jac(X2(3*(X1_sorted(Esi(i,2),2))-2,1),X2(3*(X1_sorted(Esi(i,2),2))-2,1)) = jac(X2(3*(X1_sorted(Esi(i,2),2))-2,1),X2(3*(X1_sorted(Esi(i,2),2))-2,1))...
                - esi*(edgeEsi(i,1)/lenEsi(i,1) + edgeEsi(i,2)/lenEsi(i,2));
            jac(X2(3*(X1_sorted(Esi(i,2),2))-2,1),X2(3*(X1_sorted(Esi(i,1),2))-2,1)) = jac(X2(3*(X1_sorted(Esi(i,2),2))-2,1),X2(3*(X1_sorted(Esi(i,1),2))-2,1))...
                + esi*edgeEsi(i,1)/lenEsi(i,1);
            jac(X2(3*(X1_sorted(Esi(i,2),2))-2,1),X2(3*(X1_sorted(Esi(i,3),2))-2,1)) = jac(X2(3*(X1_sorted(Esi(i,2),2))-2,1),X2(3*(X1_sorted(Esi(i,3),2))-2,1))...
                + esi*edgeEsi(i,2)/lenEsi(i,2);
            res(X2(3*(X1_sorted(Esi(i,2),2))-2,1),1) =  res(X2(3*(X1_sorted(Esi(i,2),2))-2,1),1)...
                + (edgeEsi(i,1)*lenEsi(i,1) + edgeEsi(i,2)*lenEsi(i,2))/4*q/e0*Nd...
                - esi*(edgeEsi(i,1)/lenEsi(i,1) + edgeEsi(i,2)/lenEsi(i,2))*phi(X1_sorted(Esi(i,2),2))...
                + esi*edgeEsi(i,1)/lenEsi(i,1)*phi(X1_sorted(Esi(i,1),2))...
                + esi*edgeEsi(i,2)/lenEsi(i,2)*phi(X1_sorted(Esi(i,3),2));
            exl = find(Esi(i,2)==V_silicon(:,1));
            saveexl(exl,1) =  saveexl(exl,1) + 0.25*(edgeEsi(i,1)*lenEsi(i,1) + edgeEsi(i,2)*lenEsi(i,2));

            jac(X2(3*(X1_sorted(Esi(i,3),2))-2,1),X2(3*(X1_sorted(Esi(i,3),2))-2,1)) = jac(X2(3*(X1_sorted(Esi(i,3),2))-2,1),X2(3*(X1_sorted(Esi(i,3),2))-2,1))...
                - esi*(edgeEsi(i,2)/lenEsi(i,2) + edgeEsi(i,3)/lenEsi(i,3));
            jac(X2(3*(X1_sorted(Esi(i,3),2))-2,1),X2(3*(X1_sorted(Esi(i,2),2))-2,1)) = jac(X2(3*(X1_sorted(Esi(i,3),2))-2,1),X2(3*(X1_sorted(Esi(i,2),2))-2,1))...
                + esi*edgeEsi(i,2)/lenEsi(i,2);
            jac(X2(3*(X1_sorted(Esi(i,3),2))-2,1),X2(3*(X1_sorted(Esi(i,1),2))-2,1)) = jac(X2(3*(X1_sorted(Esi(i,3),2))-2,1),X2(3*(X1_sorted(Esi(i,1),2))-2,1))...
                + esi*edgeEsi(i,3)/lenEsi(i,3);
            res(X2(3*(X1_sorted(Esi(i,3),2))-2,1),1) =  res(X2(3*(X1_sorted(Esi(i,3),2))-2,1),1)...
                + (edgeEsi(i,2)*lenEsi(i,2) + edgeEsi(i,3)*lenEsi(i,3))/4*q/e0*Nd...
                - esi*(edgeEsi(i,2)/lenEsi(i,2) + edgeEsi(i,3)/lenEsi(i,3))*phi(X1_sorted(Esi(i,3),2))...
                + esi*edgeEsi(i,2)/lenEsi(i,2)*phi(X1_sorted(Esi(i,2),2))...
                + esi*edgeEsi(i,3)/lenEsi(i,3)*phi(X1_sorted(Esi(i,1),2));
            exl = find(Esi(i,3)==V_silicon(:,1));
            saveexl(exl,1) =  saveexl(exl,1) + 0.25*(edgeEsi(i,2)*lenEsi(i,2) + edgeEsi(i,3)*lenEsi(i,3));

        elseif i > E2row && i <= (E2row + E3row)
            jac(X2(3*(X1_sorted(Esi(i,1),2))-2,1),X2(3*(X1_sorted(Esi(i,1),2))-2,1)) = jac(X2(3*(X1_sorted(Esi(i,1),2))-2,1),X2(3*(X1_sorted(Esi(i,1),2))-2,1))...
                - esi*(edgeEsi(i,1)/lenEsi(i,1) + edgeEsi(i,3)/lenEsi(i,3));
            jac(X2(3*(X1_sorted(Esi(i,1),2))-2,1),X2(3*(X1_sorted(Esi(i,2),2))-2,1)) = jac(X2(3*(X1_sorted(Esi(i,1),2))-2,1),X2(3*(X1_sorted(Esi(i,2),2))-2,1))...
                + esi*edgeEsi(i,1)/lenEsi(i,1);
            jac(X2(3*(X1_sorted(Esi(i,1),2))-2,1),X2(3*(X1_sorted(Esi(i,3),2))-2,1)) = jac(X2(3*(X1_sorted(Esi(i,1),2))-2,1),X2(3*(X1_sorted(Esi(i,3),2))-2,1))...
                + esi*edgeEsi(i,3)/lenEsi(i,3);
            res(X2(3*(X1_sorted(Esi(i,1),2))-2,1),1) =  res(X2(3*(X1_sorted(Esi(i,1),2))-2,1),1)...
                + (edgeEsi(i,1)*lenEsi(i,1) + edgeEsi(i,3)*lenEsi(i,3))/4*q/e0*Na...
                - esi*(edgeEsi(i,1)/lenEsi(i,1) + edgeEsi(i,3)/lenEsi(i,3))*phi(X1_sorted(Esi(i,1),2))...
                + esi*edgeEsi(i,1)/lenEsi(i,1)*phi(X1_sorted(Esi(i,2),2))...
                + esi*edgeEsi(i,3)/lenEsi(i,3)*phi(X1_sorted(Esi(i,3),2));
            exl = find(Esi(i,1)==V_silicon(:,1));
            saveexl(exl,1) =  saveexl(exl,1) + 0.25*(edgeEsi(i,1)*lenEsi(i,1) + edgeEsi(i,3)*lenEsi(i,3));

            jac(X2(3*(X1_sorted(Esi(i,2),2))-2,1),X2(3*(X1_sorted(Esi(i,2),2))-2,1)) = jac(X2(3*(X1_sorted(Esi(i,2),2))-2,1),X2(3*(X1_sorted(Esi(i,2),2))-2,1))...
                - esi*(edgeEsi(i,1)/lenEsi(i,1) + edgeEsi(i,2)/lenEsi(i,2));
            jac(X2(3*(X1_sorted(Esi(i,2),2))-2,1),X2(3*(X1_sorted(Esi(i,1),2))-2,1)) = jac(X2(3*(X1_sorted(Esi(i,2),2))-2,1),X2(3*(X1_sorted(Esi(i,1),2))-2,1))...
                + esi*edgeEsi(i,1)/lenEsi(i,1);
            jac(X2(3*(X1_sorted(Esi(i,2),2))-2,1),X2(3*(X1_sorted(Esi(i,3),2))-2,1)) = jac(X2(3*(X1_sorted(Esi(i,2),2))-2,1),X2(3*(X1_sorted(Esi(i,3),2))-2,1))...
                + esi*edgeEsi(i,2)/lenEsi(i,2);
            res(X2(3*(X1_sorted(Esi(i,2),2))-2,1),1) =  res(X2(3*(X1_sorted(Esi(i,2),2))-2,1),1)...
                + (edgeEsi(i,1)*lenEsi(i,1) + edgeEsi(i,2)*lenEsi(i,2))/4*q/e0*Na...
                - esi*(edgeEsi(i,1)/lenEsi(i,1) + edgeEsi(i,2)/lenEsi(i,2))*phi(X1_sorted(Esi(i,2),2))...
                + esi*edgeEsi(i,1)/lenEsi(i,1)*phi(X1_sorted(Esi(i,1),2))...
                + esi*edgeEsi(i,2)/lenEsi(i,2)*phi(X1_sorted(Esi(i,3),2));
            exl = find(Esi(i,2)==V_silicon(:,1));
            saveexl(exl,1) =  saveexl(exl,1) + 0.25*(edgeEsi(i,1)*lenEsi(i,1) + edgeEsi(i,2)*lenEsi(i,2));

            jac(X2(3*(X1_sorted(Esi(i,3),2))-2,1),X2(3*(X1_sorted(Esi(i,3),2))-2,1)) = jac(X2(3*(X1_sorted(Esi(i,3),2))-2,1),X2(3*(X1_sorted(Esi(i,3),2))-2,1))...
                - esi*(edgeEsi(i,2)/lenEsi(i,2) + edgeEsi(i,3)/lenEsi(i,3));
            jac(X2(3*(X1_sorted(Esi(i,3),2))-2,1),X2(3*(X1_sorted(Esi(i,2),2))-2,1)) = jac(X2(3*(X1_sorted(Esi(i,3),2))-2,1),X2(3*(X1_sorted(Esi(i,2),2))-2,1))...
                + esi*edgeEsi(i,2)/lenEsi(i,2);
            jac(X2(3*(X1_sorted(Esi(i,3),2))-2,1),X2(3*(X1_sorted(Esi(i,1),2))-2,1)) = jac(X2(3*(X1_sorted(Esi(i,3),2))-2,1),X2(3*(X1_sorted(Esi(i,1),2))-2,1))...
                + esi*edgeEsi(i,3)/lenEsi(i,3);
            res(X2(3*(X1_sorted(Esi(i,3),2))-2,1),1) =  res(X2(3*(X1_sorted(Esi(i,3),2))-2,1),1)...
                + (edgeEsi(i,2)*lenEsi(i,2) + edgeEsi(i,3)*lenEsi(i,3))/4*q/e0*Na...
                - esi*(edgeEsi(i,2)/lenEsi(i,2) + edgeEsi(i,3)/lenEsi(i,3))*phi(X1_sorted(Esi(i,3),2))...
                + esi*edgeEsi(i,2)/lenEsi(i,2)*phi(X1_sorted(Esi(i,2),2))...
                + esi*edgeEsi(i,3)/lenEsi(i,3)*phi(X1_sorted(Esi(i,1),2));
            exl = find(Esi(i,3)==V_silicon(:,1));
            saveexl(exl,1) =  saveexl(exl,1) + 0.25*(edgeEsi(i,2)*lenEsi(i,2) + edgeEsi(i,3)*lenEsi(i,3));

        elseif i > (E2row+E3row) && i <= Esirow
            jac(X2(3*(X1_sorted(Esi(i,1),2))-2,1),X2(3*(X1_sorted(Esi(i,1),2))-2,1)) = jac(X2(3*(X1_sorted(Esi(i,1),2))-2,1),X2(3*(X1_sorted(Esi(i,1),2))-2,1))...
                - esi*(edgeEsi(i,1)/lenEsi(i,1) + edgeEsi(i,3)/lenEsi(i,3));
            jac(X2(3*(X1_sorted(Esi(i,1),2))-2,1),X2(3*(X1_sorted(Esi(i,2),2))-2,1)) = jac(X2(3*(X1_sorted(Esi(i,1),2))-2,1),X2(3*(X1_sorted(Esi(i,2),2))-2,1))...
                + esi*edgeEsi(i,1)/lenEsi(i,1);
            jac(X2(3*(X1_sorted(Esi(i,1),2))-2,1),X2(3*(X1_sorted(Esi(i,3),2))-2,1)) = jac(X2(3*(X1_sorted(Esi(i,1),2))-2,1),X2(3*(X1_sorted(Esi(i,3),2))-2,1))...
                + esi*edgeEsi(i,3)/lenEsi(i,3);
            res(X2(3*(X1_sorted(Esi(i,1),2))-2,1),1) =  res(X2(3*(X1_sorted(Esi(i,1),2))-2,1),1)...
                + (edgeEsi(i,1)*lenEsi(i,1) + edgeEsi(i,3)*lenEsi(i,3))/4*q/e0*Nd...
                - esi*(edgeEsi(i,1)/lenEsi(i,1) + edgeEsi(i,3)/lenEsi(i,3))*phi(X1_sorted(Esi(i,1),2))...
                + esi*edgeEsi(i,1)/lenEsi(i,1)*phi(X1_sorted(Esi(i,2),2))...
                + esi*edgeEsi(i,3)/lenEsi(i,3)*phi(X1_sorted(Esi(i,3),2));
            exl = find(Esi(i,1)==V_silicon(:,1));
            saveexl(exl,1) =  saveexl(exl,1) + 0.25*(edgeEsi(i,1)*lenEsi(i,1) + edgeEsi(i,3)*lenEsi(i,3));

            jac(X2(3*(X1_sorted(Esi(i,2),2))-2,1),X2(3*(X1_sorted(Esi(i,2),2))-2,1)) = jac(X2(3*(X1_sorted(Esi(i,2),2))-2,1),X2(3*(X1_sorted(Esi(i,2),2))-2,1))...
                - esi*(edgeEsi(i,1)/lenEsi(i,1) + edgeEsi(i,2)/lenEsi(i,2));
            jac(X2(3*(X1_sorted(Esi(i,2),2))-2,1),X2(3*(X1_sorted(Esi(i,1),2))-2,1)) = jac(X2(3*(X1_sorted(Esi(i,2),2))-2,1),X2(3*(X1_sorted(Esi(i,1),2))-2,1))...
                + esi*edgeEsi(i,1)/lenEsi(i,1);
            jac(X2(3*(X1_sorted(Esi(i,2),2))-2,1),X2(3*(X1_sorted(Esi(i,3),2))-2,1)) = jac(X2(3*(X1_sorted(Esi(i,2),2))-2,1),X2(3*(X1_sorted(Esi(i,3),2))-2,1))...
                + esi*edgeEsi(i,2)/lenEsi(i,2);
            res(X2(3*(X1_sorted(Esi(i,2),2))-2,1),1) =  res(X2(3*(X1_sorted(Esi(i,2),2))-2,1),1)...
                + (edgeEsi(i,1)*lenEsi(i,1) + edgeEsi(i,2)*lenEsi(i,2))/4*q/e0*Nd...
                - esi*(edgeEsi(i,1)/lenEsi(i,1) + edgeEsi(i,2)/lenEsi(i,2))*phi(X1_sorted(Esi(i,2),2))...
                + esi*edgeEsi(i,1)/lenEsi(i,1)*phi(X1_sorted(Esi(i,1),2))...
                + esi*edgeEsi(i,2)/lenEsi(i,2)*phi(X1_sorted(Esi(i,3),2));
            exl = find(Esi(i,2)==V_silicon(:,1));
            saveexl(exl,1) = saveexl(exl,1) + 0.25*(edgeEsi(i,1)*lenEsi(i,1) + edgeEsi(i,2)*lenEsi(i,2));

            jac(X2(3*X1_sorted(Esi(i,3),2)-2,1),X2(3*(X1_sorted(Esi(i,3),2))-2,1)) = jac(X2(3*(X1_sorted(Esi(i,3),2))-2,1),X2(3*(X1_sorted(Esi(i,3),2))-2,1))...
                - esi*(edgeEsi(i,2)/lenEsi(i,2) + edgeEsi(i,3)/lenEsi(i,3));
            jac(X2(3*X1_sorted(Esi(i,3),2)-2,1),X2(3*(X1_sorted(Esi(i,2),2))-2,1)) = jac(X2(3*(X1_sorted(Esi(i,3),2))-2,1),X2(3*(X1_sorted(Esi(i,2),2))-2,1))...
                + esi*edgeEsi(i,2)/lenEsi(i,2);
            jac(X2(3*X1_sorted(Esi(i,3),2)-2,1),X2(3*(X1_sorted(Esi(i,1),2))-2,1)) = jac(X2(3*(X1_sorted(Esi(i,3),2))-2,1),X2(3*(X1_sorted(Esi(i,1),2))-2,1))...
                + esi*edgeEsi(i,3)/lenEsi(i,3);
            res(X2(3*X1_sorted(Esi(i,3),2)-2,1),1) =  res(X2(3*X1_sorted(Esi(i,3),2)-2,1),1)...
                + (edgeEsi(i,2)*lenEsi(i,2) + edgeEsi(i,3)*lenEsi(i,3))/4*q/e0*Nd...
                - esi*(edgeEsi(i,2)/lenEsi(i,2) + edgeEsi(i,3)/lenEsi(i,3))*phi(X1_sorted(Esi(i,3),2))...
                + esi*edgeEsi(i,2)/lenEsi(i,2)*phi(X1_sorted(Esi(i,2),2))...
                + esi*edgeEsi(i,3)/lenEsi(i,3)*phi(X1_sorted(Esi(i,1),2));
            exl = find(Esi(i,3)==V_silicon(:,1));
            saveexl(exl,1) = saveexl(exl,1) + 0.25*(edgeEsi(i,2)*lenEsi(i,2) + edgeEsi(i,3)*lenEsi(i,3));
        end
    end

    for i=1:E5row
        jac(X2(3*(X1_sorted(E5(i,1),3))-2,1),X2(3*(X1_sorted(E5(i,1),3))-2,1)) = jac(X2(3*(X1_sorted(E5(i,1),3))-2,1),X2(3*(X1_sorted(E5(i,1),3))-2,1))...
            - eox*(edgeE5(i,1)/lenE5(i,1) + edgeE5(i,3)/lenE5(i,3));
        jac(X2(3*(X1_sorted(E5(i,1),3))-2,1),X2(3*(X1_sorted(E5(i,2),3))-2,1)) = jac(X2(3*(X1_sorted(E5(i,1),3))-2,1),X2(3*(X1_sorted(E5(i,2),3))-2,1))...
            + eox*edgeE5(i,1)/lenE5(i,1);
        jac(X2(3*(X1_sorted(E5(i,1),3))-2,1),X2(3*(X1_sorted(E5(i,3),3))-2,1)) = jac(X2(3*(X1_sorted(E5(i,1),3))-2,1),X2(3*(X1_sorted(E5(i,3),3))-2,1))...
            + eox*edgeE5(i,3)/lenE5(i,3);
        res(X2(3*(X1_sorted(E5(i,1),3))-2,1),1) = res(X2(3*(X1_sorted(E5(i,1),3))-2,1),1)...
            - eox*(edgeE5(i,1)/lenE5(i,1) + edgeE5(i,3)/lenE5(i,3))*phi(X1_sorted(E5(i,1),3))...
            + eox*edgeE5(i,1)/lenE5(i,1)*phi(X1_sorted(E5(i,2),3))...
            + eox*edgeE5(i,3)/lenE5(i,3)*phi(X1_sorted(E5(i,3),3));

        jac(X2(3*(X1_sorted(E5(i,2),3))-2,1),X2(3*(X1_sorted(E5(i,2),3))-2,1)) = jac(X2(3*(X1_sorted(E5(i,2),3))-2,1),X2(3*(X1_sorted(E5(i,2),3))-2,1))...
            - eox*(edgeE5(i,1)/lenE5(i,1) + edgeE5(i,2)/lenE5(i,2));
        jac(X2(3*(X1_sorted(E5(i,2),3))-2,1),X2(3*(X1_sorted(E5(i,1),3))-2,1)) = jac(X2(3*(X1_sorted(E5(i,2),3))-2,1),X2(3*(X1_sorted(E5(i,1),3))-2,1))...
            + eox*edgeE5(i,1)/lenE5(i,1);
        jac(X2(3*(X1_sorted(E5(i,2),3))-2,1),X2(3*(X1_sorted(E5(i,3),3))-2,1)) = jac(X2(3*(X1_sorted(E5(i,2),3))-2,1),X2(3*(X1_sorted(E5(i,3),3))-2,1))...
            + eox*edgeE5(i,2)/lenE5(i,2);
        res(X2(3*(X1_sorted(E5(i,2),3))-2,1),1) = res(X2(3*(X1_sorted(E5(i,2),3))-2,1),1)...
            - eox*(edgeE5(i,1)/lenE5(i,1) + edgeE5(i,2)/lenE5(i,2))*phi(X1_sorted(E5(i,2),3))...
            + eox*edgeE5(i,1)/lenE5(i,1)*phi(X1_sorted(E5(i,1),3))...
            + eox*edgeE5(i,2)/lenE5(i,2)*phi(X1_sorted(E5(i,3),3));

        jac(X2(3*(X1_sorted(E5(i,3),3))-2,1),X2(3*(X1_sorted(E5(i,3),3))-2,1)) = jac(X2(3*(X1_sorted(E5(i,3),3))-2,1),X2(3*(X1_sorted(E5(i,3),3))-2,1))...
            - eox*(edgeE5(i,2)/lenE5(i,2) + edgeE5(i,3)/lenE5(i,3));
        jac(X2(3*(X1_sorted(E5(i,3),3))-2,1),X2(3*(X1_sorted(E5(i,2),3))-2,1)) = jac(X2(3*(X1_sorted(E5(i,3),3))-2,1),X2(3*(X1_sorted(E5(i,2),3))-2,1))...
            + eox*edgeE5(i,2)/lenE5(i,2);
        jac(X2(3*(X1_sorted(E5(i,3),3))-2,1),X2(3*(X1_sorted(E5(i,1),3))-2,1)) = jac(X2(3*(X1_sorted(E5(i,3),3))-2,1),X2(3*(X1_sorted(E5(i,1),3))-2,1))...
            + eox*edgeE5(i,3)/lenE5(i,3);
        res(X2(3*(X1_sorted(E5(i,3),3))-2,1),1) = res(X2(3*(X1_sorted(E5(i,3),3))-2,1),1)...
            - eox*(edgeE5(i,2)/lenE5(i,2) + edgeE5(i,3)/lenE5(i,3))*phi(X1_sorted(E5(i,3),3))...
            + eox*edgeE5(i,2)/lenE5(i,2)*phi(X1_sorted(E5(i,2),3))...
            + eox*edgeE5(i,3)/lenE5(i,3)*phi(X1_sorted(E5(i,1),3));
    end

    % electron charge
    for i=1:V_silicon_row
        % n-nit*exp(phi/Vt)
        jac(X2(3*V_silicon(i,2)-1),X2(3*V_silicon(i,2)-1)) = 1;
        jac(X2(3*V_silicon(i,2)-1),X2(3*V_silicon(i,2)-2)) = - nint/Vt*exp(phi(V_silicon(i,2))/Vt);
        res(X2(3*V_silicon(i,2)-1),1) = n(i) - nint*exp(phi(V_silicon(i,2))/Vt);

        % - qn
        jac(X2(3*V_silicon(i,2)-2),X2(3*V_silicon(i,2)-1)) = - saveexl(i,1)*q/e0;
        res(X2(3*V_silicon(i,2)-2),1) = res(X2(3*V_silicon(i,2)-2),1) - saveexl(i,1)*q/e0*n(i);
    end

    % hole charge
    for i=1:V_silicon_row
        % p-nit*exp(-phi/Vt)
        jac(X2(3*V_silicon(i,2)),X2(3*V_silicon(i,2))) = 1;
        jac(X2(3*V_silicon(i,2)),X2(3*V_silicon(i,2)-2)) = nint/Vt*exp(-phi(V_silicon(i,2))/Vt);
        res(X2(3*V_silicon(i,2)),1) = p(i) - nint*exp(-phi(V_silicon(i,2))/Vt);

        % + qp
        jac(X2(3*V_silicon(i,2)-2),X2(3*V_silicon(i,2))) = saveexl(i,1)*q/e0;
        res(X2(3*V_silicon(i,2)-2),1) = res(X2(3*V_silicon(i,2)-2),1) + saveexl(i,1)*q/e0*p(i);
    end

    for i=1:is1row
        minus1 = X1_sorted(is1(i),1);
        plus1 = X1_sorted(is1(i),2);
        jac(X2(3*plus1-2),:) = jac(X2(3*plus1-2),:) + jac(X2(3*minus1-2),:);
        res(X2(3*plus1-2)) = res(X2(3*plus1-2)) + res(X2(3*minus1-2));
        jac(X2(3*minus1-2),:) = 0;
        res(X2(3*minus1-2),1) = 0;
        jac(X2(3*minus1-2), X2(3*minus1-2)) = -1;
        jac(X2(3*minus1-2), X2(3*plus1-2)) = 1;
    end

    for i=1:is2row
        minus1 = X1_sorted(is2(i),3);
        plus1 = X1_sorted(is2(i),2);
        jac(X2(3*plus1-2),:) = jac(X2(3*plus1-2),:) + jac(X2(3*minus1-2),:);
        res(X2(3*plus1-2),1) = res(X2(3*plus1-2),1) + res(X2(3*minus1-2),1);
        jac(X2(3*minus1-2),:) = 0;
        res(X2(3*minus1-2),1) = 0;
        jac(X2(3*minus1-2), X2(3*minus1-2)) = -1;
        jac(X2(3*minus1-2), X2(3*plus1-2)) = 1;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Dirichlet boundary condition for 1st iteration

    % Top gate
    for i=1:Ncon1
        for j=1:concol
            fcon = X1_sorted(contact(i,j),1);
            jac(X2(3*fcon-2),:) = 0;
            jac(X2(3*fcon-2), X2(3*fcon-2)) = 1;
            res(X2(3*fcon-2),1) = phi(fcon,1) - 0.33374;
        end
    end

    % Bottom gate
    for i=Ncon1+1:(Ncon1+Ncon2)
        for j=1:concol
            fcon = X1_sorted(contact(i,j),3);
            jac(X2(3*fcon-2),:) = 0;
            jac(X2(3*fcon-2), X2(3*fcon-2)) = 1;
            res(X2(3*fcon-2),1) = phi(fcon) - 0.33374;
        end
    end

    % source & drain
    for i=(Ncon1+Ncon2)+1:(Ncon1+Ncon2+Ncon3+Ncon4)
        for j=1:concol
            fcon = X1_sorted(contact(i,j),2);
            fcon_np = find(contact(i,j)==V_silicon(:,1));

            jac(X2(3*fcon-2),:) = 0;
            jac(X2(3*fcon-1),:) = 0;
            jac(X2(3*fcon),:) = 0;

            jac(X2(3*fcon-2), X2(3*fcon-2)) = 1;
            jac(X2(3*fcon-1), X2(3*fcon-1)) = 1;
            jac(X2(3*fcon), X2(3*fcon)) = 1;

            if i<=(Ncon1+Ncon2+Ncon3) % source
                res(X2(3*fcon-2),1) = phi(fcon) - Vt*log(Nd/nint);
                res(X2(3*fcon-1),1) = n(fcon_np) - Nd;
                res(X2(3*fcon),1) = p(fcon_np) - nint*nint/Nd;

            elseif i>(Ncon1+Ncon2+Ncon3) &&  i<=(Ncon1+Ncon2+Ncon3+Ncon4) % drain
                res(X2(3*fcon-2),1) = phi(fcon) - Vt*log(Nd/nint);
                res(X2(3*fcon-1),1) = n(fcon_np) - Nd;
                res(X2(3*fcon),1) = p(fcon_np) - nint*nint/Nd;
            end
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% scaled part
    cvector = zeros(size(jac,1),1);
    for i=1:size(jac,1)
        if rem(X2_sorted(i),3) == 1
            cvector(i) = Vt;
        elseif rem(X2_sorted(i),3) == 2
            cvector(i) = Nd;
        elseif rem(X2_sorted(i),3) == 0
            cvector(i) = Nd;
        end
    end

    cmatrix = spdiags(cvector,0,size(X2_sorted,1),size(X2_sorted,1));
    jac_scaled = jac * cmatrix;
    rvector = 1./sum(abs(jac_scaled),2);
    rmatrix = spdiags(rvector,0,size(X2_sorted,1),size(X2_sorted,1));
    jac_scaled = rmatrix * jac_scaled;
    res_scaled = rmatrix * res;
    upd_scaled = jac_scaled \ (-res_scaled);
    upd = cmatrix * upd_scaled;

    countphi1=1;
    countn1=1;
    countp1=1;

    maxupdate = 0;

    for i=1:size(X2_sorted,1)
        if rem(X2_sorted(i),3) == 1
            phi(countphi1) = phi(countphi1) + upd(i);

            if maxupdate < max(abs(upd(i)))
                maxupdate = max(abs(upd(i)));
            end
            countphi1 = countphi1 + 1;

        elseif rem(X2_sorted(i),3) == 2
            n(countn1) = n(countn1) + upd(i);
            countn1 = countn1 + 1;

        elseif rem(X2_sorted(i),3) == 0
            p(countp1) = p(countp1) + upd(i);
            countp1 = countp1 + 1;
        end
    end

    save_phi_upd(it,1) = maxupdate;

    if maxupdate < 1e-15
        break
    end

end

before_lamping_phi = phi(:,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Double gate MOSFET gate lamping %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
final_Vgate = 0.6;
gatestep = 0.05;
Ngatestep = round(final_Vgate/gatestep) + 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
save_Vg = zeros(Ngatestep,1);
for step=1:Ngatestep
    Vg = gatestep*(step-1);
    
    save_Vg(step,1) = Vg;

    for it=1:Niter

        Jaco = zeros(size(X2_sorted,1)+13,size(X2_sorted,1)+13);
        Res = zeros(size(X2_sorted,1)+13,1);

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

        saveexl = zeros(V_silicon_row,1);

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
                exl = find(Esi(i,1)==V_silicon(:,1));
                saveexl(exl,1) = saveexl(exl,1) + 0.25*(edgeEsi(i,1)*lenEsi(i,1) + edgeEsi(i,3)*lenEsi(i,3));

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
                exl = find(Esi(i,2)==V_silicon(:,1));
                saveexl(exl,1) = saveexl(exl,1) + 0.25*(edgeEsi(i,1)*lenEsi(i,1) + edgeEsi(i,2)*lenEsi(i,2));

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
                exl = find(Esi(i,3)==V_silicon(:,1));
                saveexl(exl,1) = saveexl(exl,1) + 0.25*(edgeEsi(i,2)*lenEsi(i,2) + edgeEsi(i,3)*lenEsi(i,3));

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
                exl = find(Esi(i,1)==V_silicon(:,1));
                saveexl(exl,1) = saveexl(exl,1) + 0.25*(edgeEsi(i,1)*lenEsi(i,1) + edgeEsi(i,3)*lenEsi(i,3));

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
                exl = find(Esi(i,2)==V_silicon(:,1));
                saveexl(exl,1) = saveexl(exl,1) + 0.25*(edgeEsi(i,1)*lenEsi(i,1) + edgeEsi(i,2)*lenEsi(i,2));

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
                exl = find(Esi(i,3)==V_silicon(:,1));
                saveexl(exl,1) = saveexl(exl,1) + 0.25*(edgeEsi(i,2)*lenEsi(i,2) + edgeEsi(i,3)*lenEsi(i,3));

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
                exl = find(Esi(i,1)==V_silicon(:,1));
                saveexl(exl,1) = saveexl(exl,1) + 0.25*(edgeEsi(i,1)*lenEsi(i,1) + edgeEsi(i,3)*lenEsi(i,3));

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
                exl = find(Esi(i,2)==V_silicon(:,1));
                saveexl(exl,1) = saveexl(exl,1) + 0.25*(edgeEsi(i,1)*lenEsi(i,1) + edgeEsi(i,2)*lenEsi(i,2));

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
                exl = find(Esi(i,3)==V_silicon(:,1));
                saveexl(exl,1) = saveexl(exl,1) + 0.25*(edgeEsi(i,2)*lenEsi(i,2) + edgeEsi(i,3)*lenEsi(i,3));
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

        %  Poisson charge
        for i=1:V_silicon_row
            % - qn
            Jaco(X2(3*V_silicon(i,2)-2),X2(3*V_silicon(i,2)-1)) = Jaco(X2(3*V_silicon(i,2)-2),X2(3*V_silicon(i,2)-1)) - saveexl(i,1)*q/e0;
            Res(X2(3*V_silicon(i,2)-2),1) = Res(X2(3*V_silicon(i,2)-2),1) - saveexl(i,1)*q/e0*n(i);

            % + qp
            Jaco(X2(3*V_silicon(i,2)-2),X2(3*V_silicon(i,2))) =  Jaco(X2(3*V_silicon(i,2)-2),X2(3*V_silicon(i,2))) + saveexl(i,1)*q/e0;
            Res(X2(3*V_silicon(i,2)-2),1) = Res(X2(3*V_silicon(i,2)-2),1) + saveexl(i,1)*q/e0*p(i);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% electron DD
        for i=1:Esirow

            f1 = find(V_silicon(:,1)==Esi(i,1));
            f2 = find(V_silicon(:,1)==Esi(i,2));
            f3 = find(V_silicon(:,1)==Esi(i,3));

            %%% 1st element
            % n part
            Jaco(X2(3*X1_sorted(Esi(i,1),2)-1,1), X2(3*X1_sorted(Esi(i,1),2)-1,1)) = Jaco(X2(3*X1_sorted(Esi(i,1),2)-1,1), X2(3*X1_sorted(Esi(i,1),2)-1,1))...
                - ncoeff*Vt*edgeEsi(i,1)/lenEsi(i,1)*bern(phi(X1_sorted(Esi(i,2),2)),phi(X1_sorted(Esi(i,1),2)))...
                - ncoeff*Vt*edgeEsi(i,3)/lenEsi(i,3)*bern(phi(X1_sorted(Esi(i,3),2)),phi(X1_sorted(Esi(i,1),2)));
            Jaco(X2(3*X1_sorted(Esi(i,1),2)-1,1), X2(3*X1_sorted(Esi(i,2),2)-1,1)) = Jaco(X2(3*X1_sorted(Esi(i,1),2)-1,1), X2(3*X1_sorted(Esi(i,2),2)-1,1))...
                + ncoeff*Vt*edgeEsi(i,1)/lenEsi(i,1)*bern(phi(X1_sorted(Esi(i,1),2)),phi(X1_sorted(Esi(i,2),2)));
            Jaco(X2(3*X1_sorted(Esi(i,1),2)-1,1), X2(3*X1_sorted(Esi(i,3),2)-1,1)) = Jaco(X2(3*X1_sorted(Esi(i,1),2)-1,1), X2(3*X1_sorted(Esi(i,3),2)-1,1))...
                + ncoeff*Vt*edgeEsi(i,3)/lenEsi(i,3)*bern(phi(X1_sorted(Esi(i,1),2)),phi(X1_sorted(Esi(i,3),2)));

            % phi part
            Jaco(X2(3*X1_sorted(Esi(i,1),2)-1,1), X2(3*X1_sorted(Esi(i,1),2)-2,1)) = Jaco(X2(3*X1_sorted(Esi(i,1),2)-1,1), X2(3*X1_sorted(Esi(i,1),2)-2,1))...
                - ncoeff*edgeEsi(i,1)/lenEsi(i,1)*dbern(phi(X1_sorted(Esi(i,1),2)), phi(X1_sorted(Esi(i,2),2)))*n(f2)...
                - ncoeff*edgeEsi(i,1)/lenEsi(i,1)*dbern(phi(X1_sorted(Esi(i,2),2)), phi(X1_sorted(Esi(i,1),2)))*n(f1)...
                - ncoeff*edgeEsi(i,3)/lenEsi(i,3)*dbern(phi(X1_sorted(Esi(i,1),2)), phi(X1_sorted(Esi(i,3),2)))*n(f3)...
                - ncoeff*edgeEsi(i,3)/lenEsi(i,3)*dbern(phi(X1_sorted(Esi(i,3),2)), phi(X1_sorted(Esi(i,1),2)))*n(f1);
            Jaco(X2(3*X1_sorted(Esi(i,1),2)-1,1), X2(3*X1_sorted(Esi(i,2),2)-2,1)) = Jaco(X2(3*X1_sorted(Esi(i,1),2)-1,1), X2(3*X1_sorted(Esi(i,2),2)-2,1))...
                + ncoeff*edgeEsi(i,1)/lenEsi(i,1)*dbern(phi(X1_sorted(Esi(i,1),2)), phi(X1_sorted(Esi(i,2),2)))*n(f2)...
                + ncoeff*edgeEsi(i,1)/lenEsi(i,1)*dbern(phi(X1_sorted(Esi(i,2),2)), phi(X1_sorted(Esi(i,1),2)))*n(f1);
            Jaco(X2(3*X1_sorted(Esi(i,1),2)-1,1), X2(3*X1_sorted(Esi(i,3),2)-2,1)) = Jaco(X2(3*X1_sorted(Esi(i,1),2)-1,1), X2(3*X1_sorted(Esi(i,3),2)-2,1))...
                + ncoeff*edgeEsi(i,3)/lenEsi(i,3)*dbern(phi(X1_sorted(Esi(i,1),2)), phi(X1_sorted(Esi(i,3),2)))*n(f3)...
                + ncoeff*edgeEsi(i,3)/lenEsi(i,3)*dbern(phi(X1_sorted(Esi(i,3),2)), phi(X1_sorted(Esi(i,1),2)))*n(f1);

            Res(X2(3*X1_sorted(Esi(i,1),2)-1,1),1) = Res(X2(3*X1_sorted(Esi(i,1),2)-1,1),1)...
                + ncoeff*Vt*edgeEsi(i,1)/lenEsi(i,1)*bern(phi(X1_sorted(Esi(i,1),2)),phi(X1_sorted(Esi(i,2),2)))*n(f2)...
                - ncoeff*Vt*edgeEsi(i,1)/lenEsi(i,1)*bern(phi(X1_sorted(Esi(i,2),2)),phi(X1_sorted(Esi(i,1),2)))*n(f1)...
                + ncoeff*Vt*edgeEsi(i,3)/lenEsi(i,3)*bern(phi(X1_sorted(Esi(i,1),2)),phi(X1_sorted(Esi(i,3),2)))*n(f3)...
                - ncoeff*Vt*edgeEsi(i,3)/lenEsi(i,3)*bern(phi(X1_sorted(Esi(i,3),2)),phi(X1_sorted(Esi(i,1),2)))*n(f1);

            %%% 2nd element
            % n part
            Jaco(X2(3*X1_sorted(Esi(i,2),2)-1,1), X2(3*X1_sorted(Esi(i,2),2)-1,1)) = Jaco(X2(3*X1_sorted(Esi(i,2),2)-1,1), X2(3*X1_sorted(Esi(i,2),2)-1,1))...
                - ncoeff*Vt*edgeEsi(i,2)/lenEsi(i,2)*bern(phi(X1_sorted(Esi(i,3),2)), phi(X1_sorted(Esi(i,2),2)))...
                - ncoeff*Vt*edgeEsi(i,1)/lenEsi(i,1)*bern(phi(X1_sorted(Esi(i,1),2)), phi(X1_sorted(Esi(i,2),2)));
            Jaco(X2(3*X1_sorted(Esi(i,2),2)-1,1), X2(3*X1_sorted(Esi(i,3),2)-1,1)) = Jaco(X2(3*X1_sorted(Esi(i,2),2)-1,1), X2(3*X1_sorted(Esi(i,3),2)-1,1))...
                + ncoeff*Vt*edgeEsi(i,2)/lenEsi(i,2)*bern(phi(X1_sorted(Esi(i,2),2)), phi(X1_sorted(Esi(i,3),2)));
            Jaco(X2(3*X1_sorted(Esi(i,2),2)-1,1), X2(3*X1_sorted(Esi(i,1),2)-1,1)) = Jaco(X2(3*X1_sorted(Esi(i,2),2)-1,1), X2(3*X1_sorted(Esi(i,1),2)-1,1))...
                + ncoeff*Vt*edgeEsi(i,1)/lenEsi(i,1)*bern(phi(X1_sorted(Esi(i,2),2)), phi(X1_sorted(Esi(i,1),2)));

            % phi part
            Jaco(X2(3*X1_sorted(Esi(i,2),2)-1,1), X2(3*X1_sorted(Esi(i,2),2)-2,1)) = Jaco(X2(3*X1_sorted(Esi(i,2),2)-1,1), X2(3*X1_sorted(Esi(i,2),2)-2,1))...
                - ncoeff*edgeEsi(i,2)/lenEsi(i,2)*dbern(phi(X1_sorted(Esi(i,2),2)), phi(X1_sorted(Esi(i,3),2)))*n(f3)...
                - ncoeff*edgeEsi(i,2)/lenEsi(i,2)*dbern(phi(X1_sorted(Esi(i,3),2)), phi(X1_sorted(Esi(i,2),2)))*n(f2)...
                - ncoeff*edgeEsi(i,1)/lenEsi(i,1)*dbern(phi(X1_sorted(Esi(i,2),2)), phi(X1_sorted(Esi(i,1),2)))*n(f1)...
                - ncoeff*edgeEsi(i,1)/lenEsi(i,1)*dbern(phi(X1_sorted(Esi(i,1),2)), phi(X1_sorted(Esi(i,2),2)))*n(f2);
            Jaco(X2(3*X1_sorted(Esi(i,2),2)-1,1), X2(3*X1_sorted(Esi(i,3),2)-2,1)) = Jaco(X2(3*X1_sorted(Esi(i,2),2)-1,1), X2(3*X1_sorted(Esi(i,3),2)-2,1))...
                + ncoeff*edgeEsi(i,2)/lenEsi(i,2)*dbern(phi(X1_sorted(Esi(i,2),2)), phi(X1_sorted(Esi(i,3),2)))*n(f3)...
                + ncoeff*edgeEsi(i,2)/lenEsi(i,2)*dbern(phi(X1_sorted(Esi(i,3),2)), phi(X1_sorted(Esi(i,2),2)))*n(f2);
            Jaco(X2(3*X1_sorted(Esi(i,2),2)-1,1), X2(3*X1_sorted(Esi(i,1),2)-2,1)) = Jaco(X2(3*X1_sorted(Esi(i,2),2)-1,1), X2(3*X1_sorted(Esi(i,1),2)-2,1))...
                + ncoeff*edgeEsi(i,1)/lenEsi(i,1)*dbern(phi(X1_sorted(Esi(i,2),2)), phi(X1_sorted(Esi(i,1),2)))*n(f1)...
                + ncoeff*edgeEsi(i,1)/lenEsi(i,1)*dbern(phi(X1_sorted(Esi(i,1),2)), phi(X1_sorted(Esi(i,2),2)))*n(f2);

            Res(X2(3*X1_sorted(Esi(i,2),2)-1,1),1) = Res(X2(3*X1_sorted(Esi(i,2),2)-1,1),1)...
                + ncoeff*Vt*edgeEsi(i,2)/lenEsi(i,2)*bern(phi(X1_sorted(Esi(i,2),2)),phi(X1_sorted(Esi(i,3),2)))*n(f3)...
                - ncoeff*Vt*edgeEsi(i,2)/lenEsi(i,2)*bern(phi(X1_sorted(Esi(i,3),2)),phi(X1_sorted(Esi(i,2),2)))*n(f2)...
                + ncoeff*Vt*edgeEsi(i,1)/lenEsi(i,1)*bern(phi(X1_sorted(Esi(i,2),2)),phi(X1_sorted(Esi(i,1),2)))*n(f1)...
                - ncoeff*Vt*edgeEsi(i,1)/lenEsi(i,1)*bern(phi(X1_sorted(Esi(i,1),2)),phi(X1_sorted(Esi(i,2),2)))*n(f2);

            %%% 3rd element
            % n part
            Jaco(X2(3*X1_sorted(Esi(i,3),2)-1,1), X2(3*X1_sorted(Esi(i,3),2)-1,1)) = Jaco(X2(3*X1_sorted(Esi(i,3),2)-1,1), X2(3*X1_sorted(Esi(i,3),2)-1,1))...
                - ncoeff*Vt*edgeEsi(i,3)/lenEsi(i,3)*bern(phi(X1_sorted(Esi(i,1),2)),phi(X1_sorted(Esi(i,3),2)))...
                - ncoeff*Vt*edgeEsi(i,2)/lenEsi(i,2)*bern(phi(X1_sorted(Esi(i,2),2)),phi(X1_sorted(Esi(i,3),2)));
            Jaco(X2(3*X1_sorted(Esi(i,3),2)-1,1), X2(3*X1_sorted(Esi(i,1),2)-1,1)) = Jaco(X2(3*X1_sorted(Esi(i,3),2)-1,1), X2(3*X1_sorted(Esi(i,1),2)-1,1))...
                + ncoeff*Vt*edgeEsi(i,3)/lenEsi(i,3)*bern(phi(X1_sorted(Esi(i,3),2)),phi(X1_sorted(Esi(i,1),2)));
            Jaco(X2(3*X1_sorted(Esi(i,3),2)-1,1), X2(3*X1_sorted(Esi(i,2),2)-1,1)) = Jaco(X2(3*X1_sorted(Esi(i,3),2)-1,1), X2(3*X1_sorted(Esi(i,2),2)-1,1))...
                + ncoeff*Vt*edgeEsi(i,2)/lenEsi(i,2)*bern(phi(X1_sorted(Esi(i,3),2)),phi(X1_sorted(Esi(i,2),2)));

            % phi part
            Jaco(X2(3*X1_sorted(Esi(i,3),2)-1,1), X2(3*X1_sorted(Esi(i,3),2)-2,1)) = Jaco(X2(3*X1_sorted(Esi(i,3),2)-1,1), X2(3*X1_sorted(Esi(i,3),2)-2,1))...
                - ncoeff*edgeEsi(i,3)/lenEsi(i,3)*dbern(phi(X1_sorted(Esi(i,3),2)), phi(X1_sorted(Esi(i,1),2)))*n(f1)...
                - ncoeff*edgeEsi(i,3)/lenEsi(i,3)*dbern(phi(X1_sorted(Esi(i,1),2)), phi(X1_sorted(Esi(i,3),2)))*n(f3)...
                - ncoeff*edgeEsi(i,2)/lenEsi(i,2)*dbern(phi(X1_sorted(Esi(i,3),2)), phi(X1_sorted(Esi(i,2),2)))*n(f2)...
                - ncoeff*edgeEsi(i,2)/lenEsi(i,2)*dbern(phi(X1_sorted(Esi(i,2),2)), phi(X1_sorted(Esi(i,3),2)))*n(f3);
            Jaco(X2(3*X1_sorted(Esi(i,3),2)-1,1), X2(3*X1_sorted(Esi(i,1),2)-2,1)) = Jaco(X2(3*X1_sorted(Esi(i,3),2)-1,1), X2(3*X1_sorted(Esi(i,1),2)-2,1))...
                + ncoeff*edgeEsi(i,3)/lenEsi(i,3)*dbern(phi(X1_sorted(Esi(i,3),2)), phi(X1_sorted(Esi(i,1),2)))*n(f1)...
                + ncoeff*edgeEsi(i,3)/lenEsi(i,3)*dbern(phi(X1_sorted(Esi(i,1),2)), phi(X1_sorted(Esi(i,3),2)))*n(f3);
            Jaco(X2(3*X1_sorted(Esi(i,3),2)-1,1), X2(3*X1_sorted(Esi(i,2),2)-2,1)) = Jaco(X2(3*X1_sorted(Esi(i,3),2)-1,1), X2(3*X1_sorted(Esi(i,2),2)-2,1))...
                + ncoeff*edgeEsi(i,2)/lenEsi(i,2)*dbern(phi(X1_sorted(Esi(i,3),2)), phi(X1_sorted(Esi(i,2),2)))*n(f2)...
                + ncoeff*edgeEsi(i,2)/lenEsi(i,2)*dbern(phi(X1_sorted(Esi(i,2),2)), phi(X1_sorted(Esi(i,3),2)))*n(f3);

            Res(X2(3*X1_sorted(Esi(i,3),2)-1,1),1) = Res(X2(3*X1_sorted(Esi(i,3),2)-1,1),1)...
                + ncoeff*Vt*edgeEsi(i,3)/lenEsi(i,3)*bern(phi(X1_sorted(Esi(i,3),2)),phi(X1_sorted(Esi(i,1),2)))*n(f1)...
                - ncoeff*Vt*edgeEsi(i,3)/lenEsi(i,3)*bern(phi(X1_sorted(Esi(i,1),2)),phi(X1_sorted(Esi(i,3),2)))*n(f3)...
                + ncoeff*Vt*edgeEsi(i,2)/lenEsi(i,2)*bern(phi(X1_sorted(Esi(i,3),2)),phi(X1_sorted(Esi(i,2),2)))*n(f2)...
                - ncoeff*Vt*edgeEsi(i,2)/lenEsi(i,2)*bern(phi(X1_sorted(Esi(i,2),2)),phi(X1_sorted(Esi(i,3),2)))*n(f3);
        end

        %%% hole DD
        for i=1:Esirow

            f1 = find(V_silicon(:,1)==Esi(i,1));
            f2 = find(V_silicon(:,1)==Esi(i,2));
            f3 = find(V_silicon(:,1)==Esi(i,3));

            %%% 1st element
            % p part
            Jaco(X2(3*X1_sorted(Esi(i,1),2),1), X2(3*X1_sorted(Esi(i,1),2),1)) = Jaco(X2(3*X1_sorted(Esi(i,1),2),1), X2(3*X1_sorted(Esi(i,1),2),1))...
                + pcoeff*Vt*edgeEsi(i,1)/lenEsi(i,1)*bern(phi(X1_sorted(Esi(i,1),2)),phi(X1_sorted(Esi(i,2),2)))...
                + pcoeff*Vt*edgeEsi(i,3)/lenEsi(i,3)*bern(phi(X1_sorted(Esi(i,1),2)),phi(X1_sorted(Esi(i,3),2)));
            Jaco(X2(3*X1_sorted(Esi(i,1),2),1), X2(3*X1_sorted(Esi(i,2),2),1)) = Jaco(X2(3*X1_sorted(Esi(i,1),2),1), X2(3*X1_sorted(Esi(i,2),2),1))...
                - pcoeff*Vt*edgeEsi(i,1)/lenEsi(i,1)*bern(phi(X1_sorted(Esi(i,2),2)),phi(X1_sorted(Esi(i,1),2)));
            Jaco(X2(3*X1_sorted(Esi(i,1),2),1), X2(3*X1_sorted(Esi(i,3),2),1)) = Jaco(X2(3*X1_sorted(Esi(i,1),2),1), X2(3*X1_sorted(Esi(i,3),2),1))...
                - pcoeff*Vt*edgeEsi(i,3)/lenEsi(i,3)*bern(phi(X1_sorted(Esi(i,3),2)),phi(X1_sorted(Esi(i,1),2)));

            % phi part
            Jaco(X2(3*X1_sorted(Esi(i,1),2),1), X2(3*X1_sorted(Esi(i,1),2)-2,1)) = Jaco(X2(3*X1_sorted(Esi(i,1),2),1), X2(3*X1_sorted(Esi(i,1),2)-2,1))...
                - pcoeff*edgeEsi(i,1)/lenEsi(i,1)*dbern(phi(X1_sorted(Esi(i,2),2)),phi(X1_sorted(Esi(i,1),2)))*p(f2)...
                - pcoeff*edgeEsi(i,1)/lenEsi(i,1)*dbern(phi(X1_sorted(Esi(i,1),2)),phi(X1_sorted(Esi(i,2),2)))*p(f1)...
                - pcoeff*edgeEsi(i,3)/lenEsi(i,3)*dbern(phi(X1_sorted(Esi(i,3),2)),phi(X1_sorted(Esi(i,1),2)))*p(f3)...
                - pcoeff*edgeEsi(i,3)/lenEsi(i,3)*dbern(phi(X1_sorted(Esi(i,1),2)),phi(X1_sorted(Esi(i,3),2)))*p(f1);
            Jaco(X2(3*X1_sorted(Esi(i,1),2),1), X2(3*X1_sorted(Esi(i,2),2)-2,1)) = Jaco(X2(3*X1_sorted(Esi(i,1),2),1), X2(3*X1_sorted(Esi(i,2),2)-2,1))...
                + pcoeff*edgeEsi(i,1)/lenEsi(i,1)*dbern(phi(X1_sorted(Esi(i,2),2)),phi(X1_sorted(Esi(i,1),2)))*p(f2)...
                + pcoeff*edgeEsi(i,1)/lenEsi(i,1)*dbern(phi(X1_sorted(Esi(i,1),2)),phi(X1_sorted(Esi(i,2),2)))*p(f1);
            Jaco(X2(3*X1_sorted(Esi(i,1),2),1), X2(3*X1_sorted(Esi(i,3),2)-2,1)) = Jaco(X2(3*X1_sorted(Esi(i,1),2),1), X2(3*X1_sorted(Esi(i,3),2)-2,1))...
                + pcoeff*edgeEsi(i,3)/lenEsi(i,3)*dbern(phi(X1_sorted(Esi(i,3),2)),phi(X1_sorted(Esi(i,1),2)))*p(f3)...
                + pcoeff*edgeEsi(i,3)/lenEsi(i,3)*dbern(phi(X1_sorted(Esi(i,1),2)),phi(X1_sorted(Esi(i,3),2)))*p(f1);

            Res(X2(3*X1_sorted(Esi(i,1),2),1),1) = Res(X2(3*X1_sorted(Esi(i,1),2),1),1)...
                - pcoeff*Vt*edgeEsi(i,1)/lenEsi(i,1)*bern(phi(X1_sorted(Esi(i,2),2)),phi(X1_sorted(Esi(i,1),2)))*p(f2)...
                + pcoeff*Vt*edgeEsi(i,1)/lenEsi(i,1)*bern(phi(X1_sorted(Esi(i,1),2)),phi(X1_sorted(Esi(i,2),2)))*p(f1)...
                - pcoeff*Vt*edgeEsi(i,3)/lenEsi(i,3)*bern(phi(X1_sorted(Esi(i,3),2)),phi(X1_sorted(Esi(i,1),2)))*p(f3)...
                + pcoeff*Vt*edgeEsi(i,3)/lenEsi(i,3)*bern(phi(X1_sorted(Esi(i,1),2)),phi(X1_sorted(Esi(i,3),2)))*p(f1);

            %%% 2nd element
            % p part
            Jaco(X2(3*X1_sorted(Esi(i,2),2),1), X2(3*X1_sorted(Esi(i,2),2),1)) = Jaco(X2(3*X1_sorted(Esi(i,2),2),1), X2(3*X1_sorted(Esi(i,2),2),1))...
                + pcoeff*Vt*edgeEsi(i,2)/lenEsi(i,2)*bern(phi(X1_sorted(Esi(i,2),2)),phi(X1_sorted(Esi(i,3),2)))...
                + pcoeff*Vt*edgeEsi(i,1)/lenEsi(i,1)*bern(phi(X1_sorted(Esi(i,2),2)),phi(X1_sorted(Esi(i,1),2)));
            Jaco(X2(3*X1_sorted(Esi(i,2),2),1), X2(3*X1_sorted(Esi(i,3),2),1)) = Jaco(X2(3*X1_sorted(Esi(i,2),2),1), X2(3*X1_sorted(Esi(i,3),2),1))...
                - pcoeff*Vt*edgeEsi(i,2)/lenEsi(i,2)*bern(phi(X1_sorted(Esi(i,3),2)),phi(X1_sorted(Esi(i,2),2)));
            Jaco(X2(3*X1_sorted(Esi(i,2),2),1), X2(3*X1_sorted(Esi(i,1),2),1)) 
