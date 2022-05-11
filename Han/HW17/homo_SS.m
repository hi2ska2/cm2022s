clc
clear
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% homogeneous %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
width = 1e-6;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
q = 1.6021e-19;
e0 = 8.8542e-12;
esi = 11.7;
nint = 1e16;
Kb = 1.38e-23;
Temp = 300;
Vt = Kb*Temp/q;
mun = 518e-4;
mup = 250e-4;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Doping density
Nd = 2e23;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Iteration number
Niter1 = 20;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
V = importdata("vertex.txt"); Vrow = size(V,1);
E = importdata("element.txt"); Erow = size(E,1);
E2 = importdata("element_region2.txt"); E2row = size(E2,1);
E3 = importdata("element_region3.txt"); E3row = size(E3,1);
E4 = importdata("element_region4.txt"); E4row = size(E4,1);
Esi = vertcat(E2, E3, E4); Esirow = size(Esi,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% edge contact information
contact = importdata("contact.txt");
[conrow, concol] = size(contact);
unicon = unique(contact);
Ncon1 = 10;
Ncon2 = 10;

% Anode contact analysis
conin = zeros(2*Ncon2,1);
ccc=1;
for i=Ncon1+1:Ncon1+Ncon2
    for j=1:concol
        conin(ccc,1) = contact(i,j);
        ccc = ccc+1;
    end
end
uniconin = unique(conin);

% number of variable
Nvariable = 3;

% number of vertex
VE = unique(E); VErow= size(VE,1);
VE2 = unique(E2); VE2row = size(VE2, 1);
VE3 = unique(E3); VE3row = size(VE3, 1);
VE4 = unique(E4); VE4row = size(VE4, 1);

V_silicon = zeros(size(unique(Esi),1),2);
V_silicon_row = size(V_silicon, 1);
V_silicon(:,1) = unique(Esi);

%%% initial potential indexing
X1 = zeros(V_silicon_row,2);
X1(:,1) = V_silicon(:,1);
X1(:,2) = V_silicon(:,1);

X1_sorted = X1;

%%% Initial potential re-indexing
V_silicon(:,2) = X1(:,2);

%%% Total Variable indexing
X2 = zeros(3*(size(X1,1)),1);
count2 = 0;
index2 = 0;
for ivertex = 1:V_silicon_row
    for ivariable = 1:3
        count2 = count2 + 1;
        index2 = index2 + 1;
        X2(count2,1) = index2;
    end
end
X2_sorted = X2;

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

%%% edge per region
edgeE2=zeros(E2row,3); lenE2=zeros(E2row,3);
edgeE3=zeros(E3row,3); lenE3=zeros(E3row,3);
edgeE4=zeros(E4row,3); lenE4=zeros(E4row,3);

bb=1; cc=1; dd=1;

for i=1:Erow
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
end

edgeEsi = vertcat(edgeE2, edgeE3, edgeE4);
lenEsi = vertcat(lenE2, lenE3, lenE4);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
phi = zeros(size(X1,1),1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Initial potential_edit by euqation : Vt*ln(Ndop/nint)
for i=1:VE2row
    phi(X1_sorted(VE2(i),2),1) = Vt*log(Nd/nint);
end

for i=1:VE3row
    phi(X1_sorted(VE3(i),2),1) = Vt*log(Nd/nint);
end

for i=1:VE4row
    phi(X1_sorted(VE4(i),2),1) = Vt*log(Nd/nint);
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
%%% jacobian, residue matrix (non-linear Possion equation)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
save_phi_upd = zeros(Niter1,1);
save_upd = zeros(size(X2_sorted,1),1);

for it=1:Niter1

    jac = sparse(size(X2_sorted,1),size(X2_sorted,1));
    res = zeros(size(X2_sorted,1),1);

    saveexl = zeros(size(X1,1),1);

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
            saveexl(X1_sorted(Esi(i,1),2)) = saveexl(X1_sorted(Esi(i,1),2)) + (edgeEsi(i,1)*lenEsi(i,1) + edgeEsi(i,3)*lenEsi(i,3))/4;

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
            saveexl(X1_sorted(Esi(i,2),2)) = saveexl(X1_sorted(Esi(i,2),2)) + (edgeEsi(i,1)*lenEsi(i,1) + edgeEsi(i,2)*lenEsi(i,2))/4;

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
            saveexl(X1_sorted(Esi(i,3),2)) = saveexl(X1_sorted(Esi(i,3),2)) + (edgeEsi(i,2)*lenEsi(i,2) + edgeEsi(i,3)*lenEsi(i,3))/4;

        elseif i > E2row && i <= (E2row + E3row)
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
            saveexl(X1_sorted(Esi(i,1),2)) = saveexl(X1_sorted(Esi(i,1),2)) + (edgeEsi(i,1)*lenEsi(i,1) + edgeEsi(i,3)*lenEsi(i,3))/4;

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
            saveexl(X1_sorted(Esi(i,2),2)) = saveexl(X1_sorted(Esi(i,2),2)) + (edgeEsi(i,1)*lenEsi(i,1) + edgeEsi(i,2)*lenEsi(i,2))/4;

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
            saveexl(X1_sorted(Esi(i,3),2)) = saveexl(X1_sorted(Esi(i,3),2)) + (edgeEsi(i,2)*lenEsi(i,2) + edgeEsi(i,3)*lenEsi(i,3))/4;

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
            saveexl(X1_sorted(Esi(i,1),2)) = saveexl(X1_sorted(Esi(i,1),2)) + (edgeEsi(i,1)*lenEsi(i,1) + edgeEsi(i,3)*lenEsi(i,3))/4;

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
            saveexl(X1_sorted(Esi(i,2),2)) = saveexl(X1_sorted(Esi(i,2),2)) + (edgeEsi(i,1)*lenEsi(i,1) + edgeEsi(i,2)*lenEsi(i,2))/4;

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
            saveexl(X1_sorted(Esi(i,3),2)) = saveexl(X1_sorted(Esi(i,3),2)) + (edgeEsi(i,2)*lenEsi(i,2) + edgeEsi(i,3)*lenEsi(i,3))/4;
        end
    end

    %%%  for electron
    for i=1:V_silicon_row
        %%% n-nit*exp(phi/Vt)
        jac(X2(3*V_silicon(i,2)-1),X2(3*V_silicon(i,2)-1)) = 1;
        jac(X2(3*V_silicon(i,2)-1),X2(3*V_silicon(i,2)-2)) = - nint/Vt*exp(phi(V_silicon(i,2))/Vt);
        res(X2(3*V_silicon(i,2)-1),1) = n(i) - nint*exp(phi(V_silicon(i,2))/Vt);
        %%% - qn
        jac(X2(3*V_silicon(i,2)-2),X2(3*V_silicon(i,2)-1)) = - saveexl(V_silicon(i,2))*q/e0;
        res(X2(3*V_silicon(i,2)-2),1) = res(X2(3*V_silicon(i,2)-2),1) - saveexl(V_silicon(i,2))*q/e0*n(i);
    end

    %%%  for hole
    for i=1:V_silicon_row
        %%% p-nit*exp(-phi/Vt)
        jac(X2(3*V_silicon(i,2)),X2(3*V_silicon(i,2))) = 1;
        jac(X2(3*V_silicon(i,2)),X2(3*V_silicon(i,2)-2)) = nint/Vt*exp(-phi(V_silicon(i,2))/Vt);
        res(X2(3*V_silicon(i,2)),1) = p(i) - nint*exp(-phi(V_silicon(i,2))/Vt);
        %%% + qp
        jac(X2(3*V_silicon(i,2)-2),X2(3*V_silicon(i,2))) = saveexl(V_silicon(i,2))*q/e0;
        res(X2(3*V_silicon(i,2)-2),1) = res(X2(3*V_silicon(i,2)-2),1) + saveexl(V_silicon(i,2))*q/e0*p(i);
    end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Dirichlet boundary condition for 1st iteration

    % ground & input
    for i=1:(Ncon1+Ncon2)
        for j=1:concol
            fcon = X1_sorted(contact(i,j),2);
            fcon_np = find(contact(i,j)==V_silicon(:,1));

            jac(X2(3*fcon-2),:) = 0;
            jac(X2(3*fcon-1),:) = 0;
            jac(X2(3*fcon),:) = 0;

            jac(X2(3*fcon-2), X2(3*fcon-2)) = 1;
            jac(X2(3*fcon-1), X2(3*fcon-1)) = 1;
            jac(X2(3*fcon), X2(3*fcon)) = 1;

            % Cathode
            if i<=Ncon1
                res(X2(3*fcon-2),1) = phi(fcon) - Vt*log(Nd/nint);
                res(X2(3*fcon-1),1) = n(fcon_np) - Nd;
                res(X2(3*fcon),1) = p(fcon_np) - nint*nint/Nd;
                % Anode
            elseif i>(Ncon1) &&  i<=(Ncon1+Ncon2)
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

save_phi_update = zeros(Niter2,1);
save_update = zeros(size(X2_sorted,1),1);

for step=1:51

    Vex = (step-1)*0.01;

    for it=1:Niter2

        Jaco = sparse(size(X2_sorted,1),size(X2_sorted,1));
        Res = zeros(size(X2_sorted,1),1);

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
        %%% Dirichlet boundary condition for final iteration
        for i=1:(Ncon1+Ncon2)
            for j=1:concol
                fcon = X1_sorted(contact(i,j),2);
                fcon_np = find(contact(i,j)==V_silicon(:,1));

                Jaco(X2(3*fcon-2),:) = 0;
                Jaco(X2(3*fcon-1),:) = 0;
                Jaco(X2(3*fcon),:) = 0;
                Jaco(X2(3*fcon-2), X2(3*fcon-2)) = 1;
                Jaco(X2(3*fcon-1), X2(3*fcon-1)) = 1;
                Jaco(X2(3*fcon), X2(3*fcon)) = 1;

                % cathode
                if i<=Ncon1
                    Res(X2(3*fcon-2),1) = phi(fcon) - Vt*log(Nd/nint);
                    Res(X2(3*fcon-1),1) = n(fcon_np) - Nd;
                    Res(X2(3*fcon),1) = p(fcon_np) - nint*nint/Nd;
                    % Anode
                elseif i>(Ncon1) &&  i<=(Ncon1+Ncon2)
                    Res(X2(3*fcon-2),1) = phi(fcon) - Vt*log(Nd/nint)-Vex;
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

        save_update(:,it) = update(:,1);
        save_phi_update(it,1) = maxupdate;

        if maxupdate < 1e-10
            break
        end
    end
end
n_DC = n(:,1);
p_DC = p(:,1);
phi_DC = phi(:,1);
