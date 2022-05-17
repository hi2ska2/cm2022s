clc
clear
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% homogeneous SS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% jacobian, residue matrix (Green-function)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

jacA = sparse(size(X2_sorted,1),size(X2_sorted,1));

saveexl = zeros(size(X1,1),1);

for i=1:Esirow
    if i <= E2row
        jacA(X2(3*(X1_sorted(Esi(i,1),2))-2,1),X2(3*(X1_sorted(Esi(i,1),2))-2,1)) = jacA(X2(3*(X1_sorted(Esi(i,1),2))-2,1),X2(3*(X1_sorted(Esi(i,1),2))-2,1))...
            - esi*(edgeEsi(i,1)/lenEsi(i,1) + edgeEsi(i,3)/lenEsi(i,3));
        jacA(X2(3*(X1_sorted(Esi(i,1),2))-2,1),X2(3*(X1_sorted(Esi(i,2),2))-2,1)) = jacA(X2(3*(X1_sorted(Esi(i,1),2))-2,1),X2(3*(X1_sorted(Esi(i,2),2))-2,1))...
            + esi*edgeEsi(i,1)/lenEsi(i,1);
        jacA(X2(3*(X1_sorted(Esi(i,1),2))-2,1),X2(3*(X1_sorted(Esi(i,3),2))-2,1)) = jacA(X2(3*(X1_sorted(Esi(i,1),2))-2,1),X2(3*(X1_sorted(Esi(i,3),2))-2,1))...
            + esi*edgeEsi(i,3)/lenEsi(i,3);

        jacA(X2(3*(X1_sorted(Esi(i,2),2))-2,1),X2(3*(X1_sorted(Esi(i,2),2))-2,1)) = jacA(X2(3*(X1_sorted(Esi(i,2),2))-2,1),X2(3*(X1_sorted(Esi(i,2),2))-2,1))...
            - esi*(edgeEsi(i,1)/lenEsi(i,1) + edgeEsi(i,2)/lenEsi(i,2));
        jacA(X2(3*(X1_sorted(Esi(i,2),2))-2,1),X2(3*(X1_sorted(Esi(i,1),2))-2,1)) = jacA(X2(3*(X1_sorted(Esi(i,2),2))-2,1),X2(3*(X1_sorted(Esi(i,1),2))-2,1))...
            + esi*edgeEsi(i,1)/lenEsi(i,1);
        jacA(X2(3*(X1_sorted(Esi(i,2),2))-2,1),X2(3*(X1_sorted(Esi(i,3),2))-2,1)) = jacA(X2(3*(X1_sorted(Esi(i,2),2))-2,1),X2(3*(X1_sorted(Esi(i,3),2))-2,1))...
            + esi*edgeEsi(i,2)/lenEsi(i,2);

        jacA(X2(3*(X1_sorted(Esi(i,3),2))-2,1),X2(3*(X1_sorted(Esi(i,3),2))-2,1)) = jacA(X2(3*(X1_sorted(Esi(i,3),2))-2,1),X2(3*(X1_sorted(Esi(i,3),2))-2,1))...
            - esi*(edgeEsi(i,2)/lenEsi(i,2) + edgeEsi(i,3)/lenEsi(i,3));
        jacA(X2(3*(X1_sorted(Esi(i,3),2))-2,1),X2(3*(X1_sorted(Esi(i,2),2))-2,1)) = jacA(X2(3*(X1_sorted(Esi(i,3),2))-2,1),X2(3*(X1_sorted(Esi(i,2),2))-2,1))...
            + esi*edgeEsi(i,2)/lenEsi(i,2);
        jacA(X2(3*(X1_sorted(Esi(i,3),2))-2,1),X2(3*(X1_sorted(Esi(i,1),2))-2,1)) = jacA(X2(3*(X1_sorted(Esi(i,3),2))-2,1),X2(3*(X1_sorted(Esi(i,1),2))-2,1))...
            + esi*edgeEsi(i,3)/lenEsi(i,3);

    elseif i > E2row && i <= (E2row + E3row)
        jacA(X2(3*(X1_sorted(Esi(i,1),2))-2,1),X2(3*(X1_sorted(Esi(i,1),2))-2,1)) = jacA(X2(3*(X1_sorted(Esi(i,1),2))-2,1),X2(3*(X1_sorted(Esi(i,1),2))-2,1))...
            - esi*(edgeEsi(i,1)/lenEsi(i,1) + edgeEsi(i,3)/lenEsi(i,3));
        jacA(X2(3*(X1_sorted(Esi(i,1),2))-2,1),X2(3*(X1_sorted(Esi(i,2),2))-2,1)) = jacA(X2(3*(X1_sorted(Esi(i,1),2))-2,1),X2(3*(X1_sorted(Esi(i,2),2))-2,1))...
            + esi*edgeEsi(i,1)/lenEsi(i,1);
        jacA(X2(3*(X1_sorted(Esi(i,1),2))-2,1),X2(3*(X1_sorted(Esi(i,3),2))-2,1)) = jacA(X2(3*(X1_sorted(Esi(i,1),2))-2,1),X2(3*(X1_sorted(Esi(i,3),2))-2,1))...
            + esi*edgeEsi(i,3)/lenEsi(i,3);

        jacA(X2(3*(X1_sorted(Esi(i,2),2))-2,1),X2(3*(X1_sorted(Esi(i,2),2))-2,1)) = jacA(X2(3*(X1_sorted(Esi(i,2),2))-2,1),X2(3*(X1_sorted(Esi(i,2),2))-2,1))...
            - esi*(edgeEsi(i,1)/lenEsi(i,1) + edgeEsi(i,2)/lenEsi(i,2));
        jacA(X2(3*(X1_sorted(Esi(i,2),2))-2,1),X2(3*(X1_sorted(Esi(i,1),2))-2,1)) = jacA(X2(3*(X1_sorted(Esi(i,2),2))-2,1),X2(3*(X1_sorted(Esi(i,1),2))-2,1))...
            + esi*edgeEsi(i,1)/lenEsi(i,1);
        jacA(X2(3*(X1_sorted(Esi(i,2),2))-2,1),X2(3*(X1_sorted(Esi(i,3),2))-2,1)) = jacA(X2(3*(X1_sorted(Esi(i,2),2))-2,1),X2(3*(X1_sorted(Esi(i,3),2))-2,1))...
            + esi*edgeEsi(i,2)/lenEsi(i,2);

        jacA(X2(3*(X1_sorted(Esi(i,3),2))-2,1),X2(3*(X1_sorted(Esi(i,3),2))-2,1)) = jacA(X2(3*(X1_sorted(Esi(i,3),2))-2,1),X2(3*(X1_sorted(Esi(i,3),2))-2,1))...
            - esi*(edgeEsi(i,2)/lenEsi(i,2) + edgeEsi(i,3)/lenEsi(i,3));
        jacA(X2(3*(X1_sorted(Esi(i,3),2))-2,1),X2(3*(X1_sorted(Esi(i,2),2))-2,1)) = jacA(X2(3*(X1_sorted(Esi(i,3),2))-2,1),X2(3*(X1_sorted(Esi(i,2),2))-2,1))...
            + esi*edgeEsi(i,2)/lenEsi(i,2);
        jacA(X2(3*(X1_sorted(Esi(i,3),2))-2,1),X2(3*(X1_sorted(Esi(i,1),2))-2,1)) = jacA(X2(3*(X1_sorted(Esi(i,3),2))-2,1),X2(3*(X1_sorted(Esi(i,1),2))-2,1))...
            + esi*edgeEsi(i,3)/lenEsi(i,3);

    elseif i > (E2row+E3row) && i <= Esirow
        jacA(X2(3*(X1_sorted(Esi(i,1),2))-2,1),X2(3*(X1_sorted(Esi(i,1),2))-2,1)) = jacA(X2(3*(X1_sorted(Esi(i,1),2))-2,1),X2(3*(X1_sorted(Esi(i,1),2))-2,1))...
            - esi*(edgeEsi(i,1)/lenEsi(i,1) + edgeEsi(i,3)/lenEsi(i,3));
        jacA(X2(3*(X1_sorted(Esi(i,1),2))-2,1),X2(3*(X1_sorted(Esi(i,2),2))-2,1)) = jacA(X2(3*(X1_sorted(Esi(i,1),2))-2,1),X2(3*(X1_sorted(Esi(i,2),2))-2,1))...
            + esi*edgeEsi(i,1)/lenEsi(i,1);
        jacA(X2(3*(X1_sorted(Esi(i,1),2))-2,1),X2(3*(X1_sorted(Esi(i,3),2))-2,1)) = jacA(X2(3*(X1_sorted(Esi(i,1),2))-2,1),X2(3*(X1_sorted(Esi(i,3),2))-2,1))...
            + esi*edgeEsi(i,3)/lenEsi(i,3);

        jacA(X2(3*(X1_sorted(Esi(i,2),2))-2,1),X2(3*(X1_sorted(Esi(i,2),2))-2,1)) = jacA(X2(3*(X1_sorted(Esi(i,2),2))-2,1),X2(3*(X1_sorted(Esi(i,2),2))-2,1))...
            - esi*(edgeEsi(i,1)/lenEsi(i,1) + edgeEsi(i,2)/lenEsi(i,2));
        jacA(X2(3*(X1_sorted(Esi(i,2),2))-2,1),X2(3*(X1_sorted(Esi(i,1),2))-2,1)) = jacA(X2(3*(X1_sorted(Esi(i,2),2))-2,1),X2(3*(X1_sorted(Esi(i,1),2))-2,1))...
            + esi*edgeEsi(i,1)/lenEsi(i,1);
        jacA(X2(3*(X1_sorted(Esi(i,2),2))-2,1),X2(3*(X1_sorted(Esi(i,3),2))-2,1)) = jacA(X2(3*(X1_sorted(Esi(i,2),2))-2,1),X2(3*(X1_sorted(Esi(i,3),2))-2,1))...
            + esi*edgeEsi(i,2)/lenEsi(i,2);

        jacA(X2(3*X1_sorted(Esi(i,3),2)-2,1),X2(3*(X1_sorted(Esi(i,3),2))-2,1)) = jacA(X2(3*(X1_sorted(Esi(i,3),2))-2,1),X2(3*(X1_sorted(Esi(i,3),2))-2,1))...
            - esi*(edgeEsi(i,2)/lenEsi(i,2) + edgeEsi(i,3)/lenEsi(i,3));
        jacA(X2(3*X1_sorted(Esi(i,3),2)-2,1),X2(3*(X1_sorted(Esi(i,2),2))-2,1)) = jacA(X2(3*(X1_sorted(Esi(i,3),2))-2,1),X2(3*(X1_sorted(Esi(i,2),2))-2,1))...
            + esi*edgeEsi(i,2)/lenEsi(i,2);
        jacA(X2(3*X1_sorted(Esi(i,3),2)-2,1),X2(3*(X1_sorted(Esi(i,1),2))-2,1)) = jacA(X2(3*(X1_sorted(Esi(i,3),2))-2,1),X2(3*(X1_sorted(Esi(i,1),2))-2,1))...
            + esi*edgeEsi(i,3)/lenEsi(i,3);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Dirichlet boundary condition for 1st iteration
resb = zeros(size(X2_sorted,1),1);
for i=1:(Ncon1+Ncon2)
    for j=1:concol
        fcon = X1_sorted(contact(i,j),2);
        fcon_np = find(contact(i,j)==V_silicon(:,1));

        jacA(X2(3*fcon-2),:) = 0;
        jacA(X2(3*fcon-2), X2(3*fcon-2)) = 1;

        if i<=Ncon1 % Cathode
            resb(X2(3*fcon-2),1) = 0;

        elseif i>(Ncon1) &&  i<=(Ncon1+Ncon2) % Anode
            resb(X2(3*fcon-2),1) = 0;

        end

    end
end

dirac = 316;
for i=1:45
    resb(3*dirac-2,1) = 1;
    dirac = dirac + 1;
end

jacA2 = jacA(1:3:size(jacA,1),1:3:size(jacA,1));
resb2 = resb(1:3:size(resb,1),1);

phi = jacA2 \ resb2;

patch('Faces',E,'Vertices',V, 'FaceVertexCData',phi, 'EdgeColor','k','FaceColor','interp','LineWidth',1, 'Marker','o');
patch('Faces',E2,'Vertices',V, 'FaceVertexCData',phi,'EdgeColor','black','FaceColor','interp','LineWidth',1, 'Marker','o');
patch('Faces',E3,'Vertices',V, 'FaceVertexCData',phi,'EdgeColor','black','FaceColor','interp','LineWidth',1, 'Marker','o');
patch('Faces',E4,'Vertices',V, 'FaceVertexCData',phi,'EdgeColor','black','FaceColor','interp','LineWidth',1, 'Marker','o');
colorbar
