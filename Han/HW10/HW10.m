%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% HW10
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Define
q = 1.6e-19;
e0=8.85e-12;
nint = 1e16;
Ndop = 1e24;
Kb=1.38*10^-23; T=300; Vt=(Kb*T)/q;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
V = importdata("vertex.txt"); Vrow = size(V,1);
F = importdata("element.txt"); Frow = size(F,1);
F1 = importdata("element_region1.txt"); F1row = size(F1,1);
F2 = importdata("element_region2.txt"); F2row = size(F2,1);
F3 = importdata("element_region3.txt"); F3row = size(F3,1);
contact = importdata("contact.txt"); [conrow, concol] = size(contact);
unicon = unique(contact);
Ntcon = 5; % number of top contact edge
tcon = 0.33374;
Nbcon = 5; % bottom contact edge
bcon = 0.33374;
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
aaa_VF1 = zeros(VF1row,2);
aa=1;
for i=1:VF1row
    aaa_VF1(i,1)=VF1(i);
    aaa_VF1(i,2)=aa;
    aa=aa+1;
end
aaa_VF2 = zeros(VF2row,2);
bb=VF1row+1;
for i=1:VF2row
    aaa_VF2(i,1)=VF2(i);
    aaa_VF2(i,2)=bb;
    bb=bb+1;
end
aaa_VF3 = zeros(VF3row,2);
cc=VF1row+VF2row+1;
for i=1:VF3row
    aaa_VF3(i,1)=VF3(i);
    aaa_VF3(i,2)=cc;
    cc=cc+1;
end
bb1_VF2 = zeros(VF2row,1);
bb1 = VF1row+VF2row+VF3row+1;
for i=1:VF2row
    bb1_VF2(i,1)=VF2(i);
    bb1_VF2(i,2)=bb1;
    bb1=bb1+1;
end
bb2_VF2 = zeros(VF2row,1);
bb2 = VF1row+VF2row+VF3row+VF2row+1;
for i=1:VF2row
    bb2_VF2(i,1)=VF2(i);
    bb2_VF2(i,2)=bb2;
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
A = zeros(VFTrow,VFTrow);
b = zeros(VFTrow,1);
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% initial value calculation
for i=1:F1row
    for j=1:3
        fa = find(F1(i,1)==aaa_VF1(:,1));
        fb = find(F1(i,2)==aaa_VF1(:,1));
        fc = find(F1(i,3)==aaa_VF1(:,1));
        if j==1
            A(aaa_VF1(fa,2),aaa_VF1(fa,2)) = A(aaa_VF1(fa,2),aaa_VF1(fa,2)) - e1*edgeF1(i,1)/lenF1(i,1) - e1*edgeF1(i,3)/lenF1(i,3);
            A(aaa_VF1(fa,2),aaa_VF1(fb,2)) = A(aaa_VF1(fa,2),aaa_VF1(fb,2)) + e1*edgeF1(i,1)/lenF1(i,1);
            A(aaa_VF1(fa,2),aaa_VF1(fc,2)) = A(aaa_VF1(fa,2),aaa_VF1(fc,2)) + e1*edgeF1(i,3)/lenF1(i,3);
        elseif j==2
            A(aaa_VF1(fb,2),aaa_VF1(fb,2)) = A(aaa_VF1(fb,2),aaa_VF1(fb,2)) - e1*edgeF1(i,2)/lenF1(i,2) - e1*edgeF1(i,1)/lenF1(i,1);
            A(aaa_VF1(fb,2),aaa_VF1(fc,2)) = A(aaa_VF1(fb,2),aaa_VF1(fc,2)) + e1*edgeF1(i,2)/lenF1(i,2);
            A(aaa_VF1(fb,2),aaa_VF1(fa,2)) = A(aaa_VF1(fb,2),aaa_VF1(fa,2)) + e1*edgeF1(i,1)/lenF1(i,1);
        elseif j==3
            A(aaa_VF1(fc,2),aaa_VF1(fc,2)) = A(aaa_VF1(fc,2),aaa_VF1(fc,2)) - e1*edgeF1(i,3)/lenF1(i,3) - e1*edgeF1(i,2)/lenF1(i,2);
            A(aaa_VF1(fc,2),aaa_VF1(fa,2)) = A(aaa_VF1(fc,2),aaa_VF1(fa,2)) + e1*edgeF1(i,3)/lenF1(i,3);
            A(aaa_VF1(fc,2),aaa_VF1(fb,2)) = A(aaa_VF1(fc,2),aaa_VF1(fb,2)) + e1*edgeF1(i,2)/lenF1(i,2);
        end
    end
end

for i=1:F2row
    for j=1:3
        fa = find(F2(i,1)==aaa_VF2(:,1));
        fb = find(F2(i,2)==aaa_VF2(:,1));
        fc = find(F2(i,3)==aaa_VF2(:,1));
        if j==1
            A(aaa_VF2(fa,2),aaa_VF2(fa,2)) = A(aaa_VF2(fa,2),aaa_VF2(fa,2)) - (e2*edgeF2(i,1)/lenF2(i,1) + e2*edgeF2(i,3)/lenF2(i,3));
            A(aaa_VF2(fa,2),aaa_VF2(fb,2)) = A(aaa_VF2(fa,2),aaa_VF2(fb,2)) + e2*edgeF2(i,1)/lenF2(i,1);
            A(aaa_VF2(fa,2),aaa_VF2(fc,2)) = A(aaa_VF2(fa,2),aaa_VF2(fc,2)) + e2*edgeF2(i,3)/lenF2(i,3);
            b(aaa_VF2(fa,2),1) = b(aaa_VF2(fa,2),1) + (edgeF2(i,1)*lenF2(i,1) + edgeF2(i,3)*lenF2(i,3))/4*q/e0*Ndop;
        elseif j==2
            A(aaa_VF2(fb,2),aaa_VF2(fb,2)) = A(aaa_VF2(fb,2),aaa_VF2(fb,2)) - (e2*edgeF2(i,2)/lenF2(i,2) + e2*edgeF2(i,1)/lenF2(i,1));
            A(aaa_VF2(fb,2),aaa_VF2(fc,2)) = A(aaa_VF2(fb,2),aaa_VF2(fc,2)) + e2*edgeF2(i,2)/lenF2(i,2);
            A(aaa_VF2(fb,2),aaa_VF2(fa,2)) = A(aaa_VF2(fb,2),aaa_VF2(fa,2)) + e2*edgeF2(i,1)/lenF2(i,1);
            b(aaa_VF2(fb,2),1) = b(aaa_VF2(fb,2),1) + (edgeF2(i,2)*lenF2(i,2) + edgeF2(i,1)*lenF2(i,1))/4*q/e0*Ndop;
        elseif j==3
            A(aaa_VF2(fc,2),aaa_VF2(fc,2)) = A(aaa_VF2(fc,2),aaa_VF2(fc,2)) - (e2*edgeF2(i,3)/lenF2(i,3) + e2*edgeF2(i,2)/lenF2(i,2));
            A(aaa_VF2(fc,2),aaa_VF2(fa,2)) = A(aaa_VF2(fc,2),aaa_VF2(fa,2)) + e2*edgeF2(i,3)/lenF2(i,3);
            A(aaa_VF2(fc,2),aaa_VF2(fb,2)) = A(aaa_VF2(fc,2),aaa_VF2(fb,2)) + e2*edgeF2(i,2)/lenF2(i,2);
            b(aaa_VF2(fc,2),1) = b(aaa_VF2(fc,2),1) + (edgeF2(i,3)*lenF2(i,3) + edgeF2(i,2)*lenF2(i,2))/4*q/e0*Ndop;
        end
    end
end

for i=1:F3row
    for j=1:3
        fa = find(F3(i,1)==aaa_VF3(:,1));
        fb = find(F3(i,2)==aaa_VF3(:,1));
        fc = find(F3(i,3)==aaa_VF3(:,1));
        if j==1
            A(aaa_VF3(fa,2),aaa_VF3(fa,2)) = A(aaa_VF3(fa,2),aaa_VF3(fa,2)) - (e3*edgeF3(i,1)/lenF3(i,1) + e3*edgeF3(i,3)/lenF3(i,3));
            A(aaa_VF3(fa,2),aaa_VF3(fb,2)) = A(aaa_VF3(fa,2),aaa_VF3(fb,2)) + e3*edgeF3(i,1)/lenF3(i,1);
            A(aaa_VF3(fa,2),aaa_VF3(fc,2)) = A(aaa_VF3(fa,2),aaa_VF3(fc,2)) + e3*edgeF3(i,3)/lenF3(i,3);
        elseif j==2
            A(aaa_VF3(fb,2),aaa_VF3(fb,2)) = A(aaa_VF3(fb,2),aaa_VF3(fb,2)) - (e3*edgeF3(i,2)/lenF3(i,2) + e3*edgeF3(i,1)/lenF3(i,1));
            A(aaa_VF3(fb,2),aaa_VF3(fc,2)) = A(aaa_VF3(fb,2),aaa_VF3(fc,2)) + e3*edgeF3(i,2)/lenF3(i,2);
            A(aaa_VF3(fb,2),aaa_VF3(fa,2)) = A(aaa_VF3(fb,2),aaa_VF3(fa,2)) + e3*edgeF3(i,1)/lenF3(i,1);
        elseif j==3
            A(aaa_VF3(fc,2),aaa_VF3(fc,2)) = A(aaa_VF3(fc,2),aaa_VF3(fc,2)) - (e3*edgeF3(i,3)/lenF3(i,3) + e3*edgeF3(i,2)/lenF3(i,2));
            A(aaa_VF3(fc,2),aaa_VF3(fa,2)) = A(aaa_VF3(fc,2),aaa_VF3(fa,2)) + e3*edgeF3(i,3)/lenF3(i,3);
            A(aaa_VF3(fc,2),aaa_VF3(fb,2)) = A(aaa_VF3(fc,2),aaa_VF3(fb,2)) + e3*edgeF3(i,2)/lenF3(i,2);
        end
    end
end

for i=1:conrow
    for j=1:concol
        fcon=find(contact(i,j)==VFT);
        A(fcon,:) = 0;
        A(fcon,fcon) = 1;
    end
end

for ii=1:is1row
    minus1 = find(is1(ii)==aaa_VF1(:,1));
    plus1 = find(is1(ii)==aaa_VF2(:,1));
    A(aaa_VF2(plus1,2),:) = A(aaa_VF2(plus1,2),:) + A(aaa_VF1(minus1,2),:);
    A(aaa_VF1(minus1,2),:)=0;
    A(aaa_VF1(minus1,2),aaa_VF1(minus1,2))=-1;
    A(aaa_VF1(minus1,2),aaa_VF2(plus1,2))=1;
end

for ii=1:is2row
    minus1 = find(is2(ii)==aaa_VF3(:,1));
    plus1 = find(is2(ii)==aaa_VF2(:,1));
    A(aaa_VF2(plus1,2),:) = A(aaa_VF2(plus1,2),:) + A(aaa_VF3(minus1,2),:);
    A(aaa_VF3(minus1,2),:)=0;
    A(aaa_VF3(minus1,2),aaa_VF3(minus1,2))=-1;
    A(aaa_VF3(minus1,2),aaa_VF2(plus1,2))=1;
end

%%% edge contact
for i=1:conrow
    for j=1:concol
        fcon=find(contact(i,j)==VFT);
        if i<=Ntcon
            b(fcon) = tcon;
        elseif i>Ntcon &&  i<=Ntcon+Nbcon
            b(fcon) = bcon;
        end
    end
end

%%% initial potential
phi = A\b;

%%% initial electron
n = zeros(VF2row,1);
for i=1:VF2row
    n(i,1) = nint*exp(phi(aaa_VF2(i,2))/Vt);
end

%%% initial hole
p = zeros(VF2row,1);
for i=1:VF2row
    p(i) = nint*exp(-phi(aaa_VF2(i,2))/Vt);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% jac, res matrix
jac = zeros(VFTrow+VF2row+VF2row,VFTrow+VF2row+VF2row);
jacrow = size(jac,1);
res = zeros(jacrow,1);
for i=1:F1row
    for j=1:3
        fa = find(F1(i,1)==aaa_VF1(:,1));
        fb = find(F1(i,2)==aaa_VF1(:,1));
        fc = find(F1(i,3)==aaa_VF1(:,1));
        if j==1
            jac(aaa_VF1(fa,2),aaa_VF1(fa,2)) = jac(aaa_VF1(fa,2),aaa_VF1(fa,2)) - (e1*edgeF1(i,1)/lenF1(i,1) + e1*edgeF1(i,3)/lenF1(i,3));
            jac(aaa_VF1(fa,2),aaa_VF1(fb,2)) = jac(aaa_VF1(fa,2),aaa_VF1(fb,2)) + e1*edgeF1(i,1)/lenF1(i,1);
            jac(aaa_VF1(fa,2),aaa_VF1(fc,2)) = jac(aaa_VF1(fa,2),aaa_VF1(fc,2)) + e1*edgeF1(i,3)/lenF1(i,3);
            res(aaa_VF1(fa,2),1) = res(aaa_VF1(fa,2),1) - (e1*edgeF1(i,1)/lenF1(i,1) + e1*edgeF1(i,3)/lenF1(i,3))*phi(aaa_VF1(fa,2))...
                + e1*edgeF1(i,1)/lenF1(i,1)*phi(aaa_VF1(fb,2)) + e1*edgeF1(i,3)/lenF1(i,3)*phi(aaa_VF1(fc,2));

        elseif j==2
            jac(aaa_VF1(fb,2),aaa_VF1(fb,2)) = jac(aaa_VF1(fb,2),aaa_VF1(fb,2)) - (e1*edgeF1(i,2)/lenF1(i,2) + e1*edgeF1(i,1)/lenF1(i,1));
            jac(aaa_VF1(fb,2),aaa_VF1(fc,2)) = jac(aaa_VF1(fb,2),aaa_VF1(fc,2)) + e1*edgeF1(i,2)/lenF1(i,2);
            jac(aaa_VF1(fb,2),aaa_VF1(fa,2)) = jac(aaa_VF1(fb,2),aaa_VF1(fa,2)) + e1*edgeF1(i,1)/lenF1(i,1);
            res(aaa_VF1(fb,2),1) = res(aaa_VF1(fb,2),1) - (e1*edgeF1(i,2)/lenF1(i,2) + e1*edgeF1(i,1)/lenF1(i,1))*phi(aaa_VF1(fb,2))...
                + e1*edgeF1(i,2)/lenF1(i,2)*phi(aaa_VF1(fc,2)) + e1*edgeF1(i,1)/lenF1(i,1)*phi(aaa_VF1(fa,2));

        elseif j==3
            jac(aaa_VF1(fc,2),aaa_VF1(fc,2)) = jac(aaa_VF1(fc,2),aaa_VF1(fc,2)) - (e1*edgeF1(i,3)/lenF1(i,3) + e1*edgeF1(i,2)/lenF1(i,2));
            jac(aaa_VF1(fc,2),aaa_VF1(fa,2)) = jac(aaa_VF1(fc,2),aaa_VF1(fa,2)) + e1*edgeF1(i,3)/lenF1(i,3);
            jac(aaa_VF1(fc,2),aaa_VF1(fb,2)) = jac(aaa_VF1(fc,2),aaa_VF1(fb,2)) + e1*edgeF1(i,2)/lenF1(i,2);
            res(aaa_VF1(fc,2),1) = res(aaa_VF1(fc,2),1) - (e1*edgeF1(i,3)/lenF1(i,3) + e1*edgeF1(i,2)/lenF1(i,2))*phi(aaa_VF1(fc,2))...
                + e1*edgeF1(i,3)/lenF1(i,3)*phi(aaa_VF1(fa,2)) + e1*edgeF1(i,2)/lenF1(i,2)*phi(aaa_VF1(fb,2));
        end
    end
end

saveexl = zeros(VFTrow,1);

for i=1:F2row
    for j=1:3
        fa = find(F2(i,1)==aaa_VF2(:,1));
        fb = find(F2(i,2)==aaa_VF2(:,1));
        fc = find(F2(i,3)==aaa_VF2(:,1));
        if j==1
            jac(aaa_VF2(fa,2),aaa_VF2(fa,2)) = jac(aaa_VF2(fa,2),aaa_VF2(fa,2)) - (e2*edgeF2(i,1)/lenF2(i,1) + e2*edgeF2(i,3)/lenF2(i,3));
            jac(aaa_VF2(fa,2),aaa_VF2(fb,2)) = jac(aaa_VF2(fa,2),aaa_VF2(fb,2)) + e2*edgeF2(i,1)/lenF2(i,1);
            jac(aaa_VF2(fa,2),aaa_VF2(fc,2)) = jac(aaa_VF2(fa,2),aaa_VF2(fc,2)) + e2*edgeF2(i,3)/lenF2(i,3);

            res(aaa_VF2(fa,2),1) = res(aaa_VF2(fa,2),1) + (edgeF2(i,1)*lenF2(i,1) ...
                + edgeF2(i,3)*lenF2(i,3))/4*q/e0*Ndop...
                - (e2*edgeF2(i,1)/lenF2(i,1) + e2*edgeF2(i,3)/lenF2(i,3))*phi(aaa_VF2(fa,2))...
                + e2*edgeF2(i,1)/lenF2(i,1)*phi(aaa_VF2(fb,2)) + e2*edgeF2(i,3)/lenF2(i,3)*phi(aaa_VF2(fc,2));
            saveexl(aaa_VF2(fa,2),1) = saveexl(aaa_VF2(fa,2),1) + (edgeF2(i,1)*lenF2(i,1) + edgeF2(i,3)*lenF2(i,3))/4;

        elseif j==2
            jac(aaa_VF2(fb,2),aaa_VF2(fb,2)) = jac(aaa_VF2(fb,2),aaa_VF2(fb,2)) - (e2*edgeF2(i,2)/lenF2(i,2) + e2*edgeF2(i,1)/lenF2(i,1));
            jac(aaa_VF2(fb,2),aaa_VF2(fc,2)) = jac(aaa_VF2(fb,2),aaa_VF2(fc,2)) + e2*edgeF2(i,2)/lenF2(i,2);
            jac(aaa_VF2(fb,2),aaa_VF2(fa,2)) = jac(aaa_VF2(fb,2),aaa_VF2(fa,2)) + e2*edgeF2(i,1)/lenF2(i,1);

            res(aaa_VF2(fb,2),1) = res(aaa_VF2(fb,2),1) ...
                + (edgeF2(i,2)*lenF2(i,2) + edgeF2(i,1)*lenF2(i,1))/4*q/e0*Ndop...
                - (e2*edgeF2(i,2)/lenF2(i,2) + e2*edgeF2(i,1)/lenF2(i,1))*phi(aaa_VF2(fb,2))...
                + e2*edgeF2(i,2)/lenF2(i,2)*phi(aaa_VF2(fc,2)) + e2*edgeF2(i,1)/lenF2(i,1)*phi(aaa_VF2(fa,2));
            saveexl(aaa_VF2(fb,2),1) = saveexl(aaa_VF2(fb,2),1) + (edgeF2(i,2)*lenF2(i,2) + edgeF2(i,1)*lenF2(i,1))/4;

        elseif j==3
            jac(aaa_VF2(fc,2),aaa_VF2(fc,2)) = jac(aaa_VF2(fc,2),aaa_VF2(fc,2)) - (e2*edgeF2(i,3)/lenF2(i,3) + e2*edgeF2(i,2)/lenF2(i,2));
            jac(aaa_VF2(fc,2),aaa_VF2(fa,2)) = jac(aaa_VF2(fc,2),aaa_VF2(fa,2)) + e2*edgeF2(i,3)/lenF2(i,3);
            jac(aaa_VF2(fc,2),aaa_VF2(fb,2)) = jac(aaa_VF2(fc,2),aaa_VF2(fb,2)) + e2*edgeF2(i,2)/lenF2(i,2);

            res(aaa_VF2(fc,2),1) = res(aaa_VF2(fc,2),1) + (edgeF2(i,3)*lenF2(i,3) ...
                + edgeF2(i,2)*lenF2(i,2))/4*q/e0*Ndop...
                - (e2*edgeF2(i,3)/lenF2(i,3) + e2*edgeF2(i,2)/lenF2(i,2))*phi(aaa_VF2(fc,2))...
                + e2*edgeF2(i,3)/lenF2(i,3)*phi(aaa_VF2(fa,2)) + e2*edgeF2(i,2)/lenF2(i,2)*phi(aaa_VF2(fb,2));
            saveexl(aaa_VF2(fc,2),1) = saveexl(aaa_VF2(fc,2),1) + (edgeF2(i,3)*lenF2(i,3) + edgeF2(i,2)*lenF2(i,2))/4;
        end
    end
end

for i=1:F3row
    for j=1:3
        fa = find(F3(i,1)==aaa_VF3(:,1));
        fb = find(F3(i,2)==aaa_VF3(:,1));
        fc = find(F3(i,3)==aaa_VF3(:,1));
        if j==1
            jac(aaa_VF3(fa,2),aaa_VF3(fa,2)) = jac(aaa_VF3(fa,2),aaa_VF3(fa,2)) - (e3*edgeF3(i,1)/lenF3(i,1) + e3*edgeF3(i,3)/lenF3(i,3));
            jac(aaa_VF3(fa,2),aaa_VF3(fb,2)) = jac(aaa_VF3(fa,2),aaa_VF3(fb,2)) + e3*edgeF3(i,1)/lenF3(i,1);
            jac(aaa_VF3(fa,2),aaa_VF3(fc,2)) = jac(aaa_VF3(fa,2),aaa_VF3(fc,2)) + e3*edgeF3(i,3)/lenF3(i,3);
            res(aaa_VF3(fa,2),1) = res(aaa_VF3(fa,2),1) - (e3*edgeF3(i,1)/lenF3(i,1) + e3*edgeF3(i,3)/lenF3(i,3))*phi(aaa_VF3(fa,2))...
                + e3*edgeF3(i,1)/lenF3(i,1)*phi(aaa_VF3(fb,2)) + e3*edgeF3(i,3)/lenF3(i,3)*phi(aaa_VF3(fc,2));

        elseif j==2
            jac(aaa_VF3(fb,2),aaa_VF3(fb,2)) = jac(aaa_VF3(fb,2),aaa_VF3(fb,2)) - (e3*edgeF3(i,2)/lenF3(i,2) + e3*edgeF3(i,1)/lenF3(i,1));
            jac(aaa_VF3(fb,2),aaa_VF3(fc,2)) = jac(aaa_VF3(fb,2),aaa_VF3(fc,2)) + e3*edgeF3(i,2)/lenF3(i,2);
            jac(aaa_VF3(fb,2),aaa_VF3(fa,2)) = jac(aaa_VF3(fb,2),aaa_VF3(fa,2)) + e3*edgeF3(i,1)/lenF3(i,1);
            res(aaa_VF3(fb,2),1) = res(aaa_VF3(fb,2),1) - (e3*edgeF3(i,2)/lenF3(i,2) + e3*edgeF3(i,1)/lenF3(i,1))*phi(aaa_VF3(fb,2))...
                + e3*edgeF3(i,2)/lenF3(i,2)*phi(aaa_VF3(fc,2)) + e3*edgeF3(i,1)/lenF3(i,1)*phi(aaa_VF3(fa,2));

        elseif j==3
            jac(aaa_VF3(fc,2),aaa_VF3(fc,2)) = jac(aaa_VF3(fc,2),aaa_VF3(fc,2)) - (e3*edgeF3(i,3)/lenF3(i,3) + e3*edgeF3(i,2)/lenF3(i,2));
            jac(aaa_VF3(fc,2),aaa_VF3(fa,2)) = jac(aaa_VF3(fc,2),aaa_VF3(fa,2)) + e3*edgeF3(i,3)/lenF3(i,3);
            jac(aaa_VF3(fc,2),aaa_VF3(fb,2)) = jac(aaa_VF3(fc,2),aaa_VF3(fb,2)) + e3*edgeF3(i,2)/lenF3(i,2);
            res(aaa_VF3(fc,2),1) = res(aaa_VF3(fc,2),1) - (e3*edgeF3(i,3)/lenF3(i,3) + e3*edgeF3(i,2)/lenF3(i,2))*phi(aaa_VF3(fc,2))...
                + e3*edgeF3(i,3)/lenF3(i,3)*phi(aaa_VF3(fa,2)) + e3*edgeF3(i,2)/lenF3(i,2)*phi(aaa_VF3(fb,2));
        end
    end
end

for i=1:VF2row
    jac(bb1_VF2(i,2),bb1_VF2(i,2)) = 1;
    jac(bb1_VF2(i,2),aaa_VF2(i,2)) = - nint/Vt*exp(phi(aaa_VF2(i,2))/Vt);
    res(bb1_VF2(i,2),1) = n(i) - nint*exp(phi(aaa_VF2(i,2),1)/Vt);

    jac(aaa_VF2(i,2),bb1_VF2(i,2)) = - saveexl(aaa_VF2(i,2),1)*q/e0;
    res(aaa_VF2(i,2),1) = res(aaa_VF2(i,2),1) - saveexl(aaa_VF2(i,2))*q/e0*n(i);
end

for i=1:VF2row
    jac(bb2_VF2(i,2),bb2_VF2(i,2)) = 1;
    jac(bb2_VF2(i,2),aaa_VF2(i,2)) = nint/Vt*exp(-phi(aaa_VF2(i,2))/Vt);
    res(bb2_VF2(i,2),1) = p(i) - nint*exp(-phi(aaa_VF2(i,2),1)/Vt);

    jac(aaa_VF2(i,2),bb2_VF2(i,2)) = saveexl(aaa_VF2(i,2),1)*q/e0;
    res(aaa_VF2(i,2),1) = res(aaa_VF2(i,2),1) + saveexl(aaa_VF2(i,2))*q/e0*p(i);
end

for i=1:conrow
    for j=1:concol
        fcon = find(contact(i,j)==VFT);
        jac(fcon,:)=0;
        jac(fcon,fcon)=1;
        if i<=Ntcon
            res(fcon) = phi(fcon)-tcon;
        elseif i>Ntcon && i<=Ntcon+Nbcon
            res(fcon) = phi(fcon)-bcon;
        end
    end
end

for ii=1:is1row
    minus1 = find(is1(ii)==aaa_VF1(:,1));
    plus1 = find(is1(ii)==aaa_VF2(:,1));
    jac(aaa_VF2(plus1,2),:) = jac(aaa_VF2(plus1,2),:) + jac(aaa_VF1(minus1,2),:);
    jac(aaa_VF1(minus1,2),:)=0;
    jac(aaa_VF1(minus1,2),aaa_VF1(minus1,2))=-1;
    jac(aaa_VF1(minus1,2),aaa_VF2(plus1,2))=1;
end

for ii=1:is2row
    minus1 = find(is2(ii)==aaa_VF3(:,1));
    plus1 = find(is2(ii)==aaa_VF2(:,1));
    jac(aaa_VF2(plus1,2),:) = jac(aaa_VF2(plus1,2),:) + jac(aaa_VF3(minus1,2),:);
    jac(aaa_VF3(minus1,2),:)=0;
    jac(aaa_VF3(minus1,2),aaa_VF3(minus1,2))=-1;
    jac(aaa_VF3(minus1,2),aaa_VF2(plus1,2))=1;
end

update = jac \(-res);

for i=1:jacrow
    if i <= VFTrow
        phi(i) = phi(i) + update(i);
    elseif i > VFTrow && i <= (VFTrow + VF2row)
        n(i-VFTrow) = n(i-VFTrow) + update(i);
    elseif i > (VFTrow + VF2row) && i <= jacrow
        p(i-(VFTrow + VF2row)) = p(i-(VFTrow + VF2row)) + update(i);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pot = zeros(VFrow,1); %%% potential for total regions
% for i=1:size(pot,1)
%     pot(i)=phi_sorted(i);
% end
%
% elec = zeros(VF2row,1); %%% electron for silicon region
% for i=1:size(elec,1)
%     elec(i)=phi(VF2(i));
% end
%
% hole = zeros(VF2row,1); %%% hole for silicon region
% for i=1:size(hole,1)
%     hole(i)=phi(VF2(i));
% end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%  solution X matrix
% X = zeros(VFrow+VF2row+VF2row,1);
% for i=1:size(X,1)
%     if i <= VFrow
%         X(i) = pot(i);
%     elseif i > VFrow && i<=VFrow+VF2row
%         X(i) = elec(i-VFrow);
%     elseif i > VFrow+VF2row && i<=VFrow+VF2row+VF2row
%         X(i) = hole(i-(VFrow+VF2row));
%     end
% end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%% scan part
% while 1
%     var = input('variable?\n ▶ pot\n ▶ elec\n ▶ hole\n : ', 's');
%     reg = input('region?\n (1)oxide\n (2)silicon\n : ');
%     if (strcmp(var, 'elec') && reg == 1) || (strcmp(var, 'hole') && reg == 2)
%         fprintf('Error!\n');
%         continue
%     end
%     ver = input('vertex number? : ');
%     %%%%% print part
%     if strcmp(var, 'pot') && size(find(ver==VF),1)==1
%         xindex = find(ver==VF);
%         fprintf('X[%d] = %f\n',xindex, X(xindex));
%     elseif strcmp(var, 'elec') && size(find(ver==VF2),1)==1
%         xindex = VFrow+find(ver==VF2);
%         fprintf('X[%d] = %f\n',xindex, X(xindex));
%     elseif strcmp(var, 'hole') && size(find(ver==VF2),1)==1
%         xindex = VFrow+VF2row+find(ver==VF2);
%         fprintf('X[%d] = %f\n',xindex, X(xindex));
%     else
%         fprintf('Error !\n');
%     end
%     break
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% visualize part
% figure
% patch('Faces',F1,'Vertices',V, 'FaceVertexCData',phi, 'EdgeColor','blue','FaceColor','interp','LineWidth',1, 'Marker','o');
% patch('Faces',F2,'Vertices',V, 'FaceVertexCData',phi, 'EdgeColor','green','FaceColor','interp','LineWidth',1, 'Marker','o');
% patch('Faces',F3,'Vertices',V, 'FaceVertexCData',phi, 'EdgeColor','blue','FaceColor','interp','LineWidth',1, 'Marker','o');
% title('Φ visualizing structure')
% patch('Faces',interedge1,'Vertices',V, 'EdgeColor','magenta','FaceColor','none','LineWidth', 3, 'Marker','o');
% patch('Faces',interedge2,'Vertices',V, 'EdgeColor','magenta','FaceColor','none','LineWidth', 3, 'Marker','o');
% patch('Faces',interedge3,'Vertices',V, 'EdgeColor','magenta','FaceColor','none','LineWidth', 3, 'Marker','o');
% colorbar
