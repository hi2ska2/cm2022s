% due to : 2022.03.10. (Thur) / HW2

% 2D Triangle Laplace equation

% Vertex information

P = zeros(3,2);

for ii = 1:3
    if ii == 1
        P(ii,1) = 0;
        P(ii,2) = 0;
    elseif ii == 2
        P(ii,1) = 6;
        P(ii,2) = 0;
    elseif ii == 3
        P(ii,1) = 4;
        P(ii,2) = 4;
    end
end

% Geometry information

L=zeros(3,1);

for ii = 1:3
    if ii == 1 % between 1 and 2
        L(ii,1) = sqrt(abs((P(ii+1,1)-P(ii,1))^2+(P(ii+1,2)-P(ii,2))^2));
    elseif ii == 2 % between 2 and 3
        L(ii,1) = sqrt(abs((P(ii+1,1)-P(ii,1))^2+(P(ii+1,2)-P(ii,2))^2));
    elseif ii ==3 % between 3 and 1
        L(ii,1) = sqrt(abs((P(ii,1)-P(ii-2,1))^2+(P(ii,2)-P(ii-2,2))^2));
    end
end

% Vector
a12 = [P(2,1)-P(1,1), P(2,2)-P(1,2)];
b13 = [P(3,1)-P(1,1), P(3,2)-P(1,2)];
a21=-a12;
b23 = [P(3,1)-P(2,1), P(3,2)-P(2,2)];
a31=-b13;
b32=-b23;

Area = (1/2)*abs(a12(1,1)*b13(1,2)-a12(1,2)*b13(1,1));
R = (L(1,1)*L(2,1)*L(3,1))/(4*Area);

A=zeros(3,1);

A(1,1) = sqrt(R^2-(L(1,1)/2)^2);
A(2,1) = sqrt(R^2-(L(2,1)/2)^2);
A(3,1) = sqrt(R^2-(L(3,1)/2)^2);

% Potential 

Jaco = zeros(3,3);
res = zeros(3,1);
res(3,1) = 1 ;

for ii = 1:3
    if ii == 1
%         Jaco(ii,ii) = -(A(ii,1)/L(ii,1)+A(ii+2,1)/L(ii+2,1));
%         Jaco(ii,ii+1) = (A(ii,1)/L(ii,1));
%         Jaco(ii,ii+2) = (A(ii+2,1)/L(ii+2,1));
        Jaco(ii,ii) = 1;
        
    elseif ii == 2
        Jaco(ii,ii-1) = (A(ii-1,1)/L(ii-1,1));
        Jaco(ii,ii) = -(A(ii-1,1)/L(ii-1,1)+A(ii,1)/L(ii,1));
        Jaco(ii,ii+1) = (A(ii,1)/L(ii,1));
        
    elseif ii == 3
%         Jaco(ii,ii-2) = (A(ii-2,1)/L(ii-2,1));
%         Jaco(ii,ii-1) = (A(ii-1,1)/L(ii-1,1));
%         Jaco(ii,ii) = -(A(ii,1)/L(ii,1)+A(ii-1,1)/L(ii-1,1));
        Jaco(ii,ii) = 1;

    end
end

phi=Jaco\res;