% written by Seong-Min, Han

P = zeros(3,2);
P(1,:) = [0 0];
P(2,:) = [1 0];
P(3,:) = [1 1]; % case1) 
% P(3,:) = [0.5 0.5*sqrt(3)]; % case2) An equilateral triangle
% P(3,:) = [3 4]; % case3)
% P(3,:) = [2/3 3/2]; case4) % acute triangle

v21 = P(2,:)-P(1,:);
v31 = P(3,:)-P(1,:);
norm21 = norm(P(2,:)-P(1,:));
norm31 = norm(P(3,:)-P(1,:));
area = sqrt(norm21*norm21 * norm31*norm31- dot(v21,v31)*dot(v21,v31))/2;

L = zeros(3,1);
L(1,1) = norm(P(2,:)-P(1,:));
L(2,1) = norm(P(3,:)-P(2,:));
L(3,1) = norm(P(3,:)-P(1,:));

R = (L(1,1)*L(2,1)*L(3,1)) / (4*area);

edge = zeros(3,1);
for i=1:3
    edge(i,1) = sqrt(R*R - (L(i)/2)*(L(i)/2));
end

jac = zeros(3,3);
res = zeros(3,1);
res(3,1) = 1;
for i=1:3
    if i==1 || i==3
        jac(i,i) = 1;
    else
        jac(i,i-1) = edge(i-1)/L(i-1)/2;
        jac(i,i) = -edge(i-1)/L(i-1)/2 - edge(i)/L(i)/2;
        jac(i,i+1) = edge(i)/L(i)/2;
    end
end
phi = jac\ res;

disp(phi(2))
x=1:1:3;
plot(x, phi, '-.r*')
title('2D Laplacian equation solver')
xlabel('vertex')
ylabel('potential')
