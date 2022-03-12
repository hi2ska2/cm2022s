%%% HW3_Vertex maker: use num circle
%%% written by Seong-Min, Han
N = 4; % 가중치
numc = 4; % circle의 개수
numv = 2*numc*numc+2*numc+1; % Vertex 개수
Vertex = zeros(numv, 2);

R = zeros(numc,1); % 각 원의 반지름
for i=2:numc
    R(1,1) = 3;
    R(i,1) = R(1,1) * i;
end

for i=1:numc
    for j=1:1:i*N
        Vertex(2*(i*i-i+1)-1+j,1) = R(i,1) * cos((j-1)*pi/(2*i));
        Vertex(2*(i*i-i+1)-1+j,2) = R(i,1) * sin((j-1)*pi/(2*i));
    end
end

writematrix(Vertex,'Vertex.txt');
type Vertex.txt
