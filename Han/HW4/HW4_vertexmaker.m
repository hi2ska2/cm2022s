%%% HW4_Vertex maker
%%% written by Seong-Min, Han
numc = 3;
numv = 2*numc*numc+6*numc+1;
vertex = zeros(numv, 2);

R = zeros(numc,1); % 각 원의 반지름
R(1,1) = 1;
R(2,1) = R(1,1) * 2;
R(3,1) = R(1,1) * 2.5;

for i=1:numc
    for j=1:(i*4+4)
        vertex(2*(i*i+i-1)-1+j,1) = R(i,1) * cospi((j-1)/(2*(i+1)));
        vertex(2*(i*i+i-1)-1+j,2) = R(i,1) * sinpi((j-1)/(2*(i+1)));
    end
end

writematrix(vertex,'vertex.txt');
type vertex.txt
