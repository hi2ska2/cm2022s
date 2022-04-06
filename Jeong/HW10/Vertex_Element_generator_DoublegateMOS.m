clear; close all; clc;

% Make a Dublegate Mosfet (Vertex/Element)

% Variable
a=[5 3]; % size
N=0.5; % 점 사이의 간격
RegionX_ox1=[0 1];  % 각 Region 범위 지정
RegionY_ox1=[1 2];
RegionX_si=[1 4];
RegionY_si=[1 2];
RegionX_ox2=[4 5];
RegionY_ox2=[1 2];

%%%%%%%%%%%%%%%%%%%%%%% Vertex %%%%%%%%%%%%%%%%%%%%%%%
x=(0:N:a(1,1))';
y=(0:N:a(1,2))';
point_tmp=zeros(length(x)*length(y),2);
for ii=1:length(y)
    for j=1:length(x)
        for m=1:2
            v=[x(j) y(ii)];
            k=length(x)*ii-(length(x)-1)+(j-1);
            point_tmp(k,m)=v(:,m);
        end
    end
end

% 각 Region별로 vertex 저장
n=1;
for ii=1:length(point_tmp)     % ox1
    x_tmp=point_tmp(ii,1); y_tmp=point_tmp(ii,2);
    if x_tmp>=RegionX_ox1(1,1) && x_tmp<=RegionX_ox1(1,2) && y_tmp>=RegionY_ox1(1,1) && y_tmp<=RegionY_ox1(1,2)
        vertex_ox1(n,:)=[point_tmp(ii,:) 0];
        n=n+1;
    end
end
n=1;
for ii=1:length(point_tmp)     % ox2
    x_tmp=point_tmp(ii,1); y_tmp=point_tmp(ii,2);
    if x_tmp>=RegionX_ox2(1,1) && x_tmp<=RegionX_ox2(1,2) && y_tmp>=RegionY_ox2(1,1) && y_tmp<=RegionY_ox2(1,2)
        vertex_ox2(n,:)=[point_tmp(ii,:) 0];
        n=n+1;
    end
end
n=1;
for ii=1:length(point_tmp)     % si
    x_tmp=point_tmp(ii,1); y_tmp=point_tmp(ii,2);
    if x_tmp>=RegionX_si(1,1) && x_tmp<=RegionX_si(1,2) && y_tmp>=RegionY_si(1,1) && y_tmp<=RegionY_si(1,2)
        vertex_si(n,:)=[point_tmp(ii,:) 0];
        n=n+1;
    end
end


writematrix(vertex_ox1,'Vertex_ox1.txt'); % vertex 파일 출력
writematrix(vertex_ox2,'Vertex_ox2.txt'); % vertex 파일 출력
writematrix(vertex_si,'Vertex_si.txt'); % vertex 파일 출력
%%%%%%%%%%%%%%%%%%%%%%% Vertex end %%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%% Element %%%%%%%%%%%%%%%%%%%%%%%

% Region1(ox1)
x1=unique(sort(vertex_ox1(:,1)));
y1=unique(sort(vertex_ox1(:,2)));
n=1;
for ii=1:length(y1)-1
    for m=1:length(x1)-1
        Element_ox1(n,:)=[length(x1)*(ii-1)+m length(x1)*(ii)+m+1 length(x1)*(ii)+m];
        Element_ox1(n+1,:)=[length(x1)*(ii-1)+m length(x1)*(ii-1)+m+1 length(x1)*(ii)+m+1];
        n=n+2;
    end
end

% Region1(si)
x2=unique(sort(vertex_si(:,1)));
y2=unique(sort(vertex_si(:,2)));
n=1;
for ii=1:length(y2)-1
    for m=1:length(x2)-1
        Element_si(n,:)=[length(x2)*(ii-1)+m length(x2)*(ii)+m+1 length(x2)*(ii)+m];
        Element_si(n+1,:)=[length(x2)*(ii-1)+m length(x2)*(ii-1)+m+1 length(x2)*(ii)+m+1];
        n=n+2;
    end
end

% Region1(ox2)
x3=unique(sort(vertex_ox2(:,1)));
y3=unique(sort(vertex_ox2(:,2)));
n=1;
for ii=1:length(y3)-1
    for m=1:length(x3)-1
        Element_ox2(n,:)=[length(x3)*(ii-1)+m length(x3)*(ii)+m+1 length(x3)*(ii)+m];
        Element_ox2(n+1,:)=[length(x3)*(ii-1)+m length(x3)*(ii-1)+m+1 length(x3)*(ii)+m+1];
        n=n+2;
    end
end

writematrix(Element_ox1,'Element_ox1.txt'); % element 파일 출력
writematrix(Element_si,'Element_si.txt'); % element 파일 출력
writematrix(Element_ox2,'Element_ox2.txt'); % element 파일 출력
%%%%%%%%%%%%%%%%%%%%%%% Element end %%%%%%%%%%%%%%%%%%%%%%%

% Visualization
figure(1); plot(point_tmp(:,1),point_tmp(:,2), '*');
figure(2); plot(vertex_ox1(:,1), vertex_ox1(:,2), '*');
axis([0 a(1,1) 0 a(1,2)])
figure(3); plot(vertex_ox2(:,1), vertex_ox2(:,2), '*');
axis([0 a(1,1) 0 a(1,2)])
figure(4); plot(vertex_si(:,1), vertex_si(:,2), '*');
axis([0 a(1,1) 0 a(1,2)])
figure(5);
patch('Faces',Element_ox1,'Vertices',vertex_ox1, 'EdgeColor','black','FaceColor','none','LineWidth',2);
xlabel('X');
ylabel('Y');
title('Structure')
axis([0 a(1,1) 0 a(1,2)])
figure(6);
patch('Faces',Element_si,'Vertices',vertex_si, 'EdgeColor','black','FaceColor','none','LineWidth',2);
xlabel('X');
ylabel('Y');
title('Structure')
axis([0 a(1,1) 0 a(1,2)])
figure(7);
patch('Faces',Element_ox2,'Vertices',vertex_ox2, 'EdgeColor','black','FaceColor','none','LineWidth',2);
xlabel('X');
ylabel('Y');
title('Structure')
axis([0 a(1,1) 0 a(1,2)])