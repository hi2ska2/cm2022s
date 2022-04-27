function vertex = vertex_doublegate(x1, y1, x2, y2, dx, dy)
% double-gate의 형태를 가진 모델의 vertex를 찍기 위해 작성했다.
%
% 직사각형 형태의 범위를 일정 간격으로 나누어 vertex를 생성하는 함수이다.
% 입력 인수는 총 6개이며 (x1, y1, x2, y2, dx, dy) 순서로 입력하면 vertex가 작성이 된다.
% vertex의 배치(index)는 좌측 하단에서 우측 상단으로 올라간다.
% 먼저 y를 오름차순으로 정렬 한 후 같은 y내에서 x를 오름차순으로 정렬한다.
% 추가적으로 Vertex.txt파일을 생성한다.
%
% 사용예시는 다음과 같다.
%
% vertex = vertex_doublegate(x1, y1, x2, y2, dx, dy)

    format long

     nm=1e-9;
     
     x1=x1*nm; y1=y1*nm; x2=x2*nm; y2=y2*nm; dx=dx*nm; dy=dy*nm;

    x=(0:dx:x2-x1)'+x1;
    y=(0:dy:y2-y1)'+y1;
    vertex=zeros(length(x)*length(y),2);
    for ii=1:length(y)
        for j=1:length(x)
            for m=1:2
                v=[x(j) y(ii)];
                k=length(x)*ii-(length(x)-1)+(j-1);
                vertex(k,m)=v(:,m);
            end
        end
    end
    
    vertex=[vertex zeros(length(vertex),1)];
    writematrix(vertex,'Vertex.txt'); % vertex 파일 출력
end