function Element = element(x1, y1, x2, y2)
% 직사각형 형태의 Region을 가질때, 이 Region 내 element를 작성했다.
%
% 직사각형 형태의 범위를 일정 간격으로 나누어 element를 생성하는 함수이다.
% 입력 인수는 총 4개이며 (x1, y1, x2, y2) 순서로 입력하면 element가 작성이 된다.
%
% 사용예시는 다음과 같다.
%
% Element = element(x1, y1, x2, y2)

  
    nm=1e-9;
    
    x1=x1*nm; y1=y1*nm; x2=x2*nm; y2=y2*nm;

    norm=1e-20;

    vertex = importdata("Vertex.txt");
    n=1;
    for ii=1:length(vertex)
        x_tmp=vertex(ii,1); y_tmp=vertex(ii,2);
        if x_tmp-x1>=-norm && x_tmp-x2<=norm && y_tmp-y1>=-norm && y_tmp-y2<=norm
            vertex_tmp(n,:)=vertex(ii,:);
            table(n,:)=[n ii]; %[re-index original-index]
            n=n+1;
        end
    end
    
    % vetex_tmp를 기반으로 한 element 작성
    x_unique=unique(sort(vertex_tmp(:,1)));
    y_unique=unique(sort(vertex_tmp(:,2)));
    n=1;
    for ii=1:length(y_unique)-1
        for m=1:length(x_unique)-1
            Element_tmp(n,:)=[length(x_unique)*(ii-1)+m length(x_unique)*(ii)+m+1 length(x_unique)*(ii)+m];
            Element_tmp(n+1,:)=[length(x_unique)*(ii-1)+m length(x_unique)*(ii-1)+m+1 length(x_unique)*(ii)+m+1];
            n=n+2;
        end
    end
    
    % Tabel을 활용해서 기존 vertex의 index에 맞게 element를 다시 구성
    Element=zeros(size(Element_tmp,1),size(Element_tmp,2));
    for ii=1:size(Element_tmp,1)
        for j=1:3
            Element(ii,j)=table(Element_tmp(ii,j),2);
        end
    end
end