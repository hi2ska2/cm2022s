function result = Contact(x1, y1, x2, y2, potential)
% 사용 범위 안에 들면 Contact으로 설정되게 했다.

  
    nm=1e-9;
    
    x1=x1*nm; y1=y1*nm; x2=x2*nm; y2=y2*nm;
    
    norm=1e-20;

    vertex = importdata("Vertex.txt");
    n=1;
    contact=[];
    for ii=1:length(vertex)
        x_tmp=vertex(ii,1); y_tmp=vertex(ii,2);
        if x_tmp-x1>=-norm && x_tmp-x2<=norm && y_tmp-y1>=-norm && y_tmp-y2<=norm
            contact(1,n)=ii;
            n=n+1;
        end
    end
    result=[contact potential];