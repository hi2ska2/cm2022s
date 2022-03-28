clear; close all; clc;

%% Vertex/Element import %%
% vertex 값 불러오기
point = importdata("Vertex.txt");
point_tmp=zeros(length(point),1);
point = [point point_tmp];

% element 값 불러오기
F_1 = importdata("element_1.txt");
F_2 = importdata("element_2.txt");
F_3 = importdata("element_3.txt");
Element = [F_1; F_2; F_3];

%% Find interface %%
% 각 Region의 edge 정렬
n=1;
for ii=1:length(F_1) % Si region
   for j=1:3
      edge_1(n,:)= [F_1(ii,2) F_1(ii,1)];
      edge_1(n+1,:)= [F_1(ii,3) F_1(ii,2)];
      edge_1(n+2,:)= [F_1(ii,1) F_1(ii,3)];
      n=n+3;
   end
end
n=1;
for ii=1:length(F_2) % ox region
   for j=1:3
      edge_2(n,:)= [F_2(ii,1) F_2(ii,2)];
      edge_2(n+1,:)= [F_2(ii,2) F_2(ii,3)];
      edge_2(n+2,:)= [F_2(ii,3) F_2(ii,1)];
      n=n+3;
   end
end
n=1;
for ii=1:length(F_3) % Si region
   for j=1:3
      edge_3(n,:)= [F_3(ii,2) F_3(ii,1)];
      edge_3(n+1,:)= [F_3(ii,3) F_3(ii,2)];
      edge_3(n+2,:)= [F_3(ii,1) F_3(ii,3)];
      n=n+3;
   end
end

% 중복 edge 정보 제거
edge_1=sort(edge_1, 2); % edge_si의 각 행을 오름차순으로 정렬(중복값 제거 위해)
edge_1=unique(edge_1,'rows');
edge_2=sort(edge_2, 2); % edge_ox의 각 행을 오름차순으로 정렬(중복값 제거 위해)
edge_2=unique(edge_2,'rows');
edge_3=sort(edge_3, 2); % edge_si의 각 행을 오름차순으로 정렬(중복값 제거 위해)
edge_3=unique(edge_3,'rows');


% 각 region에서 동일한 edge가 있는지 확인, Interface
edge_interface1 = intersect(edge_1,edge_2,'rows');
edge_interface2 = intersect(edge_1,edge_3,'rows');
edge_interface3 = intersect(edge_2,edge_3,'rows');
edge_interface=[edge_interface1; edge_interface2; edge_interface3];

%% Find the solution index %%    <--- HW7 part.
%각 region을 구성하는 vertex의 index를 확인하는 과정
Vertex_region1 = unique(sort(reshape(F_1,[],1))); 
Vertex_region2 = unique(sort(reshape(F_2,[],1)));
Vertex_region3 = unique(sort(reshape(F_3,[],1)));

%각 Reigon 내의 vertex 갯수 저장
Number_of_vertex=[length(Vertex_region1); length(Vertex_region2); length(Vertex_region3)];

% 임의의 Variable 지정
aaa_tmp = rand(sum(Number_of_vertex),1);
bbb_tmp = rand(Number_of_vertex(2,1),1);
ccc_tmp = rand(Number_of_vertex(1,1)+Number_of_vertex(3,1),1);

% Variable, region 정보 저장
Variable_region=[1 1 1; % Variable 'aaa' is in Region1,2,3(all) 
                 0 1 0; % Variable 'bbb' is in Region2(si) 
                 1 0 1]; % Variable 'ccc' is in Region1,3(ox) 
% row: Variable. ex)aaa,bbb,ccc), columm: Region

% Solution vector x 생성
aaa=NaN(sum(Number_of_vertex),1);
bbb=NaN(sum(Number_of_vertex),1);
ccc=NaN(sum(Number_of_vertex),1);
%aaa
aaa=aaa_tmp; 
%bbb
for ii=Number_of_vertex(1,1)+1:Number_of_vertex(1,1)+Number_of_vertex(2,1)
    j=ii-Number_of_vertex(1,1);
    bbb(ii,1)=bbb_tmp(j,1);
end
%ccc
for ii=1:Number_of_vertex(1,1)
    ccc(ii,1)=ccc_tmp(ii,1);
end
for ii=Number_of_vertex(1,1)+Number_of_vertex(2,1)+1:sum(Number_of_vertex)
    j=ii-Number_of_vertex(2,1);
    ccc(ii,1)=ccc_tmp(j,1);
end
x_tmp=[aaa bbb ccc];

%x를 행백터로 변환
n=1;
for ii=1:sum(Number_of_vertex)
    for j=1:size(Variable_region,1)
        x(n,1)=x_tmp(ii,j);
        n=n+1;
    end
end
%Nan값 제거
x=x(~isnan(x));

Repeat = "y";
while Repeat=="y"
    %Input is [variable region vetex]
    prompt = 'Please enter the input. ex) [Variable Region Vetex]\n: ';
    information=input(prompt);

    %입력한 값이 올바른 값인지 검증
    while information(1,1)>size(Variable_region,1) || information(1,2)>size(Variable_region,2)
        fprintf('\n%% Please enter the corret information %%\n\n')
        %Input is [variable region vetex]
        prompt = '\n Please enter the input. ex) [Variable Region Vetex]\n: ';
        information=input(prompt);
    end
    while Variable_region(information(1,1), information(1,2))==0 || information(1,3) >= Number_of_vertex(information(1,2),1)
        fprintf('\n%% Please enter the corret information %%\n\n')
        %Input is [variable region vetex]
        prompt = '\n Please enter the input. ex) [Variable Region Vetex]\n: ';
        information=input(prompt);

    end

    % Find index
    if information(1,2)==1
        if information(1,1)==1
            index= sum(Variable_region(:,1))*information(1,3)-1;
        elseif information(1,1)==3
            index= sum(Variable_region(:,1))*information(1,3);
        end
    elseif information(1,2)==2
        if information(1,1)==1
            index= sum(Variable_region(:,2))*information(1,3)-1;
        elseif information(1,1)==2
            index= sum(Variable_region(:,2))*information(1,3);
        end
        index=index+Number_of_vertex(1,1)*sum(Variable_region(:,1));
    elseif information(1,2)==3
        if information(1,1)==1
            index= sum(Variable_region(:,3))*information(1,3)-1;
        elseif information(1,1)==3
            index= sum(Variable_region(:,3))*information(1,3);
        end
        index=index+Number_of_vertex(1,1)*sum(Variable_region(:,1))+Number_of_vertex(2,1)*sum(Variable_region(:,2));
    end

    % 찾은 index 표시
    fprintf('\n Index : %d\n x(%d) : %f \n\n',index, index, x(index))

    % 다시 index를 찾을것인지?
    prompt = 'Do you want more? Y/N [Y]: ';
    Repeat = input(prompt,'s');
    fprintf('----------------------------------------------------\n')
end
%% Visualize %%
% % figure(1);
% % patch('Faces',Element(:,1:3),'Vertices',point(:,[1 2]), 'EdgeColor','black','FaceColor','none','LineWidth',0.5);
% % xlabel('X');
% % ylabel('Y');
% % title('Structure')
% 
% figure(2)
% patch('Faces',F_1,'Vertices',point, 'EdgeColor','Blue','FaceColor','none','LineWidth',2)
% hold on
% patch('Faces',F_2,'Vertices',point, 'EdgeColor','Red','FaceColor','none','LineWidth',2)
% hold on
% patch('Faces',F_3,'Vertices',point, 'EdgeColor','Yellow','FaceColor','none','LineWidth',2)
% hold on
% patch('Faces',edge_interface,'Vertices',point, 'EdgeColor','Green','FaceColor','none','LineWidth',5)
% hold off
% title('Structure')
% xlabel('X');
% ylabel('Y');
% axis([0 5 0 5]);
% % caxis([0 2]);