clear; clc; close all;
load('./data/homework25_dd.mat')
mkdir('./data/homework25_Vd')
%% Ramping Drain Voltage
Iteration_Vdd=100000;

% Jacobian matrix / res vector
solution_vector_Vdd(:,1)=solution_vector_DD(:,1);
solution_vector_Vdd_saved(:,1)=solution_vector_Vdd(:,1);
count_Vdd_current=1;
for bias_Vd=1:11 %%% Vd는 1V까지 올림.
    Vdd=0.1*bias_Vd-0.1;
    for Newton_Vdd=1:Iteration_Vdd
        fprintf("Ramping Drain Voltage, Vd=%.2f, Newton_Vd=%d\n" , Vdd, Newton_Vdd)
        clearvars Jaco_Vdd; clearvars res_Vdd;
        %     Jaco_Vdd=sparse(zeros(size(solution_vector_Vdd,1),size(solution_vector_Vdd,1))); res_Vdd=zeros(size(solution_vector_Vdd,1),1);
        Jaco_Vdd=zeros(size(solution_vector_Vdd,1),size(solution_vector_Vdd,1)); res_Vdd=zeros(size(solution_vector_Vdd,1),1); %% 희소행렬 아님. Jaco matrix 확인용.

        % Jaco/res %
        for ii=1:length(Element)
            % Control_Volume 계산
            Control_Volume=zeros(3,1);
            Control_Volume(1,1)=(0.25)*L(ii,3)*edge(ii,3)+(0.25)*L(ii,1)*edge(ii,1);
            Control_Volume(2,1)=(0.25)*L(ii,1)*edge(ii,1)+(0.25)*L(ii,2)*edge(ii,2);
            Control_Volume(3,1)=(0.25)*L(ii,2)*edge(ii,2)+(0.25)*L(ii,3)*edge(ii,3);

            % 현재 region의 eps_r 불러오기
            eps_now=eps(Table_element_region(ii,1),1);

            if Table_element_region(ii,1)==2 || Table_element_region(ii,1)==3 || Table_element_region(ii,1)==4 %%% Region 변경되면 꼭 수정되어야함!!
                % index : 현재 element에 속해있는 vertex가 resdue에 몇번째 있는가 저장
                index=zeros(9,1);
                n=1;
                for j=1:3
                    index(n,1)=find(Table_Jaco(:,1)==Table_element_region(ii,1) & Table_Jaco(:,2)==Element(ii,j) & Table_Jaco(:,3)==1);
                    index(n+1,1)=index(n,1)+1;
                    index(n+2,1)=index(n,1)+2;
                    n=n+3;
                end

                % 식 간단히 표기하기 위해 미리 계산.
                x12_Vdd=solution_vector_Vdd(index(1,1),1)-solution_vector_Vdd(index(4,1),1); x21_Vdd=-x12_Vdd;
                x23_Vdd=solution_vector_Vdd(index(4,1),1)-solution_vector_Vdd(index(7,1),1); x32_Vdd=-x23_Vdd;
                x13_Vdd=solution_vector_Vdd(index(1,1),1)-solution_vector_Vdd(index(7,1),1); x31_Vdd=-x13_Vdd;
                n1_Vdd=solution_vector_Vdd(index(1,1)+1,1); n2_Vdd=solution_vector_Vdd(index(4,1)+1,1); n3_Vdd=solution_vector_Vdd(index(7,1)+1,1);
                p1_Vdd=solution_vector_Vdd(index(1,1)+2,1); p2_Vdd=solution_vector_Vdd(index(4,1)+2,1); p3_Vdd=solution_vector_Vdd(index(7,1)+2,1);

                % Jaco %
                Jaco_tmp_Vdd=zeros(9,9);
                Jaco_tmp_Vdd(1,1)=eps_now*(-edge(ii,1)/L(ii,1)-edge(ii,3)/L(ii,3));
                Jaco_tmp_Vdd(1,2)=-coeff*Control_Volume(1,1);
                Jaco_tmp_Vdd(1,3)=coeff*Control_Volume(1,1);
                Jaco_tmp_Vdd(1,4)=eps_now*edge(ii,1)/L(ii,1);
                Jaco_tmp_Vdd(1,7)=eps_now*edge(ii,3)/L(ii,3);

                Jaco_tmp_Vdd(2,1)=coeff_Jn*((edge(ii,1)/L(ii,1))*(n2_Vdd/(-V_T)*Ber_d(x21_Vdd/V_T)-n1_Vdd/(V_T)*Ber_d(-x21_Vdd/V_T))+((edge(ii,3)/L(ii,3))*(n3_Vdd/(-V_T)*Ber_d(x31_Vdd/V_T)-n1_Vdd/(V_T)*Ber_d(-x31_Vdd/V_T))));
                Jaco_tmp_Vdd(2,2)=coeff_Jn*((edge(ii,1)/L(ii,1))*(-Ber(-x21_Vdd/V_T))+((edge(ii,3)/L(ii,3))*(-Ber(-x31_Vdd/V_T))));
                Jaco_tmp_Vdd(2,4)=coeff_Jn*((edge(ii,1)/L(ii,1))*(n2_Vdd/(V_T)*Ber_d(x21_Vdd/V_T)-n1_Vdd/(-V_T)*Ber_d(-x21_Vdd/V_T)));
                Jaco_tmp_Vdd(2,5)=coeff_Jn*((edge(ii,1)/L(ii,1))*(Ber(x21_Vdd/V_T)));
                Jaco_tmp_Vdd(2,7)=coeff_Jn*((edge(ii,3)/L(ii,3))*(n3_Vdd/(V_T)*Ber_d(x31_Vdd/V_T)-n1_Vdd/(-V_T)*Ber_d(-x31_Vdd/V_T)));
                Jaco_tmp_Vdd(2,8)=coeff_Jn*((edge(ii,3)/L(ii,3))*(Ber(x31_Vdd/V_T)));

                Jaco_tmp_Vdd(3,1)=coeff_Jp*((edge(ii,1)/L(ii,1))*(p2_Vdd/(V_T)*Ber_d(-x21_Vdd/V_T)-p1_Vdd/(-V_T)*Ber_d(x21_Vdd/V_T))+((edge(ii,3)/L(ii,3))*(p3_Vdd/(V_T)*Ber_d(-x31_Vdd/V_T)-p1_Vdd/(-V_T)*Ber_d(x31_Vdd/V_T))));
                Jaco_tmp_Vdd(3,3)=coeff_Jp*((edge(ii,1)/L(ii,1))*(-Ber(x21_Vdd/V_T))+(edge(ii,3)/L(ii,3))*(-Ber(x31_Vdd/V_T)));
                Jaco_tmp_Vdd(3,4)=coeff_Jp*((edge(ii,1)/L(ii,1))*(p2_Vdd/(-V_T)*Ber_d(-x21_Vdd/V_T)-p1_Vdd/(V_T)*Ber_d(x21_Vdd/V_T)));
                Jaco_tmp_Vdd(3,6)=coeff_Jp*((edge(ii,1)/L(ii,1))*(Ber(-x21_Vdd/V_T)));
                Jaco_tmp_Vdd(3,7)=coeff_Jp*((edge(ii,3)/L(ii,3))*(p3_Vdd/(-V_T)*Ber_d(-x31_Vdd/V_T)-p1_Vdd/(V_T)*Ber_d(x31_Vdd/V_T)));
                Jaco_tmp_Vdd(3,9)=coeff_Jp*((edge(ii,3)/L(ii,3))*(Ber(-x31_Vdd/V_T)));

                Jaco_tmp_Vdd(4,1)=eps_now*edge(ii,1)/L(ii,1);
                Jaco_tmp_Vdd(4,4)=eps_now*(-edge(ii,1)/L(ii,1)-edge(ii,2)/L(ii,2));
                Jaco_tmp_Vdd(4,5)=-coeff*Control_Volume(2,1);
                Jaco_tmp_Vdd(4,6)=coeff*Control_Volume(2,1);
                Jaco_tmp_Vdd(4,7)=eps_now*edge(ii,2)/L(ii,2);

                Jaco_tmp_Vdd(5,1)=coeff_Jn*((edge(ii,1)/L(ii,1))*(n1_Vdd/(V_T)*Ber_d(x12_Vdd/V_T)-n2_Vdd/(-V_T)*Ber_d(-x12_Vdd/V_T)));
                Jaco_tmp_Vdd(5,2)=coeff_Jn*((edge(ii,1)/L(ii,1))*(Ber(x12_Vdd/V_T)));
                Jaco_tmp_Vdd(5,4)=coeff_Jn*((edge(ii,2)/L(ii,2))*(n3_Vdd/(-V_T)*Ber_d(x32_Vdd/V_T)-n2_Vdd/(V_T)*Ber_d(-x32_Vdd/V_T))+((edge(ii,1)/L(ii,1))*(n1_Vdd/(-V_T)*Ber_d(x12_Vdd/V_T)-n2_Vdd/(V_T)*Ber_d(-x12_Vdd/V_T))));
                Jaco_tmp_Vdd(5,5)=coeff_Jn*((edge(ii,2)/L(ii,2))*(-Ber(-x32_Vdd/V_T))+((edge(ii,1)/L(ii,1))*(-Ber(-x12_Vdd/V_T))));
                Jaco_tmp_Vdd(5,7)=coeff_Jn*((edge(ii,2)/L(ii,2))*(n3_Vdd/(V_T)*Ber_d(x32_Vdd/V_T)-n2_Vdd/(-V_T)*Ber_d(-x32_Vdd/V_T)));
                Jaco_tmp_Vdd(5,8)=coeff_Jn*((edge(ii,2)/L(ii,2))*(Ber(x32_Vdd/V_T)));

                Jaco_tmp_Vdd(6,1)=coeff_Jp*((edge(ii,1)/L(ii,1))*(p1_Vdd/(-V_T)*Ber_d(-x12_Vdd/V_T)-p2_Vdd/(V_T)*Ber_d(x12_Vdd/V_T)));
                Jaco_tmp_Vdd(6,3)=coeff_Jp*((edge(ii,1)/L(ii,1))*(Ber(-x12_Vdd/V_T)));
                Jaco_tmp_Vdd(6,4)=coeff_Jp*((edge(ii,2)/L(ii,2))*(p3_Vdd/(V_T)*Ber_d(-x32_Vdd/V_T)-p2_Vdd/(-V_T)*Ber_d(x32_Vdd/V_T))+((edge(ii,1)/L(ii,1))*(p1_Vdd/(V_T)*Ber_d(-x12_Vdd/V_T)-p2_Vdd/(-V_T)*Ber_d(x12_Vdd/V_T))));
                Jaco_tmp_Vdd(6,6)=coeff_Jp*((edge(ii,2)/L(ii,2))*(-Ber(x32_Vdd/V_T))+(edge(ii,1)/L(ii,1))*(-Ber(x12_Vdd/V_T)));
                Jaco_tmp_Vdd(6,7)=coeff_Jp*((edge(ii,2)/L(ii,2))*(p3_Vdd/(-V_T)*Ber_d(-x32_Vdd/V_T)-p2_Vdd/(V_T)*Ber_d(x32_Vdd/V_T)));
                Jaco_tmp_Vdd(6,9)=coeff_Jp*((edge(ii,2)/L(ii,2))*(Ber(-x32_Vdd/V_T)));

                Jaco_tmp_Vdd(7,1)=eps_now*edge(ii,3)/L(ii,3);
                Jaco_tmp_Vdd(7,4)=eps_now*edge(ii,2)/L(ii,2);
                Jaco_tmp_Vdd(7,7)=eps_now*(-edge(ii,2)/L(ii,2)-edge(ii,3)/L(ii,3));
                Jaco_tmp_Vdd(7,8)=-coeff*Control_Volume(3,1);
                Jaco_tmp_Vdd(7,9)=coeff*Control_Volume(3,1);

                Jaco_tmp_Vdd(8,1)=coeff_Jn*((edge(ii,3)/L(ii,3))*(n1_Vdd/(V_T)*Ber_d(x13_Vdd/V_T)-n3_Vdd/(-V_T)*Ber_d(-x13_Vdd/V_T)));
                Jaco_tmp_Vdd(8,2)=coeff_Jn*((edge(ii,3)/L(ii,3))*(Ber(x13_Vdd/V_T)));
                Jaco_tmp_Vdd(8,4)=coeff_Jn*((edge(ii,2)/L(ii,2))*(n2_Vdd/(V_T)*Ber_d(x23_Vdd/V_T)-n3_Vdd/(-V_T)*Ber_d(-x23_Vdd/V_T)));
                Jaco_tmp_Vdd(8,5)=coeff_Jn*((edge(ii,2)/L(ii,2))*(Ber(x23_Vdd/V_T)));
                Jaco_tmp_Vdd(8,7)=coeff_Jn*((edge(ii,3)/L(ii,3))*(n1_Vdd/(-V_T)*Ber_d(x13_Vdd/V_T)-n3_Vdd/(V_T)*Ber_d(-x13_Vdd/V_T))+((edge(ii,2)/L(ii,2))*(n2_Vdd/(-V_T)*Ber_d(x23_Vdd/V_T)-n3_Vdd/(V_T)*Ber_d(-x23_Vdd/V_T))));
                Jaco_tmp_Vdd(8,8)=coeff_Jn*((edge(ii,3)/L(ii,3))*(-Ber(-x13_Vdd/V_T))+((edge(ii,2)/L(ii,2))*(-Ber(-x23_Vdd/V_T))));

                Jaco_tmp_Vdd(9,1)=coeff_Jp*((edge(ii,3)/L(ii,3))*(p1_Vdd/(-V_T)*Ber_d(-x13_Vdd/V_T)-p3_Vdd/(V_T)*Ber_d(x13_Vdd/V_T)));
                Jaco_tmp_Vdd(9,3)=coeff_Jp*((edge(ii,3)/L(ii,3))*(Ber(-x13_Vdd/V_T)));
                Jaco_tmp_Vdd(9,4)=coeff_Jp*((edge(ii,2)/L(ii,2))*(p2_Vdd/(-V_T)*Ber_d(-x23_Vdd/V_T)-p3_Vdd/(V_T)*Ber_d(x23_Vdd/V_T)));
                Jaco_tmp_Vdd(9,6)=coeff_Jp*((edge(ii,2)/L(ii,2))*(Ber(-x23_Vdd/V_T)));
                Jaco_tmp_Vdd(9,7)=coeff_Jp*((edge(ii,3)/L(ii,3))*(p1_Vdd/(V_T)*Ber_d(-x13_Vdd/V_T)-p3_Vdd/(-V_T)*Ber_d(x13_Vdd/V_T))+((edge(ii,2)/L(ii,2))*(p2_Vdd/(V_T)*Ber_d(-x23_Vdd/V_T)-p3_Vdd/(-V_T)*Ber_d(x23_Vdd/V_T))));
                Jaco_tmp_Vdd(9,9)=coeff_Jp*((edge(ii,3)/L(ii,3))*(-Ber(x13_Vdd/V_T))+(edge(ii,2)/L(ii,2))*(-Ber(x23_Vdd/V_T)));


                for j=1:9
                    for k=1:9
                        Jaco_Vdd(index(j,1),index(k,1))=Jaco_Vdd(index(j,1),index(k,1))+Jaco_tmp_Vdd(j,k);
                    end
                end

                % res vector
                res_tmp_Vdd=zeros(9,1);
                res_potential_tmp1_Vdd(1,1)=eps_now*(edge(ii,1)/L(ii,1)*(x21_Vdd)+edge(ii,3)/L(ii,3)*(x31_Vdd));
                res_potential_tmp1_Vdd(2,1)=eps_now*(edge(ii,2)/L(ii,2)*(x32_Vdd)+edge(ii,1)/L(ii,1)*(x12_Vdd));
                res_potential_tmp1_Vdd(3,1)=eps_now*(edge(ii,3)/L(ii,3)*(x13_Vdd)+edge(ii,2)/L(ii,2)*(x23_Vdd));
                for j=1:3
                    res_potential_tmp2_Vdd(j,1)=(N_dop(Table_element_region(ii,1),1)-solution_vector_Vdd(index(3*j-1,1),1)+solution_vector_Vdd(index(3*j,1),1))*coeff*Control_Volume(j,1);
                end

                res_potential_Vdd=res_potential_tmp1_Vdd+res_potential_tmp2_Vdd;

                % Jn
                res_Jn=zeros(3,1);
                res_Jn(1,1)=coeff_Jn*((edge(ii,1)/L(ii,1))*(n2_Vdd*Ber(x21_Vdd/V_T)-n1_Vdd*Ber(-x21_Vdd/V_T))+(edge(ii,3)/L(ii,3))*(n3_Vdd*Ber(x31_Vdd/V_T)-n1_Vdd*Ber(-x31_Vdd/V_T)));
                res_Jn(2,1)=coeff_Jn*((edge(ii,2)/L(ii,2))*(n3_Vdd*Ber(x32_Vdd/V_T)-n2_Vdd*Ber(-x32_Vdd/V_T))+(edge(ii,1)/L(ii,1))*(n1_Vdd*Ber(x12_Vdd/V_T)-n2_Vdd*Ber(-x12_Vdd/V_T)));
                res_Jn(3,1)=coeff_Jn*((edge(ii,3)/L(ii,3))*(n1_Vdd*Ber(x13_Vdd/V_T)-n3_Vdd*Ber(-x13_Vdd/V_T))+(edge(ii,2)/L(ii,2))*(n2_Vdd*Ber(x23_Vdd/V_T)-n3_Vdd*Ber(-x23_Vdd/V_T)));

                % Jp
                res_Jp=zeros(3,1);
                res_Jp(1,1)=coeff_Jp*((edge(ii,1)/L(ii,1))*(p2_Vdd*Ber(-x21_Vdd/V_T)-p1_Vdd*Ber(x21_Vdd/V_T))+(edge(ii,3)/L(ii,3))*(p3_Vdd*Ber(-x31_Vdd/V_T)-p1_Vdd*Ber(x31_Vdd/V_T)));
                res_Jp(2,1)=coeff_Jp*((edge(ii,2)/L(ii,2))*(p3_Vdd*Ber(-x32_Vdd/V_T)-p2_Vdd*Ber(x32_Vdd/V_T))+(edge(ii,1)/L(ii,1))*(p1_Vdd*Ber(-x12_Vdd/V_T)-p2_Vdd*Ber(x12_Vdd/V_T)));
                res_Jp(3,1)=coeff_Jp*((edge(ii,3)/L(ii,3))*(p1_Vdd*Ber(-x13_Vdd/V_T)-p3_Vdd*Ber(x13_Vdd/V_T))+(edge(ii,2)/L(ii,2))*(p2_Vdd*Ber(-x23_Vdd/V_T)-p3_Vdd*Ber(x23_Vdd/V_T)));

                n=1;
                for j=1:3
                    res_tmp_Vdd(n,1)=res_potential_Vdd(j,1);
                    res_tmp_Vdd(n+1,1)=res_Jn(j,1);
                    res_tmp_Vdd(n+2,1)=res_Jp(j,1);
                    n=n+3;
                end

                for j=1:9
                    res_Vdd(index(j,1),1)=res_Vdd(index(j,1),1)+res_tmp_Vdd(j,1);
                end

            else
                % index : 현재 element에 속해있는 vertex가 resdue에 몇번째 있는가 저장
                index=zeros(3,1);
                n=1;
                for j=1:3
                    index(j,1)=find(Table_Jaco(:,1)==Table_element_region(ii,1) & Table_Jaco(:,2)==Element(ii,j) & Table_Jaco(:,3)==1);
                end

                % Jaco_potential matrix
                Jaco_tmp_Vdd=zeros(3,3);
                Jaco_tmp_Vdd=[(-edge(ii,1)/L(ii,1)-edge(ii,3)/L(ii,3)) edge(ii,1)/L(ii,1) edge(ii,3)/L(ii,3);
                    edge(ii,1)/L(ii,1) -edge(ii,1)/L(ii,1)-edge(ii,2)/L(ii,2) edge(ii,2)/L(ii,2);
                    edge(ii,3)/L(ii,3) edge(ii,2)/L(ii,2) -edge(ii,2)/L(ii,2)-edge(ii,3)/L(ii,3)];
                Jaco_tmp_Vdd=eps_now*Jaco_tmp_Vdd;
                for j=1:3
                    for k=1:3
                        Jaco_Vdd(index(j,1),index(k,1))=Jaco_Vdd(index(j,1),index(k,1))+Jaco_tmp_Vdd(j,k);
                    end
                end

                % res vector
                res_tmp_Vdd=zeros(3,1);
                res_tmp_Vdd=eps_now.*[(-edge(ii,1)/L(ii,1)-edge(ii,3)/L(ii,3))*solution_vector_Vdd(index(1,1),1)+edge(ii,1)/L(ii,1)*solution_vector_Vdd(index(2,1),1)+edge(ii,3)/L(ii,3)*solution_vector_Vdd(index(3,1),1);
                    edge(ii,1)/L(ii,1)*solution_vector_Vdd(index(1,1),1)+(-edge(ii,1)/L(ii,1)-edge(ii,2)/L(ii,2))*solution_vector_Vdd(index(2,1),1)+edge(ii,2)/L(ii,2)*solution_vector_Vdd(index(3,1),1);
                    edge(ii,3)/L(ii,3)*solution_vector_Vdd(index(1,1),1)+edge(ii,2)/L(ii,2)*solution_vector_Vdd(index(2,1),1)+(-edge(ii,2)/L(ii,2)-edge(ii,3)/L(ii,3))*solution_vector_Vdd(index(3,1),1)];
                for j=1:3
                    res_Vdd(index(j,1),1)=res_Vdd(index(j,1),1)+res_tmp_Vdd(j,1);
                end
            end
        end

        % ox 의 값을 si에 더해준 후 ox값은 1,-1로 변경하기
        for ii=[1 5]
            for j=[3 2 4]
                clearvars Vertex_interface_tmp
                Vertex_interface_tmp=Vertex_interface{ii,j};
                for k=1:size(Vertex_interface_tmp,1)
                    index_ox=find(Table_Jaco(:,1)==ii & Table_Jaco(:,2)==Vertex_interface_tmp(k,1) & Table_Jaco(:,3)==1);
                    index_si=find(Table_Jaco(:,1)==j & Table_Jaco(:,2)==Vertex_interface_tmp(k,1) & Table_Jaco(:,3)==1);
                    if length(sort(nonzeros(Jaco_Vdd(index_ox,:))))>2
                        Jaco_Vdd(index_si,:)=Jaco_Vdd(index_si,:)+Jaco_Vdd(index_ox,:);
                        Jaco_Vdd(index_ox,:)=0; Jaco_Vdd(index_ox,index_ox)=1; Jaco_Vdd(index_ox,index_si)=-1;
                        res_Vdd(index_si,1)=res_Vdd(index_si,1)+res_Vdd(index_ox,1);
                        res_Vdd(index_ox,1)=0;
                    elseif ~isequal(sort(nonzeros(Jaco_Vdd(index_ox,:))),[-1; 1]) % 점 3개가 접해있는 지역에서 중복적으로 합쳐지는것을 막기위해 추가함. jaco의 ox-channel행이 -1과1이라면 이미 다른행과 합쳐진 것이므로 pass.
                        Jaco_Vdd(index_si,:)=Jaco_Vdd(index_si,:)+Jaco_Vdd(index_ox,:);
                        Jaco_Vdd(index_ox,:)=0; Jaco_Vdd(index_ox,index_ox)=1; Jaco_Vdd(index_ox,index_si)=-1;
                        res_Vdd(index_si,1)=res_Vdd(index_si,1)+res_Vdd(index_ox,1);
                        res_Vdd(index_ox,1)=0;
                    end
                end
            end
        end
        % source/drain 값을 channel에 더해주기. Si부분이니까 electron/hole density도 모두 바꾸어주어야함.
        for ii=[2 4]
            for j=3
                clearvars Vertex_interface_tmp
                Vertex_interface_tmp=Vertex_interface{ii,j};
                for k=1:size(Vertex_interface_tmp,1)
                    % index 설정
                    index_source_drain_potential=find(Table_Jaco(:,1)==ii & Table_Jaco(:,2)==Vertex_interface_tmp(k,1) & Table_Jaco(:,3)==1);
                    index_source_drain_electron=index_source_drain_potential+1;
                    index_source_drain_hole=index_source_drain_potential+2;
                    index_channel_potential=find(Table_Jaco(:,1)==j & Table_Jaco(:,2)==Vertex_interface_tmp(k,1) & Table_Jaco(:,3)==1);
                    index_channel_electron=index_channel_potential+1;
                    index_channel_hole=index_channel_potential+2;
                    % Jaco 변경, potential/electron/hole 순서
                    Jaco_Vdd(index_channel_potential,:)=Jaco_Vdd(index_channel_potential,:)+Jaco_Vdd(index_source_drain_potential,:);
                    Jaco_Vdd(index_channel_electron,:)=Jaco_Vdd(index_channel_electron,:)+Jaco_Vdd(index_source_drain_electron,:);
                    Jaco_Vdd(index_channel_hole,:)=Jaco_Vdd(index_channel_hole,:)+Jaco_Vdd(index_source_drain_hole,:);
                    Jaco_Vdd(index_source_drain_potential,:)=0; Jaco_Vdd(index_source_drain_potential,index_source_drain_potential)=1; Jaco_Vdd(index_source_drain_potential,index_channel_potential)=-1;
                    Jaco_Vdd(index_source_drain_electron,:)=0; Jaco_Vdd(index_source_drain_electron,index_source_drain_electron)=1; Jaco_Vdd(index_source_drain_electron,index_channel_electron)=-1;
                    Jaco_Vdd(index_source_drain_hole,:)=0; Jaco_Vdd(index_source_drain_hole,index_source_drain_hole)=1; Jaco_Vdd(index_source_drain_hole,index_channel_hole)=-1;

                    % res 변경, potential/electron/hole 순서
                    res_Vdd(index_channel_potential,1)=res_Vdd(index_channel_potential,1)+res_Vdd(index_source_drain_potential,1);
                    res_Vdd(index_channel_electron,1)=res_Vdd(index_channel_electron,1)+res_Vdd(index_source_drain_electron,1);
                    res_Vdd(index_channel_hole,1)=res_Vdd(index_channel_hole,1)+res_Vdd(index_source_drain_hole,1);
                    res_Vdd(index_source_drain_potential,1)=0;
                    res_Vdd(index_source_drain_electron,1)=0;
                    res_Vdd(index_source_drain_hole,1)=0;
                end
            end
        end


        % Boundary condition %
        % Dirichlet_BC
        for ii=1:size(Dirichlet_BC,1)
            BC_index=find(Table_Jaco(:,2)==Dirichlet_BC(ii,1) & Table_Jaco(:,3)==1);
            Jaco_Vdd(BC_index,:)=0;
            Jaco_Vdd(BC_index,BC_index)=1;
            Jaco_Vdd(BC_index,ind_Vg)=-1;
            res_Vdd(BC_index,1)=solution_vector_Vdd(BC_index,1)-V_gate_workfunction-V_gate;
        end

        % Source_BC
        for ii=1:size(Source_BC,1)
            BC_index=find(Table_Jaco(:,1)==2 & Table_Jaco(:,2)==Source_BC(ii,1) & Table_Jaco(:,3)==1);
            Jaco_Vdd(BC_index,:)=0; Jaco_Vdd(BC_index+1,:)=0; Jaco_Vdd(BC_index+2,:)=0;
            Jaco_Vdd(BC_index,BC_index)=1; Jaco_Vdd(BC_index+1,BC_index+1)=1; Jaco_Vdd(BC_index+2,BC_index+2)=1;
            Jaco_Vdd(BC_index,ind_Vs)=-1;
            res_Vdd(BC_index,1)=solution_vector_Vdd(BC_index,1)-V_T*log(Nd/n_int)-solution_vector_Vdd(ind_Vs,1);
            res_Vdd(BC_index+1,1)=solution_vector_Vdd(BC_index+1,1)-Nd;
            res_Vdd(BC_index+2,1)=solution_vector_Vdd(BC_index+2,1)-n_int^2/Nd;
        end

        % Darin_BC
        for ii=1:size(Drain_BC,1)
            BC_index=find(Table_Jaco(:,1)==4 & Table_Jaco(:,2)==Drain_BC(ii,1) & Table_Jaco(:,3)==1);

            % for mixed simulation
            Jaco_Vdd(ind_Id,:)=Jaco_Vdd(ind_Id,:)-Jaco_Vdd(BC_index+1,:)*width;
            Jaco_Vdd(ind_Id,:)=Jaco_Vdd(ind_Id,:)-Jaco_Vdd(BC_index+2,:)*width;
            res_Vdd(ind_Id,:)=res_Vdd(ind_Id,:)-res_Vdd(BC_index+1,:)*width;
            res_Vdd(ind_Id,:)=res_Vdd(ind_Id,:)-res_Vdd(BC_index+2,:)*width;

            % Boundary condition
            Jaco_Vdd(BC_index,:)=0; Jaco_Vdd(BC_index+1,:)=0; Jaco_Vdd(BC_index+2,:)=0;
            Jaco_Vdd(BC_index,BC_index)=1; Jaco_Vdd(BC_index+1,BC_index+1)=1; Jaco_Vdd(BC_index+2,BC_index+2)=1;
            Jaco_Vdd(BC_index,ind_Vd)=-1;
            res_Vdd(BC_index,1)=solution_vector_Vdd(BC_index,1)-V_T*log(Nd/n_int)-solution_vector_Vdd(ind_Vd,1);
            res_Vdd(BC_index+1,1)=solution_vector_Vdd(BC_index+1,1)-Nd;
            res_Vdd(BC_index+2,1)=solution_vector_Vdd(BC_index+2,1)-n_int^2/Nd;
        end

        %% mixed mode simulation
        % Id
        Jaco_Vdd(ind_Id,ind_Id)=1;
        res_Vdd(ind_Id,1)=solution_vector_Vdd(ind_Id,1)+res_Vdd(ind_Id,1);

        % Is
        Jaco_Vdd(ind_Is,ind_Is)=1; Jaco_Vdd(ind_Is,ind_Id)=1; Jaco_Vdd(ind_Is,ind_Ig)=1;
        res_Vdd(ind_Is,1)=0;

        % Ig
        Jaco_Vdd(ind_Ig,ind_Ig)=1;
        res_Vdd(ind_Ig,1)=solution_vector_Vdd(ind_Ig,1)-0;

        % I1
        Jaco_Vdd(ind_I1,ind_I1)=1; Jaco_Vdd(ind_I1,ind_I2)=1;
        res_Vdd(ind_I1,1)=0;

        % I2
        Jaco_Vdd(ind_I2,ind_I2)=1; Jaco_Vdd(ind_I2,ind_Vout)=1/R1; Jaco_Vdd(ind_I2,ind_V2)=-1/R1;
        res_Vdd(ind_I2,1)=solution_vector_Vdd(ind_I2,1)-(solution_vector_Vdd(ind_V2,1)-solution_vector_Vdd(ind_Vout,1))/R1;

        % Vd
        Jaco_Vdd(ind_Vd,ind_Vd)=1; Jaco_Vdd(ind_Vd,ind_Vout)=-1;
        res_Vdd(ind_Vd,1)=0;

        % Vs
        Jaco_Vdd(ind_Vs,ind_Vs)=1;
        res_Vdd(ind_Vs,1)=solution_vector_Vdd(ind_Vs,1)-0;

        % Vg
        Jaco_Vdd(ind_Vg,ind_Vg)=1;
        res_Vdd(ind_Vg,1)=solution_vector_Vdd(ind_Vg,1)-V_gate;

        % V1
        Jaco_Vdd(ind_V1,ind_V1)=1; Jaco_Vdd(ind_V1,ind_Vout)=-1;
        res_Vdd(ind_V1,1)=0;

        % V2
        Jaco_Vdd(ind_V2,ind_V2)=1;
        res_Vdd(ind_V2,1)=solution_vector_Vdd(ind_V2,1)-Vdd;

        % Vout
        Jaco_Vdd(ind_Vout,ind_Id)=1; Jaco_Vdd(ind_Vout,ind_I1)=1;
        res_Vdd(ind_Vout,1)=0;


        % Scaling %
        Table_Jaco_mixed_tmp=zeros(11,3);
        Table_Jaco_mixed_tmp(1:5,3)=4; Table_Jaco_mixed_tmp(6:11,3)=5;
        Table_Jaco_mixed=[Table_Jaco; Table_Jaco_mixed_tmp];
        Cvector_Vdd=zeros(size(solution_vector_Vdd,1),1); v=[V_T; max(abs(N_dop)); max(abs(N_dop)); 1; 1];
        for ii=1:size(solution_vector_Vdd,1)
            Cvector_Vdd(ii,1)=v(Table_Jaco_mixed(ii,3),1);
        end
        Jaco_Vdd=sparse(Jaco_Vdd);
        Cmatrix_Vdd=spdiags(Cvector_Vdd,0,size(solution_vector_Vdd,1),size(solution_vector_Vdd,1));
        Jaco_scaled_Vdd=Jaco_Vdd*Cmatrix_Vdd;
        Rvector_Vdd=1./sum(abs(Jaco_scaled_Vdd),2);
        Rmatrix_Vdd=spdiags(Rvector_Vdd,0,size(solution_vector_Vdd,1),size(solution_vector_Vdd,1));
        Jaco_scaled_Vdd=Rmatrix_Vdd*Jaco_scaled_Vdd;
        res_scaled_Vdd=Rmatrix_Vdd*res_Vdd;
        update_scaled_Vdd=Jaco_scaled_Vdd\(-res_scaled_Vdd);
        update_Vdd(:,Newton_Vdd)=Cmatrix_Vdd*update_scaled_Vdd;

        solution_vector_Vdd(:,1)=solution_vector_Vdd(:,1)+update_Vdd(:,Newton_Vdd);
        solution_vector_Vdd_saved(:,Newton_Vdd+1)=solution_vector_Vdd(:,1);

        % break
        update_Vdd_phi_for_break=update_Vdd(find(Table_Jaco(:,1)==2, 1 , 'first'):3:find(Table_Jaco(:,1)==4, 1 , 'last'),Newton_Vdd);
        fprintf("Max Updatevector=%d\n\n" , max(abs(update_Vdd_phi_for_break)))
        if max(abs(update_Vdd_phi_for_break))<1e-15
            break;
        end

        %     % break
        %     update_Vdd_phi_for_break=update_Vdd(size(solution_vector_poisson,1)+1,Newton_Vdd);
        %     fprintf("Max Updatevector=%d\n\n" , max(abs(update_Vdd_phi_for_break)))
        %     if max(abs(update_Vdd_phi_for_break))<1e-10
        %         break;
        %     end

    end
    Result_Vdd(:,count_Vdd_current)=solution_vector_Vdd(size(solution_vector_poisson,1)+1:size(solution_vector_Vdd,1));


    % Save
    FILENAME = sprintf('./data/homework25_Vd/homework25_Vd_%.2fV.mat' , Vdd);
    save(FILENAME);
    count_Vdd_current=count_Vdd_current+1;
end

%% Visualize - Vd
Visual_solution_vector_Vd=zeros(size(Vertex,1),3);
for ii=1:size(Table_Jaco,1)
    Visual_solution_vector_Vd(Table_Jaco(ii,2),Table_Jaco(ii,3))=solution_vector_Vdd_saved(ii,Newton_Vdd);
end

figure(1) % 점으로 3D
plot3(Vertex(:,1),Vertex(:,2),Visual_solution_vector_Vd(:,1),'*')
figure(2) % 점으로 3D
plot3(Vertex(:,1),Vertex(:,2),Visual_solution_vector_Vd(:,2),'*')
figure(3) % 점으로 3D
plot3(Vertex(:,1),Vertex(:,2),Visual_solution_vector_Vd(:,3),'*')

% figure(2); % mesh 모양
% patch('Faces',Element,'Vertices',Vertex(:,[1 2]),'EdgeColor','black','FaceColor','none','LineWidth',1);
% title('Mesh')
% xlim([0 50*1e-9])
% xlabel('X');
% ylabel('Y');
% hold on
% patch('Faces',edge_interface,'Vertices',Vertex(:,[1 2]),'EdgeColor','green','FaceColor','none','LineWidth',3);
% hold on
% patch('Faces',Element_si_source,'Vertices',Vertex(:,[1 2]), 'EdgeColor','black','FaceColor','yellow','LineWidth',0.5)
% hold on
% patch('Faces',Element_si_drain,'Vertices',Vertex(:,[1 2]), 'EdgeColor','black','FaceColor','yellow','LineWidth',0.5)
% hold on
% patch('Faces',Element_si_channel,'Vertices',Vertex(:,[1 2]), 'EdgeColor','black','FaceColor','red','LineWidth',0.5)
% hold on
% patch('Faces',Element_ox1,'Vertices',Vertex(:,[1 2]), 'EdgeColor','black','FaceColor','blue','LineWidth',0.5)
% hold on
% patch('Faces',Element_ox2,'Vertices',Vertex(:,[1 2]), 'EdgeColor','black','FaceColor','blue','LineWidth',0.5)
% hold on
% patch('Faces',Element_ox2,'Vertices',Vertex(:,[1 2]), 'EdgeColor','black','FaceColor','blue','LineWidth',0.5)
% hold on
% patch('Faces',Element_source_drain,'Vertices',Vertex(:,[1 2]), 'EdgeColor','m','FaceColor','none','LineWidth',5)
% hold on
% patch('Faces',Element_Dirichlet_BC,'Vertices',Vertex(:,[1 2]), 'EdgeColor','g','FaceColor','none','LineWidth',5)
% hold off
%
% figure(3); % potential
% patch('Faces',Element,'Vertices',Vertex(:,[1 2]), 'FaceVertexCData',Visual_solution_vector_Vd(:,1), 'EdgeColor','black','FaceColor','interp');
% title('Initial potential')
% hold on
% patch('Faces',edge_interface,'Vertices',Vertex(:,[1 2]),'EdgeColor','black','FaceColor','none','LineWidth',5);
% hold on
% patch('Faces',Element_source_drain,'Vertices',Vertex(:,[1 2]), 'EdgeColor','m','FaceColor','none','LineWidth',5)
% hold on
% patch('Faces',Element_Dirichlet_BC,'Vertices',Vertex(:,[1 2]), 'EdgeColor','g','FaceColor','none','LineWidth',5)
% hold off
% xlim([0 50*1e-9])
% xlabel('X');
% ylabel('Y');
% colorbar
%
% figure(4); % electron
% patch('Faces',Element,'Vertices',Vertex(:,[1 2]), 'FaceVertexCData',Visual_solution_vector_Vd(:,2), 'EdgeColor','black','FaceColor','interp');
% title('elctron')
% xlim([0 50*1e-9])
% xlabel('X');
% ylabel('Y');
% colorbar
%
% figure(5); % hole
% patch('Faces',Element,'Vertices',Vertex(:,[1 2]), 'FaceVertexCData',Visual_solution_vector_Vd(:,3), 'EdgeColor','black','FaceColor','interp');
% title('hole')
% xlim([0 50*1e-9])
% xlabel('X');
% ylabel('Y');
% colorbar