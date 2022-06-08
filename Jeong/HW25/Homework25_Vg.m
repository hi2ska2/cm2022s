clear; clc; close all;


% Load할 Vd 지정
Vd=0.50;
FILENAME = sprintf('./data/homework25_Vd/homework25_Vd_%.2fV.mat' , Vd);
load(FILENAME)
mkdir('./data/homework25_V_gate')
%% Ramping Gate Voltage

Iteration_Vg=100000;

count_Vg_current=1;
for bias_Vg=1:29
    V_gate=0.05*bias_Vg-0.05;

    % Jacobian matrix / res vector
    solution_vector_Vg(:,1)=solution_vector_Vdd(:,1);
    solution_vector_Vg_saved(:,1)=solution_vector_Vg(:,1);
    for Newton_Vg=1:Iteration_Vg
        fprintf("Ramping Gate Voltage, Vg=%.5f, Newton_Vg=%d\n" , V_gate, Newton_Vg)
        clearvars Jaco_Vg; clearvars res_Vg;
        %     Jaco_Vg=sparse(zeros(size(solution_vector_Vg,1),size(solution_vector_Vg,1))); res_Vg=zeros(size(solution_vector_Vg,1),1);
        Jaco_Vg=zeros(size(solution_vector_Vg,1),size(solution_vector_Vg,1)); res_Vg=zeros(size(solution_vector_Vg,1),1); %% 희소행렬 아님. Jaco matrix 확인용.

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
                x12_Vg=solution_vector_Vg(index(1,1),1)-solution_vector_Vg(index(4,1),1); x21_Vg=-x12_Vg;
                x23_Vg=solution_vector_Vg(index(4,1),1)-solution_vector_Vg(index(7,1),1); x32_Vg=-x23_Vg;
                x13_Vg=solution_vector_Vg(index(1,1),1)-solution_vector_Vg(index(7,1),1); x31_Vg=-x13_Vg;
                n1_Vg=solution_vector_Vg(index(1,1)+1,1); n2_Vg=solution_vector_Vg(index(4,1)+1,1); n3_Vg=solution_vector_Vg(index(7,1)+1,1);
                p1_Vg=solution_vector_Vg(index(1,1)+2,1); p2_Vg=solution_vector_Vg(index(4,1)+2,1); p3_Vg=solution_vector_Vg(index(7,1)+2,1);

                % Jaco %
                Jaco_tmp_Vg=zeros(9,9);
                Jaco_tmp_Vg(1,1)=eps_now*(-edge(ii,1)/L(ii,1)-edge(ii,3)/L(ii,3));
                Jaco_tmp_Vg(1,2)=-coeff*Control_Volume(1,1);
                Jaco_tmp_Vg(1,3)=coeff*Control_Volume(1,1);
                Jaco_tmp_Vg(1,4)=eps_now*edge(ii,1)/L(ii,1);
                Jaco_tmp_Vg(1,7)=eps_now*edge(ii,3)/L(ii,3);

                Jaco_tmp_Vg(2,1)=coeff_Jn*((edge(ii,1)/L(ii,1))*(n2_Vg/(-V_T)*Ber_d(x21_Vg/V_T)-n1_Vg/(V_T)*Ber_d(-x21_Vg/V_T))+((edge(ii,3)/L(ii,3))*(n3_Vg/(-V_T)*Ber_d(x31_Vg/V_T)-n1_Vg/(V_T)*Ber_d(-x31_Vg/V_T))));
                Jaco_tmp_Vg(2,2)=coeff_Jn*((edge(ii,1)/L(ii,1))*(-Ber(-x21_Vg/V_T))+((edge(ii,3)/L(ii,3))*(-Ber(-x31_Vg/V_T))));
                Jaco_tmp_Vg(2,4)=coeff_Jn*((edge(ii,1)/L(ii,1))*(n2_Vg/(V_T)*Ber_d(x21_Vg/V_T)-n1_Vg/(-V_T)*Ber_d(-x21_Vg/V_T)));
                Jaco_tmp_Vg(2,5)=coeff_Jn*((edge(ii,1)/L(ii,1))*(Ber(x21_Vg/V_T)));
                Jaco_tmp_Vg(2,7)=coeff_Jn*((edge(ii,3)/L(ii,3))*(n3_Vg/(V_T)*Ber_d(x31_Vg/V_T)-n1_Vg/(-V_T)*Ber_d(-x31_Vg/V_T)));
                Jaco_tmp_Vg(2,8)=coeff_Jn*((edge(ii,3)/L(ii,3))*(Ber(x31_Vg/V_T)));

                Jaco_tmp_Vg(3,1)=coeff_Jp*((edge(ii,1)/L(ii,1))*(p2_Vg/(V_T)*Ber_d(-x21_Vg/V_T)-p1_Vg/(-V_T)*Ber_d(x21_Vg/V_T))+((edge(ii,3)/L(ii,3))*(p3_Vg/(V_T)*Ber_d(-x31_Vg/V_T)-p1_Vg/(-V_T)*Ber_d(x31_Vg/V_T))));
                Jaco_tmp_Vg(3,3)=coeff_Jp*((edge(ii,1)/L(ii,1))*(-Ber(x21_Vg/V_T))+(edge(ii,3)/L(ii,3))*(-Ber(x31_Vg/V_T)));
                Jaco_tmp_Vg(3,4)=coeff_Jp*((edge(ii,1)/L(ii,1))*(p2_Vg/(-V_T)*Ber_d(-x21_Vg/V_T)-p1_Vg/(V_T)*Ber_d(x21_Vg/V_T)));
                Jaco_tmp_Vg(3,6)=coeff_Jp*((edge(ii,1)/L(ii,1))*(Ber(-x21_Vg/V_T)));
                Jaco_tmp_Vg(3,7)=coeff_Jp*((edge(ii,3)/L(ii,3))*(p3_Vg/(-V_T)*Ber_d(-x31_Vg/V_T)-p1_Vg/(V_T)*Ber_d(x31_Vg/V_T)));
                Jaco_tmp_Vg(3,9)=coeff_Jp*((edge(ii,3)/L(ii,3))*(Ber(-x31_Vg/V_T)));

                Jaco_tmp_Vg(4,1)=eps_now*edge(ii,1)/L(ii,1);
                Jaco_tmp_Vg(4,4)=eps_now*(-edge(ii,1)/L(ii,1)-edge(ii,2)/L(ii,2));
                Jaco_tmp_Vg(4,5)=-coeff*Control_Volume(2,1);
                Jaco_tmp_Vg(4,6)=coeff*Control_Volume(2,1);
                Jaco_tmp_Vg(4,7)=eps_now*edge(ii,2)/L(ii,2);

                Jaco_tmp_Vg(5,1)=coeff_Jn*((edge(ii,1)/L(ii,1))*(n1_Vg/(V_T)*Ber_d(x12_Vg/V_T)-n2_Vg/(-V_T)*Ber_d(-x12_Vg/V_T)));
                Jaco_tmp_Vg(5,2)=coeff_Jn*((edge(ii,1)/L(ii,1))*(Ber(x12_Vg/V_T)));
                Jaco_tmp_Vg(5,4)=coeff_Jn*((edge(ii,2)/L(ii,2))*(n3_Vg/(-V_T)*Ber_d(x32_Vg/V_T)-n2_Vg/(V_T)*Ber_d(-x32_Vg/V_T))+((edge(ii,1)/L(ii,1))*(n1_Vg/(-V_T)*Ber_d(x12_Vg/V_T)-n2_Vg/(V_T)*Ber_d(-x12_Vg/V_T))));
                Jaco_tmp_Vg(5,5)=coeff_Jn*((edge(ii,2)/L(ii,2))*(-Ber(-x32_Vg/V_T))+((edge(ii,1)/L(ii,1))*(-Ber(-x12_Vg/V_T))));
                Jaco_tmp_Vg(5,7)=coeff_Jn*((edge(ii,2)/L(ii,2))*(n3_Vg/(V_T)*Ber_d(x32_Vg/V_T)-n2_Vg/(-V_T)*Ber_d(-x32_Vg/V_T)));
                Jaco_tmp_Vg(5,8)=coeff_Jn*((edge(ii,2)/L(ii,2))*(Ber(x32_Vg/V_T)));

                Jaco_tmp_Vg(6,1)=coeff_Jp*((edge(ii,1)/L(ii,1))*(p1_Vg/(-V_T)*Ber_d(-x12_Vg/V_T)-p2_Vg/(V_T)*Ber_d(x12_Vg/V_T)));
                Jaco_tmp_Vg(6,3)=coeff_Jp*((edge(ii,1)/L(ii,1))*(Ber(-x12_Vg/V_T)));
                Jaco_tmp_Vg(6,4)=coeff_Jp*((edge(ii,2)/L(ii,2))*(p3_Vg/(V_T)*Ber_d(-x32_Vg/V_T)-p2_Vg/(-V_T)*Ber_d(x32_Vg/V_T))+((edge(ii,1)/L(ii,1))*(p1_Vg/(V_T)*Ber_d(-x12_Vg/V_T)-p2_Vg/(-V_T)*Ber_d(x12_Vg/V_T))));
                Jaco_tmp_Vg(6,6)=coeff_Jp*((edge(ii,2)/L(ii,2))*(-Ber(x32_Vg/V_T))+(edge(ii,1)/L(ii,1))*(-Ber(x12_Vg/V_T)));
                Jaco_tmp_Vg(6,7)=coeff_Jp*((edge(ii,2)/L(ii,2))*(p3_Vg/(-V_T)*Ber_d(-x32_Vg/V_T)-p2_Vg/(V_T)*Ber_d(x32_Vg/V_T)));
                Jaco_tmp_Vg(6,9)=coeff_Jp*((edge(ii,2)/L(ii,2))*(Ber(-x32_Vg/V_T)));

                Jaco_tmp_Vg(7,1)=eps_now*edge(ii,3)/L(ii,3);
                Jaco_tmp_Vg(7,4)=eps_now*edge(ii,2)/L(ii,2);
                Jaco_tmp_Vg(7,7)=eps_now*(-edge(ii,2)/L(ii,2)-edge(ii,3)/L(ii,3));
                Jaco_tmp_Vg(7,8)=-coeff*Control_Volume(3,1);
                Jaco_tmp_Vg(7,9)=coeff*Control_Volume(3,1);

                Jaco_tmp_Vg(8,1)=coeff_Jn*((edge(ii,3)/L(ii,3))*(n1_Vg/(V_T)*Ber_d(x13_Vg/V_T)-n3_Vg/(-V_T)*Ber_d(-x13_Vg/V_T)));
                Jaco_tmp_Vg(8,2)=coeff_Jn*((edge(ii,3)/L(ii,3))*(Ber(x13_Vg/V_T)));
                Jaco_tmp_Vg(8,4)=coeff_Jn*((edge(ii,2)/L(ii,2))*(n2_Vg/(V_T)*Ber_d(x23_Vg/V_T)-n3_Vg/(-V_T)*Ber_d(-x23_Vg/V_T)));
                Jaco_tmp_Vg(8,5)=coeff_Jn*((edge(ii,2)/L(ii,2))*(Ber(x23_Vg/V_T)));
                Jaco_tmp_Vg(8,7)=coeff_Jn*((edge(ii,3)/L(ii,3))*(n1_Vg/(-V_T)*Ber_d(x13_Vg/V_T)-n3_Vg/(V_T)*Ber_d(-x13_Vg/V_T))+((edge(ii,2)/L(ii,2))*(n2_Vg/(-V_T)*Ber_d(x23_Vg/V_T)-n3_Vg/(V_T)*Ber_d(-x23_Vg/V_T))));
                Jaco_tmp_Vg(8,8)=coeff_Jn*((edge(ii,3)/L(ii,3))*(-Ber(-x13_Vg/V_T))+((edge(ii,2)/L(ii,2))*(-Ber(-x23_Vg/V_T))));

                Jaco_tmp_Vg(9,1)=coeff_Jp*((edge(ii,3)/L(ii,3))*(p1_Vg/(-V_T)*Ber_d(-x13_Vg/V_T)-p3_Vg/(V_T)*Ber_d(x13_Vg/V_T)));
                Jaco_tmp_Vg(9,3)=coeff_Jp*((edge(ii,3)/L(ii,3))*(Ber(-x13_Vg/V_T)));
                Jaco_tmp_Vg(9,4)=coeff_Jp*((edge(ii,2)/L(ii,2))*(p2_Vg/(-V_T)*Ber_d(-x23_Vg/V_T)-p3_Vg/(V_T)*Ber_d(x23_Vg/V_T)));
                Jaco_tmp_Vg(9,6)=coeff_Jp*((edge(ii,2)/L(ii,2))*(Ber(-x23_Vg/V_T)));
                Jaco_tmp_Vg(9,7)=coeff_Jp*((edge(ii,3)/L(ii,3))*(p1_Vg/(V_T)*Ber_d(-x13_Vg/V_T)-p3_Vg/(-V_T)*Ber_d(x13_Vg/V_T))+((edge(ii,2)/L(ii,2))*(p2_Vg/(V_T)*Ber_d(-x23_Vg/V_T)-p3_Vg/(-V_T)*Ber_d(x23_Vg/V_T))));
                Jaco_tmp_Vg(9,9)=coeff_Jp*((edge(ii,3)/L(ii,3))*(-Ber(x13_Vg/V_T))+(edge(ii,2)/L(ii,2))*(-Ber(x23_Vg/V_T)));


                for j=1:9
                    for k=1:9
                        Jaco_Vg(index(j,1),index(k,1))=Jaco_Vg(index(j,1),index(k,1))+Jaco_tmp_Vg(j,k);
                    end
                end

                % res vector
                res_tmp_Vg=zeros(9,1);
                res_potential_tmp1_Vg(1,1)=eps_now*(edge(ii,1)/L(ii,1)*(x21_Vg)+edge(ii,3)/L(ii,3)*(x31_Vg));
                res_potential_tmp1_Vg(2,1)=eps_now*(edge(ii,2)/L(ii,2)*(x32_Vg)+edge(ii,1)/L(ii,1)*(x12_Vg));
                res_potential_tmp1_Vg(3,1)=eps_now*(edge(ii,3)/L(ii,3)*(x13_Vg)+edge(ii,2)/L(ii,2)*(x23_Vg));
                for j=1:3
                    res_potential_tmp2_Vg(j,1)=(N_dop(Table_element_region(ii,1),1)-solution_vector_Vg(index(3*j-1,1),1)+solution_vector_Vg(index(3*j,1),1))*coeff*Control_Volume(j,1);
                end

                res_potential_Vg=res_potential_tmp1_Vg+res_potential_tmp2_Vg;

                % Jn
                res_Jn=zeros(3,1);
                res_Jn(1,1)=coeff_Jn*((edge(ii,1)/L(ii,1))*(n2_Vg*Ber(x21_Vg/V_T)-n1_Vg*Ber(-x21_Vg/V_T))+(edge(ii,3)/L(ii,3))*(n3_Vg*Ber(x31_Vg/V_T)-n1_Vg*Ber(-x31_Vg/V_T)));
                res_Jn(2,1)=coeff_Jn*((edge(ii,2)/L(ii,2))*(n3_Vg*Ber(x32_Vg/V_T)-n2_Vg*Ber(-x32_Vg/V_T))+(edge(ii,1)/L(ii,1))*(n1_Vg*Ber(x12_Vg/V_T)-n2_Vg*Ber(-x12_Vg/V_T)));
                res_Jn(3,1)=coeff_Jn*((edge(ii,3)/L(ii,3))*(n1_Vg*Ber(x13_Vg/V_T)-n3_Vg*Ber(-x13_Vg/V_T))+(edge(ii,2)/L(ii,2))*(n2_Vg*Ber(x23_Vg/V_T)-n3_Vg*Ber(-x23_Vg/V_T)));

                % Jp
                res_Jp=zeros(3,1);
                res_Jp(1,1)=coeff_Jp*((edge(ii,1)/L(ii,1))*(p2_Vg*Ber(-x21_Vg/V_T)-p1_Vg*Ber(x21_Vg/V_T))+(edge(ii,3)/L(ii,3))*(p3_Vg*Ber(-x31_Vg/V_T)-p1_Vg*Ber(x31_Vg/V_T)));
                res_Jp(2,1)=coeff_Jp*((edge(ii,2)/L(ii,2))*(p3_Vg*Ber(-x32_Vg/V_T)-p2_Vg*Ber(x32_Vg/V_T))+(edge(ii,1)/L(ii,1))*(p1_Vg*Ber(-x12_Vg/V_T)-p2_Vg*Ber(x12_Vg/V_T)));
                res_Jp(3,1)=coeff_Jp*((edge(ii,3)/L(ii,3))*(p1_Vg*Ber(-x13_Vg/V_T)-p3_Vg*Ber(x13_Vg/V_T))+(edge(ii,2)/L(ii,2))*(p2_Vg*Ber(-x23_Vg/V_T)-p3_Vg*Ber(x23_Vg/V_T)));

                n=1;
                for j=1:3
                    res_tmp_Vg(n,1)=res_potential_Vg(j,1);
                    res_tmp_Vg(n+1,1)=res_Jn(j,1);
                    res_tmp_Vg(n+2,1)=res_Jp(j,1);
                    n=n+3;
                end

                for j=1:9
                    res_Vg(index(j,1),1)=res_Vg(index(j,1),1)+res_tmp_Vg(j,1);
                end

            else
                % index : 현재 element에 속해있는 vertex가 resdue에 몇번째 있는가 저장
                index=zeros(3,1);
                n=1;
                for j=1:3
                    index(j,1)=find(Table_Jaco(:,1)==Table_element_region(ii,1) & Table_Jaco(:,2)==Element(ii,j) & Table_Jaco(:,3)==1);
                end

                % Jaco_potential matrix
                Jaco_tmp_Vg=zeros(3,3);
                Jaco_tmp_Vg=[(-edge(ii,1)/L(ii,1)-edge(ii,3)/L(ii,3)) edge(ii,1)/L(ii,1) edge(ii,3)/L(ii,3);
                    edge(ii,1)/L(ii,1) -edge(ii,1)/L(ii,1)-edge(ii,2)/L(ii,2) edge(ii,2)/L(ii,2);
                    edge(ii,3)/L(ii,3) edge(ii,2)/L(ii,2) -edge(ii,2)/L(ii,2)-edge(ii,3)/L(ii,3)];
                Jaco_tmp_Vg=eps_now*Jaco_tmp_Vg;
                for j=1:3
                    for k=1:3
                        Jaco_Vg(index(j,1),index(k,1))=Jaco_Vg(index(j,1),index(k,1))+Jaco_tmp_Vg(j,k);
                    end
                end

                % res vector
                res_tmp_Vg=zeros(3,1);
                res_tmp_Vg=eps_now.*[(-edge(ii,1)/L(ii,1)-edge(ii,3)/L(ii,3))*solution_vector_Vg(index(1,1),1)+edge(ii,1)/L(ii,1)*solution_vector_Vg(index(2,1),1)+edge(ii,3)/L(ii,3)*solution_vector_Vg(index(3,1),1);
                    edge(ii,1)/L(ii,1)*solution_vector_Vg(index(1,1),1)+(-edge(ii,1)/L(ii,1)-edge(ii,2)/L(ii,2))*solution_vector_Vg(index(2,1),1)+edge(ii,2)/L(ii,2)*solution_vector_Vg(index(3,1),1);
                    edge(ii,3)/L(ii,3)*solution_vector_Vg(index(1,1),1)+edge(ii,2)/L(ii,2)*solution_vector_Vg(index(2,1),1)+(-edge(ii,2)/L(ii,2)-edge(ii,3)/L(ii,3))*solution_vector_Vg(index(3,1),1)];
                for j=1:3
                    res_Vg(index(j,1),1)=res_Vg(index(j,1),1)+res_tmp_Vg(j,1);
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
                    if length(sort(nonzeros(Jaco_Vg(index_ox,:))))>2
                        Jaco_Vg(index_si,:)=Jaco_Vg(index_si,:)+Jaco_Vg(index_ox,:);
                        Jaco_Vg(index_ox,:)=0; Jaco_Vg(index_ox,index_ox)=1; Jaco_Vg(index_ox,index_si)=-1;
                        res_Vg(index_si,1)=res_Vg(index_si,1)+res_Vg(index_ox,1);
                        res_Vg(index_ox,1)=0;
                    elseif ~isequal(sort(nonzeros(Jaco_Vg(index_ox,:))),[-1; 1]) % 점 3개가 접해있는 지역에서 중복적으로 합쳐지는것을 막기위해 추가함. jaco의 ox-channel행이 -1과1이라면 이미 다른행과 합쳐진 것이므로 pass.
                        Jaco_Vg(index_si,:)=Jaco_Vg(index_si,:)+Jaco_Vg(index_ox,:);
                        Jaco_Vg(index_ox,:)=0; Jaco_Vg(index_ox,index_ox)=1; Jaco_Vg(index_ox,index_si)=-1;
                        res_Vg(index_si,1)=res_Vg(index_si,1)+res_Vg(index_ox,1);
                        res_Vg(index_ox,1)=0;
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
                    Jaco_Vg(index_channel_potential,:)=Jaco_Vg(index_channel_potential,:)+Jaco_Vg(index_source_drain_potential,:);
                    Jaco_Vg(index_channel_electron,:)=Jaco_Vg(index_channel_electron,:)+Jaco_Vg(index_source_drain_electron,:);
                    Jaco_Vg(index_channel_hole,:)=Jaco_Vg(index_channel_hole,:)+Jaco_Vg(index_source_drain_hole,:);
                    Jaco_Vg(index_source_drain_potential,:)=0; Jaco_Vg(index_source_drain_potential,index_source_drain_potential)=1; Jaco_Vg(index_source_drain_potential,index_channel_potential)=-1;
                    Jaco_Vg(index_source_drain_electron,:)=0; Jaco_Vg(index_source_drain_electron,index_source_drain_electron)=1; Jaco_Vg(index_source_drain_electron,index_channel_electron)=-1;
                    Jaco_Vg(index_source_drain_hole,:)=0; Jaco_Vg(index_source_drain_hole,index_source_drain_hole)=1; Jaco_Vg(index_source_drain_hole,index_channel_hole)=-1;

                    % res 변경, potential/electron/hole 순서
                    res_Vg(index_channel_potential,1)=res_Vg(index_channel_potential,1)+res_Vg(index_source_drain_potential,1);
                    res_Vg(index_channel_electron,1)=res_Vg(index_channel_electron,1)+res_Vg(index_source_drain_electron,1);
                    res_Vg(index_channel_hole,1)=res_Vg(index_channel_hole,1)+res_Vg(index_source_drain_hole,1);
                    res_Vg(index_source_drain_potential,1)=0;
                    res_Vg(index_source_drain_electron,1)=0;
                    res_Vg(index_source_drain_hole,1)=0;
                end
            end
        end


        % Boundary condition %
        % Dirichlet_BC
        for ii=1:size(Dirichlet_BC,1)
            BC_index=find(Table_Jaco(:,2)==Dirichlet_BC(ii,1) & Table_Jaco(:,3)==1);
            Jaco_Vg(BC_index,:)=0;
            Jaco_Vg(BC_index,BC_index)=1;
            Jaco_Vg(BC_index,ind_Vg)=-1;
            res_Vg(BC_index,1)=solution_vector_Vg(BC_index,1)-V_gate_workfunction-V_gate;
        end

        % Source_BC
        for ii=1:size(Source_BC,1)
            BC_index=find(Table_Jaco(:,1)==2 & Table_Jaco(:,2)==Source_BC(ii,1) & Table_Jaco(:,3)==1);
            Jaco_Vg(BC_index,:)=0; Jaco_Vg(BC_index+1,:)=0; Jaco_Vg(BC_index+2,:)=0;
            Jaco_Vg(BC_index,BC_index)=1; Jaco_Vg(BC_index+1,BC_index+1)=1; Jaco_Vg(BC_index+2,BC_index+2)=1;
            Jaco_Vg(BC_index,ind_Vs)=-1;
            res_Vg(BC_index,1)=solution_vector_Vg(BC_index,1)-V_T*log(Nd/n_int)-solution_vector_Vg(ind_Vs,1);
            res_Vg(BC_index+1,1)=solution_vector_Vg(BC_index+1,1)-Nd;
            res_Vg(BC_index+2,1)=solution_vector_Vg(BC_index+2,1)-n_int^2/Nd;
        end

        % Darin_BC
        for ii=1:size(Drain_BC,1)
            BC_index=find(Table_Jaco(:,1)==4 & Table_Jaco(:,2)==Drain_BC(ii,1) & Table_Jaco(:,3)==1);

            % for mixed simulation
            Jaco_Vg(ind_Id,:)=Jaco_Vg(ind_Id,:)-Jaco_Vg(BC_index+1,:)*width;
            Jaco_Vg(ind_Id,:)=Jaco_Vg(ind_Id,:)-Jaco_Vg(BC_index+2,:)*width;
            res_Vg(ind_Id,:)=res_Vg(ind_Id,:)-res_Vg(BC_index+1,:)*width;
            res_Vg(ind_Id,:)=res_Vg(ind_Id,:)-res_Vg(BC_index+2,:)*width;

            % Boundary condition
            Jaco_Vg(BC_index,:)=0; Jaco_Vg(BC_index+1,:)=0; Jaco_Vg(BC_index+2,:)=0;
            Jaco_Vg(BC_index,BC_index)=1; Jaco_Vg(BC_index+1,BC_index+1)=1; Jaco_Vg(BC_index+2,BC_index+2)=1;
            Jaco_Vg(BC_index,ind_Vd)=-1;
            res_Vg(BC_index,1)=solution_vector_Vg(BC_index,1)-V_T*log(Nd/n_int)-solution_vector_Vg(ind_Vd,1);
            res_Vg(BC_index+1,1)=solution_vector_Vg(BC_index+1,1)-Nd;
            res_Vg(BC_index+2,1)=solution_vector_Vg(BC_index+2,1)-n_int^2/Nd;
        end

        %% mixed mode simulation
        % Id
        Jaco_Vg(ind_Id,ind_Id)=1;
        res_Vg(ind_Id,1)=solution_vector_Vg(ind_Id,1)+res_Vg(ind_Id,1);

        % Is
        Jaco_Vg(ind_Is,ind_Is)=1; Jaco_Vg(ind_Is,ind_Id)=1; Jaco_Vg(ind_Is,ind_Ig)=1;
        res_Vg(ind_Is,1)=0;

        % Ig
        Jaco_Vg(ind_Ig,ind_Ig)=1;
        res_Vg(ind_Ig,1)=solution_vector_Vg(ind_Ig,1)-0;

        % I1
        Jaco_Vg(ind_I1,ind_I1)=1; Jaco_Vg(ind_I1,ind_I2)=1;
        res_Vg(ind_I1,1)=0;

        % I2
        Jaco_Vg(ind_I2,ind_I2)=1; Jaco_Vg(ind_I2,ind_Vout)=1/R1; Jaco_Vg(ind_I2,ind_V2)=-1/R1;
        res_Vg(ind_I2,1)=solution_vector_Vg(ind_I2,1)-(solution_vector_Vg(ind_V2,1)-solution_vector_Vg(ind_Vout,1))/R1;

        % Vd
        Jaco_Vg(ind_Vd,ind_Vd)=1; Jaco_Vg(ind_Vd,ind_Vout)=-1;
        res_Vg(ind_Vd,1)=0;

        % Vs
        Jaco_Vg(ind_Vs,ind_Vs)=1;
        res_Vg(ind_Vs,1)=solution_vector_Vg(ind_Vs,1)-0;

        % Vg
        Jaco_Vg(ind_Vg,ind_Vg)=1;
        res_Vg(ind_Vg,1)=solution_vector_Vg(ind_Vg,1)-V_gate;

        % V1
        Jaco_Vg(ind_V1,ind_V1)=1; Jaco_Vg(ind_V1,ind_Vout)=-1;
        res_Vg(ind_V1,1)=0;

        % V2
        Jaco_Vg(ind_V2,ind_V2)=1;
        res_Vg(ind_V2,1)=solution_vector_Vg(ind_V2,1)-Vdd;

        % Vout
        Jaco_Vg(ind_Vout,ind_Id)=1; Jaco_Vg(ind_Vout,ind_I1)=1;
        res_Vg(ind_Vout,1)=0;


        % Scaling %
        Table_Jaco_mixed_tmp=zeros(11,3);
        Table_Jaco_mixed_tmp(1:5,3)=4; Table_Jaco_mixed_tmp(6:11,3)=5;
        Table_Jaco_mixed=[Table_Jaco; Table_Jaco_mixed_tmp];
        Cvector_Vg=zeros(size(solution_vector_Vg,1),1); v=[V_T; max(abs(N_dop)); max(abs(N_dop)); 1; 1];
        for ii=1:size(solution_vector_Vg,1)
            Cvector_Vg(ii,1)=v(Table_Jaco_mixed(ii,3),1);
        end
        Jaco_Vg=sparse(Jaco_Vg);
        Cmatrix_Vg=spdiags(Cvector_Vg,0,size(solution_vector_Vg,1),size(solution_vector_Vg,1));
        Jaco_scaled_Vg=Jaco_Vg*Cmatrix_Vg;
        Rvector_Vg=1./sum(abs(Jaco_scaled_Vg),2);
        Rmatrix_Vg=spdiags(Rvector_Vg,0,size(solution_vector_Vg,1),size(solution_vector_Vg,1));
        Jaco_scaled_Vg=Rmatrix_Vg*Jaco_scaled_Vg;
        res_scaled_Vg=Rmatrix_Vg*res_Vg;
        update_scaled_Vg=Jaco_scaled_Vg\(-res_scaled_Vg);
        update_Vg(:,Newton_Vg)=Cmatrix_Vg*update_scaled_Vg;

        solution_vector_Vg(:,1)=solution_vector_Vg(:,1)+update_Vg(:,Newton_Vg);
        solution_vector_Vg_saved(:,Newton_Vg+1)=solution_vector_Vg(:,1);

        % break
        update_Vg_phi_for_break=update_Vg(find(Table_Jaco(:,1)==2, 1 , 'first'):3:find(Table_Jaco(:,1)==4, 1 , 'last'),Newton_Vg);
        fprintf("Max Updatevector=%d\n\n" , max(abs(update_Vg_phi_for_break)))
        if max(abs(update_Vg_phi_for_break))<1e-15
            break;
        end

        %     % break
        %     update_Vg_phi_for_break=update_Vg(size(solution_vector_poisson,1)+1,Newton_Vg);
        %     fprintf("Max Updatevector=%d\n\n" , max(abs(update_Vg_phi_for_break)))
        %     if max(abs(update_Vg_phi_for_break))<1e-10
        %         break;
        %     end

    end
    Result_Vg(:,count_Vg_current)=solution_vector_Vg(size(solution_vector_poisson,1)+1:size(solution_vector_Vg,1));


    %% save
    FILENAME = sprintf('./data/homework25_V_gate/homework25_Vd_%.2fV_V_gate_%.5fV.mat' ,Vd, V_gate);
    save(FILENAME);
    count_Vg_current=count_Vg_current+1;
end

%% visualize -Vg
%     while str==y || str==Y
%         close all
%         prompt = 'Please enter the time? ';
%         time_input=input(prompt);
% 
%         for ii=1:length(time)
%             if abs(time(ii,1)-time_input)<=1e-15
%                 TIME=ii;
%             end
%         end
% 
% 
%         Visual_solution_vector_Vg=zeros(size(Vertex,1),3);
%         for ii=1:size(Table_Jaco,1)
%             Visual_solution_vector_Vg(Table_Jaco(ii,2),Table_Jaco(ii,3))=solution_vector_Vg_saved(ii,TIME);
%         end
% 
%         figure(1); % hole
%         patch('Faces',Element,'Vertices',Vertex(:,[1 2]), 'FaceVertexCData',Visual_solution_vector_Vg(:,3), 'EdgeColor','black','FaceColor','interp');
%         title('hole')
%         xlabel('X');
%         ylabel('Y');
%         colorbar
% 
%         figure(2); % electron
%         patch('Faces',Element,'Vertices',Vertex(:,[1 2]), 'FaceVertexCData',Visual_solution_vector_Vg(:,2), 'EdgeColor','black','FaceColor','interp');
%         title('elctron')
%         xlabel('X');
%         ylabel('Y');
%         colorbar
% 
%         figure(3); % potential
%         patch('Faces',Element,'Vertices',Vertex(:,[1 2]), 'FaceVertexCData',Visual_solution_vector_Vg(:,1), 'EdgeColor','black','FaceColor','interp');
%         title('Initial potential')
%         hold on
%         patch('Faces',edge_interface,'Vertices',Vertex(:,[1 2]),'EdgeColor','black','FaceColor','none','LineWidth',5);
%         hold on
%         patch('Faces',Element_source_drain,'Vertices',Vertex(:,[1 2]), 'EdgeColor','m','FaceColor','none','LineWidth',5)
%         hold on
%         patch('Faces',Element_Dirichlet_BC,'Vertices',Vertex(:,[1 2]), 'EdgeColor','g','FaceColor','none','LineWidth',5)
%         hold off
%         xlabel('X');
%         ylabel('Y');
%         colorbar
% 
%         figure(4); % mesh 모양
%         patch('Faces',Element,'Vertices',Vertex(:,[1 2]),'EdgeColor','black','FaceColor','none','LineWidth',1);
%         title('Mesh')
%         xlabel('X');
%         ylabel('Y');
%         hold on
%         patch('Faces',edge_interface,'Vertices',Vertex(:,[1 2]),'EdgeColor','green','FaceColor','none','LineWidth',3);
%         hold on
%         patch('Faces',Element_si_source,'Vertices',Vertex(:,[1 2]), 'EdgeColor','black','FaceColor','yellow','LineWidth',0.5)
%         hold on
%         patch('Faces',Element_si_drain,'Vertices',Vertex(:,[1 2]), 'EdgeColor','black','FaceColor','yellow','LineWidth',0.5)
%         hold on
%         patch('Faces',Element_si_channel,'Vertices',Vertex(:,[1 2]), 'EdgeColor','black','FaceColor','red','LineWidth',0.5)
%         hold on
%         patch('Faces',Element_ox1,'Vertices',Vertex(:,[1 2]), 'EdgeColor','black','FaceColor','blue','LineWidth',0.5)
%         hold on
%         patch('Faces',Element_ox2,'Vertices',Vertex(:,[1 2]), 'EdgeColor','black','FaceColor','blue','LineWidth',0.5)
%         hold on
%         patch('Faces',Element_ox2,'Vertices',Vertex(:,[1 2]), 'EdgeColor','black','FaceColor','blue','LineWidth',0.5)
%         hold on
%         patch('Faces',Element_source_drain,'Vertices',Vertex(:,[1 2]), 'EdgeColor','m','FaceColor','none','LineWidth',5)
%         hold on
%         patch('Faces',Element_Dirichlet_BC,'Vertices',Vertex(:,[1 2]), 'EdgeColor','g','FaceColor','none','LineWidth',5)
%         hold off
% 
%         figure(5) % 점으로 3D
%         plot3(Vertex(:,1),Vertex(:,2),Visual_solution_vector_Vg(:,1),'*')
% 
%         figure(6);
%         subplot(1,2,1)
%         plot(time,Vg,'-o')
%         title('Gate Voltage vs time')
%         xlabel('Time [ns]');
%         ylabel('Gate Voltage [V]');
% 
%         subplot(1,2,2)
%         plot(time,I,'-o')
%         title('Drain Current vs time')
%         xlabel('Time [ns]');
%         ylabel('Drain Current [V]');
% 
% 
%         figure(8) % current
%         plot(time, I )
%         title('Drain Current vs time')
%         xlabel('Time [s]');
%         ylabel('Drain Current [V]');
%         xlim=([0.8*1e-10 1.2*1e-10]);
% 
%         figure(7); % gate
%         plot(time, Vg)
%         title('Gate Voltage vs time')
%         xlabel('Time [s]');
%         ylabel('Gate Voltage [V]');
% 
%         prompt = 'Do you want more? Y/N [Y]: ';
%         str = input(prompt,'s');
%         if isempty(str)
%             str = 'Y';
%         end
% 
%     end