clear; clc; close all;
load('./data/Homework16_dd.mat')
mkdir('./data/Homework16_Vd')
%% Ramping Drain Voltage
Iteration_Vd=40;

% Jacobian matrix / res vector
solution_vector_Vd(:,1)=solution_vector_DD(:,1);
solution_vector_Vd_saved(:,1)=solution_vector_Vd(:,1);
count_Vd_current=1;
for bias_Vd=1:11 %%% Vd는 1V까지 올림.
    Vd=0.1*bias_Vd-0.1;

    if bias_Vd>=2
        clearvars('-except', 'Iteration_Vd','bias_Vd','Vd')

        Vd_load=0.1*(bias_Vd-1)-0.1;
        FILENAME = sprintf('./data/Homework16_Vd/Homework16_Vd_%.2fV.mat',Vd_load);
        load(FILENAME)
        Vd=Vd_load+0.1;
        clearvars Vd_load
    end
    for Newton_Vd=1:Iteration_Vd
        fprintf("Ramping Drain Voltage, Vd=%.2f, Newton_Vd=%d\n" , Vd, Newton_Vd)
        clearvars Jaco_Vd; clearvars res_Vd;
        %         Jaco_Vd=sparse(zeros(size(solution_vector_Vd,1),size(solution_vector_Vd,1))); res_Vd=zeros(size(solution_vector_Vd,1),1);
        Jaco_Vd=zeros(size(solution_vector_Vd,1),size(solution_vector_Vd,1)); res_Vd=zeros(size(solution_vector_Vd,1),1); %% 희소행렬 아님. Jaco matrix 확인용.

        % Jaco/res %
        for ii=1:length(Element)
            % Control_Volume 계산
            Control_Volume=zeros(3,1);
            Control_Volume(1,1)=(0.25)*L(ii,3)*edge(ii,3)+(0.25)*L(ii,1)*edge(ii,1);
            Control_Volume(2,1)=(0.25)*L(ii,1)*edge(ii,1)+(0.25)*L(ii,2)*edge(ii,2);
            Control_Volume(3,1)=(0.25)*L(ii,2)*edge(ii,2)+(0.25)*L(ii,3)*edge(ii,3);

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
            x12=solution_vector_Vd(index(1,1),1)-solution_vector_Vd(index(4,1),1); x21=-x12;
            x23=solution_vector_Vd(index(4,1),1)-solution_vector_Vd(index(7,1),1); x32=-x23;
            x13=solution_vector_Vd(index(1,1),1)-solution_vector_Vd(index(7,1),1); x31=-x13;
            n1=solution_vector_Vd(index(1,1)+1,1); n2=solution_vector_Vd(index(4,1)+1,1); n3=solution_vector_Vd(index(7,1)+1,1);
            p1=solution_vector_Vd(index(1,1)+2,1); p2=solution_vector_Vd(index(4,1)+2,1); p3=solution_vector_Vd(index(7,1)+2,1);

                % Jaco %
                Jaco_tmp_Vd=zeros(9,9);
                Jaco_tmp_Vd(1,1)=eps(Table_element_region(ii,1),1)*(-edge(ii,1)/L(ii,1)-edge(ii,3)/L(ii,3));
                Jaco_tmp_Vd(1,2)=-coeff*Control_Volume(1,1);
                Jaco_tmp_Vd(1,3)=coeff*Control_Volume(1,1);
                Jaco_tmp_Vd(1,4)=eps(Table_element_region(ii,1),1)*edge(ii,1)/L(ii,1);
                Jaco_tmp_Vd(1,7)=eps(Table_element_region(ii,1),1)*edge(ii,3)/L(ii,3);

                Jaco_tmp_Vd(2,1)=coeff_Jn*((edge(ii,1)/L(ii,1))*(n2/(-V_T)*Ber_d(x21/V_T)-n1/(V_T)*Ber_d(-x21/V_T))+((edge(ii,3)/L(ii,3))*(n3/(-V_T)*Ber_d(x31/V_T)-n1/(V_T)*Ber_d(-x31/V_T))));
                Jaco_tmp_Vd(2,2)=coeff_Jn*((edge(ii,1)/L(ii,1))*(-Ber(-x21/V_T))+((edge(ii,3)/L(ii,3))*(-Ber(-x31/V_T))));
                Jaco_tmp_Vd(2,4)=coeff_Jn*((edge(ii,1)/L(ii,1))*(n2/(V_T)*Ber_d(x21/V_T)-n1/(-V_T)*Ber_d(-x21/V_T)));
                Jaco_tmp_Vd(2,5)=coeff_Jn*((edge(ii,1)/L(ii,1))*(Ber(x21/V_T)));
                Jaco_tmp_Vd(2,7)=coeff_Jn*((edge(ii,3)/L(ii,3))*(n3/(V_T)*Ber_d(x31/V_T)-n1/(-V_T)*Ber_d(-x31/V_T)));
                Jaco_tmp_Vd(2,8)=coeff_Jn*((edge(ii,3)/L(ii,3))*(Ber(x31/V_T)));

                Jaco_tmp_Vd(3,1)=coeff_Jp*((edge(ii,1)/L(ii,1))*(p2/(V_T)*Ber_d(-x21/V_T)-p1/(-V_T)*Ber_d(x21/V_T))+((edge(ii,3)/L(ii,3))*(p3/(V_T)*Ber_d(-x31/V_T)-p1/(-V_T)*Ber_d(x31/V_T))));
                Jaco_tmp_Vd(3,3)=coeff_Jp*((edge(ii,1)/L(ii,1))*(-Ber(x21/V_T))+(edge(ii,3)/L(ii,3))*(-Ber(x31/V_T)));
                Jaco_tmp_Vd(3,4)=coeff_Jp*((edge(ii,1)/L(ii,1))*(p2/(-V_T)*Ber_d(-x21/V_T)-p1/(V_T)*Ber_d(x21/V_T)));
                Jaco_tmp_Vd(3,6)=coeff_Jp*((edge(ii,1)/L(ii,1))*(Ber(-x21/V_T)));
                Jaco_tmp_Vd(3,7)=coeff_Jp*((edge(ii,3)/L(ii,3))*(p3/(-V_T)*Ber_d(-x31/V_T)-p1/(V_T)*Ber_d(x31/V_T)));
                Jaco_tmp_Vd(3,9)=coeff_Jp*((edge(ii,3)/L(ii,3))*(Ber(-x31/V_T)));

                Jaco_tmp_Vd(4,1)=eps(Table_element_region(ii,1),1)*edge(ii,1)/L(ii,1);
                Jaco_tmp_Vd(4,4)=eps(Table_element_region(ii,1),1)*(-edge(ii,1)/L(ii,1)-edge(ii,2)/L(ii,2));
                Jaco_tmp_Vd(4,5)=-coeff*Control_Volume(2,1);
                Jaco_tmp_Vd(4,6)=coeff*Control_Volume(2,1);
                Jaco_tmp_Vd(4,7)=eps(Table_element_region(ii,1),1)*edge(ii,2)/L(ii,2);

                Jaco_tmp_Vd(5,1)=coeff_Jn*((edge(ii,1)/L(ii,1))*(n1/(V_T)*Ber_d(x12/V_T)-n2/(-V_T)*Ber_d(-x12/V_T)));
                Jaco_tmp_Vd(5,2)=coeff_Jn*((edge(ii,1)/L(ii,1))*(Ber(x12/V_T)));
                Jaco_tmp_Vd(5,4)=coeff_Jn*((edge(ii,2)/L(ii,2))*(n3/(-V_T)*Ber_d(x32/V_T)-n2/(V_T)*Ber_d(-x32/V_T))+((edge(ii,1)/L(ii,1))*(n1/(-V_T)*Ber_d(x12/V_T)-n2/(V_T)*Ber_d(-x12/V_T))));
                Jaco_tmp_Vd(5,5)=coeff_Jn*((edge(ii,2)/L(ii,2))*(-Ber(-x32/V_T))+((edge(ii,1)/L(ii,1))*(-Ber(-x12/V_T))));
                Jaco_tmp_Vd(5,7)=coeff_Jn*((edge(ii,2)/L(ii,2))*(n3/(V_T)*Ber_d(x32/V_T)-n2/(-V_T)*Ber_d(-x32/V_T)));
                Jaco_tmp_Vd(5,8)=coeff_Jn*((edge(ii,2)/L(ii,2))*(Ber(x32/V_T)));

                Jaco_tmp_Vd(6,1)=coeff_Jp*((edge(ii,1)/L(ii,1))*(p1/(-V_T)*Ber_d(-x12/V_T)-p2/(V_T)*Ber_d(x12/V_T)));
                Jaco_tmp_Vd(6,3)=coeff_Jp*((edge(ii,1)/L(ii,1))*(Ber(-x12/V_T)));
                Jaco_tmp_Vd(6,4)=coeff_Jp*((edge(ii,2)/L(ii,2))*(p3/(V_T)*Ber_d(-x32/V_T)-p2/(-V_T)*Ber_d(x32/V_T))+((edge(ii,1)/L(ii,1))*(p1/(V_T)*Ber_d(-x12/V_T)-p2/(-V_T)*Ber_d(x12/V_T))));
                Jaco_tmp_Vd(6,6)=coeff_Jp*((edge(ii,2)/L(ii,2))*(-Ber(x32/V_T))+(edge(ii,1)/L(ii,1))*(-Ber(x12/V_T)));
                Jaco_tmp_Vd(6,7)=coeff_Jp*((edge(ii,2)/L(ii,2))*(p3/(-V_T)*Ber_d(-x32/V_T)-p2/(V_T)*Ber_d(x32/V_T)));
                Jaco_tmp_Vd(6,9)=coeff_Jp*((edge(ii,2)/L(ii,2))*(Ber(-x32/V_T)));

                Jaco_tmp_Vd(7,1)=eps(Table_element_region(ii,1),1)*edge(ii,3)/L(ii,3);
                Jaco_tmp_Vd(7,4)=eps(Table_element_region(ii,1),1)*edge(ii,2)/L(ii,2);
                Jaco_tmp_Vd(7,7)=eps(Table_element_region(ii,1),1)*(-edge(ii,2)/L(ii,2)-edge(ii,3)/L(ii,3));
                Jaco_tmp_Vd(7,8)=-coeff*Control_Volume(3,1);
                Jaco_tmp_Vd(7,9)=coeff*Control_Volume(3,1);

                Jaco_tmp_Vd(8,1)=coeff_Jn*((edge(ii,3)/L(ii,3))*(n1/(V_T)*Ber_d(x13/V_T)-n3/(-V_T)*Ber_d(-x13/V_T)));
                Jaco_tmp_Vd(8,2)=coeff_Jn*((edge(ii,3)/L(ii,3))*(Ber(x13/V_T)));
                Jaco_tmp_Vd(8,4)=coeff_Jn*((edge(ii,2)/L(ii,2))*(n2/(V_T)*Ber_d(x23/V_T)-n3/(-V_T)*Ber_d(-x23/V_T)));
                Jaco_tmp_Vd(8,5)=coeff_Jn*((edge(ii,2)/L(ii,2))*(Ber(x23/V_T)));
                Jaco_tmp_Vd(8,7)=coeff_Jn*((edge(ii,3)/L(ii,3))*(n1/(-V_T)*Ber_d(x13/V_T)-n3/(V_T)*Ber_d(-x13/V_T))+((edge(ii,2)/L(ii,2))*(n2/(-V_T)*Ber_d(x23/V_T)-n3/(V_T)*Ber_d(-x23/V_T))));
                Jaco_tmp_Vd(8,8)=coeff_Jn*((edge(ii,3)/L(ii,3))*(-Ber(-x13/V_T))+((edge(ii,2)/L(ii,2))*(-Ber(-x23/V_T))));

                Jaco_tmp_Vd(9,1)=coeff_Jp*((edge(ii,3)/L(ii,3))*(p1/(-V_T)*Ber_d(-x13/V_T)-p3/(V_T)*Ber_d(x13/V_T)));
                Jaco_tmp_Vd(9,3)=coeff_Jp*((edge(ii,3)/L(ii,3))*(Ber(-x13/V_T)));
                Jaco_tmp_Vd(9,4)=coeff_Jp*((edge(ii,2)/L(ii,2))*(p2/(-V_T)*Ber_d(-x23/V_T)-p3/(V_T)*Ber_d(x23/V_T)));
                Jaco_tmp_Vd(9,6)=coeff_Jp*((edge(ii,2)/L(ii,2))*(Ber(-x23/V_T)));
                Jaco_tmp_Vd(9,7)=coeff_Jp*((edge(ii,3)/L(ii,3))*(p1/(V_T)*Ber_d(-x13/V_T)-p3/(-V_T)*Ber_d(x13/V_T))+((edge(ii,2)/L(ii,2))*(p2/(V_T)*Ber_d(-x23/V_T)-p3/(-V_T)*Ber_d(x23/V_T))));
                Jaco_tmp_Vd(9,9)=coeff_Jp*((edge(ii,3)/L(ii,3))*(-Ber(x13/V_T))+(edge(ii,2)/L(ii,2))*(-Ber(x23/V_T)));


                for j=1:9
                    for k=1:9
                        Jaco_Vd(index(j,1),index(k,1))=Jaco_Vd(index(j,1),index(k,1))+Jaco_tmp_Vd(j,k);
                    end
                end

                % res vector
                res_tmp_Vd=zeros(9,1);
                res_potential_tmp1_Vd(1,1)=eps(Table_element_region(ii,1),1)*(edge(ii,1)/L(ii,1)*(solution_vector_Vd(index(4,1),1)-solution_vector_Vd(index(1,1),1))+edge(ii,3)/L(ii,3)*(solution_vector_Vd(index(7,1),1)-solution_vector_Vd(index(1,1),1)));
                res_potential_tmp1_Vd(2,1)=eps(Table_element_region(ii,1),1)*(edge(ii,2)/L(ii,2)*(solution_vector_Vd(index(7,1),1)-solution_vector_Vd(index(4,1),1))+edge(ii,1)/L(ii,1)*(solution_vector_Vd(index(1,1),1)-solution_vector_Vd(index(4,1),1)));
                res_potential_tmp1_Vd(3,1)=eps(Table_element_region(ii,1),1)*(edge(ii,3)/L(ii,3)*(solution_vector_Vd(index(1,1),1)-solution_vector_Vd(index(7,1),1))+edge(ii,2)/L(ii,2)*(solution_vector_Vd(index(4,1),1)-solution_vector_Vd(index(7,1),1)));
                for j=1:3
                    res_potential_tmp2_Vd(j,1)=(N_dop(Table_element_region(ii,1),1)-solution_vector_Vd(index(3*j-1,1),1)+solution_vector_Vd(index(3*j,1),1))*coeff*Control_Volume(j,1);
                end

                res_potential_Vd=res_potential_tmp1_Vd+res_potential_tmp2_Vd;

                % Jn
                res_Jn=zeros(3,1);
                res_Jn(1,1)=coeff_Jn*((edge(ii,1)/L(ii,1))*(n2*Ber(x21/V_T)-n1*Ber(-x21/V_T))+(edge(ii,3)/L(ii,3))*(n3*Ber(x31/V_T)-n1*Ber(-x31/V_T)));
                res_Jn(2,1)=coeff_Jn*((edge(ii,2)/L(ii,2))*(n3*Ber(x32/V_T)-n2*Ber(-x32/V_T))+(edge(ii,1)/L(ii,1))*(n1*Ber(x12/V_T)-n2*Ber(-x12/V_T)));
                res_Jn(3,1)=coeff_Jn*((edge(ii,3)/L(ii,3))*(n1*Ber(x13/V_T)-n3*Ber(-x13/V_T))+(edge(ii,2)/L(ii,2))*(n2*Ber(x23/V_T)-n3*Ber(-x23/V_T)));

                % Jp
                res_Jp=zeros(3,1);
                res_Jp(1,1)=coeff_Jp*((edge(ii,1)/L(ii,1))*(p2*Ber(-x21/V_T)-p1*Ber(x21/V_T))+(edge(ii,3)/L(ii,3))*(p3*Ber(-x31/V_T)-p1*Ber(x31/V_T)));
                res_Jp(2,1)=coeff_Jp*((edge(ii,2)/L(ii,2))*(p3*Ber(-x32/V_T)-p2*Ber(x32/V_T))+(edge(ii,1)/L(ii,1))*(p1*Ber(-x12/V_T)-p2*Ber(x12/V_T)));
                res_Jp(3,1)=coeff_Jp*((edge(ii,3)/L(ii,3))*(p1*Ber(-x13/V_T)-p3*Ber(x13/V_T))+(edge(ii,2)/L(ii,2))*(p2*Ber(-x23/V_T)-p3*Ber(x23/V_T)));

                n=1;
                for j=1:3
                    res_tmp_Vd(n,1)=res_potential_Vd(j,1);
                    res_tmp_Vd(n+1,1)=res_Jn(j,1);
                    res_tmp_Vd(n+2,1)=res_Jp(j,1);
                    n=n+3;
                end

                for j=1:9
                    res_Vd(index(j,1),1)=res_Vd(index(j,1),1)+res_tmp_Vd(j,1);
                end
        end

        % Boundary condition %
        % cathode_BC
        for ii=1:size(anode_BC,1)
            BC_index=find(Table_Jaco(:,1)==1 & Table_Jaco(:,2)==anode_BC(ii,1) & Table_Jaco(:,3)==1);
            Jaco_Vd(BC_index,:)=0; Jaco_Vd(BC_index+1,:)=0; Jaco_Vd(BC_index+2,:)=0;
            Jaco_Vd(BC_index,BC_index)=1; Jaco_Vd(BC_index+1,BC_index+1)=1; Jaco_Vd(BC_index+2,BC_index+2)=1;
            res_Vd(BC_index,1)=solution_vector_Vd(BC_index,1)-V_T*log(Nd/n_int);
            res_Vd(BC_index+1,1)=solution_vector_Vd(BC_index+1,1)-Nd;
            res_Vd(BC_index+2,1)=solution_vector_Vd(BC_index+2,1)-n_int^2/Nd;
        end

        % cathode_BC
        for ii=1:size(cathode_BC,1)
            BC_index=find(Table_Jaco(:,1)==1 & Table_Jaco(:,2)==cathode_BC(ii,1) & Table_Jaco(:,3)==1);
            Jaco_Vd(BC_index,:)=0; Jaco_Vd(BC_index+1,:)=0; Jaco_Vd(BC_index+2,:)=0;
            Jaco_Vd(BC_index,BC_index)=1; Jaco_Vd(BC_index+1,BC_index+1)=1; Jaco_Vd(BC_index+2,BC_index+2)=1;
            res_Vd(BC_index,1)=solution_vector_Vd(BC_index,1)-V_T*log(Nd/n_int)-Vd;
            res_Vd(BC_index+1,1)=solution_vector_Vd(BC_index+1,1)-Nd;
            res_Vd(BC_index+2,1)=solution_vector_Vd(BC_index+2,1)-n_int^2/Nd;
        end

        % Scaling %
        Cvector_Vd=zeros(size(solution_vector_Vd,1),1); v=[V_T; max(abs(N_dop)); max(abs(N_dop))];
        for ii=1:size(solution_vector_Vd,1)
            Cvector_Vd(ii,1)=v(Table_Jaco(ii,3),1);
        end
        Jaco_Vd=sparse(Jaco_Vd);
        Cmatrix_Vd=spdiags(Cvector_Vd,0,size(solution_vector_Vd,1),size(solution_vector_Vd,1));
        Jaco_scaled_Vd=Jaco_Vd*Cmatrix_Vd;
        Rvector_Vd=1./sum(abs(Jaco_scaled_Vd),2);
        Rmatrix_Vd=spdiags(Rvector_Vd,0,size(solution_vector_Vd,1),size(solution_vector_Vd,1));
        Jaco_scaled_Vd=Rmatrix_Vd*Jaco_scaled_Vd;
        res_scaled_Vd=Rmatrix_Vd*res_Vd;
        update_scaled_Vd=Jaco_scaled_Vd\(-res_scaled_Vd);
        update_Vd(:,Newton_Vd)=Cmatrix_Vd*update_scaled_Vd;

        solution_vector_Vd(:,1)=solution_vector_Vd(:,1)+update_Vd(:,Newton_Vd);
        solution_vector_Vd_saved(:,Newton_Vd+1)=solution_vector_Vd(:,1);

        % break
        update_Vd_phi_for_break=update_Vd(find(Table_Jaco(:,1)==1, 1 , 'first'):3:find(Table_Jaco(:,1)==1, 1 , 'last'),Newton_Vd);
        fprintf("Max Updatevector=%d\n\n" , max(abs(update_Vd_phi_for_break)))
        if max(abs(update_Vd_phi_for_break))<1e-15
            break;
        end

    end
    update_Vd_phi=update_Vd(find(Table_Jaco(:,1)==1, 1 , 'first'):3:find(Table_Jaco(:,1)==1, 1 , 'last'),:);

    % Current %
    I_n_Vd(count_Vd_current,1)=0; I_p_Vd(count_Vd_current,1)=0;
    for ii=1:size(cathode_BC,1)
        vertex_current=cathode_BC(ii,1);
        [row,col]=find(Element_si==vertex_current);
        for k=1:size(row,1)

            index=zeros(9,1);
            n=1;
            for j=1:3
                index(n,1)=find(Table_Jaco(:,1)==1 & Table_Jaco(:,2)==Element_si(row(k,1),j) & Table_Jaco(:,3)==1);
                index(n+1,1)=index(n,1)+1;
                index(n+2,1)=index(n,1)+2;
                n=n+3;
            end

            x12=solution_vector_Vd(index(1,1),1)-solution_vector_Vd(index(4,1),1); x21=-x12;
            x23=solution_vector_Vd(index(4,1),1)-solution_vector_Vd(index(7,1),1); x32=-x23;
            x13=solution_vector_Vd(index(1,1),1)-solution_vector_Vd(index(7,1),1); x31=-x13;
            n1=solution_vector_Vd(index(1,1)+1,1); n2=solution_vector_Vd(index(4,1)+1,1); n3=solution_vector_Vd(index(7,1)+1,1);
            p1=solution_vector_Vd(index(1,1)+2,1); p2=solution_vector_Vd(index(4,1)+2,1); p3=solution_vector_Vd(index(7,1)+2,1);


            index_element=find(Table_element_region(:,1)==1 & Table_element_region(:,2)==row(k,1));
            if col==1
                I_n_Vd_tmp=coeff_Jn*((edge(index_element,1)/L(index_element,1))*(n2*Ber(x21/V_T)-n1*Ber(-x21/V_T))+(edge(index_element,3)/L(index_element,3))*(n3*Ber(x31/V_T)-n1*Ber(-x31/V_T)));
                I_p_Vd_tmp=coeff_Jp*((edge(index_element,1)/L(index_element,1))*(p2*Ber(-x21/V_T)-p1*Ber(x21/V_T))+(edge(index_element,3)/L(index_element,3))*(p3*Ber(-x31/V_T)-p1*Ber(x31/V_T)));
            elseif col==2
                I_n_Vd_tmp=coeff_Jn*((edge(index_element,2)/L(index_element,2))*(n3*Ber(x32/V_T)-n2*Ber(-x32/V_T))+(edge(index_element,1)/L(index_element,1))*(n1*Ber(x12/V_T)-n2*Ber(-x12/V_T)));
                I_p_Vd_tmp=coeff_Jp*((edge(index_element,2)/L(index_element,2))*(p3*Ber(-x32/V_T)-p2*Ber(x32/V_T))+(edge(index_element,1)/L(index_element,1))*(p1*Ber(-x12/V_T)-p2*Ber(x12/V_T)));
            elseif col==3
                I_n_Vd_tmp=coeff_Jn*((edge(index_element,3)/L(index_element,3))*(n1*Ber(x13/V_T)-n3*Ber(-x13/V_T))+(edge(index_element,2)/L(index_element,2))*(n2*Ber(x23/V_T)-n3*Ber(-x23/V_T)));
                I_p_Vd_tmp=coeff_Jp*((edge(index_element,3)/L(index_element,3))*(p1*Ber(-x13/V_T)-p3*Ber(x13/V_T))+(edge(index_element,2)/L(index_element,2))*(p2*Ber(-x23/V_T)-p3*Ber(x23/V_T)));
            end
            I_n_Vd(count_Vd_current,1)=I_n_Vd(count_Vd_current,1)+I_n_Vd_tmp;
            I_p_Vd(count_Vd_current,1)=I_p_Vd(count_Vd_current,1)+I_p_Vd_tmp;

        end
    end
    I_Vd(count_Vd_current,1) = (I_n_Vd(count_Vd_current,1)+I_p_Vd(count_Vd_current,1))*1e-6;
    count_Vd_current=count_Vd_current+1;

    % Save
    FILENAME = sprintf('./data/Homework16_Vd/Homework16_Vd_%.2fV.mat' , Vd);
    save(FILENAME);

    FILENAME = sprintf('./data/Homework16_Vd/Homework16_Vd_%.2fV_current.mat',Vd);
    save(FILENAME, 'I_Vd');

end

%% Visualize - Vd
Visual_solution_vector_Vd=zeros(size(Vertex,1),3);
for ii=1:size(Table_Jaco,1)
    Visual_solution_vector_Vd(Table_Jaco(ii,2),Table_Jaco(ii,3))=solution_vector_Vd_saved(ii,Newton_Vd);
end

figure(1) % 점으로 3D
plot3(Vertex(:,1),Vertex(:,2),Visual_solution_vector_Vd(:,1),'*')