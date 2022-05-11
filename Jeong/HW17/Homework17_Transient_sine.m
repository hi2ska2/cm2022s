clear
load('matlab.mat')

%% Transient
Iteration_Vg=40;
f=1E14;
period=2;
time_step=100;
Amp=1e-3;

T=1/f;
delta_t=T/time_step;
time=transpose(0:delta_t:period*T);

% Jacobian matrix / res vector
solution_vector_Vg(:,1)=solution_vector_Vd(:,1);
solution_vector_Vg_saved(:,1)=solution_vector_Vg(:,1);

count_Vg_current=1;
for time_for=1:length(time) %%% Vg는 1V까지 올림.
    Vg(time_for,1)=Amp*sin(2*pi*f*time(time_for));
    solution_vector_Vg_old=solution_vector_Vg;

    for Newton_Vg=1:Iteration_Vg
        fprintf("Ramping Drain Voltage, Vg=%.6f, Newton_Vg=%d\n" , Vg(time_for,1), Newton_Vg)
        %         Jaco_Vg=sparse(zeros(size(solution_vector_Vg,1),size(solution_vector_Vg,1))); res_Vg=zeros(size(solution_vector_Vg,1),1);
        Jaco_Vg=zeros(size(solution_vector_Vg,1),size(solution_vector_Vg,1)); res_Vg=zeros(size(solution_vector_Vg,1),1); %% 희소행렬 아님. Jaco matrix 확인용.

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
            x12=solution_vector_Vg(index(1,1),1)-solution_vector_Vg(index(4,1),1); x21=-x12;
            x23=solution_vector_Vg(index(4,1),1)-solution_vector_Vg(index(7,1),1); x32=-x23;
            x13=solution_vector_Vg(index(1,1),1)-solution_vector_Vg(index(7,1),1); x31=-x13;
            n1=solution_vector_Vg(index(1,1)+1,1); n2=solution_vector_Vg(index(4,1)+1,1); n3=solution_vector_Vg(index(7,1)+1,1);
            p1=solution_vector_Vg(index(1,1)+2,1); p2=solution_vector_Vg(index(4,1)+2,1); p3=solution_vector_Vg(index(7,1)+2,1);
            % old solution vector
            n1_old=solution_vector_Vg_old(index(1,1)+1,1); n2_old=solution_vector_Vg_old(index(4,1)+1,1); n3_old=solution_vector_Vg_old(index(7,1)+1,1);
            p1_old=solution_vector_Vg_old(index(1,1)+2,1); p2_old=solution_vector_Vg_old(index(4,1)+2,1); p3_old=solution_vector_Vg_old(index(7,1)+2,1);

                % Jaco %
                Jaco_tmp_Vg=zeros(9,9);
                Jaco_tmp_Vg(1,1)=eps(Table_element_region(ii,1),1)*(-edge(ii,1)/L(ii,1)-edge(ii,3)/L(ii,3));
                Jaco_tmp_Vg(1,2)=-coeff*Control_Volume(1,1);
                Jaco_tmp_Vg(1,3)=coeff*Control_Volume(1,1);
                Jaco_tmp_Vg(1,4)=eps(Table_element_region(ii,1),1)*edge(ii,1)/L(ii,1);
                Jaco_tmp_Vg(1,7)=eps(Table_element_region(ii,1),1)*edge(ii,3)/L(ii,3);

                Jaco_tmp_Vg(2,1)=coeff_Jn*((edge(ii,1)/L(ii,1))*(n2/(-V_T)*Ber_d(x21/V_T)-n1/(V_T)*Ber_d(-x21/V_T))+((edge(ii,3)/L(ii,3))*(n3/(-V_T)*Ber_d(x31/V_T)-n1/(V_T)*Ber_d(-x31/V_T))));
                Jaco_tmp_Vg(2,2)=coeff_Jn*((edge(ii,1)/L(ii,1))*(-Ber(-x21/V_T))+((edge(ii,3)/L(ii,3))*(-Ber(-x31/V_T))));
                Jaco_tmp_Vg(2,4)=coeff_Jn*((edge(ii,1)/L(ii,1))*(n2/(V_T)*Ber_d(x21/V_T)-n1/(-V_T)*Ber_d(-x21/V_T)));
                Jaco_tmp_Vg(2,5)=coeff_Jn*((edge(ii,1)/L(ii,1))*(Ber(x21/V_T)));
                Jaco_tmp_Vg(2,7)=coeff_Jn*((edge(ii,3)/L(ii,3))*(n3/(V_T)*Ber_d(x31/V_T)-n1/(-V_T)*Ber_d(-x31/V_T)));
                Jaco_tmp_Vg(2,8)=coeff_Jn*((edge(ii,3)/L(ii,3))*(Ber(x31/V_T)));
                % electron_1_DD_transient
                Jaco_tmp_Vg(2,2)=Jaco_tmp_Vg(2,2)-q*Control_Volume(1,1)/delta_t;

                Jaco_tmp_Vg(3,1)=coeff_Jp*((edge(ii,1)/L(ii,1))*(p2/(V_T)*Ber_d(-x21/V_T)-p1/(-V_T)*Ber_d(x21/V_T))+((edge(ii,3)/L(ii,3))*(p3/(V_T)*Ber_d(-x31/V_T)-p1/(-V_T)*Ber_d(x31/V_T))));
                Jaco_tmp_Vg(3,3)=coeff_Jp*((edge(ii,1)/L(ii,1))*(-Ber(x21/V_T))+(edge(ii,3)/L(ii,3))*(-Ber(x31/V_T)));
                Jaco_tmp_Vg(3,4)=coeff_Jp*((edge(ii,1)/L(ii,1))*(p2/(-V_T)*Ber_d(-x21/V_T)-p1/(V_T)*Ber_d(x21/V_T)));
                Jaco_tmp_Vg(3,6)=coeff_Jp*((edge(ii,1)/L(ii,1))*(Ber(-x21/V_T)));
                Jaco_tmp_Vg(3,7)=coeff_Jp*((edge(ii,3)/L(ii,3))*(p3/(-V_T)*Ber_d(-x31/V_T)-p1/(V_T)*Ber_d(x31/V_T)));
                Jaco_tmp_Vg(3,9)=coeff_Jp*((edge(ii,3)/L(ii,3))*(Ber(-x31/V_T)));
                % hole_1_transient
                Jaco_tmp_Vg(3,3)=Jaco_tmp_Vg(3,3)+q*Control_Volume(1,1)/delta_t;

                Jaco_tmp_Vg(4,1)=eps(Table_element_region(ii,1),1)*edge(ii,1)/L(ii,1);
                Jaco_tmp_Vg(4,4)=eps(Table_element_region(ii,1),1)*(-edge(ii,1)/L(ii,1)-edge(ii,2)/L(ii,2));
                Jaco_tmp_Vg(4,5)=-coeff*Control_Volume(2,1);
                Jaco_tmp_Vg(4,6)=coeff*Control_Volume(2,1);
                Jaco_tmp_Vg(4,7)=eps(Table_element_region(ii,1),1)*edge(ii,2)/L(ii,2);

                Jaco_tmp_Vg(5,1)=coeff_Jn*((edge(ii,1)/L(ii,1))*(n1/(V_T)*Ber_d(x12/V_T)-n2/(-V_T)*Ber_d(-x12/V_T)));
                Jaco_tmp_Vg(5,2)=coeff_Jn*((edge(ii,1)/L(ii,1))*(Ber(x12/V_T)));
                Jaco_tmp_Vg(5,4)=coeff_Jn*((edge(ii,2)/L(ii,2))*(n3/(-V_T)*Ber_d(x32/V_T)-n2/(V_T)*Ber_d(-x32/V_T))+((edge(ii,1)/L(ii,1))*(n1/(-V_T)*Ber_d(x12/V_T)-n2/(V_T)*Ber_d(-x12/V_T))));
                Jaco_tmp_Vg(5,5)=coeff_Jn*((edge(ii,2)/L(ii,2))*(-Ber(-x32/V_T))+((edge(ii,1)/L(ii,1))*(-Ber(-x12/V_T))));
                Jaco_tmp_Vg(5,7)=coeff_Jn*((edge(ii,2)/L(ii,2))*(n3/(V_T)*Ber_d(x32/V_T)-n2/(-V_T)*Ber_d(-x32/V_T)));
                Jaco_tmp_Vg(5,8)=coeff_Jn*((edge(ii,2)/L(ii,2))*(Ber(x32/V_T)));
                % electron_2_DD_transient
                Jaco_tmp_Vg(5,5)=Jaco_tmp_Vg(5,5)-q*Control_Volume(2,1)/delta_t;

                Jaco_tmp_Vg(6,1)=coeff_Jp*((edge(ii,1)/L(ii,1))*(p1/(-V_T)*Ber_d(-x12/V_T)-p2/(V_T)*Ber_d(x12/V_T)));
                Jaco_tmp_Vg(6,3)=coeff_Jp*((edge(ii,1)/L(ii,1))*(Ber(-x12/V_T)));
                Jaco_tmp_Vg(6,4)=coeff_Jp*((edge(ii,2)/L(ii,2))*(p3/(V_T)*Ber_d(-x32/V_T)-p2/(-V_T)*Ber_d(x32/V_T))+((edge(ii,1)/L(ii,1))*(p1/(V_T)*Ber_d(-x12/V_T)-p2/(-V_T)*Ber_d(x12/V_T))));
                Jaco_tmp_Vg(6,6)=coeff_Jp*((edge(ii,2)/L(ii,2))*(-Ber(x32/V_T))+(edge(ii,1)/L(ii,1))*(-Ber(x12/V_T)));
                Jaco_tmp_Vg(6,7)=coeff_Jp*((edge(ii,2)/L(ii,2))*(p3/(-V_T)*Ber_d(-x32/V_T)-p2/(V_T)*Ber_d(x32/V_T)));
                Jaco_tmp_Vg(6,9)=coeff_Jp*((edge(ii,2)/L(ii,2))*(Ber(-x32/V_T)));
                % hole_2_transient
                Jaco_tmp_Vg(6,6)=Jaco_tmp_Vg(6,6)+q*Control_Volume(2,1)/delta_t;

                Jaco_tmp_Vg(7,1)=eps(Table_element_region(ii,1),1)*edge(ii,3)/L(ii,3);
                Jaco_tmp_Vg(7,4)=eps(Table_element_region(ii,1),1)*edge(ii,2)/L(ii,2);
                Jaco_tmp_Vg(7,7)=eps(Table_element_region(ii,1),1)*(-edge(ii,2)/L(ii,2)-edge(ii,3)/L(ii,3));
                Jaco_tmp_Vg(7,8)=-coeff*Control_Volume(3,1);
                Jaco_tmp_Vg(7,9)=coeff*Control_Volume(3,1);

                Jaco_tmp_Vg(8,1)=coeff_Jn*((edge(ii,3)/L(ii,3))*(n1/(V_T)*Ber_d(x13/V_T)-n3/(-V_T)*Ber_d(-x13/V_T)));
                Jaco_tmp_Vg(8,2)=coeff_Jn*((edge(ii,3)/L(ii,3))*(Ber(x13/V_T)));
                Jaco_tmp_Vg(8,4)=coeff_Jn*((edge(ii,2)/L(ii,2))*(n2/(V_T)*Ber_d(x23/V_T)-n3/(-V_T)*Ber_d(-x23/V_T)));
                Jaco_tmp_Vg(8,5)=coeff_Jn*((edge(ii,2)/L(ii,2))*(Ber(x23/V_T)));
                Jaco_tmp_Vg(8,7)=coeff_Jn*((edge(ii,3)/L(ii,3))*(n1/(-V_T)*Ber_d(x13/V_T)-n3/(V_T)*Ber_d(-x13/V_T))+((edge(ii,2)/L(ii,2))*(n2/(-V_T)*Ber_d(x23/V_T)-n3/(V_T)*Ber_d(-x23/V_T))));
                Jaco_tmp_Vg(8,8)=coeff_Jn*((edge(ii,3)/L(ii,3))*(-Ber(-x13/V_T))+((edge(ii,2)/L(ii,2))*(-Ber(-x23/V_T))));
                % electron_1_DD_transient
                Jaco_tmp_Vg(8,8)=Jaco_tmp_Vg(8,8)-q*Control_Volume(3,1)/delta_t;

                Jaco_tmp_Vg(9,1)=coeff_Jp*((edge(ii,3)/L(ii,3))*(p1/(-V_T)*Ber_d(-x13/V_T)-p3/(V_T)*Ber_d(x13/V_T)));
                Jaco_tmp_Vg(9,3)=coeff_Jp*((edge(ii,3)/L(ii,3))*(Ber(-x13/V_T)));
                Jaco_tmp_Vg(9,4)=coeff_Jp*((edge(ii,2)/L(ii,2))*(p2/(-V_T)*Ber_d(-x23/V_T)-p3/(V_T)*Ber_d(x23/V_T)));
                Jaco_tmp_Vg(9,6)=coeff_Jp*((edge(ii,2)/L(ii,2))*(Ber(-x23/V_T)));
                Jaco_tmp_Vg(9,7)=coeff_Jp*((edge(ii,3)/L(ii,3))*(p1/(V_T)*Ber_d(-x13/V_T)-p3/(-V_T)*Ber_d(x13/V_T))+((edge(ii,2)/L(ii,2))*(p2/(V_T)*Ber_d(-x23/V_T)-p3/(-V_T)*Ber_d(x23/V_T))));
                Jaco_tmp_Vg(9,9)=coeff_Jp*((edge(ii,3)/L(ii,3))*(-Ber(x13/V_T))+(edge(ii,2)/L(ii,2))*(-Ber(x23/V_T)));
                % hole_3_transient
                Jaco_tmp_Vg(9,9)=Jaco_tmp_Vg(9,9)+q*Control_Volume(3,1)/delta_t;


                for j=1:9
                    for k=1:9
                        Jaco_Vg(index(j,1),index(k,1))=Jaco_Vg(index(j,1),index(k,1))+Jaco_tmp_Vg(j,k);
                    end
                end

                % res vector
                res_tmp_Vg=zeros(9,1);
                res_potential_tmp1_Vg(1,1)=eps(Table_element_region(ii,1),1)*(edge(ii,1)/L(ii,1)*(solution_vector_Vg(index(4,1),1)-solution_vector_Vg(index(1,1),1))+edge(ii,3)/L(ii,3)*(solution_vector_Vg(index(7,1),1)-solution_vector_Vg(index(1,1),1)));
                res_potential_tmp1_Vg(2,1)=eps(Table_element_region(ii,1),1)*(edge(ii,2)/L(ii,2)*(solution_vector_Vg(index(7,1),1)-solution_vector_Vg(index(4,1),1))+edge(ii,1)/L(ii,1)*(solution_vector_Vg(index(1,1),1)-solution_vector_Vg(index(4,1),1)));
                res_potential_tmp1_Vg(3,1)=eps(Table_element_region(ii,1),1)*(edge(ii,3)/L(ii,3)*(solution_vector_Vg(index(1,1),1)-solution_vector_Vg(index(7,1),1))+edge(ii,2)/L(ii,2)*(solution_vector_Vg(index(4,1),1)-solution_vector_Vg(index(7,1),1)));
                for j=1:3
                    res_potential_tmp2_Vg(j,1)=(N_dop(Table_element_region(ii,1),1)-solution_vector_Vg(index(3*j-1,1),1)+solution_vector_Vg(index(3*j,1),1))*coeff*Control_Volume(j,1);
                end

                res_potential_Vg=res_potential_tmp1_Vg+res_potential_tmp2_Vg;

                % Jn
                res_Jn=zeros(3,1);
                res_Jn(1,1)=coeff_Jn*((edge(ii,1)/L(ii,1))*(n2*Ber(x21/V_T)-n1*Ber(-x21/V_T))+(edge(ii,3)/L(ii,3))*(n3*Ber(x31/V_T)-n1*Ber(-x31/V_T)));
                res_Jn(2,1)=coeff_Jn*((edge(ii,2)/L(ii,2))*(n3*Ber(x32/V_T)-n2*Ber(-x32/V_T))+(edge(ii,1)/L(ii,1))*(n1*Ber(x12/V_T)-n2*Ber(-x12/V_T)));
                res_Jn(3,1)=coeff_Jn*((edge(ii,3)/L(ii,3))*(n1*Ber(x13/V_T)-n3*Ber(-x13/V_T))+(edge(ii,2)/L(ii,2))*(n2*Ber(x23/V_T)-n3*Ber(-x23/V_T)));
                % transient
                res_Jn(1,1)=res_Jn(1,1)-q*Control_Volume(1,1)*(n1-n1_old)/delta_t;
                res_Jn(2,1)=res_Jn(2,1)-q*Control_Volume(2,1)*(n2-n2_old)/delta_t;
                res_Jn(3,1)=res_Jn(3,1)-q*Control_Volume(3,1)*(n3-n3_old)/delta_t;

                % Jp
                res_Jp=zeros(3,1);
                res_Jp(1,1)=coeff_Jp*((edge(ii,1)/L(ii,1))*(p2*Ber(-x21/V_T)-p1*Ber(x21/V_T))+(edge(ii,3)/L(ii,3))*(p3*Ber(-x31/V_T)-p1*Ber(x31/V_T)));
                res_Jp(2,1)=coeff_Jp*((edge(ii,2)/L(ii,2))*(p3*Ber(-x32/V_T)-p2*Ber(x32/V_T))+(edge(ii,1)/L(ii,1))*(p1*Ber(-x12/V_T)-p2*Ber(x12/V_T)));
                res_Jp(3,1)=coeff_Jp*((edge(ii,3)/L(ii,3))*(p1*Ber(-x13/V_T)-p3*Ber(x13/V_T))+(edge(ii,2)/L(ii,2))*(p2*Ber(-x23/V_T)-p3*Ber(x23/V_T)));
                % transient
                res_Jp(1,1)=res_Jp(1,1)+q*Control_Volume(1,1)*(p1-p1_old)/delta_t;
                res_Jp(2,1)=res_Jp(2,1)+q*Control_Volume(2,1)*(p2-p2_old)/delta_t;
                res_Jp(3,1)=res_Jp(3,1)+q*Control_Volume(3,1)*(p3-p3_old)/delta_t;


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
        end

        % Boundary condition %
        % cathode_BC
        for ii=1:size(anode_BC,1)
            BC_index=find(Table_Jaco(:,1)==1 & Table_Jaco(:,2)==anode_BC(ii,1) & Table_Jaco(:,3)==1);
            Jaco_Vg(BC_index,:)=0; Jaco_Vg(BC_index+1,:)=0; Jaco_Vg(BC_index+2,:)=0;
            Jaco_Vg(BC_index,BC_index)=1; Jaco_Vg(BC_index+1,BC_index+1)=1; Jaco_Vg(BC_index+2,BC_index+2)=1;
            res_Vg(BC_index,1)=solution_vector_Vg(BC_index,1)-V_T*log(Nd/n_int);
            res_Vg(BC_index+1,1)=solution_vector_Vg(BC_index+1,1)-Nd;
            res_Vg(BC_index+2,1)=solution_vector_Vg(BC_index+2,1)-n_int^2/Nd;
        end

        % cathode_BC
        for ii=1:size(cathode_BC,1)
            BC_index=find(Table_Jaco(:,1)==1 & Table_Jaco(:,2)==cathode_BC(ii,1) & Table_Jaco(:,3)==1);
            Jaco_Vg(BC_index,:)=0; Jaco_Vg(BC_index+1,:)=0; Jaco_Vg(BC_index+2,:)=0;
            Jaco_Vg(BC_index,BC_index)=1; Jaco_Vg(BC_index+1,BC_index+1)=1; Jaco_Vg(BC_index+2,BC_index+2)=1;
            res_Vg(BC_index,1)=solution_vector_Vg(BC_index,1)-V_T*log(Nd/n_int)-Vd-Vg(time_for,1);
            res_Vg(BC_index+1,1)=solution_vector_Vg(BC_index+1,1)-Nd;
            res_Vg(BC_index+2,1)=solution_vector_Vg(BC_index+2,1)-n_int^2/Nd;
        end

        % Scaling %
        Cvector_Vg=zeros(size(solution_vector_Vg,1),1); v=[V_T; max(abs(N_dop)); max(abs(N_dop))];
        for ii=1:size(solution_vector_Vg,1)
            Cvector_Vg(ii,1)=v(Table_Jaco(ii,3),1);
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
        update_Vg_phi_for_break=update_Vg(find(Table_Jaco(:,1)==1, 1 , 'first'):3:find(Table_Jaco(:,1)==1, 1 , 'last'),Newton_Vg);
        update_Vg_elec_for_break=update_Vg(find(Table_Jaco(:,1)==1, 1 , 'first')+1:3:find(Table_Jaco(:,1)==1, 1 , 'last')+1,Newton_Vg);
        update_Vg_hole_for_break=update_Vg(find(Table_Jaco(:,1)==1, 1 , 'first')+2:3:find(Table_Jaco(:,1)==1, 1 , 'last')+2,Newton_Vg);
        fprintf("Max Updatevector=%d\n\n" , max(abs(update_Vg_phi_for_break)))
        if max(abs(update_Vg_phi_for_break))<1e-15 && max(abs(update_Vg_elec_for_break))<1e10 && max(abs(update_Vg_hole_for_break))<1e-2
            break;
        end

    end
    update_Vg_phi=update_Vg(find(Table_Jaco(:,1)==1, 1 , 'first'):3:find(Table_Jaco(:,1)==1, 1 , 'last'),:);

    % Current %
    I_n_Vg(count_Vg_current,1)=0; I_p_Vg(count_Vg_current,1)=0; I_d_Vg(count_Vg_current,1)=0;
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

            x12=solution_vector_Vg(index(1,1),1)-solution_vector_Vg(index(4,1),1); x21=-x12;
            x23=solution_vector_Vg(index(4,1),1)-solution_vector_Vg(index(7,1),1); x32=-x23;
            x13=solution_vector_Vg(index(1,1),1)-solution_vector_Vg(index(7,1),1); x31=-x13;
            n1=solution_vector_Vg(index(1,1)+1,1); n2=solution_vector_Vg(index(4,1)+1,1); n3=solution_vector_Vg(index(7,1)+1,1);
            p1=solution_vector_Vg(index(1,1)+2,1); p2=solution_vector_Vg(index(4,1)+2,1); p3=solution_vector_Vg(index(7,1)+2,1);

            x12_old=solution_vector_Vg_old(index(1,1),1)-solution_vector_Vg_old(index(4,1),1); x21_old=-x12_old;
            x23_old=solution_vector_Vg_old(index(4,1),1)-solution_vector_Vg_old(index(7,1),1); x32_old=-x23_old;
            x13_old=solution_vector_Vg_old(index(1,1),1)-solution_vector_Vg_old(index(7,1),1); x31_old=-x13_old;


            index_element=find(Table_element_region(:,1)==1 & Table_element_region(:,2)==row(k,1));
            eps_si=11.7;
            if col(k,1)==1
                I_n_Vg_tmp=coeff_Jn*((edge(index_element,1)/L(index_element,1))*(n2*Ber(x21/V_T)-n1*Ber(-x21/V_T))+(edge(index_element,3)/L(index_element,3))*(n3*Ber(x31/V_T)-n1*Ber(-x31/V_T)));
                I_p_Vg_tmp=coeff_Jp*((edge(index_element,1)/L(index_element,1))*(p2*Ber(-x21/V_T)-p1*Ber(x21/V_T))+(edge(index_element,3)/L(index_element,3))*(p3*Ber(-x31/V_T)-p1*Ber(x31/V_T)));
                I_d_tmp=-eps_si*eps0*(x21-x21_old)/delta_t/L(index_element,1)*edge(index_element,1)-eps_si*eps0*(x31-x31_old)/delta_t/L(index_element,3)*edge(index_element,3);
            elseif col(k,1)==2
                I_n_Vg_tmp=coeff_Jn*((edge(index_element,2)/L(index_element,2))*(n3*Ber(x32/V_T)-n2*Ber(-x32/V_T))+(edge(index_element,1)/L(index_element,1))*(n1*Ber(x12/V_T)-n2*Ber(-x12/V_T)));
                I_p_Vg_tmp=coeff_Jp*((edge(index_element,2)/L(index_element,2))*(p3*Ber(-x32/V_T)-p2*Ber(x32/V_T))+(edge(index_element,1)/L(index_element,1))*(p1*Ber(-x12/V_T)-p2*Ber(x12/V_T)));
                I_d_tmp=-eps_si*eps0*(x32-x32_old)/delta_t/L(index_element,2)*edge(index_element,2)-eps_si*eps0*(x12-x12_old)/delta_t/L(index_element,1)*edge(index_element,1);
            elseif col(k,1)==3
                I_n_Vg_tmp=coeff_Jn*((edge(index_element,3)/L(index_element,3))*(n1*Ber(x13/V_T)-n3*Ber(-x13/V_T))+(edge(index_element,2)/L(index_element,2))*(n2*Ber(x23/V_T)-n3*Ber(-x23/V_T)));
                I_p_Vg_tmp=coeff_Jp*((edge(index_element,3)/L(index_element,3))*(p1*Ber(-x13/V_T)-p3*Ber(x13/V_T))+(edge(index_element,2)/L(index_element,2))*(p2*Ber(-x23/V_T)-p3*Ber(x23/V_T)));
                I_d_tmp=-eps_si*eps0*(x23-x23_old)/delta_t/L(index_element,2)*edge(index_element,2)-eps_si*eps0*(x13-x13_old)/delta_t/L(index_element,3)*edge(index_element,3);
            end
            I_n_Vg(count_Vg_current,1)=I_n_Vg(count_Vg_current,1)+I_n_Vg_tmp*1e-6;
            I_p_Vg(count_Vg_current,1)=I_p_Vg(count_Vg_current,1)+I_p_Vg_tmp*1e-6;
            I_d_Vg(count_Vg_current,1)=I_d_Vg(count_Vg_current,1)+I_d_tmp*1e-6;


        end
    end
    I_Vg(count_Vg_current,1) = (I_n_Vg(count_Vg_current,1)+I_p_Vg(count_Vg_current,1)+I_d_Vg(count_Vg_current,1));
    count_Vg_current=count_Vg_current+1;

%     % Save
%     FILENAME = sprintf('./data/Homework17_Vg/Homework17_Vg_%.2fV.mat' , Vg);
%     save(FILENAME);
% 
%     FILENAME = sprintf('./data/Homework17_Vg/Homework17_Vg_%.2fV_current.mat',Vg);
%     save(FILENAME, 'I_Vg');

end
max(I_d_Vg)/Amp
%% Visualize - Vg
Visual_solution_vector_Vg=zeros(size(Vertex,1),3);
for ii=1:size(Table_Jaco,1)
    Visual_solution_vector_Vg(Table_Jaco(ii,2),Table_Jaco(ii,3))=solution_vector_Vg_saved(ii,Newton_Vg);
end

% figure(1) % 점으로 3D
% plot3(Vertex(:,1),Vertex(:,2),Visual_solution_vector_Vg(:,1),'*')
% figure(2) % 점으로 3D
% plot3(Vertex(:,1),Vertex(:,2),Visual_solution_vector_Vg(:,2),'*')
% figure(3) % 점으로 3D
% plot3(Vertex(:,1),Vertex(:,2),Visual_solution_vector_Vg(:,3),'*')