clear; clc; close all;


% Load할 Vd 지정
Vd=0.5;
FILENAME = sprintf('./data/Homework16_Vd/Homework16_Vd_%.2fV.mat' , Vd);
load(FILENAME)

% 제작한 pulse function을 활용했다.
% [Voltage_save, time_save] = pulse(T, Time_step, Amp, Duty_cycle, Slope, Period)

T=1e-23; % [ns], 주기
Time_step=100;
Amp=2;
period=2;

DIRNAME = datestr(now,'yyyymmdd_HHMMSS');
Make_DIRNAME = sprintf('./data/Homework16_Transient/%s_T_%.2e_sine' ,DIRNAME, T);
mkdir(Make_DIRNAME);


%% Transient simulation
Iteration_Vg=10000;
delta_t=T/Time_step;
time=0:delta_t:T*period;
% count_save=1;


solution_vector_Vg(:,1)=solution_vector_Vd(:,1);
solution_vector_Vg_old(:,1)=solution_vector_Vd(:,1);
solution_vector_Vg_saved(:,1)=solution_vector_Vg(:,1);

Vg=zeros(length(time),1);
I=zeros(length(time),1);
Newton_Vg_save=zeros(length(time),1);
I_n=zeros(length(time),1); I_p=zeros(length(time),1); I_d=zeros(length(time),1);

% Jacobian matrix / res vector
for TIME_for=1:length(time)
    TIME_save=TIME_for;
    t=time(TIME_save);
    Vg_tmp=Amp*sin(2*pi/T*t);
    Vg(TIME_save,1)=Vg_tmp;

    for Newton_Vg=1:Iteration_Vg
        fprintf("##Transient simulation##\nRamping Gate Voltage, Time=%.2es -> Vg=%.5fV, Newton_Vg=%d\n" , t, Vg_tmp, Newton_Vg)
        clearvars Jaco_Vg; clearvars res_Vg;
        update_Vg=zeros(size(solution_vector_Vg,1),Iteration_Vg);
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

                % 식 간단히 표기하기 위해 미리 계산 %
                % present solution vector
                x12=solution_vector_Vg(index(1,1),1)-solution_vector_Vg(index(4,1),1); x21=-x12;
                x23=solution_vector_Vg(index(4,1),1)-solution_vector_Vg(index(7,1),1); x32=-x23;
                x13=solution_vector_Vg(index(1,1),1)-solution_vector_Vg(index(7,1),1); x31=-x13;
                n1=solution_vector_Vg(index(1,1)+1,1); n2=solution_vector_Vg(index(4,1)+1,1); n3=solution_vector_Vg(index(7,1)+1,1);
                p1=solution_vector_Vg(index(1,1)+2,1); p2=solution_vector_Vg(index(4,1)+2,1); p3=solution_vector_Vg(index(7,1)+2,1);
                % old solution vector
                n1_old=solution_vector_Vg_old(index(1,1)+1,1); n2_old=solution_vector_Vg_old(index(4,1)+1,1); n3_old=solution_vector_Vg_old(index(7,1)+1,1);
                p1_old=solution_vector_Vg_old(index(1,1)+2,1); p2_old=solution_vector_Vg_old(index(4,1)+2,1); p3_old=solution_vector_Vg_old(index(7,1)+2,1);


                % Jaco %
                % potential_1
                Jaco_tmp_Vg=zeros(9,9);
                Jaco_tmp_Vg(1,1)=eps(Table_element_region(ii,1),1)*(-edge(ii,1)/L(ii,1)-edge(ii,3)/L(ii,3));
                Jaco_tmp_Vg(1,2)=-coeff*Control_Volume(1,1);
                Jaco_tmp_Vg(1,3)=coeff*Control_Volume(1,1);
                Jaco_tmp_Vg(1,4)=eps(Table_element_region(ii,1),1)*edge(ii,1)/L(ii,1);
                Jaco_tmp_Vg(1,7)=eps(Table_element_region(ii,1),1)*edge(ii,3)/L(ii,3);
                
                % electron_1_DD
                Jaco_tmp_Vg(2,1)=coeff_Jn*((edge(ii,1)/L(ii,1))*(n2/(-V_T)*Ber_d(x21/V_T)-n1/(V_T)*Ber_d(-x21/V_T))+((edge(ii,3)/L(ii,3))*(n3/(-V_T)*Ber_d(x31/V_T)-n1/(V_T)*Ber_d(-x31/V_T))));
                Jaco_tmp_Vg(2,2)=coeff_Jn*((edge(ii,1)/L(ii,1))*(-Ber(-x21/V_T))+((edge(ii,3)/L(ii,3))*(-Ber(-x31/V_T))));
                Jaco_tmp_Vg(2,4)=coeff_Jn*((edge(ii,1)/L(ii,1))*(n2/(V_T)*Ber_d(x21/V_T)-n1/(-V_T)*Ber_d(-x21/V_T)));
                Jaco_tmp_Vg(2,5)=coeff_Jn*((edge(ii,1)/L(ii,1))*(Ber(x21/V_T)));
                Jaco_tmp_Vg(2,7)=coeff_Jn*((edge(ii,3)/L(ii,3))*(n3/(V_T)*Ber_d(x31/V_T)-n1/(-V_T)*Ber_d(-x31/V_T)));
                Jaco_tmp_Vg(2,8)=coeff_Jn*((edge(ii,3)/L(ii,3))*(Ber(x31/V_T)));
                % electron_1_DD_transient
                Jaco_tmp_Vg(2,2)=Jaco_tmp_Vg(2,2)-q*Control_Volume(1,1)/delta_t;

                % hole_1_DD
                Jaco_tmp_Vg(3,1)=coeff_Jp*((edge(ii,1)/L(ii,1))*(p2/(V_T)*Ber_d(-x21/V_T)-p1/(-V_T)*Ber_d(x21/V_T))+((edge(ii,3)/L(ii,3))*(p3/(V_T)*Ber_d(-x31/V_T)-p1/(-V_T)*Ber_d(x31/V_T))));
                Jaco_tmp_Vg(3,3)=coeff_Jp*((edge(ii,1)/L(ii,1))*(-Ber(x21/V_T))+(edge(ii,3)/L(ii,3))*(-Ber(x31/V_T)));
                Jaco_tmp_Vg(3,4)=coeff_Jp*((edge(ii,1)/L(ii,1))*(p2/(-V_T)*Ber_d(-x21/V_T)-p1/(V_T)*Ber_d(x21/V_T)));
                Jaco_tmp_Vg(3,6)=coeff_Jp*((edge(ii,1)/L(ii,1))*(Ber(-x21/V_T)));
                Jaco_tmp_Vg(3,7)=coeff_Jp*((edge(ii,3)/L(ii,3))*(p3/(-V_T)*Ber_d(-x31/V_T)-p1/(V_T)*Ber_d(x31/V_T)));
                Jaco_tmp_Vg(3,9)=coeff_Jp*((edge(ii,3)/L(ii,3))*(Ber(-x31/V_T)));
                % hole_1_transient
                Jaco_tmp_Vg(3,3)=Jaco_tmp_Vg(3,3)+q*Control_Volume(1,1)/delta_t;

                % potential_2
                Jaco_tmp_Vg(4,1)=eps(Table_element_region(ii,1),1)*edge(ii,1)/L(ii,1);
                Jaco_tmp_Vg(4,4)=eps(Table_element_region(ii,1),1)*(-edge(ii,1)/L(ii,1)-edge(ii,2)/L(ii,2));
                Jaco_tmp_Vg(4,5)=-coeff*Control_Volume(2,1);
                Jaco_tmp_Vg(4,6)=coeff*Control_Volume(2,1);
                Jaco_tmp_Vg(4,7)=eps(Table_element_region(ii,1),1)*edge(ii,2)/L(ii,2);

                % electron_2_DD
                Jaco_tmp_Vg(5,1)=coeff_Jn*((edge(ii,1)/L(ii,1))*(n1/(V_T)*Ber_d(x12/V_T)-n2/(-V_T)*Ber_d(-x12/V_T)));
                Jaco_tmp_Vg(5,2)=coeff_Jn*((edge(ii,1)/L(ii,1))*(Ber(x12/V_T)));
                Jaco_tmp_Vg(5,4)=coeff_Jn*((edge(ii,2)/L(ii,2))*(n3/(-V_T)*Ber_d(x32/V_T)-n2/(V_T)*Ber_d(-x32/V_T))+((edge(ii,1)/L(ii,1))*(n1/(-V_T)*Ber_d(x12/V_T)-n2/(V_T)*Ber_d(-x12/V_T))));
                Jaco_tmp_Vg(5,5)=coeff_Jn*((edge(ii,2)/L(ii,2))*(-Ber(-x32/V_T))+((edge(ii,1)/L(ii,1))*(-Ber(-x12/V_T))));
                Jaco_tmp_Vg(5,7)=coeff_Jn*((edge(ii,2)/L(ii,2))*(n3/(V_T)*Ber_d(x32/V_T)-n2/(-V_T)*Ber_d(-x32/V_T)));
                Jaco_tmp_Vg(5,8)=coeff_Jn*((edge(ii,2)/L(ii,2))*(Ber(x32/V_T)));
                % electron_2_DD_transient
                Jaco_tmp_Vg(5,5)=Jaco_tmp_Vg(5,5)-q*Control_Volume(2,1)/delta_t;

                % hole_2_DD
                Jaco_tmp_Vg(6,1)=coeff_Jp*((edge(ii,1)/L(ii,1))*(p1/(-V_T)*Ber_d(-x12/V_T)-p2/(V_T)*Ber_d(x12/V_T)));
                Jaco_tmp_Vg(6,3)=coeff_Jp*((edge(ii,1)/L(ii,1))*(Ber(-x12/V_T)));
                Jaco_tmp_Vg(6,4)=coeff_Jp*((edge(ii,2)/L(ii,2))*(p3/(V_T)*Ber_d(-x32/V_T)-p2/(-V_T)*Ber_d(x32/V_T))+((edge(ii,1)/L(ii,1))*(p1/(V_T)*Ber_d(-x12/V_T)-p2/(-V_T)*Ber_d(x12/V_T))));
                Jaco_tmp_Vg(6,6)=coeff_Jp*((edge(ii,2)/L(ii,2))*(-Ber(x32/V_T))+(edge(ii,1)/L(ii,1))*(-Ber(x12/V_T)));
                Jaco_tmp_Vg(6,7)=coeff_Jp*((edge(ii,2)/L(ii,2))*(p3/(-V_T)*Ber_d(-x32/V_T)-p2/(V_T)*Ber_d(x32/V_T)));
                Jaco_tmp_Vg(6,9)=coeff_Jp*((edge(ii,2)/L(ii,2))*(Ber(-x32/V_T)));
                % hole_2_transient
                Jaco_tmp_Vg(6,6)=Jaco_tmp_Vg(6,6)+q*Control_Volume(2,1)/delta_t;

                % potential_3
                Jaco_tmp_Vg(7,1)=eps(Table_element_region(ii,1),1)*edge(ii,3)/L(ii,3);
                Jaco_tmp_Vg(7,4)=eps(Table_element_region(ii,1),1)*edge(ii,2)/L(ii,2);
                Jaco_tmp_Vg(7,7)=eps(Table_element_region(ii,1),1)*(-edge(ii,2)/L(ii,2)-edge(ii,3)/L(ii,3));
                Jaco_tmp_Vg(7,8)=-coeff*Control_Volume(3,1);
                Jaco_tmp_Vg(7,9)=coeff*Control_Volume(3,1);

                % electron_3_DD
                Jaco_tmp_Vg(8,1)=coeff_Jn*((edge(ii,3)/L(ii,3))*(n1/(V_T)*Ber_d(x13/V_T)-n3/(-V_T)*Ber_d(-x13/V_T)));
                Jaco_tmp_Vg(8,2)=coeff_Jn*((edge(ii,3)/L(ii,3))*(Ber(x13/V_T)));
                Jaco_tmp_Vg(8,4)=coeff_Jn*((edge(ii,2)/L(ii,2))*(n2/(V_T)*Ber_d(x23/V_T)-n3/(-V_T)*Ber_d(-x23/V_T)));
                Jaco_tmp_Vg(8,5)=coeff_Jn*((edge(ii,2)/L(ii,2))*(Ber(x23/V_T)));
                Jaco_tmp_Vg(8,7)=coeff_Jn*((edge(ii,3)/L(ii,3))*(n1/(-V_T)*Ber_d(x13/V_T)-n3/(V_T)*Ber_d(-x13/V_T))+((edge(ii,2)/L(ii,2))*(n2/(-V_T)*Ber_d(x23/V_T)-n3/(V_T)*Ber_d(-x23/V_T))));
                Jaco_tmp_Vg(8,8)=coeff_Jn*((edge(ii,3)/L(ii,3))*(-Ber(-x13/V_T))+((edge(ii,2)/L(ii,2))*(-Ber(-x23/V_T))));
                % electron_1_DD_transient
                Jaco_tmp_Vg(8,8)=Jaco_tmp_Vg(8,8)-q*Control_Volume(3,1)/delta_t;

                % hole_3_DD
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

                % Res vector %
                % potnetial
                res_tmp_Vg=zeros(9,1);
                res_potential_tmp1_Vg(1,1)=eps(Table_element_region(ii,1),1)*(edge(ii,1)/L(ii,1)*(solution_vector_Vg(index(4,1),1)-solution_vector_Vg(index(1,1),1))+edge(ii,3)/L(ii,3)*(solution_vector_Vg(index(7,1),1)-solution_vector_Vg(index(1,1),1)));
                res_potential_tmp1_Vg(2,1)=eps(Table_element_region(ii,1),1)*(edge(ii,2)/L(ii,2)*(solution_vector_Vg(index(7,1),1)-solution_vector_Vg(index(4,1),1))+edge(ii,1)/L(ii,1)*(solution_vector_Vg(index(1,1),1)-solution_vector_Vg(index(4,1),1)));
                res_potential_tmp1_Vg(3,1)=eps(Table_element_region(ii,1),1)*(edge(ii,3)/L(ii,3)*(solution_vector_Vg(index(1,1),1)-solution_vector_Vg(index(7,1),1))+edge(ii,2)/L(ii,2)*(solution_vector_Vg(index(4,1),1)-solution_vector_Vg(index(7,1),1)));
                res_potential_tmp2_Vg=zeros(3,1);
                for j=1:3
                    res_potential_tmp2_Vg(j,1)=(N_dop(Table_element_region(ii,1),1)-solution_vector_Vg(index(3*j-1,1),1)+solution_vector_Vg(index(3*j,1),1))*coeff*Control_Volume(j,1);
                end
                res_potential_Vg=res_potential_tmp1_Vg+res_potential_tmp2_Vg;

                % Jn %
                % DD
                res_Jn=zeros(3,1);
                res_Jn(1,1)=coeff_Jn*((edge(ii,1)/L(ii,1))*(n2*Ber(x21/V_T)-n1*Ber(-x21/V_T))+(edge(ii,3)/L(ii,3))*(n3*Ber(x31/V_T)-n1*Ber(-x31/V_T)));
                res_Jn(2,1)=coeff_Jn*((edge(ii,2)/L(ii,2))*(n3*Ber(x32/V_T)-n2*Ber(-x32/V_T))+(edge(ii,1)/L(ii,1))*(n1*Ber(x12/V_T)-n2*Ber(-x12/V_T)));
                res_Jn(3,1)=coeff_Jn*((edge(ii,3)/L(ii,3))*(n1*Ber(x13/V_T)-n3*Ber(-x13/V_T))+(edge(ii,2)/L(ii,2))*(n2*Ber(x23/V_T)-n3*Ber(-x23/V_T)));
                % transient
                res_Jn(1,1)=res_Jn(1,1)-q*Control_Volume(1,1)*(n1-n1_old)/delta_t;
                res_Jn(2,1)=res_Jn(2,1)-q*Control_Volume(2,1)*(n2-n2_old)/delta_t;
                res_Jn(3,1)=res_Jn(3,1)-q*Control_Volume(3,1)*(n3-n3_old)/delta_t;


                % Jp %
                % DD
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
            res_Vg(BC_index,1)=solution_vector_Vg(BC_index,1)-V_T*log(Nd/n_int)-Vd-Vg_tmp;
            res_Vg(BC_index+1,1)=solution_vector_Vg(BC_index+1,1)-Nd;
            res_Vg(BC_index+2,1)=solution_vector_Vg(BC_index+2,1)-n_int^2/Nd;
        end

        % Scaling %
        Cvector_Vg=zeros(size(solution_vector_Vg,1),1); v=[V_T; max(abs(N_dop)); max(abs(N_dop))];
        for ii=1:size(solution_vector_Vg,1)
            Cvector_Vg(ii,1)=v(Table_Jaco(ii,3),1);
        end
%         Jaco_Vg=sparse(Jaco_Vg);
        Cmatrix_Vg=spdiags(Cvector_Vg,0,size(solution_vector_Vg,1),size(solution_vector_Vg,1));
        Jaco_scaled_Vg=Jaco_Vg*Cmatrix_Vg;
        Rvector_Vg=1./sum(abs(Jaco_scaled_Vg),2);
        Rmatrix_Vg=spdiags(Rvector_Vg,0,size(solution_vector_Vg,1),size(solution_vector_Vg,1));
        Jaco_scaled_Vg=Rmatrix_Vg*Jaco_scaled_Vg;
        res_scaled_Vg=Rmatrix_Vg*res_Vg;
        update_scaled_Vg=Jaco_scaled_Vg\(-res_scaled_Vg);
        update_Vg(:,Newton_Vg)=Cmatrix_Vg*update_scaled_Vg;

        % update solution vector
        solution_vector_Vg(:,1)=solution_vector_Vg(:,1)+update_Vg(:,Newton_Vg);

        % break
        update_Vg_phi_for_break=update_Vg(find(Table_Jaco(:,1)==1, 1 , 'first'):3:find(Table_Jaco(:,1)==1, 1 , 'last'),Newton_Vg);
        fprintf("Max Updatevector=%d\n\n" , max(abs(update_Vg_phi_for_break)))
        if max(abs(update_Vg_phi_for_break))<1e-15
            break;
        end

    end

    % save solution vector
    solution_vector_Vg_saved(:,TIME_for)=solution_vector_Vg(:,1);
    if TIME_save>=2
        solution_vector_Vg_old(:,1)=solution_vector_Vg_saved(:,TIME_for-1);
    end
    Newton_Vg_save(TIME_for,1)=Newton_Vg;

    % Current %
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

            Control_Volume=zeros(3,1);
            Control_Volume(1,1)=(0.25)*L(ii,3)*edge(ii,3)+(0.25)*L(ii,1)*edge(ii,1);
            Control_Volume(2,1)=(0.25)*L(ii,1)*edge(ii,1)+(0.25)*L(ii,2)*edge(ii,2);
            Control_Volume(3,1)=(0.25)*L(ii,2)*edge(ii,2)+(0.25)*L(ii,3)*edge(ii,3);


            index_element=find(Table_element_region(:,1)==1 & Table_element_region(:,2)==row(k,1));
            eps_si=11.7;
            if col==1
                I_n_tmp=coeff_Jn*((edge(index_element,1)/L(index_element,1))*(n2*Ber(x21/V_T)-n1*Ber(-x21/V_T))+(edge(index_element,3)/L(index_element,3))*(n3*Ber(x31/V_T)-n1*Ber(-x31/V_T)));
                I_p_tmp=coeff_Jp*((edge(index_element,1)/L(index_element,1))*(p2*Ber(-x21/V_T)-p1*Ber(x21/V_T))+(edge(index_element,3)/L(index_element,3))*(p3*Ber(-x31/V_T)-p1*Ber(x31/V_T)));
                I_d_tmp=-eps_si*eps0*(x21-x21_old)/delta_t/L(index_element,1)*Control_Volume(1,1)-eps_si*eps0*(x31-x31_old)/delta_t/L(index_element,3)*Control_Volume(3,1);
            elseif col==2
                I_n_tmp=coeff_Jn*((edge(index_element,2)/L(index_element,2))*(n3*Ber(x32/V_T)-n2*Ber(-x32/V_T))+(edge(index_element,1)/L(index_element,1))*(n1*Ber(x12/V_T)-n2*Ber(-x12/V_T)));
                I_p_tmp=coeff_Jp*((edge(index_element,2)/L(index_element,2))*(p3*Ber(-x32/V_T)-p2*Ber(x32/V_T))+(edge(index_element,1)/L(index_element,1))*(p1*Ber(-x12/V_T)-p2*Ber(x12/V_T)));
                I_d_tmp=-eps_si*eps0*(x32-x32_old)/delta_t/L(index_element,2)*Control_Volume(2,1)-eps_si*eps0*(x12-x12_old)/delta_t/L(index_element,1)*Control_Volume(1,1);
            elseif col==3
                I_n_tmp=coeff_Jn*((edge(index_element,3)/L(index_element,3))*(n1*Ber(x13/V_T)-n3*Ber(-x13/V_T))+(edge(index_element,2)/L(index_element,2))*(n2*Ber(x23/V_T)-n3*Ber(-x23/V_T)));
                I_p_tmp=coeff_Jp*((edge(index_element,3)/L(index_element,3))*(p1*Ber(-x13/V_T)-p3*Ber(x13/V_T))+(edge(index_element,2)/L(index_element,2))*(p2*Ber(-x23/V_T)-p3*Ber(x23/V_T)));
                I_d_tmp=-eps_si*eps0*(x23-x23_old)/delta_t/L(index_element,2)*Control_Volume(2,1)-eps_si*eps0*(x13-x13_old)/delta_t/L(index_element,3)*Control_Volume(3,1);
            end
            I_n(TIME_for,1)=I_n(TIME_for,1)+I_n_tmp;
            I_p(TIME_for,1)=I_p(TIME_for,1)+I_p_tmp;
            I_d(TIME_for,1)=I_d(TIME_for,1)+I_d_tmp;

        end
    end
    I(TIME_for,1) = (I_n(TIME_for,1)+I_p(TIME_for,1)+I_d(TIME_for,1))*1e-6;
    update_Vg_phi=update_Vg(find(Table_Jaco(:,1)==1, 1 , 'first'):3:find(Table_Jaco(:,1)==1, 1 , 'last'),:);

%     % Save
%     if TIME_for==1 % 10번마다 1번씩 저장하도록 설정.
%         FILENAME = sprintf('./data/Homework16_Transient/%s/Homework16_Vd_%.2fV_T_%.2es_Amp_%.2fV_Dutycycle_%.2f_Slope_%.2f.mat',DIRNAME,Vd,T,Amp,Duty_cycle,Slope);
%         save(FILENAME);
% 
%         FILENAME = sprintf('./data/Homework16_Transient/%s/Homework16_Vd_%.2fV_T_%.2es_Amp_%.2fV_Dutycycle_%.2f_Slope_%.2f_current.mat',DIRNAME,Vd,T,Amp,Duty_cycle,Slope);
%         save(FILENAME, 'I');
%     elseif count_save==10
%         FILENAME = sprintf('./data/Homework16_Transient/%s/Homework16_Vd_%.2fV_T_%.2es_Amp_%.2fV_Dutycycle_%.2f_Slope_%.2f.mat',DIRNAME,Vd,T,Amp,Duty_cycle,Slope);
%         save(FILENAME);
% 
%         FILENAME = sprintf('./data/Homework16_Transient/%s/Homework16_Vd_%.2fV_T_%.2es_Amp_%.2fV_Dutycycle_%.2f_Slope_%.2f_current.mat',DIRNAME,Vd,T,Amp,Duty_cycle,Slope);
%         save(FILENAME, 'I');
%         count_save=0;
%     end
%     count_save=count_save+1;
end
%% save

FILENAME = sprintf('%s/Homework16_Vd_%.2fV_T_%.2es_Amp_%.2fV.mat',Make_DIRNAME,Vd,T,Amp);
save(FILENAME);

FILENAME = sprintf('%s/Homework16_Vd_%.2fV_T_%.2es_Amp_%.2fV_current.mat',Make_DIRNAME,Vd,T,Amp);
save(FILENAME, 'I', 'time','Vg');
count_save=0;

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
%         patch('Faces',Element_si,'Vertices',Vertex(:,[1 2]), 'EdgeColor','black','FaceColor','yellow','LineWidth',0.5)
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