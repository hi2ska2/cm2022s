clear; clc; close all;


% Load할 Vd 지정
load('matlab.mat')

% 제작한 pulse function을 활용했다.
% [Voltage_save, time_save] = pulse(T, Time_step, Amp, Duty_cycle, Slope, Period)

f=1e6;
T=1/f; % [ns], 주기
Time_step=100;
Amp=1e-3;
period=2;

coeff=q/eps0;
coeff_dJn=-q*mobility_n;
coeff_dJp=-q*mobility_p;

%% initial value
solution_vector_Sg=zeros(length(Table_Jaco),1);
for ii=1:length(Table_Jaco)
    if Table_Jaco(ii,3)==1 && (Table_Jaco(ii,1)==1 || Table_Jaco(ii,1)==5)
        solution_vector_Sg(ii,1)=0;
    elseif Table_Jaco(ii,3)==1 && (Table_Jaco(ii,1)==2 || Table_Jaco(ii,1)==3 || Table_Jaco(ii,1)==4)
        solution_vector_Sg(ii,1)=0;
    elseif Table_Jaco(ii,3)==2 || Table_Jaco(ii,3)==3
        solution_vector_Sg(ii,1)=0;
    end
end

solution_vector_Sg_saved(:,1)=solution_vector_Sg(:,1);

%% Small signal
Iteration_Sg=100000;
delta_t=T/Time_step;
time=0:delta_t:T*period;

% Jacobian matrix / res vector
count_Sg_current=1;
solution_vector_Sg_old=solution_vector_Sg;
for TIME_for=1:length(time)
    TIME_save=TIME_for;
    Sg(TIME_for,1)=Amp*sin(2*pi*f*time(TIME_for));
%     solution_vector_Sg_old=solution_vector_Sg;

    update_Sg=zeros(size(solution_vector_Sg,1),Iteration_Sg);
    for Newton_Sg=1:Iteration_Sg
        fprintf("Input Small Signal, Frequency=%.3e, Sg=%.6f, Newton_Sg=%d\n" , f, Sg(TIME_for,1), Newton_Sg)
        %         Jaco_Sg=sparse(zeros(size(solution_vector_Sg,1),size(solution_vector_Sg,1))); res_Sg=zeros(size(solution_vector_Sg,1),1);
        Jaco_Sg=zeros(size(solution_vector_Sg,1),size(solution_vector_Sg,1)); res_Sg=zeros(size(solution_vector_Sg,1),1); %% 희소행렬 아님. Jaco matrix 확인용.

        % Jaco/res %
        for ii=1:length(Element)
            % Control_Volume 계산
            Control_Volume=zeros(3,1);
            Control_Volume(1,1)=(0.25)*L(ii,3)*edge(ii,3)+(0.25)*L(ii,1)*edge(ii,1);
            Control_Volume(2,1)=(0.25)*L(ii,1)*edge(ii,1)+(0.25)*L(ii,2)*edge(ii,2);
            Control_Volume(3,1)=(0.25)*L(ii,2)*edge(ii,2)+(0.25)*L(ii,3)*edge(ii,3);

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

                % potential, elec, hole 계산 %
                %%% DC %%%
                % potential
                x1_dc=solution_vector_Vg(index(1,1),1);
                x2_dc=solution_vector_Vg(index(4,1),1);
                x3_dc=solution_vector_Vg(index(7,1),1);
                x12_dc=x1_dc-x2_dc; x21_dc=-x12_dc;
                x23_dc=x2_dc-x3_dc; x32_dc=-x23_dc;
                x31_dc=x3_dc-x1_dc; x13_dc=-x31_dc;
                x12_dc_avr=(x1_dc+x2_dc)/2; x23_dc_avr=(x2_dc+x3_dc)/2; x31_dc_avr=(x3_dc+x1_dc)/2;
                x21_dc_avr=x12_dc_avr; x32_dc_avr=x23_dc_avr; x13_dc_avr=x31_dc_avr;

                % electron
                n1_dc=solution_vector_Vg(index(1,1)+1,1);
                n2_dc=solution_vector_Vg(index(4,1)+1,1);
                n3_dc=solution_vector_Vg(index(7,1)+1,1);
                n12_dc=n1_dc-n2_dc; n21_dc=-n12_dc;
                n23_dc=n2_dc-n3_dc; n32_dc=-n23_dc;
                n31_dc=n3_dc-n1_dc; n13_dc=-n31_dc;
                n12_dc_avr=(n1_dc+n2_dc)/2; n23_dc_avr=(n2_dc+n3_dc)/2; n31_dc_avr=(n3_dc+n1_dc)/2;
                n21_dc_avr=n12_dc_avr; n32_dc_avr=n23_dc_avr; n13_dc_avr=n31_dc_avr;

                % hole
                p1_dc=solution_vector_Vg(index(1,1)+2,1);
                p2_dc=solution_vector_Vg(index(4,1)+2,1);
                p3_dc=solution_vector_Vg(index(7,1)+2,1);
                p12_dc=p1_dc-p2_dc; p21_dc=-p12_dc;
                p23_dc= p2_dc-p3_dc; p32_dc=-p23_dc;
                p31_dc= p3_dc-p1_dc; p13_dc=-p31_dc;
                p12_dc_avr=(p1_dc+p2_dc)/2; p23_dc_avr=(p2_dc+p3_dc)/2; p31_dc_avr=(p3_dc+p1_dc)/2;
                p21_dc_avr=p12_dc_avr; p32_dc_avr=p23_dc_avr; p13_dc_avr=p31_dc_avr;

                %%% small signal %%%
                % potential
                x1_delta=solution_vector_Sg(index(1,1),1);
                x2_delta=solution_vector_Sg(index(4,1),1);
                x3_delta=solution_vector_Sg(index(7,1),1);
                x12_delta=x1_delta-x2_delta; x21_delta=-x12_delta;
                x23_delta=x2_delta-x3_delta; x32_delta=-x23_delta;
                x31_delta=x3_delta-x1_delta; x13_delta=-x31_delta;
                x12_delta_avr=(x1_delta+x2_delta)/2; x23_delta_avr=(x2_delta+x3_delta)/2; x31_delta_avr=(x3_delta+x1_delta)/2;
                x21_delta_avr=x12_delta_avr; x32_delta_avr=x23_delta_avr; x13_delta_avr=x31_delta_avr;

                % electron
                n1_delta=solution_vector_Sg(index(1,1)+1,1);
                n2_delta=solution_vector_Sg(index(4,1)+1,1);
                n3_delta=solution_vector_Sg(index(7,1)+1,1);
                n12_delta=n1_delta-n2_delta; n21_delta=-n12_delta;
                n23_delta=n2_delta-n3_delta; n32_delta=-n23_delta;
                n31_delta=n3_delta-n1_delta; n13_delta=-n31_delta;
                n12_delta_avr=(n1_delta+n2_delta)/2; n23_delta_avr=(n2_delta+n3_delta)/2; n31_delta_avr=(n3_delta+n1_delta)/2;
                n21_delta_avr=n12_delta_avr; n32_delta_avr=n23_delta_avr; n13_delta_avr=n31_delta_avr;

                % hole
                p1_delta=solution_vector_Sg(index(1,1)+2,1);
                p2_delta=solution_vector_Sg(index(4,1)+2,1);
                p3_delta=solution_vector_Sg(index(7,1)+2,1);
                p12_delta=p1_delta-p2_delta; p21_delta=-p12_delta;
                p23_delta= p2_delta-p3_delta; p32_delta=-p23_delta;
                p31_delta= p3_delta-p1_delta; p13_delta=-p31_delta;
                p12_delta_avr=(p1_delta+p2_delta)/2; p23_delta_avr=(p2_delta+p3_delta)/2; p31_delta_avr=(p3_delta+p1_delta)/2;
                p21_delta_avr=p12_delta_avr; p32_delta_avr=p23_delta_avr; p13_delta_avr=p31_delta_avr;

                %%% old %%%
                % potential
                x1_delta_old=solution_vector_Sg_old(index(1,1),1);
                x2_delta_old=solution_vector_Sg_old(index(4,1),1);
                x3_delta_old=solution_vector_Sg_old(index(7,1),1);
                x12_delta_old=x1_delta_old-x2_delta_old; x21_delta_old=-x12_delta_old;
                x23_delta_old=x2_delta_old-x3_delta_old; x32_delta_old=-x23_delta_old;
                x31_delta_old=x3_delta_old-x1_delta_old; x13_delta_old=-x31_delta_old;
                x12_delta_old_avr=(x1_delta_old+x2_delta_old)/2; x23_delta_old_avr=(x2_delta_old+x3_delta_old)/2; x31_delta_old_avr=(x3_delta_old+x1_delta_old)/2;
                x21_delta_old_avr=x12_delta_old_avr; x32_delta_old_avr=x23_delta_old_avr; x13_delta_old_avr=x31_delta_old_avr;

                % electron
                n1_delta_old=solution_vector_Sg_old(index(1,1)+1,1);
                n2_delta_old=solution_vector_Sg_old(index(4,1)+1,1);
                n3_delta_old=solution_vector_Sg_old(index(7,1)+1,1);
                n12_delta_old=n1_delta_old-n2_delta_old; n21_delta_old=-n12_delta_old;
                n23_delta_old=n2_delta_old-n3_delta_old; n32_delta_old=-n23_delta_old;
                n31_delta_old=n3_delta_old-n1_delta_old; n13_delta_old=-n31_delta_old;
                n12_delta_old_avr=(n1_delta_old+n2_delta_old)/2; n23_delta_old_avr=(n2_delta_old+n3_delta_old)/2; n31_delta_old_avr=(n3_delta_old+n1_delta_old)/2;
                n21_delta_old_avr=n12_delta_old_avr; n32_delta_old_avr=n23_delta_old_avr; n13_delta_old_avr=n31_delta_old_avr;

                % hole
                p1_delta_old=solution_vector_Sg_old(index(1,1)+2,1);
                p2_delta_old=solution_vector_Sg_old(index(4,1)+2,1);
                p3_delta_old=solution_vector_Sg_old(index(7,1)+2,1);
                p12_delta_old=p1_delta_old-p2_delta_old; p21_delta_old=-p12_delta_old;
                p23_delta_old= p2_delta_old-p3_delta_old; p32_delta_old=-p23_delta_old;
                p31_delta_old= p3_delta_old-p1_delta_old; p13_delta_old=-p31_delta_old;
                p12_delta_old_avr=(p1_delta_old+p2_delta_old)/2; p23_delta_old_avr=(p2_delta_old+p3_delta_old)/2; p31_delta_old_avr=(p3_delta_old+p1_delta_old)/2;
                p21_delta_old_avr=p12_delta_old_avr; p32_delta_old_avr=p23_delta_old_avr; p13_delta_old_avr=p31_delta_old_avr;


                %% res vector
                res_tmp_Sg=zeros(9,1);

                % Potential-DD
                res_potential_Sg=zeros(3,1);
                res_potential_Sg(1,1)=res_potential_Sg(1,1)+eps(Table_element_region(ii,1),1)*(edge(ii,1)/L(ii,1)*(x2_delta-x1_delta)+edge(ii,3)/L(ii,3)*(x3_delta-x1_delta));
                res_potential_Sg(1,1)=res_potential_Sg(1,1)+(-n1_delta+p1_delta)*coeff*Control_Volume(1,1);
                res_potential_Sg(2,1)=res_potential_Sg(2,1)+eps(Table_element_region(ii,1),1)*(edge(ii,2)/L(ii,2)*(x3_delta-x2_delta)+edge(ii,1)/L(ii,1)*(x1_delta-x2_delta));
                res_potential_Sg(2,1)=res_potential_Sg(2,1)+(-n2_delta+p2_delta)*coeff*Control_Volume(2,1);
                res_potential_Sg(3,1)=res_potential_Sg(3,1)+eps(Table_element_region(ii,1),1)*(edge(ii,3)/L(ii,3)*(x1_delta-x3_delta)+edge(ii,2)/L(ii,2)*(x2_delta-x3_delta));
                res_potential_Sg(3,1)=res_potential_Sg(3,1)+(-n3_delta+p3_delta)*coeff*Control_Volume(3,1);

                % Jn-DD
                res_Jn_Sg=zeros(3,1);
                res_Jn_Sg(1,1)=res_Jn_Sg(1,1)+coeff_dJn*((edge(ii,1)/L(ii,1))*((n21_dc_avr*x21_delta)+(n21_delta_avr*x21_dc)-V_T*n21_delta));
                res_Jn_Sg(1,1)=res_Jn_Sg(1,1)+coeff_dJn*((edge(ii,3)/L(ii,3))*((n31_dc_avr*x31_delta)+(n31_delta_avr*x31_dc)-V_T*n31_delta));
                res_Jn_Sg(2,1)=res_Jn_Sg(2,1)+coeff_dJn*((edge(ii,2)/L(ii,2))*((n32_dc_avr*x32_delta)+(n32_delta_avr*x32_dc)-V_T*n32_delta));
                res_Jn_Sg(2,1)=res_Jn_Sg(2,1)+coeff_dJn*((edge(ii,1)/L(ii,1))*((n12_dc_avr*x12_delta)+(n12_delta_avr*x12_dc)-V_T*n12_delta));
                res_Jn_Sg(3,1)=res_Jn_Sg(3,1)+coeff_dJn*((edge(ii,3)/L(ii,3))*((n13_dc_avr*x13_delta)+(n13_delta_avr*x13_dc)-V_T*n13_delta));
                res_Jn_Sg(3,1)=res_Jn_Sg(3,1)+coeff_dJn*((edge(ii,2)/L(ii,2))*((n23_dc_avr*x23_delta)+(n23_delta_avr*x23_dc)-V_T*n23_delta));
                % Jn-transient
                res_Jn_Sg(1,1)=res_Jn_Sg(1,1)-q*Control_Volume(1,1)*(n1_delta-n1_delta_old)/delta_t;
                res_Jn_Sg(2,1)=res_Jn_Sg(2,1)-q*Control_Volume(2,1)*(n2_delta-n2_delta_old)/delta_t;
                res_Jn_Sg(3,1)=res_Jn_Sg(3,1)-q*Control_Volume(3,1)*(n3_delta-n3_delta_old)/delta_t;

                % Jp-DD
                res_Jp_Sg=zeros(3,1);
                res_Jp_Sg(1,1)=res_Jp_Sg(1,1)+coeff_dJp*((edge(ii,1)/L(ii,1))*((p21_dc_avr*x21_delta)+(p21_delta_avr*x21_dc)+V_T*p21_delta));
                res_Jp_Sg(1,1)=res_Jp_Sg(1,1)+coeff_dJp*((edge(ii,3)/L(ii,3))*((p31_dc_avr*x31_delta)+(p31_delta_avr*x31_dc)+V_T*p31_delta));
                res_Jp_Sg(2,1)=res_Jp_Sg(2,1)+coeff_dJp*((edge(ii,2)/L(ii,2))*((p32_dc_avr*x32_delta)+(p32_delta_avr*x32_dc)+V_T*p32_delta));
                res_Jp_Sg(2,1)=res_Jp_Sg(2,1)+coeff_dJp*((edge(ii,1)/L(ii,1))*((p12_dc_avr*x12_delta)+(p12_delta_avr*x12_dc)+V_T*p12_delta));
                res_Jp_Sg(3,1)=res_Jp_Sg(3,1)+coeff_dJp*((edge(ii,3)/L(ii,3))*((p13_dc_avr*x13_delta)+(p13_delta_avr*x13_dc)+V_T*p13_delta));
                res_Jp_Sg(3,1)=res_Jp_Sg(3,1)+coeff_dJp*((edge(ii,2)/L(ii,2))*((p23_dc_avr*x23_delta)+(p23_delta_avr*x23_dc)+V_T*p23_delta));
                % Jp-transient
                res_Jp_Sg(1,1)=res_Jp_Sg(1,1)+q*Control_Volume(1,1)*(p1_delta-p1_delta_old)/delta_t;
                res_Jp_Sg(2,1)=res_Jp_Sg(2,1)+q*Control_Volume(2,1)*(p2_delta-p2_delta_old)/delta_t;
                res_Jp_Sg(3,1)=res_Jp_Sg(3,1)+q*Control_Volume(3,1)*(p3_delta-p3_delta_old)/delta_t;


                n=1;
                for j=1:3
                    res_tmp_Sg(n,1)=res_potential_Sg(j,1);
                    res_tmp_Sg(n+1,1)=res_Jn_Sg(j,1);
                    res_tmp_Sg(n+2,1)=res_Jp_Sg(j,1);
                    n=n+3;
                end

                for j=1:9
                    res_Sg(index(j,1),1)=res_Sg(index(j,1),1)+res_tmp_Sg(j,1);
                end


                %% Jacobian matrix
                Jaco_tmp_Sg=zeros(9,9);

                % node1 - potential
                Jaco_tmp_Sg(1,1)=eps(Table_element_region(ii,1),1)*(-edge(ii,1)/L(ii,1)-edge(ii,3)/L(ii,3));
                Jaco_tmp_Sg(1,2)=-coeff*Control_Volume(1,1);
                Jaco_tmp_Sg(1,3)=coeff*Control_Volume(1,1);
                Jaco_tmp_Sg(1,4)=eps(Table_element_region(ii,1),1)*edge(ii,1)/L(ii,1);
                Jaco_tmp_Sg(1,7)=eps(Table_element_region(ii,1),1)*edge(ii,3)/L(ii,3);

                % node1 - electron
                Jaco_tmp_Sg(2,1)=coeff_dJn*((edge(ii,1)/L(ii,1))*(-n21_dc_avr)+(edge(ii,3)/L(ii,3))*(-n31_dc_avr));
                Jaco_tmp_Sg(2,2)=coeff_dJn*((edge(ii,1)/L(ii,1))*(0.5*x21_dc+V_T)+(edge(ii,3)/L(ii,3))*(0.5*x31_dc+V_T));
                Jaco_tmp_Sg(2,4)=coeff_dJn*((edge(ii,1)/L(ii,1))*(n21_dc_avr));
                Jaco_tmp_Sg(2,5)=coeff_dJn*((edge(ii,1)/L(ii,1))*(0.5*x21_dc-V_T));
                Jaco_tmp_Sg(2,7)=coeff_dJn*((edge(ii,3)/L(ii,3))*(n31_dc_avr));
                Jaco_tmp_Sg(2,8)=coeff_dJn*((edge(ii,3)/L(ii,3))*(0.5*x31_dc-V_T));
                % node1 - electron_transient
                Jaco_tmp_Sg(2,2)=Jaco_tmp_Sg(2,2)-q*Control_Volume(1,1)/delta_t;

                % node1 - hole
                Jaco_tmp_Sg(3,1)=coeff_dJp*((edge(ii,1)/L(ii,1))*(-n21_dc_avr)+(edge(ii,3)/L(ii,3))*(-n31_dc_avr));
                Jaco_tmp_Sg(3,3)=coeff_dJp*((edge(ii,1)/L(ii,1))*(0.5*x21_dc-V_T)+(edge(ii,3)/L(ii,3))*(0.5*x31_dc-V_T));
                Jaco_tmp_Sg(3,4)=coeff_dJp*((edge(ii,1)/L(ii,1))*(n21_dc_avr));
                Jaco_tmp_Sg(3,6)=coeff_dJp*((edge(ii,1)/L(ii,1))*(0.5*x21_dc+V_T));
                Jaco_tmp_Sg(3,7)=coeff_dJp*((edge(ii,3)/L(ii,3))*(n31_dc_avr));
                Jaco_tmp_Sg(3,9)=coeff_dJp*((edge(ii,3)/L(ii,3))*(0.5*x31_dc+V_T));
                % node1 - hole_transient
                Jaco_tmp_Sg(3,3)=Jaco_tmp_Sg(3,3)+q*Control_Volume(1,1)/delta_t;

                % node2 - potential
                Jaco_tmp_Sg(4,1)=eps(Table_element_region(ii,1),1)*edge(ii,1)/L(ii,1);
                Jaco_tmp_Sg(4,4)=eps(Table_element_region(ii,1),1)*(-edge(ii,1)/L(ii,1)-edge(ii,2)/L(ii,2));
                Jaco_tmp_Sg(4,5)=-coeff*Control_Volume(2,1);
                Jaco_tmp_Sg(4,6)=coeff*Control_Volume(2,1);
                Jaco_tmp_Sg(4,7)=eps(Table_element_region(ii,1),1)*edge(ii,2)/L(ii,2);

                % node2 - electron
                Jaco_tmp_Sg(5,1)=coeff_dJn*((edge(ii,1)/L(ii,1))*(n12_dc_avr));
                Jaco_tmp_Sg(5,2)=coeff_dJn*((edge(ii,1)/L(ii,1))*(0.5*x12_dc-V_T));
                Jaco_tmp_Sg(5,4)=coeff_dJn*((edge(ii,2)/L(ii,2))*(-n32_dc_avr)+(edge(ii,1)/L(ii,1))*(-n12_dc_avr));
                Jaco_tmp_Sg(5,5)=coeff_dJn*((edge(ii,2)/L(ii,2))*(0.5*x32_dc+V_T)+(edge(ii,1)/L(ii,1))*(0.5*x12_dc+V_T));
                Jaco_tmp_Sg(5,7)=coeff_dJn*((edge(ii,2)/L(ii,2))*(n32_dc_avr));
                Jaco_tmp_Sg(5,8)=coeff_dJn*((edge(ii,2)/L(ii,2))*(0.5*x32_dc-V_T));
                % node2 - electron_transient
                Jaco_tmp_Sg(5,5)=Jaco_tmp_Sg(5,5)-q*Control_Volume(2,1)/delta_t;

                % node2 - hole
                Jaco_tmp_Sg(6,1)=coeff_dJp*((edge(ii,1)/L(ii,1))*(n12_dc_avr));
                Jaco_tmp_Sg(6,3)=coeff_dJp*((edge(ii,1)/L(ii,1))*(0.5*x12_dc+V_T));
                Jaco_tmp_Sg(6,4)=coeff_dJp*((edge(ii,2)/L(ii,2))*(-n32_dc_avr)+(edge(ii,1)/L(ii,1))*(-n12_dc_avr));
                Jaco_tmp_Sg(6,6)=coeff_dJp*((edge(ii,2)/L(ii,2))*(0.5*x32_dc-V_T)+(edge(ii,1)/L(ii,1))*(0.5*x12_dc-V_T));
                Jaco_tmp_Sg(6,7)=coeff_dJp*((edge(ii,2)/L(ii,2))*(n32_dc_avr));
                Jaco_tmp_Sg(6,9)=coeff_dJp*((edge(ii,2)/L(ii,2))*(0.5*x32_dc+V_T));
                % node2 - hole_transient
                Jaco_tmp_Sg(6,6)=Jaco_tmp_Sg(6,6)+q*Control_Volume(2,1)/delta_t;

                % node3 - potential
                Jaco_tmp_Sg(7,1)=eps(Table_element_region(ii,1),1)*edge(ii,3)/L(ii,3);
                Jaco_tmp_Sg(7,4)=eps(Table_element_region(ii,1),1)*edge(ii,2)/L(ii,2);
                Jaco_tmp_Sg(7,7)=eps(Table_element_region(ii,1),1)*(-edge(ii,2)/L(ii,2)-edge(ii,3)/L(ii,3));
                Jaco_tmp_Sg(7,8)=-coeff*Control_Volume(3,1);
                Jaco_tmp_Sg(7,9)=coeff*Control_Volume(3,1);

                % node3 - electron
                Jaco_tmp_Sg(8,1)=coeff_dJn*((edge(ii,3)/L(ii,3))*(n13_dc_avr));
                Jaco_tmp_Sg(8,2)=coeff_dJn*((edge(ii,3)/L(ii,3))*(0.5*x13_dc-V_T));
                Jaco_tmp_Sg(8,4)=coeff_dJn*((edge(ii,2)/L(ii,2))*(n23_dc_avr));
                Jaco_tmp_Sg(8,5)=coeff_dJn*((edge(ii,2)/L(ii,2))*(0.5*x23_dc-V_T));
                Jaco_tmp_Sg(8,7)=coeff_dJn*((edge(ii,3)/L(ii,3))*(-n13_dc_avr)+(edge(ii,2)/L(ii,2))*(-n23_dc_avr));
                Jaco_tmp_Sg(8,8)=coeff_dJn*((edge(ii,3)/L(ii,3))*(0.5*x13_dc+V_T)+(edge(ii,2)/L(ii,2))*(0.5*x23_dc+V_T));
                % node3 - electron_transient
                Jaco_tmp_Sg(8,8)=Jaco_tmp_Sg(8,8)-q*Control_Volume(3,1)/delta_t;

                % node3 - hole
                Jaco_tmp_Sg(9,1)=coeff_dJp*((edge(ii,3)/L(ii,3))*(n13_dc_avr));
                Jaco_tmp_Sg(9,3)=coeff_dJp*((edge(ii,3)/L(ii,3))*(0.5*x13_dc+V_T));
                Jaco_tmp_Sg(9,4)=coeff_dJp*((edge(ii,2)/L(ii,2))*(n23_dc_avr));
                Jaco_tmp_Sg(9,6)=coeff_dJp*((edge(ii,2)/L(ii,2))*(0.5*x23_dc+V_T));
                Jaco_tmp_Sg(9,7)=coeff_dJp*((edge(ii,3)/L(ii,3))*(-n13_dc_avr)+(edge(ii,2)/L(ii,2))*(-n23_dc_avr));
                Jaco_tmp_Sg(9,9)=coeff_dJp*((edge(ii,3)/L(ii,3))*(0.5*x13_dc-V_T)+(edge(ii,2)/L(ii,2))*(0.5*x23_dc-V_T));
                % node3 - hole_transient
                Jaco_tmp_Sg(9,9)=Jaco_tmp_Sg(9,9)+q*Control_Volume(3,1)/delta_t;


                for j=1:9
                    for k=1:9
                        Jaco_Sg(index(j,1),index(k,1))=Jaco_Sg(index(j,1),index(k,1))+Jaco_tmp_Sg(j,k);
                    end
                end


            else
                % index : 현재 element에 속해있는 vertex가 resdue에 몇번째 있는가 저장
                index=zeros(3,1);
                n=1;
                for j=1:3
                    index(j,1)=find(Table_Jaco(:,1)==Table_element_region(ii,1) & Table_Jaco(:,2)==Element(ii,j) & Table_Jaco(:,3)==1);
                end

                %%% small signal %%%
                % potential
                x1_delta=solution_vector_Sg(index(1,1),1);
                x2_delta=solution_vector_Sg(index(2,1),1);
                x3_delta=solution_vector_Sg(index(3,1),1);

                %% res vector
                res_tmp_Sg(1,1)=eps_now*((-edge(ii,1)/L(ii,1)-edge(ii,3)/L(ii,3))*x1_delta+edge(ii,1)/L(ii,1)*x2_delta+edge(ii,3)/L(ii,3)*x3_delta);
                res_tmp_Sg(2,1)=eps_now*(edge(ii,1)/L(ii,1)*x1_delta+(-edge(ii,1)/L(ii,1)-edge(ii,2)/L(ii,2))*x2_delta+edge(ii,2)/L(ii,2)*x3_delta);
                res_tmp_Sg(3,1)=eps_now*(edge(ii,3)/L(ii,3)*x1_delta+edge(ii,2)/L(ii,2)*x2_delta+(-edge(ii,2)/L(ii,2)-edge(ii,3)/L(ii,3))*x3_delta);
                for j=1:3
                    res_Sg(index(j,1),1)=res_Sg(index(j,1),1)+res_tmp_Sg(j,1);
                end

                %% Jaco_potential matrix
                Jaco_tmp_Sg=zeros(3,3);

                Jaco_tmp_Sg(1,1)=eps_now*(-edge(ii,1)/L(ii,1)-edge(ii,3)/L(ii,3));
                Jaco_tmp_Sg(1,2)=eps_now*(edge(ii,1)/L(ii,1));
                Jaco_tmp_Sg(1,3)=eps_now*(edge(ii,3)/L(ii,3));
                Jaco_tmp_Sg(2,1)=eps_now*(edge(ii,1)/L(ii,1));
                Jaco_tmp_Sg(2,2)=eps_now*(-edge(ii,1)/L(ii,1)-edge(ii,2)/L(ii,2));
                Jaco_tmp_Sg(2,3)=eps_now*(edge(ii,2)/L(ii,2));
                Jaco_tmp_Sg(3,1)=eps_now*(edge(ii,3)/L(ii,3));
                Jaco_tmp_Sg(3,2)=eps_now*(edge(ii,2)/L(ii,2));
                Jaco_tmp_Sg(3,3)=eps_now*(-edge(ii,2)/L(ii,2)-edge(ii,3)/L(ii,3));

                for j=1:3
                    for k=1:3
                        Jaco_Sg(index(j,1),index(k,1))=Jaco_Sg(index(j,1),index(k,1))+Jaco_tmp_Sg(j,k);
                    end
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
                    if length(sort(nonzeros(Jaco_Sg(index_ox,:))))>2
                        Jaco_Sg(index_si,:)=Jaco_Sg(index_si,:)+Jaco_Sg(index_ox,:);
                        Jaco_Sg(index_ox,:)=0; Jaco_Sg(index_ox,index_ox)=1; Jaco_Sg(index_ox,index_si)=-1;
                        res_Sg(index_si,1)=res_Sg(index_si,1)+res_Sg(index_ox,1);
                        res_Sg(index_ox,1)=0;
                    elseif ~isequal(sort(nonzeros(Jaco_Sg(index_ox,:))),[-1; 1]) % 점 3개가 접해있는 지역에서 중복적으로 합쳐지는것을 막기위해 추가함. jaco의 ox-channel행이 -1과1이라면 이미 다른행과 합쳐진 것이므로 pass.
                        Jaco_Sg(index_si,:)=Jaco_Sg(index_si,:)+Jaco_Sg(index_ox,:);
                        Jaco_Sg(index_ox,:)=0; Jaco_Sg(index_ox,index_ox)=1; Jaco_Sg(index_ox,index_si)=-1;
                        res_Sg(index_si,1)=res_Sg(index_si,1)+res_Sg(index_ox,1);
                        res_Sg(index_ox,1)=0;
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
                    Jaco_Sg(index_channel_potential,:)=Jaco_Sg(index_channel_potential,:)+Jaco_Sg(index_source_drain_potential,:);
                    Jaco_Sg(index_channel_electron,:)=Jaco_Sg(index_channel_electron,:)+Jaco_Sg(index_source_drain_electron,:);
                    Jaco_Sg(index_channel_hole,:)=Jaco_Sg(index_channel_hole,:)+Jaco_Sg(index_source_drain_hole,:);
                    Jaco_Sg(index_source_drain_potential,:)=0; Jaco_Sg(index_source_drain_potential,index_source_drain_potential)=1; Jaco_Sg(index_source_drain_potential,index_channel_potential)=-1;
                    Jaco_Sg(index_source_drain_electron,:)=0; Jaco_Sg(index_source_drain_electron,index_source_drain_electron)=1; Jaco_Sg(index_source_drain_electron,index_channel_electron)=-1;
                    Jaco_Sg(index_source_drain_hole,:)=0; Jaco_Sg(index_source_drain_hole,index_source_drain_hole)=1; Jaco_Sg(index_source_drain_hole,index_channel_hole)=-1;

                    % res 변경, potential/electron/hole 순서
                    res_Sg(index_channel_potential,1)=res_Sg(index_channel_potential,1)+res_Sg(index_source_drain_potential,1);
                    res_Sg(index_channel_electron,1)=res_Sg(index_channel_electron,1)+res_Sg(index_source_drain_electron,1);
                    res_Sg(index_channel_hole,1)=res_Sg(index_channel_hole,1)+res_Sg(index_source_drain_hole,1);
                    res_Sg(index_source_drain_potential,1)=0;
                    res_Sg(index_source_drain_electron,1)=0;
                    res_Sg(index_source_drain_hole,1)=0;
                end
            end
        end


        % Boundary condition %
        % Dirichlet_BC
        for ii=1:size(Dirichlet_BC,1)
            BC_index=find(Table_Jaco(:,2)==Dirichlet_BC(ii,1) & Table_Jaco(:,3)==1);
            Jaco_Sg(BC_index,:)=0;
            Jaco_Sg(BC_index,BC_index)=1;
            res_Sg(BC_index,1)=solution_vector_Sg(BC_index,1)-Sg(TIME_for,1);
        end

        % Source_BC
        for ii=1:size(Source_BC,1)
            BC_index=find(Table_Jaco(:,1)==2 & Table_Jaco(:,2)==Source_BC(ii,1) & Table_Jaco(:,3)==1);
            Jaco_Sg(BC_index,:)=0; Jaco_Sg(BC_index+1,:)=0; Jaco_Sg(BC_index+2,:)=0;
            Jaco_Sg(BC_index,BC_index)=1; Jaco_Sg(BC_index+1,BC_index+1)=1; Jaco_Sg(BC_index+2,BC_index+2)=1;
            res_Sg(BC_index,1)=solution_vector_Sg(BC_index,1)-0;
            res_Sg(BC_index+1,1)=solution_vector_Sg(BC_index,1)-0;
            res_Sg(BC_index+2,1)=solution_vector_Sg(BC_index,1)-0;
        end

        % Darin_BC
        for ii=1:size(Drain_BC,1)
            BC_index=find(Table_Jaco(:,1)==4 & Table_Jaco(:,2)==Drain_BC(ii,1) & Table_Jaco(:,3)==1);
            Jaco_Sg(BC_index,:)=0; Jaco_Sg(BC_index+1,:)=0; Jaco_Sg(BC_index+2,:)=0;
            Jaco_Sg(BC_index,BC_index)=1; Jaco_Sg(BC_index+1,BC_index+1)=1; Jaco_Sg(BC_index+2,BC_index+2)=1;
            res_Sg(BC_index,1)=solution_vector_Sg(BC_index,1)-0;
            res_Sg(BC_index+1,1)=solution_vector_Sg(BC_index,1)-0;
            res_Sg(BC_index+2,1)=solution_vector_Sg(BC_index,1)-0;
        end

        % Scaling %
        Cvector_Sg=zeros(size(solution_vector_Sg,1),1); v=[V_T; max(abs(N_dop)); max(abs(N_dop))];
        for ii=1:size(solution_vector_Sg,1)
            Cvector_Sg(ii,1)=v(Table_Jaco(ii,3),1);
        end
        Jaco_Sg=sparse(Jaco_Sg);
        Cmatrix_Sg=spdiags(Cvector_Sg,0,size(solution_vector_Sg,1),size(solution_vector_Sg,1));
        Jaco_scaled_Sg=Jaco_Sg*Cmatrix_Sg;
        Rvector_Sg=1./sum(abs(Jaco_scaled_Sg),2);
        Rmatrix_Sg=spdiags(Rvector_Sg,0,size(solution_vector_Sg,1),size(solution_vector_Sg,1));
        Jaco_scaled_Sg=Rmatrix_Sg*Jaco_scaled_Sg;
        res_scaled_Sg=Rmatrix_Sg*res_Sg;
        update_scaled_Sg=Jaco_scaled_Sg\(-res_scaled_Sg);
        update_Sg(:,Newton_Sg)=Cmatrix_Sg*update_scaled_Sg;

        % update solution vector
        solution_vector_Sg(:,1)=solution_vector_Sg(:,1)+update_Sg(:,Newton_Sg);

        % break
        update_Sg_phi_for_break=update_Sg(find(Table_Jaco(:,1)==2, 1 , 'first'):3:find(Table_Jaco(:,1)==4, 1 , 'last'),Newton_Sg);
        fprintf("Max Updatevector=%d\n\n" , max(abs(update_Sg_phi_for_break)))
        if max(abs(update_Sg_phi_for_break))<1e-10
            break;
        end

    end
    solution_vector_Sg_saved(:,Newton_Sg)=solution_vector_Sg(:,1);

    % save solution vector
    solution_vector_Sg_saved(:,TIME_for)=solution_vector_Vg(:,1);
    if TIME_save>=2
        solution_vector_Sg_old(:,1)=solution_vector_Sg_saved(:,TIME_for-1);
    end
    Newton_Vg_save(TIME_for,1)=Newton_Vg;

    %% Current %%
    I_n_Sg(count_Sg_current,1)=0; I_p_Sg(count_Sg_current,1)=0; I_d_Sg(count_Sg_current,1)=0;
    for ii=1:size(Drain_BC,1)
        vertex_current=Drain_BC(ii,1);
        [row,col]=find(Element_si_drain==vertex_current);
        for k=1:size(row,1)

            index=zeros(9,1);
            n=1;
            for j=1:3
                index(n,1)=find(Table_Jaco(:,1)==4 & Table_Jaco(:,2)==Element_si_drain(row(k,1),j) & Table_Jaco(:,3)==1);
                index(n+1,1)=index(n,1)+1;
                index(n+2,1)=index(n,1)+2;
                n=n+3;
            end

            % potential, elec, hole 계산 %
            %%% DC %%%
            % potential
            x1_dc=solution_vector_Vg(index(1,1),1);
            x2_dc=solution_vector_Vg(index(4,1),1);
            x3_dc=solution_vector_Vg(index(7,1),1);
            x12_dc=x1_dc-x2_dc; x21_dc=-x12_dc;
            x23_dc=x2_dc-x3_dc; x32_dc=-x23_dc;
            x31_dc=x3_dc-x1_dc; x13_dc=-x31_dc;
            x12_dc_avr=(x1_dc+x2_dc)/2; x23_dc_avr=(x2_dc+x3_dc)/2; x31_dc_avr=(x3_dc+x1_dc)/2;
            x21_dc_avr=x12_dc_avr; x32_dc_avr=x23_dc_avr; x13_dc_avr=x31_dc_avr;

            % electron
            n1_dc=solution_vector_Vg(index(1,1)+1,1);
            n2_dc=solution_vector_Vg(index(4,1)+1,1);
            n3_dc=solution_vector_Vg(index(7,1)+1,1);
            n12_dc=n1_dc-n2_dc; n21_dc=-n12_dc;
            n23_dc=n2_dc-n3_dc; n32_dc=-n23_dc;
            n31_dc=n3_dc-n1_dc; n13_dc=-n31_dc;
            n12_dc_avr=(n1_dc+n2_dc)/2; n23_dc_avr=(n2_dc+n3_dc)/2; n31_dc_avr=(n3_dc+n1_dc)/2;
            n21_dc_avr=n12_dc_avr; n32_dc_avr=n23_dc_avr; n13_dc_avr=n31_dc_avr;

            % hole
            p1_dc=solution_vector_Vg(index(1,1)+2,1);
            p2_dc=solution_vector_Vg(index(4,1)+2,1);
            p3_dc=solution_vector_Vg(index(7,1)+2,1);
            p12_dc=p1_dc-p2_dc; p21_dc=-p12_dc;
            p23_dc= p2_dc-p3_dc; p32_dc=-p23_dc;
            p31_dc= p3_dc-p1_dc; p13_dc=-p31_dc;
            p12_dc_avr=(p1_dc+p2_dc)/2; p23_dc_avr=(p2_dc+p3_dc)/2; p31_dc_avr=(p3_dc+p1_dc)/2;
            p21_dc_avr=p12_dc_avr; p32_dc_avr=p23_dc_avr; p13_dc_avr=p31_dc_avr;

            %%% small signal %%%
            % potential
            x1_delta=solution_vector_Sg(index(1,1),1);
            x2_delta=solution_vector_Sg(index(4,1),1);
            x3_delta=solution_vector_Sg(index(7,1),1);
            x12_delta=x1_delta-x2_delta; x21_delta=-x12_delta;
            x23_delta=x2_delta-x3_delta; x32_delta=-x23_delta;
            x31_delta=x3_delta-x1_delta; x13_delta=-x31_delta;
            x12_delta_avr=(x1_delta+x2_delta)/2; x23_delta_avr=(x2_delta+x3_delta)/2; x31_delta_avr=(x3_delta+x1_delta)/2;
            x21_delta_avr=x12_delta_avr; x32_delta_avr=x23_delta_avr; x13_delta_avr=x31_delta_avr;

            % electron
            n1_delta=solution_vector_Sg(index(1,1)+1,1);
            n2_delta=solution_vector_Sg(index(4,1)+1,1);
            n3_delta=solution_vector_Sg(index(7,1)+1,1);
            n12_delta=n1_delta-n2_delta; n21_delta=-n12_delta;
            n23_delta=n2_delta-n3_delta; n32_delta=-n23_delta;
            n31_delta=n3_delta-n1_delta; n13_delta=-n31_delta;
            n12_delta_avr=(n1_delta+n2_delta)/2; n23_delta_avr=(n2_delta+n3_delta)/2; n31_delta_avr=(n3_delta+n1_delta)/2;
            n21_delta_avr=n12_delta_avr; n32_delta_avr=n23_delta_avr; n13_delta_avr=n31_delta_avr;

            % hole
            p1_delta=solution_vector_Sg(index(1,1)+2,1);
            p2_delta=solution_vector_Sg(index(4,1)+2,1);
            p3_delta=solution_vector_Sg(index(7,1)+2,1);
            p12_delta=p1_delta-p2_delta; p21_delta=-p12_delta;
            p23_delta= p2_delta-p3_delta; p32_delta=-p23_delta;
            p31_delta= p3_delta-p1_delta; p13_delta=-p31_delta;
            p12_delta_avr=(p1_delta+p2_delta)/2; p23_delta_avr=(p2_delta+p3_delta)/2; p31_delta_avr=(p3_delta+p1_delta)/2;
            p21_delta_avr=p12_delta_avr; p32_delta_avr=p23_delta_avr; p13_delta_avr=p31_delta_avr;

            %%% old %%%
            % potential
            x1_delta_old=solution_vector_Sg_old(index(1,1),1);
            x2_delta_old=solution_vector_Sg_old(index(4,1),1);
            x3_delta_old=solution_vector_Sg_old(index(7,1),1);
            x12_delta_old=x1_delta_old-x2_delta_old; x21_delta_old=-x12_delta_old;
            x23_delta_old=x2_delta_old-x3_delta_old; x32_delta_old=-x23_delta_old;
            x31_delta_old=x3_delta_old-x1_delta_old; x13_delta_old=-x31_delta_old;
            x12_delta_old_avr=(x1_delta_old+x2_delta_old)/2; x23_delta_old_avr=(x2_delta_old+x3_delta_old)/2; x31_delta_old_avr=(x3_delta_old+x1_delta_old)/2;
            x21_delta_old_avr=x12_delta_old_avr; x32_delta_old_avr=x23_delta_old_avr; x13_delta_old_avr=x31_delta_old_avr;

            % electron
            n1_delta_old=solution_vector_Sg_old(index(1,1)+1,1);
            n2_delta_old=solution_vector_Sg_old(index(4,1)+1,1);
            n3_delta_old=solution_vector_Sg_old(index(7,1)+1,1);
            n12_delta_old=n1_delta_old-n2_delta_old; n21_delta_old=-n12_delta_old;
            n23_delta_old=n2_delta_old-n3_delta_old; n32_delta_old=-n23_delta_old;
            n31_delta_old=n3_delta_old-n1_delta_old; n13_delta_old=-n31_delta_old;
            n12_delta_old_avr=(n1_delta_old+n2_delta_old)/2; n23_delta_old_avr=(n2_delta_old+n3_delta_old)/2; n31_delta_old_avr=(n3_delta_old+n1_delta_old)/2;
            n21_delta_old_avr=n12_delta_old_avr; n32_delta_old_avr=n23_delta_old_avr; n13_delta_old_avr=n31_delta_old_avr;

            % hole
            p1_delta_old=solution_vector_Sg_old(index(1,1)+2,1);
            p2_delta_old=solution_vector_Sg_old(index(4,1)+2,1);
            p3_delta_old=solution_vector_Sg_old(index(7,1)+2,1);
            p12_delta_old=p1_delta_old-p2_delta_old; p21_delta_old=-p12_delta_old;
            p23_delta_old= p2_delta_old-p3_delta_old; p32_delta_old=-p23_delta_old;
            p31_delta_old= p3_delta_old-p1_delta_old; p13_delta_old=-p31_delta_old;
            p12_delta_old_avr=(p1_delta_old+p2_delta_old)/2; p23_delta_old_avr=(p2_delta_old+p3_delta_old)/2; p31_delta_old_avr=(p3_delta_old+p1_delta_old)/2;
            p21_delta_old_avr=p12_delta_old_avr; p32_delta_old_avr=p23_delta_old_avr; p13_delta_old_avr=p31_delta_old_avr;


            index_element=find(Table_element_region(:,1)==4 & Table_element_region(:,2)==row(k,1));

            % Control_Volume 계산
            Control_Volume=zeros(3,1);
            Control_Volume(1,1)=(0.25)*L(index_element,3)*edge(index_element,3)+(0.25)*L(index_element,1)*edge(index_element,1);
            Control_Volume(2,1)=(0.25)*L(index_element,1)*edge(index_element,1)+(0.25)*L(index_element,2)*edge(index_element,2);
            Control_Volume(3,1)=(0.25)*L(index_element,2)*edge(index_element,2)+(0.25)*L(index_element,3)*edge(index_element,3);

            eps_si=11.7; I_n_Sg_tmp=0; I_p_Sg_tmp=0; I_d_tmp=0;
            if col(k,1)==1
                I_n_Sg_tmp=I_n_Sg_tmp+coeff_dJn*((edge(ii,1)/L(ii,1))*((n21_dc_avr*x21_delta)+(n21_delta_avr*x21_dc)-V_T*n21_delta));
                I_n_Sg_tmp=I_n_Sg_tmp+coeff_dJn*((edge(ii,3)/L(ii,3))*((n31_dc_avr*x31_delta)+(n31_delta_avr*x31_dc)-V_T*n31_delta));
                I_p_Sg_tmp=I_p_Sg_tmp+coeff_dJp*((edge(ii,1)/L(ii,1))*((p21_dc_avr*x21_delta)+(p21_delta_avr*x21_dc)+V_T*p21_delta));
                I_p_Sg_tmp=I_p_Sg_tmp+coeff_dJp*((edge(ii,3)/L(ii,3))*((p31_dc_avr*x31_delta)+(p31_delta_avr*x31_dc)+V_T*p31_delta));
                I_d_tmp=I_d_tmp-eps_si*eps0*(x21_delta-x21_delta_old)/delta_t/L(index_element,1)*edge(index_element,1);
                I_d_tmp=I_d_tmp-eps_si*eps0*(x31_delta-x3_delta1_old)/delta_t/L(index_element,3)*edge(index_element,3);
            elseif col(k,1)==2
                I_n_Sg_tmp=I_n_Sg_tmp+coeff_dJn*((edge(index_element,2)/L(index_element,2))*((n32_dc_avr*x32_delta)+(n32_delta_avr*x32_dc)-V_T*n32_delta));
                I_n_Sg_tmp=I_n_Sg_tmp+coeff_dJn*((edge(index_element,1)/L(index_element,1))*((n12_dc_avr*x12_delta)+(n12_delta_avr*x12_dc)-V_T*n12_delta));
                I_p_Sg_tmp=I_p_Sg_tmp+coeff_dJp*((edge(index_element,2)/L(index_element,2))*((p32_dc_avr*x32_delta)+(p32_delta_avr*x32_dc)+V_T*p32_delta));
                I_p_Sg_tmp=I_p_Sg_tmp+coeff_dJp*((edge(index_element,1)/L(index_element,1))*((p12_dc_avr*x12_delta)+(p12_delta_avr*x12_dc)+V_T*p12_delta));
                I_d_tmp=I_d_tmp-eps_si*eps0*(x32_delta-x32_delta_old)/delta_t/L(index_element,2)*edge(index_element,2);
                I_d_tmp=I_d_tmp-eps_si*eps0*(x12_delta-x12_delta_old)/delta_t/L(index_element,1)*edge(index_element,1);
            elseif col(k,1)==3
                I_n_Sg_tmp=I_n_Sg_tmp+coeff_dJn*((edge(index_element,3)/L(index_element,3))*((n13_dc_avr*x13_delta)+(n13_delta_avr*x13_dc)-V_T*n13_delta));
                I_n_Sg_tmp=I_n_Sg_tmp+coeff_dJn*((edge(index_element,2)/L(index_element,2))*((n23_dc_avr*x23_delta)+(n23_delta_avr*x23_dc)-V_T*n23_delta));
                I_p_Sg_tmp=I_p_Sg_tmp+coeff_dJp*((edge(index_element,3)/L(index_element,3))*((p13_dc_avr*x13_delta)+(p13_delta_avr*x13_dc)+V_T*p13_delta));
                I_p_Sg_tmp=I_p_Sg_tmp+coeff_dJp*((edge(index_element,2)/L(index_element,2))*((p23_dc_avr*x23_delta)+(p23_delta_avr*x23_dc)+V_T*p23_delta));
                I_d_tmp=I_d_tmp-eps_si*eps0*(x13_delta-x13_delta_old)/delta_t/L(index_element,3)*edge(index_element,3);
                I_d_tmp=I_d_tmp-eps_si*eps0*(x23_delta-x23_delta_old)/delta_t/L(index_element,2)*edge(index_element,2);
            end
            I_n_Sg(count_Sg_current,1)=I_n_Sg(count_Sg_current,1)+I_n_Sg_tmp*width;
            I_p_Sg(count_Sg_current,1)=I_p_Sg(count_Sg_current,1)+I_p_Sg_tmp*width;
            I_d_Sg(count_Sg_current,1)=I_d_Sg(count_Sg_current,1)+I_d_tmp*width;


        end
    end
    I_Sg(count_Sg_current,1) = (I_n_Sg(count_Sg_current,1)+I_p_Sg(count_Sg_current,1)+I_d_Sg(count_Sg_current,1));
    count_Sg_current=count_Sg_current+1;
end