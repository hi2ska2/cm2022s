clear; clc; close all
load('Ramping_0.5V_data.mat')

%% AC simulation

freq=logspace(10,14, 20)';
for freq_for=1:length(freq)

    clearvars -except freq freq_for I_AC Admittance Admittance_real Admittance_imag
    load('Ramping_0.5V_data.mat')

    coeff=q/eps0;
    coeff_dJn=-q*mobility_n;
    coeff_dJp=-q*mobility_p;


    f=freq(freq_for,1);
    Amp=1e-3;
    width=1e-6; % m
    w=2*pi*f;

    % initial value
    solution_vector_AC=zeros(length(solution_vector_Vd),1);

    % 여기서 나오는 Jaco, res는 사실 AC에서는 A matrix, b vector이다.(oneshot이므로)
    Jaco_AC=zeros(size(solution_vector_AC,1),size(solution_vector_AC,1)); res_AC=zeros(size(solution_vector_AC,1),1);
    % Jaco/res 작성 %
    for ii=1:length(Element)
        %% 각종 변수 사전 계산
        % index : 현재 element에 속해있는 vertex가 resdue에 몇번째 있는가 저장
        index=zeros(9,1);
        n=1;
        for j=1:3
            index(n,1)=find(Table_Jaco(:,1)==Table_element_region(ii,1) & Table_Jaco(:,2)==Element(ii,j) & Table_Jaco(:,3)==1);
            index(n+1,1)=index(n,1)+1;
            index(n+2,1)=index(n,1)+2;
            n=n+3;
        end

        % Control_Volume 계산
        Control_Volume=zeros(3,1);
        Control_Volume(1,1)=(0.25)*L(ii,3)*edge(ii,3)+(0.25)*L(ii,1)*edge(ii,1);
        Control_Volume(2,1)=(0.25)*L(ii,1)*edge(ii,1)+(0.25)*L(ii,2)*edge(ii,2);
        Control_Volume(3,1)=(0.25)*L(ii,2)*edge(ii,2)+(0.25)*L(ii,3)*edge(ii,3);

        % potential, elec, hole 계산 %
        %%% DC %%%
        % potential
        x1_dc=solution_vector_Vd(index(1,1),1);
        x2_dc=solution_vector_Vd(index(4,1),1);
        x3_dc=solution_vector_Vd(index(7,1),1);
        x12_dc=x1_dc-x2_dc; x21_dc=-x12_dc;
        x23_dc=x2_dc-x3_dc; x32_dc=-x23_dc;
        x31_dc=x3_dc-x1_dc; x13_dc=-x31_dc;
        x12_dc_avr=(x1_dc+x2_dc)/2; x23_dc_avr=(x2_dc+x3_dc)/2; x31_dc_avr=(x3_dc+x1_dc)/2;
        x21_dc_avr=x12_dc_avr; x32_dc_avr=x23_dc_avr; x13_dc_avr=x31_dc_avr;

        % electron
        n1_dc=solution_vector_Vd(index(1,1)+1,1);
        n2_dc=solution_vector_Vd(index(4,1)+1,1);
        n3_dc=solution_vector_Vd(index(7,1)+1,1);
        n12_dc=n1_dc-n2_dc; n21_dc=-n12_dc;
        n23_dc=n2_dc-n3_dc; n32_dc=-n23_dc;
        n31_dc=n3_dc-n1_dc; n13_dc=-n31_dc;
        n12_dc_avr=(n1_dc+n2_dc)/2; n23_dc_avr=(n2_dc+n3_dc)/2; n31_dc_avr=(n3_dc+n1_dc)/2;
        n21_dc_avr=n12_dc_avr; n32_dc_avr=n23_dc_avr; n13_dc_avr=n31_dc_avr;

        % hole
        p1_dc=solution_vector_Vd(index(1,1)+2,1);
        p2_dc=solution_vector_Vd(index(4,1)+2,1);
        p3_dc=solution_vector_Vd(index(7,1)+2,1);
        p12_dc=p1_dc-p2_dc; p21_dc=-p12_dc;
        p23_dc= p2_dc-p3_dc; p32_dc=-p23_dc;
        p31_dc= p3_dc-p1_dc; p13_dc=-p31_dc;
        p12_dc_avr=(p1_dc+p2_dc)/2; p23_dc_avr=(p2_dc+p3_dc)/2; p31_dc_avr=(p3_dc+p1_dc)/2;
        p21_dc_avr=p12_dc_avr; p32_dc_avr=p23_dc_avr; p13_dc_avr=p31_dc_avr;

        %% res vector
        res_tmp_AC=zeros(9,1); % oneshot에 푸는것임. 사실 res가 아니라 b vector.

        %% Jacobian matrix
        Jaco_tmp_AC=zeros(9,9);

        % node1 - potential
        Jaco_tmp_AC(1,1)=eps(Table_element_region(ii,1),1)*(-edge(ii,1)/L(ii,1)-edge(ii,3)/L(ii,3));
        Jaco_tmp_AC(1,2)=-coeff*Control_Volume(1,1);
        Jaco_tmp_AC(1,3)=coeff*Control_Volume(1,1);
        Jaco_tmp_AC(1,4)=eps(Table_element_region(ii,1),1)*edge(ii,1)/L(ii,1);
        Jaco_tmp_AC(1,7)=eps(Table_element_region(ii,1),1)*edge(ii,3)/L(ii,3);

        % node1 - electron
        Jaco_tmp_AC(2,1)=coeff_dJn*((edge(ii,1)/L(ii,1))*(-n21_dc_avr)+(edge(ii,3)/L(ii,3))*(-n31_dc_avr));
        Jaco_tmp_AC(2,2)=coeff_dJn*((edge(ii,1)/L(ii,1))*(0.5*x21_dc+V_T)+(edge(ii,3)/L(ii,3))*(0.5*x31_dc+V_T));
        Jaco_tmp_AC(2,4)=coeff_dJn*((edge(ii,1)/L(ii,1))*(n21_dc_avr));
        Jaco_tmp_AC(2,5)=coeff_dJn*((edge(ii,1)/L(ii,1))*(0.5*x21_dc-V_T));
        Jaco_tmp_AC(2,7)=coeff_dJn*((edge(ii,3)/L(ii,3))*(n31_dc_avr));
        Jaco_tmp_AC(2,8)=coeff_dJn*((edge(ii,3)/L(ii,3))*(0.5*x31_dc-V_T));
        % node1 - electron_transient
        Jaco_tmp_AC(2,2)=Jaco_tmp_AC(2,2)-q*Control_Volume(1,1)*1i*w;

        % node1 - hole
        Jaco_tmp_AC(3,1)=coeff_dJp*((edge(ii,1)/L(ii,1))*(-n21_dc_avr)+(edge(ii,3)/L(ii,3))*(-n31_dc_avr));
        Jaco_tmp_AC(3,3)=coeff_dJp*((edge(ii,1)/L(ii,1))*(0.5*x21_dc-V_T)+(edge(ii,3)/L(ii,3))*(0.5*x31_dc-V_T));
        Jaco_tmp_AC(3,4)=coeff_dJp*((edge(ii,1)/L(ii,1))*(n21_dc_avr));
        Jaco_tmp_AC(3,6)=coeff_dJp*((edge(ii,1)/L(ii,1))*(0.5*x21_dc+V_T));
        Jaco_tmp_AC(3,7)=coeff_dJp*((edge(ii,3)/L(ii,3))*(n31_dc_avr));
        Jaco_tmp_AC(3,9)=coeff_dJp*((edge(ii,3)/L(ii,3))*(0.5*x31_dc+V_T));
        % node1 - hole_transient
        Jaco_tmp_AC(3,3)=Jaco_tmp_AC(3,3)+q*Control_Volume(1,1)*1i*w;

        % node2 - potential
        Jaco_tmp_AC(4,1)=eps(Table_element_region(ii,1),1)*edge(ii,1)/L(ii,1);
        Jaco_tmp_AC(4,4)=eps(Table_element_region(ii,1),1)*(-edge(ii,1)/L(ii,1)-edge(ii,2)/L(ii,2));
        Jaco_tmp_AC(4,5)=-coeff*Control_Volume(2,1);
        Jaco_tmp_AC(4,6)=coeff*Control_Volume(2,1);
        Jaco_tmp_AC(4,7)=eps(Table_element_region(ii,1),1)*edge(ii,2)/L(ii,2);

        % node2 - electron
        Jaco_tmp_AC(5,1)=coeff_dJn*((edge(ii,1)/L(ii,1))*(n12_dc_avr));
        Jaco_tmp_AC(5,2)=coeff_dJn*((edge(ii,1)/L(ii,1))*(0.5*x12_dc-V_T));
        Jaco_tmp_AC(5,4)=coeff_dJn*((edge(ii,2)/L(ii,2))*(-n32_dc_avr)+(edge(ii,1)/L(ii,1))*(-n12_dc_avr));
        Jaco_tmp_AC(5,5)=coeff_dJn*((edge(ii,2)/L(ii,2))*(0.5*x32_dc+V_T)+(edge(ii,1)/L(ii,1))*(0.5*x12_dc+V_T));
        Jaco_tmp_AC(5,7)=coeff_dJn*((edge(ii,2)/L(ii,2))*(n32_dc_avr));
        Jaco_tmp_AC(5,8)=coeff_dJn*((edge(ii,2)/L(ii,2))*(0.5*x32_dc-V_T));
        % node2 - electron_transient
        Jaco_tmp_AC(5,5)=Jaco_tmp_AC(5,5)-q*Control_Volume(2,1)*1i*w;

        % node2 - hole
        Jaco_tmp_AC(6,1)=coeff_dJp*((edge(ii,1)/L(ii,1))*(n12_dc_avr));
        Jaco_tmp_AC(6,3)=coeff_dJp*((edge(ii,1)/L(ii,1))*(0.5*x12_dc+V_T));
        Jaco_tmp_AC(6,4)=coeff_dJp*((edge(ii,2)/L(ii,2))*(-n32_dc_avr)+(edge(ii,1)/L(ii,1))*(-n12_dc_avr));
        Jaco_tmp_AC(6,6)=coeff_dJp*((edge(ii,2)/L(ii,2))*(0.5*x32_dc-V_T)+(edge(ii,1)/L(ii,1))*(0.5*x12_dc-V_T));
        Jaco_tmp_AC(6,7)=coeff_dJp*((edge(ii,2)/L(ii,2))*(n32_dc_avr));
        Jaco_tmp_AC(6,9)=coeff_dJp*((edge(ii,2)/L(ii,2))*(0.5*x32_dc+V_T));
        % node2 - hole_transient
        Jaco_tmp_AC(6,6)=Jaco_tmp_AC(6,6)+q*Control_Volume(2,1)*1i*w;

        % node3 - potential
        Jaco_tmp_AC(7,1)=eps(Table_element_region(ii,1),1)*edge(ii,3)/L(ii,3);
        Jaco_tmp_AC(7,4)=eps(Table_element_region(ii,1),1)*edge(ii,2)/L(ii,2);
        Jaco_tmp_AC(7,7)=eps(Table_element_region(ii,1),1)*(-edge(ii,2)/L(ii,2)-edge(ii,3)/L(ii,3));
        Jaco_tmp_AC(7,8)=-coeff*Control_Volume(3,1);
        Jaco_tmp_AC(7,9)=coeff*Control_Volume(3,1);

        % node3 - electron
        Jaco_tmp_AC(8,1)=coeff_dJn*((edge(ii,3)/L(ii,3))*(n13_dc_avr));
        Jaco_tmp_AC(8,2)=coeff_dJn*((edge(ii,3)/L(ii,3))*(0.5*x13_dc-V_T));
        Jaco_tmp_AC(8,4)=coeff_dJn*((edge(ii,2)/L(ii,2))*(n23_dc_avr));
        Jaco_tmp_AC(8,5)=coeff_dJn*((edge(ii,2)/L(ii,2))*(0.5*x23_dc-V_T));
        Jaco_tmp_AC(8,7)=coeff_dJn*((edge(ii,3)/L(ii,3))*(-n13_dc_avr)+(edge(ii,2)/L(ii,2))*(-n23_dc_avr));
        Jaco_tmp_AC(8,8)=coeff_dJn*((edge(ii,3)/L(ii,3))*(0.5*x13_dc+V_T)+(edge(ii,2)/L(ii,2))*(0.5*x23_dc+V_T));
        % node3 - electron_transient
        Jaco_tmp_AC(8,8)=Jaco_tmp_AC(8,8)-q*Control_Volume(3,1)*1i*w;

        % node3 - hole
        Jaco_tmp_AC(9,1)=coeff_dJp*((edge(ii,3)/L(ii,3))*(n13_dc_avr));
        Jaco_tmp_AC(9,3)=coeff_dJp*((edge(ii,3)/L(ii,3))*(0.5*x13_dc+V_T));
        Jaco_tmp_AC(9,4)=coeff_dJp*((edge(ii,2)/L(ii,2))*(n23_dc_avr));
        Jaco_tmp_AC(9,6)=coeff_dJp*((edge(ii,2)/L(ii,2))*(0.5*x23_dc+V_T));
        Jaco_tmp_AC(9,7)=coeff_dJp*((edge(ii,3)/L(ii,3))*(-n13_dc_avr)+(edge(ii,2)/L(ii,2))*(-n23_dc_avr));
        Jaco_tmp_AC(9,9)=coeff_dJp*((edge(ii,3)/L(ii,3))*(0.5*x13_dc-V_T)+(edge(ii,2)/L(ii,2))*(0.5*x23_dc-V_T));
        % node3 - hole_transient
        Jaco_tmp_AC(9,9)=Jaco_tmp_AC(9,9)+q*Control_Volume(3,1)*1i*w;


        for j=1:9
            for k=1:9
                Jaco_AC(index(j,1),index(k,1))=Jaco_AC(index(j,1),index(k,1))+Jaco_tmp_AC(j,k);
            end
        end
    end

    % Boundary condition %
    % cathode_BC
    for ii=1:size(anode_BC,1)
        BC_index=find(Table_Jaco(:,1)==1 & Table_Jaco(:,2)==anode_BC(ii,1) & Table_Jaco(:,3)==1);
        Jaco_AC(BC_index,:)=0; Jaco_AC(BC_index+1,:)=0; Jaco_AC(BC_index+2,:)=0;
        Jaco_AC(BC_index,BC_index)=1; Jaco_AC(BC_index+1,BC_index+1)=1; Jaco_AC(BC_index+2,BC_index+2)=1;
        res_AC(BC_index,1)=solution_vector_AC(BC_index,1)-0;
        res_AC(BC_index+1,1)=solution_vector_AC(BC_index+1,1)-0;
        res_AC(BC_index+2,1)=solution_vector_AC(BC_index+2,1)-0;
    end

    % cathode_BC
    for ii=1:size(cathode_BC,1)
        BC_index=find(Table_Jaco(:,1)==1 & Table_Jaco(:,2)==cathode_BC(ii,1) & Table_Jaco(:,3)==1);
        Jaco_AC(BC_index,:)=0; Jaco_AC(BC_index+1,:)=0; Jaco_AC(BC_index+2,:)=0;
        Jaco_AC(BC_index,BC_index)=1; Jaco_AC(BC_index+1,BC_index+1)=1; Jaco_AC(BC_index+2,BC_index+2)=1;
        res_AC(BC_index,1)=solution_vector_AC(BC_index,1)-1i*Amp;
        res_AC(BC_index+1,1)=solution_vector_AC(BC_index+1,1)-0;
        res_AC(BC_index+2,1)=solution_vector_AC(BC_index+2,1)-0;
    end

    % Scaling %
    Cvector_AC=zeros(size(solution_vector_AC,1),1); v=[V_T; max(abs(N_dop)); max(abs(N_dop))];
    for ii=1:size(solution_vector_AC,1)
        Cvector_AC(ii,1)=v(Table_Jaco(ii,3),1);
    end
    Jaco_AC=sparse(Jaco_AC);
    Cmatrix_AC=spdiags(Cvector_AC,0,size(solution_vector_AC,1),size(solution_vector_AC,1));
    Jaco_scaled_AC=Jaco_AC*Cmatrix_AC;
    Rvector_AC=1./sum(abs(Jaco_scaled_AC),2);
    Rmatrix_AC=spdiags(Rvector_AC,0,size(solution_vector_AC,1),size(solution_vector_AC,1));
    Jaco_scaled_AC=Rmatrix_AC*Jaco_scaled_AC;
    res_scaled_AC=Rmatrix_AC*res_AC;
    update_scaled_AC=Jaco_scaled_AC\(-res_scaled_AC);
    update_AC(:,1)=Cmatrix_AC*update_scaled_AC;

    solution_vector_AC(:,1)=update_AC(:,1);


    %% Current %%
    I_n_AC(1,1)=0; I_p_AC(1,1)=0; I_d_AC(1,1)=0;
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

            % potential, elec, hole 계산 %
            %%% DC %%%
            % potential
            x1_dc=solution_vector_Vd(index(1,1),1);
            x2_dc=solution_vector_Vd(index(4,1),1);
            x3_dc=solution_vector_Vd(index(7,1),1);
            x12_dc=x1_dc-x2_dc; x21_dc=-x12_dc;
            x23_dc=x2_dc-x3_dc; x32_dc=-x23_dc;
            x31_dc=x3_dc-x1_dc; x13_dc=-x31_dc;
            x12_dc_avr=(x1_dc+x2_dc)/2; x23_dc_avr=(x2_dc+x3_dc)/2; x31_dc_avr=(x3_dc+x1_dc)/2;
            x21_dc_avr=x12_dc_avr; x32_dc_avr=x23_dc_avr; x13_dc_avr=x31_dc_avr;

            % electron
            n1_dc=solution_vector_Vd(index(1,1)+1,1);
            n2_dc=solution_vector_Vd(index(4,1)+1,1);
            n3_dc=solution_vector_Vd(index(7,1)+1,1);
            n12_dc=n1_dc-n2_dc; n21_dc=-n12_dc;
            n23_dc=n2_dc-n3_dc; n32_dc=-n23_dc;
            n31_dc=n3_dc-n1_dc; n13_dc=-n31_dc;
            n12_dc_avr=(n1_dc+n2_dc)/2; n23_dc_avr=(n2_dc+n3_dc)/2; n31_dc_avr=(n3_dc+n1_dc)/2;
            n21_dc_avr=n12_dc_avr; n32_dc_avr=n23_dc_avr; n13_dc_avr=n31_dc_avr;

            % hole
            p1_dc=solution_vector_Vd(index(1,1)+2,1);
            p2_dc=solution_vector_Vd(index(4,1)+2,1);
            p3_dc=solution_vector_Vd(index(7,1)+2,1);
            p12_dc=p1_dc-p2_dc; p21_dc=-p12_dc;
            p23_dc= p2_dc-p3_dc; p32_dc=-p23_dc;
            p31_dc= p3_dc-p1_dc; p13_dc=-p31_dc;
            p12_dc_avr=(p1_dc+p2_dc)/2; p23_dc_avr=(p2_dc+p3_dc)/2; p31_dc_avr=(p3_dc+p1_dc)/2;
            p21_dc_avr=p12_dc_avr; p32_dc_avr=p23_dc_avr; p13_dc_avr=p31_dc_avr;

            %%% small signal %%%
            % potential
            x1_delta=solution_vector_AC(index(1,1),1);
            x2_delta=solution_vector_AC(index(4,1),1);
            x3_delta=solution_vector_AC(index(7,1),1);
            x12_delta=x1_delta-x2_delta; x21_delta=-x12_delta;
            x23_delta=x2_delta-x3_delta; x32_delta=-x23_delta;
            x31_delta=x3_delta-x1_delta; x13_delta=-x31_delta;
            x12_delta_avr=(x1_delta+x2_delta)/2; x23_delta_avr=(x2_delta+x3_delta)/2; x31_delta_avr=(x3_delta+x1_delta)/2;
            x21_delta_avr=x12_delta_avr; x32_delta_avr=x23_delta_avr; x13_delta_avr=x31_delta_avr;

            % electron
            n1_delta=solution_vector_AC(index(1,1)+1,1);
            n2_delta=solution_vector_AC(index(4,1)+1,1);
            n3_delta=solution_vector_AC(index(7,1)+1,1);
            n12_delta=n1_delta-n2_delta; n21_delta=-n12_delta;
            n23_delta=n2_delta-n3_delta; n32_delta=-n23_delta;
            n31_delta=n3_delta-n1_delta; n13_delta=-n31_delta;
            n12_delta_avr=(n1_delta+n2_delta)/2; n23_delta_avr=(n2_delta+n3_delta)/2; n31_delta_avr=(n3_delta+n1_delta)/2;
            n21_delta_avr=n12_delta_avr; n32_delta_avr=n23_delta_avr; n13_delta_avr=n31_delta_avr;

            % hole
            p1_delta=solution_vector_AC(index(1,1)+2,1);
            p2_delta=solution_vector_AC(index(4,1)+2,1);
            p3_delta=solution_vector_AC(index(7,1)+2,1);
            p12_delta=p1_delta-p2_delta; p21_delta=-p12_delta;
            p23_delta= p2_delta-p3_delta; p32_delta=-p23_delta;
            p31_delta= p3_delta-p1_delta; p13_delta=-p31_delta;
            p12_delta_avr=(p1_delta+p2_delta)/2; p23_delta_avr=(p2_delta+p3_delta)/2; p31_delta_avr=(p3_delta+p1_delta)/2;
            p21_delta_avr=p12_delta_avr; p32_delta_avr=p23_delta_avr; p13_delta_avr=p31_delta_avr;


            index_element=find(Table_element_region(:,1)==1 & Table_element_region(:,2)==row(k,1));
            eps_si=11.7; I_n_AC_tmp=0; I_p_AC_tmp=0; I_d_tmp=0;
            if col(k,1)==1
                I_n_AC_tmp=I_n_AC_tmp+coeff_dJn*((edge(ii,1)/L(ii,1))*((n21_dc_avr*x21_delta)+(n21_delta_avr*x21_dc)-V_T*n21_delta));
                I_n_AC_tmp=I_n_AC_tmp+coeff_dJn*((edge(ii,3)/L(ii,3))*((n31_dc_avr*x31_delta)+(n31_delta_avr*x31_dc)-V_T*n31_delta));
                I_p_AC_tmp=I_p_AC_tmp+coeff_dJp*((edge(ii,1)/L(ii,1))*((p21_dc_avr*x21_delta)+(p21_delta_avr*x21_dc)+V_T*p21_delta));
                I_p_AC_tmp=I_p_AC_tmp+coeff_dJp*((edge(ii,3)/L(ii,3))*((p31_dc_avr*x31_delta)+(p31_delta_avr*x31_dc)+V_T*p31_delta));
                I_d_tmp=I_d_tmp-eps_si*eps0*(x21_delta)*1i*w/L(index_element,1)*edge(index_element,1);
                I_d_tmp=I_d_tmp-eps_si*eps0*(x31_delta)*1i*w/L(index_element,3)*edge(index_element,3);
            elseif col(k,1)==2
                I_n_AC_tmp=I_n_AC_tmp+coeff_dJn*((edge(index_element,2)/L(index_element,2))*((n32_dc_avr*x32_delta)+(n32_delta_avr*x32_dc)-V_T*n32_delta));
                I_n_AC_tmp=I_n_AC_tmp+coeff_dJn*((edge(index_element,1)/L(index_element,1))*((n12_dc_avr*x12_delta)+(n12_delta_avr*x12_dc)-V_T*n12_delta));
                I_p_AC_tmp=I_p_AC_tmp+coeff_dJp*((edge(index_element,2)/L(index_element,2))*((p32_dc_avr*x32_delta)+(p32_delta_avr*x32_dc)+V_T*p32_delta));
                I_p_AC_tmp=I_p_AC_tmp+coeff_dJp*((edge(index_element,1)/L(index_element,1))*((p12_dc_avr*x12_delta)+(p12_delta_avr*x12_dc)+V_T*p12_delta));
                I_d_tmp=I_d_tmp-eps_si*eps0*(x32_delta)*1i*w/L(index_element,2)*edge(index_element,2);
                I_d_tmp=I_d_tmp-eps_si*eps0*(x12_delta)*1i*w/L(index_element,1)*edge(index_element,1);
            elseif col(k,1)==3
                I_n_AC_tmp=I_n_AC_tmp+coeff_dJn*((edge(index_element,3)/L(index_element,3))*((n13_dc_avr*x13_delta)+(n13_delta_avr*x13_dc)-V_T*n13_delta));
                I_n_AC_tmp=I_n_AC_tmp+coeff_dJn*((edge(index_element,2)/L(index_element,2))*((n23_dc_avr*x23_delta)+(n23_delta_avr*x23_dc)-V_T*n23_delta));
                I_p_AC_tmp=I_p_AC_tmp+coeff_dJp*((edge(index_element,3)/L(index_element,3))*((p13_dc_avr*x13_delta)+(p13_delta_avr*x13_dc)+V_T*p13_delta));
                I_p_AC_tmp=I_p_AC_tmp+coeff_dJp*((edge(index_element,2)/L(index_element,2))*((p23_dc_avr*x23_delta)+(p23_delta_avr*x23_dc)+V_T*p23_delta));
                I_d_tmp=I_d_tmp-eps_si*eps0*(x13_delta)*1i*w/L(index_element,3)*edge(index_element,3);
                I_d_tmp=I_d_tmp-eps_si*eps0*(x23_delta)*1i*w/L(index_element,2)*edge(index_element,2);
            end
            I_n_AC(1,1)=I_n_AC(1,1)+I_n_AC_tmp*width;
            I_p_AC(1,1)=I_p_AC(1,1)+I_p_AC_tmp*width;
            I_d_AC(1,1)=I_d_AC(1,1)+I_d_tmp*width;


        end
    end
    I_AC(1,freq_for) = (I_n_AC(1,1)+I_p_AC(1,1)+I_d_AC(1,1));

    Admittance(freq_for,1)=I_AC(1,freq_for)/(Amp*1i);
    Admittance_real(freq_for,1)=real(Admittance(freq_for,1));
    Admittance_imag(freq_for,1)=imag(Admittance(freq_for,1));

end
