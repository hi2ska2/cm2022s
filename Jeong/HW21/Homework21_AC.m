clear; clc; close all
load('./data/homework21_dd.mat')

%% AC simulation

freq=1e9; %logspace(10,14, 20)';
for freq_for=1:length(freq)

    clearvars -except freq freq_for I_AC_cathode Admittance Admittance_real Admittance_imag
    load('./data/homework21_dd.mat')

    coeff=q/eps0;
    coeff_dJn=-q*mobility_n;
    coeff_dJp=-q*mobility_p;


    f=freq(freq_for,1);
    Amp=1e-3;
    width=1e-6; % m
    w=2*pi*f;

    % initial value
    solution_vector_AC=zeros(length(solution_vector_DD),length(solution_vector_DD));

    % 여기서 나오는 Jaco, res는 사실 AC에서는 A matrix, b vector이다.(oneshot이므로)
    Jaco_AC=zeros(size(solution_vector_AC,1),size(solution_vector_AC,1)); res_AC=zeros(size(solution_vector_AC,1),size(solution_vector_AC,1));
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
        x1_dc=solution_vector_DD(index(1,1),1);
        x2_dc=solution_vector_DD(index(4,1),1);
        x3_dc=solution_vector_DD(index(7,1),1);
        x12_dc=x1_dc-x2_dc; x21_dc=-x12_dc;
        x23_dc=x2_dc-x3_dc; x32_dc=-x23_dc;
        x31_dc=x3_dc-x1_dc; x13_dc=-x31_dc;
        x12_dc_avr=(x1_dc+x2_dc)/2; x23_dc_avr=(x2_dc+x3_dc)/2; x31_dc_avr=(x3_dc+x1_dc)/2;
        x21_dc_avr=x12_dc_avr; x32_dc_avr=x23_dc_avr; x13_dc_avr=x31_dc_avr;

        % electron
        n1_dc=solution_vector_DD(index(1,1)+1,1);
        n2_dc=solution_vector_DD(index(4,1)+1,1);
        n3_dc=solution_vector_DD(index(7,1)+1,1);
        n12_dc=n1_dc-n2_dc; n21_dc=-n12_dc;
        n23_dc=n2_dc-n3_dc; n32_dc=-n23_dc;
        n31_dc=n3_dc-n1_dc; n13_dc=-n31_dc;
        n12_dc_avr=(n1_dc+n2_dc)/2; n23_dc_avr=(n2_dc+n3_dc)/2; n31_dc_avr=(n3_dc+n1_dc)/2;
        n21_dc_avr=n12_dc_avr; n32_dc_avr=n23_dc_avr; n13_dc_avr=n31_dc_avr;

        % hole
        p1_dc=solution_vector_DD(index(1,1)+2,1);
        p2_dc=solution_vector_DD(index(4,1)+2,1);
        p3_dc=solution_vector_DD(index(7,1)+2,1);
        p12_dc=p1_dc-p2_dc; p21_dc=-p12_dc;
        p23_dc= p2_dc-p3_dc; p32_dc=-p23_dc;
        p31_dc= p3_dc-p1_dc; p13_dc=-p31_dc;
        p12_dc_avr=(p1_dc+p2_dc)/2; p23_dc_avr=(p2_dc+p3_dc)/2; p31_dc_avr=(p3_dc+p1_dc)/2;
        p21_dc_avr=p12_dc_avr; p32_dc_avr=p23_dc_avr; p13_dc_avr=p31_dc_avr;

        %% res vector
        res_AC=diag(ones(length(solution_vector_AC),1));

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
        Jaco_tmp_AC(3,1)=coeff_dJp*((edge(ii,1)/L(ii,1))*(-p21_dc_avr)+(edge(ii,3)/L(ii,3))*(-p31_dc_avr));
        Jaco_tmp_AC(3,3)=coeff_dJp*((edge(ii,1)/L(ii,1))*(0.5*x21_dc-V_T)+(edge(ii,3)/L(ii,3))*(0.5*x31_dc-V_T));
        Jaco_tmp_AC(3,4)=coeff_dJp*((edge(ii,1)/L(ii,1))*(p21_dc_avr));
        Jaco_tmp_AC(3,6)=coeff_dJp*((edge(ii,1)/L(ii,1))*(0.5*x21_dc+V_T));
        Jaco_tmp_AC(3,7)=coeff_dJp*((edge(ii,3)/L(ii,3))*(p31_dc_avr));
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
        Jaco_tmp_AC(6,1)=coeff_dJp*((edge(ii,1)/L(ii,1))*(p12_dc_avr));
        Jaco_tmp_AC(6,3)=coeff_dJp*((edge(ii,1)/L(ii,1))*(0.5*x12_dc+V_T));
        Jaco_tmp_AC(6,4)=coeff_dJp*((edge(ii,2)/L(ii,2))*(-p32_dc_avr)+(edge(ii,1)/L(ii,1))*(-p12_dc_avr));
        Jaco_tmp_AC(6,6)=coeff_dJp*((edge(ii,2)/L(ii,2))*(0.5*x32_dc-V_T)+(edge(ii,1)/L(ii,1))*(0.5*x12_dc-V_T));
        Jaco_tmp_AC(6,7)=coeff_dJp*((edge(ii,2)/L(ii,2))*(p32_dc_avr));
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
        Jaco_tmp_AC(9,1)=coeff_dJp*((edge(ii,3)/L(ii,3))*(p13_dc_avr));
        Jaco_tmp_AC(9,3)=coeff_dJp*((edge(ii,3)/L(ii,3))*(0.5*x13_dc+V_T));
        Jaco_tmp_AC(9,4)=coeff_dJp*((edge(ii,2)/L(ii,2))*(p23_dc_avr));
        Jaco_tmp_AC(9,6)=coeff_dJp*((edge(ii,2)/L(ii,2))*(0.5*x23_dc+V_T));
        Jaco_tmp_AC(9,7)=coeff_dJp*((edge(ii,3)/L(ii,3))*(-p13_dc_avr)+(edge(ii,2)/L(ii,2))*(-p23_dc_avr));
        Jaco_tmp_AC(9,9)=coeff_dJp*((edge(ii,3)/L(ii,3))*(0.5*x13_dc-V_T)+(edge(ii,2)/L(ii,2))*(0.5*x23_dc-V_T));
        % node3 - hole_transient
        Jaco_tmp_AC(9,9)=Jaco_tmp_AC(9,9)+q*Control_Volume(3,1)*1i*w;


        for j=1:9
            for k=1:9
                Jaco_AC(index(j,1),index(k,1))=Jaco_AC(index(j,1),index(k,1))+Jaco_tmp_AC(j,k);
            end
        end
    end

    %% Boundary condition %
    % cathode_BC
    for ii=1:size(anode_BC,1)
        BC_index=find(Table_Jaco(:,1)==1 & Table_Jaco(:,2)==anode_BC(ii,1) & Table_Jaco(:,3)==1);
        Jaco_AC(BC_index,:)=0; Jaco_AC(BC_index+1,:)=0; Jaco_AC(BC_index+2,:)=0;
        Jaco_AC(BC_index,BC_index)=1; Jaco_AC(BC_index+1,BC_index+1)=1; Jaco_AC(BC_index+2,BC_index+2)=1;
        res_AC(BC_index,1)=0;
        res_AC(BC_index+1,1)=0;
        res_AC(BC_index+2,1)=0;
    end

    % cathode_BC
    for ii=1:size(cathode_BC,1)
        BC_index=find(Table_Jaco(:,1)==1 & Table_Jaco(:,2)==cathode_BC(ii,1) & Table_Jaco(:,3)==1);
        Jaco_AC(BC_index,:)=0; Jaco_AC(BC_index+1,:)=0; Jaco_AC(BC_index+2,:)=0;
        Jaco_AC(BC_index,BC_index)=1; Jaco_AC(BC_index+1,BC_index+1)=1; Jaco_AC(BC_index+2,BC_index+2)=1;
        res_AC(BC_index,1)=0;
        res_AC(BC_index+1,1)=0;
        res_AC(BC_index+2,1)=0;
    end

    % Scaling %
    Cvector_AC=zeros(size(solution_vector_AC,1),1); v=[V_T; max(abs(N_dop)); max(abs(N_dop))];
    for ii=1:size(solution_vector_AC,1)
        Cvector_AC(ii,1)=v(Table_Jaco(ii,3),1);
    end
    Jaco_AC=sparse(Jaco_AC);
    res_AC=sparse(res_AC);
    Cmatrix_AC=spdiags(Cvector_AC,0,size(solution_vector_AC,1),size(solution_vector_AC,1));
    Jaco_scaled_AC=Jaco_AC*Cmatrix_AC;
    Rvector_AC=1./sum(abs(Jaco_scaled_AC),2);
    Rmatrix_AC=spdiags(Rvector_AC,0,size(solution_vector_AC,1),size(solution_vector_AC,1));
    Jaco_scaled_AC=Rmatrix_AC*Jaco_scaled_AC;
    res_scaled_AC=Rmatrix_AC*res_AC;
    update_scaled_AC=Jaco_scaled_AC\(res_scaled_AC);
    update_AC(:,:)=Cmatrix_AC*update_scaled_AC;

    solution_vector_AC(:,:)=update_AC(:,:);
    solution_vector_phi_AC=solution_vector_AC(:,1:3:length(solution_vector_AC));
    solution_vector_elec_AC=solution_vector_AC(:,2:3:length(solution_vector_AC));
    solution_vector_hole_AC=solution_vector_AC(:,3:3:length(solution_vector_AC));

    FILENAME = sprintf('./data/homework21_AC_freq_%.2EHz.mat' , freq);
    save(FILENAME)
end