clear; close all; clc;

freq=1.00E+09;
w=2*pi*freq;
FILENAME = sprintf('./data/homework21_AC_freq_%.2EHz.mat' , freq);
load(FILENAME)

%% Current %%
excitation_node=522; perturbed="elec";

if perturbed=="potential"
    solution_vector_AC_for_current=solution_vector_phi_AC(:,excitation_node);
elseif perturbed=="elec"
    solution_vector_AC_for_current=solution_vector_elec_AC(:,excitation_node);
elseif perturbed=="hole"
    solution_vector_AC_for_current=solution_vector_hole_AC(:,excitation_node);
end

% Cathode current
J_AC_cathode(1,freq_for)=0; J_n_AC_cathode(1,freq_for)=0; J_p_AC_cathode(1,freq_for)=0; J_d_AC_cathode(1,freq_for)=0;
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

        %%% small signal %%%
        % potential
        x1_delta=solution_vector_AC_for_current(index(1,1),freq_for);
        x2_delta=solution_vector_AC_for_current(index(4,1),freq_for);
        x3_delta=solution_vector_AC_for_current(index(7,1),freq_for);
        x12_delta=x1_delta-x2_delta; x21_delta=-x12_delta;
        x23_delta=x2_delta-x3_delta; x32_delta=-x23_delta;
        x31_delta=x3_delta-x1_delta; x13_delta=-x31_delta;
        x12_delta_avr=(x1_delta+x2_delta)/2; x23_delta_avr=(x2_delta+x3_delta)/2; x31_delta_avr=(x3_delta+x1_delta)/2;
        x21_delta_avr=x12_delta_avr; x32_delta_avr=x23_delta_avr; x13_delta_avr=x31_delta_avr;

        % electron
        n1_delta=solution_vector_AC_for_current(index(1,1)+1,freq_for);
        n2_delta=solution_vector_AC_for_current(index(4,1)+1,freq_for);
        n3_delta=solution_vector_AC_for_current(index(7,1)+1,freq_for);
        n12_delta=n1_delta-n2_delta; n21_delta=-n12_delta;
        n23_delta=n2_delta-n3_delta; n32_delta=-n23_delta;
        n31_delta=n3_delta-n1_delta; n13_delta=-n31_delta;
        n12_delta_avr=(n1_delta+n2_delta)/2; n23_delta_avr=(n2_delta+n3_delta)/2; n31_delta_avr=(n3_delta+n1_delta)/2;
        n21_delta_avr=n12_delta_avr; n32_delta_avr=n23_delta_avr; n13_delta_avr=n31_delta_avr;

        % hole
        p1_delta=solution_vector_AC_for_current(index(1,1)+2,freq_for);
        p2_delta=solution_vector_AC_for_current(index(4,1)+2,freq_for);
        p3_delta=solution_vector_AC_for_current(index(7,1)+2,freq_for);
        p12_delta=p1_delta-p2_delta; p21_delta=-p12_delta;
        p23_delta= p2_delta-p3_delta; p32_delta=-p23_delta;
        p31_delta= p3_delta-p1_delta; p13_delta=-p31_delta;
        p12_delta_avr=(p1_delta+p2_delta)/2; p23_delta_avr=(p2_delta+p3_delta)/2; p31_delta_avr=(p3_delta+p1_delta)/2;
        p21_delta_avr=p12_delta_avr; p32_delta_avr=p23_delta_avr; p13_delta_avr=p31_delta_avr;


        index_element=find(Table_element_region(:,1)==1 & Table_element_region(:,2)==row(k,1));
        eps_si=11.7; J_n_AC_tmp=0; J_p_AC_tmp=0; J_d_AC_tmp=0;
        if col(k,1)==1
            J_n_AC_tmp=J_n_AC_tmp+coeff_dJn*((edge(ii,1)/L(ii,1))*((n21_dc_avr*x21_delta)+(n21_delta_avr*x21_dc)-V_T*n21_delta));
            J_n_AC_tmp=J_n_AC_tmp+coeff_dJn*((edge(ii,3)/L(ii,3))*((n31_dc_avr*x31_delta)+(n31_delta_avr*x31_dc)-V_T*n31_delta));
            J_p_AC_tmp=J_p_AC_tmp+coeff_dJp*((edge(ii,1)/L(ii,1))*((p21_dc_avr*x21_delta)+(p21_delta_avr*x21_dc)+V_T*p21_delta));
            J_p_AC_tmp=J_p_AC_tmp+coeff_dJp*((edge(ii,3)/L(ii,3))*((p31_dc_avr*x31_delta)+(p31_delta_avr*x31_dc)+V_T*p31_delta));
            J_d_AC_tmp=J_d_AC_tmp-eps_si*eps0*(x21_delta)*1i*w/L(index_element,1)*edge(index_element,1);
            J_d_AC_tmp=J_d_AC_tmp-eps_si*eps0*(x31_delta)*1i*w/L(index_element,3)*edge(index_element,3);
        elseif col(k,1)==2
            J_n_AC_tmp=J_n_AC_tmp+coeff_dJn*((edge(index_element,2)/L(index_element,2))*((n32_dc_avr*x32_delta)+(n32_delta_avr*x32_dc)-V_T*n32_delta));
            J_n_AC_tmp=J_n_AC_tmp+coeff_dJn*((edge(index_element,1)/L(index_element,1))*((n12_dc_avr*x12_delta)+(n12_delta_avr*x12_dc)-V_T*n12_delta));
            J_p_AC_tmp=J_p_AC_tmp+coeff_dJp*((edge(index_element,2)/L(index_element,2))*((p32_dc_avr*x32_delta)+(p32_delta_avr*x32_dc)+V_T*p32_delta));
            J_p_AC_tmp=J_p_AC_tmp+coeff_dJp*((edge(index_element,1)/L(index_element,1))*((p12_dc_avr*x12_delta)+(p12_delta_avr*x12_dc)+V_T*p12_delta));
            J_d_AC_tmp=J_d_AC_tmp-eps_si*eps0*(x32_delta)*1i*w/L(index_element,2)*edge(index_element,2);
            J_d_AC_tmp=J_d_AC_tmp-eps_si*eps0*(x12_delta)*1i*w/L(index_element,1)*edge(index_element,1);
        elseif col(k,1)==3
            J_n_AC_tmp=J_n_AC_tmp+coeff_dJn*((edge(index_element,3)/L(index_element,3))*((n13_dc_avr*x13_delta)+(n13_delta_avr*x13_dc)-V_T*n13_delta));
            J_n_AC_tmp=J_n_AC_tmp+coeff_dJn*((edge(index_element,2)/L(index_element,2))*((n23_dc_avr*x23_delta)+(n23_delta_avr*x23_dc)-V_T*n23_delta));
            J_p_AC_tmp=J_p_AC_tmp+coeff_dJp*((edge(index_element,3)/L(index_element,3))*((p13_dc_avr*x13_delta)+(p13_delta_avr*x13_dc)+V_T*p13_delta));
            J_p_AC_tmp=J_p_AC_tmp+coeff_dJp*((edge(index_element,2)/L(index_element,2))*((p23_dc_avr*x23_delta)+(p23_delta_avr*x23_dc)+V_T*p23_delta));
            J_d_AC_tmp=J_d_AC_tmp-eps_si*eps0*(x13_delta)*1i*w/L(index_element,3)*edge(index_element,3);
            J_d_AC_tmp=J_d_AC_tmp-eps_si*eps0*(x23_delta)*1i*w/L(index_element,2)*edge(index_element,2);
        end
        J_n_AC_cathode(1,freq_for)=J_n_AC_cathode(1,freq_for)+J_n_AC_tmp;
        J_p_AC_cathode(1,freq_for)=J_p_AC_cathode(1,freq_for)+J_p_AC_tmp;
        J_d_AC_cathode(1,freq_for)=J_d_AC_cathode(1,freq_for)+J_d_AC_tmp;
    end
end
J_AC_cathode(1,freq_for) = (J_n_AC_cathode(1,freq_for)+J_p_AC_cathode(1,freq_for)+J_d_AC_cathode(1,freq_for));

J_AC_anode(1,freq_for)=0; J_n_AC_anode(1,freq_for)=0; J_p_AC_anode(1,freq_for)=0; J_d_AC_anode(1,freq_for)=0;
for ii=1:size(anode_BC,1)
    vertex_current=anode_BC(ii,1);
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

        %%% small signal %%%
        % potential
        x1_delta=solution_vector_AC_for_current(index(1,1),freq_for);
        x2_delta=solution_vector_AC_for_current(index(4,1),freq_for);
        x3_delta=solution_vector_AC_for_current(index(7,1),freq_for);
        x12_delta=x1_delta-x2_delta; x21_delta=-x12_delta;
        x23_delta=x2_delta-x3_delta; x32_delta=-x23_delta;
        x31_delta=x3_delta-x1_delta; x13_delta=-x31_delta;
        x12_delta_avr=(x1_delta+x2_delta)/2; x23_delta_avr=(x2_delta+x3_delta)/2; x31_delta_avr=(x3_delta+x1_delta)/2;
        x21_delta_avr=x12_delta_avr; x32_delta_avr=x23_delta_avr; x13_delta_avr=x31_delta_avr;

        % electron
        n1_delta=solution_vector_AC_for_current(index(1,1)+1,freq_for);
        n2_delta=solution_vector_AC_for_current(index(4,1)+1,freq_for);
        n3_delta=solution_vector_AC_for_current(index(7,1)+1,freq_for);
        n12_delta=n1_delta-n2_delta; n21_delta=-n12_delta;
        n23_delta=n2_delta-n3_delta; n32_delta=-n23_delta;
        n31_delta=n3_delta-n1_delta; n13_delta=-n31_delta;
        n12_delta_avr=(n1_delta+n2_delta)/2; n23_delta_avr=(n2_delta+n3_delta)/2; n31_delta_avr=(n3_delta+n1_delta)/2;
        n21_delta_avr=n12_delta_avr; n32_delta_avr=n23_delta_avr; n13_delta_avr=n31_delta_avr;

        % hole
        p1_delta=solution_vector_AC_for_current(index(1,1)+2,freq_for);
        p2_delta=solution_vector_AC_for_current(index(4,1)+2,freq_for);
        p3_delta=solution_vector_AC_for_current(index(7,1)+2,freq_for);
        p12_delta=p1_delta-p2_delta; p21_delta=-p12_delta;
        p23_delta= p2_delta-p3_delta; p32_delta=-p23_delta;
        p31_delta= p3_delta-p1_delta; p13_delta=-p31_delta;
        p12_delta_avr=(p1_delta+p2_delta)/2; p23_delta_avr=(p2_delta+p3_delta)/2; p31_delta_avr=(p3_delta+p1_delta)/2;
        p21_delta_avr=p12_delta_avr; p32_delta_avr=p23_delta_avr; p13_delta_avr=p31_delta_avr;


        index_element=find(Table_element_region(:,1)==1 & Table_element_region(:,2)==row(k,1));
        eps_si=11.7; J_n_AC_tmp=0; J_p_AC_tmp=0; J_d_AC_tmp=0;
        if col(k,1)==1
            J_n_AC_tmp=J_n_AC_tmp+coeff_dJn*((edge(ii,1)/L(ii,1))*((n21_dc_avr*x21_delta)+(n21_delta_avr*x21_dc)-V_T*n21_delta));
            J_n_AC_tmp=J_n_AC_tmp+coeff_dJn*((edge(ii,3)/L(ii,3))*((n31_dc_avr*x31_delta)+(n31_delta_avr*x31_dc)-V_T*n31_delta));
            J_p_AC_tmp=J_p_AC_tmp+coeff_dJp*((edge(ii,1)/L(ii,1))*((p21_dc_avr*x21_delta)+(p21_delta_avr*x21_dc)+V_T*p21_delta));
            J_p_AC_tmp=J_p_AC_tmp+coeff_dJp*((edge(ii,3)/L(ii,3))*((p31_dc_avr*x31_delta)+(p31_delta_avr*x31_dc)+V_T*p31_delta));
            J_d_AC_tmp=J_d_AC_tmp-eps_si*eps0*(x21_delta)*1i*w/L(index_element,1)*edge(index_element,1);
            J_d_AC_tmp=J_d_AC_tmp-eps_si*eps0*(x31_delta)*1i*w/L(index_element,3)*edge(index_element,3);
        elseif col(k,1)==2
            J_n_AC_tmp=J_n_AC_tmp+coeff_dJn*((edge(index_element,2)/L(index_element,2))*((n32_dc_avr*x32_delta)+(n32_delta_avr*x32_dc)-V_T*n32_delta));
            J_n_AC_tmp=J_n_AC_tmp+coeff_dJn*((edge(index_element,1)/L(index_element,1))*((n12_dc_avr*x12_delta)+(n12_delta_avr*x12_dc)-V_T*n12_delta));
            J_p_AC_tmp=J_p_AC_tmp+coeff_dJp*((edge(index_element,2)/L(index_element,2))*((p32_dc_avr*x32_delta)+(p32_delta_avr*x32_dc)+V_T*p32_delta));
            J_p_AC_tmp=J_p_AC_tmp+coeff_dJp*((edge(index_element,1)/L(index_element,1))*((p12_dc_avr*x12_delta)+(p12_delta_avr*x12_dc)+V_T*p12_delta));
            J_d_AC_tmp=J_d_AC_tmp-eps_si*eps0*(x32_delta)*1i*w/L(index_element,2)*edge(index_element,2);
            J_d_AC_tmp=J_d_AC_tmp-eps_si*eps0*(x12_delta)*1i*w/L(index_element,1)*edge(index_element,1);
        elseif col(k,1)==3
            J_n_AC_tmp=J_n_AC_tmp+coeff_dJn*((edge(index_element,3)/L(index_element,3))*((n13_dc_avr*x13_delta)+(n13_delta_avr*x13_dc)-V_T*n13_delta));
            J_n_AC_tmp=J_n_AC_tmp+coeff_dJn*((edge(index_element,2)/L(index_element,2))*((n23_dc_avr*x23_delta)+(n23_delta_avr*x23_dc)-V_T*n23_delta));
            J_p_AC_tmp=J_p_AC_tmp+coeff_dJp*((edge(index_element,3)/L(index_element,3))*((p13_dc_avr*x13_delta)+(p13_delta_avr*x13_dc)+V_T*p13_delta));
            J_p_AC_tmp=J_p_AC_tmp+coeff_dJp*((edge(index_element,2)/L(index_element,2))*((p23_dc_avr*x23_delta)+(p23_delta_avr*x23_dc)+V_T*p23_delta));
            J_d_AC_tmp=J_d_AC_tmp-eps_si*eps0*(x13_delta)*1i*w/L(index_element,3)*edge(index_element,3);
            J_d_AC_tmp=J_d_AC_tmp-eps_si*eps0*(x23_delta)*1i*w/L(index_element,2)*edge(index_element,2);
        end
        J_n_AC_anode(1,freq_for)=J_n_AC_anode(1,freq_for)+J_n_AC_tmp;
        J_p_AC_anode(1,freq_for)=J_p_AC_anode(1,freq_for)+J_p_AC_tmp;
        J_d_AC_anode(1,freq_for)=J_d_AC_anode(1,freq_for)+J_d_AC_tmp;
    end
end
J_AC_anode(1,freq_for) = (J_n_AC_anode(1,freq_for)+J_p_AC_anode(1,freq_for)+J_d_AC_anode(1,freq_for));
sum=(J_d_AC_cathode+J_d_AC_anode);

%% 검증 (delta_n+delta+p)*q*Control_volume=delta_q인것을 확인
delta_n=0; delta_p=0;
for ii=1:length(Element)
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

    for j=1:3
        delta_n=delta_n+Control_Volume(j,1)*solution_vector_AC_for_current(index(3*j-1,1),freq_for);
        delta_p=delta_p+Control_Volume(j,1)*solution_vector_AC_for_current(index(3*j,1),freq_for);
    end
end
prove=(delta_n+delta_p)*1i*w*q;