clear; clc;
load("Ramping_information.mat")

coeff_dJn=-q*m_n;
coeff_dJp=-q*m_p;
i=sqrt(-1);

dphi=zeros(VR2,VR2);
delec=nint*exp(dphi);
dhole=nint*exp(dphi);


%%%%% Green's function %%%%%

freq = 1e+12;
Amp = 1e-3; 

Jn=zeros(1,1);
Jp=zeros(1,1);
Jd=zeros(1,1);

Jaco_ac = sparse(T_row,T_row);
res_ac = eye(T_row,T_row);

for ii=1:T_row

    if Table(ii,1) == 1     % Region1
        K = find(E_R1 == Table(ii,2));
        K_row = size(K,1);

        for rr = 1:K_row

            if Table(ii,3) == Potential

                if K(rr,1) <= R_row1

                    jj = 1;
                    V3 = find(Table(1:VR1,2)==E_R1(K(rr,1),jj+2)) ;
                    V2 = find(Table(1:VR1,2)==E_R1(K(rr,1),jj+1));
                    V1 = find(Table(1:VR1,2)==E_R1(K(rr,1),jj));

                    Jaco_ac(ii, V3) = Jaco_ac(ii, V3)+ eox*(Area(K(rr,1),jj+2)/L(K(rr,1),jj+2));
                    Jaco_ac(ii, V2) = Jaco_ac(ii, V2)+eox*(Area(K(rr,1),jj)/L(K(rr,1),jj));
                    Jaco_ac(ii, V1) = Jaco_ac(ii, V1)-eox*(Area(K(rr,1),jj+2)/L(K(rr,1),jj+2))-eox*(Area(K(rr,1),jj)/L(K(rr,1),jj));

                elseif K(rr,1) > R_row1 && K(rr,1) <= 2*R_row1
                    jj = 2;
                    K(rr,1) = K(rr,1)-R_row1;

                    V3 = find(Table(1:VR1,2)==E_R1(K(rr,1),jj+1)) ;
                    V2 = find(Table(1:VR1,2)==E_R1(K(rr,1),jj));
                    V1 = find(Table(1:VR1,2)==E_R1(K(rr,1),jj-1));

                    Jaco_ac(ii, V3) = Jaco_ac(ii, V3)+eox*(Area(K(rr,1),jj)/L(K(rr,1),jj));
                    Jaco_ac(ii, V2) =  Jaco_ac(ii, V2)-eox*(Area(K(rr,1),jj)/L(K(rr,1),jj))-eox*(Area(K(rr,1),jj-1)/L(K(rr,1),jj-1));
                    Jaco_ac(ii, V1) = Jaco_ac(ii, V1)+eox*(Area(K(rr,1),jj-1)/L(K(rr,1),jj-1));

                elseif K(rr,1) > 2*R_row1 && K(rr,1) <= 3*R_row1
                    jj = 3;
                    K(rr,1) = K(rr,1)-2*R_row1;

                    V3 = find(Table(1:VR1,2)==E_R1(K(rr,1),jj-2)) ;
                    V2 = find(Table(1:VR1,2)==E_R1(K(rr,1),jj-1));
                    V1 = find(Table(1:VR1,2)==E_R1(K(rr,1),jj));

                    Jaco_ac(ii, V3) = Jaco_ac(ii, V3)+eox*(Area(K(rr,1),jj)/L(K(rr,1),jj));
                    Jaco_ac(ii, V2) = Jaco_ac(ii, V2)+eox*(Area(K(rr,1),jj-1)/L(K(rr,1),jj-1));
                    Jaco_ac(ii, V1) = Jaco_ac(ii, V1)-eox*(Area(K(rr,1),jj)/L(K(rr,1),jj))-eox*(Area(K(rr,1),jj-1)/L(K(rr,1),jj-1));

                end
            end
        end

    elseif Table(ii,1) == 2     %region2

        K = find(E_R2 == Table(ii,2));
        K_row = size(K,1);

        if Table(ii,3) == Potential

            for rr = 1:K_row

                if K(rr,1) <= R_row2
                    jj = 1;

                    V3 = find(Table(1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj+2)) ;
                    V2 = find(Table(1+1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj+1));
                    V1 = find(Table(1+1+1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj));

                    Jaco_ac(ii, 3*(V3)+(VR1-2)) = Jaco_ac(ii, 3*(V3)+(VR1-2))+ esi*(Area(K(rr,1)+R_row1,jj+2)/L(K(rr,1)+R_row1,jj+2));
                    Jaco_ac(ii, 3*(V2)+(VR1-2)) = Jaco_ac(ii, 3*(V2)+(VR1-2))+ esi*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj));
                    Jaco_ac(ii, 3*(V1)+(VR1-2)) = Jaco_ac(ii, 3*(V1)+(VR1-2))- esi*(Area(K(rr,1)+R_row1,jj+2)/L(K(rr,1)+R_row1,jj+2))-esi*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj));

                    Jaco_ac(ii, 3*(V1)+(VR1-1)) = Jaco_ac(ii, 3*(V1)+(VR1-1)) - coeff*(Area(K(rr,1)+R_row1,jj)*L(K(rr,1)+R_row1,jj)+Area(K(rr,1)+R_row1,jj+2)*L(K(rr,1)+R_row1,jj+2));
                    Jaco_ac(ii, 3*(V1)+(VR1)) = Jaco_ac(ii, 3*(V1)+(VR1)) + coeff*(Area(K(rr,1)+R_row1,jj)*L(K(rr,1)+R_row1,jj)+Area(K(rr,1)+R_row1,jj+2)*L(K(rr,1)+R_row1,jj+2));

                elseif K(rr,1) > R_row2 && K(rr,1) <= 2*R_row2
                    jj = 2;
                    K(rr,1) = K(rr,1)-R_row2;

                    V3 = find(Table(1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj+1)) ;
                    V2 = find(Table(1+1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj));
                    V1 = find(Table(1+1+1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj-1));

                    Jaco_ac(ii, 3*(V3)+(VR1-2)) = Jaco_ac(ii, 3*(V3)+(VR1-2)) + esi*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj));
                    Jaco_ac(ii, 3*(V2)+(VR1-2)) =  Jaco_ac(ii, 3*(V2)+(VR1-2)) - esi*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj))-esi*(Area(K(rr,1)+R_row1,jj-1)/L(K(rr,1)+R_row1,jj-1));
                    Jaco_ac(ii, 3*(V1)+(VR1-2)) = Jaco_ac(ii, 3*(V1)+(VR1-2)) + esi*(Area(K(rr,1)+R_row1,jj-1)/L(K(rr,1)+R_row1,jj-1));

                    Jaco_ac(ii, 3*(V2)+(VR1-1)) = Jaco_ac(ii, 3*(V2)+(VR1-1)) - coeff*(Area(K(rr,1)+R_row1,jj)*L(K(rr,1)+R_row1,jj)+Area(K(rr,1)+R_row1,jj-1)*L(K(rr,1)+R_row1,jj-1));
                    Jaco_ac(ii, 3*(V2)+(VR1)) = Jaco_ac(ii, 3*(V2)+(VR1)) + coeff*(Area(K(rr,1)+R_row1,jj)*L(K(rr,1)+R_row1,jj)+Area(K(rr,1)+R_row1,jj-1)*L(K(rr,1)+R_row1,jj-1));

                elseif K(rr,1) > 2*R_row2 && K(rr,1) <= 3*R_row2
                    jj = 3;
                    K(rr,1) = K(rr,1)-2*R_row2;

                    V3 = find(Table(1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj-2)) ;
                    V2 = find(Table(1+1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj-1));
                    V1 = find(Table(1+1+1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj));

                    Jaco_ac(ii, 3*(V3)+(VR1-2)) = Jaco_ac(ii, 3*(V3)+(VR1-2)) + esi*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj));
                    Jaco_ac(ii, 3*(V2)+(VR1-2)) = Jaco_ac(ii, 3*(V2)+(VR1-2)) + esi*(Area(K(rr,1)+R_row1,jj-1)/L(K(rr,1)+R_row1,jj-1));
                    Jaco_ac(ii, 3*(V1)+(VR1-2)) = Jaco_ac(ii, 3*(V1)+(VR1-2)) - esi*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj))-esi*(Area(K(rr,1)+R_row1,jj-1)/L(K(rr,1)+R_row1,jj-1));

                    Jaco_ac(ii, 3*(V1)+(VR1-1)) = Jaco_ac(ii, 3*(V1)+(VR1-1)) - coeff*(Area(K(rr,1)+R_row1,jj)*L(K(rr,1)+R_row1,jj)+Area(K(rr,1)+R_row1,jj-1)*L(K(rr,1)+R_row1,jj-1));
                    Jaco_ac(ii, 3*(V1)+(VR1)) = Jaco_ac(ii, 3*(V1)+(VR1)) + coeff*(Area(K(rr,1)+R_row1,jj)*L(K(rr,1)+R_row1,jj)+Area(K(rr,1)+R_row1,jj-1)*L(K(rr,1)+R_row1,jj-1));
                end
            end

        elseif Table(ii,3) == eDensity

            for rr = 1:K_row
                if K(rr,1) <= R_row2
                    jj = 1;

                    V3 = find(Table(1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj+2)) ;  % node 3
                    V2 = find(Table(1+1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj+1)); % node 2
                    V1 = find(Table(2+1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj));   % node 1

                    elecDC_av1 = 0.5*(elec_DC(VR1+V2,1)+elec_DC(VR1+V1,1));
                    delec_av1 = 0.5*(delec(VR1+V2,1)+delec(VR1+V1,1));
                    ddphi1 = dphi(VR1+V2,1)-dphi(VR1+V1,1);
                    ddelec1 = delec(VR1+V2,1)-delec(VR1+V1,1);
                    dphi_DC1 = phi_DC(VR1+V2,1)-phi_DC(VR1+V1,1);

                    elecDC_av2 = 0.5*(elec_DC(VR1+V3,1)+elec_DC(VR1+V1,1));
                    delec_av2 = 0.5*(delec(VR1+V3,1)+delec(VR1+V1,1));
                    ddphi2 = dphi(VR1+V3,1)-dphi(VR1+V1,1);
                    ddelec2 = delec(VR1+V3,1)-delec(VR1+V1,1);
                    dphi_DC2 = phi_DC(VR1+V3,1)-phi_DC(VR1+V1,1);

                    % eDensity
                    Jaco_ac(ii, 3*(V1)+(VR1-1)) = Jaco_ac(ii, 3*(V1)+(VR1-1)) + coeff_dJn*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj))*(0.5*dphi_DC1+Thermal_V);
                    Jaco_ac(ii, 3*(V1)+(VR1-1)) = Jaco_ac(ii, 3*(V1)+(VR1-1)) + coeff_dJn*(Area(K(rr,1)+R_row1,jj+2)/L(K(rr,1)+R_row1,jj+2))*(0.5*dphi_DC2+Thermal_V);

                    Jaco_ac(ii, 3*(V2)+(VR1-1)) = Jaco_ac(ii, 3*(V2)+(VR1-1)) + coeff_dJn*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj))*(0.5*dphi_DC1-Thermal_V);
                    Jaco_ac(ii, 3*(V3)+(VR1-1)) = Jaco_ac(ii, 3*(V3)+(VR1-1)) + coeff_dJn*(Area(K(rr,1)+R_row1,jj+2)/L(K(rr,1)+R_row1,jj+2))*(0.5*dphi_DC2-Thermal_V);

                    % AC
                    Jaco_ac(ii,3*(V1)+(VR1-1)) = Jaco_ac(ii,3*(V1)+(VR1-1))  - (q/4)*(2*pi*freq*i)*(Area(K(rr,1)+R_row1,jj)*L(K(rr,1)+R_row1,jj)+Area(K(rr,1)+R_row1,jj+2)*L(K(rr,1)+R_row1,jj+2));

                    % Potential
                    Jaco_ac(ii, 3*(V1)+(VR1-2)) = Jaco_ac(ii, 3*(V1)+(VR1-2)) + coeff_dJn*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj))*(-elecDC_av1);
                    Jaco_ac(ii, 3*(V1)+(VR1-2)) = Jaco_ac(ii, 3*(V1)+(VR1-2)) + coeff_dJn*(Area(K(rr,1)+R_row1,jj+2)/L(K(rr,1)+R_row1,jj+2))*(-elecDC_av2);

                    Jaco_ac(ii, 3*(V2)+(VR1-2)) = Jaco_ac(ii, 3*(V2)+(VR1-2)) + coeff_dJn*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj))*elecDC_av1;
                    Jaco_ac(ii, 3*(V3)+(VR1-2)) = Jaco_ac(ii, 3*(V3)+(VR1-2)) + coeff_dJn*(Area(K(rr,1)+R_row1,jj+2)/L(K(rr,1)+R_row1,jj+2))*elecDC_av2;


                elseif K(rr,1) > R_row2 && K(rr,1) <= 2*R_row2
                    jj = 2;
                    K(rr,1) = K(rr,1)-R_row2;

                    V3 = find(Table(1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj+1)) ; % node 3
                    V2 = find(Table(1+1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj));  % node 2
                    V1 = find(Table(2+1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj-1)); % node 1

                    elecDC_av1 = 0.5*(elec_DC(VR1+V1,1)+elec_DC(VR1+V2,1));
                    delec_av1 = 0.5*(delec(VR1+V1,1)+delec(VR1+V2,1));
                    ddphi1 = dphi(VR1+V1,1)-dphi(VR1+V2,1);
                    ddelec1 = delec(VR1+V1,1)-delec(VR1+V2,1);
                    dphi_DC1 = phi_DC(VR1+V1,1)-phi_DC(VR1+V2,1);

                    elecDC_av2 = 0.5*(elec_DC(VR1+V3,1)+elec_DC(VR1+V2,1));
                    delec_av2 = 0.5*(delec(VR1+V3,1)+delec(VR1+V2,1));
                    ddphi2 = dphi(VR1+V3,1)-dphi(VR1+V2,1);
                    ddelec2 = delec(VR1+V3,1)-delec(VR1+V2,1);
                    dphi_DC2 = phi_DC(VR1+V3,1)-phi_DC(VR1+V2,1);

                    % eDensity
                    Jaco_ac(ii, 3*(V2)+(VR1-1)) = Jaco_ac(ii, 3*(V2)+(VR1-1)) + coeff_dJn*(Area(K(rr,1)+R_row1,jj-1)/L(K(rr,1)+R_row1,jj-1))*(0.5*dphi_DC1+Thermal_V);
                    Jaco_ac(ii, 3*(V2)+(VR1-1)) = Jaco_ac(ii, 3*(V2)+(VR1-1)) + coeff_dJn*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj))*(0.5*dphi_DC2+Thermal_V);

                    Jaco_ac(ii, 3*(V1)+(VR1-1)) = Jaco_ac(ii, 3*(V1)+(VR1-1)) + coeff_dJn*(Area(K(rr,1)+R_row1,jj-1)/L(K(rr,1)+R_row1,jj-1))*(0.5*dphi_DC1-Thermal_V);
                    Jaco_ac(ii, 3*(V3)+(VR1-1)) = Jaco_ac(ii, 3*(V3)+(VR1-1)) + coeff_dJn*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj))*(0.5*dphi_DC2-Thermal_V);

                    % AC
                    Jaco_ac(ii,3*(V2)+(VR1-1)) = Jaco_ac(ii,3*(V2)+(VR1-1)) - (q/4)*(2*pi*freq*i)*(Area(K(rr,1)+R_row1,jj-1)*L(K(rr,1)+R_row1,jj-1)+Area(K(rr,1)+R_row1,jj)*L(K(rr,1)+R_row1,jj));

                    % Potential
                    Jaco_ac(ii, 3*(V2)+(VR1-2)) = Jaco_ac(ii, 3*(V2)+(VR1-2)) + coeff_dJn*(Area(K(rr,1)+R_row1,jj-1)/L(K(rr,1)+R_row1,jj-1))*(-elecDC_av1);
                    Jaco_ac(ii, 3*(V2)+(VR1-2)) = Jaco_ac(ii, 3*(V2)+(VR1-2)) + coeff_dJn*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj))*(-elecDC_av2);

                    Jaco_ac(ii, 3*(V1)+(VR1-2)) = Jaco_ac(ii, 3*(V1)+(VR1-2)) + coeff_dJn*(Area(K(rr,1)+R_row1,jj-1)/L(K(rr,1)+R_row1,jj-1))*elecDC_av1;
                    Jaco_ac(ii, 3*(V3)+(VR1-2)) = Jaco_ac(ii, 3*(V3)+(VR1-2)) + coeff_dJn*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj))*elecDC_av2;

                elseif K(rr,1) > 2*R_row2 && K(rr,1) <= 3*R_row2

                    jj = 3;
                    K(rr,1) = K(rr,1)-2*R_row2;

                    V3 = find(Table(2+1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj)); % node 3
                    V2 = find(Table(1+1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj-1)); % node 2
                    V1 = find(Table(1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj-2)) ;  % node 1

                    elecDC_av1 = 0.5*(elec_DC(VR1+V1,1)+elec_DC(VR1+V3,1));
                    delec_av1 = 0.5*(delec(VR1+V1,1)+delec(VR1+V3,1));
                    ddphi1 = dphi(VR1+V1,1)-dphi(VR1+V3,1);
                    ddelec1 = delec(VR1+V1,1)-delec(VR1+V3,1);
                    dphi_DC1 = phi_DC(VR1+V1,1)-phi_DC(VR1+V3,1);

                    elecDC_av2 = 0.5*(elec_DC(VR1+V2,1)+elec_DC(VR1+V3,1));
                    delec_av2 = 0.5*(delec(VR1+V2,1)+delec(VR1+V3,1));
                    ddphi2 = dphi(VR1+V2,1)-dphi(VR1+V3,1);
                    ddelec2 = delec(VR1+V2,1)-delec(VR1+V3,1);
                    dphi_DC2 = phi_DC(VR1+V2,1)-phi_DC(VR1+V3,1);


                    % eDensity
                    Jaco_ac(ii, 3*(V1)+(VR1-1)) = Jaco_ac(ii, 3*(V1)+(VR1-1)) + coeff_dJn*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj))*(0.5*dphi_DC1-Thermal_V);
                    Jaco_ac(ii, 3*(V2)+(VR1-1)) = Jaco_ac(ii, 3*(V2)+(VR1-1)) + coeff_dJn*(Area(K(rr,1)+R_row1,jj-1)/L(K(rr,1)+R_row1,jj-1))*(0.5*dphi_DC2-Thermal_V);

                    Jaco_ac(ii, 3*(V3)+(VR1-1)) = Jaco_ac(ii, 3*(V3)+(VR1-1)) + coeff_dJn*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj))*(0.5*dphi_DC1+Thermal_V);
                    Jaco_ac(ii, 3*(V3)+(VR1-1)) = Jaco_ac(ii, 3*(V3)+(VR1-1)) + coeff_dJn*(Area(K(rr,1)+R_row1,jj-1)/L(K(rr,1)+R_row1,jj-1))*(0.5*dphi_DC2+Thermal_V);

                    % AC
                    Jaco_ac(ii,3*(V3)+(VR1-1)) = Jaco_ac(ii,3*(V3)+(VR1-1)) - (q/4)*(2*pi*freq*i)*(Area(K(rr,1)+R_row1,jj)*L(K(rr,1)+R_row1,jj)+Area(K(rr,1)+R_row1,jj-1)*L(K(rr,1)+R_row1,jj-1));

                    % Potential
                    Jaco_ac(ii, 3*(V3)+(VR1-2)) = Jaco_ac(ii, 3*(V3)+(VR1-2)) + coeff_dJn*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj))*(-elecDC_av1);
                    Jaco_ac(ii, 3*(V3)+(VR1-2)) = Jaco_ac(ii, 3*(V3)+(VR1-2)) + coeff_dJn*(Area(K(rr,1)+R_row1,jj-1)/L(K(rr,1)+R_row1,jj-1))*(-elecDC_av2);

                    Jaco_ac(ii, 3*(V2)+(VR1-2)) = Jaco_ac(ii, 3*(V2)+(VR1-2)) + coeff_dJn*(Area(K(rr,1)+R_row1,jj-1)/L(K(rr,1)+R_row1,jj-1))*(elecDC_av2);
                    Jaco_ac(ii, 3*(V1)+(VR1-2)) = Jaco_ac(ii, 3*(V1)+(VR1-2)) + coeff_dJn*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj))*(elecDC_av1);
                end
            end


        elseif Table(ii,3) == hDensity

            for rr = 1:K_row

                if K(rr,1) <= R_row2
                    jj = 1;

                    V3 = find(Table(1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj+2)) ;  % node 3
                    V2 = find(Table(1+1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj+1));  % node 2
                    V1 = find(Table(2+1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj));  % node 1

                    holeDC_av1 = 0.5*(hole_DC(VR1+V2,1)+hole_DC(VR1+V1,1));
                    dhole_av1 = 0.5*(dhole(VR1+V2,1)+dhole(VR1+V1,1));
                    ddphi1 = dphi(VR1+V2,1)-dphi(VR1+V1,1);
                    ddhole1 = dhole(VR1+V2,1)-dhole(VR1+V1,1);
                    dphi_DC1 = phi_DC(VR1+V2,1)-phi_DC(VR1+V1,1);

                    holeDC_av2 = 0.5*(hole_DC(VR1+V3,1)+hole_DC(VR1+V1,1));
                    dhole_av2 = 0.5*(dhole(VR1+V3,1)+dhole(VR1+V1,1));
                    ddphi2 = dphi(VR1+V3,1)-dphi(VR1+V1,1);
                    ddhole2 = dhole(VR1+V3,1)-dhole(VR1+V1,1);
                    dphi_DC2 = phi_DC(VR1+V3,1)-phi_DC(VR1+V1,1);


                    % hDensity
                    Jaco_ac(ii, 3*(V1)+(VR1)) = Jaco_ac(ii, 3*(V1)+(VR1)) + coeff_dJp*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj))*(0.5*dphi_DC1-Thermal_V);
                    Jaco_ac(ii, 3*(V1)+(VR1)) = Jaco_ac(ii, 3*(V1)+(VR1)) + coeff_dJp*(Area(K(rr,1)+R_row1,jj+2)/L(K(rr,1)+R_row1,jj+2))*(0.5*dphi_DC2-Thermal_V);

                    Jaco_ac(ii, 3*(V2)+(VR1)) = Jaco_ac(ii, 3*(V2)+(VR1)) + coeff_dJp*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj))*(0.5*dphi_DC1+Thermal_V);
                    Jaco_ac(ii, 3*(V3)+(VR1)) = Jaco_ac(ii, 3*(V3)+(VR1)) + coeff_dJp*(Area(K(rr,1)+R_row1,jj+2)/L(K(rr,1)+R_row1,jj+2))*(0.5*dphi_DC2+Thermal_V);

                    % AC
                    Jaco_ac(ii,3*(V1)+(VR1)) = Jaco_ac(ii,3*(V1)+(VR1)) + (q/4)*(2*pi*freq*i)*(Area(K(rr,1)+R_row1,jj)*L(K(rr,1)+R_row1,jj)+Area(K(rr,1)+R_row1,jj+2)*L(K(rr,1)+R_row1,jj+2));

                    % Potential
                    Jaco_ac(ii, 3*(V1)+(VR1-2)) = Jaco_ac(ii, 3*(V1)+(VR1-2)) + coeff_dJp*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj))*(-holeDC_av1);
                    Jaco_ac(ii, 3*(V1)+(VR1-2)) = Jaco_ac(ii, 3*(V1)+(VR1-2)) + coeff_dJp*(Area(K(rr,1)+R_row1,jj+2)/L(K(rr,1)+R_row1,jj+2))*(-holeDC_av2);

                    Jaco_ac(ii, 3*(V2)+(VR1-2)) = Jaco_ac(ii, 3*(V2)+(VR1-2)) + coeff_dJp*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj))*holeDC_av1;
                    Jaco_ac(ii, 3*(V3)+(VR1-2)) = Jaco_ac(ii, 3*(V3)+(VR1-2)) + coeff_dJp*(Area(K(rr,1)+R_row1,jj+2)/L(K(rr,1)+R_row1,jj+2))*holeDC_av2;


                elseif K(rr,1) > R_row2 && K(rr,1) <= 2*R_row2
                    jj = 2;
                    K(rr,1) = K(rr,1)-R_row2;

                    V3 = find(Table(1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj+1)) ;  % node 3
                    V2 = find(Table(1+1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj));  % node 2
                    V1 = find(Table(2+1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj-1));  % node 1

                    holeDC_av1 = 0.5*(hole_DC(VR1+V1,1)+hole_DC(VR1+V2,1));
                    dhole_av1 = 0.5*(dhole(VR1+V1,1)+dhole(VR1+V2,1));
                    ddphi1 = dphi(VR1+V1,1)-dphi(VR1+V2,1);
                    ddhole1 = dhole(VR1+V1,1)-dhole(VR1+V2,1);
                    dphi_DC1 = phi_DC(VR1+V1,1)-phi_DC(VR1+V2,1);

                    holeDC_av2 = 0.5*(hole_DC(VR1+V3,1)+hole_DC(VR1+V2,1));
                    dhole_av2 = 0.5*(dhole(VR1+V3,1)+dhole(VR1+V2,1));
                    ddphi2 = dphi(VR1+V3,1)-dphi(VR1+V2,1);
                    ddhole2 = dhole(VR1+V3,1)-dhole(VR1+V2,1);
                    dphi_DC2 = phi_DC(VR1+V3,1)-phi_DC(VR1+V2,1);


                    % hDensity
                    Jaco_ac(ii, 3*(V1)+(VR1)) = Jaco_ac(ii, 3*(V1)+(VR1)) + coeff_dJp*(Area(K(rr,1)+R_row1,jj-1)/L(K(rr,1)+R_row1,jj-1))*(0.5*dphi_DC1+Thermal_V);
                    Jaco_ac(ii, 3*(V3)+(VR1)) = Jaco_ac(ii, 3*(V3)+(VR1)) + coeff_dJp*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj))*(0.5*dphi_DC2+Thermal_V);

                    Jaco_ac(ii, 3*(V2)+(VR1)) = Jaco_ac(ii, 3*(V2)+(VR1)) + coeff_dJp*(Area(K(rr,1)+R_row1,jj-1)/L(K(rr,1)+R_row1,jj-1))*(0.5*dphi_DC1-Thermal_V);
                    Jaco_ac(ii, 3*(V2)+(VR1)) = Jaco_ac(ii, 3*(V2)+(VR1)) + coeff_dJp*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj))*(0.5*dphi_DC2-Thermal_V);

                    % AC
                    Jaco_ac(ii,3*(V2)+(VR1)) = Jaco_ac(ii,3*(V2)+(VR1)) + (q/4)*(2*pi*freq*i)*(Area(K(rr,1)+R_row1,jj-1)*L(K(rr,1)+R_row1,jj-1)+Area(K(rr,1)+R_row1,jj)*L(K(rr,1)+R_row1,jj));

                    % Potential
                    Jaco_ac(ii, 3*(V2)+(VR1-2)) = Jaco_ac(ii, 3*(V2)+(VR1-2)) + coeff_dJp*(Area(K(rr,1)+R_row1,jj-1)/L(K(rr,1)+R_row1,jj-1))*(-holeDC_av1);
                    Jaco_ac(ii, 3*(V2)+(VR1-2)) = Jaco_ac(ii, 3*(V2)+(VR1-2)) + coeff_dJp*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj))*(-holeDC_av2);

                    Jaco_ac(ii, 3*(V1)+(VR1-2)) = Jaco_ac(ii, 3*(V1)+(VR1-2)) + coeff_dJp*(Area(K(rr,1)+R_row1,jj-1)/L(K(rr,1)+R_row1,jj-1))*(holeDC_av1);
                    Jaco_ac(ii, 3*(V3)+(VR1-2)) = Jaco_ac(ii, 3*(V3)+(VR1-2)) + coeff_dJp*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj))*(holeDC_av2);


                elseif K(rr,1) > 2*R_row2 && K(rr,1) <= 3*R_row2
                    jj = 3;
                    K(rr,1) = K(rr,1)-2*R_row2;

                    V3 = find(Table(2+1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj));  % node 3
                    V2 = find(Table(1+1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj-1));  % node 2
                    V1 = find(Table(1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj-2)) ;  % node 1

                    holeDC_av1 = 0.5*(hole_DC(VR1+V1,1)+hole_DC(VR1+V3,1));
                    dhole_av1 = 0.5*(dhole(VR1+V1,1)+dhole(VR1+V3,1));
                    ddphi1 = dphi(VR1+V1,1)-dphi(VR1+V3,1);
                    ddhole1 = dhole(VR1+V1,1)-dhole(VR1+V3,1);
                    dphi_DC1 = phi_DC(VR1+V1,1)-phi_DC(VR1+V3,1);

                    holeDC_av2 = 0.5*(hole_DC(VR1+V2,1)+hole_DC(VR1+V3,1));
                    dhole_av2 = 0.5*(dhole(VR1+V2,1)+dhole(VR1+V3,1));
                    ddphi2 = dphi(VR1+V2,1)-dphi(VR1+V3,1);
                    ddhole2 = dhole(VR1+V2,1)-dhole(VR1+V3,1);
                    dphi_DC2 = phi_DC(VR1+V2,1)-phi_DC(VR1+V3,1);

                    %hDensity
                    Jaco_ac(ii, 3*(V1)+(VR1)) = Jaco_ac(ii, 3*(V1)+(VR1)) + coeff_dJp*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj))*(0.5*dphi_DC1+Thermal_V);
                    Jaco_ac(ii, 3*(V2)+(VR1)) = Jaco_ac(ii, 3*(V2)+(VR1)) + coeff_dJp*(Area(K(rr,1)+R_row1,jj-1)/L(K(rr,1)+R_row1,jj-1))*(0.5*dphi_DC2+Thermal_V);

                    Jaco_ac(ii, 3*(V3)+(VR1)) = Jaco_ac(ii, 3*(V3)+(VR1)) + coeff_dJp*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj))*(0.5*dphi_DC1-Thermal_V);
                    Jaco_ac(ii, 3*(V3)+(VR1)) = Jaco_ac(ii, 3*(V3)+(VR1)) + coeff_dJp*(Area(K(rr,1)+R_row1,jj-1)/L(K(rr,1)+R_row1,jj-1))*(0.5*dphi_DC2-Thermal_V);

                    % AC
                    Jaco_ac(ii,3*(V3)+(VR1)) = Jaco_ac(ii,3*(V3)+(VR1))  + (q/4)*(2*pi*freq*i)*(Area(K(rr,1)+R_row1,jj)*L(K(rr,1)+R_row1,jj)+Area(K(rr,1)+R_row1,jj-1)*L(K(rr,1)+R_row1,jj-1));

                    %Potential
                    Jaco_ac(ii, 3*(V3)+(VR1-2)) = Jaco_ac(ii, 3*(V3)+(VR1-2)) + coeff_dJp*(Area(K(rr,1)+R_row1,jj-1)/L(K(rr,1)+R_row1,jj-1))*(-holeDC_av1);
                    Jaco_ac(ii, 3*(V3)+(VR1-2)) = Jaco_ac(ii, 3*(V3)+(VR1-2)) + coeff_dJp*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj))*(-holeDC_av2);

                    Jaco_ac(ii, 3*(V2)+(VR1-2)) = Jaco_ac(ii, 3*(V2)+(VR1-2)) + coeff_dJp*(Area(K(rr,1)+R_row1,jj-1)/L(K(rr,1)+R_row1,jj-1))*holeDC_av1;
                    Jaco_ac(ii, 3*(V1)+(VR1-2)) = Jaco_ac(ii, 3*(V1)+(VR1-2)) + coeff_dJp*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj))*holeDC_av2;
                end
            end
        end

    elseif Table(ii,1) == 3  % Region3

        K = find(E_R3 == Table(ii,2));
        K_row = size(K,1);

        for rr = 1:K_row
            if Table(ii,3) == Potential

                if K(rr,1) <= R_row3
                    jj = 1;

                    V3 = find(Table(1+VR1+3*VR2:T1,2)==E_R3(K(rr,1),jj+2)) ;
                    V2 = find(Table(1+VR1+3*VR2:T1,2)==E_R3(K(rr,1),jj+1));
                    V1 = find(Table(1+VR1+3*VR2:T1,2)==E_R3(K(rr,1),jj));

                    Jaco_ac(ii, VR1+3*VR2+V3) = Jaco_ac(ii, VR1+3*VR2+V3)+ eox*(Area(K(rr,1)+R_row2+R_row1,jj+2)/L(K(rr,1)+R_row2+R_row1,jj+2));
                    Jaco_ac(ii, VR1+3*VR2+V2) = Jaco_ac(ii, VR1+3*VR2+V2)+eox*(Area(K(rr,1)+R_row2+R_row1,jj)/L(K(rr,1)+R_row2+R_row1,jj));
                    Jaco_ac(ii, VR1+3*VR2+V1) = Jaco_ac(ii, VR1+3*VR2+V1)-eox*(Area(K(rr,1)+R_row2+R_row1,jj+2)/L(K(rr,1)+R_row2+R_row1,jj+2))-eox*(Area(K(rr,1)+R_row2+R_row1,jj)/L(K(rr,1)+R_row2+R_row1,jj));

                elseif K(rr,1) > R_row3 && K(rr,1) <= 2*R_row3

                    jj = 2;
                    K(rr,1) = K(rr,1)-R_row3;

                    V3 = find(Table(1+VR1+3*VR2:T1,2)==E_R3(K(rr,1),jj+1)) ;
                    V2 = find(Table(1+VR1+3*VR2:T1,2)==E_R3(K(rr,1),jj));
                    V1 = find(Table(1+VR1+3*VR2:T1,2)==E_R3(K(rr,1),jj-1));

                    Jaco_ac(ii, VR1+3*VR2+V3) = Jaco_ac(ii, VR1+3*VR2+V3)+eox*(Area(K(rr,1)+R_row2+R_row1,jj)/L(K(rr,1)+R_row2+R_row1,jj));
                    Jaco_ac(ii, VR1+3*VR2+V2) =  Jaco_ac(ii, VR1+3*VR2+V2)-eox*(Area(K(rr,1)+R_row2+R_row1,jj)/L(K(rr,1)+R_row2+R_row1,jj))-eox*(Area(K(rr,1)+R_row2+R_row1,jj-1)/L(K(rr,1)+R_row2+R_row1,jj-1));
                    Jaco_ac(ii, VR1+3*VR2+V1) = Jaco_ac(ii, VR1+3*VR2+V1)+eox*(Area(K(rr,1)+R_row2+R_row1,jj-1)/L(K(rr,1)+R_row2+R_row1,jj-1));

                elseif K(rr,1) > 2*R_row3 && K(rr,1) <= 3*R_row3

                    jj = 3;
                    K(rr,1) = K(rr,1)-2*R_row3;

                    V3 = find(Table(1+VR1+3*VR2:T1,2)==E_R3(K(rr,1),jj-2)) ;
                    V2 = find(Table(1+VR1+3*VR2:T1,2)==E_R3(K(rr,1),jj-1));
                    V1 = find(Table(1+VR1+3*VR2:T1,2)==E_R3(K(rr,1),jj));

                    Jaco_ac(ii, VR1+3*VR2+V3) = Jaco_ac(ii, VR1+3*VR2+V3) + eox*(Area(K(rr,1)+R_row2+R_row1,jj)/L(K(rr,1)+R_row2+R_row1,jj));
                    Jaco_ac(ii, VR1+3*VR2+V2) = Jaco_ac(ii, VR1+3*VR2+V2) + eox*(Area(K(rr,1)+R_row2+R_row1,jj-1)/L(K(rr,1)+R_row2+R_row1,jj-1));
                    Jaco_ac(ii, VR1+3*VR2+V1) = Jaco_ac(ii, VR1+3*VR2+V1) - eox*(Area(K(rr,1)+R_row2+R_row1,jj)/L(K(rr,1)+R_row2+R_row1,jj))-eox*(Area(K(rr,1)+R_row2+R_row1,jj-1)/L(K(rr,1)+R_row2+R_row1,jj-1));

                end
            end
        end
    end
end

%%%% Anode, Cathode %%%%
q_row = size(ry1,1);
w_row = size(ry2,1);

for rr = 1:q_row  % Cathode

    V = find(Table(1+VR1:3:VR1+3*VR2,2)==ry1(rr,1));
    Jaco_ac(3*V+(VR1-2),:) = 0;
    Jaco_ac(3*V+(VR1-2),3*V+(VR1-2)) = 1;
    Jaco_ac(3*V+(VR1-1),:) = 0;
    Jaco_ac(3*V+(VR1-1),3*V+(VR1-1)) = 1;
    Jaco_ac(3*V+(VR1),:) = 0;
    Jaco_ac(3*V+(VR1),3*V+(VR1)) = 1;

    res_ac(3*V+(VR1-2),:) = dphi(VR1+V,1)-0;
    res_ac(3*V+(VR1-1),:) = delec(VR1+V,1)-0;
    res_ac(3*V+(VR1),:) = dhole(VR1+V,1)-0;
end

for rr = 1:w_row  % Anode

    V = find(Table(1+VR1:3:VR1+3*VR2,2)==ry2(rr,1));

    Jaco_ac(3*V+(VR1-2),:) = 0;
    Jaco_ac(3*V+(VR1-2),3*V+(VR1-2)) = 1;
    Jaco_ac(3*V+(VR1-1),:) = 0;
    Jaco_ac(3*V+(VR1-1),3*V+(VR1-1)) = 1;
    Jaco_ac(3*V+(VR1),:) = 0;
    Jaco_ac(3*V+(VR1),3*V+(VR1)) = 1;

    res_ac(3*V+(VR1-2),:) = dphi(VR1+V,1)-Amp*i;
    res_ac(3*V+(VR1-1),:) = delec(VR1+V,1)-0;
    res_ac(3*V+(VR1),:) = dhole(VR1+V,1)-0;
end

%     scaling
Cvector = zeros(T_row,1);
Cvector(1:VR1,1) = Thermal_V;
Cvector(VR1+3*VR2:T_row,1) = Thermal_V;
Cvector(1+VR1:3:VR1+3*VR2,1) = Thermal_V;
Cvector(2+VR1:3:VR1+3*VR2,1) = -Na;
Cvector(3+VR1:3:VR1+3*VR2,1) = -Na;

Cmatrix = spdiags(Cvector,0,T_row,T_row);
Jaco_scaled = Jaco_ac * Cmatrix;
Rvector = 1./sum(abs(Jaco_scaled),2);
Rmatrix = spdiags(Rvector,0,T_row,T_row);
Jaco_scaled = Rmatrix* Jaco_scaled;
res_scaled = Rmatrix *res_ac;
sol_scaled = Jaco_scaled \ (res_scaled);
Sol_perturbed_AC(:,:) = Cmatrix* sol_scaled;

G_phi = Sol_perturbed_AC(:,1:3:T_row);  % Potential perturbation

G_elec = Sol_perturbed_AC(:,2:3:T_row);  % Electron perturbation

G_hole = Sol_perturbed_AC(:,3:3:T_row);  % Hole perturbation

% %%%%%% Current Calculation %%%%%
% w_row = size(ry2,1);
% 
% for tt=1:w_row
% 
%     K = find(E_R2==ry2(tt,1));
%     K_row = size(K,1);
% 
%     for rr = 1:K_row
% 
%         if K(rr,1) <= R_row2
%             jj = 1;
% 
%             V3 = find(Table(1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj+2)) ; % node 3
%             V2 = find(Table(1+1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj+1)); % node 2
%             V1 = find(Table(2+1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj)); %node 1
% 
%             elecDC_av1 = 0.5*(elec_DC(VR1+V2,1)+elec_DC(VR1+V1,1));
%             delec_av1 = 0.5*(delec(VR1+V2,1)+delec(VR1+V1,1));
%             ddelec1 = delec(VR1+V2,1)-delec(VR1+V1,1);
% 
%             elecDC_av2 = 0.5*(elec_DC(VR1+V3,1)+elec_DC(VR1+V1,1));
%             delec_av2 = 0.5*(delec(VR1+V3,1)+delec(VR1+V1,1));
%             ddelec2 = delec(VR1+V3,1)-delec(VR1+V1,1);
% 
%             holeDC_av1 = 0.5*(hole_DC(VR1+V2,1)+hole_DC(VR1+V1,1));
%             dhole_av1 = 0.5*(dhole(VR1+V2,1)+dhole(VR1+V1,1));
%             ddphi1 = dphi(VR1+V2,1)-dphi(VR1+V1,1);
%             ddhole1 = dhole(VR1+V2,1)-dhole(VR1+V1,1);
%             dphi_DC1 = phi_DC(VR1+V2,1)-phi_DC(VR1+V1,1);
% 
%             holeDC_av2 = 0.5*(hole_DC(VR1+V3,1)+hole_DC(VR1+V1,1));
%             dhole_av2 = 0.5*(dhole(VR1+V3,1)+dhole(VR1+V1,1));
%             ddphi2 = dphi(VR1+V3,1)-dphi(VR1+V1,1);
%             ddhole2 = dhole(VR1+V3,1)-dhole(VR1+V1,1);
%             dphi_DC2 = phi_DC(VR1+V3,1)-phi_DC(VR1+V1,1);
% 
%             Jn(1,:) = Jn(1,:) + coeff_dJn*((elecDC_av1*ddphi1+delec_av1*dphi_DC1)-Thermal_V*ddelec1)*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj));
%             Jn(1,:) = Jn(1,:) + coeff_dJn*((elecDC_av2*ddphi2+delec_av2*dphi_DC2)-Thermal_V*ddelec2)*(Area(K(rr,1)+R_row1,jj+2)/L(K(rr,1)+R_row1,jj+2));
%             Jp(1,:) = Jp(1,:) + coeff_dJp*((holeDC_av1*ddphi1+dhole_av1*dphi_DC1)+Thermal_V*ddhole1)*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj));
%             Jp(1,:) = Jp(1,:) + coeff_dJp*((holeDC_av2*ddphi2+dhole_av2*dphi_DC2)+Thermal_V*ddhole2)*(Area(K(rr,1)+R_row1,jj+2)/L(K(rr,1)+R_row1,jj+2));
%             Jd(1,:) = Jd(1,:) - esi*e0*(2*pi*freq*i)*(dphi(VR1+V3,1)-dphi(VR1+V1,1))*Area(K(rr,1)+R_row1,jj+2)/L(K(rr,1)+R_row1,jj+2);
%             Jd(1,:) = Jd(1,:) - esi*e0*(2*pi*freq*i)*(dphi(VR1+V2,1)-dphi(VR1+V1,1))*Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj);
% 
%         elseif K(rr,1) > R_row2 && K(rr,1) <= 2*R_row2
% 
%             jj = 2;
%             K(rr,1) = K(rr,1)-R_row2;
% 
%             V3 = find(Table(1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj+1)) ; % node 3
%             V2 = find(Table(1+1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj)); % node 2
%             V1 = find(Table(2+1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj-1)); %node 1
% 
%             elecDC_av1 = 0.5*(elec_DC(VR1+V1,1)+elec_DC(VR1+V2,1));
%             delec_av1 = 0.5*(delec(VR1+V1,1)+delec(VR1+V2,1));
%             ddelec1 = delec(VR1+V1,1)-delec(VR1+V2,1);
% 
%             elecDC_av2 = 0.5*(elec_DC(VR1+V3,1)+elec_DC(VR1+V2,1));
%             delec_av2 = 0.5*(delec(VR1+V3,1)+delec(VR1+V2,1));
%             ddelec2 = delec(VR1+V3,1)-delec(VR1+V2,1);
% 
%             holeDC_av1 = 0.5*(hole_DC(VR1+V1,1)+hole_DC(VR1+V2,1));
%             dhole_av1 = 0.5*(dhole(VR1+V1,1)+dhole(VR1+V2,1));
%             ddphi1 = dphi(VR1+V1,1)-dphi(VR1+V2,1);
%             ddhole1 = dhole(VR1+V1,1)-dhole(VR1+V2,1);
%             dphi_DC1 = phi_DC(VR1+V1,1)-phi_DC(VR1+V2,1);
% 
%             holeDC_av2 = 0.5*(hole_DC(VR1+V3,1)+hole_DC(VR1+V2,1));
%             dhole_av2 = 0.5*(dhole(VR1+V3,1)+dhole(VR1+V2,1));
%             ddphi2 = dphi(VR1+V3,1)-dphi(VR1+V2,1);
%             ddhole2 = dhole(VR1+V3,1)-dhole(VR1+V2,1);
%             dphi_DC2 = phi_DC(VR1+V3,1)-phi_DC(VR1+V2,1);
% 
%             Jn(1,:) = Jn(1,:) + coeff_dJn*((Area(K(rr,1)+R_row1,jj-1)/L(K(rr,1)+R_row1,jj-1))*((elecDC_av1*ddphi1+delec_av1*dphi_DC1)-Thermal_V*ddelec1));
%             Jn(1,:) = Jn(1,:) + coeff_dJn*((Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj))*((elecDC_av2*ddphi2+delec_av2*dphi_DC2)-Thermal_V*ddelec2));
%             Jp(1,:) = Jp(1,:) + coeff_dJp*((Area(K(rr,1)+R_row1,jj-1)/L(K(rr,1)+R_row1,jj-1))*((holeDC_av1*ddphi1+dhole_av1*dphi_DC1)+Thermal_V*ddhole1));
%             Jp(1,:) = Jp(1,:) + coeff_dJp*((Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj))*((holeDC_av2*ddphi2+dhole_av2*dphi_DC2)+Thermal_V*ddhole2));
%             Jd(1,:) = Jd(1,:) - esi*e0*(2*pi*freq*i)*(dphi(VR1+V1,1)-dphi(VR1+V2,1))*Area(K(rr,1)+R_row1,jj-1)/L(K(rr,1)+R_row1,jj-1);
%             Jd(1,:) = Jd(1,:) - esi*e0*(2*pi*freq*i)*(dphi(VR1+V3,1)-dphi(VR1+V2,1))*Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj);
% 
%         elseif K(rr,1) > 2*R_row2 && K(rr,1) <= 3*R_row2
% 
%             jj = 3;
%             K(rr,1) = K(rr,1)-2*R_row2;
% 
%             V3 = find(Table(1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj)) ; % node 3
%             V2 = find(Table(1+1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj-1)); % node 2
%             V1 = find(Table(2+1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj-2)); %node 1
% 
%             elecDC_av1 = 0.5*(elec_DC(VR1+V1,1)+elec_DC(VR1+V3,1));
%             delec_av1 = 0.5*(delec(VR1+V1,1)+delec(VR1+V3,1));
%             ddelec1 = delec(VR1+V1,1)-delec(VR1+V3,1);
% 
%             elecDC_av2 = 0.5*(elec_DC(VR1+V2,1)+elec_DC(VR1+V3,1));
%             delec_av2 = 0.5*(delec(VR1+V2,1)+delec(VR1+V3,1));
%             ddelec2 = delec(VR1+V2,1)-delec(VR1+V3,1);
% 
%             holeDC_av1 = 0.5*(hole_DC(VR1+V1,1)+hole_DC(VR1+V3,1));
%             dhole_av1 = 0.5*(dhole(VR1+V1,1)+dhole(VR1+V3,1));
%             ddphi1 = dphi(VR1+V1,1)-dphi(VR1+V3,1);
%             ddhole1 = dhole(VR1+V1,1)-dhole(VR1+V3,1);
%             dphi_DC1 = phi_DC(VR1+V1,1)-phi_DC(VR1+V3,1);
% 
%             holeDC_av2 = 0.5*(hole_DC(VR1+V2,1)+hole_DC(VR1+V3,1));
%             dhole_av2 = 0.5*(dhole(VR1+V2,1)+dhole(VR1+V3,1));
%             ddphi2 = dphi(VR1+V2,1)-dphi(VR1+V3,1);
%             ddhole2 = dhole(VR1+V2,1)-dhole(VR1+V3,1);
%             dphi_DC2 = phi_DC(VR1+V2,1)-phi_DC(VR1+V3,1);
% 
%             Jn(1,:) = Jn(1,:) + coeff_dJn*((Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj))*((elecDC_av1*ddphi1+delec_av1*dphi_DC1)-Thermal_V*ddelec1));
%             Jn(1,:) = Jn(1,:) + coeff_dJn*((Area(K(rr,1)+R_row1,jj-1)/L(K(rr,1)+R_row1,jj-1))*((elecDC_av2*ddphi2+delec_av2*dphi_DC2)-Thermal_V*ddelec2));
%             Jp(1,:) = Jp(1,:) + coeff_dJp*((Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj))*((holeDC_av1*ddphi1+dhole_av1*dphi_DC1)+Thermal_V*ddhole1));
%             Jp(1,:) = Jp(1,:) + coeff_dJp*((Area(K(rr,1)+R_row1,jj-1)/L(K(rr,1)+R_row1,jj-1))*((holeDC_av2*ddphi2+dhole_av2*dphi_DC2)+Thermal_V*ddhole2));
%             Jd(1,:) = Jd(1,:) - esi*e0*(2*pi*freq*i)*(dphi(VR1+V1,1)-dphi(VR1+V3,1))*Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj);
%             Jd(1,:) = Jd(1,:) - esi*e0*(2*pi*freq*i)*(dphi(VR1+V2,1)-dphi(VR1+V3,1))*Area(K(rr,1)+R_row1,jj-1)/L(K(rr,1)+R_row1,jj-1);
%         end
%         I = (Jn+Jp+Jd)*1e-6;
%     end
% end
