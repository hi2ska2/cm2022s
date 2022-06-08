clear; close all; clc;
load('Double_gate_equil.mat')

%%% due to : 2022.06.09. / HW25

%%%%% HW 25 %%%%%%%
%%% Multi-Mode Simulation

%%%%%% Drift-Diffusion %%%%%%%%
parameter = 13;
T1 = VR1+3*VR2+VR3;
Table = zeros(T1+parameter,3);

T_row = size(Table,1);

for iRegion = 1:4
    if iRegion ==1  % Oxide , only potential

        Table(1:VR1,1) = 1;

        for iVertex =  1:VR1
            Table(iVertex,2) = V_R1(iVertex,1);
            Table(iVertex,3) = Potential;   % only Potential
        end

    elseif iRegion == 2     % Silicon, Potential, eDensity, hDensity

        Table(1+VR1:VR1+3*VR2,1) = 2;

        for iVertex =  1:VR2

            Table(3*iVertex+(VR1-2),2) = V_R2(iVertex,1);
            Table(3*iVertex+(VR1-1),2) = V_R2(iVertex,1);
            Table(3*iVertex+(VR1),2) = V_R2(iVertex,1);

            Table(3*iVertex+(VR1-2),3) = Potential;
            Table(3*iVertex+(VR1-1),3) = eDensity;
            Table(3*iVertex+(VR1),3) = hDensity;
        end

    elseif iRegion ==3  % Oxide , only potential

        Table(1+VR1+3*VR2:T1,1) = 3;
        
        for iVertex = 1:VR3
            Table(VR1+3*VR2+iVertex,2) = V_R3(iVertex,1);
            Table(VR1+3*VR2+iVertex,3) = Potential;
        end

    elseif iRegion == 4  % Oxide , only potential

        Table(1+T1:T1+parameter,1) = 4;

        for iVertex = 1:parameter
            Table(1+T1:T1+parameter,3) = 4;
        end
    end
end

for Newton = 1

    Jaco = zeros(T_row,T_row);
    res = zeros(T_row,1);

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

                        res(ii,1) = res(ii,1) + eox*(phi(V2,1)-phi(V1,1))*(Area(K(rr,1),jj)/L(K(rr,1),jj)) + eox*(phi(V3,1)-phi(V1,1))*(Area(K(rr,1),jj+2)/L(K(rr,1),jj+2));
                        Jaco(ii, V3) = Jaco(ii, V3)+ eox*(Area(K(rr,1),jj+2)/L(K(rr,1),jj+2));
                        Jaco(ii, V2) = Jaco(ii, V2)+eox*(Area(K(rr,1),jj)/L(K(rr,1),jj));
                        Jaco(ii, V1) = Jaco(ii, V1)-eox*(Area(K(rr,1),jj+2)/L(K(rr,1),jj+2))-eox*(Area(K(rr,1),jj)/L(K(rr,1),jj));

                    elseif K(rr,1) > R_row1 && K(rr,1) <= 2*R_row1
                        jj = 2;
                        K(rr,1) = K(rr,1)-R_row1;

                        V3 = find(Table(1:VR1,2)==E_R1(K(rr,1),jj+1)) ;
                        V2 = find(Table(1:VR1,2)==E_R1(K(rr,1),jj));
                        V1 = find(Table(1:VR1,2)==E_R1(K(rr,1),jj-1));

                        res(ii,1) = res(ii,1) + eox*(phi(V1,1)-phi(V2,1))*(Area(K(rr,1),jj-1)/L(K(rr,1),jj-1)) + eox*(phi(V3,1)-phi(V2,1))*(Area(K(rr,1),jj)/L(K(rr,1),jj));
                        Jaco(ii, V3) = Jaco(ii, V3)+eox*(Area(K(rr,1),jj)/L(K(rr,1),jj));
                        Jaco(ii, V2) =  Jaco(ii, V2)-eox*(Area(K(rr,1),jj)/L(K(rr,1),jj))-eox*(Area(K(rr,1),jj-1)/L(K(rr,1),jj-1));
                        Jaco(ii, V1) = Jaco(ii, V1)+eox*(Area(K(rr,1),jj-1)/L(K(rr,1),jj-1));

                    elseif K(rr,1) > 2*R_row1 && K(rr,1) <= 3*R_row1
                        jj = 3;
                        K(rr,1) = K(rr,1)-2*R_row1;

                        V3 = find(Table(1:VR1,2)==E_R1(K(rr,1),jj-2)) ;
                        V2 = find(Table(1:VR1,2)==E_R1(K(rr,1),jj-1));
                        V1 = find(Table(1:VR1,2)==E_R1(K(rr,1),jj));

                        res(ii,1) = res(ii,1) + eox*(phi(V3,1)-phi(V1,1))*(Area(K(rr,1),jj)/L(K(rr,1),jj)) + eox*(phi(V2,1)-phi(V1,1))*(Area(K(rr,1),jj-1)/L(K(rr,1),jj-1));
                        Jaco(ii, V3) = Jaco(ii, V3)+eox*(Area(K(rr,1),jj)/L(K(rr,1),jj));
                        Jaco(ii, V2) = Jaco(ii, V2)+eox*(Area(K(rr,1),jj-1)/L(K(rr,1),jj-1));
                        Jaco(ii, V1) = Jaco(ii, V1)-eox*(Area(K(rr,1),jj)/L(K(rr,1),jj))-eox*(Area(K(rr,1),jj-1)/L(K(rr,1),jj-1));

                    end
                end
            end

        elseif Table(ii,1) == 2     %region2

            for kk = 1:rr_row
                if Table(ii,2) >= iy1(kk,1) && Table(ii,2) <= iy2(kk,1)

                    K = find(E_R2 == Table(ii,2));
                    K_row = size(K,1);


                    if Table(ii,3) == Potential

                        for rr = 1:K_row

                            if K(rr,1) <= R_row2
                                jj = 1;

                                V3 = find(Table(1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj+2)) ;
                                V2 = find(Table(1+1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj+1));
                                V1 = find(Table(1+1+1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj));

                                res(ii,1) = res(ii,1) + esi*(phi(VR1+V2,1)-phi(VR1+V1,1))*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj)) + esi*(phi(VR1+V3,1)-phi(VR1+V1,1))*(Area(K(rr,1)+R_row1,jj+2)/L(K(rr,1)+R_row1,jj+2));
                                res(ii,1) = res(ii,1) + coeff*(Area(K(rr,1)+R_row1,jj)*L(K(rr,1)+R_row1,jj)+Area(K(rr,1)+R_row1,jj+2)*L(K(rr,1)+R_row1,jj+2))*(-Na-elec(VR1+V1,1)+hole(VR1+V1,1));

                                Jaco(ii, 3*(V3)+(VR1-2)) = Jaco(ii, 3*(V3)+(VR1-2))+ esi*(Area(K(rr,1)+R_row1,jj+2)/L(K(rr,1)+R_row1,jj+2));
                                Jaco(ii, 3*(V2)+(VR1-2)) = Jaco(ii, 3*(V2)+(VR1-2))+esi*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj));
                                Jaco(ii, 3*(V1)+(VR1-2)) = Jaco(ii, 3*(V1)+(VR1-2))-esi*(Area(K(rr,1)+R_row1,jj+2)/L(K(rr,1)+R_row1,jj+2))-esi*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj));

                                Jaco(ii, 3*(V1)+(VR1-1)) = Jaco(ii, 3*(V1)+(VR1-1)) - coeff*(Area(K(rr,1)+R_row1,jj)*L(K(rr,1)+R_row1,jj)+Area(K(rr,1)+R_row1,jj+2)*L(K(rr,1)+R_row1,jj+2));
                                Jaco(ii, 3*(V1)+(VR1)) = Jaco(ii, 3*(V1)+(VR1)) + coeff*(Area(K(rr,1)+R_row1,jj)*L(K(rr,1)+R_row1,jj)+Area(K(rr,1)+R_row1,jj+2)*L(K(rr,1)+R_row1,jj+2));

                            elseif K(rr,1) > R_row2 && K(rr,1) <= 2*R_row2
                                jj = 2;
                                K(rr,1) = K(rr,1)-R_row2;

                                V3 = find(Table(1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj+1)) ;
                                V2 = find(Table(1+1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj));
                                V1 = find(Table(1+1+1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj-1));

                                res(ii,1) = res(ii,1) + esi*(phi(VR1+V1,1)-phi(VR1+V2,1))*(Area(K(rr,1)+R_row1,jj-1)/L(K(rr,1)+R_row1,jj-1)) + esi*(phi(VR1+V3,1)-phi(VR1+V2,1))*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj));
                                res(ii,1) = res(ii,1) + coeff*(Area(K(rr,1)+R_row1,jj)*L(K(rr,1)+R_row1,jj)+Area(K(rr,1)+R_row1,jj-1)*L(K(rr,1)+R_row1,jj-1))*(-Na-elec(VR1+V2,1)+hole(VR1+V2,1));

                                Jaco(ii, 3*(V3)+(VR1-2)) = Jaco(ii, 3*(V3)+(VR1-2))+esi*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj));
                                Jaco(ii, 3*(V2)+(VR1-2)) =  Jaco(ii, 3*(V2)+(VR1-2))-esi*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj))-esi*(Area(K(rr,1)+R_row1,jj-1)/L(K(rr,1)+R_row1,jj-1));
                                Jaco(ii, 3*(V1)+(VR1-2)) = Jaco(ii, 3*(V1)+(VR1-2))+esi*(Area(K(rr,1)+R_row1,jj-1)/L(K(rr,1)+R_row1,jj-1));

                                Jaco(ii, 3*(V2)+(VR1-1)) = Jaco(ii, 3*(V2)+(VR1-1)) - coeff*(Area(K(rr,1)+R_row1,jj)*L(K(rr,1)+R_row1,jj)+Area(K(rr,1)+R_row1,jj-1)*L(K(rr,1)+R_row1,jj-1));
                                Jaco(ii, 3*(V2)+(VR1)) = Jaco(ii, 3*(V2)+(VR1)) + coeff*(Area(K(rr,1)+R_row1,jj)*L(K(rr,1)+R_row1,jj)+Area(K(rr,1)+R_row1,jj-1)*L(K(rr,1)+R_row1,jj-1));

                            elseif K(rr,1) > 2*R_row2 && K(rr,1) <= 3*R_row2
                                jj = 3;
                                K(rr,1) = K(rr,1)-2*R_row2;

                                V3 = find(Table(1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj-2)) ;
                                V2 = find(Table(1+1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj-1));
                                V1 = find(Table(1+1+1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj));

                                res(ii,1) = res(ii,1) + esi*(phi(VR1+V3,1)-phi(VR1+V1,1))*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj)) + esi*(phi(VR1+V2,1)-phi(VR1+V1,1))*(Area(K(rr,1)+R_row1,jj-1)/L(K(rr,1)+R_row1,jj-1));
                                res(ii,1) = res(ii,1) + coeff*(Area(K(rr,1)+R_row1,jj)*L(K(rr,1)+R_row1,jj)+Area(K(rr,1)+R_row1,jj-1)*L(K(rr,1)+R_row1,jj-1))*(-Na-elec(VR1+V1,1)+hole(VR1+V1,1));

                                Jaco(ii, 3*(V3)+(VR1-2)) = Jaco(ii, 3*(V3)+(VR1-2))+esi*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj));
                                Jaco(ii, 3*(V2)+(VR1-2)) = Jaco(ii, 3*(V2)+(VR1-2))+esi*(Area(K(rr,1)+R_row1,jj-1)/L(K(rr,1)+R_row1,jj-1));
                                Jaco(ii, 3*(V1)+(VR1-2)) = Jaco(ii, 3*(V1)+(VR1-2))-esi*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj))-esi*(Area(K(rr,1)+R_row1,jj-1)/L(K(rr,1)+R_row1,jj-1));

                                Jaco(ii, 3*(V1)+(VR1-1)) = Jaco(ii, 3*(V1)+(VR1-1)) - coeff*(Area(K(rr,1)+R_row1,jj)*L(K(rr,1)+R_row1,jj)+Area(K(rr,1)+R_row1,jj-1)*L(K(rr,1)+R_row1,jj-1));
                                Jaco(ii, 3*(V1)+(VR1)) = Jaco(ii, 3*(V1)+(VR1)) + coeff*(Area(K(rr,1)+R_row1,jj)*L(K(rr,1)+R_row1,jj)+Area(K(rr,1)+R_row1,jj-1)*L(K(rr,1)+R_row1,jj-1));
                            end
                        end

                    elseif Table(ii,3) == eDensity

                        for rr = 1:K_row
                            if K(rr,1) <= R_row2
                                jj = 1;

                                V3 = find(Table(1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj+2)) ;  % node 3
                                V2 = find(Table(1+1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj+1)); % node 2
                                V1 = find(Table(2+1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj));   % node 1

                                x21 = (phi(VR1+V2,1)-phi(VR1+V1,1))/Thermal_V;
                                x12 = -x21;
                                x31 = (phi(VR1+V3,1)-phi(VR1+V1,1))/Thermal_V;
                                x13 = -x31;

                                res(ii,1) = res(ii,1) + coeff_Jn*((Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj))*(elec(VR1+V2,1)*Ber(x21) - elec(VR1+V1,1)*Ber(x12)));
                                res(ii,1) = res(ii,1) + coeff_Jn*((Area(K(rr,1)+R_row1,jj+2)/L(K(rr,1)+R_row1,jj+2))*(elec(VR1+V3,1)*Ber(x31) - elec(VR1+V1,1)*Ber(x13)));

                                % eDensity
                                Jaco(ii, 3*(V1)+(VR1-1)) = Jaco(ii, 3*(V1)+(VR1-1)) + coeff_Jn*-Ber(x12)*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj));
                                Jaco(ii, 3*(V1)+(VR1-1)) = Jaco(ii, 3*(V1)+(VR1-1)) + coeff_Jn*(-Ber(x13))*(Area(K(rr,1)+R_row1,jj+2)/L(K(rr,1)+R_row1,jj+2));

                                Jaco(ii, 3*(V2)+(VR1-1)) = Jaco(ii, 3*(V2)+(VR1-1)) + coeff_Jn*Ber(x21)*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj));
                                Jaco(ii, 3*(V3)+(VR1-1)) = Jaco(ii, 3*(V3)+(VR1-1)) + coeff_Jn*Ber(x31)*(Area(K(rr,1)+R_row1,jj+2)/L(K(rr,1)+R_row1,jj+2));

                                % Potential
                                Jaco(ii, 3*(V1)+(VR1-2)) = Jaco(ii, 3*(V1)+(VR1-2)) + coeff_Jn/Thermal_V*(-elec(VR1+V2,1)*Ber_d(x21)-elec(VR1+V1,1)*Ber_d(x12))*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj));
                                Jaco(ii, 3*(V1)+(VR1-2)) = Jaco(ii, 3*(V1)+(VR1-2)) + coeff_Jn/Thermal_V*(-elec(VR1+V3,1)*Ber_d(x31)-elec(VR1+V1,1)*Ber_d(x13))*(Area(K(rr,1)+R_row1,jj+2)/L(K(rr,1)+R_row1,jj+2));

                                Jaco(ii, 3*(V2)+(VR1-2)) = Jaco(ii, 3*(V2)+(VR1-2)) + coeff_Jn/Thermal_V*(elec(VR1+V2,1)*Ber_d(x21)+elec(VR1+V1,1)*Ber_d(x12))*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj));
                                Jaco(ii, 3*(V3)+(VR1-2)) = Jaco(ii, 3*(V3)+(VR1-2)) + coeff_Jn/Thermal_V*(elec(VR1+V3,1)*Ber_d(x31)+elec(VR1+V1,1)*Ber_d(x13))*(Area(K(rr,1)+R_row1,jj+2)/L(K(rr,1)+R_row1,jj+2));


                            elseif K(rr,1) > R_row2 && K(rr,1) <= 2*R_row2
                                jj = 2;
                                K(rr,1) = K(rr,1)-R_row2;

                                V3 = find(Table(1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj+1)) ; % node 3
                                V2 = find(Table(1+1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj));  % node 2
                                V1 = find(Table(2+1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj-1)); % node 1

                                x12 = (phi(VR1+V1,1)-phi(VR1+V2,1))/Thermal_V;
                                x21 = -x12;
                                x32 = (phi(VR1+V3,1)-phi(VR1+V2,1))/Thermal_V;
                                x23 = -x32;

                                res(ii,1) =  res(ii,1) + coeff_Jn*((Area(K(rr,1)+R_row1,jj-1)/L(K(rr,1)+R_row1,jj-1))*(elec(VR1+V1,1)*Ber(x12) - elec(VR1+V2,1)*Ber(x21)));
                                res(ii,1) =  res(ii,1) + coeff_Jn*((Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj))*(elec(VR1+V3,1)*Ber(x32) - elec(VR1+V2,1)*Ber(x23)));

                                % eDensity
                                Jaco(ii, 3*(V2)+(VR1-1)) = Jaco(ii, 3*(V2)+(VR1-1)) + coeff_Jn*(-Ber(x21))*(Area(K(rr,1)+R_row1,jj-1)/L(K(rr,1)+R_row1,jj-1));
                                Jaco(ii, 3*(V2)+(VR1-1)) = Jaco(ii, 3*(V2)+(VR1-1)) + coeff_Jn*(-Ber(x23))*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj));

                                Jaco(ii, 3*(V1)+(VR1-1)) = Jaco(ii, 3*(V1)+(VR1-1)) + coeff_Jn*Ber(x12)*(Area(K(rr,1)+R_row1,jj-1)/L(K(rr,1)+R_row1,jj-1));
                                Jaco(ii, 3*(V3)+(VR1-1)) = Jaco(ii, 3*(V3)+(VR1-1)) + coeff_Jn*Ber(x32)*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj));

                                % Potential
                                Jaco(ii, 3*(V2)+(VR1-2)) = Jaco(ii, 3*(V2)+(VR1-2)) + coeff_Jn/Thermal_V*(-elec(VR1+V1,1)*Ber_d(x12)-elec(VR1+V2,1)*Ber_d(x21))*(Area(K(rr,1)+R_row1,jj-1)/L(K(rr,1)+R_row1,jj-1));
                                Jaco(ii, 3*(V2)+(VR1-2)) = Jaco(ii, 3*(V2)+(VR1-2)) + coeff_Jn/Thermal_V*(-elec(VR1+V3,1)*Ber_d(x32)-elec(VR1+V2,1)*Ber_d(x23))*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj));

                                Jaco(ii, 3*(V1)+(VR1-2)) = Jaco(ii, 3*(V1)+(VR1-2)) + coeff_Jn/Thermal_V*(elec(VR1+V1,1)*Ber_d(x12)+elec(VR1+V2,1)*Ber_d(x21))*(Area(K(rr,1)+R_row1,jj-1)/L(K(rr,1)+R_row1,jj-1));
                                Jaco(ii, 3*(V3)+(VR1-2)) = Jaco(ii, 3*(V3)+(VR1-2)) + coeff_Jn/Thermal_V*(elec(VR1+V3,1)*Ber_d(x32)+elec(VR1+V2,1)*Ber_d(x23))*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj));

                            elseif K(rr,1) > 2*R_row2 && K(rr,1) <= 3*R_row2
                                jj = 3;
                                K(rr,1) = K(rr,1)-2*R_row2;

                                V3 = find(Table(2+1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj)); % node 3
                                V2 = find(Table(1+1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj-1)); % node 2
                                V1 = find(Table(1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj-2)) ;  % node 1

                                x13 = (phi(VR1+V1,1)-phi(VR1+V3,1))/Thermal_V;
                                x31 = -x13;
                                x23 = (phi(VR1+V2,1)-phi(VR1+V3,1))/Thermal_V;
                                x32 = -x23;

                                res(ii,1) = res(ii,1) + coeff_Jn*(elec(VR1+V1,1)*Ber(x13) - elec(VR1+V3,1)*Ber(x31))*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj));
                                res(ii,1) = res(ii,1) + coeff_Jn*(elec(VR1+V2,1)*Ber(x23) - elec(VR1+V3,1)*Ber(x32))*(Area(K(rr,1)+R_row1,jj-1)/L(K(rr,1)+R_row1,jj-1));

                                % eDensity
                                Jaco(ii, 3*(V1)+(VR1-1)) = Jaco(ii, 3*(V1)+(VR1-1)) + coeff_Jn*Ber(x13)*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj));
                                Jaco(ii, 3*(V2)+(VR1-1)) = Jaco(ii, 3*(V2)+(VR1-1)) + coeff_Jn*Ber(x23)*(Area(K(rr,1)+R_row1,jj-1)/L(K(rr,1)+R_row1,jj-1));

                                Jaco(ii, 3*(V3)+(VR1-1)) = Jaco(ii, 3*(V3)+(VR1-1)) + coeff_Jn*(-Ber(x31))*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj));
                                Jaco(ii, 3*(V3)+(VR1-1)) = Jaco(ii, 3*(V3)+(VR1-1)) + coeff_Jn*(-Ber(x32))*(Area(K(rr,1)+R_row1,jj-1)/L(K(rr,1)+R_row1,jj-1));

                                % Potential
                                Jaco(ii, 3*(V3)+(VR1-2)) = Jaco(ii, 3*(V3)+(VR1-2)) + coeff_Jn/Thermal_V*(-elec(VR1+V1,1)*Ber_d(x13)-elec(VR1+V3,1)*Ber_d(x31))*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj));
                                Jaco(ii, 3*(V3)+(VR1-2)) = Jaco(ii, 3*(V3)+(VR1-2)) + coeff_Jn/Thermal_V*(-elec(VR1+V2,1)*Ber_d(x23)-elec(VR1+V3,1)*Ber_d(x32))*(Area(K(rr,1)+R_row1,jj-1)/L(K(rr,1)+R_row1,jj-1));

                                Jaco(ii, 3*(V2)+(VR1-2)) = Jaco(ii, 3*(V2)+(VR1-2)) + coeff_Jn/Thermal_V*(elec(VR1+V2,1)*Ber_d(x23)+elec(VR1+V3,1)*Ber_d(x32))*(Area(K(rr,1)+R_row1,jj-1)/L(K(rr,1)+R_row1,jj-1));
                                Jaco(ii, 3*(V1)+(VR1-2)) = Jaco(ii, 3*(V1)+(VR1-2)) + coeff_Jn/Thermal_V*(elec(VR1+V1,1)*Ber_d(x13)+elec(VR1+V3,1)*Ber_d(x31))*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj));
                            end
                        end


                    elseif Table(ii,3) == hDensity

                        for rr = 1:K_row

                            if K(rr,1) <= R_row2
                                jj = 1;

                                V3 = find(Table(1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj+2)) ;  % node 3
                                V2 = find(Table(1+1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj+1));  % node 2
                                V1 = find(Table(2+1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj));  % node 1

                                x21 = (phi(VR1+V2,1)-phi(VR1+V1,1))/Thermal_V;
                                x12 = -x21;
                                x31 = (phi(VR1+V3,1)-phi(VR1+V1,1))/Thermal_V;
                                x13 = -x31;

                                res(ii,1) = res(ii,1) + coeff_Jp*(hole(VR1+V2,1)*Ber(x12) - hole(VR1+V1,1)*Ber(x21))*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj));
                                res(ii,1) = res(ii,1) + coeff_Jp*(hole(VR1+V3,1)*Ber(x13) - hole(VR1+V1,1)*Ber(x31))*(Area(K(rr,1)+R_row1,jj+2)/L(K(rr,1)+R_row1,jj+2));

                                %hDensity
                                Jaco(ii, 3*(V1)+(VR1)) = Jaco(ii, 3*(V1)+(VR1)) + coeff_Jp*(-Ber(x21))*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj));
                                Jaco(ii, 3*(V1)+(VR1)) = Jaco(ii, 3*(V1)+(VR1)) + coeff_Jp*(-Ber(x31))*(Area(K(rr,1)+R_row1,jj+2)/L(K(rr,1)+R_row1,jj+2));

                                Jaco(ii, 3*(V2)+(VR1)) = Jaco(ii, 3*(V2)+(VR1)) + coeff_Jp*Ber(x12)*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj));
                                Jaco(ii, 3*(V3)+(VR1)) = Jaco(ii, 3*(V3)+(VR1)) + coeff_Jp*Ber(x13)*(Area(K(rr,1)+R_row1,jj+2)/L(K(rr,1)+R_row1,jj+2));

                                %Potential
                                Jaco(ii, 3*(V1)+(VR1-2)) = Jaco(ii, 3*(V1)+(VR1-2)) + coeff_Jp/Thermal_V*(hole(VR1+V2,1)*Ber_d(x12)+hole(VR1+V1,1)*Ber_d(x21))*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj));
                                Jaco(ii, 3*(V1)+(VR1-2)) = Jaco(ii, 3*(V1)+(VR1-2)) + coeff_Jp/Thermal_V*(hole(VR1+V3,1)*Ber_d(x13)+hole(VR1+V1,1)*Ber_d(x31))*(Area(K(rr,1)+R_row1,jj+2)/L(K(rr,1)+R_row1,jj+2));

                                Jaco(ii, 3*(V2)+(VR1-2)) = Jaco(ii, 3*(V2)+(VR1-2)) + coeff_Jp/Thermal_V*(-hole(VR1+V2,1)*Ber_d(x12)-hole(VR1+V1,1)*Ber_d(x21))*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj));
                                Jaco(ii, 3*(V3)+(VR1-2)) = Jaco(ii, 3*(V3)+(VR1-2)) + coeff_Jp/Thermal_V*(-hole(VR1+V3,1)*Ber_d(x13)-hole(VR1+V1,1)*Ber_d(x31))*(Area(K(rr,1)+R_row1,jj+2)/L(K(rr,1)+R_row1,jj+2));


                            elseif K(rr,1) > R_row2 && K(rr,1) <= 2*R_row2
                                jj = 2;
                                K(rr,1) = K(rr,1)-R_row2;

                                V3 = find(Table(1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj+1)) ;  % node 3
                                V2 = find(Table(1+1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj));  % node 2
                                V1 = find(Table(2+1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj-1));  % node 1

                                x12 = (phi(VR1+V1,1)-phi(VR1+V2,1))/Thermal_V;
                                x21 = -x12;
                                x32 = (phi(VR1+V3,1)-phi(VR1+V2,1))/Thermal_V;
                                x23 = -x32;

                                res(ii,1) =  res(ii,1) + coeff_Jp*((Area(K(rr,1)+R_row1,jj-1)/L(K(rr,1)+R_row1,jj-1))*(hole(VR1+V1,1)*Ber(x21) - hole(VR1+V2,1)*Ber(x12)));
                                res(ii,1) =  res(ii,1) + coeff_Jp*((Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj))*(hole(VR1+V3,1)*Ber(x23) - hole(VR1+V2,1)*Ber(x32)));

                                % hDensity
                                Jaco(ii, 3*(V1)+(VR1)) = Jaco(ii, 3*(V1)+(VR1)) + coeff_Jp*Ber(x21)*(Area(K(rr,1)+R_row1,jj-1)/L(K(rr,1)+R_row1,jj-1));
                                Jaco(ii, 3*(V3)+(VR1)) = Jaco(ii, 3*(V3)+(VR1)) + coeff_Jp*Ber(x23)*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj));

                                Jaco(ii, 3*(V2)+(VR1)) = Jaco(ii, 3*(V2)+(VR1)) + coeff_Jp*(-Ber(x12))*(Area(K(rr,1)+R_row1,jj-1)/L(K(rr,1)+R_row1,jj-1));
                                Jaco(ii, 3*(V2)+(VR1)) = Jaco(ii, 3*(V2)+(VR1)) + coeff_Jp*(-Ber(x32))*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj));


                                %Potential
                                Jaco(ii, 3*(V2)+(VR1-2)) = Jaco(ii, 3*(V2)+(VR1-2)) + coeff_Jp/Thermal_V*(hole(VR1+V1,1)*Ber_d(x21)+hole(VR1+V2,1)*Ber_d(x12))*(Area(K(rr,1)+R_row1,jj-1)/L(K(rr,1)+R_row1,jj-1));
                                Jaco(ii, 3*(V2)+(VR1-2)) = Jaco(ii, 3*(V2)+(VR1-2)) + coeff_Jp/Thermal_V*(hole(VR1+V3,1)*Ber_d(x23)+hole(VR1+V2,1)*Ber_d(x32))*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj));

                                Jaco(ii, 3*(V1)+(VR1-2)) = Jaco(ii, 3*(V1)+(VR1-2)) + coeff_Jp/Thermal_V*(-hole(VR1+V1,1)*Ber_d(x21)-hole(VR1+V2,1)*Ber_d(x12))*(Area(K(rr,1)+R_row1,jj-1)/L(K(rr,1)+R_row1,jj-1));
                                Jaco(ii, 3*(V3)+(VR1-2)) = Jaco(ii, 3*(V3)+(VR1-2)) + coeff_Jp/Thermal_V*(-hole(VR1+V3,1)*Ber_d(x23)-hole(VR1+V2,1)*Ber_d(x32))*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj));


                            elseif K(rr,1) > 2*R_row2 && K(rr,1) <= 3*R_row2
                                jj = 3;
                                K(rr,1) = K(rr,1)-2*R_row2;

                                V3 = find(Table(2+1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj));  % node 3
                                V2 = find(Table(1+1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj-1));  % node 2
                                V1 = find(Table(1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj-2)) ;  % node 1

                                x13 = (phi(VR1+V1,1)-phi(VR1+V3,1))/Thermal_V;
                                x31 = -x13;
                                x23 = (phi(VR1+V2,1)-phi(VR1+V3,1))/Thermal_V;
                                x32 = -x23;

                                res(ii,1) = res(ii,1) + coeff_Jp*((Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj))*(hole(VR1+V1,1)*Ber(x31) - hole(VR1+V3,1)*Ber(x13)));
                                res(ii,1) = res(ii,1) + coeff_Jp*((Area(K(rr,1)+R_row1,jj-1)/L(K(rr,1)+R_row1,jj-1))*(hole(VR1+V2,1)*Ber(x32) - hole(VR1+V3,1)*Ber(x23)));

                                %hDensity
                                Jaco(ii, 3*(V1)+(VR1)) = Jaco(ii, 3*(V1)+(VR1)) + coeff_Jp*Ber(x31)*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj));
                                Jaco(ii, 3*(V2)+(VR1)) = Jaco(ii, 3*(V2)+(VR1)) + coeff_Jp*Ber(x32)*(Area(K(rr,1)+R_row1,jj-1)/L(K(rr,1)+R_row1,jj-1));

                                Jaco(ii, 3*(V3)+(VR1)) = Jaco(ii, 3*(V3)+(VR1)) + coeff_Jp*(-Ber(x13))*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj));
                                Jaco(ii, 3*(V3)+(VR1)) = Jaco(ii, 3*(V3)+(VR1)) + coeff_Jp*(-Ber(x23))*(Area(K(rr,1)+R_row1,jj-1)/L(K(rr,1)+R_row1,jj-1));

                                %Potential
                                Jaco(ii, 3*(V3)+(VR1-2)) = Jaco(ii, 3*(V3)+(VR1-2)) + coeff_Jp/Thermal_V*(Area(K(rr,1)+R_row1,jj-1)/L(K(rr,1)+R_row1,jj-1))*(hole(VR1+V2,1)*Ber_d(x32)+hole(VR1+V3,1)*Ber_d(x23));
                                Jaco(ii, 3*(V3)+(VR1-2)) = Jaco(ii, 3*(V3)+(VR1-2)) + coeff_Jp/Thermal_V*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj))*(hole(VR1+V1,1)*Ber_d(x31)+hole(VR1+V3,1)*Ber_d(x13));

                                Jaco(ii, 3*(V2)+(VR1-2)) = Jaco(ii, 3*(V2)+(VR1-2)) + coeff_Jp/Thermal_V*(Area(K(rr,1)+R_row1,jj-1)/L(K(rr,1)+R_row1,jj-1))*(-hole(VR1+V2,1)*Ber_d(x32)-hole(VR1+V3,1)*Ber_d(x23));
                                Jaco(ii, 3*(V1)+(VR1-2)) = Jaco(ii, 3*(V1)+(VR1-2)) + coeff_Jp/Thermal_V*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj))*(-hole(VR1+V1,1)*Ber_d(x31)-hole(VR1+V3,1)*Ber_d(x13));
                            end
                        end

                    end
                end

                if (ry1(kk,1) <= Table(ii,2) && Table(ii,2) <= iy1(kk,1)) || (ry2(kk,1) >= Table(ii,2) && Table(ii,2) >= iy2(kk,1))

                    K = find(E_R2 == Table(ii,2));
                    K_row = size(K,1);
                    if Table(ii,3) == Potential

                        for rr = 1:K_row

                            if K(rr,1) <= R_row2
                                jj = 1;

                                V3 = find(Table(1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj+2)) ;
                                V2 = find(Table(1+1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj+1));
                                V1 = find(Table(1+1+1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj));

                                res(ii,1) = res(ii,1) + esi*(phi(VR1+V2,1)-phi(VR1+V1,1))*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj)) + esi*(phi(VR1+V3,1)-phi(VR1+V1,1))*(Area(K(rr,1)+R_row1,jj+2)/L(K(rr,1)+R_row1,jj+2));
                                res(ii,1) = res(ii,1) + coeff*(Area(K(rr,1)+R_row1,jj)*L(K(rr,1)+R_row1,jj)+Area(K(rr,1)+R_row1,jj+2)*L(K(rr,1)+R_row1,jj+2))*(Nd-elec(VR1+V1,1)+hole(VR1+V1,1));

                                Jaco(ii, 3*(V3)+(VR1-2)) = Jaco(ii, 3*(V3)+(VR1-2))+ esi*(Area(K(rr,1)+R_row1,jj+2)/L(K(rr,1)+R_row1,jj+2));
                                Jaco(ii, 3*(V2)+(VR1-2)) = Jaco(ii, 3*(V2)+(VR1-2))+esi*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj));
                                Jaco(ii, 3*(V1)+(VR1-2)) = Jaco(ii, 3*(V1)+(VR1-2))-esi*(Area(K(rr,1)+R_row1,jj+2)/L(K(rr,1)+R_row1,jj+2))-esi*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj));

                                Jaco(ii, 3*(V1)+(VR1-1)) = Jaco(ii, 3*(V1)+(VR1-1)) - coeff*(Area(K(rr,1)+R_row1,jj)*L(K(rr,1)+R_row1,jj)+Area(K(rr,1)+R_row1,jj+2)*L(K(rr,1)+R_row1,jj+2));
                                Jaco(ii, 3*(V1)+(VR1)) = Jaco(ii, 3*(V1)+(VR1)) + coeff*(Area(K(rr,1)+R_row1,jj)*L(K(rr,1)+R_row1,jj)+Area(K(rr,1)+R_row1,jj+2)*L(K(rr,1)+R_row1,jj+2));

                            elseif K(rr,1) > R_row2 && K(rr,1) <= 2*R_row2
                                jj = 2;
                                K(rr,1) = K(rr,1)-R_row2;

                                V3 = find(Table(1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj+1)) ;
                                V2 = find(Table(1+1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj));
                                V1 = find(Table(1+1+1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj-1));

                                res(ii,1) = res(ii,1) + esi*(phi(VR1+V1,1)-phi(VR1+V2,1))*(Area(K(rr,1)+R_row1,jj-1)/L(K(rr,1)+R_row1,jj-1)) + esi*(phi(VR1+V3,1)-phi(VR1+V2,1))*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj));
                                res(ii,1) = res(ii,1) + coeff*(Area(K(rr,1)+R_row1,jj)*L(K(rr,1)+R_row1,jj)+Area(K(rr,1)+R_row1,jj-1)*L(K(rr,1)+R_row1,jj-1))*(Nd-elec(VR1+V2,1)+hole(VR1+V2,1));

                                Jaco(ii, 3*(V3)+(VR1-2)) = Jaco(ii, 3*(V3)+(VR1-2))+esi*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj));
                                Jaco(ii, 3*(V2)+(VR1-2)) =  Jaco(ii, 3*(V2)+(VR1-2))-esi*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj))-esi*(Area(K(rr,1)+R_row1,jj-1)/L(K(rr,1)+R_row1,jj-1));
                                Jaco(ii, 3*(V1)+(VR1-2)) = Jaco(ii, 3*(V1)+(VR1-2))+esi*(Area(K(rr,1)+R_row1,jj-1)/L(K(rr,1)+R_row1,jj-1));

                                Jaco(ii, 3*(V2)+(VR1-1)) = Jaco(ii, 3*(V2)+(VR1-1)) - coeff*(Area(K(rr,1)+R_row1,jj)*L(K(rr,1)+R_row1,jj)+Area(K(rr,1)+R_row1,jj-1)*L(K(rr,1)+R_row1,jj-1));
                                Jaco(ii, 3*(V2)+(VR1)) = Jaco(ii, 3*(V2)+(VR1)) + coeff*(Area(K(rr,1)+R_row1,jj)*L(K(rr,1)+R_row1,jj)+Area(K(rr,1)+R_row1,jj-1)*L(K(rr,1)+R_row1,jj-1));

                            elseif K(rr,1) > 2*R_row2 && K(rr,1) <= 3*R_row2
                                jj = 3;
                                K(rr,1) = K(rr,1)-2*R_row2;

                                V3 = find(Table(1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj-2)) ;
                                V2 = find(Table(1+1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj-1));
                                V1 = find(Table(1+1+1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj));

                                res(ii,1) = res(ii,1) + esi*(phi(VR1+V3,1)-phi(VR1+V1,1))*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj)) + esi*(phi(VR1+V2,1)-phi(VR1+V1,1))*(Area(K(rr,1)+R_row1,jj-1)/L(K(rr,1)+R_row1,jj-1));
                                res(ii,1) = res(ii,1) + coeff*(Area(K(rr,1)+R_row1,jj)*L(K(rr,1)+R_row1,jj)+Area(K(rr,1)+R_row1,jj-1)*L(K(rr,1)+R_row1,jj-1))*(Nd-elec(VR1+V1,1)+hole(VR1+V1,1));

                                Jaco(ii, 3*(V3)+(VR1-2)) = Jaco(ii, 3*(V3)+(VR1-2))+esi*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj));
                                Jaco(ii, 3*(V2)+(VR1-2)) = Jaco(ii, 3*(V2)+(VR1-2))+esi*(Area(K(rr,1)+R_row1,jj-1)/L(K(rr,1)+R_row1,jj-1));
                                Jaco(ii, 3*(V1)+(VR1-2)) = Jaco(ii, 3*(V1)+(VR1-2))-esi*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj))-esi*(Area(K(rr,1)+R_row1,jj-1)/L(K(rr,1)+R_row1,jj-1));

                                Jaco(ii, 3*(V1)+(VR1-1)) = Jaco(ii, 3*(V1)+(VR1-1)) - coeff*(Area(K(rr,1)+R_row1,jj)*L(K(rr,1)+R_row1,jj)+Area(K(rr,1)+R_row1,jj-1)*L(K(rr,1)+R_row1,jj-1));
                                Jaco(ii, 3*(V1)+(VR1)) = Jaco(ii, 3*(V1)+(VR1)) + coeff*(Area(K(rr,1)+R_row1,jj)*L(K(rr,1)+R_row1,jj)+Area(K(rr,1)+R_row1,jj-1)*L(K(rr,1)+R_row1,jj-1));
                            end
                        end

                    elseif Table(ii,3) == eDensity

                        for rr = 1:K_row
                            if K(rr,1) <= R_row2
                                jj = 1;

                                V3 = find(Table(1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj+2)) ; % node 3
                                V2 = find(Table(1+1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj+1)); % node 2
                                V1 = find(Table(2+1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj)); %node 1

                                x21 = (phi(VR1+V2,1)-phi(VR1+V1,1))/Thermal_V;
                                x12 = -x21;
                                x31 = (phi(VR1+V3,1)-phi(VR1+V1,1))/Thermal_V;
                                x13 = -x31;

                                res(ii,1) = res(ii,1) + coeff_Jn*(elec(VR1+V2,1)*Ber(x21) - elec(VR1+V1,1)*Ber(x12))*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj));
                                res(ii,1) = res(ii,1) + coeff_Jn*(elec(VR1+V3,1)*Ber(x31) - elec(VR1+V1,1)*Ber(x13))*(Area(K(rr,1)+R_row1,jj+2)/L(K(rr,1)+R_row1,jj+2));

                                % eDensity
                                Jaco(ii, 3*(V1)+(VR1-1)) = Jaco(ii, 3*(V1)+(VR1-1)) + coeff_Jn*-Ber(x12)*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj));
                                Jaco(ii, 3*(V1)+(VR1-1)) = Jaco(ii, 3*(V1)+(VR1-1)) + coeff_Jn*(-Ber(x13))*(Area(K(rr,1)+R_row1,jj+2)/L(K(rr,1)+R_row1,jj+2));

                                Jaco(ii, 3*(V2)+(VR1-1)) = Jaco(ii, 3*(V2)+(VR1-1)) + coeff_Jn*Ber(x21)*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj));
                                Jaco(ii, 3*(V3)+(VR1-1)) = Jaco(ii, 3*(V3)+(VR1-1)) + coeff_Jn*Ber(x31)*(Area(K(rr,1)+R_row1,jj+2)/L(K(rr,1)+R_row1,jj+2));

                                % Potential
                                Jaco(ii, 3*(V1)+(VR1-2)) = Jaco(ii, 3*(V1)+(VR1-2)) + coeff_Jn/Thermal_V*(-elec(VR1+V2,1)*Ber_d(x21)-elec(VR1+V1,1)*Ber_d(x12))*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj));
                                Jaco(ii, 3*(V1)+(VR1-2)) = Jaco(ii, 3*(V1)+(VR1-2)) + coeff_Jn/Thermal_V*(-elec(VR1+V3,1)*Ber_d(x31)-elec(VR1+V1,1)*Ber_d(x13))*(Area(K(rr,1)+R_row1,jj+2)/L(K(rr,1)+R_row1,jj+2));

                                Jaco(ii, 3*(V2)+(VR1-2)) = Jaco(ii, 3*(V2)+(VR1-2)) + coeff_Jn/Thermal_V*(elec(VR1+V2,1)*Ber_d(x21)+elec(VR1+V1,1)*Ber_d(x12))*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj));
                                Jaco(ii, 3*(V3)+(VR1-2)) = Jaco(ii, 3*(V3)+(VR1-2)) + coeff_Jn/Thermal_V*(elec(VR1+V3,1)*Ber_d(x31)+elec(VR1+V1,1)*Ber_d(x13))*(Area(K(rr,1)+R_row1,jj+2)/L(K(rr,1)+R_row1,jj+2));

                            elseif K(rr,1) > R_row2 && K(rr,1) <= 2*R_row2
                                jj = 2;
                                K(rr,1) = K(rr,1)-R_row2;

                                V3 = find(Table(1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj+1)) ; % node 3
                                V2 = find(Table(1+1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj)); % node 2
                                V1 = find(Table(1+1+1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj-1)); % node 1

                                x12 = (phi(VR1+V1,1)-phi(VR1+V2,1))/Thermal_V;
                                x21 = -x12;
                                x32 = (phi(VR1+V3,1)-phi(VR1+V2,1))/Thermal_V;
                                x23 = -x32;

                                res(ii,1) =  res(ii,1) + coeff_Jn*((Area(K(rr,1)+R_row1,jj-1)/L(K(rr,1)+R_row1,jj-1))*(elec(VR1+V1,1)*Ber(x12) - elec(VR1+V2,1)*Ber(x21)));
                                res(ii,1) =  res(ii,1) + coeff_Jn*((Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj))*(elec(VR1+V3,1)*Ber(x32) - elec(VR1+V2,1)*Ber(x23)));

                                % eDensity
                                Jaco(ii, 3*(V2)+(VR1-1)) = Jaco(ii, 3*(V2)+(VR1-1)) + coeff_Jn*(-Ber(x21))*(Area(K(rr,1)+R_row1,jj-1)/L(K(rr,1)+R_row1,jj-1));
                                Jaco(ii, 3*(V2)+(VR1-1)) = Jaco(ii, 3*(V2)+(VR1-1)) + coeff_Jn*(-Ber(x23))*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj));

                                Jaco(ii, 3*(V1)+(VR1-1)) = Jaco(ii, 3*(V1)+(VR1-1)) + coeff_Jn*Ber(x12)*(Area(K(rr,1)+R_row1,jj-1)/L(K(rr,1)+R_row1,jj-1));
                                Jaco(ii, 3*(V3)+(VR1-1)) = Jaco(ii, 3*(V3)+(VR1-1)) + coeff_Jn*Ber(x32)*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj));

                                % Potential
                                Jaco(ii, 3*(V2)+(VR1-2)) = Jaco(ii, 3*(V2)+(VR1-2)) + coeff_Jn/Thermal_V*(-elec(VR1+V1,1)*Ber_d(x12)-elec(VR1+V2,1)*Ber_d(x21))*(Area(K(rr,1)+R_row1,jj-1)/L(K(rr,1)+R_row1,jj-1));
                                Jaco(ii, 3*(V2)+(VR1-2)) = Jaco(ii, 3*(V2)+(VR1-2)) + coeff_Jn/Thermal_V*(-elec(VR1+V3,1)*Ber_d(x32)-elec(VR1+V2,1)*Ber_d(x23))*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj));

                                Jaco(ii, 3*(V1)+(VR1-2)) = Jaco(ii, 3*(V1)+(VR1-2)) + coeff_Jn/Thermal_V*(elec(VR1+V1,1)*Ber_d(x12)+elec(VR1+V2,1)*Ber_d(x21))*(Area(K(rr,1)+R_row1,jj-1)/L(K(rr,1)+R_row1,jj-1));
                                Jaco(ii, 3*(V3)+(VR1-2)) = Jaco(ii, 3*(V3)+(VR1-2)) + coeff_Jn/Thermal_V*(elec(VR1+V3,1)*Ber_d(x32)+elec(VR1+V2,1)*Ber_d(x23))*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj));

                            elseif K(rr,1) > 2*R_row2 && K(rr,1) <= 3*R_row2
                                jj = 3;
                                K(rr,1) = K(rr,1)-2*R_row2;

                                V1 = find(Table(1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj-2)) ; % node 1
                                V2 = find(Table(1+1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj-1)); % node 2
                                V3 = find(Table(1+1+1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj)); % node 3

                                x13 = (phi(VR1+V1,1)-phi(VR1+V3,1))/Thermal_V;
                                x31 = -x13;
                                x23 = (phi(VR1+V2,1)-phi(VR1+V3,1))/Thermal_V;
                                x32 = -x23;

                                res(ii,1) = res(ii,1) + coeff_Jn*((Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj))*(elec(VR1+V1,1)*Ber(x13) - elec(VR1+V3,1)*Ber(x31)));
                                res(ii,1) = res(ii,1) + coeff_Jn*((Area(K(rr,1)+R_row1,jj-1)/L(K(rr,1)+R_row1,jj-1))*(elec(VR1+V2,1)*Ber(x23) - elec(VR1+V3,1)*Ber(x32)));

                                % eDensity
                                Jaco(ii, 3*(V1)+(VR1-1)) = Jaco(ii, 3*(V1)+(VR1-1)) + coeff_Jn*Ber(x13)*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj));
                                Jaco(ii, 3*(V2)+(VR1-1)) = Jaco(ii, 3*(V2)+(VR1-1)) + coeff_Jn*Ber(x23)*(Area(K(rr,1)+R_row1,jj-1)/L(K(rr,1)+R_row1,jj-1));

                                Jaco(ii, 3*(V3)+(VR1-1)) = Jaco(ii, 3*(V3)+(VR1-1)) + coeff_Jn*(-Ber(x31))*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj));
                                Jaco(ii, 3*(V3)+(VR1-1)) = Jaco(ii, 3*(V3)+(VR1-1)) + coeff_Jn*(-Ber(x32))*(Area(K(rr,1)+R_row1,jj-1)/L(K(rr,1)+R_row1,jj-1));

                                % Potential
                                Jaco(ii, 3*(V3)+(VR1-2)) = Jaco(ii, 3*(V3)+(VR1-2)) + coeff_Jn/Thermal_V*(-elec(VR1+V1,1)*Ber_d(x13)-elec(VR1+V3,1)*Ber_d(x31))*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj));
                                Jaco(ii, 3*(V3)+(VR1-2)) = Jaco(ii, 3*(V3)+(VR1-2)) + coeff_Jn/Thermal_V*(-elec(VR1+V2,1)*Ber_d(x23)-elec(VR1+V3,1)*Ber_d(x32))*(Area(K(rr,1)+R_row1,jj-1)/L(K(rr,1)+R_row1,jj-1));

                                Jaco(ii, 3*(V2)+(VR1-2)) = Jaco(ii, 3*(V2)+(VR1-2)) + coeff_Jn/Thermal_V*(elec(VR1+V2,1)*Ber_d(x23)+elec(VR1+V3,1)*Ber_d(x32))*(Area(K(rr,1)+R_row1,jj-1)/L(K(rr,1)+R_row1,jj-1));
                                Jaco(ii, 3*(V1)+(VR1-2)) = Jaco(ii, 3*(V1)+(VR1-2)) + coeff_Jn/Thermal_V*(elec(VR1+V1,1)*Ber_d(x13)+elec(VR1+V3,1)*Ber_d(x31))*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj));

                            end
                        end


                    elseif Table(ii,3) == hDensity

                        for rr = 1:K_row

                            if K(rr,1) <= R_row2
                                jj = 1;

                                V3 = find(Table(1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj+2)) ;  % node 3
                                V2 = find(Table(1+1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj+1));  % node 2
                                V1 = find(Table(2+1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj));  % node 1

                                x21 = (phi(VR1+V2,1)-phi(VR1+V1,1))/Thermal_V;
                                x12 = -x21;
                                x31 = (phi(VR1+V3,1)-phi(VR1+V1,1))/Thermal_V;
                                x13 = -x31;

                                res(ii,1) = res(ii,1) + coeff_Jp*(hole(VR1+V2,1)*Ber(x12) - hole(VR1+V1,1)*Ber(x21))*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj));
                                res(ii,1) = res(ii,1) + coeff_Jp*(hole(VR1+V3,1)*Ber(x13) - hole(VR1+V1,1)*Ber(x31))*(Area(K(rr,1)+R_row1,jj+2)/L(K(rr,1)+R_row1,jj+2));

                                %hDensity
                                Jaco(ii, 3*(V1)+(VR1)) = Jaco(ii, 3*(V1)+(VR1)) + coeff_Jp*(-Ber(x21))*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj));
                                Jaco(ii, 3*(V1)+(VR1)) = Jaco(ii, 3*(V1)+(VR1)) + coeff_Jp*(-Ber(x31))*(Area(K(rr,1)+R_row1,jj+2)/L(K(rr,1)+R_row1,jj+2));

                                Jaco(ii, 3*(V2)+(VR1)) = Jaco(ii, 3*(V2)+(VR1)) + coeff_Jp*Ber(x12)*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj));
                                Jaco(ii, 3*(V3)+(VR1)) = Jaco(ii, 3*(V3)+(VR1)) + coeff_Jp*Ber(x13)*(Area(K(rr,1)+R_row1,jj+2)/L(K(rr,1)+R_row1,jj+2));

                                %Potential
                                Jaco(ii, 3*(V1)+(VR1-2)) = Jaco(ii, 3*(V1)+(VR1-2)) + coeff_Jp/Thermal_V*(hole(VR1+V2,1)*Ber_d(x12)+hole(VR1+V1,1)*Ber_d(x21))*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj));
                                Jaco(ii, 3*(V1)+(VR1-2)) = Jaco(ii, 3*(V1)+(VR1-2)) + coeff_Jp/Thermal_V*(hole(VR1+V3,1)*Ber_d(x13)+hole(VR1+V1,1)*Ber_d(x31))*(Area(K(rr,1)+R_row1,jj+2)/L(K(rr,1)+R_row1,jj+2));

                                Jaco(ii, 3*(V2)+(VR1-2)) = Jaco(ii, 3*(V2)+(VR1-2)) + coeff_Jp/Thermal_V*(-hole(VR1+V2,1)*Ber_d(x12)-hole(VR1+V1,1)*Ber_d(x21))*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj));
                                Jaco(ii, 3*(V3)+(VR1-2)) = Jaco(ii, 3*(V3)+(VR1-2)) + coeff_Jp/Thermal_V*(-hole(VR1+V3,1)*Ber_d(x13)-hole(VR1+V1,1)*Ber_d(x31))*(Area(K(rr,1)+R_row1,jj+2)/L(K(rr,1)+R_row1,jj+2));


                            elseif K(rr,1) > R_row2 && K(rr,1) <= 2*R_row2
                                jj = 2;
                                K(rr,1) = K(rr,1)-R_row2;

                                V3 = find(Table(1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj+1)) ;  % node 3
                                V2 = find(Table(1+1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj));  % node 2
                                V1 = find(Table(2+1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj-1));  % node 1

                                x12 = (phi(VR1+V1,1)-phi(VR1+V2,1))/Thermal_V;
                                x21 = -x12;
                                x32 = (phi(VR1+V3,1)-phi(VR1+V2,1))/Thermal_V;
                                x23 = -x32;

                                res(ii,1) =  res(ii,1) + coeff_Jp*((Area(K(rr,1)+R_row1,jj-1)/L(K(rr,1)+R_row1,jj-1))*(hole(VR1+V1,1)*Ber(x21) - hole(VR1+V2,1)*Ber(x12)));
                                res(ii,1) =  res(ii,1) + coeff_Jp*((Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj))*(hole(VR1+V3,1)*Ber(x23) - hole(VR1+V2,1)*Ber(x32)));

                                % hDensity
                                Jaco(ii, 3*(V1)+(VR1)) = Jaco(ii, 3*(V1)+(VR1)) + coeff_Jp*Ber(x21)*(Area(K(rr,1)+R_row1,jj-1)/L(K(rr,1)+R_row1,jj-1));
                                Jaco(ii, 3*(V3)+(VR1)) = Jaco(ii, 3*(V3)+(VR1)) + coeff_Jp*Ber(x23)*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj));

                                Jaco(ii, 3*(V2)+(VR1)) = Jaco(ii, 3*(V2)+(VR1)) + coeff_Jp*(-Ber(x12))*(Area(K(rr,1)+R_row1,jj-1)/L(K(rr,1)+R_row1,jj-1));
                                Jaco(ii, 3*(V2)+(VR1)) = Jaco(ii, 3*(V2)+(VR1)) + coeff_Jp*(-Ber(x32))*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj));


                                %Potential
                                Jaco(ii, 3*(V2)+(VR1-2)) = Jaco(ii, 3*(V2)+(VR1-2)) + coeff_Jp/Thermal_V*(hole(VR1+V1,1)*Ber_d(x21)+hole(VR1+V2,1)*Ber_d(x12))*(Area(K(rr,1)+R_row1,jj-1)/L(K(rr,1)+R_row1,jj-1));
                                Jaco(ii, 3*(V2)+(VR1-2)) = Jaco(ii, 3*(V2)+(VR1-2)) + coeff_Jp/Thermal_V*(hole(VR1+V3,1)*Ber_d(x23)+hole(VR1+V2,1)*Ber_d(x32))*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj));

                                Jaco(ii, 3*(V1)+(VR1-2)) = Jaco(ii, 3*(V1)+(VR1-2)) + coeff_Jp/Thermal_V*(-hole(VR1+V1,1)*Ber_d(x21)-hole(VR1+V2,1)*Ber_d(x12))*(Area(K(rr,1)+R_row1,jj-1)/L(K(rr,1)+R_row1,jj-1));
                                Jaco(ii, 3*(V3)+(VR1-2)) = Jaco(ii, 3*(V3)+(VR1-2)) + coeff_Jp/Thermal_V*(-hole(VR1+V3,1)*Ber_d(x23)-hole(VR1+V2,1)*Ber_d(x32))*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj));


                            elseif K(rr,1) > 2*R_row2 && K(rr,1) <= 3*R_row2
                                jj = 3;
                                K(rr,1) = K(rr,1)-2*R_row2;

                                V3 = find(Table(2+1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj));  % node 3
                                V2 = find(Table(1+1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj-1));  % node 2
                                V1 = find(Table(1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj-2)) ;  % node 1

                                x13 = (phi(VR1+V1,1)-phi(VR1+V3,1))/Thermal_V;
                                x31 = -x13;
                                x23 = (phi(VR1+V2,1)-phi(VR1+V3,1))/Thermal_V;
                                x32 = -x23;

                                res(ii,1) = res(ii,1) + coeff_Jp*((Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj))*(hole(VR1+V1,1)*Ber(x31) - hole(VR1+V3,1)*Ber(x13)));
                                res(ii,1) = res(ii,1) + coeff_Jp*((Area(K(rr,1)+R_row1,jj-1)/L(K(rr,1)+R_row1,jj-1))*(hole(VR1+V2,1)*Ber(x32) - hole(VR1+V3,1)*Ber(x23)));

                                %hDensity
                                Jaco(ii, 3*(V1)+(VR1)) = Jaco(ii, 3*(V1)+(VR1)) + coeff_Jp*Ber(x31)*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj));
                                Jaco(ii, 3*(V2)+(VR1)) = Jaco(ii, 3*(V2)+(VR1)) + coeff_Jp*Ber(x32)*(Area(K(rr,1)+R_row1,jj-1)/L(K(rr,1)+R_row1,jj-1));

                                Jaco(ii, 3*(V3)+(VR1)) = Jaco(ii, 3*(V3)+(VR1)) + coeff_Jp*(-Ber(x13))*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj));
                                Jaco(ii, 3*(V3)+(VR1)) = Jaco(ii, 3*(V3)+(VR1)) + coeff_Jp*(-Ber(x23))*(Area(K(rr,1)+R_row1,jj-1)/L(K(rr,1)+R_row1,jj-1));

                                %Potential
                                Jaco(ii, 3*(V3)+(VR1-2)) = Jaco(ii, 3*(V3)+(VR1-2)) + coeff_Jp/Thermal_V*(Area(K(rr,1)+R_row1,jj-1)/L(K(rr,1)+R_row1,jj-1))*(hole(VR1+V2,1)*Ber_d(x32)+hole(VR1+V3,1)*Ber_d(x23));
                                Jaco(ii, 3*(V3)+(VR1-2)) = Jaco(ii, 3*(V3)+(VR1-2)) + coeff_Jp/Thermal_V*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj))*(hole(VR1+V1,1)*Ber_d(x31)+hole(VR1+V3,1)*Ber_d(x13));

                                Jaco(ii, 3*(V2)+(VR1-2)) = Jaco(ii, 3*(V2)+(VR1-2)) + coeff_Jp/Thermal_V*(Area(K(rr,1)+R_row1,jj-1)/L(K(rr,1)+R_row1,jj-1))*(-hole(VR1+V2,1)*Ber_d(x32)-hole(VR1+V3,1)*Ber_d(x23));
                                Jaco(ii, 3*(V1)+(VR1-2)) = Jaco(ii, 3*(V1)+(VR1-2)) + coeff_Jp/Thermal_V*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj))*(-hole(VR1+V1,1)*Ber_d(x31)-hole(VR1+V3,1)*Ber_d(x13));
                            end
                        end

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

                        res(ii,1) = res(ii,1) + eox*(phi(VR1+VR2+V2,1)-phi(VR1+VR2+V1,1))*(Area(K(rr,1)+R_row2+R_row1,jj)/L(K(rr,1)+R_row2+R_row1,jj)) + eox*(phi(VR1+VR2+V3,1)-phi(VR1+VR2+V1,1))*(Area(K(rr,1)+R_row2+R_row1,jj+2)/L(K(rr,1)+R_row2+R_row1,jj+2));

                        Jaco(ii, VR1+3*VR2+V3) = Jaco(ii, VR1+3*VR2+V3)+ eox*(Area(K(rr,1)+R_row2+R_row1,jj+2)/L(K(rr,1)+R_row2+R_row1,jj+2));
                        Jaco(ii, VR1+3*VR2+V2) = Jaco(ii, VR1+3*VR2+V2)+eox*(Area(K(rr,1)+R_row2+R_row1,jj)/L(K(rr,1)+R_row2+R_row1,jj));
                        Jaco(ii, VR1+3*VR2+V1) = Jaco(ii, VR1+3*VR2+V1)-eox*(Area(K(rr,1)+R_row2+R_row1,jj+2)/L(K(rr,1)+R_row2+R_row1,jj+2))-eox*(Area(K(rr,1)+R_row2+R_row1,jj)/L(K(rr,1)+R_row2+R_row1,jj));

                    elseif K(rr,1) > R_row3 && K(rr,1) <= 2*R_row3

                        jj = 2;
                        K(rr,1) = K(rr,1)-R_row3;

                        V3 = find(Table(1+VR1+3*VR2:T1,2)==E_R3(K(rr,1),jj+1)) ;
                        V2 = find(Table(1+VR1+3*VR2:T1,2)==E_R3(K(rr,1),jj));
                        V1 = find(Table(1+VR1+3*VR2:T1,2)==E_R3(K(rr,1),jj-1));

                        res(ii,1) = res(ii,1) + eox*(phi(VR1+VR2+V1,1)-phi(VR1+VR2+V2,1))*(Area(K(rr,1)+R_row2+R_row1,jj-1)/L(K(rr,1)+R_row2+R_row1,jj-1)) + eox*(phi(VR1+VR2+V3,1)-phi(VR1+VR2+V2,1))*(Area(K(rr,1)+R_row2+R_row1,jj)/L(K(rr,1)+R_row2+R_row1,jj));

                        Jaco(ii, VR1+3*VR2+V3) = Jaco(ii, VR1+3*VR2+V3)+eox*(Area(K(rr,1)+R_row2+R_row1,jj)/L(K(rr,1)+R_row2+R_row1,jj));
                        Jaco(ii, VR1+3*VR2+V2) =  Jaco(ii, VR1+3*VR2+V2)-eox*(Area(K(rr,1)+R_row2+R_row1,jj)/L(K(rr,1)+R_row2+R_row1,jj))-eox*(Area(K(rr,1)+R_row2+R_row1,jj-1)/L(K(rr,1)+R_row2+R_row1,jj-1));
                        Jaco(ii, VR1+3*VR2+V1) = Jaco(ii, VR1+3*VR2+V1)+eox*(Area(K(rr,1)+R_row2+R_row1,jj-1)/L(K(rr,1)+R_row2+R_row1,jj-1));


                    elseif K(rr,1) > 2*R_row3 && K(rr,1) <= 3*R_row3

                        jj = 3;
                        K(rr,1) = K(rr,1)-2*R_row3;

                        V3 = find(Table(1+VR1+3*VR2:T1,2)==E_R3(K(rr,1),jj-2)) ;
                        V2 = find(Table(1+VR1+3*VR2:T1,2)==E_R3(K(rr,1),jj-1));
                        V1 = find(Table(1+VR1+3*VR2:T1,2)==E_R3(K(rr,1),jj));

                        res(ii,1) = res(ii,1) + eox*(phi(VR1+VR2+V3,1)-phi(VR1+VR2+V1,1))*(Area(K(rr,1)+R_row2+R_row1,jj)/L(K(rr,1)+R_row2+R_row1,jj)) + eox*(phi(VR1+VR2+V2,1)-phi(VR1+VR2+V1,1))*(Area(K(rr,1)+R_row2+R_row1,jj-1)/L(K(rr,1)+R_row2+R_row1,jj-1));

                        Jaco(ii, VR1+3*VR2+V3) = Jaco(ii, VR1+3*VR2+V3) + eox*(Area(K(rr,1)+R_row2+R_row1,jj)/L(K(rr,1)+R_row2+R_row1,jj));
                        Jaco(ii, VR1+3*VR2+V2) = Jaco(ii, VR1+3*VR2+V2) + eox*(Area(K(rr,1)+R_row2+R_row1,jj-1)/L(K(rr,1)+R_row2+R_row1,jj-1));
                        Jaco(ii, VR1+3*VR2+V1) = Jaco(ii, VR1+3*VR2+V1) - eox*(Area(K(rr,1)+R_row2+R_row1,jj)/L(K(rr,1)+R_row2+R_row1,jj))-eox*(Area(K(rr,1)+R_row2+R_row1,jj-1)/L(K(rr,1)+R_row2+R_row1,jj-1));

                    end
                end
            end

        end
    end

    %    Boundary Condition

    for rr = 1:i_row1

        Vi11 = find(Table(1:VR1,2)==interface1(rr,1));
        Vi21 = find(Table(1+VR1:3:VR1+3*VR2,2)==interface1(rr,1));

        Jaco(3*Vi21+(VR1-2),:)=Jaco(3*Vi21+(VR1-2),:)+Jaco(Vi11,:);
        Jaco(Vi11,:)=0; Jaco(Vi11,Vi11)=1; Jaco(Vi11,3*Vi21+(VR1-2))=-1;

        res(3*Vi21+(VR1-2),1)=res(3*Vi21+(VR1-2),1)+res(Vi11,1);
        res(Vi11,1)=0;
    end

    for rr = 1:i_row2

        Vi12 = find(Table(1+VR1+3*VR2:T1,2)==interface2(rr,1));
        Vi22 = find(Table(1+VR1:3:VR1+3*VR2,2)==interface2(rr,1));

        Jaco(3*Vi22+(VR1-2),:)=Jaco(3*Vi22+(VR1-2),:)+Jaco(VR1+3*VR2+Vi12,:);
        Jaco(VR1+3*VR2+Vi12,:)=0; Jaco(VR1+3*VR2+Vi12,VR1+3*VR2+Vi12)=1; Jaco(VR1+3*VR2+Vi12,3*Vi22+(VR1-2))=-1;

        res(3*Vi22+(VR1-2),1)=res(3*Vi22+(VR1-2),1)+res(VR1+3*VR2+Vi12,1);
        res(VR1+3*VR2+Vi12,1)=0;

    end

    %%%% Dirichlet Boundary Condition %%%%%

    for n=1:C_col
        for i =1: C_row
            if i <= 20 % Bottom Gate
                V = find(Table(1:VR1,2)==Contact(i,n));

                Jaco(V,:) =0;
                Jaco(V,V) = 1;
                Jaco(V,T1+7) = -1;
                res(V,1)=phi(V,1)-0.33374;

            else % Top Gate

                V = find(Table(1+VR1+3*VR2:T_row,2)==Contact(i,n));

                Jaco(VR1+3*VR2+V,:) =0;
                Jaco(VR1+3*VR2+V,VR1+3*VR2+V) = 1;
                Jaco(VR1+3*VR2+V,T1+8) = -1;
                res(VR1+3*VR2+V,1)= phi(VR1+VR2+V,1)-0.33374;
            end
        end
    end

    %%%% Source, Drain %%%%
    q_row = size(ry1,1);
    w_row = size(ry2,1);

    for rr = 1:q_row  % Source

        V = find(Table(1+VR1:3:VR1+3*VR2,2)==ry1(rr,1));

        Jaco(3*V+(VR1-2),:) = 0;
        Jaco(3*V+(VR1-2),3*V+(VR1-2)) = 1;
        Jaco(3*V+(VR1-2),T1+5) = -1;
        Jaco(3*V+(VR1-1),:) = 0;
        Jaco(3*V+(VR1-1),3*V+(VR1-1)) = 1;
        Jaco(3*V+(VR1),:) = 0;
        Jaco(3*V+(VR1),3*V+(VR1)) = 1;

        res(3*V+(VR1-2),1) = phi(VR1+V,1)-Thermal_V*log(Nd/nint);
        res(3*V+(VR1-1),1) = elec(VR1+V,1)-Nd;
        res(3*V+(VR1),1) = hole(VR1+V,1)-(nint^2/Nd);

    end

    for rr = 1:w_row  % Drain

        V = find(Table(1+VR1:3:VR1+3*VR2,2)==ry2(rr,1));

        Jaco(T1+2,:)=Jaco(T1+2,:)+Jaco(3*V+(VR1-1),:)+Jaco(3*V+(VR1),:);
        Jaco(T1+2,:)=-1*Jaco(T1+2,:);

        res(T1+2,1)= res(T1+2,1)+res(3*V+(VR1-1),1)+res(3*V+(VR1),1);
        res(T1+2,1)=-1*res(T1+2,1);

        Jaco(3*V+(VR1-2),:) = 0;
        Jaco(3*V+(VR1-2),3*V+(VR1-2)) = 1;
        Jaco(3*V+(VR1-2),T1+6) = -1;
        Jaco(3*V+(VR1-1),:) = 0;
        Jaco(3*V+(VR1-1),3*V+(VR1-1)) = 1;
        Jaco(3*V+(VR1),:) = 0;
        Jaco(3*V+(VR1),3*V+(VR1)) = 1;

        res(3*V+(VR1-2),1) = phi(VR1+V,1)-Thermal_V*log(Nd/nint);
        res(3*V+(VR1-1),1) = elec(VR1+V,1)-Nd;
        res(3*V+(VR1),1) = hole(VR1+V,1)-(nint^2/Nd);

    end

    %%%%% Multi-Mode Parameter %%%%%%% 
    R=1000; %%% ohm
    Vout=0;
    Vdd=1;

    % Is
    Jaco(T1+1,T1+1) = 1; Jaco(T1+1,T1+2) = 1;
    res(T1+1,1)=0; % equation
    %Id
    Jaco(T1+2,T1+2) = 1;
    % Ig1
    Jaco(T1+3,T1+3) = 1; res(T1+3,1)=0; % equation Ig1-0
    % Ig2
    Jaco(T1+4,T1+4) = 1; res(T1+4,1)=0; % equation Ig2-0
    % Vs
    Jaco(T1+5,T1+5) = 1; res(T1+5,1)=0; 
    % Vd
    Jaco(T1+6,T1+6) = 1; Jaco(T1+6,T1+13) = -1; res(T1+6,1)= 0;
    % Vg1
    Jaco(T1+7,T1+7) = 1; res(T1+7,1)=0;
    % Vg2
    Jaco(T1+8,T1+8) = 1; res(T1+8,1)=0;
    % I1
    Jaco(T1+9,T1+9) = 1; res(T1+9,1)=0;
    % I2
    Jaco(T1+10,T1+10) = 1; Jaco(T1+10,T1+11) = 1/R ; Jaco(T1+10,T1+12) = -1/R; res(T1+10,1)=0; %%
    % V1
    Jaco(T1+11,T1+11) = 1; res(T1+12,1)=0-Vout;
    % V2
    Jaco(T1+12,T1+12) = 1; res(T1+12,1)=0-1;
    % Vout
    Jaco(T1+13,T1+9) = 1; Jaco(T1+13,T1+2) = 1; Jaco(T1+13,T1+13) = 0; res(T1+parameter,1)=0;

    %     scaling
    Cvector = zeros(T_row,1);
    Cvector(1:VR1,1) = Thermal_V;
    Cvector(VR1+3*VR2:T_row,1) = Thermal_V;
    Cvector(1+VR1:3:VR1+3*VR2,1) = Thermal_V;
    Cvector(2+VR1:3:VR1+3*VR2,1) = Nd;
    Cvector(3+VR1:3:VR1+3*VR2,1) = Nd;
    Cvector(1+T1:T1+parameter,1) = Thermal_V;

    Cmatrix = spdiags(Cvector,0,T_row,T_row-parameter);
    Jaco_scaled = Jaco * Cmatrix;
    Rvector = 1./sum(abs(Jaco_scaled),2);
    Rmatrix = spdiags(Rvector,0,T_row,T_row);
    Jaco_scaled = Rmatrix* Jaco_scaled;
    res_scaled = Rmatrix *res;
    update_scaled = Jaco_scaled \ (-res_scaled);
    update_DD_vector_Gate(:,Newton) = Cmatrix* update_scaled;

    phi(1:VR1,1) = phi(1:VR1,1) + update_DD_vector_Gate(1:VR1,Newton);
    phi(1+VR1:VR2+VR1,1) = phi(1+VR1:VR2+VR1,1) + update_DD_vector_Gate(1+VR1:3:VR1+3*VR2,Newton);
    phi(1+VR1+VR2:T,1) = phi(1+VR1+VR2:T,1) + update_DD_vector_Gate(1+VR1+3*VR2:T_row-parameter,Newton);

    elec(1+VR1:1:VR2+VR1,1) = elec(1+VR1:1:VR2+VR1,1) + update_DD_vector_Gate(2+VR1:3:VR1+3*VR2,Newton);
    hole(1+VR1:1:VR2+VR1,1) = hole(1+VR1:1:VR2+VR1,1) + update_DD_vector_Gate(3+VR1:3:VR1+3*VR2,Newton);

    update_DD_poisson_Gate(:,Newton) = abs(update_DD_vector_Gate(1+VR1:3:VR1+3*VR2,Newton));

    if norm(abs(update_DD_vector_Gate(1+VR1:3:VR1+3*VR2,Newton)),inf) < 1e-15
        break;
    end

end
