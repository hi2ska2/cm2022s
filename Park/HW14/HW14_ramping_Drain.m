%%%%%%%% Drain 1V ramping %%%%%%%%
% In=zeros(11,1);
% Ip=zeros(11,1);
% Vg=1.0;

for bias = 0:10
    Vd = bias*0.1;
    
    for Newton = 1:20
    
        Jaco = sparse(T_row,T_row);
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
                if i <= 10 % Bottom Gate
                    V = find(Table(1:VR1,2)==Contact(i,n));

                    Jaco(V,:) =0;
                    Jaco(V,V) = 1;
                    res(V,1)=phi(V,1)-0.33374;

                else % Top Gate

                    V = find(Table(1+VR1+3*VR2:T_row,2)==Contact(i,n));

                    Jaco(VR1+3*VR2+V,:) =0;
                    Jaco(VR1+3*VR2+V,VR1+3*VR2+V) = 1;
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

            Jaco(3*V+(VR1-2),:) = 0;
            Jaco(3*V+(VR1-2),3*V+(VR1-2)) = 1;
            Jaco(3*V+(VR1-1),:) = 0;
            Jaco(3*V+(VR1-1),3*V+(VR1-1)) = 1;
            Jaco(3*V+(VR1),:) = 0;
            Jaco(3*V+(VR1),3*V+(VR1)) = 1;

            res(3*V+(VR1-2),1) = phi(VR1+V,1)-Thermal_V*log(Nd/nint)-Vd;
            res(3*V+(VR1-1),1) = elec(VR1+V,1)-Nd;
            res(3*V+(VR1),1) = hole(VR1+V,1)-(nint^2/Nd);

        end

        %     scaling
        Cvector = zeros(T_row,1);
        Cvector(1:VR1,1) = Thermal_V;
        Cvector(VR1+3*VR2:T_row,1) = Thermal_V;
        Cvector(1+VR1:3:VR1+3*VR2,1) = Thermal_V;
        Cvector(2+VR1:3:VR1+3*VR2,1) = Nd;
        Cvector(3+VR1:3:VR1+3*VR2,1) = Nd;

        Cmatrix = spdiags(Cvector,0,T_row,T_row);
        Jaco_scaled = Jaco * Cmatrix;
        Rvector = 1./sum(abs(Jaco_scaled),2);
        Rmatrix = spdiags(Rvector,0,T_row,T_row);
        Jaco_scaled = Rmatrix* Jaco_scaled;
        res_scaled = Rmatrix *res;
        update_scaled = Jaco_scaled \ (-res_scaled);
        update_DD_vector_Drain(:,Newton) = Cmatrix* update_scaled;

        phi(1:VR1,1) = phi(1:VR1,1) + update_DD_vector_Drain(1:VR1,Newton);
        phi(1+VR1:VR2+VR1,1) = phi(1+VR1:VR2+VR1,1) + update_DD_vector_Drain(1+VR1:3:VR1+3*VR2,Newton);
        phi(1+VR1+VR2:T,1) = phi(1+VR1+VR2:T,1) + update_DD_vector_Drain(1+VR1+3*VR2:T_row,Newton);

        elec(1+VR1:1:VR2+VR1,1) = elec(1+VR1:1:VR2+VR1,1) + update_DD_vector_Drain(2+VR1:3:VR1+3*VR2,Newton);
        hole(1+VR1:1:VR2+VR1,1) = hole(1+VR1:1:VR2+VR1,1) + update_DD_vector_Drain(3+VR1:3:VR1+3*VR2,Newton);

        update_DD_Phi_Drain(:,Newton) = abs(update_DD_vector_Drain(1+VR1:3:VR1+3*VR2,Newton));

        if norm(abs(update_DD_vector_Drain(1+VR1:3:VR1+3*VR2,Newton)),inf) < 1e-15
            break;
        end

    end

%     %%%%%% Current Calculation %%%%%
%     w_row = size(ry2,1);
%     Jn = zeros(bias,1);
%     Jp = zeros(bias,1);
% 
%     for tt=1:w_row
% 
%         K = find(E_R2==ry2(tt,1));
%         K_row = size(K,1);
% 
%         for rr = 1:K_row
% 
%             if K(rr,1) <= R_row2
%                 jj = 1;
% 
%                 V3 = E_R2(K(rr,1),jj+2) ; % node 3
%                 V2 = E_R2(K(rr,1),jj+1); % node 2
%                 V1 = E_R2(K(rr,1),jj); %node 1
% 
%                 x21 = (phi(VR1+V2,1)-phi(VR1+V1,1))/Thermal_V;
%                 x12 = -x21;
%                 x31 = (phi(VR1+V3,1)-phi(VR1+V1,1))/Thermal_V;
%                 x13 = -x31;
% 
%                 Jn(bias,1) = Jn(bias,1) + coeff_Jn*(elec(VR1+V2,1)*Ber(x21) - elec(VR1+V1,1)*Ber(x12))*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj));
%                 Jn(bias,1) = Jn(bias,1) + coeff_Jn*(elec(VR1+V3,1)*Ber(x31) - elec(VR1+V1,1)*Ber(x13))*(Area(K(rr,1)+R_row1,jj+2)/L(K(rr,1)+R_row1,jj+2));
%                 Jp(bias,1) = Jp(bias,1) + coeff_Jp*(hole(VR1+V2,1)*Ber(x12) - hole(VR1+V1,1)*Ber(x21))*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj));
%                 Jp(bias,1) = Jp(bias,1) + coeff_Jp*(hole(VR1+V3,1)*Ber(x13) - hole(VR1+V1,1)*Ber(x31))*(Area(K(rr,1)+R_row1,jj+2)/L(K(rr,1)+R_row1,jj+2));
% 
%             elseif K(rr,1) > R_row2 && K(rr,1) <= 2*R_row2
% 
%                 jj = 2;
%                 K(rr,1) = K(rr,1)-R_row2;
% 
%                 V3 = E_R2(K(rr,1),jj+1) ; % node 3
%                 V2 = E_R2(K(rr,1),jj); % node 2
%                 V1 = E_R2(K(rr,1),jj-1); % node 1
% 
%                 x12 = (phi(VR1+V1,1)-phi(VR1+V2,1))/Thermal_V;
%                 x21 = -x12;
%                 x32 = (phi(VR1+V3,1)-phi(VR1+V2,1))/Thermal_V;
%                 x23 = -x32;
% 
%                 Jn(bias,1) =  Jn(bias,1) + coeff_Jn*((Area(K(rr,1)+R_row1,jj-1)/L(K(rr,1)+R_row1,jj-1))*(elec(VR1+V1,1)*Ber(x12) - elec(VR1+V2,1)*Ber(x21)));
%                 Jn(bias,1) =  Jn(bias,1) + coeff_Jn*((Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj))*(elec(VR1+V3,1)*Ber(x32) - elec(VR1+V2,1)*Ber(x23)));
%                 Jp(bias,1) =  Jp(bias,1) + coeff_Jp*((Area(K(rr,1)+R_row1,jj-1)/L(K(rr,1)+R_row1,jj-1))*(hole(VR1+V1,1)*Ber(x21) - hole(VR1+V2,1)*Ber(x12)));
%                 Jp(bias,1) =  Jp(bias,1) + coeff_Jp*((Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj))*(hole(VR1+V3,1)*Ber(x23) - hole(VR1+V2,1)*Ber(x32)));
% 
%             elseif K(rr,1) > 2*R_row2 && K(rr,1) <= 3*R_row2
% 
%                 jj = 3;
%                 K(rr,1) = K(rr,1)-2*R_row2;
% 
%                 V1 = E_R2(K(rr,1),jj-2) ; % node 1
%                 V2 = E_R2(K(rr,1),jj-1); % node 2
%                 V3 = E_R2(K(rr,1),jj); % node 3
% 
%                 x13 = (phi(VR1+V1,1)-phi(VR1+V3,1))/Thermal_V;
%                 x31 = -x13;
%                 x23 = (phi(VR1+V2,1)-phi(VR1+V3,1))/Thermal_V;
%                 x32 = -x23;
% 
%                 Jn(bias,1) = Jn(bias,1) + coeff_Jn*((Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj))*(elec(VR1+V1,1)*Ber(x13) - elec(VR1+V3,1)*Ber(x31)));
%                 Jn(bias,1) = Jn(bias,1) + coeff_Jn*((Area(K(rr,1)+R_row1,jj-1)/L(K(rr,1)+R_row1,jj-1))*(elec(VR1+V2,1)*Ber(x23) - elec(VR1+V3,1)*Ber(x32)));
%                 Jp(bias,1) = Jp(bias,1) + coeff_Jp*((Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj))*(hole(VR1+V1,1)*Ber(x31) - hole(VR1+V3,1)*Ber(x13)));
%                 Jp(bias,1) = Jp(bias,1) + coeff_Jp*((Area(K(rr,1)+R_row1,jj-1)/L(K(rr,1)+R_row1,jj-1))*(hole(VR1+V2,1)*Ber(x32) - hole(VR1+V3,1)*Ber(x23)));
%             end
% 
%             I = (Jn+Jp)/1e-6;
%         end
%     end
end



