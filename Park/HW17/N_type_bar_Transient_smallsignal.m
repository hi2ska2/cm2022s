coeff_dJn=-q*m_n;
coeff_dJp=-q*m_p;

dphi=zeros(size(phi));
delec=nint*exp(dphi);
dhole=nint*exp(dphi);

%%%%% Transient input information %%%%%

time_x = 1e-12;

Time = 1*time_x; % Total time

freq = 1/Time;
deltat = time_x/100;
Step = 2*round(Time/deltat);
t = 0:deltat:Step*deltat;
Amp=1e-3;

Jn=zeros(Step+1,1);
Jp=zeros(Step+1,1);
Jd=zeros(Step+1,1);

for T_step = 1:Step+1

    dphi_old = dphi(:,1);
    delec_old= delec(:,1);
    dhole_old = dhole(:,1);

    for Newton = 1:20

        Jaco = sparse(T_row,T_row);
        res = zeros(T_row,1);

        for ii=1:T_row

            if Table(ii,1) == 1     % Region1
                K = find(E_R1 == Table(ii,2));
                K_row = size(K,1);

%                 for rr = 1:K_row
% 
%                     if Table(ii,3) == Potential
% 
%                         if K(rr,1) <= R_row1
% 
%                             jj = 1;
%                             V3 = find(Table(1:VR1,2)==E_R1(K(rr,1),jj+2)) ;
%                             V2 = find(Table(1:VR1,2)==E_R1(K(rr,1),jj+1));
%                             V1 = find(Table(1:VR1,2)==E_R1(K(rr,1),jj));
% 
%                             res(ii,1) = res(ii,1) + eox*(phi(V2,1)-phi(V1,1))*(Area(K(rr,1),jj)/L(K(rr,1),jj)) + eox*(phi(V3,1)-phi(V1,1))*(Area(K(rr,1),jj+2)/L(K(rr,1),jj+2));
%                             Jaco(ii, V3) = Jaco(ii, V3)+ eox*(Area(K(rr,1),jj+2)/L(K(rr,1),jj+2));
%                             Jaco(ii, V2) = Jaco(ii, V2)+eox*(Area(K(rr,1),jj)/L(K(rr,1),jj));
%                             Jaco(ii, V1) = Jaco(ii, V1)-eox*(Area(K(rr,1),jj+2)/L(K(rr,1),jj+2))-eox*(Area(K(rr,1),jj)/L(K(rr,1),jj));
% 
%                         elseif K(rr,1) > R_row1 && K(rr,1) <= 2*R_row1
%                             jj = 2;
%                             K(rr,1) = K(rr,1)-R_row1;
% 
%                             V3 = find(Table(1:VR1,2)==E_R1(K(rr,1),jj+1)) ;
%                             V2 = find(Table(1:VR1,2)==E_R1(K(rr,1),jj));
%                             V1 = find(Table(1:VR1,2)==E_R1(K(rr,1),jj-1));
% 
%                             res(ii,1) = res(ii,1) + eox*(phi(V1,1)-phi(V2,1))*(Area(K(rr,1),jj-1)/L(K(rr,1),jj-1)) + eox*(phi(V3,1)-phi(V2,1))*(Area(K(rr,1),jj)/L(K(rr,1),jj));
%                             Jaco(ii, V3) = Jaco(ii, V3)+eox*(Area(K(rr,1),jj)/L(K(rr,1),jj));
%                             Jaco(ii, V2) =  Jaco(ii, V2)-eox*(Area(K(rr,1),jj)/L(K(rr,1),jj))-eox*(Area(K(rr,1),jj-1)/L(K(rr,1),jj-1));
%                             Jaco(ii, V1) = Jaco(ii, V1)+eox*(Area(K(rr,1),jj-1)/L(K(rr,1),jj-1));
% 
%                         elseif K(rr,1) > 2*R_row1 && K(rr,1) <= 3*R_row1
%                             jj = 3;
%                             K(rr,1) = K(rr,1)-2*R_row1;
% 
%                             V3 = find(Table(1:VR1,2)==E_R1(K(rr,1),jj-2)) ;
%                             V2 = find(Table(1:VR1,2)==E_R1(K(rr,1),jj-1));
%                             V1 = find(Table(1:VR1,2)==E_R1(K(rr,1),jj));
% 
%                             res(ii,1) = res(ii,1) + eox*(phi(V3,1)-phi(V1,1))*(Area(K(rr,1),jj)/L(K(rr,1),jj)) + eox*(phi(V2,1)-phi(V1,1))*(Area(K(rr,1),jj-1)/L(K(rr,1),jj-1));
%                             Jaco(ii, V3) = Jaco(ii, V3)+eox*(Area(K(rr,1),jj)/L(K(rr,1),jj));
%                             Jaco(ii, V2) = Jaco(ii, V2)+eox*(Area(K(rr,1),jj-1)/L(K(rr,1),jj-1));
%                             Jaco(ii, V1) = Jaco(ii, V1)-eox*(Area(K(rr,1),jj)/L(K(rr,1),jj))-eox*(Area(K(rr,1),jj-1)/L(K(rr,1),jj-1));
% 
%                         end
%                     end
%                 end

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

                            res(ii,1) = res(ii,1) + esi*(dphi(VR1+V2,1)-dphi(VR1+V1,1))*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj)) + esi*(dphi(VR1+V3,1)-dphi(VR1+V1,1))*(Area(K(rr,1)+R_row1,jj+2)/L(K(rr,1)+R_row1,jj+2));
                            res(ii,1) = res(ii,1) + coeff*(Area(K(rr,1)+R_row1,jj)*L(K(rr,1)+R_row1,jj)+Area(K(rr,1)+R_row1,jj+2)*L(K(rr,1)+R_row1,jj+2))*(Na1-elec(VR1+V1,1)+hole(VR1+V1,1));

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

                            res(ii,1) = res(ii,1) + esi*(dphi(VR1+V1,1)-dphi(VR1+V2,1))*(Area(K(rr,1)+R_row1,jj-1)/L(K(rr,1)+R_row1,jj-1)) + esi*(dphi(VR1+V3,1)-dphi(VR1+V2,1))*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj));
                            res(ii,1) = res(ii,1) + coeff*(Area(K(rr,1)+R_row1,jj)*L(K(rr,1)+R_row1,jj)+Area(K(rr,1)+R_row1,jj-1)*L(K(rr,1)+R_row1,jj-1))*(Na1-elec(VR1+V2,1)+hole(VR1+V2,1));

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

                            res(ii,1) = res(ii,1) + esi*(dphi(VR1+V3,1)-dphi(VR1+V1,1))*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj)) + esi*(dphi(VR1+V2,1)-dphi(VR1+V1,1))*(Area(K(rr,1)+R_row1,jj-1)/L(K(rr,1)+R_row1,jj-1));
                            res(ii,1) = res(ii,1) + coeff*(Area(K(rr,1)+R_row1,jj)*L(K(rr,1)+R_row1,jj)+Area(K(rr,1)+R_row1,jj-1)*L(K(rr,1)+R_row1,jj-1))*(Na1-elec(VR1+V1,1)+hole(VR1+V1,1));

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

                            res(ii,1) = res(ii,1) + coeff_dJn*((Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj))*((elecDC_av1*ddphi1+delec_av1*dphi_DC1)-Thermal_V*ddelec1));
                            res(ii,1) = res(ii,1) + coeff_dJn*((Area(K(rr,1)+R_row1,jj+2)/L(K(rr,1)+R_row1,jj+2))*((elecDC_av2*ddphi2+delec_av2*dphi_DC2)-Thermal_V*ddelec2));

                            % eDensity
                            Jaco(ii, 3*(V1)+(VR1-1)) = Jaco(ii, 3*(V1)+(VR1-1)) + coeff_dJn*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj))*(0.5*dphi_DC1+Thermal_V);
                            Jaco(ii, 3*(V1)+(VR1-1)) = Jaco(ii, 3*(V1)+(VR1-1)) + coeff_dJn*(Area(K(rr,1)+R_row1,jj+2)/L(K(rr,1)+R_row1,jj+2))*(0.5*dphi_DC2+Thermal_V);

                            Jaco(ii, 3*(V2)+(VR1-1)) = Jaco(ii, 3*(V2)+(VR1-1)) + coeff_dJn*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj))*(0.5*dphi_DC1-Thermal_V);
                            Jaco(ii, 3*(V3)+(VR1-1)) = Jaco(ii, 3*(V3)+(VR1-1)) + coeff_dJn*(Area(K(rr,1)+R_row1,jj+2)/L(K(rr,1)+R_row1,jj+2))*(0.5*dphi_DC2-Thermal_V);

                            % transient
                            res(ii,1) = res(ii,1) - (q/4)/deltat*(Area(K(rr,1)+R_row1,jj)*L(K(rr,1)+R_row1,jj)+Area(K(rr,1)+R_row1,jj+2)*L(K(rr,1)+R_row1,jj+2))*(delec(VR1+V1,1)-delec_old(VR1+V1,1));

                            Jaco(ii,3*(V1)+(VR1-1)) = Jaco(ii,3*(V1)+(VR1-1)) - (q/4)/deltat*(Area(K(rr,1)+R_row1,jj)*L(K(rr,1)+R_row1,jj)+Area(K(rr,1)+R_row1,jj+2)*L(K(rr,1)+R_row1,jj+2));

                            % Potential
                            Jaco(ii, 3*(V1)+(VR1-2)) = Jaco(ii, 3*(V1)+(VR1-2)) + coeff_dJn*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj))*(-elecDC_av1);
                            Jaco(ii, 3*(V1)+(VR1-2)) = Jaco(ii, 3*(V1)+(VR1-2)) + coeff_dJn*(Area(K(rr,1)+R_row1,jj+2)/L(K(rr,1)+R_row1,jj+2))*(-elecDC_av2);

                            Jaco(ii, 3*(V2)+(VR1-2)) = Jaco(ii, 3*(V2)+(VR1-2)) + coeff_dJn*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj))*elecDC_av1;
                            Jaco(ii, 3*(V3)+(VR1-2)) = Jaco(ii, 3*(V3)+(VR1-2)) + coeff_dJn*(Area(K(rr,1)+R_row1,jj+2)/L(K(rr,1)+R_row1,jj+2))*elecDC_av2;


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

                            res(ii,1) =  res(ii,1) + coeff_dJn*((Area(K(rr,1)+R_row1,jj-1)/L(K(rr,1)+R_row1,jj-1))*((elecDC_av1*ddphi1+delec_av1*dphi_DC1)-Thermal_V*ddelec1));
                            res(ii,1) =  res(ii,1) + coeff_dJn*((Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj))*((elecDC_av2*ddphi2+delec_av2*dphi_DC2)-Thermal_V*ddelec2));

                            % eDensity
                            Jaco(ii, 3*(V2)+(VR1-1)) = Jaco(ii, 3*(V2)+(VR1-1)) + coeff_dJn*(Area(K(rr,1)+R_row1,jj-1)/L(K(rr,1)+R_row1,jj-1))*(0.5*dphi_DC1+Thermal_V);
                            Jaco(ii, 3*(V2)+(VR1-1)) = Jaco(ii, 3*(V2)+(VR1-1)) + coeff_dJn*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj))*(0.5*dphi_DC2+Thermal_V);

                            Jaco(ii, 3*(V1)+(VR1-1)) = Jaco(ii, 3*(V1)+(VR1-1)) + coeff_dJn*(Area(K(rr,1)+R_row1,jj-1)/L(K(rr,1)+R_row1,jj-1))*(0.5*dphi_DC1-Thermal_V);
                            Jaco(ii, 3*(V3)+(VR1-1)) = Jaco(ii, 3*(V3)+(VR1-1)) + coeff_dJn*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj))*(0.5*dphi_DC2-Thermal_V);

                            % transient
                            res(ii,1) = res(ii,1) - (q/4)/deltat*(Area(K(rr,1)+R_row1,jj-1)*L(K(rr,1)+R_row1,jj-1)+Area(K(rr,1)+R_row1,jj)*L(K(rr,1)+R_row1,jj))*(delec(VR1+V2,1)-delec_old(VR1+V2,1));

                            Jaco(ii,3*(V2)+(VR1-1)) = Jaco(ii,3*(V2)+(VR1-1)) - (q/4)/deltat*(Area(K(rr,1)+R_row1,jj-1)*L(K(rr,1)+R_row1,jj-1)+Area(K(rr,1)+R_row1,jj)*L(K(rr,1)+R_row1,jj));

                            % Potential
                            Jaco(ii, 3*(V2)+(VR1-2)) = Jaco(ii, 3*(V2)+(VR1-2)) + coeff_dJn*(Area(K(rr,1)+R_row1,jj-1)/L(K(rr,1)+R_row1,jj-1))*(-elecDC_av1);
                            Jaco(ii, 3*(V2)+(VR1-2)) = Jaco(ii, 3*(V2)+(VR1-2)) + coeff_dJn*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj))*(-elecDC_av2);

                            Jaco(ii, 3*(V1)+(VR1-2)) = Jaco(ii, 3*(V1)+(VR1-2)) + coeff_dJn*(Area(K(rr,1)+R_row1,jj-1)/L(K(rr,1)+R_row1,jj-1))*elecDC_av1;
                            Jaco(ii, 3*(V3)+(VR1-2)) = Jaco(ii, 3*(V3)+(VR1-2)) + coeff_dJn*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj))*elecDC_av2;

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

                            res(ii,1) = res(ii,1) + coeff_dJn*((Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj))*((elecDC_av1*ddphi1+delec_av1*dphi_DC1)-Thermal_V*ddelec1));
                            res(ii,1) = res(ii,1) + coeff_dJn*((Area(K(rr,1)+R_row1,jj-1)/L(K(rr,1)+R_row1,jj-1))*((elecDC_av2*ddphi2+delec_av2*dphi_DC2)-Thermal_V*ddelec2));

                            % eDensity
                            Jaco(ii, 3*(V1)+(VR1-1)) = Jaco(ii, 3*(V1)+(VR1-1)) + coeff_dJn*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj))*(0.5*dphi_DC1-Thermal_V);
                            Jaco(ii, 3*(V2)+(VR1-1)) = Jaco(ii, 3*(V2)+(VR1-1)) + coeff_dJn*(Area(K(rr,1)+R_row1,jj-1)/L(K(rr,1)+R_row1,jj-1))*(0.5*dphi_DC2-Thermal_V);

                            Jaco(ii, 3*(V3)+(VR1-1)) = Jaco(ii, 3*(V3)+(VR1-1)) + coeff_dJn*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj))*(0.5*dphi_DC1+Thermal_V);
                            Jaco(ii, 3*(V3)+(VR1-1)) = Jaco(ii, 3*(V3)+(VR1-1)) + coeff_dJn*(Area(K(rr,1)+R_row1,jj-1)/L(K(rr,1)+R_row1,jj-1))*(0.5*dphi_DC2+Thermal_V);

                            % transient
                            res(ii,1) = res(ii,1) - (q/4)/deltat*(Area(K(rr,1)+R_row1,jj)*L(K(rr,1)+R_row1,jj)+Area(K(rr,1)+R_row1,jj-1)*L(K(rr,1)+R_row1,jj-1))*(delec(VR1+V3,1)-delec_old(VR1+V3,1));

                            Jaco(ii,3*(V3)+(VR1-1)) = Jaco(ii,3*(V3)+(VR1-1)) - (q/4)/deltat*(Area(K(rr,1)+R_row1,jj)*L(K(rr,1)+R_row1,jj)+Area(K(rr,1)+R_row1,jj-1)*L(K(rr,1)+R_row1,jj-1));

                            % Potential
                            Jaco(ii, 3*(V3)+(VR1-2)) = Jaco(ii, 3*(V3)+(VR1-2)) + coeff_dJn*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj))*(-elecDC_av1);
                            Jaco(ii, 3*(V3)+(VR1-2)) = Jaco(ii, 3*(V3)+(VR1-2)) + coeff_dJn*(Area(K(rr,1)+R_row1,jj-1)/L(K(rr,1)+R_row1,jj-1))*(-elecDC_av2);

                            Jaco(ii, 3*(V2)+(VR1-2)) = Jaco(ii, 3*(V2)+(VR1-2)) + coeff_dJn*(Area(K(rr,1)+R_row1,jj-1)/L(K(rr,1)+R_row1,jj-1))*(elecDC_av2);
                            Jaco(ii, 3*(V1)+(VR1-2)) = Jaco(ii, 3*(V1)+(VR1-2)) + coeff_dJn*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj))*(elecDC_av1);
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

                            res(ii,1) = res(ii,1) + coeff_dJp*((Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj))*((holeDC_av1*ddphi1+dhole_av1*dphi_DC1)+Thermal_V*ddhole1));
                            res(ii,1) = res(ii,1) + coeff_dJp*((Area(K(rr,1)+R_row1,jj+2)/L(K(rr,1)+R_row1,jj+2))*((holeDC_av2*ddphi2+dhole_av2*dphi_DC2)+Thermal_V*ddhole2));

                            %hDensity
                            Jaco(ii, 3*(V1)+(VR1)) = Jaco(ii, 3*(V1)+(VR1)) + coeff_dJp*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj))*(0.5*dphi_DC1-Thermal_V);
                            Jaco(ii, 3*(V1)+(VR1)) = Jaco(ii, 3*(V1)+(VR1)) + coeff_dJp*(Area(K(rr,1)+R_row1,jj+2)/L(K(rr,1)+R_row1,jj+2))*(0.5*dphi_DC2-Thermal_V);

                            Jaco(ii, 3*(V2)+(VR1)) = Jaco(ii, 3*(V2)+(VR1)) + coeff_dJp*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj))*(0.5*dphi_DC1+Thermal_V);
                            Jaco(ii, 3*(V3)+(VR1)) = Jaco(ii, 3*(V3)+(VR1)) + coeff_dJp*(Area(K(rr,1)+R_row1,jj+2)/L(K(rr,1)+R_row1,jj+2))*(0.5*dphi_DC2+Thermal_V);

                            % transient
                            res(ii,1) = res(ii,1) + (q/4)/deltat*(Area(K(rr,1)+R_row1,jj)*L(K(rr,1)+R_row1,jj)+Area(K(rr,1)+R_row1,jj+2)*L(K(rr,1)+R_row1,jj+2))*(dhole(VR1+V1,1)-dhole_old(VR1+V1,1));

                            Jaco(ii,3*(V1)+(VR1)) = Jaco(ii,3*(V1)+(VR1)) + (q/4)/deltat*(Area(K(rr,1)+R_row1,jj)*L(K(rr,1)+R_row1,jj)+Area(K(rr,1)+R_row1,jj+2)*L(K(rr,1)+R_row1,jj+2));

                            %Potential
                            Jaco(ii, 3*(V1)+(VR1-2)) = Jaco(ii, 3*(V1)+(VR1-2)) + coeff_dJp*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj))*(-holeDC_av1);
                            Jaco(ii, 3*(V1)+(VR1-2)) = Jaco(ii, 3*(V1)+(VR1-2)) + coeff_dJp*(Area(K(rr,1)+R_row1,jj+2)/L(K(rr,1)+R_row1,jj+2))*(-holeDC_av2);

                            Jaco(ii, 3*(V2)+(VR1-2)) = Jaco(ii, 3*(V2)+(VR1-2)) + coeff_dJp*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj))*holeDC_av1;
                            Jaco(ii, 3*(V3)+(VR1-2)) = Jaco(ii, 3*(V3)+(VR1-2)) + coeff_dJp*(Area(K(rr,1)+R_row1,jj+2)/L(K(rr,1)+R_row1,jj+2))*holeDC_av2;


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

                            res(ii,1) =  res(ii,1) + coeff_dJp*((Area(K(rr,1)+R_row1,jj-1)/L(K(rr,1)+R_row1,jj-1))*((holeDC_av1*ddphi1+dhole_av1*dphi_DC1)+Thermal_V*ddhole1));
                            res(ii,1) =  res(ii,1) + coeff_dJp*((Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj))*((holeDC_av2*ddphi2+dhole_av2*dphi_DC2)+Thermal_V*ddhole2));

                            % hDensity
                            Jaco(ii, 3*(V1)+(VR1)) = Jaco(ii, 3*(V1)+(VR1)) + coeff_dJp*(Area(K(rr,1)+R_row1,jj-1)/L(K(rr,1)+R_row1,jj-1))*(0.5*dphi_DC1+Thermal_V);
                            Jaco(ii, 3*(V3)+(VR1)) = Jaco(ii, 3*(V3)+(VR1)) + coeff_dJp*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj))*(0.5*dphi_DC2+Thermal_V);

                            Jaco(ii, 3*(V2)+(VR1)) = Jaco(ii, 3*(V2)+(VR1)) + coeff_dJp*(Area(K(rr,1)+R_row1,jj-1)/L(K(rr,1)+R_row1,jj-1))*(0.5*dphi_DC1-Thermal_V);
                            Jaco(ii, 3*(V2)+(VR1)) = Jaco(ii, 3*(V2)+(VR1)) + coeff_dJp*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj))*(0.5*dphi_DC2-Thermal_V);

                            % transient
                            res(ii,1) = res(ii,1) + (q/4)/deltat*(Area(K(rr,1)+R_row1,jj-1)*L(K(rr,1)+R_row1,jj-1)+Area(K(rr,1)+R_row1,jj)*L(K(rr,1)+R_row1,jj))*(dhole(VR1+V2,1)-dhole_old(VR1+V2,1));

                            Jaco(ii,3*(V2)+(VR1)) = Jaco(ii,3*(V2)+(VR1)) + (q/4)/deltat*(Area(K(rr,1)+R_row1,jj-1)*L(K(rr,1)+R_row1,jj-1)+Area(K(rr,1)+R_row1,jj)*L(K(rr,1)+R_row1,jj));

                            %Potential
                            Jaco(ii, 3*(V2)+(VR1-2)) = Jaco(ii, 3*(V2)+(VR1-2)) + coeff_dJp*(Area(K(rr,1)+R_row1,jj-1)/L(K(rr,1)+R_row1,jj-1))*(-holeDC_av1);
                            Jaco(ii, 3*(V2)+(VR1-2)) = Jaco(ii, 3*(V2)+(VR1-2)) + coeff_dJp*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj))*(-holeDC_av2);

                            Jaco(ii, 3*(V1)+(VR1-2)) = Jaco(ii, 3*(V1)+(VR1-2)) + coeff_dJp*(Area(K(rr,1)+R_row1,jj-1)/L(K(rr,1)+R_row1,jj-1))*(holeDC_av1);
                            Jaco(ii, 3*(V3)+(VR1-2)) = Jaco(ii, 3*(V3)+(VR1-2)) + coeff_dJp*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj))*(holeDC_av2);


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

                            res(ii,1) = res(ii,1) + coeff_dJp*((Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj))*((holeDC_av1*ddphi1+dhole_av1*dphi_DC1)+Thermal_V*ddhole1));
                            res(ii,1) = res(ii,1) + coeff_dJp*((Area(K(rr,1)+R_row1,jj-1)/L(K(rr,1)+R_row1,jj-1))*((holeDC_av2*ddphi2+dhole_av2*dphi_DC2)+Thermal_V*ddhole2));

                            %hDensity
                            Jaco(ii, 3*(V1)+(VR1)) = Jaco(ii, 3*(V1)+(VR1)) + coeff_dJp*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj))*(0.5*dphi_DC1+Thermal_V);
                            Jaco(ii, 3*(V2)+(VR1)) = Jaco(ii, 3*(V2)+(VR1)) + coeff_dJp*(Area(K(rr,1)+R_row1,jj-1)/L(K(rr,1)+R_row1,jj-1))*(0.5*dphi_DC2+Thermal_V);

                            Jaco(ii, 3*(V3)+(VR1)) = Jaco(ii, 3*(V3)+(VR1)) + coeff_dJp*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj))*(0.5*dphi_DC1-Thermal_V);
                            Jaco(ii, 3*(V3)+(VR1)) = Jaco(ii, 3*(V3)+(VR1)) + coeff_dJp*(Area(K(rr,1)+R_row1,jj-1)/L(K(rr,1)+R_row1,jj-1))*(0.5*dphi_DC2-Thermal_V);

                            % transient
                            res(ii,1) = res(ii,1) + (q/4)/deltat*(Area(K(rr,1)+R_row1,jj)*L(K(rr,1)+R_row1,jj)+Area(K(rr,1)+R_row1,jj-1)*L(K(rr,1)+R_row1,jj-1))*(dhole(VR1+V3,1)-dhole_old(VR1+V3,1));

                            Jaco(ii,3*(V3)+(VR1)) = Jaco(ii,3*(V3)+(VR1)) + (q/4)/deltat*(Area(K(rr,1)+R_row1,jj)*L(K(rr,1)+R_row1,jj)+Area(K(rr,1)+R_row1,jj-1)*L(K(rr,1)+R_row1,jj-1));

                            %Potential
                            Jaco(ii, 3*(V3)+(VR1-2)) = Jaco(ii, 3*(V3)+(VR1-2)) + coeff_dJp*(Area(K(rr,1)+R_row1,jj-1)/L(K(rr,1)+R_row1,jj-1))*(-holeDC_av1);
                            Jaco(ii, 3*(V3)+(VR1-2)) = Jaco(ii, 3*(V3)+(VR1-2)) + coeff_dJp*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj))*(-holeDC_av2);

                            Jaco(ii, 3*(V2)+(VR1-2)) = Jaco(ii, 3*(V2)+(VR1-2)) + coeff_dJp*(Area(K(rr,1)+R_row1,jj-1)/L(K(rr,1)+R_row1,jj-1))*holeDC_av1;
                            Jaco(ii, 3*(V1)+(VR1-2)) = Jaco(ii, 3*(V1)+(VR1-2)) + coeff_dJp*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj))*holeDC_av2;
                        end
                    end
                end

            elseif Table(ii,1) == 3  % Region3

%                 K = find(E_R3 == Table(ii,2));
%                 K_row = size(K,1);
% 
%                 for rr = 1:K_row
%                     if Table(ii,3) == Potential
% 
%                         if K(rr,1) <= R_row3
%                             jj = 1;
% 
%                             V3 = find(Table(1+VR1+3*VR2:T1,2)==E_R3(K(rr,1),jj+2)) ;
%                             V2 = find(Table(1+VR1+3*VR2:T1,2)==E_R3(K(rr,1),jj+1));
%                             V1 = find(Table(1+VR1+3*VR2:T1,2)==E_R3(K(rr,1),jj));
% 
%                             res(ii,1) = res(ii,1) + eox*(phi(VR1+VR2+V2,1)-phi(VR1+VR2+V1,1))*(Area(K(rr,1)+R_row2+R_row1,jj)/L(K(rr,1)+R_row2+R_row1,jj)) + eox*(phi(VR1+VR2+V3,1)-phi(VR1+VR2+V1,1))*(Area(K(rr,1)+R_row2+R_row1,jj+2)/L(K(rr,1)+R_row2+R_row1,jj+2));
% 
%                             Jaco(ii, VR1+3*VR2+V3) = Jaco(ii, VR1+3*VR2+V3)+ eox*(Area(K(rr,1)+R_row2+R_row1,jj+2)/L(K(rr,1)+R_row2+R_row1,jj+2));
%                             Jaco(ii, VR1+3*VR2+V2) = Jaco(ii, VR1+3*VR2+V2)+eox*(Area(K(rr,1)+R_row2+R_row1,jj)/L(K(rr,1)+R_row2+R_row1,jj));
%                             Jaco(ii, VR1+3*VR2+V1) = Jaco(ii, VR1+3*VR2+V1)-eox*(Area(K(rr,1)+R_row2+R_row1,jj+2)/L(K(rr,1)+R_row2+R_row1,jj+2))-eox*(Area(K(rr,1)+R_row2+R_row1,jj)/L(K(rr,1)+R_row2+R_row1,jj));
% 
%                         elseif K(rr,1) > R_row3 && K(rr,1) <= 2*R_row3
% 
%                             jj = 2;
%                             K(rr,1) = K(rr,1)-R_row3;
% 
%                             V3 = find(Table(1+VR1+3*VR2:T1,2)==E_R3(K(rr,1),jj+1)) ;
%                             V2 = find(Table(1+VR1+3*VR2:T1,2)==E_R3(K(rr,1),jj));
%                             V1 = find(Table(1+VR1+3*VR2:T1,2)==E_R3(K(rr,1),jj-1));
% 
%                             res(ii,1) = res(ii,1) + eox*(phi(VR1+VR2+V1,1)-phi(VR1+VR2+V2,1))*(Area(K(rr,1)+R_row2+R_row1,jj-1)/L(K(rr,1)+R_row2+R_row1,jj-1)) + eox*(phi(VR1+VR2+V3,1)-phi(VR1+VR2+V2,1))*(Area(K(rr,1)+R_row2+R_row1,jj)/L(K(rr,1)+R_row2+R_row1,jj));
% 
%                             Jaco(ii, VR1+3*VR2+V3) = Jaco(ii, VR1+3*VR2+V3)+eox*(Area(K(rr,1)+R_row2+R_row1,jj)/L(K(rr,1)+R_row2+R_row1,jj));
%                             Jaco(ii, VR1+3*VR2+V2) =  Jaco(ii, VR1+3*VR2+V2)-eox*(Area(K(rr,1)+R_row2+R_row1,jj)/L(K(rr,1)+R_row2+R_row1,jj))-eox*(Area(K(rr,1)+R_row2+R_row1,jj-1)/L(K(rr,1)+R_row2+R_row1,jj-1));
%                             Jaco(ii, VR1+3*VR2+V1) = Jaco(ii, VR1+3*VR2+V1)+eox*(Area(K(rr,1)+R_row2+R_row1,jj-1)/L(K(rr,1)+R_row2+R_row1,jj-1));
% 
%                         elseif K(rr,1) > 2*R_row3 && K(rr,1) <= 3*R_row3
% 
%                             jj = 3;
%                             K(rr,1) = K(rr,1)-2*R_row3;
% 
%                             V3 = find(Table(1+VR1+3*VR2:T1,2)==E_R3(K(rr,1),jj-2)) ;
%                             V2 = find(Table(1+VR1+3*VR2:T1,2)==E_R3(K(rr,1),jj-1));
%                             V1 = find(Table(1+VR1+3*VR2:T1,2)==E_R3(K(rr,1),jj));
% 
%                             res(ii,1) = res(ii,1) + eox*(phi(VR1+VR2+V3,1)-phi(VR1+VR2+V1,1))*(Area(K(rr,1)+R_row2+R_row1,jj)/L(K(rr,1)+R_row2+R_row1,jj)) + eox*(phi(VR1+VR2+V2,1)-phi(VR1+VR2+V1,1))*(Area(K(rr,1)+R_row2+R_row1,jj-1)/L(K(rr,1)+R_row2+R_row1,jj-1));
% 
%                             Jaco(ii, VR1+3*VR2+V3) = Jaco(ii, VR1+3*VR2+V3) + eox*(Area(K(rr,1)+R_row2+R_row1,jj)/L(K(rr,1)+R_row2+R_row1,jj));
%                             Jaco(ii, VR1+3*VR2+V2) = Jaco(ii, VR1+3*VR2+V2) + eox*(Area(K(rr,1)+R_row2+R_row1,jj-1)/L(K(rr,1)+R_row2+R_row1,jj-1));
%                             Jaco(ii, VR1+3*VR2+V1) = Jaco(ii, VR1+3*VR2+V1) - eox*(Area(K(rr,1)+R_row2+R_row1,jj)/L(K(rr,1)+R_row2+R_row1,jj))-eox*(Area(K(rr,1)+R_row2+R_row1,jj-1)/L(K(rr,1)+R_row2+R_row1,jj-1));
% 
%                         end
%                     end
%                 end
            end
        end

        %%%% Anode, Cathode %%%%
        q_row = size(ry1,1);
        w_row = size(ry2,1);

        for rr = 1:q_row  % Cathode

            V = find(Table(1+VR1:3:VR1+3*VR2,2)==ry1(rr,1));
            Jaco(3*V+(VR1-2),:) = 0;
            Jaco(3*V+(VR1-2),3*V+(VR1-2)) = 1;
            Jaco(3*V+(VR1-1),:) = 0;
            Jaco(3*V+(VR1-1),3*V+(VR1-1)) = 1;
            Jaco(3*V+(VR1),:) = 0;
            Jaco(3*V+(VR1),3*V+(VR1)) = 1;

            res(3*V+(VR1-2),1) = dphi(VR1+V,1);%-Thermal_V*log(-Na/nint);
            res(3*V+(VR1-1),1) = delec(VR1+V,1);
            res(3*V+(VR1),1) = dhole(VR1+V,1);

        end

        for rr = 1:w_row  % Anode

            V = find(Table(1+VR1:3:VR1+3*VR2,2)==ry2(rr,1));

            Jaco(3*V+(VR1-2),:) = 0;
            Jaco(3*V+(VR1-2),3*V+(VR1-2)) = 1;
            Jaco(3*V+(VR1-1),:) = 0;
            Jaco(3*V+(VR1-1),3*V+(VR1-1)) = 1;
            Jaco(3*V+(VR1),:) = 0;
            Jaco(3*V+(VR1),3*V+(VR1)) = 1;

            res(3*V+(VR1-2),1) = dphi(VR1+V,1)-Amp*sin(2*pi*freq*t(:,T_step));
            res(3*V+(VR1-1),1) = delec(VR1+V,1);
            res(3*V+(VR1),1) = dhole(VR1+V,1);

        end

        %     scaling
        Cvector = zeros(T_row,1);
        Cvector(1:VR1,1) = Thermal_V;
        Cvector(VR1+3*VR2:T_row,1) = Thermal_V;
        Cvector(1+VR1:3:VR1+3*VR2,1) = Thermal_V;
        Cvector(2+VR1:3:VR1+3*VR2,1) = -Na;
        Cvector(3+VR1:3:VR1+3*VR2,1) = -Na;

        Cmatrix = spdiags(Cvector,0,T_row,T_row);
        Jaco_scaled = Jaco * Cmatrix;
        Rvector = 1./sum(abs(Jaco_scaled),2);
        Rmatrix = spdiags(Rvector,0,T_row,T_row);
        Jaco_scaled = Rmatrix* Jaco_scaled;
        res_scaled = Rmatrix *res;
        update_scaled = Jaco_scaled \ (-res_scaled);
        update_DD_vector_tran(:,Newton) = Cmatrix* update_scaled;

        dphi(1:VR1,1) = dphi(1:VR1,1) + update_DD_vector_tran(1:VR1,Newton);
        dphi(1+VR1:VR2+VR1,1) = dphi(1+VR1:VR2+VR1,1) + update_DD_vector_tran(1+VR1:3:VR1+3*VR2,Newton);
        dphi(1+VR1+VR2:T,1) = dphi(1+VR1+VR2:T,1) + update_DD_vector_tran(1+VR1+3*VR2:T_row,Newton);

        delec(1+VR1:1:VR2+VR1,1) = delec(1+VR1:1:VR2+VR1,1) + update_DD_vector_tran(2+VR1:3:VR1+3*VR2,Newton);
        dhole(1+VR1:1:VR2+VR1,1) = dhole(1+VR1:1:VR2+VR1,1) + update_DD_vector_tran(3+VR1:3:VR1+3*VR2,Newton);

        update_DD_poisson_tran(:,Newton) = abs(update_DD_vector_tran(1+VR1:3:VR1+3*VR2,Newton));

        if norm(abs(update_DD_poisson_tran(:,Newton)),inf) < 1e-15
            break;
        end
    end

    %%%%%% Current Calculation %%%%%
    w_row = size(ry2,1);

    for tt=1:w_row

        K = find(E_R2==ry2(tt,1));
        K_row = size(K,1);

        for rr = 1:K_row

            if K(rr,1) <= R_row2
                jj = 1;

                V3 = find(Table(1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj+2)) ; % node 3
                V2 = find(Table(1+1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj+1)); % node 2
                V1 = find(Table(2+1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj)); %node 1

                elecDC_av1 = 0.5*(elec_DC(VR1+V2,1)+elec_DC(VR1+V1,1));
                delec_av1 = 0.5*(delec(VR1+V2,1)+delec(VR1+V1,1));
                ddelec1 = delec(VR1+V2,1)-delec(VR1+V1,1);

                elecDC_av2 = 0.5*(elec_DC(VR1+V3,1)+elec_DC(VR1+V1,1));
                delec_av2 = 0.5*(delec(VR1+V3,1)+delec(VR1+V1,1));
                ddelec2 = delec(VR1+V3,1)-delec(VR1+V1,1);

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

                Jn(T_step,1) = Jn(T_step,1) + coeff_dJn*((elecDC_av1*ddphi1+delec_av1*dphi_DC1)-Thermal_V*ddelec1)*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj));
                Jn(T_step,1) = Jn(T_step,1) + coeff_dJn*((elecDC_av2*ddphi2+delec_av2*dphi_DC2)-Thermal_V*ddelec2)*(Area(K(rr,1)+R_row1,jj+2)/L(K(rr,1)+R_row1,jj+2));
                Jp(T_step,1) = Jp(T_step,1) + coeff_dJp*((holeDC_av1*ddphi1+dhole_av1*dphi_DC1)+Thermal_V*ddhole1)*(Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj));
                Jp(T_step,1) = Jp(T_step,1) + coeff_dJp*((holeDC_av2*ddphi2+dhole_av2*dphi_DC2)+Thermal_V*ddhole2)*(Area(K(rr,1)+R_row1,jj+2)/L(K(rr,1)+R_row1,jj+2));
                Jd(T_step,1) = Jd(T_step,1) - esi*e0*((dphi(VR1+V3,1)-dphi(VR1+V1,1))-(dphi_old(VR1+V3,1)-dphi_old(VR1+V1,1)))*Area(K(rr,1)+R_row1,jj+2)/L(K(rr,1)+R_row1,jj+2);
                Jd(T_step,1) = Jd(T_step,1) - esi*e0/deltat*((dphi(VR1+V2,1)-dphi(VR1+V1,1))-(dphi_old(VR1+V2,1)-dphi_old(VR1+V1,1)))*Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj);

            elseif K(rr,1) > R_row2 && K(rr,1) <= 2*R_row2

                jj = 2;
                K(rr,1) = K(rr,1)-R_row2;

                V3 = find(Table(1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj+1)) ; % node 3
                V2 = find(Table(1+1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj)); % node 2
                V1 = find(Table(2+1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj-1)); %node 1

                elecDC_av1 = 0.5*(elec_DC(VR1+V1,1)+elec_DC(VR1+V2,1));
                delec_av1 = 0.5*(delec(VR1+V1,1)+delec(VR1+V2,1));
                ddelec1 = delec(VR1+V1,1)-delec(VR1+V2,1);

                elecDC_av2 = 0.5*(elec_DC(VR1+V3,1)+elec_DC(VR1+V2,1));
                delec_av2 = 0.5*(delec(VR1+V3,1)+delec(VR1+V2,1));
                ddelec2 = delec(VR1+V3,1)-delec(VR1+V2,1);

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

                Jn(T_step,1) = Jn(T_step,1) + coeff_dJn*((Area(K(rr,1)+R_row1,jj-1)/L(K(rr,1)+R_row1,jj-1))*((elecDC_av1*ddphi1+delec_av1*dphi_DC1)-Thermal_V*ddelec1));
                Jn(T_step,1) = Jn(T_step,1) + coeff_dJn*((Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj))*((elecDC_av2*ddphi2+delec_av2*dphi_DC2)-Thermal_V*ddelec2));
                Jp(T_step,1) = Jp(T_step,1) + coeff_dJp*((Area(K(rr,1)+R_row1,jj-1)/L(K(rr,1)+R_row1,jj-1))*((holeDC_av1*ddphi1+dhole_av1*dphi_DC1)+Thermal_V*ddhole1));
                Jp(T_step,1) = Jp(T_step,1) + coeff_dJp*((Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj))*((holeDC_av2*ddphi2+dhole_av2*dphi_DC2)+Thermal_V*ddhole2));
                Jd(T_step,1) = Jd(T_step,1) - esi*e0/deltat*((dphi(VR1+V1,1)-dphi(VR1+V2,1))-(dphi_old(VR1+V1,1)-dphi_old(VR1+V2,1)))*Area(K(rr,1)+R_row1,jj-1)/L(K(rr,1)+R_row1,jj-1);
                Jd(T_step,1) = Jd(T_step,1) - esi*e0/deltat*((dphi(VR1+V3,1)-dphi(VR1+V2,1))-(dphi_old(VR1+V3,1)-dphi_old(VR1+V2,1)))*Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj);

            elseif K(rr,1) > 2*R_row2 && K(rr,1) <= 3*R_row2

                jj = 3;
                K(rr,1) = K(rr,1)-2*R_row2;

                V3 = find(Table(1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj)) ; % node 3
                V2 = find(Table(1+1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj-1)); % node 2
                V1 = find(Table(2+1+VR1:3:VR1+3*VR2,2)==E_R2(K(rr,1),jj-2)); %node 1

                elecDC_av1 = 0.5*(elec_DC(VR1+V1,1)+elec_DC(VR1+V3,1));
                delec_av1 = 0.5*(delec(VR1+V1,1)+delec(VR1+V3,1));
                ddelec1 = delec(VR1+V1,1)-delec(VR1+V3,1);

                elecDC_av2 = 0.5*(elec_DC(VR1+V2,1)+elec_DC(VR1+V3,1));
                delec_av2 = 0.5*(delec(VR1+V2,1)+delec(VR1+V3,1));
                ddelec2 = delec(VR1+V2,1)-delec(VR1+V3,1);

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

                Jn(T_step,1) = Jn(T_step,1) + coeff_dJn*((Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj))*((elecDC_av1*ddphi1+delec_av1*dphi_DC1)-Thermal_V*ddelec1));
                Jn(T_step,1) = Jn(T_step,1) + coeff_dJn*((Area(K(rr,1)+R_row1,jj-1)/L(K(rr,1)+R_row1,jj-1))*((elecDC_av2*ddphi2+delec_av2*dphi_DC2)-Thermal_V*ddelec2));
                Jp(T_step,1) = Jp(T_step,1) + coeff_dJp*((Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj))*((holeDC_av1*ddphi1+dhole_av1*dphi_DC1)+Thermal_V*ddhole1));
                Jp(T_step,1) = Jp(T_step,1) + coeff_dJp*((Area(K(rr,1)+R_row1,jj-1)/L(K(rr,1)+R_row1,jj-1))*((holeDC_av2*ddphi2+dhole_av2*dphi_DC2)+Thermal_V*ddhole2));
                Jd(T_step,1) = Jd(T_step,1) - esi*e0/deltat*((dphi(VR1+V1,1)-dphi(VR1+V3,1))-(dphi_old(VR1+V1,1)-dphi_old(VR1+V3,1)))*Area(K(rr,1)+R_row1,jj)/L(K(rr,1)+R_row1,jj);
                Jd(T_step,1) = Jd(T_step,1) - esi*e0/deltat*((dphi(VR1+V2,1)-dphi(VR1+V3,1))-(dphi_old(VR1+V2,1)-dphi_old(VR1+V3,1)))*Area(K(rr,1)+R_row1,jj-1)/L(K(rr,1)+R_row1,jj-1);

            end
            I = (Jn+Jp+Jd)*1e-6;
        end
    end
    save_Vg(T_step,1) = Amp*sin(2*pi*freq*t(:,T_step));
end