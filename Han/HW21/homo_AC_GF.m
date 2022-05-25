clear
load("homo_ramping_information.mat")
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% homo AC Green Function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
freq = 1e12;
Vamp = 1e-3;
inum=sqrt(-1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ncoeff = q*mun;
pcoeff = q*mup;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dphi = zeros(size(phi));
dn = zeros(size(n));
dp = zeros(size(p));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
In_anode = zeros(1,1);
Ip_anode = zeros(1,1);
Id_anode = zeros(1,1);

In_cathode = zeros(1,1);
Ip_cathode = zeros(1,1);
Id_cathode = zeros(1,1);

Jaco = sparse(size(X2_sorted,1),size(X2_sorted,1));

saveexl = zeros(size(X1,1),1);

for i=1:Esirow
    if i <= E2row
        Jaco(X2(3*(X1_sorted(Esi(i,1),2))-2,1), X2(3*(X1_sorted(Esi(i,1),2))-2,1)) = Jaco(X2(3*(X1_sorted(Esi(i,1),2))-2,1), X2(3*(X1_sorted(Esi(i,1),2))-2,1))...
            - esi*(edgeEsi(i,1)/lenEsi(i,1) + edgeEsi(i,3)/lenEsi(i,3));
        Jaco(X2(3*(X1_sorted(Esi(i,1),2))-2,1), X2(3*(X1_sorted(Esi(i,2),2))-2,1)) = Jaco(X2(3*(X1_sorted(Esi(i,1),2))-2,1), X2(3*(X1_sorted(Esi(i,2),2))-2,1))...
            + esi*edgeEsi(i,1)/lenEsi(i,1);
        Jaco(X2(3*(X1_sorted(Esi(i,1),2))-2,1), X2(3*(X1_sorted(Esi(i,3),2))-2,1)) = Jaco(X2(3*(X1_sorted(Esi(i,1),2))-2,1), X2(3*(X1_sorted(Esi(i,3),2))-2,1))...
            + esi*edgeEsi(i,3)/lenEsi(i,3);
%         Res(X2(3*(X1_sorted(Esi(i,1),2))-2,1),1) =  Res(X2(3*(X1_sorted(Esi(i,1),2))-2,1),1)...
%             + esi*edgeEsi(i,1)/lenEsi(i,1)*(dphi(X1_sorted(Esi(i,2),2))-dphi(X1_sorted(Esi(i,1),2)))...
%             + esi*edgeEsi(i,3)/lenEsi(i,3)*(dphi(X1_sorted(Esi(i,3),2))-dphi(X1_sorted(Esi(i,1),2)));
        saveexl(X1_sorted(Esi(i,1),2)) = saveexl(X1_sorted(Esi(i,1),2)) + 0.25*(edgeEsi(i,1)*lenEsi(i,1) + edgeEsi(i,3)*lenEsi(i,3));

        Jaco(X2(3*(X1_sorted(Esi(i,2),2))-2,1), X2(3*(X1_sorted(Esi(i,2),2))-2,1)) = Jaco(X2(3*(X1_sorted(Esi(i,2),2))-2,1), X2(3*(X1_sorted(Esi(i,2),2))-2,1))...
            - esi*(edgeEsi(i,1)/lenEsi(i,1) + edgeEsi(i,2)/lenEsi(i,2));
        Jaco(X2(3*(X1_sorted(Esi(i,2),2))-2,1), X2(3*(X1_sorted(Esi(i,1),2))-2,1)) = Jaco(X2(3*(X1_sorted(Esi(i,2),2))-2,1), X2(3*(X1_sorted(Esi(i,1),2))-2,1))...
            + esi*edgeEsi(i,1)/lenEsi(i,1);
        Jaco(X2(3*(X1_sorted(Esi(i,2),2))-2,1), X2(3*(X1_sorted(Esi(i,3),2))-2,1)) = Jaco(X2(3*(X1_sorted(Esi(i,2),2))-2,1), X2(3*(X1_sorted(Esi(i,3),2))-2,1))...
            + esi*edgeEsi(i,2)/lenEsi(i,2);
%         Res(X2(3*(X1_sorted(Esi(i,2),2))-2,1),1) =  Res(X2(3*(X1_sorted(Esi(i,2),2))-2,1),1)...
%             + esi*edgeEsi(i,1)/lenEsi(i,1)*(dphi(X1_sorted(Esi(i,1),2))-dphi(X1_sorted(Esi(i,2),2)))...
%             + esi*edgeEsi(i,2)/lenEsi(i,2)*(dphi(X1_sorted(Esi(i,3),2))-dphi(X1_sorted(Esi(i,2),2)));
        saveexl(X1_sorted(Esi(i,2),2)) = saveexl(X1_sorted(Esi(i,2),2)) + 0.25*(edgeEsi(i,1)*lenEsi(i,1) + edgeEsi(i,2)*lenEsi(i,2));

        Jaco(X2(3*(X1_sorted(Esi(i,3),2))-2,1), X2(3*(X1_sorted(Esi(i,3),2))-2,1)) = Jaco(X2(3*(X1_sorted(Esi(i,3),2))-2,1), X2(3*(X1_sorted(Esi(i,3),2))-2,1))...
            - esi*(edgeEsi(i,2)/lenEsi(i,2) + edgeEsi(i,3)/lenEsi(i,3));
        Jaco(X2(3*(X1_sorted(Esi(i,3),2))-2,1), X2(3*(X1_sorted(Esi(i,2),2))-2,1)) = Jaco(X2(3*(X1_sorted(Esi(i,3),2))-2,1), X2(3*(X1_sorted(Esi(i,2),2))-2,1))...
            + esi*edgeEsi(i,2)/lenEsi(i,2);
        Jaco(X2(3*(X1_sorted(Esi(i,3),2))-2,1), X2(3*(X1_sorted(Esi(i,1),2))-2,1)) = Jaco(X2(3*(X1_sorted(Esi(i,3),2))-2,1), X2(3*(X1_sorted(Esi(i,1),2))-2,1))...
            + esi*edgeEsi(i,3)/lenEsi(i,3);
%         Res(X2(3*(X1_sorted(Esi(i,3),2))-2,1),1) =  Res(X2(3*(X1_sorted(Esi(i,3),2))-2,1),1)...
%             + esi*edgeEsi(i,2)/lenEsi(i,2)*(dphi(X1_sorted(Esi(i,2),2))-dphi(X1_sorted(Esi(i,3),2)))...
%             + esi*edgeEsi(i,3)/lenEsi(i,3)*(dphi(X1_sorted(Esi(i,1),2))-dphi(X1_sorted(Esi(i,3),2)));
        saveexl(X1_sorted(Esi(i,3),2)) = saveexl(X1_sorted(Esi(i,3),2)) + 0.25*(edgeEsi(i,2)*lenEsi(i,2) + edgeEsi(i,3)*lenEsi(i,3));

    elseif i > E2row && i <= (E2row + E3row)
        Jaco(X2(3*(X1_sorted(Esi(i,1),2))-2,1), X2(3*(X1_sorted(Esi(i,1),2))-2,1)) = Jaco(X2(3*(X1_sorted(Esi(i,1),2))-2,1), X2(3*(X1_sorted(Esi(i,1),2))-2,1))...
            - esi*(edgeEsi(i,1)/lenEsi(i,1) + edgeEsi(i,3)/lenEsi(i,3));
        Jaco(X2(3*(X1_sorted(Esi(i,1),2))-2,1), X2(3*(X1_sorted(Esi(i,2),2))-2,1)) = Jaco(X2(3*(X1_sorted(Esi(i,1),2))-2,1), X2(3*(X1_sorted(Esi(i,2),2))-2,1))...
            + esi*edgeEsi(i,1)/lenEsi(i,1);
        Jaco(X2(3*(X1_sorted(Esi(i,1),2))-2,1), X2(3*(X1_sorted(Esi(i,3),2))-2,1)) = Jaco(X2(3*(X1_sorted(Esi(i,1),2))-2,1), X2(3*(X1_sorted(Esi(i,3),2))-2,1))...
            + esi*edgeEsi(i,3)/lenEsi(i,3);
%         Res(X2(3*(X1_sorted(Esi(i,1),2))-2,1),1) =  Res(X2(3*(X1_sorted(Esi(i,1),2))-2,1),1)...
%             + esi*edgeEsi(i,1)/lenEsi(i,1)*(dphi(X1_sorted(Esi(i,2),2))-dphi(X1_sorted(Esi(i,1),2)))...
%             + esi*edgeEsi(i,3)/lenEsi(i,3)*(dphi(X1_sorted(Esi(i,3),2))-dphi(X1_sorted(Esi(i,1),2)));
        saveexl(X1_sorted(Esi(i,1),2)) = saveexl(X1_sorted(Esi(i,1),2)) + 0.25*(edgeEsi(i,1)*lenEsi(i,1) + edgeEsi(i,3)*lenEsi(i,3));

        Jaco(X2(3*(X1_sorted(Esi(i,2),2))-2,1), X2(3*(X1_sorted(Esi(i,2),2))-2,1)) = Jaco(X2(3*(X1_sorted(Esi(i,2),2))-2,1), X2(3*(X1_sorted(Esi(i,2),2))-2,1))...
            - esi*(edgeEsi(i,1)/lenEsi(i,1) + edgeEsi(i,2)/lenEsi(i,2));
        Jaco(X2(3*(X1_sorted(Esi(i,2),2))-2,1), X2(3*(X1_sorted(Esi(i,1),2))-2,1)) = Jaco(X2(3*(X1_sorted(Esi(i,2),2))-2,1), X2(3*(X1_sorted(Esi(i,1),2))-2,1))...
            + esi*edgeEsi(i,1)/lenEsi(i,1);
        Jaco(X2(3*(X1_sorted(Esi(i,2),2))-2,1), X2(3*(X1_sorted(Esi(i,3),2))-2,1)) = Jaco(X2(3*(X1_sorted(Esi(i,2),2))-2,1), X2(3*(X1_sorted(Esi(i,3),2))-2,1))...
            + esi*edgeEsi(i,2)/lenEsi(i,2);
%         Res(X2(3*(X1_sorted(Esi(i,2),2))-2,1),1) =  Res(X2(3*(X1_sorted(Esi(i,2),2))-2,1),1)...
%             + esi*edgeEsi(i,1)/lenEsi(i,1)*(dphi(X1_sorted(Esi(i,1),2))-dphi(X1_sorted(Esi(i,2),2)))...
%             + esi*edgeEsi(i,2)/lenEsi(i,2)*(dphi(X1_sorted(Esi(i,3),2))-dphi(X1_sorted(Esi(i,2),2)));
        saveexl(X1_sorted(Esi(i,2),2)) = saveexl(X1_sorted(Esi(i,2),2)) + 0.25*(edgeEsi(i,1)*lenEsi(i,1) + edgeEsi(i,2)*lenEsi(i,2));

        Jaco(X2(3*(X1_sorted(Esi(i,3),2))-2,1), X2(3*(X1_sorted(Esi(i,3),2))-2,1)) = Jaco(X2(3*(X1_sorted(Esi(i,3),2))-2,1), X2(3*(X1_sorted(Esi(i,3),2))-2,1))...
            - esi*(edgeEsi(i,2)/lenEsi(i,2) + edgeEsi(i,3)/lenEsi(i,3));
        Jaco(X2(3*(X1_sorted(Esi(i,3),2))-2,1), X2(3*(X1_sorted(Esi(i,2),2))-2,1)) = Jaco(X2(3*(X1_sorted(Esi(i,3),2))-2,1), X2(3*(X1_sorted(Esi(i,2),2))-2,1))...
            + esi*edgeEsi(i,2)/lenEsi(i,2);
        Jaco(X2(3*(X1_sorted(Esi(i,3),2))-2,1), X2(3*(X1_sorted(Esi(i,1),2))-2,1)) = Jaco(X2(3*(X1_sorted(Esi(i,3),2))-2,1), X2(3*(X1_sorted(Esi(i,1),2))-2,1))...
            + esi*edgeEsi(i,3)/lenEsi(i,3);
%         Res(X2(3*(X1_sorted(Esi(i,3),2))-2,1),1) =  Res(X2(3*(X1_sorted(Esi(i,3),2))-2,1),1)...
%             + esi*edgeEsi(i,2)/lenEsi(i,2)*(dphi(X1_sorted(Esi(i,2),2))-dphi(X1_sorted(Esi(i,3),2)))...
%             + esi*edgeEsi(i,3)/lenEsi(i,3)*(dphi(X1_sorted(Esi(i,1),2))-dphi(X1_sorted(Esi(i,3),2)));
        saveexl(X1_sorted(Esi(i,3),2)) = saveexl(X1_sorted(Esi(i,3),2)) + 0.25*(edgeEsi(i,2)*lenEsi(i,2) + edgeEsi(i,3)*lenEsi(i,3));

    elseif i > (E2row+E3row) && i <= Esirow
        Jaco(X2(3*(X1_sorted(Esi(i,1),2))-2,1), X2(3*(X1_sorted(Esi(i,1),2))-2,1)) = Jaco(X2(3*(X1_sorted(Esi(i,1),2))-2,1), X2(3*(X1_sorted(Esi(i,1),2))-2,1))...
            - esi*(edgeEsi(i,1)/lenEsi(i,1) + edgeEsi(i,3)/lenEsi(i,3));
        Jaco(X2(3*(X1_sorted(Esi(i,1),2))-2,1), X2(3*(X1_sorted(Esi(i,2),2))-2,1)) = Jaco(X2(3*(X1_sorted(Esi(i,1),2))-2,1), X2(3*(X1_sorted(Esi(i,2),2))-2,1))...
            + esi*edgeEsi(i,1)/lenEsi(i,1);
        Jaco(X2(3*(X1_sorted(Esi(i,1),2))-2,1), X2(3*(X1_sorted(Esi(i,3),2))-2,1)) = Jaco(X2(3*(X1_sorted(Esi(i,1),2))-2,1), X2(3*(X1_sorted(Esi(i,3),2))-2,1))...
            + esi*edgeEsi(i,3)/lenEsi(i,3);
%         Res(X2(3*(X1_sorted(Esi(i,1),2))-2,1),1) =  Res(X2(3*(X1_sorted(Esi(i,1),2))-2,1),1)...
%             + esi*edgeEsi(i,1)/lenEsi(i,1)*(dphi(X1_sorted(Esi(i,2),2))-dphi(X1_sorted(Esi(i,1),2)))...
%             + esi*edgeEsi(i,3)/lenEsi(i,3)*(dphi(X1_sorted(Esi(i,3),2))-dphi(X1_sorted(Esi(i,1),2)));
        saveexl(X1_sorted(Esi(i,1),2)) = saveexl(X1_sorted(Esi(i,1),2)) + 0.25*(edgeEsi(i,1)*lenEsi(i,1) + edgeEsi(i,3)*lenEsi(i,3));

        Jaco(X2(3*(X1_sorted(Esi(i,2),2))-2,1), X2(3*(X1_sorted(Esi(i,2),2))-2,1)) = Jaco(X2(3*(X1_sorted(Esi(i,2),2))-2,1), X2(3*(X1_sorted(Esi(i,2),2))-2,1))...
            - esi*(edgeEsi(i,1)/lenEsi(i,1) + edgeEsi(i,2)/lenEsi(i,2));
        Jaco(X2(3*(X1_sorted(Esi(i,2),2))-2,1), X2(3*(X1_sorted(Esi(i,1),2))-2,1)) = Jaco(X2(3*(X1_sorted(Esi(i,2),2))-2,1), X2(3*(X1_sorted(Esi(i,1),2))-2,1))...
            + esi*edgeEsi(i,1)/lenEsi(i,1);
        Jaco(X2(3*(X1_sorted(Esi(i,2),2))-2,1), X2(3*(X1_sorted(Esi(i,3),2))-2,1)) = Jaco(X2(3*(X1_sorted(Esi(i,2),2))-2,1), X2(3*(X1_sorted(Esi(i,3),2))-2,1))...
            + esi*edgeEsi(i,2)/lenEsi(i,2);
%         Res(X2(3*(X1_sorted(Esi(i,2),2))-2,1),1) =  Res(X2(3*(X1_sorted(Esi(i,2),2))-2,1),1)...
%             + esi*edgeEsi(i,1)/lenEsi(i,1)*(dphi(X1_sorted(Esi(i,1),2))-dphi(X1_sorted(Esi(i,2),2)))...
%             + esi*edgeEsi(i,2)/lenEsi(i,2)*(dphi(X1_sorted(Esi(i,3),2))-dphi(X1_sorted(Esi(i,2),2)));
        saveexl(X1_sorted(Esi(i,2),2)) = saveexl(X1_sorted(Esi(i,2),2)) + 0.25*(edgeEsi(i,1)*lenEsi(i,1) + edgeEsi(i,2)*lenEsi(i,2));

        Jaco(X2(3*(X1_sorted(Esi(i,3),2))-2,1), X2(3*(X1_sorted(Esi(i,3),2))-2,1)) = Jaco(X2(3*(X1_sorted(Esi(i,3),2))-2,1), X2(3*(X1_sorted(Esi(i,3),2))-2,1))...
            - esi*(edgeEsi(i,2)/lenEsi(i,2) + edgeEsi(i,3)/lenEsi(i,3));
        Jaco(X2(3*(X1_sorted(Esi(i,3),2))-2,1), X2(3*(X1_sorted(Esi(i,2),2))-2,1)) = Jaco(X2(3*(X1_sorted(Esi(i,3),2))-2,1), X2(3*(X1_sorted(Esi(i,2),2))-2,1))...
            + esi*edgeEsi(i,2)/lenEsi(i,2);
        Jaco(X2(3*(X1_sorted(Esi(i,3),2))-2,1), X2(3*(X1_sorted(Esi(i,1),2))-2,1)) = Jaco(X2(3*(X1_sorted(Esi(i,3),2))-2,1), X2(3*(X1_sorted(Esi(i,1),2))-2,1))...
            + esi*edgeEsi(i,3)/lenEsi(i,3);
%         Res(X2(3*(X1_sorted(Esi(i,3),2))-2,1),1) =  Res(X2(3*(X1_sorted(Esi(i,3),2))-2,1),1)...
%             + esi*edgeEsi(i,2)/lenEsi(i,2)*(dphi(X1_sorted(Esi(i,2),2))-dphi(X1_sorted(Esi(i,3),2)))...
%             + esi*edgeEsi(i,3)/lenEsi(i,3)*(dphi(X1_sorted(Esi(i,1),2))-dphi(X1_sorted(Esi(i,3),2)));
        saveexl(X1_sorted(Esi(i,3),2)) = saveexl(X1_sorted(Esi(i,3),2)) + 0.25*(edgeEsi(i,2)*lenEsi(i,2) + edgeEsi(i,3)*lenEsi(i,3));

    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Possion charge
for i=1:V_silicon_row
    % - qn
    Jaco(X2(3*V_silicon(i,2)-2),X2(3*V_silicon(i,2)-1)) = Jaco(X2(3*V_silicon(i,2)-2),X2(3*V_silicon(i,2)-1)) - saveexl(V_silicon(i,2))*q/e0;
%     Res(X2(3*V_silicon(i,2)-2),1) = Res(X2(3*V_silicon(i,2)-2),1) - saveexl(V_silicon(i,2))*q/e0*dn(i);

    % + qp
    Jaco(X2(3*V_silicon(i,2)-2),X2(3*V_silicon(i,2))) =  Jaco(X2(3*V_silicon(i,2)-2),X2(3*V_silicon(i,2))) + saveexl(V_silicon(i,2))*q/e0;
%     Res(X2(3*V_silicon(i,2)-2),1) = Res(X2(3*V_silicon(i,2)-2),1) + saveexl(V_silicon(i,2))*q/e0*dp(i);
end

%%%  electron DD
for i=1:Esirow

    f1 = find(V_silicon(:,1)==Esi(i,1));
    f2 = find(V_silicon(:,1)==Esi(i,2));
    f3 = find(V_silicon(:,1)==Esi(i,3));

    for j = 1:3

        if j == 1 %%% 1st element

            n_DC_av1 = 0.5*(n_DC(f2)+n_DC(f1));
            n_DC_av2 = 0.5*(n_DC(f3)+n_DC(f1));
            dn_av1 = 0.5*(dn(f2)+dn(f1));
            dn_av2 = 0.5*(dn(f3)+dn(f1));

            ddphi1 = dphi(X1_sorted(Esi(i,2),2))-dphi(X1_sorted(Esi(i,1),2));
            ddphi2 = dphi(X1_sorted(Esi(i,3),2))-dphi(X1_sorted(Esi(i,1),2));
            dphi_DC1 = phi_DC(X1_sorted(Esi(i,2),2))-phi_DC(X1_sorted(Esi(i,1),2));
            dphi_DC2 = phi_DC(X1_sorted(Esi(i,3),2))-phi_DC(X1_sorted(Esi(i,1),2));

            ddn1 = dn(f2)-dn(f1);
            ddn2 = dn(f3)-dn(f1);

            % n part
            Jaco(X2(3*X1_sorted(Esi(i,1),2)-1,1), X2(3*X1_sorted(Esi(i,1),2)-1,1)) = Jaco(X2(3*X1_sorted(Esi(i,1),2)-1,1), X2(3*X1_sorted(Esi(i,1),2)-1,1))...
                - ncoeff*edgeEsi(i,1)/lenEsi(i,1)*(0.5*dphi_DC1 + Vt)...
                - ncoeff*edgeEsi(i,3)/lenEsi(i,3)*(0.5*dphi_DC2 + Vt);

            Jaco(X2(3*X1_sorted(Esi(i,1),2)-1,1), X2(3*X1_sorted(Esi(i,2),2)-1,1)) = Jaco(X2(3*X1_sorted(Esi(i,1),2)-1,1), X2(3*X1_sorted(Esi(i,2),2)-1,1))...
                - ncoeff*edgeEsi(i,1)/lenEsi(i,1)*(0.5*dphi_DC1 - Vt);
            Jaco(X2(3*X1_sorted(Esi(i,1),2)-1,1), X2(3*X1_sorted(Esi(i,3),2)-1,1)) = Jaco(X2(3*X1_sorted(Esi(i,1),2)-1,1), X2(3*X1_sorted(Esi(i,3),2)-1,1))...
                - ncoeff*edgeEsi(i,3)/lenEsi(i,3)*(0.5*dphi_DC2 - Vt);

            % phi part
            Jaco(X2(3*X1_sorted(Esi(i,1),2)-1,1), X2(3*X1_sorted(Esi(i,1),2)-2,1)) = Jaco(X2(3*X1_sorted(Esi(i,1),2)-1,1), X2(3*X1_sorted(Esi(i,1),2)-2,1))...
                - ncoeff*edgeEsi(i,1)/lenEsi(i,1)*(-n_DC_av1)...
                - ncoeff*edgeEsi(i,3)/lenEsi(i,3)*(-n_DC_av2);

            Jaco(X2(3*X1_sorted(Esi(i,1),2)-1,1), X2(3*X1_sorted(Esi(i,2),2)-2,1)) = Jaco(X2(3*X1_sorted(Esi(i,1),2)-1,1), X2(3*X1_sorted(Esi(i,2),2)-2,1))...
                - ncoeff*edgeEsi(i,1)/lenEsi(i,1)*n_DC_av1;
            Jaco(X2(3*X1_sorted(Esi(i,1),2)-1,1), X2(3*X1_sorted(Esi(i,3),2)-2,1)) = Jaco(X2(3*X1_sorted(Esi(i,1),2)-1,1), X2(3*X1_sorted(Esi(i,3),2)-2,1))...
                - ncoeff*edgeEsi(i,3)/lenEsi(i,3)*n_DC_av2;

            % n AC part
            Jaco(X2(3*X1_sorted(Esi(i,1),2)-1,1), X2(3*X1_sorted(Esi(i,1),2)-1,1)) = Jaco(X2(3*X1_sorted(Esi(i,1),2)-1,1), X2(3*X1_sorted(Esi(i,1),2)-1,1))...
                - (2*pi*freq*inum)*0.25*q*(lenEsi(i,1)*edgeEsi(i,1) + lenEsi(i,3)*edgeEsi(i,3));

        elseif j == 2 %%% 2nd element

            n_DC_av1 = 0.5*(n_DC(f3)+n_DC(f2));
            n_DC_av2 = 0.5*(n_DC(f1)+n_DC(f2));
            dn_av1 = 0.5*(dn(f3)+dn(f2));
            dn_av2 = 0.5*(dn(f1)+dn(f2));

            ddphi1 = dphi(X1_sorted(Esi(i,3),2))-dphi(X1_sorted(Esi(i,2),2));
            ddphi2 = dphi(X1_sorted(Esi(i,1),2))-dphi(X1_sorted(Esi(i,2),2));
            dphi_DC1 = phi_DC(X1_sorted(Esi(i,3),2))-phi_DC(X1_sorted(Esi(i,2),2));
            dphi_DC2 = phi_DC(X1_sorted(Esi(i,1),2))-phi_DC(X1_sorted(Esi(i,2),2));

            ddn1 = dn(f3)-dn(f2);
            ddn2 = dn(f1)-dn(f2);

            % n part
            Jaco(X2(3*X1_sorted(Esi(i,2),2)-1,1), X2(3*X1_sorted(Esi(i,2),2)-1,1)) = Jaco(X2(3*X1_sorted(Esi(i,2),2)-1,1), X2(3*X1_sorted(Esi(i,2),2)-1,1))...
                - ncoeff*edgeEsi(i,2)/lenEsi(i,2)*(0.5*dphi_DC1 + Vt)...
                - ncoeff*edgeEsi(i,1)/lenEsi(i,1)*(0.5*dphi_DC2 + Vt);
            Jaco(X2(3*X1_sorted(Esi(i,2),2)-1,1), X2(3*X1_sorted(Esi(i,3),2)-1,1)) = Jaco(X2(3*X1_sorted(Esi(i,2),2)-1,1), X2(3*X1_sorted(Esi(i,3),2)-1,1))...
                - ncoeff*edgeEsi(i,2)/lenEsi(i,2)*(0.5*dphi_DC1 - Vt);
            Jaco(X2(3*X1_sorted(Esi(i,2),2)-1,1), X2(3*X1_sorted(Esi(i,1),2)-1,1)) = Jaco(X2(3*X1_sorted(Esi(i,2),2)-1,1), X2(3*X1_sorted(Esi(i,1),2)-1,1))...
                - ncoeff*edgeEsi(i,1)/lenEsi(i,1)*(0.5*dphi_DC2 - Vt);

            % phi part
            Jaco(X2(3*X1_sorted(Esi(i,2),2)-1,1), X2(3*X1_sorted(Esi(i,2),2)-2,1)) = Jaco(X2(3*X1_sorted(Esi(i,2),2)-1,1), X2(3*X1_sorted(Esi(i,2),2)-2,1))...
                - ncoeff*edgeEsi(i,2)/lenEsi(i,2)*(-n_DC_av1)...
                - ncoeff*edgeEsi(i,1)/lenEsi(i,1)*(-n_DC_av2);
            Jaco(X2(3*X1_sorted(Esi(i,2),2)-1,1), X2(3*X1_sorted(Esi(i,3),2)-2,1)) = Jaco(X2(3*X1_sorted(Esi(i,2),2)-1,1), X2(3*X1_sorted(Esi(i,3),2)-2,1))...
                - ncoeff*edgeEsi(i,2)/lenEsi(i,2)*n_DC_av1;
            Jaco(X2(3*X1_sorted(Esi(i,2),2)-1,1), X2(3*X1_sorted(Esi(i,1),2)-2,1)) = Jaco(X2(3*X1_sorted(Esi(i,2),2)-1,1), X2(3*X1_sorted(Esi(i,1),2)-2,1))...
                - ncoeff*edgeEsi(i,1)/lenEsi(i,1)*n_DC_av2;

            % n AC part
            Jaco(X2(3*X1_sorted(Esi(i,2),2)-1,1), X2(3*X1_sorted(Esi(i,2),2)-1,1)) = Jaco(X2(3*X1_sorted(Esi(i,2),2)-1,1), X2(3*X1_sorted(Esi(i,2),2)-1,1))...
                - (2*pi*freq*inum)*0.25*q*(lenEsi(i,2)*edgeEsi(i,2) + lenEsi(i,1)*edgeEsi(i,1));

        elseif j == 3  %%% 3rd element

            n_DC_av1 = 0.5*(n_DC(f1)+n_DC(f3));
            n_DC_av2 = 0.5*(n_DC(f2)+n_DC(f3));
            dn_av1 = 0.5*(dn(f1)+dn(f3));
            dn_av2 = 0.5*(dn(f2)+dn(f3));

            ddphi1 = dphi(X1_sorted(Esi(i,1),2))-dphi(X1_sorted(Esi(i,3),2));
            ddphi2 = dphi(X1_sorted(Esi(i,2),2))-dphi(X1_sorted(Esi(i,3),2));
            dphi_DC1 = phi_DC(X1_sorted(Esi(i,1),2))-phi_DC(X1_sorted(Esi(i,3),2));
            dphi_DC2 = phi_DC(X1_sorted(Esi(i,2),2))-phi_DC(X1_sorted(Esi(i,3),2));

            ddn1 = dn(f1)-dn(f3);
            ddn2 = dn(f2)-dn(f3);

            % n part
            Jaco(X2(3*X1_sorted(Esi(i,3),2)-1,1), X2(3*X1_sorted(Esi(i,3),2)-1,1)) = Jaco(X2(3*X1_sorted(Esi(i,3),2)-1,1), X2(3*X1_sorted(Esi(i,3),2)-1,1))...
                - ncoeff*edgeEsi(i,3)/lenEsi(i,3)*(0.5*dphi_DC1 + Vt)...
                - ncoeff*edgeEsi(i,2)/lenEsi(i,2)*(0.5*dphi_DC2 + Vt);
            Jaco(X2(3*X1_sorted(Esi(i,3),2)-1,1), X2(3*X1_sorted(Esi(i,1),2)-1,1)) = Jaco(X2(3*X1_sorted(Esi(i,3),2)-1,1), X2(3*X1_sorted(Esi(i,1),2)-1,1))...
                - ncoeff*edgeEsi(i,3)/lenEsi(i,3)*(0.5*dphi_DC1 - Vt);
            Jaco(X2(3*X1_sorted(Esi(i,3),2)-1,1), X2(3*X1_sorted(Esi(i,2),2)-1,1)) = Jaco(X2(3*X1_sorted(Esi(i,3),2)-1,1), X2(3*X1_sorted(Esi(i,2),2)-1,1))...
                - ncoeff*edgeEsi(i,2)/lenEsi(i,2)*(0.5*dphi_DC2 - Vt);

            % phi part
            Jaco(X2(3*X1_sorted(Esi(i,3),2)-1,1), X2(3*X1_sorted(Esi(i,3),2)-2,1)) = Jaco(X2(3*X1_sorted(Esi(i,3),2)-1,1), X2(3*X1_sorted(Esi(i,3),2)-2,1))...
                - ncoeff*edgeEsi(i,3)/lenEsi(i,3)*(-n_DC_av1)...
                - ncoeff*edgeEsi(i,2)/lenEsi(i,2)*(-n_DC_av2);
            Jaco(X2(3*X1_sorted(Esi(i,3),2)-1,1), X2(3*X1_sorted(Esi(i,1),2)-2,1)) = Jaco(X2(3*X1_sorted(Esi(i,3),2)-1,1), X2(3*X1_sorted(Esi(i,1),2)-2,1))...
                - ncoeff*edgeEsi(i,3)/lenEsi(i,3)*n_DC_av1;
            Jaco(X2(3*X1_sorted(Esi(i,3),2)-1,1), X2(3*X1_sorted(Esi(i,2),2)-2,1)) = Jaco(X2(3*X1_sorted(Esi(i,3),2)-1,1), X2(3*X1_sorted(Esi(i,2),2)-2,1))...
                - ncoeff*edgeEsi(i,2)/lenEsi(i,2)*n_DC_av2;

            % n AC part
            Jaco(X2(3*X1_sorted(Esi(i,3),2)-1,1), X2(3*X1_sorted(Esi(i,3),2)-1,1)) = Jaco(X2(3*X1_sorted(Esi(i,3),2)-1,1), X2(3*X1_sorted(Esi(i,3),2)-1,1))...
                - (2*pi*freq*inum)*0.25*q*(lenEsi(i,3)*edgeEsi(i,3) + lenEsi(i,2)*edgeEsi(i,2));
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% hole DD
for i=1:Esirow

    f1 = find(V_silicon(:,1)==Esi(i,1));
    f2 = find(V_silicon(:,1)==Esi(i,2));
    f3 = find(V_silicon(:,1)==Esi(i,3));

    for j = 1:3

        if j == 1 %%% 1st element

            p_DC_av1 = 0.5*(p_DC(f2)+p_DC(f1));
            p_DC_av2 = 0.5*(p_DC(f3)+p_DC(f1));
            dp_av1 = 0.5*(dp(f2)+dp(f1));
            dp_av2 = 0.5*(dp(f3)+dp(f1));

            ddphi1 = dphi(X1_sorted(Esi(i,2),2))-dphi(X1_sorted(Esi(i,1),2));
            ddphi2 = dphi(X1_sorted(Esi(i,3),2))-dphi(X1_sorted(Esi(i,1),2));
            dphi_DC1 = phi_DC(X1_sorted(Esi(i,2),2))-phi_DC(X1_sorted(Esi(i,1),2));
            dphi_DC2 = phi_DC(X1_sorted(Esi(i,3),2))-phi_DC(X1_sorted(Esi(i,1),2));

            ddp1 = dp(f2)-dp(f1);
            ddp2 = dp(f3)-dp(f1);

            % p part
            Jaco(X2(3*X1_sorted(Esi(i,1),2),1), X2(3*X1_sorted(Esi(i,1),2),1)) = Jaco(X2(3*X1_sorted(Esi(i,1),2),1), X2(3*X1_sorted(Esi(i,1),2),1))...
                - pcoeff*edgeEsi(i,1)/lenEsi(i,1)*(0.5*dphi_DC1 - Vt)...
                - pcoeff*edgeEsi(i,3)/lenEsi(i,3)*(0.5*dphi_DC2 - Vt);
            Jaco(X2(3*X1_sorted(Esi(i,1),2),1), X2(3*X1_sorted(Esi(i,2),2),1)) = Jaco(X2(3*X1_sorted(Esi(i,1),2),1), X2(3*X1_sorted(Esi(i,2),2),1))...
                - pcoeff*edgeEsi(i,1)/lenEsi(i,1)*(0.5*dphi_DC1 + Vt);
            Jaco(X2(3*X1_sorted(Esi(i,1),2),1), X2(3*X1_sorted(Esi(i,3),2),1)) = Jaco(X2(3*X1_sorted(Esi(i,1),2),1), X2(3*X1_sorted(Esi(i,3),2),1))...
                - pcoeff*edgeEsi(i,3)/lenEsi(i,3)*(0.5*dphi_DC2 + Vt);

            % phi part
            Jaco(X2(3*X1_sorted(Esi(i,1),2),1), X2(3*X1_sorted(Esi(i,1),2)-2,1)) = Jaco(X2(3*X1_sorted(Esi(i,1),2),1), X2(3*X1_sorted(Esi(i,1),2)-2,1))...
                - pcoeff*edgeEsi(i,1)/lenEsi(i,1)*(-p_DC_av1)...
                - pcoeff*edgeEsi(i,3)/lenEsi(i,3)*(-p_DC_av2);
            Jaco(X2(3*X1_sorted(Esi(i,1),2),1), X2(3*X1_sorted(Esi(i,2),2)-2,1)) = Jaco(X2(3*X1_sorted(Esi(i,1),2),1), X2(3*X1_sorted(Esi(i,2),2)-2,1))...
                - pcoeff*edgeEsi(i,1)/lenEsi(i,1)*p_DC_av1;
            Jaco(X2(3*X1_sorted(Esi(i,1),2),1), X2(3*X1_sorted(Esi(i,3),2)-2,1)) = Jaco(X2(3*X1_sorted(Esi(i,1),2),1), X2(3*X1_sorted(Esi(i,3),2)-2,1))...
                - pcoeff*edgeEsi(i,3)/lenEsi(i,3)*p_DC_av2;

            % p AC part
            Jaco(X2(3*X1_sorted(Esi(i,1),2),1), X2(3*X1_sorted(Esi(i,1),2),1)) = Jaco(X2(3*X1_sorted(Esi(i,1),2),1), X2(3*X1_sorted(Esi(i,1),2),1))...
                + (2*pi*freq*inum)*0.25*q*(lenEsi(i,1)*edgeEsi(i,1) + lenEsi(i,3)*edgeEsi(i,3));

        elseif j == 2 %%% 2nd element

            p_DC_av1 = 0.5*(p_DC(f3)+p_DC(f2));
            p_DC_av2 = 0.5*(p_DC(f1)+p_DC(f2));
            dp_av1 = 0.5*(dp(f3)+dp(f2));
            dp_av2 = 0.5*(dp(f1)+dp(f2));

            ddphi1 = dphi(X1_sorted(Esi(i,3),2))-dphi(X1_sorted(Esi(i,2),2));
            ddphi2 = dphi(X1_sorted(Esi(i,1),2))-dphi(X1_sorted(Esi(i,2),2));
            dphi_DC1 = phi_DC(X1_sorted(Esi(i,3),2))-phi_DC(X1_sorted(Esi(i,2),2));
            dphi_DC2 = phi_DC(X1_sorted(Esi(i,1),2))-phi_DC(X1_sorted(Esi(i,2),2));

            ddp1 = dp(f3)-dp(f2);
            ddp2 = dp(f1)-dp(f2);

            % p part
            Jaco(X2(3*X1_sorted(Esi(i,2),2),1), X2(3*X1_sorted(Esi(i,2),2),1)) = Jaco(X2(3*X1_sorted(Esi(i,2),2),1), X2(3*X1_sorted(Esi(i,2),2),1))...
                - pcoeff*edgeEsi(i,2)/lenEsi(i,2)*(0.5*dphi_DC1 - Vt)...
                - pcoeff*edgeEsi(i,1)/lenEsi(i,1)*(0.5*dphi_DC2 - Vt);
            Jaco(X2(3*X1_sorted(Esi(i,2),2),1), X2(3*X1_sorted(Esi(i,3),2),1)) = Jaco(X2(3*X1_sorted(Esi(i,2),2),1), X2(3*X1_sorted(Esi(i,3),2),1))...
                - pcoeff*edgeEsi(i,2)/lenEsi(i,2)*(0.5*dphi_DC1 + Vt);
            Jaco(X2(3*X1_sorted(Esi(i,2),2),1), X2(3*X1_sorted(Esi(i,1),2),1)) = Jaco(X2(3*X1_sorted(Esi(i,2),2),1), X2(3*X1_sorted(Esi(i,1),2),1))...
                - pcoeff*edgeEsi(i,1)/lenEsi(i,1)*(0.5*dphi_DC2 + Vt);

            % phi part
            Jaco(X2(3*X1_sorted(Esi(i,2),2),1), X2(3*X1_sorted(Esi(i,2),2)-2,1)) = Jaco(X2(3*X1_sorted(Esi(i,2),2),1), X2(3*X1_sorted(Esi(i,2),2)-2,1))...
                - pcoeff*edgeEsi(i,2)/lenEsi(i,2)*(-p_DC_av1)...
                - pcoeff*edgeEsi(i,1)/lenEsi(i,1)*(-p_DC_av2);
            Jaco(X2(3*X1_sorted(Esi(i,2),2),1), X2(3*X1_sorted(Esi(i,3),2)-2,1)) = Jaco(X2(3*X1_sorted(Esi(i,2),2),1), X2(3*X1_sorted(Esi(i,3),2)-2,1))...
                - pcoeff*edgeEsi(i,2)/lenEsi(i,2)*p_DC_av1;
            Jaco(X2(3*X1_sorted(Esi(i,2),2),1), X2(3*X1_sorted(Esi(i,1),2)-2,1)) = Jaco(X2(3*X1_sorted(Esi(i,2),2),1), X2(3*X1_sorted(Esi(i,1),2)-2,1))...
                - pcoeff*edgeEsi(i,1)/lenEsi(i,1)*p_DC_av2;

            % p AC part
            Jaco(X2(3*X1_sorted(Esi(i,2),2),1), X2(3*X1_sorted(Esi(i,2),2),1)) = Jaco(X2(3*X1_sorted(Esi(i,2),2),1), X2(3*X1_sorted(Esi(i,2),2),1))...
                + (2*pi*freq*inum)*0.25*q*(lenEsi(i,2)*edgeEsi(i,2) + lenEsi(i,1)*edgeEsi(i,1));

        elseif j == 3 %%% 3rd element

            p_DC_av1 = (p_DC(f1)+p_DC(f3))/2;
            p_DC_av2 = (p_DC(f2)+p_DC(f3))/2;
            dp_av1 = (dp(f1)+dp(f3))/2;
            dp_av2 = (dp(f2)+dp(f3))/2;

            ddphi1 = dphi(X1_sorted(Esi(i,1),2))-dphi(X1_sorted(Esi(i,3),2));
            ddphi2 = dphi(X1_sorted(Esi(i,2),2))-dphi(X1_sorted(Esi(i,3),2));
            dphi_DC1 = phi_DC(X1_sorted(Esi(i,1),2))-phi_DC(X1_sorted(Esi(i,3),2));
            dphi_DC2 = phi_DC(X1_sorted(Esi(i,2),2))-phi_DC(X1_sorted(Esi(i,3),2));

            ddp1 = dp(f1)-dp(f3);
            ddp2 = dp(f2)-dp(f3);

            % p part
            Jaco(X2(3*X1_sorted(Esi(i,3),2),1), X2(3*X1_sorted(Esi(i,3),2),1)) = Jaco(X2(3*X1_sorted(Esi(i,3),2),1), X2(3*X1_sorted(Esi(i,3),2),1))...
                - pcoeff*edgeEsi(i,3)/lenEsi(i,3)*(0.5*dphi_DC1 - Vt)...
                - pcoeff*edgeEsi(i,2)/lenEsi(i,2)*(0.5*dphi_DC2 - Vt);
            Jaco(X2(3*X1_sorted(Esi(i,3),2),1), X2(3*X1_sorted(Esi(i,1),2),1)) = Jaco(X2(3*X1_sorted(Esi(i,3),2),1), X2(3*X1_sorted(Esi(i,1),2),1))...
                - pcoeff*edgeEsi(i,3)/lenEsi(i,3)*(0.5*dphi_DC1 + Vt);
            Jaco(X2(3*X1_sorted(Esi(i,3),2),1), X2(3*X1_sorted(Esi(i,2),2),1)) = Jaco(X2(3*X1_sorted(Esi(i,3),2),1), X2(3*X1_sorted(Esi(i,2),2),1))...
                - pcoeff*edgeEsi(i,2)/lenEsi(i,2)*(0.5*dphi_DC2 + Vt);

            % phi part
            Jaco(X2(3*X1_sorted(Esi(i,3),2),1), X2(3*X1_sorted(Esi(i,3),2)-2,1)) = Jaco(X2(3*X1_sorted(Esi(i,3),2),1), X2(3*X1_sorted(Esi(i,3),2)-2,1))...
                - pcoeff*edgeEsi(i,3)/lenEsi(i,3)*(-p_DC_av1)...
                - pcoeff*edgeEsi(i,2)/lenEsi(i,2)*(-p_DC_av2);
            Jaco(X2(3*X1_sorted(Esi(i,3),2),1), X2(3*X1_sorted(Esi(i,1),2)-2,1)) = Jaco(X2(3*X1_sorted(Esi(i,3),2),1), X2(3*X1_sorted(Esi(i,1),2)-2,1))...
                - pcoeff*edgeEsi(i,3)/lenEsi(i,3)*p_DC_av1;
            Jaco(X2(3*X1_sorted(Esi(i,3),2),1), X2(3*X1_sorted(Esi(i,2),2)-2,1)) = Jaco(X2(3*X1_sorted(Esi(i,3),2),1), X2(3*X1_sorted(Esi(i,2),2)-2,1))...
                - pcoeff*edgeEsi(i,2)/lenEsi(i,2)*p_DC_av2;

            % p AC part
            Jaco(X2(3*X1_sorted(Esi(i,3),2),1), X2(3*X1_sorted(Esi(i,3),2),1)) = Jaco(X2(3*X1_sorted(Esi(i,3),2),1), X2(3*X1_sorted(Esi(i,3),2),1))...
                + (2*pi*freq*inum)*0.25*q*(lenEsi(i,3)*edgeEsi(i,3) + lenEsi(i,2)*edgeEsi(i,2));
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Residue matrix
Res = zeros(size(Jaco,1),1);

%%% perturbed
% total vertex : 671
node = 336; % node

% Res(3*node-2,1) = 1; % potential disturbed 
% Res(3*node-1,1) = 1; % electron disturbed 
Res(3*node,1) = 1; % hole disturbed 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Dirichlet boundary 
for i=1:(Ncon1+Ncon2)
    for j=1:concol
        fcon = X1_sorted(contact(i,j),2);
        fcon_np = find(contact(i,j)==V_silicon(:,1));

        Jaco(X2(3*fcon-2),:) = 0;
        Jaco(X2(3*fcon-1),:) = 0;
        Jaco(X2(3*fcon),:) = 0;

        Jaco(X2(3*fcon-2), X2(3*fcon-2)) = 1;
        Jaco(X2(3*fcon-1), X2(3*fcon-1)) = 1;
        Jaco(X2(3*fcon), X2(3*fcon)) = 1;

        if i<=Ncon1 % Cathode
            Res(X2(3*fcon-2),1) = 0;
            Res(X2(3*fcon-1),1) = 0;
            Res(X2(3*fcon),1) = 0;

        elseif i>(Ncon1) &&  i<=(Ncon1+Ncon2) % Anode
            Res(X2(3*fcon-2),1) = 0;
            Res(X2(3*fcon-1),1) = 0;
            Res(X2(3*fcon),1) = 0;
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Scaled part

Cvector = zeros(size(Jaco,1),1);
for i=1:size(Jaco,1)
    if rem(X2_sorted(i),3) == 1
        Cvector(i) = Vt;
    elseif rem(X2_sorted(i),3) == 2
        Cvector(i) = Nd;
    elseif rem(X2_sorted(i),3) == 0
        Cvector(i) = Nd;
    end
end

Cmatrix = spdiags(Cvector,0,size(X2_sorted,1),size(X2_sorted,1));
Jaco_scaled = Jaco * Cmatrix;
Rvector = 1./sum(abs(Jaco_scaled),2);
Rmatrix = spdiags(Rvector,0,size(X2_sorted,1),size(X2_sorted,1));

Jaco_scaled = Rmatrix * Jaco_scaled;
Res_scaled = Rmatrix * Res;

sol_scaled = Jaco_scaled \ (Res_scaled);
sol = Cmatrix * sol_scaled;

for ii=1:VErow
    dphi(ii,1) = sol(3*ii-2,1);
    dn(ii,1) = sol(3*ii-1,1);
    dp(ii,1) = sol(3*ii,1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Anode current calculation

for i = 1:size(uniconin,1)
    for j = 1:E4row

        fcon = find(uniconin(i) == E4(j,:));
        f1 = find(V_silicon(:,1) == E4(j,1));
        f2 = find(V_silicon(:,1) == E4(j,2));
        f3 = find(V_silicon(:,1) == E4(j,3));

        if fcon == 1
            n_DC_av1 = 0.5*(n_DC(f1)+n_DC(f2));
            n_DC_av2 = 0.5*(n_DC(f1)+n_DC(f3));
            dn_av1 = 0.5*(dn(f1)+dn(f2));
            dn_av2 = 0.5*(dn(f1)+dn(f3));
            ddn1 = dn(f2)-dn(f1);
            ddn2 = dn(f3)-dn(f1);

            p_DC_av1 = 0.5*(p_DC(f1)+p_DC(f2));
            p_DC_av2 = 0.5*(p_DC(f1)+p_DC(f3));
            dp_av1 = 0.5*(dp(f1)+dp(f2));
            dp_av2 = 0.5*(dp(f1)+dp(f3));
            ddp1 = dp(f2)-dp(f1);
            ddp2 = dp(f3)-dp(f1);

            ddphi1 = dphi(X1_sorted(E4(j,2),2))-dphi(X1_sorted(E4(j,1),2));
            ddphi2 = dphi(X1_sorted(E4(j,3),2))-dphi(X1_sorted(E4(j,1),2));
            dphi_DC1 = phi_DC(X1_sorted(E4(j,2),2))-phi_DC(X1_sorted(E4(j,1),2));
            dphi_DC2 = phi_DC(X1_sorted(E4(j,3),2))-phi_DC(X1_sorted(E4(j,1),2));

            In_anode(1,1) = In_anode(1,1) ...
                - ncoeff*edgeE4(j,1)/lenE4(j,1)*(n_DC_av1*ddphi1 + dn_av1*dphi_DC1 - ddn1*Vt)*width...
                - ncoeff*edgeE4(j,3)/lenE4(j,3)*(n_DC_av2*ddphi2 + dn_av2*dphi_DC2 - ddn2*Vt)*width;
            Ip_anode(1,1) = Ip_anode(1,1) ...
                - pcoeff*edgeE4(j,1)/lenE4(j,1)*(p_DC_av1*ddphi1 + dp_av1*dphi_DC1 + ddp1*Vt)*width...
                - pcoeff*edgeE4(j,3)/lenE4(j,3)*(p_DC_av2*ddphi2 + dp_av2*dphi_DC2 + ddp2*Vt)*width;
            Id_anode(1,1) = Id_anode(1,1) ...
                - esi*e0*edgeE4(j,1)/lenE4(j,1)*(2*pi*freq*inum)*ddphi1*width...
                - esi*e0*edgeE4(j,3)/lenE4(j,3)*(2*pi*freq*inum)*ddphi2*width;

        elseif fcon == 2
            n_DC_av1 = 0.5*(n_DC(f3)+n_DC(f2));
            n_DC_av2 = 0.5*(n_DC(f1)+n_DC(f2));
            dn_av1 = 0.5*(dn(f3)+dn(f2));
            dn_av2 = 0.5*(dn(f1)+dn(f2));
            ddn1 = dn(f3)-dn(f2);
            ddn2 = dn(f1)-dn(f2);

            p_DC_av1 = 0.5*(p_DC(f2)+p_DC(f3));
            p_DC_av2 = 0.5*(p_DC(f2)+p_DC(f1));
            dp_av1 = 0.5*(dp(f2)+dp(f3));
            dp_av2 = 0.5*(dp(f2)+dp(f1));
            ddp1 = dp(f3)-dp(f2);
            ddp2 = dp(f1)-dp(f2);

            ddphi1 = dphi(X1_sorted(E4(j,3),2))-dphi(X1_sorted(E4(j,2),2));
            ddphi2 = dphi(X1_sorted(E4(j,1),2))-dphi(X1_sorted(E4(j,2),2));
            dphi_DC1 = phi_DC(X1_sorted(E4(j,3),2))-phi_DC(X1_sorted(E4(j,2),2));
            dphi_DC2 = phi_DC(X1_sorted(E4(j,1),2))-phi_DC(X1_sorted(E4(j,2),2));

            In_anode(1,1) = In_anode(1,1) ...
                - ncoeff*edgeE4(j,2)/lenE4(j,2)*(n_DC_av1*ddphi1 + dn_av1*dphi_DC1 - ddn1*Vt)*width...
                - ncoeff*edgeE4(j,1)/lenE4(j,1)*(n_DC_av2*ddphi2 + dn_av2*dphi_DC2 - ddn2*Vt)*width;
            Ip_anode(1,1) = Ip_anode(1,1) ...
                - pcoeff*edgeE4(j,2)/lenE4(j,2)*(p_DC_av1*ddphi1 + dp_av1*dphi_DC1 + ddp1*Vt)*width...
                - pcoeff*edgeE4(j,1)/lenE4(j,1)*(p_DC_av2*ddphi2 + dp_av2*dphi_DC2 + ddp2*Vt)*width;
            Id_anode(1,1) = Id_anode(1,1) ...
                - esi*e0*edgeE4(j,2)/lenE4(j,2)*(2*pi*freq*inum)*ddphi1*width...
                - esi*e0*edgeE4(j,1)/lenE4(j,1)*(2*pi*freq*inum)*ddphi2*width;

        elseif fcon == 3
            n_DC_av1 = 0.5*(n_DC(f3)+n_DC(f1));
            n_DC_av2 = 0.5*(n_DC(f3)+n_DC(f2));
            dn_av1 = 0.5*(dn(f3)+dn(f1));
            dn_av2 = 0.5*(dn(f3)+dn(f2));
            ddn1 = dn(f1)-dn(f3);
            ddn2 = dn(f2)-dn(f3);

            p_DC_av1 = 0.5*(p_DC(f3)+p_DC(f1));
            p_DC_av2 = 0.5*(p_DC(f3)+p_DC(f2));
            dp_av1 = 0.5*(dp(f3)+dp(f1));
            dp_av2 = 0.5*(dp(f3)+dp(f2));
            ddp1 = dp(f1)-dp(f3);
            ddp2 = dp(f2)-dp(f3);

            ddphi1 = dphi(X1_sorted(E4(j,1),2))-dphi(X1_sorted(E4(j,3),2));
            ddphi2 = dphi(X1_sorted(E4(j,2),2))-dphi(X1_sorted(E4(j,3),2));
            dphi_DC1 = phi_DC(X1_sorted(E4(j,1),2))-phi_DC(X1_sorted(E4(j,3),2));
            dphi_DC2 = phi_DC(X1_sorted(E4(j,2),2))-phi_DC(X1_sorted(E4(j,3),2));

            In_anode(1,1) = In_anode(1,1) ...
                - ncoeff*edgeE4(j,3)/lenE4(j,3)*(n_DC_av1*ddphi1 + dn_av1*dphi_DC1 - ddn1*Vt)*width...
                - ncoeff*edgeE4(j,2)/lenE4(j,2)*(n_DC_av2*ddphi2 + dn_av2*dphi_DC2 - ddn2*Vt)*width;
            Ip_anode(1,1) = Ip_anode(1,1) ...
                - pcoeff*edgeE4(j,3)/lenE4(j,3)*(p_DC_av1*ddphi1 + dp_av1*dphi_DC1 + ddp1*Vt)*width...
                - pcoeff*edgeE4(j,2)/lenE4(j,2)*(p_DC_av2*ddphi2 + dp_av2*dphi_DC2 + ddp2*Vt)*width;
            Id_anode(1,1) = Id_anode(1,1) ...
                - esi*e0*edgeE4(j,3)/lenE4(j,3)*(2*pi*freq*inum)*ddphi1*width...
                - esi*e0*edgeE4(j,2)/lenE4(j,2)*(2*pi*freq*inum)*ddphi2*width;
        end
    end
end

I_anode = In_anode + Ip_anode + Id_anode;

ansmat_anode = zeros(4,1);
ansmat_anode(1,1) = max(I_anode);
ansmat_anode(2,1) = max(In_anode);
ansmat_anode(3,1) = max(Ip_anode);
ansmat_anode(4,1) = max(Id_anode);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Cathode current calculation

for i = 1:size(uniconout,1)
    for j = 1:E2row

        fcon = find(uniconout(i) == E2(j,:));
        f1 = find(V_silicon(:,1) == E2(j,1));
        f2 = find(V_silicon(:,1) == E2(j,2));
        f3 = find(V_silicon(:,1) == E2(j,3));

        if fcon == 1
            n_DC_av1 = 0.5*(n_DC(f1)+n_DC(f2));
            n_DC_av2 = 0.5*(n_DC(f1)+n_DC(f3));
            dn_av1 = 0.5*(dn(f1)+dn(f2));
            dn_av2 = 0.5*(dn(f1)+dn(f3));
            ddn1 = dn(f2)-dn(f1);
            ddn2 = dn(f3)-dn(f1);

            p_DC_av1 = 0.5*(p_DC(f1)+p_DC(f2));
            p_DC_av2 = 0.5*(p_DC(f1)+p_DC(f3));
            dp_av1 = 0.5*(dp(f1)+dp(f2));
            dp_av2 = 0.5*(dp(f1)+dp(f3));
            ddp1 = dp(f2)-dp(f1);
            ddp2 = dp(f3)-dp(f1);

            ddphi1 = dphi(X1_sorted(E2(j,2),2))-dphi(X1_sorted(E2(j,1),2));
            ddphi2 = dphi(X1_sorted(E2(j,3),2))-dphi(X1_sorted(E2(j,1),2));
            dphi_DC1 = phi_DC(X1_sorted(E2(j,2),2))-phi_DC(X1_sorted(E2(j,1),2));
            dphi_DC2 = phi_DC(X1_sorted(E2(j,3),2))-phi_DC(X1_sorted(E2(j,1),2));

            In_cathode(1,1) = In_cathode(1,1) ...
                - ncoeff*edgeE2(j,1)/lenE2(j,1)*(n_DC_av1*ddphi1 + dn_av1*dphi_DC1 - ddn1*Vt)*width...
                - ncoeff*edgeE2(j,3)/lenE2(j,3)*(n_DC_av2*ddphi2 + dn_av2*dphi_DC2 - ddn2*Vt)*width;
            Ip_cathode(1,1) = Ip_cathode(1,1) ...
                - pcoeff*edgeE2(j,1)/lenE2(j,1)*(p_DC_av1*ddphi1 + dp_av1*dphi_DC1 + ddp1*Vt)*width...
                - pcoeff*edgeE2(j,3)/lenE2(j,3)*(p_DC_av2*ddphi2 + dp_av2*dphi_DC2 + ddp2*Vt)*width;
            Id_cathode(1,1) = Id_cathode(1,1) ...
                - esi*e0*edgeE2(j,1)/lenE2(j,1)*(2*pi*freq*inum)*ddphi1*width...
                - esi*e0*edgeE2(j,3)/lenE2(j,3)*(2*pi*freq*inum)*ddphi2*width;

        elseif fcon == 2
            n_DC_av1 = 0.5*(n_DC(f3)+n_DC(f2));
            n_DC_av2 = 0.5*(n_DC(f1)+n_DC(f2));
            dn_av1 = 0.5*(dn(f3)+dn(f2));
            dn_av2 = 0.5*(dn(f1)+dn(f2));
            ddn1 = dn(f3)-dn(f2);
            ddn2 = dn(f1)-dn(f2);

            p_DC_av1 = 0.5*(p_DC(f2)+p_DC(f3));
            p_DC_av2 = 0.5*(p_DC(f2)+p_DC(f1));
            dp_av1 = 0.5*(dp(f2)+dp(f3));
            dp_av2 = 0.5*(dp(f2)+dp(f1));
            ddp1 = dp(f3)-dp(f2);
            ddp2 = dp(f1)-dp(f2);

            ddphi1 = dphi(X1_sorted(E2(j,3),2))-dphi(X1_sorted(E2(j,2),2));
            ddphi2 = dphi(X1_sorted(E2(j,1),2))-dphi(X1_sorted(E2(j,2),2));
            dphi_DC1 = phi_DC(X1_sorted(E2(j,3),2))-phi_DC(X1_sorted(E2(j,2),2));
            dphi_DC2 = phi_DC(X1_sorted(E2(j,1),2))-phi_DC(X1_sorted(E2(j,2),2));

            In_cathode(1,1) = In_cathode(1,1) ...
                - ncoeff*edgeE2(j,2)/lenE2(j,2)*(n_DC_av1*ddphi1 + dn_av1*dphi_DC1 - ddn1*Vt)*width...
                - ncoeff*edgeE2(j,1)/lenE2(j,1)*(n_DC_av2*ddphi2 + dn_av2*dphi_DC2 - ddn2*Vt)*width;
            Ip_cathode(1,1) = Ip_cathode(1,1) ...
                - pcoeff*edgeE2(j,2)/lenE2(j,2)*(p_DC_av1*ddphi1 + dp_av1*dphi_DC1 + ddp1*Vt)*width...
                - pcoeff*edgeE2(j,1)/lenE2(j,1)*(p_DC_av2*ddphi2 + dp_av2*dphi_DC2 + ddp2*Vt)*width;
            Id_cathode(1,1) = Id_cathode(1,1) ...
                - esi*e0*edgeE2(j,2)/lenE2(j,2)*(2*pi*freq*inum)*ddphi1*width...
                - esi*e0*edgeE2(j,1)/lenE2(j,1)*(2*pi*freq*inum)*ddphi2*width;

        elseif fcon == 3
            n_DC_av1 = 0.5*(n_DC(f3)+n_DC(f1));
            n_DC_av2 = 0.5*(n_DC(f3)+n_DC(f2));
            dn_av1 = 0.5*(dn(f3)+dn(f1));
            dn_av2 = 0.5*(dn(f3)+dn(f2));
            ddn1 = dn(f1)-dn(f3);
            ddn2 = dn(f2)-dn(f3);

            p_DC_av1 = 0.5*(p_DC(f3)+p_DC(f1));
            p_DC_av2 = 0.5*(p_DC(f3)+p_DC(f2));
            dp_av1 = 0.5*(dp(f3)+dp(f1));
            dp_av2 = 0.5*(dp(f3)+dp(f2));
            ddp1 = dp(f1)-dp(f3);
            ddp2 = dp(f2)-dp(f3);

            ddphi1 = dphi(X1_sorted(E2(j,1),2))-dphi(X1_sorted(E2(j,3),2));
            ddphi2 = dphi(X1_sorted(E2(j,2),2))-dphi(X1_sorted(E2(j,3),2));
            dphi_DC1 = phi_DC(X1_sorted(E2(j,1),2))-phi_DC(X1_sorted(E2(j,3),2));
            dphi_DC2 = phi_DC(X1_sorted(E2(j,2),2))-phi_DC(X1_sorted(E2(j,3),2));

            In_cathode(1,1) = In_cathode(1,1) ...
                - ncoeff*edgeE2(j,3)/lenE2(j,3)*(n_DC_av1*ddphi1 + dn_av1*dphi_DC1 - ddn1*Vt)*width...
                - ncoeff*edgeE2(j,2)/lenE2(j,2)*(n_DC_av2*ddphi2 + dn_av2*dphi_DC2 - ddn2*Vt)*width;
            Ip_cathode(1,1) = Ip_cathode(1,1) ...
                - pcoeff*edgeE2(j,3)/lenE2(j,3)*(p_DC_av1*ddphi1 + dp_av1*dphi_DC1 + ddp1*Vt)*width...
                - pcoeff*edgeE2(j,2)/lenE2(j,2)*(p_DC_av2*ddphi2 + dp_av2*dphi_DC2 + ddp2*Vt)*width;
            Id_cathode(1,1) = Id_cathode(1,1) ...
                - esi*e0*edgeE2(j,3)/lenE2(j,3)*(2*pi*freq*inum)*ddphi1*width...
                - esi*e0*edgeE2(j,2)/lenE2(j,2)*(2*pi*freq*inum)*ddphi2*width;
        end
    end
end

I_cathode = In_cathode + Ip_cathode + Id_cathode;

ansmat_cathode = zeros(4,1);
ansmat_cathode(1,1) = max(I_cathode);
ansmat_cathode(2,1) = max(In_cathode);
ansmat_cathode(3,1) = max(Ip_cathode);
ansmat_cathode(4,1) = max(Id_cathode);

ans_I_sum = ansmat_cathode(1,1) + ansmat_anode(1,1);

ans_total_charge = 0;
for i=1:size(VE,1)
    ans_total_charge = ans_total_charge - q*(-dn(i,1)+dp(i,1))*saveexl(i,1)*inum*(2*pi*freq);
end

Jd_anode = Id_anode / width;
Jd_cathode = Id_cathode / width;
ans_Jd_sum = Jd_anode + Jd_cathode;
