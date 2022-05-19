clear
load("ramping_information.mat")

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  homo AC %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% transient edit part
freq = 1e12;
cycle = 2;
T = 1/freq;
Vamp = 1e-3;
Nstep = 100;
iN=sqrt(-1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ncoeff = q*mun;
pcoeff = q*mup;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dphi = zeros(size(phi));
dn = nint*exp(dphi);
dp = nint*exp(dphi);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
save_phi_update = zeros(20,1);
save_update = zeros(size(X2_sorted,1),1);

In = zeros(1,1);
Ip = zeros(1,1);
Id = zeros(1,1);

Jaco_ac = sparse(size(X2_sorted,1),size(X2_sorted,1));
Res_ac = zeros(size(X2_sorted,1),1);

saveexl = zeros(size(X1,1),1);

for i=1:Esirow
    if i <= E2row
        Jaco_ac(X2(3*(X1_sorted(Esi(i,1),2))-2,1), X2(3*(X1_sorted(Esi(i,1),2))-2,1)) = Jaco_ac(X2(3*(X1_sorted(Esi(i,1),2))-2,1), X2(3*(X1_sorted(Esi(i,1),2))-2,1))...
            - esi*(edgeEsi(i,1)/lenEsi(i,1) + edgeEsi(i,3)/lenEsi(i,3));
        Jaco_ac(X2(3*(X1_sorted(Esi(i,1),2))-2,1), X2(3*(X1_sorted(Esi(i,2),2))-2,1)) = Jaco_ac(X2(3*(X1_sorted(Esi(i,1),2))-2,1), X2(3*(X1_sorted(Esi(i,2),2))-2,1))...
            + esi*edgeEsi(i,1)/lenEsi(i,1);
        Jaco_ac(X2(3*(X1_sorted(Esi(i,1),2))-2,1), X2(3*(X1_sorted(Esi(i,3),2))-2,1)) = Jaco_ac(X2(3*(X1_sorted(Esi(i,1),2))-2,1), X2(3*(X1_sorted(Esi(i,3),2))-2,1))...
            + esi*edgeEsi(i,3)/lenEsi(i,3);
        Res_ac(X2(3*(X1_sorted(Esi(i,1),2))-2,1),1) =  Res_ac(X2(3*(X1_sorted(Esi(i,1),2))-2,1),1)...
            + esi*edgeEsi(i,1)/lenEsi(i,1)*(dphi(X1_sorted(Esi(i,2),2))-dphi(X1_sorted(Esi(i,1),2)))...
            + esi*edgeEsi(i,3)/lenEsi(i,3)*(dphi(X1_sorted(Esi(i,3),2))-dphi(X1_sorted(Esi(i,1),2)));
        saveexl(X1_sorted(Esi(i,1),2)) = saveexl(X1_sorted(Esi(i,1),2)) + (edgeEsi(i,1)*lenEsi(i,1) + edgeEsi(i,3)*lenEsi(i,3))/4;

        Jaco_ac(X2(3*(X1_sorted(Esi(i,2),2))-2,1), X2(3*(X1_sorted(Esi(i,2),2))-2,1)) = Jaco_ac(X2(3*(X1_sorted(Esi(i,2),2))-2,1), X2(3*(X1_sorted(Esi(i,2),2))-2,1))...
            - esi*(edgeEsi(i,1)/lenEsi(i,1) + edgeEsi(i,2)/lenEsi(i,2));
        Jaco_ac(X2(3*(X1_sorted(Esi(i,2),2))-2,1), X2(3*(X1_sorted(Esi(i,1),2))-2,1)) = Jaco_ac(X2(3*(X1_sorted(Esi(i,2),2))-2,1), X2(3*(X1_sorted(Esi(i,1),2))-2,1))...
            + esi*edgeEsi(i,1)/lenEsi(i,1);
        Jaco_ac(X2(3*(X1_sorted(Esi(i,2),2))-2,1), X2(3*(X1_sorted(Esi(i,3),2))-2,1)) = Jaco_ac(X2(3*(X1_sorted(Esi(i,2),2))-2,1), X2(3*(X1_sorted(Esi(i,3),2))-2,1))...
            + esi*edgeEsi(i,2)/lenEsi(i,2);
        Res_ac(X2(3*(X1_sorted(Esi(i,2),2))-2,1),1) =  Res_ac(X2(3*(X1_sorted(Esi(i,2),2))-2,1),1)...
            + esi*edgeEsi(i,1)/lenEsi(i,1)*(dphi(X1_sorted(Esi(i,1),2))-dphi(X1_sorted(Esi(i,2),2)))...
            + esi*edgeEsi(i,2)/lenEsi(i,2)*(dphi(X1_sorted(Esi(i,3),2))-dphi(X1_sorted(Esi(i,2),2)));
        saveexl(X1_sorted(Esi(i,2),2)) = saveexl(X1_sorted(Esi(i,2),2)) + (edgeEsi(i,1)*lenEsi(i,1) + edgeEsi(i,2)*lenEsi(i,2))/4;

        Jaco_ac(X2(3*(X1_sorted(Esi(i,3),2))-2,1), X2(3*(X1_sorted(Esi(i,3),2))-2,1)) = Jaco_ac(X2(3*(X1_sorted(Esi(i,3),2))-2,1), X2(3*(X1_sorted(Esi(i,3),2))-2,1))...
            - esi*(edgeEsi(i,2)/lenEsi(i,2) + edgeEsi(i,3)/lenEsi(i,3));
        Jaco_ac(X2(3*(X1_sorted(Esi(i,3),2))-2,1), X2(3*(X1_sorted(Esi(i,2),2))-2,1)) = Jaco_ac(X2(3*(X1_sorted(Esi(i,3),2))-2,1), X2(3*(X1_sorted(Esi(i,2),2))-2,1))...
            + esi*edgeEsi(i,2)/lenEsi(i,2);
        Jaco_ac(X2(3*(X1_sorted(Esi(i,3),2))-2,1), X2(3*(X1_sorted(Esi(i,1),2))-2,1)) = Jaco_ac(X2(3*(X1_sorted(Esi(i,3),2))-2,1), X2(3*(X1_sorted(Esi(i,1),2))-2,1))...
            + esi*edgeEsi(i,3)/lenEsi(i,3);
        Res_ac(X2(3*(X1_sorted(Esi(i,3),2))-2,1),1) =  Res_ac(X2(3*(X1_sorted(Esi(i,3),2))-2,1),1)...
            + esi*edgeEsi(i,2)/lenEsi(i,2)*(dphi(X1_sorted(Esi(i,2),2))-dphi(X1_sorted(Esi(i,3),2)))...
            + esi*edgeEsi(i,3)/lenEsi(i,3)*(dphi(X1_sorted(Esi(i,1),2))-dphi(X1_sorted(Esi(i,3),2)));
        saveexl(X1_sorted(Esi(i,3),2)) = saveexl(X1_sorted(Esi(i,3),2)) + (edgeEsi(i,2)*lenEsi(i,2) + edgeEsi(i,3)*lenEsi(i,3))/4;

    elseif i > E2row && i <= (E2row + E3row)
        Jaco_ac(X2(3*(X1_sorted(Esi(i,1),2))-2,1), X2(3*(X1_sorted(Esi(i,1),2))-2,1)) = Jaco_ac(X2(3*(X1_sorted(Esi(i,1),2))-2,1), X2(3*(X1_sorted(Esi(i,1),2))-2,1))...
            - esi*(edgeEsi(i,1)/lenEsi(i,1) + edgeEsi(i,3)/lenEsi(i,3));
        Jaco_ac(X2(3*(X1_sorted(Esi(i,1),2))-2,1), X2(3*(X1_sorted(Esi(i,2),2))-2,1)) = Jaco_ac(X2(3*(X1_sorted(Esi(i,1),2))-2,1), X2(3*(X1_sorted(Esi(i,2),2))-2,1))...
            + esi*edgeEsi(i,1)/lenEsi(i,1);
        Jaco_ac(X2(3*(X1_sorted(Esi(i,1),2))-2,1), X2(3*(X1_sorted(Esi(i,3),2))-2,1)) = Jaco_ac(X2(3*(X1_sorted(Esi(i,1),2))-2,1), X2(3*(X1_sorted(Esi(i,3),2))-2,1))...
            + esi*edgeEsi(i,3)/lenEsi(i,3);
        Res_ac(X2(3*(X1_sorted(Esi(i,1),2))-2,1),1) =  Res_ac(X2(3*(X1_sorted(Esi(i,1),2))-2,1),1)...
            + esi*edgeEsi(i,1)/lenEsi(i,1)*(dphi(X1_sorted(Esi(i,2),2))-dphi(X1_sorted(Esi(i,1),2)))...
            + esi*edgeEsi(i,3)/lenEsi(i,3)*(dphi(X1_sorted(Esi(i,3),2))-dphi(X1_sorted(Esi(i,1),2)));
        saveexl(X1_sorted(Esi(i,1),2)) = saveexl(X1_sorted(Esi(i,1),2)) + (edgeEsi(i,1)*lenEsi(i,1) + edgeEsi(i,3)*lenEsi(i,3))/4;

        Jaco_ac(X2(3*(X1_sorted(Esi(i,2),2))-2,1), X2(3*(X1_sorted(Esi(i,2),2))-2,1)) = Jaco_ac(X2(3*(X1_sorted(Esi(i,2),2))-2,1), X2(3*(X1_sorted(Esi(i,2),2))-2,1))...
            - esi*(edgeEsi(i,1)/lenEsi(i,1) + edgeEsi(i,2)/lenEsi(i,2));
        Jaco_ac(X2(3*(X1_sorted(Esi(i,2),2))-2,1), X2(3*(X1_sorted(Esi(i,1),2))-2,1)) = Jaco_ac(X2(3*(X1_sorted(Esi(i,2),2))-2,1), X2(3*(X1_sorted(Esi(i,1),2))-2,1))...
            + esi*edgeEsi(i,1)/lenEsi(i,1);
        Jaco_ac(X2(3*(X1_sorted(Esi(i,2),2))-2,1), X2(3*(X1_sorted(Esi(i,3),2))-2,1)) = Jaco_ac(X2(3*(X1_sorted(Esi(i,2),2))-2,1), X2(3*(X1_sorted(Esi(i,3),2))-2,1))...
            + esi*edgeEsi(i,2)/lenEsi(i,2);
        Res_ac(X2(3*(X1_sorted(Esi(i,2),2))-2,1),1) =  Res_ac(X2(3*(X1_sorted(Esi(i,2),2))-2,1),1)...
            + esi*edgeEsi(i,1)/lenEsi(i,1)*(dphi(X1_sorted(Esi(i,1),2))-dphi(X1_sorted(Esi(i,2),2)))...
            + esi*edgeEsi(i,2)/lenEsi(i,2)*(dphi(X1_sorted(Esi(i,3),2))-dphi(X1_sorted(Esi(i,2),2)));
        saveexl(X1_sorted(Esi(i,2),2)) = saveexl(X1_sorted(Esi(i,2),2)) + (edgeEsi(i,1)*lenEsi(i,1) + edgeEsi(i,2)*lenEsi(i,2))/4;

        Jaco_ac(X2(3*(X1_sorted(Esi(i,3),2))-2,1), X2(3*(X1_sorted(Esi(i,3),2))-2,1)) = Jaco_ac(X2(3*(X1_sorted(Esi(i,3),2))-2,1), X2(3*(X1_sorted(Esi(i,3),2))-2,1))...
            - esi*(edgeEsi(i,2)/lenEsi(i,2) + edgeEsi(i,3)/lenEsi(i,3));
        Jaco_ac(X2(3*(X1_sorted(Esi(i,3),2))-2,1), X2(3*(X1_sorted(Esi(i,2),2))-2,1)) = Jaco_ac(X2(3*(X1_sorted(Esi(i,3),2))-2,1), X2(3*(X1_sorted(Esi(i,2),2))-2,1))...
            + esi*edgeEsi(i,2)/lenEsi(i,2);
        Jaco_ac(X2(3*(X1_sorted(Esi(i,3),2))-2,1), X2(3*(X1_sorted(Esi(i,1),2))-2,1)) = Jaco_ac(X2(3*(X1_sorted(Esi(i,3),2))-2,1), X2(3*(X1_sorted(Esi(i,1),2))-2,1))...
            + esi*edgeEsi(i,3)/lenEsi(i,3);
        Res_ac(X2(3*(X1_sorted(Esi(i,3),2))-2,1),1) =  Res_ac(X2(3*(X1_sorted(Esi(i,3),2))-2,1),1)...
            + esi*edgeEsi(i,2)/lenEsi(i,2)*(dphi(X1_sorted(Esi(i,2),2))-dphi(X1_sorted(Esi(i,3),2)))...
            + esi*edgeEsi(i,3)/lenEsi(i,3)*(dphi(X1_sorted(Esi(i,1),2))-dphi(X1_sorted(Esi(i,3),2)));
        saveexl(X1_sorted(Esi(i,3),2)) = saveexl(X1_sorted(Esi(i,3),2)) + (edgeEsi(i,2)*lenEsi(i,2) + edgeEsi(i,3)*lenEsi(i,3))/4;

    elseif i > (E2row+E3row) && i <= Esirow
        Jaco_ac(X2(3*(X1_sorted(Esi(i,1),2))-2,1), X2(3*(X1_sorted(Esi(i,1),2))-2,1)) = Jaco_ac(X2(3*(X1_sorted(Esi(i,1),2))-2,1), X2(3*(X1_sorted(Esi(i,1),2))-2,1))...
            - esi*(edgeEsi(i,1)/lenEsi(i,1) + edgeEsi(i,3)/lenEsi(i,3));
        Jaco_ac(X2(3*(X1_sorted(Esi(i,1),2))-2,1), X2(3*(X1_sorted(Esi(i,2),2))-2,1)) = Jaco_ac(X2(3*(X1_sorted(Esi(i,1),2))-2,1), X2(3*(X1_sorted(Esi(i,2),2))-2,1))...
            + esi*edgeEsi(i,1)/lenEsi(i,1);
        Jaco_ac(X2(3*(X1_sorted(Esi(i,1),2))-2,1), X2(3*(X1_sorted(Esi(i,3),2))-2,1)) = Jaco_ac(X2(3*(X1_sorted(Esi(i,1),2))-2,1), X2(3*(X1_sorted(Esi(i,3),2))-2,1))...
            + esi*edgeEsi(i,3)/lenEsi(i,3);
        Res_ac(X2(3*(X1_sorted(Esi(i,1),2))-2,1),1) =  Res_ac(X2(3*(X1_sorted(Esi(i,1),2))-2,1),1)...
            + esi*edgeEsi(i,1)/lenEsi(i,1)*(dphi(X1_sorted(Esi(i,2),2))-dphi(X1_sorted(Esi(i,1),2)))...
            + esi*edgeEsi(i,3)/lenEsi(i,3)*(dphi(X1_sorted(Esi(i,3),2))-dphi(X1_sorted(Esi(i,1),2)));
        saveexl(X1_sorted(Esi(i,1),2)) = saveexl(X1_sorted(Esi(i,1),2)) + (edgeEsi(i,1)*lenEsi(i,1) + edgeEsi(i,3)*lenEsi(i,3))/4;

        Jaco_ac(X2(3*(X1_sorted(Esi(i,2),2))-2,1), X2(3*(X1_sorted(Esi(i,2),2))-2,1)) = Jaco_ac(X2(3*(X1_sorted(Esi(i,2),2))-2,1), X2(3*(X1_sorted(Esi(i,2),2))-2,1))...
            - esi*(edgeEsi(i,1)/lenEsi(i,1) + edgeEsi(i,2)/lenEsi(i,2));
        Jaco_ac(X2(3*(X1_sorted(Esi(i,2),2))-2,1), X2(3*(X1_sorted(Esi(i,1),2))-2,1)) = Jaco_ac(X2(3*(X1_sorted(Esi(i,2),2))-2,1), X2(3*(X1_sorted(Esi(i,1),2))-2,1))...
            + esi*edgeEsi(i,1)/lenEsi(i,1);
        Jaco_ac(X2(3*(X1_sorted(Esi(i,2),2))-2,1), X2(3*(X1_sorted(Esi(i,3),2))-2,1)) = Jaco_ac(X2(3*(X1_sorted(Esi(i,2),2))-2,1), X2(3*(X1_sorted(Esi(i,3),2))-2,1))...
            + esi*edgeEsi(i,2)/lenEsi(i,2);
        Res_ac(X2(3*(X1_sorted(Esi(i,2),2))-2,1),1) =  Res_ac(X2(3*(X1_sorted(Esi(i,2),2))-2,1),1)...
            + esi*edgeEsi(i,1)/lenEsi(i,1)*(dphi(X1_sorted(Esi(i,1),2))-dphi(X1_sorted(Esi(i,2),2)))...
            + esi*edgeEsi(i,2)/lenEsi(i,2)*(dphi(X1_sorted(Esi(i,3),2))-dphi(X1_sorted(Esi(i,2),2)));
        saveexl(X1_sorted(Esi(i,2),2)) = saveexl(X1_sorted(Esi(i,2),2)) + (edgeEsi(i,1)*lenEsi(i,1) + edgeEsi(i,2)*lenEsi(i,2))/4;

        Jaco_ac(X2(3*(X1_sorted(Esi(i,3),2))-2,1), X2(3*(X1_sorted(Esi(i,3),2))-2,1)) = Jaco_ac(X2(3*(X1_sorted(Esi(i,3),2))-2,1), X2(3*(X1_sorted(Esi(i,3),2))-2,1))...
            - esi*(edgeEsi(i,2)/lenEsi(i,2) + edgeEsi(i,3)/lenEsi(i,3));
        Jaco_ac(X2(3*(X1_sorted(Esi(i,3),2))-2,1), X2(3*(X1_sorted(Esi(i,2),2))-2,1)) = Jaco_ac(X2(3*(X1_sorted(Esi(i,3),2))-2,1), X2(3*(X1_sorted(Esi(i,2),2))-2,1))...
            + esi*edgeEsi(i,2)/lenEsi(i,2);
        Jaco_ac(X2(3*(X1_sorted(Esi(i,3),2))-2,1), X2(3*(X1_sorted(Esi(i,1),2))-2,1)) = Jaco_ac(X2(3*(X1_sorted(Esi(i,3),2))-2,1), X2(3*(X1_sorted(Esi(i,1),2))-2,1))...
            + esi*edgeEsi(i,3)/lenEsi(i,3);
        Res_ac(X2(3*(X1_sorted(Esi(i,3),2))-2,1),1) =  Res_ac(X2(3*(X1_sorted(Esi(i,3),2))-2,1),1)...
            + esi*edgeEsi(i,2)/lenEsi(i,2)*(dphi(X1_sorted(Esi(i,2),2))-dphi(X1_sorted(Esi(i,3),2)))...
            + esi*edgeEsi(i,3)/lenEsi(i,3)*(dphi(X1_sorted(Esi(i,1),2))-dphi(X1_sorted(Esi(i,3),2)));
        saveexl(X1_sorted(Esi(i,3),2)) = saveexl(X1_sorted(Esi(i,3),2)) + (edgeEsi(i,2)*lenEsi(i,2) + edgeEsi(i,3)*lenEsi(i,3))/4;
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  for electron
for i=1:V_silicon_row
    % - qn
    Jaco_ac(X2(3*V_silicon(i,2)-2),X2(3*V_silicon(i,2)-1)) = Jaco_ac(X2(3*V_silicon(i,2)-2),X2(3*V_silicon(i,2)-1)) - saveexl(V_silicon(i,2))*q/e0;
    Res_ac(X2(3*V_silicon(i,2)-2),1) = Res_ac(X2(3*V_silicon(i,2)-2),1) - saveexl(V_silicon(i,2))*q/e0*dn(i);

    % + qp
    Jaco_ac(X2(3*V_silicon(i,2)-2),X2(3*V_silicon(i,2))) =  Jaco_ac(X2(3*V_silicon(i,2)-2),X2(3*V_silicon(i,2))) + saveexl(V_silicon(i,2))*q/e0;
    Res_ac(X2(3*V_silicon(i,2)-2),1) = Res_ac(X2(3*V_silicon(i,2)-2),1) + saveexl(V_silicon(i,2))*q/e0*dp(i);
end

%%%  for electron DD
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
            Jaco_ac(X2(3*X1_sorted(Esi(i,1),2)-1,1), X2(3*X1_sorted(Esi(i,1),2)-1,1)) = Jaco_ac(X2(3*X1_sorted(Esi(i,1),2)-1,1), X2(3*X1_sorted(Esi(i,1),2)-1,1))...
                - ncoeff*edgeEsi(i,1)/lenEsi(i,1)*(0.5*dphi_DC1 + Vt)...
                - ncoeff*edgeEsi(i,3)/lenEsi(i,3)*(0.5*dphi_DC2 + Vt);

            Jaco_ac(X2(3*X1_sorted(Esi(i,1),2)-1,1), X2(3*X1_sorted(Esi(i,2),2)-1,1)) = Jaco_ac(X2(3*X1_sorted(Esi(i,1),2)-1,1), X2(3*X1_sorted(Esi(i,2),2)-1,1))...
                - ncoeff*edgeEsi(i,1)/lenEsi(i,1)*(0.5*dphi_DC1 - Vt);
            Jaco_ac(X2(3*X1_sorted(Esi(i,1),2)-1,1), X2(3*X1_sorted(Esi(i,3),2)-1,1)) = Jaco_ac(X2(3*X1_sorted(Esi(i,1),2)-1,1), X2(3*X1_sorted(Esi(i,3),2)-1,1))...
                - ncoeff*edgeEsi(i,3)/lenEsi(i,3)*(0.5*dphi_DC2 - Vt);

            % phi part
            Jaco_ac(X2(3*X1_sorted(Esi(i,1),2)-1,1), X2(3*X1_sorted(Esi(i,1),2)-2,1)) = Jaco_ac(X2(3*X1_sorted(Esi(i,1),2)-1,1), X2(3*X1_sorted(Esi(i,1),2)-2,1))...
                - ncoeff*edgeEsi(i,1)/lenEsi(i,1)*(-n_DC_av1)...
                - ncoeff*edgeEsi(i,3)/lenEsi(i,3)*(-n_DC_av2);

            Jaco_ac(X2(3*X1_sorted(Esi(i,1),2)-1,1), X2(3*X1_sorted(Esi(i,2),2)-2,1)) = Jaco_ac(X2(3*X1_sorted(Esi(i,1),2)-1,1), X2(3*X1_sorted(Esi(i,2),2)-2,1))...
                - ncoeff*edgeEsi(i,1)/lenEsi(i,1)*n_DC_av1;
            Jaco_ac(X2(3*X1_sorted(Esi(i,1),2)-1,1), X2(3*X1_sorted(Esi(i,3),2)-2,1)) = Jaco_ac(X2(3*X1_sorted(Esi(i,1),2)-1,1), X2(3*X1_sorted(Esi(i,3),2)-2,1))...
                - ncoeff*edgeEsi(i,3)/lenEsi(i,3)*n_DC_av2;

            % n ac part
            Jaco_ac(X2(3*X1_sorted(Esi(i,1),2)-1,1), X2(3*X1_sorted(Esi(i,1),2)-1,1)) = Jaco_ac(X2(3*X1_sorted(Esi(i,1),2)-1,1), X2(3*X1_sorted(Esi(i,1),2)-1,1))...
                - q/4*(lenEsi(i,1)*edgeEsi(i,1) + lenEsi(i,3)*edge(i,3))*(2*pi*freq*iN);

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
            Jaco_ac(X2(3*X1_sorted(Esi(i,2),2)-1,1), X2(3*X1_sorted(Esi(i,2),2)-1,1)) = Jaco_ac(X2(3*X1_sorted(Esi(i,2),2)-1,1), X2(3*X1_sorted(Esi(i,2),2)-1,1))...
                - ncoeff*edgeEsi(i,2)/lenEsi(i,2)*(0.5*dphi_DC1 + Vt)...
                - ncoeff*edgeEsi(i,1)/lenEsi(i,1)*(0.5*dphi_DC2 + Vt);
            Jaco_ac(X2(3*X1_sorted(Esi(i,2),2)-1,1), X2(3*X1_sorted(Esi(i,3),2)-1,1)) = Jaco_ac(X2(3*X1_sorted(Esi(i,2),2)-1,1), X2(3*X1_sorted(Esi(i,3),2)-1,1))...
                - ncoeff*edgeEsi(i,2)/lenEsi(i,2)*(0.5*dphi_DC1 - Vt);
            Jaco_ac(X2(3*X1_sorted(Esi(i,2),2)-1,1), X2(3*X1_sorted(Esi(i,1),2)-1,1)) = Jaco_ac(X2(3*X1_sorted(Esi(i,2),2)-1,1), X2(3*X1_sorted(Esi(i,1),2)-1,1))...
                - ncoeff*edgeEsi(i,1)/lenEsi(i,1)*(0.5*dphi_DC2 - Vt);

            % phi part
            Jaco_ac(X2(3*X1_sorted(Esi(i,2),2)-1,1), X2(3*X1_sorted(Esi(i,2),2)-2,1)) = Jaco_ac(X2(3*X1_sorted(Esi(i,2),2)-1,1), X2(3*X1_sorted(Esi(i,2),2)-2,1))...
                - ncoeff*edgeEsi(i,2)/lenEsi(i,2)*(-n_DC_av1)...
                - ncoeff*edgeEsi(i,1)/lenEsi(i,1)*(-n_DC_av2);
            Jaco_ac(X2(3*X1_sorted(Esi(i,2),2)-1,1), X2(3*X1_sorted(Esi(i,3),2)-2,1)) = Jaco_ac(X2(3*X1_sorted(Esi(i,2),2)-1,1), X2(3*X1_sorted(Esi(i,3),2)-2,1))...
                - ncoeff*edgeEsi(i,2)/lenEsi(i,2)*n_DC_av1;
            Jaco_ac(X2(3*X1_sorted(Esi(i,2),2)-1,1), X2(3*X1_sorted(Esi(i,1),2)-2,1)) = Jaco_ac(X2(3*X1_sorted(Esi(i,2),2)-1,1), X2(3*X1_sorted(Esi(i,1),2)-2,1))...
                - ncoeff*edgeEsi(i,1)/lenEsi(i,1)*n_DC_av2;

            % n ac part
            Jaco_ac(X2(3*X1_sorted(Esi(i,2),2)-1,1), X2(3*X1_sorted(Esi(i,2),2)-1,1)) = Jaco_ac(X2(3*X1_sorted(Esi(i,2),2)-1,1), X2(3*X1_sorted(Esi(i,2),2)-1,1))...
                - q/4*(lenEsi(i,2)*edgeEsi(i,2) + lenEsi(i,1)*edge(i,1))*(2*pi*freq*iN);

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
            Jaco_ac(X2(3*X1_sorted(Esi(i,3),2)-1,1), X2(3*X1_sorted(Esi(i,3),2)-1,1)) = Jaco_ac(X2(3*X1_sorted(Esi(i,3),2)-1,1), X2(3*X1_sorted(Esi(i,3),2)-1,1))...
                - ncoeff*edgeEsi(i,3)/lenEsi(i,3)*(0.5*dphi_DC1 + Vt)...
                - ncoeff*edgeEsi(i,2)/lenEsi(i,2)*(0.5*dphi_DC2 + Vt);
            Jaco_ac(X2(3*X1_sorted(Esi(i,3),2)-1,1), X2(3*X1_sorted(Esi(i,1),2)-1,1)) = Jaco_ac(X2(3*X1_sorted(Esi(i,3),2)-1,1), X2(3*X1_sorted(Esi(i,1),2)-1,1))...
                - ncoeff*edgeEsi(i,3)/lenEsi(i,3)*(0.5*dphi_DC1 - Vt);
            Jaco_ac(X2(3*X1_sorted(Esi(i,3),2)-1,1), X2(3*X1_sorted(Esi(i,2),2)-1,1)) = Jaco_ac(X2(3*X1_sorted(Esi(i,3),2)-1,1), X2(3*X1_sorted(Esi(i,2),2)-1,1))...
                - ncoeff*edgeEsi(i,2)/lenEsi(i,2)*(0.5*dphi_DC2 - Vt);

            % phi part
            Jaco_ac(X2(3*X1_sorted(Esi(i,3),2)-1,1), X2(3*X1_sorted(Esi(i,3),2)-2,1)) = Jaco_ac(X2(3*X1_sorted(Esi(i,3),2)-1,1), X2(3*X1_sorted(Esi(i,3),2)-2,1))...
                - ncoeff*edgeEsi(i,3)/lenEsi(i,3)*(-n_DC_av1)...
                - ncoeff*edgeEsi(i,2)/lenEsi(i,2)*(-n_DC_av2);
            Jaco_ac(X2(3*X1_sorted(Esi(i,3),2)-1,1), X2(3*X1_sorted(Esi(i,1),2)-2,1)) = Jaco_ac(X2(3*X1_sorted(Esi(i,3),2)-1,1), X2(3*X1_sorted(Esi(i,1),2)-2,1))...
                - ncoeff*edgeEsi(i,3)/lenEsi(i,3)*n_DC_av1;
            Jaco_ac(X2(3*X1_sorted(Esi(i,3),2)-1,1), X2(3*X1_sorted(Esi(i,2),2)-2,1)) = Jaco_ac(X2(3*X1_sorted(Esi(i,3),2)-1,1), X2(3*X1_sorted(Esi(i,2),2)-2,1))...
                - ncoeff*edgeEsi(i,2)/lenEsi(i,2)*n_DC_av2;

            % n ac part
            Jaco_ac(X2(3*X1_sorted(Esi(i,3),2)-1,1), X2(3*X1_sorted(Esi(i,3),2)-1,1)) = Jaco_ac(X2(3*X1_sorted(Esi(i,3),2)-1,1), X2(3*X1_sorted(Esi(i,3),2)-1,1))...
                - q/4*(lenEsi(i,3)*edgeEsi(i,3) + lenEsi(i,2)*edge(i,2))*(2*pi*freq*iN);
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% for hole DD
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
            Jaco_ac(X2(3*X1_sorted(Esi(i,1),2),1), X2(3*X1_sorted(Esi(i,1),2),1)) = Jaco_ac(X2(3*X1_sorted(Esi(i,1),2),1), X2(3*X1_sorted(Esi(i,1),2),1))...
                - pcoeff*edgeEsi(i,1)/lenEsi(i,1)*(0.5*dphi_DC1 - Vt)...
                - pcoeff*edgeEsi(i,3)/lenEsi(i,3)*(0.5*dphi_DC2 - Vt);
            Jaco_ac(X2(3*X1_sorted(Esi(i,1),2),1), X2(3*X1_sorted(Esi(i,2),2),1)) = Jaco_ac(X2(3*X1_sorted(Esi(i,1),2),1), X2(3*X1_sorted(Esi(i,2),2),1))...
                - pcoeff*edgeEsi(i,1)/lenEsi(i,1)*(0.5*dphi_DC1 + Vt);
            Jaco_ac(X2(3*X1_sorted(Esi(i,1),2),1), X2(3*X1_sorted(Esi(i,3),2),1)) = Jaco_ac(X2(3*X1_sorted(Esi(i,1),2),1), X2(3*X1_sorted(Esi(i,3),2),1))...
                - pcoeff*edgeEsi(i,3)/lenEsi(i,3)*(0.5*dphi_DC2 + Vt);

            % phi part
            Jaco_ac(X2(3*X1_sorted(Esi(i,1),2),1), X2(3*X1_sorted(Esi(i,1),2)-2,1)) = Jaco_ac(X2(3*X1_sorted(Esi(i,1),2),1), X2(3*X1_sorted(Esi(i,1),2)-2,1))...
                - pcoeff*edgeEsi(i,1)/lenEsi(i,1)*(-p_DC_av1)...
                - pcoeff*edgeEsi(i,3)/lenEsi(i,3)*(-p_DC_av2);
            Jaco_ac(X2(3*X1_sorted(Esi(i,1),2),1), X2(3*X1_sorted(Esi(i,2),2)-2,1)) = Jaco_ac(X2(3*X1_sorted(Esi(i,1),2),1), X2(3*X1_sorted(Esi(i,2),2)-2,1))...
                - pcoeff*edgeEsi(i,1)/lenEsi(i,1)*p_DC_av1;
            Jaco_ac(X2(3*X1_sorted(Esi(i,1),2),1), X2(3*X1_sorted(Esi(i,3),2)-2,1)) = Jaco_ac(X2(3*X1_sorted(Esi(i,1),2),1), X2(3*X1_sorted(Esi(i,3),2)-2,1))...
                - pcoeff*edgeEsi(i,3)/lenEsi(i,3)*p_DC_av2;

            % p AC part
            Jaco_ac(X2(3*X1_sorted(Esi(i,1),2),1), X2(3*X1_sorted(Esi(i,1),2),1)) = Jaco_ac(X2(3*X1_sorted(Esi(i,1),2),1), X2(3*X1_sorted(Esi(i,1),2),1))...
                + q/4*(lenEsi(i,1)*edgeEsi(i,1) + lenEsi(i,3)*edge(i,3))*(2*pi*freq*iN);

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
            Jaco_ac(X2(3*X1_sorted(Esi(i,2),2),1), X2(3*X1_sorted(Esi(i,2),2),1)) = Jaco_ac(X2(3*X1_sorted(Esi(i,2),2),1), X2(3*X1_sorted(Esi(i,2),2),1))...
                - pcoeff*edgeEsi(i,2)/lenEsi(i,2)*(0.5*dphi_DC1 - Vt)...
                - pcoeff*edgeEsi(i,1)/lenEsi(i,1)*(0.5*dphi_DC2 - Vt);
            Jaco_ac(X2(3*X1_sorted(Esi(i,2),2),1), X2(3*X1_sorted(Esi(i,3),2),1)) = Jaco_ac(X2(3*X1_sorted(Esi(i,2),2),1), X2(3*X1_sorted(Esi(i,3),2),1))...
                - pcoeff*edgeEsi(i,2)/lenEsi(i,2)*(0.5*dphi_DC1 + Vt);
            Jaco_ac(X2(3*X1_sorted(Esi(i,2),2),1), X2(3*X1_sorted(Esi(i,1),2),1)) = Jaco_ac(X2(3*X1_sorted(Esi(i,2),2),1), X2(3*X1_sorted(Esi(i,1),2),1))...
                - pcoeff*edgeEsi(i,1)/lenEsi(i,1)*(0.5*dphi_DC2 + Vt);

            % phi part
            Jaco_ac(X2(3*X1_sorted(Esi(i,2),2),1), X2(3*X1_sorted(Esi(i,2),2)-2,1)) = Jaco_ac(X2(3*X1_sorted(Esi(i,2),2),1), X2(3*X1_sorted(Esi(i,2),2)-2,1))...
                - pcoeff*edgeEsi(i,2)/lenEsi(i,2)*(-p_DC_av1)...
                - pcoeff*edgeEsi(i,1)/lenEsi(i,1)*(-p_DC_av2);
            Jaco_ac(X2(3*X1_sorted(Esi(i,2),2),1), X2(3*X1_sorted(Esi(i,3),2)-2,1)) = Jaco_ac(X2(3*X1_sorted(Esi(i,2),2),1), X2(3*X1_sorted(Esi(i,3),2)-2,1))...
                - pcoeff*edgeEsi(i,2)/lenEsi(i,2)*p_DC_av1;
            Jaco_ac(X2(3*X1_sorted(Esi(i,2),2),1), X2(3*X1_sorted(Esi(i,1),2)-2,1)) = Jaco_ac(X2(3*X1_sorted(Esi(i,2),2),1), X2(3*X1_sorted(Esi(i,1),2)-2,1))...
                - pcoeff*edgeEsi(i,1)/lenEsi(i,1)*p_DC_av2;

            % p AC part
            Jaco_ac(X2(3*X1_sorted(Esi(i,2),2),1), X2(3*X1_sorted(Esi(i,2),2),1)) = Jaco_ac(X2(3*X1_sorted(Esi(i,2),2),1), X2(3*X1_sorted(Esi(i,2),2),1))...
                + q/4*(lenEsi(i,2)*edgeEsi(i,2) + lenEsi(i,1)*edge(i,1))*(2*pi*freq*iN);

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
            Jaco_ac(X2(3*X1_sorted(Esi(i,3),2),1), X2(3*X1_sorted(Esi(i,3),2),1)) = Jaco_ac(X2(3*X1_sorted(Esi(i,3),2),1), X2(3*X1_sorted(Esi(i,3),2),1))...
                - pcoeff*edgeEsi(i,3)/lenEsi(i,3)*(0.5*dphi_DC1 - Vt)...
                - pcoeff*edgeEsi(i,2)/lenEsi(i,2)*(0.5*dphi_DC2 - Vt);
            Jaco_ac(X2(3*X1_sorted(Esi(i,3),2),1), X2(3*X1_sorted(Esi(i,1),2),1)) = Jaco_ac(X2(3*X1_sorted(Esi(i,3),2),1), X2(3*X1_sorted(Esi(i,1),2),1))...
                - pcoeff*edgeEsi(i,3)/lenEsi(i,3)*(0.5*dphi_DC1 + Vt);
            Jaco_ac(X2(3*X1_sorted(Esi(i,3),2),1), X2(3*X1_sorted(Esi(i,2),2),1)) = Jaco_ac(X2(3*X1_sorted(Esi(i,3),2),1), X2(3*X1_sorted(Esi(i,2),2),1))...
                - pcoeff*edgeEsi(i,2)/lenEsi(i,2)*(0.5*dphi_DC2 + Vt);

            % phi part
            Jaco_ac(X2(3*X1_sorted(Esi(i,3),2),1), X2(3*X1_sorted(Esi(i,3),2)-2,1)) = Jaco_ac(X2(3*X1_sorted(Esi(i,3),2),1), X2(3*X1_sorted(Esi(i,3),2)-2,1))...
                - pcoeff*edgeEsi(i,3)/lenEsi(i,3)*(-p_DC_av1)...
                - pcoeff*edgeEsi(i,2)/lenEsi(i,2)*(-p_DC_av2);
            Jaco_ac(X2(3*X1_sorted(Esi(i,3),2),1), X2(3*X1_sorted(Esi(i,1),2)-2,1)) = Jaco_ac(X2(3*X1_sorted(Esi(i,3),2),1), X2(3*X1_sorted(Esi(i,1),2)-2,1))...
                - pcoeff*edgeEsi(i,3)/lenEsi(i,3)*p_DC_av1;
            Jaco_ac(X2(3*X1_sorted(Esi(i,3),2),1), X2(3*X1_sorted(Esi(i,2),2)-2,1)) = Jaco_ac(X2(3*X1_sorted(Esi(i,3),2),1), X2(3*X1_sorted(Esi(i,2),2)-2,1))...
                - pcoeff*edgeEsi(i,2)/lenEsi(i,2)*p_DC_av2;

            % p AC part
            Jaco_ac(X2(3*X1_sorted(Esi(i,3),2),1), X2(3*X1_sorted(Esi(i,3),2),1)) = Jaco_ac(X2(3*X1_sorted(Esi(i,3),2),1), X2(3*X1_sorted(Esi(i,3),2),1))...
                + q/4*(lenEsi(i,3)*edgeEsi(i,3) + lenEsi(i,2)*edge(i,2))*(2*pi*freq*iN);
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Dirichlet boundary condition for final iteration
for i=1:(Ncon1+Ncon2)
    for j=1:concol
        fcon = X1_sorted(contact(i,j),2);
        fcon_np = find(contact(i,j)==V_silicon(:,1));

        Jaco_ac(X2(3*fcon-2),:) = 0;
        Jaco_ac(X2(3*fcon-1),:) = 0;
        Jaco_ac(X2(3*fcon),:) = 0;

        Jaco_ac(X2(3*fcon-2), X2(3*fcon-2)) = 1;
        Jaco_ac(X2(3*fcon-1), X2(3*fcon-1)) = 1;
        Jaco_ac(X2(3*fcon), X2(3*fcon)) = 1;

        if i<=Ncon1 % Cathode
            Res_ac(X2(3*fcon-2),1) = dphi(fcon)- 0;
            Res_ac(X2(3*fcon-1),1) = dn(fcon_np) - 0;
            Res_ac(X2(3*fcon),1) = dp(fcon_np) - 0;

        elseif i>(Ncon1) &&  i<=(Ncon1+Ncon2) % Anode
            Res_ac(X2(3*fcon-2),1) = dphi(fcon) - Vamp*iN;
            Res_ac(X2(3*fcon-1),1) = dn(fcon_np) - 0;
            Res_ac(X2(3*fcon),1) = dp(fcon_np) - 0;
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% scaled part
Cvector = zeros(size(Jaco_ac,1),1);
for i=1:size(Jaco_ac,1)
    if rem(X2_sorted(i),3) == 1
        Cvector(i) = Vt;
    elseif rem(X2_sorted(i),3) == 2
        Cvector(i) = Nd;
    elseif rem(X2_sorted(i),3) == 0
        Cvector(i) = Nd;
    end
end

Cmatrix = spdiags(Cvector,0,size(X2_sorted,1),size(X2_sorted,1));
Jaco_scaled = Jaco_ac * Cmatrix;
Rvector = 1./sum(abs(Jaco_scaled),2);
Rmatrix = spdiags(Rvector,0,size(X2_sorted,1),size(X2_sorted,1));

Jaco_scaled = Rmatrix * Jaco_scaled;
Res_scaled = Rmatrix * Res_ac;

update_scaled = Jaco_scaled \ (-Res_scaled);
update = Cmatrix * update_scaled;

for ii=1:VErow
    dphi(ii,1) = dphi(ii,1) + update(3*ii-2);
    dn(ii,1) = dn(ii,1) + update(3*ii-1);
    dp(ii,1) = dp(ii,1) + update(3*ii);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% current calculation

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

            In(1,1) = In(1,1) ...
                - ncoeff*edgeE4(j,1)/lenE4(j,1)*(n_DC_av1*ddphi1 + dn_av1*dphi_DC1 - ddn1*Vt)*width...
                - ncoeff*edgeE4(j,3)/lenE4(j,3)*(n_DC_av2*ddphi2 + dn_av2*dphi_DC2 - ddn2*Vt)*width;
            Ip(1,1) = Ip(1,1) ...
                - pcoeff*edgeE4(j,1)/lenE4(j,1)*(p_DC_av1*ddphi1 + dp_av1*dphi_DC1 + ddp1*Vt)*width...
                - pcoeff*edgeE4(j,3)/lenE4(j,3)*(p_DC_av2*ddphi2 + dp_av2*dphi_DC2 + ddp2*Vt)*width;
            Id(1,1) = Id(1,1) ...
                - esi*e0*edgeE4(j,1)/lenE4(j,1)*(2*pi*freq*iN)*ddphi1*width...
                - esi*e0*edgeE4(j,3)/lenE4(j,3)*(2*pi*freq*iN)*ddphi2*width;

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

            In(1,1) = In(1,1) ...
                - ncoeff*edgeE4(j,2)/lenE4(j,2)*(n_DC_av1*ddphi1 + dn_av1*dphi_DC1 - ddn1*Vt)*width...
                - ncoeff*edgeE4(j,1)/lenE4(j,1)*(n_DC_av2*ddphi2 + dn_av2*dphi_DC2 - ddn2*Vt)*width;
            Ip(1,1) = Ip(1,1) ...
                - pcoeff*edgeE4(j,2)/lenE4(j,2)*(p_DC_av1*ddphi1 + dp_av1*dphi_DC1 + ddp1*Vt)*width...
                - pcoeff*edgeE4(j,1)/lenE4(j,1)*(p_DC_av2*ddphi2 + dp_av2*dphi_DC2 + ddp2*Vt)*width;
            Id(1,1) = Id(1,1) ...
                - esi*e0*edgeE4(j,2)/lenE4(j,2)*(2*pi*freq*iN)*ddphi1*width...
                - esi*e0*edgeE4(j,1)/lenE4(j,1)*(2*pi*freq*iN)*ddphi2*width;

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

            In(1,1) = In(1,1) ...
                - ncoeff*edgeE4(j,3)/lenE4(j,3)*(n_DC_av1*ddphi1 + dn_av1*dphi_DC1 - ddn1*Vt)*width...
                - ncoeff*edgeE4(j,2)/lenE4(j,2)*(n_DC_av2*ddphi2 + dn_av2*dphi_DC2 - ddn2*Vt)*width;
            Ip(1,1) = Ip(1,1) ...
                - pcoeff*edgeE4(j,3)/lenE4(j,3)*(p_DC_av1*ddphi1 + dp_av1*dphi_DC1 + ddp1*Vt)*width...
                - pcoeff*edgeE4(j,2)/lenE4(j,2)*(p_DC_av2*ddphi2 + dp_av2*dphi_DC2 + ddp2*Vt)*width;
            Id(1,1) = Id(1,1) ...
                - esi*e0*edgeE4(j,3)/lenE4(j,3)*(2*pi*freq*iN)*ddphi1*width...
                - esi*e0*edgeE4(j,2)/lenE4(j,2)*(2*pi*freq*iN)*ddphi2*width;
        end
    end
end

I = In + Id;

ansmat = zeros(3,1);
ansmat(1,1) = max(I);
ansmat(2,1) = max(In);
ansmat(3,1) = max(Id);
