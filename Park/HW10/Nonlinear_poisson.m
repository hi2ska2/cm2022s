%%%%% Nonlinear Poisson %%%%%%

for Newton=1:10

    Jaco=zeros(index, index);
    res=zeros(index,1);

    for ii=1:R_row1
        for jj=1:3

            if jj==1
                res(E_R1(ii,jj),1) = res(E_R1(ii,jj),1) + eox*(phi(E_R1(ii,jj+1),1)-phi(E_R1(ii,jj),1))*(edge(ii,jj)/L(ii,jj)) + eox*(phi(E_R1(ii,jj+2),1)-phi(E_R1(ii,jj),1))*(edge(ii,jj+2)/L(ii,jj+2));
                Jaco(E_R1(ii,jj), E_R1(ii,jj+2)) = Jaco(E_R1(ii,jj), E_R1(ii,jj+2))+ eox*(edge(ii,jj+2)/L(ii,jj+2));
                Jaco(E_R1(ii,jj), E_R1(ii,jj+1)) = Jaco(E_R1(ii,jj), E_R1(ii,jj+1))+eox*(edge(ii,jj)/L(ii,jj));
                Jaco(E_R1(ii,jj), E_R1(ii,jj)) = Jaco(E_R1(ii,jj), E_R1(ii,jj))-eox*(edge(ii,jj+2)/L(ii,jj+2))-eox*(edge(ii,jj)/L(ii,jj));

            elseif jj ==2
                res(E_R1(ii,jj),1) = res(E_R1(ii,jj),1) + eox*(phi(E_R1(ii,jj-1),1)-phi(E_R1(ii,jj),1))*(edge(ii,jj-1)/L(ii,jj-1)) + eox*(phi(E_R1(ii,jj+1),1)-phi(E_R1(ii,jj),1))*(edge(ii,jj)/L(ii,jj));
                Jaco(E_R1(ii,jj), E_R1(ii,jj+1)) = Jaco(E_R1(ii,jj), E_R1(ii,jj+1))+eox*(edge(ii,jj)/L(ii,jj));
                Jaco(E_R1(ii,jj), E_R1(ii,jj)) =  Jaco(E_R1(ii,jj), E_R1(ii,jj))-eox*(edge(ii,jj)/L(ii,jj))-eox*(edge(ii,jj-1)/L(ii,jj-1));
                Jaco(E_R1(ii,jj), E_R1(ii,jj-1)) = Jaco(E_R1(ii,jj), E_R1(ii,jj-1))+eox*(edge(ii,jj-1)/L(ii,jj-1));

            elseif jj == 3
                res(E_R1(ii,jj),1) = res(E_R1(ii,jj),1) + eox*(phi(E_R1(ii,jj-2),1)-phi(E_R1(ii,jj),1))*(edge(ii,jj)/L(ii,jj)) + eox*(phi(E_R1(ii,jj-1),1)-phi(E_R1(ii,jj),1))*(edge(ii,jj-1)/L(ii,jj-1));
                Jaco(E_R1(ii,jj), E_R1(ii,jj-2)) = Jaco(E_R1(ii,jj), E_R1(ii,jj-2))+eox*(edge(ii,jj)/L(ii,jj));
                Jaco(E_R1(ii,jj), E_R1(ii,jj-1)) = Jaco(E_R1(ii,jj), E_R1(ii,jj-1))+eox*(edge(ii,jj-1)/L(ii,jj-1));
                Jaco(E_R1(ii,jj), E_R1(ii,jj)) = Jaco(E_R1(ii,jj), E_R1(ii,jj))-eox*(edge(ii,jj)/L(ii,jj))-eox*(edge(ii,jj-1)/L(ii,jj-1));
            end
        end
    end

    for ii=1:R_row2
        for jj=1:3
            if jj==1
                res(E_R2(ii,jj),1) = res(E_R2(ii,jj),1) + esi*(phi(E_R2(ii,jj+1),1)-phi(E_R2(ii,jj),1))*(edge(ii+R_row1,jj)/L(ii+R_row1,jj)) + esi*(phi(E_R2(ii,jj+2),1)-phi(E_R2(ii,jj),1))*(edge(ii+R_row1,jj+2)/L(ii+R_row1,jj+2));
                res(E_R2(ii,jj),1) = res(E_R2(ii,jj),1) + coeff*(edge(ii+R_row1,jj)*L(ii+R_row1,jj)+edge(ii+R_row1,jj+2)*L(ii+R_row1,jj+2))*(-Na-nint*exp(phi(E_R2(ii,jj),1)/Thermal_V)+nint*exp(-phi(E_R2(ii,jj),1)/Thermal_V));
                Jaco(E_R2(ii,jj), E_R2(ii,jj+2)) = Jaco(E_R2(ii,jj), E_R2(ii,jj+2))+ esi*(edge(ii+R_row1,jj+2)/L(ii+R_row1,jj+2));
                Jaco(E_R2(ii,jj), E_R2(ii,jj+1)) = Jaco(E_R2(ii,jj), E_R2(ii,jj+1))+esi*(edge(ii+R_row1,jj)/L(ii+R_row1,jj));
                Jaco(E_R2(ii,jj), E_R2(ii,jj)) = Jaco(E_R2(ii,jj), E_R2(ii,jj))-esi*(edge(ii+R_row1,jj+2)/L(ii+R_row1,jj+2))-esi*(edge(ii+R_row1,jj)/L(ii+R_row1,jj));
                Jaco(E_R2(ii,jj), E_R2(ii,jj)) = Jaco(E_R2(ii,jj), E_R2(ii,jj))+coeff*(edge(ii+R_row1,jj)*L(ii+R_row1,jj)+edge(ii+R_row1,jj+2)*L(ii+R_row1,jj+2))*((-nint/Thermal_V)*exp(phi(E_R2(ii,jj),1)/Thermal_V)+(-nint/Thermal_V)*exp(-phi(E_R2(ii,jj),1)/Thermal_V));

            elseif jj ==2
                res(E_R2(ii,jj),1) = res(E_R2(ii,jj),1) + esi*(phi(E_R2(ii,jj-1),1)-phi(E_R2(ii,jj),1))*(edge(ii+R_row1,jj-1)/L(ii+R_row1,jj-1)) + esi*(phi(E_R2(ii,jj+1),1)-phi(E_R2(ii,jj),1))*(edge(ii+R_row1,jj)/L(ii+R_row1,jj));
                res(E_R2(ii,jj),1) = res(E_R2(ii,jj),1) + coeff*(edge(ii+R_row1,jj)*L(ii+R_row1,jj)+edge(ii+R_row1,jj-1)*L(ii+R_row1,jj-1))*(-Na-nint*exp(phi(E_R2(ii,jj),1)/Thermal_V)+nint*exp(-phi(E_R2(ii,jj),1)/Thermal_V));
                Jaco(E_R2(ii,jj), E_R2(ii,jj+1)) = Jaco(E_R2(ii,jj), E_R2(ii,jj+1))+esi*(edge(ii+R_row1,jj)/L(ii+R_row1,jj));
                Jaco(E_R2(ii,jj), E_R2(ii,jj-1)) = Jaco(E_R2(ii,jj), E_R2(ii,jj-1))+esi*(edge(ii+R_row1,jj-1)/L(ii+R_row1,jj-1));
                Jaco(E_R2(ii,jj), E_R2(ii,jj)) =  Jaco(E_R2(ii,jj), E_R2(ii,jj))-esi*(edge(ii+R_row1,jj)/L(ii+R_row1,jj))-esi*(edge(ii+R_row1,jj-1)/L(ii+R_row1,jj-1));
                Jaco(E_R2(ii,jj), E_R2(ii,jj)) =  Jaco(E_R2(ii,jj), E_R2(ii,jj))+coeff*(edge(ii+R_row1,jj)*L(ii+R_row1,jj)+edge(ii+R_row1,jj-1)*L(ii+R_row1,jj-1))*((-nint/Thermal_V)*exp(phi(E_R2(ii,jj),1)/Thermal_V)+(-nint/Thermal_V)*exp(-phi(E_R2(ii,jj),1)/Thermal_V));

            elseif jj == 3
                res(E_R2(ii,jj),1) = res(E_R2(ii,jj),1) + esi*(phi(E_R2(ii,jj-2),1)-phi(E_R2(ii,jj),1))*(edge(ii+R_row1,jj)/L(ii+R_row1,jj)) + esi*(phi(E_R2(ii,jj-1),1)-phi(E_R2(ii,jj),1))*(edge(ii+R_row1,jj-1)/L(ii+R_row1,jj-1));
                res(E_R2(ii,jj),1) = res(E_R2(ii,jj),1) + coeff*(edge(ii+R_row1,jj)*L(ii+R_row1,jj)+edge(ii+R_row1,jj-1)*L(ii+R_row1,jj-1))*(-Na-nint*exp(phi(E_R2(ii,jj),1)/Thermal_V)+nint*exp(-phi(E_R2(ii,jj),1)/Thermal_V));
                Jaco(E_R2(ii,jj), E_R2(ii,jj-2)) = Jaco(E_R2(ii,jj), E_R2(ii,jj-2))+esi*(edge(ii+R_row1,jj)/L(ii+R_row1,jj));
                Jaco(E_R2(ii,jj), E_R2(ii,jj-1)) = Jaco(E_R2(ii,jj), E_R2(ii,jj-1))+esi*(edge(ii+R_row1,jj-1)/L(ii+R_row1,jj-1));
                Jaco(E_R2(ii,jj), E_R2(ii,jj)) = Jaco(E_R2(ii,jj), E_R2(ii,jj))-esi*(edge(ii+R_row1,jj)/L(ii+R_row1,jj))-esi*(edge(ii+R_row1,jj-1)/L(ii+R_row1,jj-1));
                Jaco(E_R2(ii,jj), E_R2(ii,jj)) = Jaco(E_R2(ii,jj), E_R2(ii,jj))+coeff*(edge(ii+R_row1,jj)*L(ii+R_row1,jj)+edge(ii+R_row1,jj-1)*L(ii+R_row1,jj-1))*((-nint/Thermal_V)*exp(phi(E_R2(ii,jj),1)/Thermal_V)+(-nint/Thermal_V)*exp(-phi(E_R2(ii,jj),1)/Thermal_V));
            end
        end
    end

    for ii=1:R_row3
        for jj=1:3

            if jj==1
                res(E_R3(ii,jj),1) = res(E_R3(ii,jj),1) + eox*(phi(E_R3(ii,jj+1),1)-phi(E_R3(ii,jj),1))*(edge(ii+R_row2+R_row1,jj)/L(ii+R_row2+R_row1,jj)) + eox*(phi(E_R3(ii,jj+2),1)-phi(E_R3(ii,jj),1))*(edge(ii+R_row2+R_row1,jj+2)/L(ii+R_row2+R_row1,jj+2));
                Jaco(E_R3(ii,jj), E_R3(ii,jj+2)) = Jaco(E_R3(ii,jj), E_R3(ii,jj+2))+ eox*(edge(ii+R_row2+R_row1,jj+2)/L(ii+R_row2+R_row1,jj+2));
                Jaco(E_R3(ii,jj), E_R3(ii,jj+1)) = Jaco(E_R3(ii,jj), E_R3(ii,jj+1))+eox*(edge(ii+R_row2+R_row1,jj)/L(ii+R_row2+R_row1,jj));
                Jaco(E_R3(ii,jj), E_R3(ii,jj)) = Jaco(E_R3(ii,jj), E_R3(ii,jj))-eox*(edge(ii+R_row2+R_row1,jj+2)/L(ii+R_row2+R_row1,jj+2))-eox*(edge(ii+R_row2+R_row1,jj)/L(ii+R_row2+R_row1,jj));

            elseif jj ==2
                res(E_R3(ii,jj),1) = res(E_R3(ii,jj),1) + eox*(phi(E_R3(ii,jj-1),1)-phi(E_R3(ii,jj),1))*(edge(ii+R_row2+R_row1,jj-1)/L(ii+R_row2+R_row1,jj-1)) + eox*(phi(E_R3(ii,jj+1),1)-phi(E_R3(ii,jj),1))*(edge(ii+R_row2+R_row1,jj)/L(ii+R_row2+R_row1,jj));
                Jaco(E_R3(ii,jj), E_R3(ii,jj+1)) = Jaco(E_R3(ii,jj), E_R3(ii,jj+1))+eox*(edge(ii+R_row2+R_row1,jj)/L(ii+R_row2+R_row1,jj));
                Jaco(E_R3(ii,jj), E_R3(ii,jj)) =  Jaco(E_R3(ii,jj), E_R3(ii,jj))-eox*(edge(ii+R_row2+R_row1,jj)/L(ii+R_row2+R_row1,jj))-eox*(edge(ii+R_row2+R_row1,jj-1)/L(ii+R_row2+R_row1,jj-1));
                Jaco(E_R3(ii,jj), E_R3(ii,jj-1)) = Jaco(E_R3(ii,jj), E_R3(ii,jj-1))+eox*(edge(ii+R_row2+R_row1,jj-1)/L(ii+R_row2+R_row1,jj-1));

            elseif jj == 3
                res(E_R3(ii,jj),1) = res(E_R3(ii,jj),1) + eox*(phi(E_R3(ii,jj-2),1)-phi(E_R3(ii,jj),1))*(edge(ii+R_row2+R_row1,jj)/L(ii+R_row2+R_row1,jj)) + eox*(phi(E_R3(ii,jj-1),1)-phi(E_R3(ii,jj),1))*(edge(ii+R_row2+R_row1,jj-1)/L(ii+R_row2+R_row1,jj-1));
                Jaco(E_R3(ii,jj), E_R3(ii,jj-2)) = Jaco(E_R3(ii,jj), E_R3(ii,jj-2))+eox*(edge(ii+R_row2+R_row1,jj)/L(ii+R_row2+R_row1,jj));
                Jaco(E_R3(ii,jj), E_R3(ii,jj-1)) = Jaco(E_R3(ii,jj), E_R3(ii,jj-1))+eox*(edge(ii+R_row2+R_row1,jj-1)/L(ii+R_row2+R_row1,jj-1));
                Jaco(E_R3(ii,jj), E_R3(ii,jj)) = Jaco(E_R3(ii,jj), E_R3(ii,jj))-eox*(edge(ii+R_row2+R_row1,jj)/L(ii+R_row2+R_row1,jj))-eox*(edge(ii+R_row2+R_row1,jj-1)/L(ii+R_row2+R_row1,jj-1));
            end
        end
    end

    %%%%% Dirichlet Boundary Condition %%%%%
    for n=1:C_col
        for i =1: C_row
            if i <= 2
                Jaco(Contact(i,n),:) =0;
                Jaco(Contact(i,n),Contact(i,n)) = 1;
                res(Contact(i,n),1)=phi(Contact(i,n),1)-0;

            elseif i > 2 && i <= 4
                Jaco(Contact(i,n),:) =0;
                Jaco(Contact(i,n),Contact(i,n)) = 1;
                res(Contact(i,n),1)= phi(Contact(i,n),1)-1;

            end
        end
    end

    delphi(:,Newton) = Jaco \ (-res);
    phi = phi+delphi(:,Newton);
end

elec=zeros(index,1);
hole=zeros(index,1);

elec(7:1:21,1) = 1e-6*nint*exp(phi(7:1:21,1)/Thermal_V);
hole(7:1:21,1) = 1e-6*nint*exp(-phi(7:1:21,1)/Thermal_V);

