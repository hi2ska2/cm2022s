time_x = 1e-12;

analytic_Jd = zeros(20,1);

for ii=1:100

    Time = time_x/(2*ii); % Total time

    freq(ii,1)=1/(Time);
    
    analytic_Jd(ii,1) = (2*pi*freq(ii,1)*esi*e0*((5e-9*1e-6)/120e-9));
end

analytic_Jn = q*m_n*2e23*((5e-9*1e-6)/120e-9);