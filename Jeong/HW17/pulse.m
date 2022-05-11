function [Voltage_save, time_save] = pulse(T, Time_step, Amp, Duty_cycle, Slope, Period)
% pulse만들어주는 함수를 작성했다.
% 사용변수는 주기 T, 주기를 나누어줄 간격인 Time_step, pulse의 최대 전압인 Max_voltage,
% pulse가 유지되는 시간인 duty_cycle, pulse의 기울기인 slope,
% 주기를 몇번 반복할 것인가인 period가 있다.
%
% 변수의 단위는 다음과 같다.
% T[s], Time_step[비율], Amp[V], Duty_cycle[%], Slope[%], Period[횟수]
%
% 기본으로 1주기당 56개의 노드로 나누어 주었다. Time_step을 늘리면 56*Time_step으로 노드가 늘어난다.
%
% 사용 예시는 다음과 같다.
% [Voltage_save, time_save] = pulse(T, Time_step, Amp, Duty_cycle, Slope, Period)
% 출력은 전압과, 시간이 출력된다.
% T=0.1ns, Amp=1V일때, Duty_cycle=10%, Slope=30%, Time_step=2000까지 그래프가 정상적으로 나옴을 확인.
%
% EX)
%
% T=0.1e-9; % [ns], 주기
% Time_step=1; %
% Amp=1; % [V] % Amp 값
% Duty_cycle=10; % [%]
% Slope=30; % [%]
% Period=2;
% [Vg,time]=pulse(T, Time_step, Amp, Duty_cycle, Slope, Period);


%% 중간에 안줄인 버전
% Duty_cycle=Duty_cycle/100; % Duty_cycle [%], 단위를 %로 변경함.
% Slope=Slope/100; % Slope [%], 단위를 %로 변경함.
%
% % 변수
% delta_t=T/Time_step; % 주기/나눌 간격 하면 delta_t 나옴.
% time_save=transpose(0:delta_t:T*Period);
%
%     count=1; count_delet1=1; count_delet2=1;
% for P=1:Period
%
%
%     clearvars time
%     if P==1
%         time=0:delta_t:T;
%     else
%         time=(P-1)*T+delta_t:delta_t:P*T;
%     end
%
%     for TIME=1:length(time) %%% 0~5us까지 관찰
%
%         t=time(1,TIME);
%
%         a=0+(P-1)*T; b=T*Duty_cycle*Slope+(P-1)*T; c=T*Duty_cycle*(1-Slope)+(P-1)*T; d=T*Duty_cycle+(P-1)*T; e=T+(P-1)*T;
%
%         % duty cycle이 slope의 중간지점에서 측정한것을 반영하여 c와 d를 다시 계산함.
%         duty_tmp1=(a+b)/2; duty_tmp2=(c+d)/2;
%         duty_tmp3=duty_tmp2-duty_tmp1;
%         duty_tmp4=T*Duty_cycle-duty_tmp3;
%         c=c+duty_tmp4; d=d+duty_tmp4;
%
%         if t>=a && t<b % region 1
%             Voltage=Amp/(T*Duty_cycle*Slope)*(t-(P-1)*T);
%
%         elseif t>=b && t<c % region 2
%             Voltage=Amp;
%
%         elseif t>=c && t<d % region 3
%             Voltage=-Amp/(T*Duty_cycle*Slope)*(t-(P-1)*T-duty_tmp4)+Amp/(T*Duty_cycle*Slope)*(T*Duty_cycle);
%
%         elseif t>=d && t<=e % region 4
%             Voltage=0;
%         end
%
%         Voltage_save(count,1)=Voltage;
%         count=count+1;
%     end
% end
%
% % rmmissing_nan
% tmp=[time_save Voltage_save];
% tmp=rmmissing(tmp);
%
% clearvars time_save Voltage_save
% time_save(:,1)=tmp(:,1);
% Voltage_save(:,1)=tmp(:,2);

%% 중간에 노드를 감소시킨 버전
% Duty_cycle=Duty_cycle/100; % Duty_cycle [%], 단위를 %로 변경함.
% Slope=Slope/100; % Slope [%], 단위를 %로 변경함.
%
% % 변수
% delta_t=T/Time_step; % 주기/나눌 간격 하면 delta_t 나옴.
% time_save=transpose(0:delta_t:T*Period);
% del1=40;
% del2=50;
% pls1=0; % 0~1사이의 값 주어야함. off되고 얼마나 더 관측할 것인지에 대한 값.
%
%     count=1; count_delet1=1; count_delet2=1;
% for P=1:Period
%
%
%     clearvars time
%     if P==1
%         time=0:delta_t:T;
%     else
%         time=(P-1)*T+delta_t:delta_t:P*T;
%     end
%
%     for TIME=1:length(time) %%% 0~5us까지 관찰
%
%         t=time(1,TIME);
%
%         a=0+(P-1)*T; b=T*Duty_cycle*Slope+(P-1)*T; c=T*Duty_cycle*(1-Slope)+(P-1)*T; d=T*Duty_cycle+(P-1)*T; e=T+(P-1)*T;
%
%         % duty cycle이 slope의 중간지점에서 측정한것을 반영하여 c와 d를 다시 계산함.
%         duty_tmp1=(a+b)/2; duty_tmp2=(c+d)/2;
%         duty_tmp3=duty_tmp2-duty_tmp1;
%         duty_tmp4=T*Duty_cycle-duty_tmp3;
%         c=c+duty_tmp4; d=d+duty_tmp4;
%
%        norm=1e-12;
%         if t>=a && t<b % region 1
%             Voltage=Amp/(T*Duty_cycle*Slope)*(t-(P-1)*T);
%
%         elseif t>=b && t<c % region 2
%             if t<=b+(c-b)*pls1 || t>=c-(c-b)*pls1
%                 Voltage=Amp;
%             else
%                 if count_delet1==1 || abs(t-(b+(c-b)*pls1))<=norm || abs(t-(c-(c-b)*pls1))<=norm
%                     Voltage=Amp;
%                     count_delet1=count_delet1+1;
%                 elseif count_delet1>=del1
%                     Voltage=nan;
%                     count_delet1=1;
%                 elseif t>=b+(c-b)*pls1 && t<=c-(c-b)*pls1 && count_delet1~=1
%                     Voltage=nan;
%                     count_delet1=count_delet1+1;
%                 end
%             end
%
%         elseif t>=c && t<d % region 3
%             Voltage=-Amp/(T*Duty_cycle*Slope)*(t-(P-1)*T-duty_tmp4)+Amp/(T*Duty_cycle*Slope)*(T*Duty_cycle);
%
%         elseif t>=d && t<=e % region 4
%             if t<=d+(d-c)*pls1 || t>=e-(d-c)*pls1
%                 Voltage=0;
%             else
%                 if count_delet2==1 || abs(t-(d+(d-c)*pls1))<=norm || abs(t-(e-(d-c)*pls1))<=norm
%                     Voltage=0;
%                     count_delet2=count_delet2+1;
%                 elseif count_delet2>=del2
%                     Voltage=nan;
%                     count_delet2=1;
%                 elseif t>=d+(d-c)*pls1 && t<=e-(d-c)*pls1 && count_delet2~=1
%                     Voltage=nan;
%                     count_delet2=count_delet2+1;
%                 end
%             end
%         end
%
%         Voltage_save(count,1)=Voltage;
%         count=count+1;
%     end
% end
%
% % rmmissing_nan
% tmp=[time_save Voltage_save];
% tmp=rmmissing(tmp);
%
% clearvars time_save Voltage_save
% time_save(:,1)=tmp(:,1);
% Voltage_save(:,1)=tmp(:,2);

%% Slope부분에 노드를 증가시킨 경우

Duty_cycle=Duty_cycle/100; % Duty_cycle [%], 단위를 %로 변경함.
Slope=Slope/100; % Slope [%], 단위를 %로 변경함.

a=0; b=T*Duty_cycle*Slope; c=T*Duty_cycle*(1-Slope); d=T*Duty_cycle; e=T;

% duty cycle이 slope의 중간지점에서 측정한것을 반영하여 c와 d를 다시 계산함.
duty_tmp1=(a+b)/2; duty_tmp2=(c+d)/2;
duty_tmp3=duty_tmp2-duty_tmp1;
duty_tmp4=T*Duty_cycle-duty_tmp3;
c=c+duty_tmp4; d=d+duty_tmp4;

t1=transpose(linspace(a,b,(20*Time_step+1)));
t2=transpose(b:(c-b)/(5*Time_step):c);
t3=transpose(linspace(c,d,(20*Time_step+1)));
t4=transpose(d:(e-d)/(10*Time_step):e-(e-d)/(10*Time_step));

time_tmp=[t1; t2; t3; t4];

time_tmp=unique(time_tmp);

count=1;
for TIME=1:length(time_tmp)

    t=time_tmp(TIME,1);

    if t>=a && t<b % region 1
        Voltage=Amp/(T*Duty_cycle*Slope)*t;

    elseif t>=b && t<c % region 2
        Voltage=Amp;

    elseif t>=c && t<d % region 3
        Voltage=-Amp/(T*Duty_cycle*Slope)*(t-d);

    elseif t>=d && t<=e % region 4
        Voltage=0;
    end

    Voltage_save_tmp(count,1)=Voltage;
    count=count+1;
end

time=[]; Voltage_save=[];
for P=1:Period
    time=[time; time_tmp+T*(P-1)];

    Voltage_save=[Voltage_save; Voltage_save_tmp];
end
time_save=[time; T*P]; Voltage_save=[Voltage_save; 0];

% %% Slope기울기 고정 (1V / 2.5ns)
% 
% Duty_cycle=Duty_cycle/100; % Duty_cycle [%], 단위를 %로 변경함.
% Slope=Slope/100; % Slope [%], 단위를 %로 변경함.
% 
% a=0; b=2.5e-9; c=T*Duty_cycle-2.5e-9; d=T*Duty_cycle; e=T;
% 
% % duty cycle이 slope의 중간지점에서 측정한것을 반영하여 c와 d를 다시 계산함.
% duty_tmp1=(a+b)/2; duty_tmp2=(c+d)/2;
% duty_tmp3=duty_tmp2-duty_tmp1;
% duty_tmp4=T*Duty_cycle-duty_tmp3;
% c=c+duty_tmp4; d=d+duty_tmp4;
% 
% t1=transpose(linspace(a,b,(20*Time_step+1)));
% t2=transpose(b:(c-b)/(5*Time_step):c);
% t3=transpose(linspace(c,d,(20*Time_step+1)));
% t4=transpose(d:(e-d)/(10*Time_step):e-(e-d)/(10*Time_step));
% 
% time_tmp=[t1; t2; t3; t4];
% 
% time_tmp=unique(time_tmp);
% 
% count=1;
% for TIME=1:length(time_tmp)
% 
%     t=time_tmp(TIME,1);
% 
%     if t>=a && t<b % region 1
%         Voltage=(1/2.5e-9)*t;
% 
%     elseif t>=b && t<c % region 2
%         Voltage=Amp;
% 
%     elseif t>=c && t<d % region 3
%         Voltage=-(1/2.5e-9)*(t-d);
% 
%     elseif t>=d && t<=e % region 4
%         Voltage=0;
%     end
% 
%     Voltage_save_tmp(count,1)=Voltage;
%     count=count+1;
% end
% 
% time=[]; Voltage_save=[];
% for P=1:Period
%     time=[time; time_tmp+T*(P-1)];
% 
%     Voltage_save=[Voltage_save; Voltage_save_tmp];
% end
% time_save=[time; T*P]; Voltage_save=[Voltage_save; 0];