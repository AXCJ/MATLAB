clc;clear all;close all;pause(0.01)


A=[0 1;-3 -2];
B=[0; 1];
C=[1 1];

Ts=0.18;
Tf=0.06;
Tw=0.01;
T_Sim=5.4;
Td_seq=0:Tf:T_Sim; % period = 0.06
Tw_seq=0:Tw:T_Sim; % period = 0.01
% white_noise=wgn(1,length(Tw_seq),1);
load white_noise.mat;
x0=[0;0];
white=white_noise(1:6:end);
plot(Td_seq,white,'o','LineWidth',3)
hold on
stairs(Td_seq,white,'b')
white=[0 white];
stairs([Td_seq 5.46],white,':','LineWidth',1.5,'color','r')
legend('u_{t}','u_{d}','u_{t-\tau}')
axis([0 5.5 -4 3])
xlabel('Time(s)','fontsize',10,'fontweight','b');
title('input','fontsize',10,'fontweight','b');
hold off

% sim
% simDelay=0.12;
simDelay=delay_time;
simin=[Tw_seq' white_noise'];

%
[n,m]=size(B);
[p,n]=size(C);
D=zeros(p,m);

%
[G,H]=c2d(A,B,Tf);
Gamma_0=(real(G^(1-simDelay/Ts))-eye(n))*inv(A)*B;
Gamma_1=(G-real(G^(1-simDelay/Ts)))*inv(A)*B;
PI_ok=real(G^(-2/3))*H;

%
[Tc,Xc,Yc]=sim('CTD_system_Based_TF2',Td_seq);
figure
plot(Tc,Yc,'LineWidth',2)


Yc=Yc';
Tc=Tc';
Td_Num=length(Td_seq)-1;
kk=0;
for ii=1:Td_Num
    if mod(ii-1,(Ts/Tf))==0
        kk=kk+1;
        Yd(:,kk) = Yc(:,ii);
        Td(:,kk) = Tc(:,ii);
    end
end


ii=0;kk=0;
temp_T=Tf;
for t=0:Tf:T_Sim
    ii=ii+1;

    if mod(ii-1,(temp_T/Tf))==0
        kk=kk+1;
        if kk==1
            xf(:,kk)=pinv(C)*(C*x0);
        elseif kk==2
            xf(:,kk)=G*xf(:,kk-1)+Gamma_0*Us(:,kk-1);
        else
            xf(:,kk)=G*xf(:,kk-1)+Gamma_0*Us(:,kk-1)+Gamma_1*Us(:,kk-2);
        end
        Us(:,kk) = white_noise(:,6*(kk-1)+1);
    end

    if mod(ii-1,(temp_T/Tf))==0
        Yf(:,kk)=C*xf(:,kk);
        th(:,kk)=(kk-1)*temp_T;
    end
    
%     disp(sprintf('Simulation Time:%2.3f sec',(ii-1)*Tf));
end

Yf=[zeros(1,3) Yf(1:88)];
for kk=1:size(Yc,2)
    e(:,kk)=Yc(:,kk)-Yf(:,kk);
    if abs(Yc(:,kk))==0
%         Yc(:,kk)=10^-6;
        error(:,kk)=(abs(e(:,kk))./abs(Yc(:,kk)))*100;
    else
        error(:,kk)=(abs(e(:,kk))./abs(Yc(:,kk)))*100;
    end
end

% ode
load ud_c
xdot=@(t,x)  A*x+B*ud_c(2,floor(t/Tf)+1);
[t,x]=ode45(xdot,Td_seq,[0 0]);
Y_ode=C(1)*x(:,1)+C(2)*x(:,2);



hold on
plot(Td,Yd,'o')
plot(th,Yf,'LineWidth',2)
plot(th,Ytp,'LineWidth',2)
plot(th,Yzz,'LineWidth',2)
plot(th,Y_ode,'LineWidth',2)
hold off
legend('Yc','Yd','Yf','Ytp','Yzz','Yode')
figure
plot(Tc,error)


figure(3)
plot(Td_seq,Yc);           %by LTI block
title('\bf\fontsize{14}\fontname{Cambria} Time-delay by LTI block')
xlabel('\bf\fontsize{10} t(時間)')
ylabel('\bf\fontsize{10} 系統輸出')
axis([-inf inf -inf inf]);
figure(4)
plot(Td_seq,Ytp);            %by TR block

title('\bf\fontsize{14}\fontname{Cambria} Time-delay by TR block')
xlabel('\bf\fontsize{10} t(時間)')
ylabel('\bf\fontsize{10}系統輸出')
axis([-inf inf -inf inf]);
figure(5)
plot(Td_seq,Yzz);           %Time-delay by simulink
title('\bf\fontsize{14}\fontname{Cambria} discrete-time system by simulink')
xlabel('\bf\fontsize{10} t(時間)')
ylabel('\bf\fontsize{10}離散系統輸出')
axis([-inf inf -inf inf]);
figure(6)
plot(Td_seq,Y_ode);                    %Time-delay by ode45
title('\bf\fontsize{14}\fontname{Cambria} discrete-system by loop')
xlabel('\bf\fontsize{10} t(時間)')
ylabel('\bf\fontsize{10}系統輸出')
axis([-inf inf -inf inf]);

%---(1)error for LTI
error_LTI=abs(Yc-Yf)./abs(Yc)*100;
%---(2)error for TR
error_TR=abs(Ytp'-Yf)./abs(Ytp')*100;
%---(3)error for ode45
error_ode=abs(Y_ode'-Yf)./abs(Y_ode')*100;
%---(4)error for loop 疊代和 simulink 疊代
error_IT=abs(Yzz'-Yf)./abs(Yzz')*100;



figure(7)
plot(Td_seq,error_LTI);          %----時延系統以LTI  Block模擬的誤差圖
title('\bf\fontsize{14}\fontname{Cambria}  Error about using LTI block')
xlabel('\bf\fontsize{10} t(時間)')
ylabel('\bf\fontsize{10} 誤差')
axis([-inf inf -inf inf]);
figure(8)
plot(Td_seq,error_TR);           %----時延系統以Transport  Delay模擬的誤差圖
title('\bf\fontsize{14}\fontname{Cambria}  Error about using TR block')
xlabel('\bf\fontsize{10} t(時間)')
ylabel('\bf\fontsize{10}誤差')
axis([-inf inf -inf inf]);
figure(9)
plot(Td_seq,error_ode);          %---時延系統以simulink模擬的誤差圖
title('\bf\fontsize{14}\fontname{Cambria}  Error about using ode45')
xlabel('\bf\fontsize{10} t(時間)')
ylabel('\bf\fontsize{10} 誤差')
axis([-inf inf -inf inf]);
figure(10)
plot(Td_seq,error_IT);          %----loop 疊代和 simulink 疊代
title('\bf\fontsize{14}\fontname{Cambria}  Error about using Z^{-1}')
xlabel('\bf\fontsize{10} t(時間)')
ylabel('\bf\fontsize{10} 誤差')
axis([-inf inf -inf inf]);




