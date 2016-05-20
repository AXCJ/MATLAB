clc;clear all;close all;pause(0.01)


A=[0 1;-3 -2];
B=[0; 1];
C=[1 1];

Ts=0.18;
Tf=0.06;
T_Sim=5.4;
Td_seq=0:Tf:T_Sim;
% r=wgn(1,length(Td_seq),1);

x0=[0;0];
load noise_r.mat
r(1,:)=[];
ii=0;kk=-3;

%ZOH
for t=0:Tf:T_Sim
    ii=ii+1;
    if mod(ii-1,(Ts/Tf)) == 0 
        
        kk=kk+3;
    end

    
    if t<T_Sim & kk<92
        rc(:,ii) = r(:,kk+1);
    end
end
options = simset('FixedStep',Tf,'MaxStep',Tf,'initialstep',Tf,'initialstate',x0);
[Tc,Xc,Yc]=sim('CTD_system_Based_TF',Td_seq,options,[Td_seq', rc']);
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

%
[n,m]=size(B);
[p,n]=size(C);
D=zeros(p,m);

%
[G,H]=c2d(A,B,Ts);
Gamma_0=(real(G^(1/3))-eye(n))*inv(A)*B;
Gamma_1=(G-real(G^(1/3)))*inv(A)*B;
PI_ok=real(G^(-2/3))*H;

ii=0;kk=0;
for t=0:Tf:T_Sim
    ii=ii+1;
    
    if mod(ii-1,(Ts/Tf))==0
        kk=kk+1;
        if kk==1
            xf(:,kk)=pinv(C)*(C*x0);
        elseif kk==2
            xf(:,kk)=G*xf(:,kk-1)+Gamma_0*Us(:,kk-1)
        else
            xf(:,kk)=G*xf(:,kk-1)+Gamma_0*Us(:,kk-1)+Gamma_1*Us(:,kk-2);
        end
        Us(:,kk) = r(:,3*(kk-1)+1);
    end
    
    if mod(ii-1,(Ts/Tf))==0
        yf(:,kk)=C*xf(:,kk);
        th(:,kk)=(kk-1)*Ts;
    end
    
    disp(sprintf('Simulation Time:%2.3f sec',(ii-1)*Tf));
end

for kk=1:size(Yd,2)
    e(:,kk)=Yd(:,kk)-yf(:,kk);
    if abs(Yd(:,kk))==0
        Yd(:,kk)=10^-6;
        error(:,kk)=(abs(e(:,kk))./abs(Yd(:,kk)))*100;
    else
        error(:,kk)=(abs(e(:,kk))./abs(Yd(:,kk)))*100;
    end
end

plot(Tc,Yc)
hold on
plot(Td,Yd,'o')
plot(th,yf)
hold off
legend('Yc','Yd','yf')
figure
plot(Td,error)




