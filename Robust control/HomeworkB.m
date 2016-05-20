% Optimal control
% HomeworkB
% Question1: 
%      Design an observer-based controller K(s) so that the
%      closed-loop system is stable.
clc;clear;close all;

num = 1;
den = [-1 1 0];
sys = tf(num,den);
[P , Z] = pzmap(sys)
if(find(P > 0))
    fprintf('The Open-loop pole is %d , so the system is unstable\n' , P)
else
    fprintf('The Open-loop system is stable\n')
end

[A,B,C,D] = tf2ss(num,den);
[n,r] = size(B);
p = size(C,1);

desirePole = [-4+10i,-4-10i];
desirePole1 = [-5+10i,-5-10i];
Kc = acker(A,B,desirePole);
H = acker(A',C',desirePole1)';

sys_ob = ss(A-B*Kc-H*C,H,-Kc,[]);
initial(sys_ob,[1 0])
title('Controller free response','fontsize',20)
M = -1/(C/(A-B*Kc)*B);
% M = 1;
Acl = [A-B*Kc B*Kc ; zeros(size(A-H*C,1),size(A-B*Kc,2))  A-H*C];
Bcl = [B*M;zeros(size(A-H*C,1),1)];

Ccl = [C zeros(size(C,1),size(Acl,2)-size(C,2))];
sys_cl = ss(Acl,Bcl,Ccl,[]);
figure
step(sys_cl)
title('Closed-loop step response','fontsize',20)
