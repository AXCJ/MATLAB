 clc,clear all;close all

%%AMB system necessary parameter
load AMB_discrete_parameter;
Ts = 1e-6;
sample = 100;
% Tf = 10e-6; %10us
Tend = Ts*sample;
Sam_t = 0:Ts:Tend;
T_length = length(Sam_t)
omega = 2*pi/Ts; 
%--------------------------------------%
%%
%Get input and output
x_initial = zeros(8,1);
ii=0;
tt=0:Ts:Tend-Ts;
for t=0:Ts:Tend-Ts
    ii=ii+1;
    u_star(:,ii) = [50*randn(1,1);
                          30*randn(1,1);
                          40*randn(1,1);
                         35*randn(1,1)];
end
pp=[tt' u_star'];
randn_1000=randn(4,T_length);
dev = [1.2 1.1 1.3 1.5];

kk = 0;
for t = 0:Ts:Tend-Ts
    kk = kk + 1;
    u_randn(:,kk) = [   dev(1)*randn_1000(1,kk);
                                           dev(2)*randn_1000(2,kk); 
                                           dev(3)*randn_1000(3,kk);
                                           dev(4)*randn_1000(4,kk);                                          
                                       ];
end
sim('zoh')
simout_zoh = simout_zoh';
u = simout_zoh  +  u_randn;
%%
%AMB system Matrix without Z
Ag = [zeros(2) [-alpha1 alpha1; alpha2 -alpha2]; [-alpha1 alpha1; -alpha2 alpha2] zeros(2,2)];
Ae = [[krp*beta1 2*krp*beta2; 2*krp*beta2 krp*beta3] zeros(2,2); zeros(2,2) [krp*beta1 2*krp*beta2 ; 2*krp*beta2 krp*beta3]];
Bi = [[kri*beta1 2*kri*beta2; 2*kri*beta2 kri*beta3] zeros(2,2);zeros(2,2) [kri*beta1 2*kri*beta2 ; 2*kri*beta2 kri*beta3]];
Ci = [[c 1;d 1] zeros(2,2); zeros(2,2) [c 1 ; d 1]];
Eg_s = [0 0 -1 -1]';
T = [[a b;1 1] zeros(2,2); zeros(2,2) [a b;1 1]];

% n = size(Ae,1);
% m = size(Bi,2);
% p = size(Ci,1);

A = [zeros(4) eye(4); Ae Ag];
B = [zeros(4); Bi];
C = [Ci*inv(T') zeros(4)];
Eg = [zeros(4,1); Eg_s];
g = Eg*9.8;

nr = size(A,1);     % No. of states 
mr = size(B,2);     % No. of inputs
pr = size(C,1);     % No. of outputs
D  = zeros(pr,mr);    % Without a direct feed-through term
x = zeros(nr,1); % Initial condition

save sys_XY.mat 'A' 'B' 'C' 'g' 'Ts' 'Tend'
%%
[G_d ,H_d]=c2d(A,B,Ts);
%Do OKID
kk = 0;
xd(:,1) = x_initial;
 for t = 0:Ts:Tend-Ts
    kk = kk + 1;
    xd(:,kk+1) = G_d*xd(:,kk) + H_d*simout_zoh(:,kk)+(G_d-eye(8))*inv(A)*g;
    y(:,kk) = C*xd(:,kk);
 end

upsilon = [simout_zoh; y];

Alpha = 3;
Beta = 3;
n = 8;
[m, l] = size(y);
r = size(simout_zoh,1);
v = [simout_zoh; y];
p = 2;
%%
index_p = p + 1;    
V_bar = [simout_zoh(:, index_p:l)];
V_bar=[];
for ii = 1:p
    V_bar = [V_bar ; upsilon(:, (index_p-ii):(l-ii) )];
end
y_bar = [y(:, index_p:l)];
Y_bar = y_bar*pinv(V_bar);% System markov parameters

% Get Observer system markov parameters
for ii = 1:p
    Yk_bar(:, :, ii) = Y_bar(:, ((ii-1)*(r+m)+1):(ii*(r+m)) );
end

for jj = 1:length(Yk_bar(1, 1, :))
    Yk1_bar(:, :, jj) = Yk_bar(:, 1:r, jj);
    Yk2_bar(:, :, jj) = -Yk_bar(:, ((r+1):(r+m)), jj);
end
% End of get Observer system markov parameter

Pk(:, :, 1) = [ (Yk1_bar(:, :, 1)  - Yk2_bar(:, :, 1)*D),  Yk2_bar(:, :, 1)];
Yk(:, :, 1) = Pk(:, 1:r, 1);
Yok(:, :, 1) = Pk(:, (r+1):(r+m), 1);

for i = 2:Alpha + Beta
    if (i <= p)
        sumTempPk = zeros(size(Pk(:, :, 1) ));
        for j = 1:(i-1)
            sumTempPk = sumTempPk + Yk2_bar(:, :, j)*Pk(:, :, i-j);
        end
        Pk(:, :, i) = [ (Yk1_bar(:, :, i)-Yk2_bar(:, :, i)*D),  Yk2_bar(:, :, i)] - sumTempPk;
    else
        sumTempPk = zeros(size(Pk(:, :, 1) ));
        for j = 1:p
            sumTempPk = sumTempPk + Yk2_bar(:, :, j)*Pk(:, :, i-j);
        end
        Pk(:, :, i) = -sumTempPk;
    end
    
    Yk(:, :, i) = Pk(:, 1:r, i);
    Yok(:, :, i) = Pk(:, (r+1):(r+m), i);
end

H =  [];
for i = 1:Alpha
    Htemp = [];
    for j = 1:(Beta+1)
        Htemp = [Htemp Pk(:, :, i+j-1)];
    end
    H = [H; Htemp];
end
H0 = H;
H0(:, (end-(r+m)+1):end) = [];
H1 = H;
H1(:, 1:(r+m) ) = [];

[A_ok, Btemp, C_ok, er, U, Sigma, V] = era_WT( H1, H0, r+m, m, n);
B_ok = Btemp(:, 1:r); % B : ID systeim => B
G_ok= Btemp(:, (r+1):(r+m) ); % G : ID systeim => G

x_o(:,1) = pinv(C_ok)*(C*x_initial);
for kk = 1:T_length-1
    y_o(:,kk) = C_ok*x_o(:,kk);
    x_o(:,kk+1) = A_ok*x_o(:,kk)  + B_ok*simout_zoh(:,kk) - G_ok*( y(:,kk) -  y_o(:,kk));
end

ii = 0;
for t = 0 : Ts : 0.04-Ts
    ii = ii + 1;
    td(:,ii) = ii;
end
Sam_t(:,end)=[];

i = 10;
figure(1)
plot(Sam_t(i:end),y(1,i:end),'*',Sam_t(i:end),y_o(1,i:end),'-')
xlabel('Time (sec)')
ylabel('x1(t)')
legend('x_{d1}','x_{d1\_ok}')
figure(2)
plot(Sam_t(i:end),y(2,i:end),'*',Sam_t(i:end),y_o(2,i:end),'.-')
xlabel('Time (sec)')
ylabel('x2(t)')
legend('x_{d2}','x_{d2\_ok}')
figure(3)
plot(Sam_t(i:end),y(3,i:end),'*',Sam_t(i:end),y_o(3,i:end),'.-')
xlabel('Time (sec)')
ylabel('y1(t)')
legend('y_{d1}','y_{d1\_ok}')
figure(4)
plot(Sam_t(i:end),y(4,i:end),'*',Sam_t(i:end),y_o(4,i:end),'.-')
xlabel('Time (sec)')
ylabel('y2(t)')
legend('y_{d2}','y_{d2\_ok}')

y_er=[];
%% ~~ error~~Start
% y_er(1,:)=(y_o(1,:)-y(1,:))./y(1,:)*100;
% y_er(2,:)=(y_o(2,:)-y(2,:))./y(2,:)*100;
% y_er(3,:)=(y_o(3,:)-y(3,:))./y(3,:)*100;
% y_er(4,:)=(y_o(4,:)-y(4,:))./y(4,:)*100;
y_er(1,:)=(y_o(1,:)-y(1,:));
y_er(2,:)=(y_o(2,:)-y(2,:));
y_er(3,:)=(y_o(3,:)-y(3,:));
y_er(4,:)=(y_o(4,:)-y(4,:));
figure(5)
plot(Sam_t(i:end),y_er(1,i:end),'-.b')
xlabel('Time (sec)')
ylabel('Relative error (%)')
legend('x_{d1\_error}')
figure(6)
plot(Sam_t(i:end),y_er(2,i:end),'-.b')
xlabel('Time (sec)')
ylabel('Relative error (%)')
legend('x_{d2\_error}')
figure(7)
plot(Sam_t(i:end),y_er(3,i:end),'-.b')
xlabel('Time (sec)')
ylabel('Relative error (%)')
legend('y_{d1\_error}')
figure(8)
plot(Sam_t(i:end),y_er(4,i:end),'-.b')
xlabel('Time (sec)')
ylabel('Relative error (%)')
legend('y_{d2\_error}')

%% ~~error~~End

%% ~~pole zero ~~Start
sys1=ss(A_ok,B_ok,C_ok,[],Ts);
[p1 z1]=pzmap(sys1)
%% ~~ploe zero ~~End

save OKID_XY.mat 'A_ok' 'B_ok' 'C_ok' 'G_ok'


