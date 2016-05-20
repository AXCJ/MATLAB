clc,clear ;close all

%% system necessary parameter
Ts = 2.5e-4; %0.25ms
Tf = 10e-6; %10us
Tend = 5e-2;  %0.5s
Sam_t = 0:Ts:Tend;
T_length = length(Sam_t);
load AMB_discrete_parameter;
omega = 2*pi/Ts; 
%-------------%
%%
x_initial = zeros(2,1);
ii=0;
for t=0:Ts:Tend-Ts
    ii=ii+1;
    u_star(:,ii) = [  0.001*cos(omega*Sam_t(1,ii)/10)];
end
tt=0: 2.5e-4:Tend-Ts;
pp=[tt' u_star'];
randn_1000=randn(1,T_length);
dev = [1.4];

kk = 0;
for t =  0:Ts:Tend-Ts
    kk = kk + 1;
    u_randn(:,kk) = [  dev(1)*randn_1000(1,kk) ];   
end
sim('zoh')
simout_zoh = simout_zoh';
u = simout_zoh  +  u_randn;
%%
%AMB system Matrix without XY
Ag = [zeros(1,1)];
Ae = [kai*beta4];
Bi = [kai*beta4];
Ci = [1];
Eg_s = [0]';
T = [1];

% n = size(Ae,1);
% m = size(Bi,2);
% p = size(Ci,1);
A = [zeros(1) eye(1); Ae Ag];
B = [zeros(1); Bi];
C = [Ci*inv(T') zeros(1)];
Eg = [zeros(1,1); Eg_s];
g = Eg*9.8;

nr = size(A,1);     % No. of states 
mr = size(B,2);     % No. of inputs
pr = size(C,1);     % No. of outputs
D  = zeros(pr,mr);    % Without a direct feed-through term
x = zeros(nr,1); % Initial condition

[G_d ,H_d]=c2d(A,B,2.5*(10^(-4)));


kk = 0;
xd(:,1) = x_initial;
Ts=0.01;
 for t = 0:2.5*10^(-4):Tend-2.5*10^(-4)
    kk = kk + 1;
    xd(:,kk+1) = G_d*xd(:,kk) + H_d*u_star(:,kk);
    y(:,kk) = C*xd(:,kk);
 end

upsilon = [simout_zoh; y];

Alpha = 3;
Beta = 3;
n = 2;
[m, l] = size(y);
r = size(u_star,1);
v = [u_star; y];
p = 2;
%%
index_p = p + 1;    
V_bar = [u_star(:, index_p:l)];
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
for kk = 1:200
    y_o(:,kk) = C_ok*x_o(:,kk);
    x_o(:,kk+1) = A_ok*x_o(:,kk)  + B_ok*u_star(:,kk) - G_ok*( y(:,kk) -  y_o(:,kk));
end

ii = 0;
for t = 0 : 2.5*10^(-4) : 0.04-2.5*10^(-4)
    ii = ii + 1;
    td(:,ii) = ii;
end

% save OKID_esys_without.mat G_ok H_ok C_ok D_ok F_ok
% figure(1)
% plot(Sam_t,)
Sam_t(:,end)=[];

figure(1)
plot(Sam_t,y(1,:),'-',Sam_t,y_o(1,:),'*')
xlabel('Time (sec)')
ylabel('z(t)')
legend('z_{d}','z_{d\_ok}')

%% ~~ error ~~Start
y_er(1,:)=(y_o(1,:)-y(1,:))./y(1,:)*100;
figure(2)
plot(Sam_t,y_er(1,:),'-.b')
xlabel('Time (sec)')
ylabel('Relative error (%)')
legend('z_{d\_error}')
%% ~~ error ~~End

%% ~~pole zero ~~Start
sys1=ss(A_ok,B_ok,C_ok,[]);
[p1 z1]=pzmap(sys1)
%% ~~ploe zero ~~End

A_ok_z = A_ok;
B_ok_z = B_ok;
C_ok_z = C_ok;
G_ok_z = G_ok;

save OKID_Z.mat 'A_ok_z' 'B_ok_z' 'C_ok_z' 'G_ok_z'