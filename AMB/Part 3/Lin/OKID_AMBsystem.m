clc,clear ;close all

%% system necessary parameter
Ts = 2.5*10^(-4);
Tend = 0.05;
Sam_t = 0:Ts:Tend-Ts;
T_length = length(Sam_t);

m = 2.56478;
L = 0.505;
rd = 0.0166; %Diameter of rotor
J = 0.04004;
Jz = 0.0006565;
kri = 80;
krp = 220000; 
kai = 40;
kap = 36000;
a = -0.16;
b = 0.19;
l = -a+b; 
c = a;
d = b;
% c = 0.263 d = l-c;

omega = 2*pi/Ts; % Speed

alpha1 = a*Jz*omega/(J*l);
alpha2 = b*Jz*omega/(J*l);
beta1 = 1/m + a*a/J;
beta2 = 1/m + a*b/J;
beta3 = 1/m + b*b/J;
beta4 = 1/m;
gamma1 = 1/m - a*c/J;
gamma2 = 1/m + b*c/J;
gamma3 = 1/m;
%--------------------------------------%
%%

x_initial = zeros(10,1);
ii=0;
for t=0:Ts:Tend-Ts
    ii=ii+1;
    u_star(:,ii) = [60*sin(omega*Sam_t(1,ii)/10) ;
                        80*cos(omega*Sam_t(1,ii)/10) ;
                        70*cos(omega*Sam_t(1,ii)/10) ;
                        60*sin(omega*Sam_t(1,ii)/10) ;
                        0.001*cos(omega*Sam_t(1,ii)/10)];
end
randn_1000=randn(5,T_length);
dev = [0 0 0 0 0];

% kk = 0;
% for t = 0:2.5*10^(-4):0.04-2.5*10^(-4)
%     kk = kk + 1;
%     u_randn = [    dev(1)*randn_1000(1,kk);
%                                 dev(2)*randn_1000(2,kk); 
%                                 dev(3)*randn_1000(3,kk);
%                                 dev(4)*randn_1000(4,kk);
%                                 dev(5)*randn_1000(5,kk);
%                                 dev(6)*randn_1000(6,kk);
%                                 dev(7)*randn_1000(7,kk);
%                                 dev(8)*randn_1000(8,kk)
%                            ];
                            
%     esys(:,kk) = -0.0073 + 13.3067*randn(1,1);
% 
%     u(:,kk) = u_star ;
% end

Ag = [zeros(2) [-alpha1 alpha1; alpha2 -alpha2] [0;0]; [-alpha1 alpha1; -alpha2 alpha2] zeros(2,3); zeros(1,5)];
Ae = [[krp*beta1 2*krp*beta2; 2*krp*beta2 krp*beta3] zeros(2,3); zeros(3,2) [krp*beta1 2*krp*beta2 0; 2*krp*beta2 krp*beta3 0; 0 0  kai*beta4]];
Bi = [[kri*beta1 2*kri*beta2; 2*kri*beta2 kri*beta3] zeros(2,3);zeros(3,2) [kri*beta1 2*kri*beta2 0; 2*kri*beta2 kri*beta3 0; 0 0 kai*beta4]];
Ci = [[c 1;d 1] zeros(2,3); zeros(3,2) [c 1 0; d 1 0; 0 0 1]];
Eg_s = [0 0 -1 -1 0]';
T = [[a b;1 1] zeros(2,3); zeros(3,2) [a b 0; 1 1 0; 0 0 1]];

% n = size(Ae,1);
% m = size(Bi,2);
% p = size(Ci,1);

A = [zeros(5) eye(5); Ae Ag];
B = [zeros(5); Bi];
C = [Ci*inv(T') zeros(5)];
Eg = [zeros(5,1); Eg_s];
g = Eg*9.8;

nr = size(A,1);     % No. of states 
mr = size(B,2);     % No. of inputs
pr = size(C,1);     % No. of outputs
D  = zeros(pr,mr);    % Without a direct feed-through term
x(:,1) = zeros(nr,1); % Initial condition

[G_d ,H_d]=c2d(A,B,2.5*(10^(-4)));


kk = 0;
xd(:,1) = x_initial;
Ts=0.01;
 for t = 0:2.5*10^(-4):Tend-2.5*10^(-4)
    kk = kk + 1;
    xd(:,kk+1) = G_d*xd(:,kk) + H_d*u_star(:,kk)+(G_d-eye(10))*inv(A)*g;
    y(:,kk) = C*xd(:,kk);
 end

upsilon = [u_star; y];

Alpha = 20;
Beta = 20;
n = 10;
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
for kk = 1:T_length
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
figure(2)
plot(Sam_t,y(1,:),'o',Sam_t,y_o(1,:),'-')
legend('y_{d}','y_{o}')
figure(3)
plot(Sam_t,y(2,:),'o',Sam_t,y_o(2,:),'.-')
legend('y_{d}','y_{o}')
figure(4)
plot(Sam_t,y(3,:),'o',Sam_t,y_o(3,:),'.-')
legend('y_{d}','y_{o}')
figure(5)
plot(Sam_t,y(4,:),'o',Sam_t,y_o(4,:),'.-')
legend('y_{d}','y_{o}')
figure(6)
plot(Sam_t,y(5,:),'o',Sam_t,y_o(5,:),'.-')
legend('y_{d}','y_{o}')