clc;clear all;close all;
% format short

load('OKID_XY.mat')
load('sys_XY.mat')
Tend = 0.03;
Sim_t = 0:Ts:Tend;
T_length = length(Sim_t);

%% == Reference  r(t)
ii = 0;
for t = 0 : Ts : Tend
    ii = ii+1;
    r(:,ii) = 0.1; %Unit(m)
end
%------------------------------%
m = 2.56478;
L = 0.505;
rd = 0.0166; %Diameter of rotor
J = 0.04004;
Jz = 0.0006565;
kri = 80;
krp = 220000;
kai = 40;
kap = 36000;
a = -0.16; b = 0.19;
l = -a+b;
c = a; d = b;
% c = 0.263 d = l-c;
xb = 0.4; xa = 0.4;
z0 = 0.5;
%--------%
ib = 0.9;
i0 = 1.1;
i0_z = 0.7;
Lc = 0.265;
Rc = 0.27;
%--------%
omega = 2400; % Speed

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

%% new system

%% new system


% Eg = [zeros(5,1); Eg_s];
% g = Eg*9.806;
nr = size(A,1);     % No. of states
mr = size(B,2);     % No. of inputs
pr = size(C,1);     % No. of outputs
D  = zeros(size(C,1),size(B,2));    % Without a direct feed-through term

G = A_ok; H = B_ok;  
x(:,1) = zeros(size(A,1),1); % Initial condition
% x = [5e-3 0 5e-3 0 0 0 0 0 0 0]'; %test
% pzmap(A,B,C,D);
Mc = ctrb(A, B);  % Controllability matrix
Mo = obsv(A, C);  % Observability matrix

Rc = rank(Mc);      % Check whether the plant is controllable
Ro = rank(Mo);      % Check whether the plant is observable

Qd = eye(pr);
Rd = 1e-15*eye(mr);
% Rd(5,5) = 1e-2;

gd = (G-eye(nr))*inv(A)*g;
Rd_bar = Rd + D'*Qd*D;
Nd = C_ok'*Qd*D;
[Kd,P] = dlqr(G,H,C_ok'*Qd*C_ok,Rd_bar,Nd);
Rd_wave = Rd_bar + H'*P*H;
P_bar = H'*P*G+Nd';
Ed =  inv(Rd_wave)*(D'+H'*inv(eye(nr)-(G-H*Kd)')*(C_ok-D*Kd)')*Qd;
Cd = inv(Rd_wave)*H'*(inv((G-H*Kd)'-eye(nr))*(G-H*Kd)'-eye(nr))*P*gd; %Ignore output noise s(k)
Cu_star = inv(Rd_wave)*(H'*inv((G-H*Kd)'-eye(nr))*Kd'+eye(mr))*Rd;
ud_star = [0 0 0 0]';
% i = [ixA ixB iyA iyB]'
% x = [beta x -alpha y]
% y = [XseA XseB YseA YseB]
Cd=0;

% Model-Following Approach

% ~~~ Perform Simulation
sys = ss(G, H, C_ok, D, Ts);

figure(1);
[yss,tss] = step(sys);
yss = sum(yss, 5);
str = ['x_{seA}(m)';'x_{seB}(m)'; 'y_{seA}(m)'; 'y_{seB}(m)'];
for i = 1 : 4
    subplot(220+i);
    plot(tss, yss(:,i));
    title('Open-loop step response');
    xlabel('Time (sec)');
    ylabel(str(i,:));
end
figure(2);
[po, zo] = pzmap(sys)
pzmap(sys);
opts_sim = simset('Solver', 'ode45', 'RelTol', 1e-7);
[T, X, U, Y] = sim('Model_Following_Destrete_controlled', Sim_t, opts_sim, [ Sim_t', [r' r' r' r']]);

sys2 = ss(G-H*Kd, H*Ed, C_ok, D, Ts);
figure(3);
[yss2,tss2] = step(sys2);
yss22 = sum(yss2, 3);
for i = 1 : 4
    subplot(220+i);
    plot(tss2, yss22(:,i));
    title('Closed-loop step response');
    xlabel('Time (sec)');
    ylabel(str(i,:));
end
figure(4)
[pc, zc] = pzmap(sys2)
pzmap(sys2);

figure(5);
plot(1e3*Y(:,1), 1e3*Y(:,3)); %1m = 1000mm
grid on;
title('Left rotor orbit');
xlabel('x_{1} - axis (mm)');
ylabel('y_{1} - axis (mm)');

figure(6);
plot(1e3*Y(:,2), 1e3*Y(:,4));
grid on;
title('Right rotor orbit');
xlabel('x_{2} - axis (mm)');
ylabel('y_{2} - axis (mm)');

figure(7);
plot(T, 1e3*Y);
% figure(7);
% subplot(211);
% plot(T, ib-U(:,1));
% title('Total currents in x_{1} - axis');
% xlabel('Time (sec)');
% ylabel('i_{b} - i_{x1} (A)');
% 
% subplot(212);
% plot(T, ib+U(:,1));
% xlabel('Time (sec)');
% ylabel('i_{b} + i_{x1} (A)');
% 
% figure(8);
% subplot(211);
% plot(T, ib-U(:,3));
% title('Total currents in y_{1} - axis');
% xlabel('Time (sec)');
% ylabel('i_{b} - i_{y1} (A)');
% 
% subplot(212);
% plot(T, ib+U(:,3));
% xlabel('Time (sec)');
% ylabel('i_{b} + i_{y1} (A)');
% figure(9);
% str = ['i_{x1} (A)'; 'i_{x2} (A)'; 'i_{y1} (A)'; 'i_{y2} (A)'];
% for i = 1 :  4
%     subplot(220+i);
%     plot(T, U(:,i));
%     title('Current');
%     xlabel('Time (sec)');
%     ylabel(str(i,:));
% end
% 
% figure(10);
% plot(T,U(:,5)+1.1);
% title('Total current in z - axis');
% xlabel('Time (sec)');
% ylabel('i_{0} + i_{z} (A)');
% 
% figure(11);
% plot(T,U(:,5));
% title('Current');
% xlabel('Time (sec)');
% ylabel('i_{z} (A)');
% 
% figure(12);
% plot(T,Y(:,5));
% title('Thrust displacement in z - axis');
% xlabel('Time (sec)');
% ylabel('z - axis (mm)');

[G1,H1]=c2d(A,B,Ts);
% eig_G=eig(G1)
% eig_Aok=eig(A_ok)
error_eig=eig(G1)-eig(A_ok)



