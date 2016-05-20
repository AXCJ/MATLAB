clc;clear all;close all;
% format short
% Tf = 0.002; %2ms
% Tend = 30; %30s
Tf = 2e-6; %2us
Tend = 0.5;  %0.5s
Sim_t = 0:Tf:Tend;
T_length = length(Sim_t);

 %% == Reference  r(t)  
ii = 0;
for t = 0 : Tf : Tend
    ii = ii+1;
    r(:,ii) = 0; %Unit(m)
   %r(:,ii) = 8e-3; % 8mm (test)
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
a = -0.16, b = 0.19;
l = -a+b; 
c = a, d = b;
% c = 0.263 d = l-c;
xb = 0.4, xa = 0.4;
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

Ag = [zeros(2) [-alpha1 alpha1; alpha2 -alpha2] [0;0]; [-alpha1 alpha1; -alpha2 alpha2] zeros(2,3); zeros(1,5)];
Ae = [[krp*beta1 2*krp*beta2; 2*krp*beta2 krp*beta3] zeros(2,3); zeros(3,2) [krp*beta1 2*krp*beta2 0; 2*krp*beta2 krp*beta3 0; 0 0  kai*beta4]];
Bi = [[kri*beta1 2*kri*beta2; 2*kri*beta2 kri*beta3] zeros(2,3);zeros(3,2) [kri*beta1 2*kri*beta2 0; 2*kri*beta2 kri*beta3 0; 0 0 kai*beta4]]
Ci = [[c 1;d 1] zeros(2,3); zeros(3,2) [c 1 0; d 1 0; 0 0 1]];
Eg_s = [0 0 -1 -1 0]';
T = [[a b;1 1] zeros(2,3); zeros(3,2) [a b 0; 1 1 0; 0 0 1]];

n = size(Ae,1);
m = size(Bi,2);
p = size(Ci,1);

A = [zeros(n) eye(n); Ae Ag];
B = [zeros(m); Bi];
C = [Ci*inv(T') zeros(5)];
Eg = [zeros(5,1); Eg_s];
g = Eg*9.8;

nr = size(A,1);     % No. of states 
mr = size(B,2);     % No. of inputs
pr = size(C,1);     % No. of outputs
D  = zeros(pr,mr);    % Without a direct feed-through term
x(:,1) = zeros(nr,1); % Initial condition
% x = [5e-3 0 5e-3 0 0 0 0 0 0 0]'; %test
pzmap(A,B,C,D);
Mc = ctrb(A, B);  % Controllability matrix 
Mo = obsv(A, C);  % Observability matrix

Rc = rank(Mc);      % Check whether the plant is controllable 
Ro = rank(Mo);      % Check whether the plant is observable  

Qc = eye(pr);  
Rc = 1e-9*eye(mr);

Rc_bar = Rc + D'*Qc*D;

[Kc,Pc] = lqr(A, B, C'*Qc*C, Rc_bar);
Ec = -inv(Rc_bar)*((C - D*Kc )*inv( A - B*Kc )*B - D )'*Qc;
Cc = inv(Rc_bar)*B'*inv((A-B*Kc)')*Pc*g;
Cuc_start = inv(Rc_bar)*(eye(mr) + B'*inv((A-B*Kc)')*Kc')*Rc;
uc_start = [0 0 0 0 0]';
% i = [ixA ixB iyA iyB]' 
% x = [beta x -alpha y]
% y = [XseA XseB YseA YseB]


% Model-Following Approach 

% ~~~ Perform Simulation
[y_ori,t_ori] = step(A, B, C, D);
sys = ss(A, B, C, D);

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
[po, zo] = pzmap(sys);
pzmap(sys);
opts_sim = simset('Solver', 'ode45', 'FixedStep', Tf, 'InitialStep', Tf, 'MaxStep', Tf, 'RelTol', 1e-7);
[T, X, U, Y] = sim('Model_Following', Sim_t, opts_sim, [ Sim_t', [r' r' r' r' r'] ]);

sys2 = ss(A-B*Kc, B*Ec, C, D);
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
[pc, zc] = pzmap(sys2);
pzmap(sys2);

figure(5);
plot(1e3*Y(:,1), 1e3*Y(:,3)); %1m = 1000mm
grid on;
title('Left rotor orbit');
%axis([-0.5 0.5 -0.5 0.5]);
xlabel('x_{1}-axis (mm)');
ylabel('y_{1}-axis (mm)');

figure(6);
plot(1e3*Y(:,2), 1e3*Y(:,4));
grid on;
title('Right rotor orbit');
xlabel('x_{2}- axis (mm)');
ylabel('y_{2}- axis (mm)');
%axis([-0.5 0.5 -0.5 0.5]);

figure(7);
plot(T,U(:,1));
str = ['i_{x1} (A)'; 'i_{x2} (A)'; 'i_{y1} (A)'; 'i_{y2} (A)'];
for i = 1 : 4
    subplot(220+i);
    plot(T,U(:,i));
    title('Control current');
    xlabel('Time (sec)');
    ylabel(str(i,:));
end
