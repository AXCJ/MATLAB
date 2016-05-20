clc;clear all;close all;
% format short
% Tf = 0.002; %2ms
% Tend = 30; %30s
Tf = 1e-3; %1ms
Tend = 10;  %5s
Sim_t = 0:Tf:Tend;

parameter = [0.2,0.8; 0.1 0.5; 0.1 1.2; 0.3 0.5; 0.3 1.2];

%% == Reference  r(t)
r = ones(1,size(Sim_t,2));

for i = 1:size(parameter,1)
    x = 0;
    x0_o = [0.5 0.5]';
    xi = parameter(i,1);
    omegaC = parameter(i,2);

    sys = tf([1],[1 2*xi*omegaC omegaC^2]);
    [A,B,C,D] = tf2ss([1],[1 2*xi*omegaC omegaC^2]);
    
    nr = size(A,1);     % No. of states
    mr = size(B,2);     % No. of inputs
    pr = size(C,1);     % No. of outputs
    
    Qc = 1e4*eye(pr);
    Rc = eye(mr);
    
    Rc_bar = Rc + D'*Qc*D;
    [F,Pc] = lqr(A, B, C'*Qc*C, Rc_bar);
    F = -F;
%     F = -place(A,B,[-30+5i -30-5i]);
    M = -inv(C*inv(A+B*F)*B);
    
    Ao=A';
    Bo=C';
    Co = B';
    Qo = 1e3*eye(size(Ao,1));
    [H,Po] = lqr(Ao, Bo, Qo, Rc_bar);
    H = -H';
    
    opts_sim = simset('Solver', 'ode45', 'FixedStep', Tf, 'InitialStep', Tf, 'MaxStep', Tf, 'RelTol', 1e-7);
    [T, X, U, Y] = sim('trackerSys', Sim_t, opts_sim, [ Sim_t', [r'] ]);
    T = T';
    Y = Y';
    
    figure(1);
    hold on
    plot(T, Y); %1m = 1000mm
    grid on;
    
    S(i) = stepinfo(Y,T);
    S(i).SettlingTime
    ssError(i) = Y(end)-1
end


% title('Left rotor orbit');
% xlabel('x_{1} - axis (mm)');
% ylabel('y_{1} - axis (mm)');