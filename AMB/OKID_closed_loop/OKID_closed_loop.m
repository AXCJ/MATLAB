clc; clear; close all;
% load SystemData.mat % load original system
sampleData = 100; % number of sample data 
deletePoint = 1; % from the first data to the data you choose won't be shown in figure
showPoint = sampleData-deletePoint; % total plot point
plotPoint = deletePoint + showPoint;

Tend = 0.5;
Ts = Tend / sampleData;
simT = 0 : Ts : Tend-Ts;
MarkovOrder = 2;

P = [2 3 4 5 6 7 8 9];
A = zp2ss(2, P ,8);
Pd = exp(-P * Ts);
load('OKIDData.mat' , 'H' , 'C', 'D')
B = H;
[G,H] = c2d(A,B,Ts);
Kd = place(G,H,Pd);
eig(G-H*Kd)
input_weight = [1 2 3 4 5 6 7 8]; % input wighting make input to be same level

kk = 0;
% input generator
for t = 1:1:sampleData
    kk = kk + 1;
    r_randn = [     input_weight(1)*randn(1,1);
                             input_weight(2)*randn(1,1);
                             input_weight(3)*randn(1,1);
                             input_weight(4)*randn(1,1);
                             input_weight(5)*randn(1,1);
                             input_weight(6)*randn(1,1);
                             input_weight(7)*randn(1,1);
                             input_weight(8)*randn(1,1)];
    r(:,kk) =1+r_randn;
end

r = [simT ;r]';
sim('sim_OKID_closed_loop',[simT])

sim('sim_OKID_open_loop',[simT])
Ud = rw';
Yd = y';
x_initial = zeros(8,1);
% OKID the system
IDset_ =  struct('MarkovOrder', [MarkovOrder], 'Alpha', 2, 'Beta', 2, 'n', 0.0001, 'MinRA', 'era');
[G_ok,H_ok,C_ok,D_ok,F_ok,Sigma, er, M, Ob_CanF] = ...
    OKID_fun_WT02_6(Ud, Yd, IDset_, []);

Qd = 1e2*eye(size(C_ok,1));
Rd = 1e-15*eye(size(B,2));

[Kd,P] = dlqr(G_ok,H_ok,C_ok'*Qd*C_ok,Rd,0);
sim('sim_Observer_based_OKID_closed_loop',[0:0.1:10])
% OKID_system output generator 
x_o(:,1) = pinv(C_ok)*(C*x_initial);
for kk = 1 : sampleData
    y_o(:,kk) = C_ok*x_o(:,kk);
    x_o(:,kk+1) = G_ok*x_o(:,kk)  + H_ok*Ud(:,kk) + F_ok*(y_o(:,kk) - Yd(:,kk)) ;
end
% error between original & OKIDed system.
error_Y = Yd - y_o;

ii = 0;
for t = 0 : 1 : sampleData-1
    ii = ii + 1;
    td(:,ii) = ii;
end

eigG = eig(G)
eigGok = eig(G_ok)

 
scrsz = get(groot,'ScreenSize');
% [Left bottom width height]
figure('Position',[9 scrsz(4)/1.7 scrsz(3)/3 scrsz(4)/3]) % Left_Top of screen
plot(td(deletePoint:plotPoint),error_Y(deletePoint:plotPoint)) 
figure('Position',[scrsz(3)/3 scrsz(4)/1.7 scrsz(3)/3 scrsz(4)/3]) % Middle_Top of screen
plot(td(deletePoint:plotPoint),Yd(deletePoint:plotPoint),td(deletePoint:plotPoint),y_o(deletePoint:plotPoint),'.-')
legend('y_{d}','y_{o}')

figure('Position',[scrsz(3)/1.5 scrsz(4)/1.7 scrsz(3)/3 scrsz(4)/3]) % Right_Top of screen
plot(td(deletePoint:plotPoint),error_Y(deletePoint:plotPoint)) 
figure('Position',[9 scrsz(4)/1.7/4 scrsz(3)/3 scrsz(4)/3]) % Left_Bottom of screen
plot(td(deletePoint:plotPoint),Yd(deletePoint:plotPoint),td(deletePoint:plotPoint),y_o(deletePoint:plotPoint),'.-') 
legend('y_{d}','y_{o}')
figure('Position',[scrsz(3)/1.5 scrsz(4)/1.7/4 scrsz(3)/3 scrsz(4)/3]) % Right_Bottom of screen
plot([0:0.1:10],Y) 
