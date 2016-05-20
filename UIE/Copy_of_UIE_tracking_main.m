clc; clear all; close all;
scrsz = get(groot,'ScreenSize'); % screen  size
RB = [scrsz(3)/1.45 scrsz(4)/13 scrsz(3)/3.3 scrsz(4)/2.5]; % Right_Bottom of screen
RT = [scrsz(3)/1.45 scrsz(4)/1.95 scrsz(3)/3.3 scrsz(4)/2.5]; % Right_Top of screen
LT = [30 scrsz(4)/1.95 scrsz(3)/3.3 scrsz(4)/2.5]; % Left_Top of screen
POS = LT; % figure position

Tbegin = 0;
Tf =1e-3;
Tend = 12;
t = [Tbegin : Tf : Tend]';
H = 50; % shift the y-axis (A+HI)
Hc = 0; % useless

A = [-31.31 0 -2.833e4 2; 0 -10.25 8001 0; 1 -1 0 0; 0 0 0 -100];
B = [28.06 0 0 1; 0 28.06 0 1]';
Bd = B;
C = [1 0 0 0;1 0 30 0];

[n,r] = size(B);
p = size(C,1);
zero = tzero(A,B,C, zeros(p, r))
n_ctrb_obsv = [n rank(ctrb(A,B)) rank(obsv(A,C))]
%%
% observer gain L design
Qo = 1e5*eye(n);
Ro = 1e0*eye(p);
L = lqr(A'+H*eye(n), C', Qo, Ro)' % x-x_hat=0
eigALC = eig((A-L*C))

% estimation gain Kd design
Qod = 1e7*eye(n);
Rod = 1e0*eye(p);
Lt = lqr((A-L*C)'+H*eye(n), C', Qod, Rod)'; % regulator
Kd = pinv(B)*Lt
eigAtLtC = eig((A-L*C)-Lt*C)

% subsystem of disturbance
Af = [-100 0; 0 -100];
Bf = [8 0; 0 8];
Cf = [12.5 0; 0 12.5];

de1 = disturb1(t);
% de2 = disturb2(t);
de3 = disturb3(t);
de = de1;

% get reference input
% ii = 0;
% for t_ = Tbegin : Tf : Tend
%     ii = ii + 1;
%     ref(ii, :) = ref_fun(t_, 1); 
% %     refStar(ii, :) = ref_fun(t_ + Tf, 1); 
% end
ref = ones(length(t),1);
% refStar = sin(2*[Tbegin+Tf : Tf : Tend+Tf]');
ref = [ref ref];
% refStar = [refStar refStar];
r_in = [t ref];

% ~~~ Controller design
Qc = 1e3*eye(p);
Rc = 1e0*eye(r);
Kc = lqr(A + Hc*eye(n), B, C'*Qc*C, Rc); % Tracker --> C'*Q*C
Ec = -inv(Rc)*B'*inv(A - B*Kc)'*C'*Qc;
% Ec = -pinv(C*inv(A - B*Kc)*B)
% ~~~ Controller design ~~~ End

eigOfAc = eig(A-B*Kc)
zero = tzero(A-B*Kc,B*Ec,C, zeros(p, p))
x_ini = [0 0 0 0];
xo_ini = [0 0 0 0];
% simulation
[tt , xx , de1_hat, y, y_hat, x] = sim('model_tracking',t);

err_de = de - de1_hat;
err_yr = y - ref;
err_y = y - y_hat;

% de1
for i = 1:2
figure('Name', ['de' num2str(i) ' vs. de_hat' num2str(i)], 'Position', POS);
plot(t, de(: , i), 'b'); hold on;  plot(t, de1_hat(: , i), 'r'); hold off;
legend({['$d_{e' num2str(i) '}$'] , ['$\hat{d}_{e' num2str(i) '}$']} , 'Interpreter','latex', 'fontsize',15)
title({['$d_{e' num2str(i) '}$ $ vs.$ $ \hat{d}_{e' num2str(i) '}$' ]} , 'Interpreter','latex', 'fontsize',15)
h=gca;
set(gca,'FontSize',15,'XTick',h.XTick,'YTick',h.YTick,'YLim',h.YLim)
end

% Estimation error
for i = 1:2
figure('Name', ['Estimation error' num2str(i)], 'Position', POS);
plot(t, err_de(: , i), 'b');
legend({['$e_{' num2str(i) '}$' ]} , 'Interpreter','latex', 'fontsize',15)
title(['$e_{' num2str(i) '}$' ] , 'Interpreter','latex', 'fontsize',15)
title({['$d_{e' num2str(i) '}$ $ -$ $ \hat{d}_{e' num2str(i) '}$' ]} , 'Interpreter','latex', 'fontsize',15)
h=gca;
set(gca,'FontSize',15,'XTick',h.XTick,'YTick',h.YTick,'YLim',h.YLim)
% axis([0 100 -0.23 0.13])
end

% y
for i = 1:p
figure('Name', ['y' num2str(i)], 'Position', POS);
plot(t, ref(: , i), 'k'); hold on; plot(t, y(: , i), 'b');   plot(t, y_hat(: , i), 'r'); hold off;
legend({['$ref_{' num2str(i) '}$'], ['$y_{' num2str(i) '}$'] , ['$\hat{y}_{' num2str(i) '}$']} , 'Interpreter','latex', 'fontsize',15)
% title({['$d_{e' num2str(i) '}$ $ vs.$ $ \hat{d}_{e' num2str(i) '}$' ]} , 'Interpreter','latex', 'fontsize',15)
h=gca;
set(gca,'FontSize',15,'XTick',h.XTick,'YTick',h.YTick,'YLim',h.YLim)
end

% error of y & ref
for i = 1:p
figure('Name', ['y' num2str(i) '- r' num2str(i)], 'Position', POS);
plot(t, err_yr(: , i), 'b');
% legend({['$y_{' num2str(i) '}$'] , ['$\hat{d}_{e' num2str(i) '}$']} , 'Interpreter','latex', 'fontsize',15)
title({['$y_{' num2str(i) '}$ $ - $ $ r_{' num2str(i) '}$' ]} , 'Interpreter','latex', 'fontsize',15)
h=gca;
set(gca,'FontSize',15,'XTick',h.XTick,'YTick',h.YTick,'YLim',h.YLim)
end

% error of y & y_hat
for i = 1:p
figure('Name', ['y' num2str(i) '- y' num2str(i) '_hat'], 'Position', POS);
plot(t, err_y(: , i), 'b');
% legend({['$y_{' num2str(i) '}$'] , ['$\hat{d}_{e' num2str(i) '}$']} , 'Interpreter','latex', 'fontsize',15)
title({['$y_{' num2str(i) '}-$ $ \hat{y}_{' num2str(i) '}$' ]} , 'Interpreter','latex', 'fontsize',15)
h=gca;
set(gca,'FontSize',15,'XTick',h.XTick,'YTick',h.YTick,'YLim',h.YLim)
end

% % state
% for i = 1:size(x,2)
% figure('Name', ['x' num2str(i)], 'Position', POS);
% plot(t, x(: , i), 'b');
% legend({['$x_{' num2str(i) '}$' ]} , 'Interpreter','latex', 'fontsize',15)
% title(['$x_{' num2str(i) '}$' ] , 'Interpreter','latex', 'fontsize',15)
% title({['$x_{' num2str(i) '}$' ]} , 'Interpreter','latex', 'fontsize',15)
% h=gca;
% set(gca,'FontSize',15,'XTick',h.XTick,'YTick',h.YTick,'YLim',h.YLim)
% % axis([0 100 -0.23 0.13])
% end