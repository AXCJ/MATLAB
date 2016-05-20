clear all;clc; close all;

Tf = 5e-5;
Ts = 0.01;
sample = Ts/Tf;
Tend = 3;
% ~~~ system setting
% A = [ 0.9   -2.1     0.3     0.5     0.9;
%       0.7    0.2     0.3       0     0.7;
%      -1.3    0.5    -1.1    -2.3    -0.2;
%      -0.3    0.8     1.7    -1.2    -0.3;
%      -3.5   -4.3    -0.7       0    -8.3];
% B = [   1   -0.4;
%      -1.7   -1.7;
%      -0.2    1.2; 
%       0.7    0.1;
%       0.9    1.4];
% C = [1  -2  -1   0   3;
%      0   2   1   0   4];
x_ini = [0.3 -0.2 -0.1 0.2]';
% ~~~ system setting ~~~ End

A = [0 0 1 0; 0 0 0 1; -3 1 -0.4 0.1; 1 -1 0.1 -0.1];
B = [0 0; 0 0; 1 0; 0 1];
Bd = B;
C = [1 0 0 0; 0 1 0 0];

[n,r] = size(B);
p = size(C,1);

% ~~~ Controller design
Qc = 5e6*eye(p);
Rc = 1e0*eye(r);
[Kc,Pc] = lqr(A, B, C'*Qc*C, Rc); % Tracker --> C'*Qc*C
Ec = -inv(Rc)*B'*inv(A - B*Kc)'*C'*Qc;
eigOfA = eig(A-B*Kc)
% ~~~ Controller design ~~~ End

Ao = A;
Bo = C';
Qo = 1e8*eye(n);
Ro = 1e0*eye(r);
[Ko,~] = lqr(A', C', Qo, Ro); % Regulator --> Qo
Ko = Ko';
eigA_LC = eig(A-Ko*C)
i = 0; k = 0;
for t_ = 0:Tf:Tend;
    i = i + 1;
    tf(i) = t_;
    if(mod(i-1,sample) == 0)
        k = k + 1;
        tf_plot(k) = t_;
    end
    ref(:, i) = ref_fun(t_)';
    refStar(:, i) = ref_fun(t_ + Tf);
end

r_in = [tf ;refStar]';


[t, x, yc, yco, xc, xco, uc] = sim('model_conti', tf);

% yc
for i = 1:2
figure('Name', ['yc' num2str(i) 'vs. yco' num2str(i) 'vs. r' num2str(i)]);
plot(tf_plot, yc(1:sample:end , i), 'r'); hold on;  plot(tf_plot, yco(1:sample:end , i), 'b'); plot(tf_plot, ref(i, 1:sample:end), 'k'); hold off;
legend({['$y_{c' num2str(i) '}$'] , ['$\hat{y}_{co' num2str(i) '}$'] , ['$r_' num2str(i) '$']} , 'Interpreter','latex', 'fontsize',25)
title({['$y_{c' num2str(i) '}$ $ vs.$ $ \hat{y}_{co' num2str(i) '}$ $vs.$ $r_' num2str(i) ' $']} , 'Interpreter','latex', 'fontsize',25)
h=gca;
set(gca,'FontSize',20,'XTick',h.XTick,'YTick',h.YTick,'YLim',h.YLim)
end

% xc
for i = 1:5
figure('Name', ['xc' num2str(i) ' vs. xco' num2str(i)]);
plot(tf_plot, xc(1:sample:end, i)); hold on; plot(tf_plot, xco(1:sample:end , i)); hold off;
legend({['$x_{c' num2str(i) '}$'] , ['$\hat{x}_{co' num2str(i) '}$']} , 'Interpreter','latex', 'fontsize',25)
title({['$x_{c' num2str(i) '}$ $ vs.$ $ \hat{x}_{co' num2str(i) '}$']} , 'Interpreter','latex', 'fontsize',25)
h=gca;
set(gca,'FontSize',20,'XTick',h.XTick,'YTick',h.YTick,'YLim',h.YLim)
end

% % uc
% for i = 1:2
% figure('Name', ['uc' num2str(i)]);
% plot(tf_plot, uc(1:sample:end, i))
% legend({['u_{c' num2str(i) '}']} ,'fontsize',25)
% title({['$u_{c' num2str(i) '}$']} , 'Interpreter','latex', 'fontsize',25)
% h=gca;
% set(gca,'FontSize',20,'XTick',h.XTick,'YTick',h.YTick,'YLim',h.YLim)
% end

for i = 1:5
figure('Name', ['xc' num2str(i) ' vs. xco' num2str(i)]);
plot(tf, xc(:, i)); hold on; plot(tf, xco(: , i)); hold off;
legend({['$x_{c' num2str(i) '}$'] , ['$\hat{x}_{co' num2str(i) '}$']} , 'Interpreter','latex', 'fontsize',25)
title({['$x_{c' num2str(i) '}$ $ vs.$ $ \hat{x}_{co' num2str(i) '}$']} , 'Interpreter','latex', 'fontsize',25)
h=gca;
set(gca,'FontSize',20,'XTick',h.XTick,'YTick',h.YTick,'YLim',h.YLim)
end

% error of y & ref
err_yr = yc - ref';
for i = 1:2
figure('Name', ['yc' num2str(i) '- r' num2str(i)]);
plot(tf_plot, err_yr(1:sample:end , i), 'b');
% legend({['$y_{' num2str(i) '}$'] , ['$\hat{d}_{e' num2str(i) '}$']} , 'Interpreter','latex', 'fontsize',15)
% title({['$d_{e' num2str(i) '}$ $ vs.$ $ \hat{d}_{e' num2str(i) '}$' ]} , 'Interpreter','latex', 'fontsize',15)
% h=gca;
% set(gca,'FontSize',20,'XTick',h.XTick,'YTick',h.YTick,'YLim',h.YLim)
end

% xco = xco(1:sample:end, :);
% save D:\MATLAB\Multivariable_Control\HW8\HW7_xco 'xco'