clear all;clc; close all;
Tf = 0.01;
Ts = 0.02;
sample = Ts/Tf;
Tend = 3;
% ~~~ system setting
A = [ 0.9   -2.1     0.3     0.5     0.9;
      0.7    0.2     0.3       0     0.7;
     -1.3    0.5    -1.1    -2.3    -0.2;
     -0.3    0.8     1.7    -1.2    -0.3;
     -3.5   -4.3    -0.7       0    -8.3];
B = [   1   -0.4;
     -1.7   -1.7;
     -0.2    1.2; 
      0.7    0.1;
      0.9    1.4];
C = [1  -2  -1   0   3;
         0   2   1   0   4];
x_ini = [ 0.3 -0.2 0 -0.1 0.2]';
% ~~~ system setting ~~~ End

[n,r] = size(B);
p = size(C,1);
% ~~~ Controller design
Qc = 5e6*eye(p);
Rc = 1e0*eye(r);
[Kc,Pc] = lqr(A, B, C'*Qc*C, Rc); % Tracker --> C'*Q*C
Ec = -inv(Rc)*B'*inv(A - B*Kc)'*C'*Qc;

[G,H] = c2d(A, B, Ts); % Convert from continuous- to discrete-time models
Kd = inv(eye(r) + Kc*H)*Kc*G;
Ed = inv(eye(r) + Kc*H)*Ec;
% ~~~ Controller design ~~~ End

Ao = A;
Bo = C';
Qo = 1e8*eye(n);
Ro = 1e0*eye(r);
[Ko,~] = lqr(A', C', Qo, Ro); % Regulator --> Qo
Ko = Ko';

Lo = (G-eye(n))/A*Ko/(eye(p)+C*(G-eye(n))/A*Ko);
G_tilde = (eye(n)-Lo*C)*G;
H_tilde = (eye(n)-Lo*C)*H;

i = 0; k = 0;
for t_ = 0 : Tf : Tend;
    i = i + 1;
    tf(i) = t_;
    ref(:, i) = ref_fun(t_)';
    if(mod(i-1, sample) == 0)
        k = k + 1;
        ts(k) = t_;
        refStar(:, k) = ref_fun(t_ + Ts); % r(KsTs+Ts)
        
        rtest(:,k) = ref_fun(t_)';
    end
    rtestHold(:,i) = rtest(:,k);
end
figure
plot(tf, ref(1,:), tf, rtestHold(1, :))
legend('r', 'HOLD')
r_in = [ts ;refStar]'; % simulation input
% the number of left-hand side arguments must be 2 (for T,X) plus number of root-level outport blocks
[tt, xx, yd, ydo, xd, xdo, ud] = sim('model_discre', tf);

% y
for i = 1:2
figure('Name', ['yd' num2str(i) ' vs. ydo' num2str(i) 'vs. r' num2str(i)]);
plot(tf, yd(: , i), 'r'); hold on;  plot(tf, ydo(: , i), 'b'); plot(tf, ref(i, :), 'k'); hold off;
legend({['$y_{d' num2str(i) '}$'] , ['$\hat{y}_{do' num2str(i) '}$'] , ['$r_' num2str(i) '$']} , 'Interpreter','latex', 'fontsize',25)
title({['$y_{d' num2str(i) '}$ $ vs.$ $ \hat{y}_{do' num2str(i) '}$ $vs.$ $r_' num2str(i) ' $']} , 'Interpreter','latex', 'fontsize',25)
h=gca;
set(gca,'FontSize',20,'XTick',h.XTick,'YTick',h.YTick,'YLim',h.YLim)
end

% xd vs. xdo
for i = 1:5
figure('Name', ['xd' num2str(i) 'vs. xdo' num2str(i)]);
plot(tf, xd(:, i)); hold on; plot(tf, xdo(: , i)); hold off;
legend({['$x_{d' num2str(i) '}$'] , ['$\hat{x}_{do' num2str(i) '}$']} , 'Interpreter','latex', 'fontsize',25)
title({['$x_{d' num2str(i) '}$ $ vs.$ $ \hat{x}_{do' num2str(i) '}$']} , 'Interpreter','latex', 'fontsize',25)
h=gca;
set(gca,'FontSize',20,'XTick',h.XTick,'YTick',h.YTick,'YLim',h.YLim)
end

% ud
for i = 1:2
figure('Name', ['ud' num2str(i)]);
plot(tf, ud(:, i))
legend({['u_{d' num2str(i) '}']} ,'fontsize',25)
title({['$u_{d' num2str(i) '}$']} , 'Interpreter','latex', 'fontsize',25)
h=gca;
set(gca,'FontSize',20,'XTick',h.XTick,'YTick',h.YTick,'YLim',h.YLim)
end

load('HW7_xco.mat')
% xdo vs xco
for i = 1:5
figure('Name', ['xdo' num2str(i) ' vs. xco' num2str(i)]);
plot(tf, xdo(:, i)); hold on; plot(tf, xco(: , i)); hold off;
legend({['$x_{do' num2str(i) '}$'] , ['$x_{co' num2str(i) '}$']} , 'Interpreter','latex', 'fontsize',25)
title({['$x_{do' num2str(i) '}$ $ vs.$ $ x_{co' num2str(i) '}$']} , 'Interpreter','latex', 'fontsize',25)
h=gca;
set(gca,'FontSize',20,'XTick',h.XTick,'YTick',h.YTick,'YLim',h.YLim)
end
