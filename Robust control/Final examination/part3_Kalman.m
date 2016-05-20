clc; clear all; close all;
load('laserAngle.mat')
Tf = 1; % sampling period
Tend = 850; % final simulation time
xi = [0.2 0.1 0.1 0.3 0.3];
wc = [0.8 0.5 1.2 0.5 1.2];
sys_idx = 1;
num = [1];

den = [1 2*xi(sys_idx)*wc(sys_idx) wc(sys_idx)^2];
[A,B,C,D] = tf2ss(num,den);
sys = ss(A,B,C,D);

x_ini = [ 0  ; 0 ];
% ~~~ System setting ~~~ End  

[n,r] = size(B);
p = size(C,1); % the size of the first dimension of C 
% n : the number of states  of plant
% r : the number of inputs  of plant
% p : the number of outputs of plant

% ~~~ Observer setting
[Ao,Bo,Co,Do] = tf2ss(num,den);

xo_ini = [ 1  ; 1 ];
% ~~~ Observer setting ~~~ End  

[no,ro] = size(Bo);
po = size(Co,1);
% no : the number of states  of observer
% ro : the number of inputs  of observer
% po : the number of outputs of observer

% ~~~ Observer design
Ro = 1e0*eye(ro);
Qo = 5e6*eye(no); % Regulator --> the size of Qo is no by no

[Ko,Po] = lqr(A', C', Qo, Ro); % Regulator --> Qo
Ko = Ko'; % (A, B ,Kc) <---> (A', C', Ko') for continuous-time system
% Design Ko such that ( A - Ko*C ) is asymptotically stable
% ~~~ Observer design ~~~ End
L=Ko;
% --- Controller design
R = 1e0*eye(ro);
Q = 5e6*eye(po); % Tracker --> the size of Q is po by po

[Kc,Pc] = lqr(Ao, Bo, Co'*Q*Co, Ro); % Tracker --> C'*Q*C
% Ec = -inv(R)*Bo'*inv(Ao - Bo*Kc)'*Co'*Q;
Ec = -inv(Co*inv(Ao-Bo*Kc)*Bo);
% --- Controller design --- end
sys_ob = ss(A-B*Kc-Ko*C,Ko,-Kc,[]);
[p,z] = pzmap(sys_ob);
x_c(:,1) = x_ini;
x_o(:,1) = xo_ini;

% ref = ones(1,length(0:Tf:Tend));
simin = [0:1:length(laserAngle(2:end))-1; laserAngle(2:end)]';

sim('model',[0:1:850-1])
% [ y_c,y_o,tf ] = solver(xi, wc, sys_idx, 1, Tend, laserAngle, 0.01);
load('y_c.mat')
yyc = y_c(1:100:end);
% fprintf('xi:%f   wc:%f  SettlingTime: %f  steady state error:%f\n',xi(sys_idx),wc(sys_idx),info.SettlingTime,y_o(end)-1)

figure
plot(simt, yout, 'b', 0:1:length(laserAngle)-1, laserAngle,'r-.', 'LineWidth',1.1)
title({'$\hat{\theta} $  $ $ vs.$ $ $ \psi$'},'Interpreter','latex', 'fontsize',30)
set(gca,'FontSize',20)
legend({'$\hat{\theta}$','$\psi$'},'Interpreter','latex', 'fontsize',30)
figure
plot(laserAngle(1:end-1)-yout')
title({'Error'},'Interpreter','latex', 'fontsize',30)
set(gca,'FontSize',20)
x23 = x(1:end-1);
y23 = y(1:end-1);
for i = 1:length(simt)
    sL = yout(i);
    s1 = angle1hat(1,i);
    s2 = angle2hat(1,i);
    x12(i) = 20*tan(s2)/(tan(sL)-tan(s2));
    y12(i) = tan(sL)*x12(i);
    x13(i) = -50*tan(s1)/(tan(sL)-tan(s1));
    y13(i) = tan(sL)*x12(i);
    arrrrea(i) = 0.5*abs(det([1 x13(i) y13(i);1 x12(i) y12(i) ;1 x23(i) y23(i)]));
    area(i) = 0.5*abs((x13(i)*y12(i)+x12(i)*y23(i)+x23(i)*y13(i))-(x12(i)*y13(i)+x23(i)*y12(i)+x13(i)*y23(i)));
end

figure
plot(area)
title({'Triangle'},'Interpreter','latex', 'fontsize',30)
% set(gca,'FontSize',20)
