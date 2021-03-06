% <<< Continuous-time system plus observer simulation >>> %
% ### Design procedure : Observer --> Controller ### %
clc; clear all; close all;

Tf = 0.01; % sampling period
Tend = 5; % final simulation time

% ~~~ System setting
num = [2 -4];
den = [1 -5 4];
[A,B,C,D] = tf2ss(num,den);
sys = ss(A,B,C,D);

x_ini = [ 1  ; -1 ];
% ~~~ System setting ~~~ End  

[n,r] = size(B);
p = size(C,1); % the size of the first dimension of C 
% n : the number of states  of plant
% r : the number of inputs  of plant
% p : the number of outputs of plant

% ~~~ Observer setting
[Ao,Bo,Co,Do] = tf2ss(num,den);

xo_ini = [ 0  ; 0 ];
% ~~~ Observer setting ~~~ End  

[no,ro] = size(Bo);
po = size(Co,1);
% no : the number of states  of observer
% ro : the number of inputs  of observer
% po : the number of outputs of observer

% ~~~ Observer design
Ro = 1e0*eye(ro);
Qo = 1e1*eye(no); % Regulator --> the size of Qo is no by no

[Ko,Po] = lqr(A', C', Qo, Ro); % Regulator --> Qo
Ko = Ko'; % (A, B ,Kc) <---> (A', C', Ko') for continuous-time system
% Design Ko such that ( A - Ko*C ) is asymptotically stable
% ~~~ Observer design ~~~ End

% --- Controller design
R = 1e0*eye(ro);
Q = 1e1*eye(po); % Tracker --> the size of Q is po by po

[Kc,Pc] = lqr(Ao, Bo, Co'*Q*Co, Ro); % Tracker --> C'*Q*C
Ec = -inv(R)*Bo'*inv(Ao - Bo*Kc)'*Co'*Q;
% --- Controller design --- end
sys_ob = ss(A-B*Kc-Ko*C,Ko,-Kc,[]);
[p,z] = pzmap(sys_ob)
x_c(:,1) = x_ini;
x_o(:,1) = xo_ini;
ii = 0;

% ~~~ Perform simulation
  for t_ = 0 : Tf : Tend
      ii = ii + 1;
      tf(:,ii) = t_; % continuous timespan
      ref_c(:,ii) = ref_fun(t_);
      u_c(:, ii) = Ec*ref_c(:,ii) - Kc*x_o(:,ii); % Tracker
      options = odeset('RelTol',1e-7);
      tspan = [ t_ , t_+Tf ];
      
      [Tc,Xc] = ode45(@(tt, xx) plant1(tt, xx, u_c(:, ii), A, B), tspan, x_c(:,ii), options);
      y_c(:, ii) = C*x_c(:, ii);
    
      [To,Xo] = ode45(@(tt, xx) observer1(tt, xx, u_c(:, ii), Ao, Bo, Co, Ko, y_c(:, ii)), tspan, x_o(:,ii), options);
      y_o(:, ii) = Co*x_o(:, ii);
     
      if(t_ >= Tend)
          break;
      end
    
      x_c(:,ii+1) = Xc(end,:)';
      x_o(:,ii+1) = Xo(end,:)';
    
  end
% ~~~ Perform simulation ~~~ End

% --- Plot fig. --- %
figure(1)
plot(tf, y_c, 'b', tf, y_o,'r-.', 'LineWidth',1.1)
grid off
title('\bf\fontsize{14}\fontname{Cambria} The time responses of y_{c}(t) and y_{o}(t)')
xlabel('\bf\fontsize{10} Time (sec.)')
ylabel('\bf\fontsize{10} Amplitude')
legend('\bf y_{c}','\bf y_{o}')

figure(2)
plot(tf, x_c(1,:), 'b-', tf, x_o(1,:),'r-.', 'LineWidth',1.5)
title('\bf\fontsize{14}\fontname{Cambria} The time responses of x_{1}(t)')
xlabel('\bf\fontsize{10} Time (sec.)')
ylabel('\bf\fontsize{10} Amplitude')
legend('\bf x_{c1}','\bf x_{o1}')

figure(3)
plot(tf, x_c(2,:), 'b-', tf, x_o(2,:),'r-.', 'LineWidth',1.5)
title('\bf\fontsize{14}\fontname{Cambria} The time responses of x_{2}(t)')
xlabel('\bf\fontsize{10} Time (sec.)')
ylabel('\bf\fontsize{10} Amplitude')
legend('\bf x_{c2}','\bf x_{o2}')

figure(4)
plot(tf, u_c, 'b-', 'LineWidth',1.5)
grid on
title('\bf\fontsize{14}\fontname{Cambria} The time response of control signal')
xlabel('\bf\fontsize{10} Time (sec.)')
ylabel('\bf\fontsize{10} Amplitude')
