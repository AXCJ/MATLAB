% <<< Digital Redesign >>> %
clc; clear all; close all;

Tf = 0.01; % sampling period
Ts = 0.5;  % sampling time
Tend = 10; % final simulation time
sample_div = fix(Ts/Tf);

% ~~~ System setting
A = [  0  ,   1   ;
     -1.2 , -0.8 ];
 
B = [ 0  ;
      1 ];
  
C = [ 0.2 , 0.5 ];

x_ini = [ 1  ; -1 ];
% ~~~ System setting ~~~ End  

[n,r] = size(B);
p = size(C,1); % the size of the first dimension of C
% n : the number of states  of plant
% r : the number of inputs  of plant
% p : the number of outputs of plant

% ~~~ Controller design
Rc = 1e0*eye(r);
Qc = 1e4*eye(p); % Tracker --> the size of Qc is p by p
[Kc,Pc] = lqr(A, B, C'*Qc*C, Rc); % Tracker --> C'*Qc*C
Ec = -inv(Rc)*B'*inv(A - B*Kc)'*C'*Qc;

[G,H] = c2d(A, B, Ts); % Convert from continuous- to discrete-time models
Kd = inv(eye(r) + Kc*H)*Kc*G;
Ed = inv(eye(r) + Kc*H)*Ec;
% ~~~ Controller design ~~~ End

x_c(:,1) = x_ini;
x_d(:,1) = x_c(:,1);
ii = 0; % continuous-time index
kk = 0; %   discrete-time index

% ~~~ Perform simulation
  for t_ = 0 : Tf : Tend
      ii = ii + 1;
      tf(:,ii) = t_; % continuous timespan
      
      ref_c(:,ii) = ref_fun(t_);   
      u_c(:, ii) = Ec*ref_c(:,ii) - Kc*x_c(:, ii); % Tracker --> the reference input isn't equal to zero
      u(:, ii) = u_c(:, ii);
      options = odeset('RelTol',1e-7);
      tspan = [ t_ , t_+Tf ]; % simulation interval
      [Tc,Xc] = ode45(@(tt, xx) plant1(tt, xx, u(:, ii), A, B), tspan, x_c(:,ii), options);
      y_c(:, ii) = C*x_c(:, ii);
   
      % ~~~  Signal sample
      if ( mod(ii-1, sample_div) == 0 )
          kk = kk + 1;
          ts(:, kk) = tf(:, ii); %  discrete timespan
          xd_s(:, kk) = x_d(:, ii);
          ref_d(:, kk) = ref_c(:,ii);
          ref_star(:, kk) = ref_fun(t_ + Ts); % ref_star(:, kk) is equal to ref_d(:,kk+1)
          u_d(:, kk) = Ed*ref_star(:,kk) - Kd*xd_s(:, kk);
      end
      % ~~~ Signal sample ~~~ End
      
      ud_zoh(:, ii) = u_d(:, kk); % actuating signal hold
    
      [Td,Xd] = ode45(@(tt, xx) plant1(tt, xx, ud_zoh(:, ii), A, B), tspan, x_d(:,ii), options);
      y_d(:, ii) = C*x_d(:, ii); 
   
      if(t_ >= Tend)
          break;
      end
    
      x_c(:,ii+1) = Xc(end,:)';
      x_d(:,ii+1) = Xd(end,:)';
      
  end
% ~~~ Perform simulation ~~~ End

% ~~~ Plot fig. ~~~ %
figure(1)
plot(tf, ref_c, '-', tf, y_c, ':', tf, y_d, '-.', 'LineWidth',1.5)
axis([0 Tend -0.5 2.2])
grid on
title('\bf\fontsize{14}\fontname{Cambria} The time responses of \Gamma(t) , y_{c}(t) and y_{d}(t)')
xlabel('\bf\fontsize{10} Time (sec.)')
ylabel('\bf\fontsize{10} Amplitude')
legend('\bf \Gamma','\bf y_{c}','\bf y_{d}')

figure(2)
plot(tf, x_c(1,:), 'b-', tf, x_d(1,:), 'r-.', 'LineWidth',1.5)
grid on
title('\bf\fontsize{14}\fontname{Cambria} The time responses of x_{1}(t)')
xlabel('\bf\fontsize{10} Time (sec.)')
ylabel('\bf\fontsize{10} Amplitude')
legend('\bf x_{c1}','\bf x_{d1}')

figure(3)
plot(tf, x_c(2,:), 'b-', tf, x_d(2,:), 'r-.', 'LineWidth',1.5)
grid on
title('\bf\fontsize{14}\fontname{Cambria} The time responses of x_{2}(t)')
xlabel('\bf\fontsize{10} Time (sec.)')
ylabel('\bf\fontsize{10} Amplitude')
legend('\bf x_{c2}','\bf x_{d2}')

figure(4)
plot(tf, u_c, 'b-', tf, ud_zoh,'r-.', 'LineWidth',1.5)
grid on
title('\bf\fontsize{14}\fontname{Cambria} The time responses of actuating signal')
xlabel('\bf\fontsize{10} Time (sec.)')
ylabel('\bf\fontsize{10} Amplitude')
legend('\bf u_{c}','\bf u_{d,zoh}')

figure(5)
plot(tf, ref_c, '-', ts, ref_d,':', ts, ref_star,'-.', 'LineWidth',1.5)
axis([0 Tend -0.5 2.2])
grid on
title('\bf\fontsize{14}\fontname{Cambria} The time responses of reference input')
xlabel('\bf\fontsize{10} Time (sec.)')
ylabel('\bf\fontsize{10} Amplitude')
legend('\bf \Gamma_{c}','\bf \Gamma_{d}','\bf \Gamma^{*}')
