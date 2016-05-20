% <<< Digital Redesign >>> %
clc;clear;close all;

scrsz = get(groot,'ScreenSize');
load('HW4_Tf_5e-5_xyu.mat')

Tf = 1e-2; % sampling period
Tend = 3; % final simulation time
Ts = 2e-2; % sampling time
sample_div = fix(Ts/Tf);


% ~~~ display
fprintf('sampling period Tf:%2.3f sec\n',Tf);
disp(sprintf('sampling time Ts:%2.3f sec',Ts));
disp(sprintf('sample_div :%2.3f',fix(Ts/Tf)));

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
p = size(C,1); % the size of the first dimension of C
% n : the number of states  of plant
% r : the number of inputs  of plant
% p : the number of outputs of plant

% ~~~ Controller design
Rc = 1e0*eye(r);
Qc = 5e6*eye(p); % Tracker --> the size of Qc is p by p
[Kc,Pc] = lqr(A, B, C'*Qc*C, Rc); % Tracker --> C'*Qc*C
Ec = -inv(Rc)*B'*inv(A - B*Kc)'*C'*Qc;

[G,H] = c2d(A, B, Ts); % Convert from continuous- to discrete-time models
Kd = inv(eye(r) + Kc*H)*Kc*G;
Ed = inv(eye(r) + Kc*H)*Ec;
% ~~~ Controller design ~~~ End

% x_c(:,1) = x_ini;
x_d(:,1) = x_ini;
ii = 0; % continuous-time index
kk = 0; %   discrete-time index
options = odeset('RelTol',1e-7);

% ~~~ Perform simulation
  for t_ = 0 : Tf : Tend
      ii = ii + 1;
      tf(:,ii) = t_; % continuous timespan
      
      ref_c(:,ii) = ref_fun(t_);   
      tspan = [ t_ , t_+Tf ]; % simulation interval
%       u_c(:, ii) = Ec*ref_c(:,ii) - Kc*x_c(:, ii); % Tracker --> the reference input isn't equal to zero
%       u(:, ii) = u_c(:, ii);
%       [Tc,Xc] = ode45(@(tt, xx) plant1(tt, xx, u(:, ii), A, B), tspan, x_c(:,ii), options);
%       y_c(:, ii) = C*x_c(:, ii);
%    
      % ~~~  Signal sample
      y_d(:, ii) = C*x_d(:, ii); 
      if ( mod(ii-1, sample_div) == 0 )
          kk = kk + 1;
          ts(:, kk) = tf(:, ii); %  discrete timespan
          xd_s(:, kk) = x_d(:, ii);
          y_ds(:, kk) = y_d(:, ii);
          ref_d(:, kk) = ref_c(:, ii);
          ref_star(:, kk) = ref_fun(t_ + Ts); % ref_star(:, kk) is equal to ref_d(:,kk+1)
          u_d(:, kk) = Ed*ref_star(:,kk) - Kd*xd_s(:, kk);
      end
      % ~~~ Signal sample ~~~ End
      ud_zoh(:, ii) = u_d(:, kk); % actuating signal hold
 
      [Td,Xd] = ode45(@(tt, xx) plant1(tt, xx, ud_zoh(:, ii), A, B), tspan, x_d(:,ii), options);
   
      if(t_ >= Tend)
          break;
      end
    
%       x_c(:,ii+1) = Xc(end,:)';
      x_d(:,ii+1) = Xd(end,:)';
      
  end
% ~~~ Perform simulation ~~~ End

% ~~~ Plot fig. ~~~ %
figure('Name', 'y_{c1}(t) vs. y_{d1}(t) vs. r_1(t) vs. y_{d1}(KT)' ,'Position',[9 scrsz(4)/1.7 scrsz(3)/3 scrsz(4)/3]) % Left_Top of screen
plot(tf, yc_s(1,:), tf, y_d(1,:), tf, ref_c(1,:), ts, y_ds(1, :),'o', 'LineWidth',1.5)
grid off
title('\bf\fontsize{30}\fontname{Cambria} y_{c1}(t) vs. y_{d1}(t) vs. r_1(t) vs. y_{d1}(KT)')
xlabel('\bf\fontsize{20} Time (sec.)')
ylabel('\bf\fontsize{20} Amplitude')
legend({'\bfy_{c1}(t)' , '\bfy_{d1}(t)' , '\bfr_1(t)' ,  '\bfy_{d1}(KT)'} , 'fontsize',15)
set(gca,'FontSize',20)

figure('Name', 'y_{c2}(t) vs. y_{d2}(t) vs. r_2(t) vs. y_{d2}(KT)' ,'Position',[scrsz(3)/3 scrsz(4)/1.7 scrsz(3)/3 scrsz(4)/3]) % Middle_Top of screen
plot(tf, yc_s(2,:), tf, y_d(2,:), tf, ref_c(2,:), ts, y_ds(2, :),'o', 'LineWidth',1.5)
grid off
title('\bf\fontsize{30}\fontname{Cambria} y_{c2}(t) vs. y_{d2}(t) vs. r_2(t) vs. y_{d2}(KT)')
xlabel('\bf\fontsize{20} Time (sec.)')
ylabel('\bf\fontsize{20} Amplitude')
legend({'\bfy_{c2}(t)' , '\bfy_{d2}(t)' , '\bfr_2(t)' ,  '\bfy_{d2}(KT)'} , 'fontsize',15)
set(gca,'FontSize',20)

for j = 5:-1:1 % plot x
    figure('Name', ['x' num2str(j)] ,'Position',[9+(j-1)*100 scrsz(4)/1.7/4 scrsz(3)/3 scrsz(4)/3]) % Left_Bottom of screen
    plot(tf, xc_s(j, :), tf, x_d(j, :), 'LineWidth',1.5)
    title(['\bf\fontsize{30\fontname{Cambria} x_{c' num2str(j) '}(t) vs. x_{d' num2str(j) '}(t)'])
    xlabel('\bf\fontsize{20} Time (sec.)')
    ylabel('\bf\fontsize{20} Amplitude')
    legend({['\bf x_{c' num2str(j) '}(t)'],['\bf x_{d' num2str(j) '}(t)'] }, 'fontsize',15)
    set(gca,'FontSize',20)
end

for j = 2:-1:1 % plot u
    figure('Name', ['u' num2str(j)] ,'Position',[scrsz(3)/1.5-(j-1)*100 scrsz(4)/1.7/4 scrsz(3)/3 scrsz(4)/3]) % Right_Bottom of screen
    plot(tf, uc_s(j, :), tf, ud_zoh(j, :), 'LineWidth',1.5)
    title(['\bf\fontsize{30\fontname{Cambria} u_{c' num2str(j) '}(t) vs. u_{d' num2str(j) '}(t)'])
    xlabel('\bf\fontsize{20} Time (sec.)')
    ylabel('\bf\fontsize{20} Amplitude')
    legend({['\bf u_{c' num2str(j) '}(t)'],['\bf u_{d' num2str(j) '}(t)'] }, 'fontsize',15)
    set(gca,'FontSize',20)
end

%% part2
% let   Kd = Kc
%        Ed = Ec
% re-design
clear all;
scrsz = get(groot,'ScreenSize');

Tf = 1e-2; % sampling period
Tend = 3; % final simulation time
Ts = 2e-2; % sampling time
sample_div = fix(Ts/Tf);

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
x_ini = [ 0.3;
             -0.2;
                0;
             -0.1;
              0.2];
% ~~~ system setting ~~~ End

[n,r] = size(B);
p = size(C,1); % the size of the first dimension of C
% n : the number of states  of plant
% r : the number of inputs  of plant
% p : the number of outputs of plant

% ~~~ Controller design
Rc = 1e0*eye(r);
Qc = 5e6*eye(p); % Tracker --> the size of Qc is p by p
[Kc,Pc] = lqr(A, B, C'*Qc*C, Rc); % Tracker --> C'*Qc*C
Ec = -inv(Rc)*B'*inv(A - B*Kc)'*C'*Qc;

% ~~~ Controller design ~~~ End

x_d(:,1) = x_ini;
ii = 0; % continuous-time index
kk = 0; %   discrete-time index
options = odeset('RelTol',1e-7);
Kd = Kc;
Ed = Ec;
  for t_ = 0 : Tf : Tend
      ii = ii + 1;
      tf(:,ii) = t_; % continuous timespan
      
      ref_c(:,ii) = ref_fun(t_);   
      tspan = [ t_ , t_+Tf ]; % simulation interval

      % ~~~  Signal sample
      y_d(:, ii) = C*x_d(:, ii); 
      if ( mod(ii-1, sample_div) == 0 )
          kk = kk + 1;
          ts(:, kk) = tf(:, ii); %  discrete timespan
          xd_s(:, kk) = x_d(:, ii);
          y_ds(:, kk) = y_d(:, ii);
          ref_d(:, kk) = ref_c(:, ii);
          ref_star(:, kk) = ref_fun(t_ + Ts); % ref_star(:, kk) is equal to ref_d(:,kk+1)
          u_d(:, kk) = Ed*ref_star(:,kk) - Kd*xd_s(:, kk);
      end
      % ~~~ Signal sample ~~~ End
      ud_zoh(:, ii) = u_d(:, kk); % actuating signal hold
    
      [Td,Xd] = ode45(@(tt, xx) plant1(tt, xx, ud_zoh(:, ii), A, B), tspan, x_d(:,ii), options);
   
      if(t_ >= Tend)
          break;
      end
    
      x_d(:,ii+1) = Xd(end,:)';
      
  end
  
% figure('Name', 'y_{d1}(t) vs. r_1(t)' ,'Position',[9 scrsz(4)/1.7 scrsz(3)/3 scrsz(4)/3]) % Left_Top of screen
% plot(tf, y_d(1,:), tf, ref_c(1,:), 'LineWidth',1.5)
% grid off
% title('\bf\fontsize{30}\fontname{Cambria} y_{d1}(t) vs. r_1(t)')
% xlabel('\bf\fontsize{20} Time (sec.)')
% ylabel('\bf\fontsize{20} Amplitude')
% legend({'\bfy_{d1}(t)' , '\bfr_1(t)'} , 'fontsize',15)
% set(gca,'FontSize',20)

for j = 2:-1:1 % plot u
    figure('Name', ['y' num2str(j)] ,'Position',[scrsz(3)/1.5-(j-1)*100 scrsz(4)/1.7 scrsz(3)/3 scrsz(4)/3]) % Right_Bottom of screen
    plot(tf, y_d(j, :), tf, ref_c(j, :), 'LineWidth',1.5)
    title(['\bf\fontsize{30\fontname{Cambria} y_{d' num2str(j) '}(t) vs. r_{' num2str(j) '}(t)'])
    xlabel('\bf\fontsize{20} Time (sec.)')
    ylabel('\bf\fontsize{20} Amplitude')
    legend({['\bf y_{c' num2str(j) '}(t)'],['\bf r_{' num2str(j) '}(t)'] }, 'fontsize',15)
    set(gca,'FontSize',20)
end