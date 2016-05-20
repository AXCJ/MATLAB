% <<<  Key : 先將對應之方塊圖畫出，那麼就能得知整個程式流程 >>> %
% ### Design procedure : Observer --> Controller ### %
clc; clear all; close all;

Tf = 0.01; % sampling period
Tend = 10; % final simulation time
Ts = 0.1;  % sampling time
sample_div = fix(Ts/Tf); 

% --- system setting
A = [  0  ,   1   ;
     -1.2 , -0.8 ];
 
B = [ 0  ;
      1 ];
  
C = [ 0.2 , 0.5 ];

x_ini = [ 1  ;
         -1 ];
% --- system setting --- end   
 
[n,r] = size(B);
p = size(C,1); % the size of the first dimension of C    
% n : the number of states  of plant
% r : the number of inputs  of plant
% p : the number of outputs of plant

% --- Observer design
Ao = [  0  ,   1   ;
      -1.2 , -0.8 ];
 
Bo = [ 0  ;
       1 ];
  
Co = [0.2 , 0.5];

xo_ini = [ 0  ;
           0 ];

[Go,Ho] = c2d(Ao, Bo, Ts); % Convert from continuous- to discrete-time models 

[no,ro] = size(Bo);
po = size(Co,1);
% no : the number of states  of observer
% ro : the number of inputs  of observer
% po : the number of outputs of observer

Ro = 1e0*eye(ro);
Qo = 1e4*eye(no); % Regulator --> the size of Qo is n by n

Ko = dlqr(Go', Co', Qo, Ro); % LQR for discrete-time system
Ko = Ko'; % (A, B ,Kc) <---> (A', C', Ko') for continuous-time system
% design Ko such that ( Go - Ko*Co ) is asymptotically stable
% --- Observer design --- end

x_d(:, 1) = x_ini;
x_o(:, 1) = xo_ini;
ii = 0; % continuous-time index
kk = 0; %   discrete-time index
for t_ = 0:Tf:Tend
    ii = ii + 1;
    tf(:,ii) = t_; % continuous-time 
    ref_c(:,ii) = ref_fun(t_);
    y_d(:, ii) = C*x_d(:, ii); 
    % ---  signal sample
    if ( mod(ii-1, sample_div) == 0 )
        kk = kk + 1;
        ts(:, kk) = tf(:, ii); %  discrete-time
        ref_d(:,kk) = ref_c(:,ii);
        y_s(:, kk) = y_d(:, ii);
        y_o(:, kk) = Co*x_o(:, kk);
        ud(:, kk) = ref_d(:,kk);
    end
    % --- signal sample --- end
    
    ud_zoh(:, ii) = ud(:, kk); % actuating signal hold
    [Td,Xd] = ode45(@(tt, xx) plant1(tt, xx, ud_zoh(:, ii), A, B),[t_,t_+Tf], x_d(:,ii));
    
    
    if(t_ >= Tend)
        break;
    end
    
    if ( mod(ii-1, sample_div) == 0 )
        x_o(:, kk+1) = Go*x_o(:, kk) + Ho*ud(:, kk) + Ko*(y_s(:, kk) - y_o(:, kk));
    end
        x_d(:, ii+1) = Xd(end, :)'; 
end

% --- Plot fig. --- %
figure(1)
plot(tf, ud_zoh, '-', tf, y_d, ':', ts, y_o,  '-.', 'LineWidth',1.5)
axis([0 Tend -0.8 2.2])
grid on
title('\bf\fontsize{14}\fontname{Cambria} The time responses of \Gamma(t) , y_{d}(t) and y_{o}(t)')
xlabel('\bf\fontsize{10} Time (sec.)')
ylabel('\bf\fontsize{10} Amplitude')
legend('\bf u_{d,zoh}','\bf y_{d}','\bf y_{o}')

figure(2)
plot(tf, x_d(1,:), 'b-', ts, x_o(1,:), 'r-.', 'LineWidth',1.5)
grid on
title('\bf\fontsize{14}\fontname{Cambria} The time responses of x_{d1}(t) and x_{o1}(t)')
xlabel('\bf\fontsize{10} Time (sec.)')
ylabel('\bf\fontsize{10} State Variables x_{1}(t)')
legend('\bf x_{d1}','\bf x_{o1}')

figure(3)
plot(tf, x_d(2,:), 'b-', ts, x_o(2,:), 'r-.', 'LineWidth',1.5)
grid on
title('\bf\fontsize{14}\fontname{Cambria} The time responses of x_{d2}(t) and x_{o2}(t)')
xlabel('\bf\fontsize{10} Time (sec.)')
ylabel('\bf\fontsize{10} State Variables x_{2}(t)')
legend('\bf x_{d2}','\bf x_{o2}')

