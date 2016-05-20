
clc;clear all;close all;
Linewidth = 1;
scrsz = get(groot,'ScreenSize');
Tf = 5e-3; % sampling period
Tend = 3; % final simulation time
Ts = 0.01; % sampling time
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
% n : the number of states of plant
% r : the number of inputs of plant
% p : the number of outputs of plant

% ~~~ Controller design
Rc = 1e0*eye(r);
Qc = 1e6*eye(p); % Tracker --> the size of Qc is p by p
[Kc,Pc] = lqr(A, B, C'*Qc*C, Rc); % Tracker --> C'*Q*C
Ec = -inv(Rc)*B'*inv(A - B*Kc)'*C'*Qc;
% ~~~ Controller design ~~~ End

xc(:,1) = x_ini;
ii = 0; % continuous-time index
kk = 0; %   discrete-time index
% ~~~ Perform simulation
for t_ = 0 : Tf : Tend
    ii = ii + 1;
    ref(: , ii) = ref_fun(t_)'; % get reference input
    ref_star(: , ii) = ref_fun(t_+Tf)'; % ref_star(:, ii) is equal to ref(: , ii+1)
    tf(: , ii) = t_; % continuous-time span
end

%%%~~~ another solution ~~~%%%
% Tf = 1e-2;
Td_seq = 0 : Tf : Tend;
% figure('Position',[9 scrsz(4)/1.7 scrsz(3)/3 scrsz(4)/3]) % Left_Top of screen
options = odeset('RelTol',1e-7);
[T,x_dott]=ode45(@(t,x) plant2(t,x,Kc,Ec,A,B,ref,Tf) , Td_seq , x_ini , options);
% Y_ode=C *x_dott';




