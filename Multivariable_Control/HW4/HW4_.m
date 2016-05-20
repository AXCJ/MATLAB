
clc;clear all;close all;
delete u.mat
Tf = 5e-5; % sampling period
Tend = 3; % final simulation time
Ts = 0.01; % sampling time

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
Qc = 1e6*eye(p); % Tracker --> the size of Qc is p by p
[Kc,Pc] = lqr(A, B, C'*Qc*C, Rc); % Tracker --> C'*Q*C
Ec = -inv(Rc)*B'*inv(A - B*Kc)'*C'*Qc;
% ~~~ Controller design ~~~ End

xc(:,1) = x_ini;
ii = 0; % continuous-time index
% ~~~ Perform simulation
for t_ = 0 : Tf : Tend+Tf
    ii = ii + 1;
    ref(: , ii) = ref_fun(t_)'; % get reference input
    if(t_ >=Tend+Tf)
        break;
    end
    tf(: , ii) = t_; % continuous-time span
end

Td_seq = 0 : Tf : Tend;


options = odeset('RelTol',1e-7);
[T,X]=ode45(@(tt,xx) plant2(tt,xx,Kc,Ec,A,B,ref) , Td_seq , x_ini , options);
Y_ode=C *X';
test_u = Ec*ref(:,1:end-1)- Kc*X';


% ~~~ Plot fig. ~~~ %
scrsz = get(groot,'ScreenSize');
%%% [Left bottom width height]

% figure('Position',[scrsz(3)/3 scrsz(4)/1.7 scrsz(3)/3 scrsz(4)/3])
% figure('Position',[9 scrsz(4)/1.7 scrsz(3)/3 scrsz(4)/3])
% figure('Position',[scrsz(3)/1.5 scrsz(4)/1.7 scrsz(3)/3 scrsz(4)/3])
% figure('Position',[9 scrsz(4)/1.7/4 scrsz(3)/3 scrsz(4)/3])

% figure('Position',[9 scrsz(4)/1.7 scrsz(3)/3 scrsz(4)/3])  % top of left of screen
% plot(tf , ref(1 , :) , tf , yc(1 , :) , '-.' , 'LineWidth' ,1 )
% legend({'r_{1}','y_{c1}'},'fontsize',15)
% xlabel('t(s)','fontsize',20)
% title('y_{c1} vs. r_{1}','fontsize',30)
% hh=gca;
% set(hh,'FontSize',20)

figure('Position',[9 scrsz(4)/1.7 scrsz(3)/3 scrsz(4)/3])  % top of left of screen
plot(Td_seq , ref(1 , :) , Td_seq , Y_ode(1 , :) , '-.' , 'LineWidth' ,1 )
legend({'r_{1}','y_{c1}'},'fontsize',15)
xlabel('t(s)','fontsize',20)
title('y_{c1} vs. r_{1}','fontsize',30)
hh=gca;
set(hh,'FontSize',20)

figure('Position',[scrsz(3)/1.5 scrsz(4)/1.7 scrsz(3)/3 scrsz(4)/3])
plot(Td_seq , ref(2 , 1:301), Td_seq , Y_ode(2 , :) , '-.' , 'LineWidth' , 1)
legend({'r_{2}','y_{c2}'},'fontsize',15)
xlabel('t(s)','fontsize',20)
title('y_{c2} vs. r_{2}','fontsize',30)
hh=gca;
set(hh,'FontSize',20)

% ylabel('y_{c2} vs. r_{c2}')
% set(get(gca,'YLabel'),'Rotation',0); % Rotate YLabel
% shift YLabel
% set(get(gca,'YLabel'),'Position',get(get(gca,'YLabel'),'Position') - [0.090 -.01 0])

load -ascii u.mat

figure('Position',[9 scrsz(4)/1.7 scrsz(3)/3 scrsz(4)/3])  % top of left of screen
plot(u(:,1),u(:,2))
hh=gca;
% set(hh,'FontSize',20)
xlabel('t(s)','fontsize',20)
title('u_{1}','fontsize',30)

figure('Position',[scrsz(3)/3 scrsz(4)/1.7 scrsz(3)/3 scrsz(4)/3]) % top of middle of screen
plot(u(:,1),u(:,3))
hh=gca;
set(hh,'FontSize',20)
xlabel('t(s)','fontsize',20)
title('u_{2}','fontsize',30)

sys = ss(A,B,C,0);
figure('Position',[scrsz(3)/3 scrsz(4)/1.7 scrsz(3)/3 scrsz(4)/3]) % top of middle of screen
h1 = pzplot(sys);
h1o = getoptions(h1); % Get options for plot.
h1o.Title.FontSize = 20;
h1o.Title.String = 'Open-loop Pole-Zero Map';
h1o.TickLabel.FontSize = 15;
setoptions(h1,h1o); % Apply options to plot.
syscl = ss(A-B*Kc,B*Ec,C,0);

figure('Position',[scrsz(3)/3 scrsz(4)/1.7 scrsz(3)/3 scrsz(4)/3]) % top of middle of screen
h2 = pzplot(syscl);
h2o = getoptions(h2); % Get options for plot.
h2o.Title.FontSize = 20;
h2o.Title.String = 'Closed-loop Pole-Zero Map';
h2o.TickLabel.FontSize = 15;
setoptions(h2,h2o); % Apply options to plot.

