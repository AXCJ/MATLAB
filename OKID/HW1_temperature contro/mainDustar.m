clc; clear; close all;
%% system model set
load OKIDData
load OKID_esys_13

Ts = 0.1;
Tend = 60.1;
Sam_t = 0:Ts:Tend;
Sam_ts = 0:Ts:Tend+Ts;
T_length = length(Sam_t);

n = size(G, 1); % state:n
[p, m] = size(D); % output: p, input: m
x_initial = [-351.2222;11.2992;8.3836;8.0217;-4.4119;7.6555;37.6572;-23.7958];
%%  Reference
kk = 0;
for t = 0 : Ts : Tend
    kk = kk + 1;
    rd(:,kk) = 1500;
end


%%  Control 
Qc = 1e4*eye(p);
Rc = 1e0*eye(m);

CTQC = C_ok'*Qc*C_ok;
[Kd1, P, e] = dlqr(G_ok, H_ok, CTQC, Rc);
E_bar = eye(n)-G_ok'+G_ok'*P*inv(eye(n)+H_ok*inv(Rc)*H_ok'*P)*H_ok*inv(Rc)*H_ok';
R_bar = H_ok'*P*H_ok+Rc;
Kd = inv(R_bar)*H_ok'*P*G_ok;
Ed = inv(R_bar)*H_ok'*inv(E_bar)*C_ok'*Qc;
Cu = inv(R_bar)*(Rc-H_ok'*inv(E_bar)*G_ok'*P*inv(eye(n)+H_ok*inv(Rc)*H_ok'*P)*H_ok);

P= 1000;
Kp = P*Ed;
C_tilde = inv(eye(p)+C_ok*H_ok*Kp)*C_ok*(G_ok-H_ok*Kd);
D_tilde = inv(eye(p)+C_ok*H_ok*Kp)*C_ok*H_ok*(Ed + Kp);

inv(eye(p) + C_ok*H_ok*Kp)*C_ok*(G_ok-H_ok*Kd)
inv(eye(p) + C_ok*H_ok*Kp)*C_ok*H_ok*(Ed+Kp)
inv(eye(p) + C_ok*H_ok*Kp)*C_ok
inv(eye(p) + C_ok*H_ok*Kp)*C_ok*H_ok*Cu

u_star = [5000;3.4321;1090;24;2.1837;130;4.4834;380];

load xo_410.mat
xo = xo_410;
yo = C_ok*xo;
xd = x_initial;
yd = C*xd;

% kk = 0;
% for t = 0 : Ts : Tend
%     kk = kk + 1;
%     esys(:,kk) = -0.0073 + 13.3067*randn(1,1);
%     eout(:,kk) = 1e-11*2.2241 + 1e-8*4.86*randn(1,1);
% end
load esys_eout
Kss = [-0.0722;0.0896;-0.0220;-0.0220;-0.0234;-0.0257;-0.0360;0.0164];

kk = 0;
%% 
for t = 0 : Ts : Tend-Ts
    kk = kk + 1;
    td(:, kk) = kk;    
        

    ud(:,kk) = -(Kd + Kp*C_tilde)*xo(:,kk) + (Ed + Kp*(eye(p)-D_tilde))*rd(:,kk) + Cu*u_star;

   
    xo(:,kk+1) = G_ok*xo(:,kk) + H_ok*ud(:,kk) + F_ok*(yo(:,kk) - yd(:,kk));
    yo(:,kk+1) = C_ok*xo(:,kk+1);

    xd(:,kk+1) = G*xd(:,kk) + H*ud(:,kk) + Kss*esys(:,kk);
        yd(:,kk+1) = inv(eye(p) + C_ok*H_ok*Kp)*C_ok*(G_ok-H_ok*Kd)*xd(:,kk) +...
                     inv(eye(p) + C_ok*H_ok*Kp)*C_ok*H_ok*(Ed+Kp)*rd(:,kk+1) +...
                     inv(eye(p) + C_ok*H_ok*Kp)*C_ok*Kss*esys(:,kk)+...
                     inv(eye(p) + C_ok*H_ok*Kp)*eout(:,kk+1)+...
                     inv(eye(p) + C_ok*H_ok*Kp)*C_ok*H_ok*Cu*u_star;


   
end

%% Plot
figure(1)
plot(td,yd(:,1:601),'.-')
legend('y_{d}')
xlim([1 600])


figure(2)
subplot(2,1,1)
plot(td,ud(1,:))
legend('u_{d1}')
line([td(1), td(end)], [5700,5700],'LineStyle','--','LineWidth',1.5,'Color',[0,0,0]);xlim([td(1),td(end)]);
line([td(1), td(end)], [4300,4300],'LineStyle','--','LineWidth',1.5,'Color',[0,0,0]);xlim([td(1),td(end)]);
axis([1 601 4000 6000])
subplot(2,1,2)
plot(td,ud(2,:))
legend('u_{d2}')
line([td(1), td(end)], [3.6,3.6],'LineStyle','--','LineWidth',1.5,'Color',[0,0,0]);xlim([td(1),td(end)]);
line([td(1), td(end)], [1.9,1.9],'LineStyle','--','LineWidth',1.5,'Color',[0,0,0]);xlim([td(1),td(end)]);
axis([1 601 1.5 4])

figure(3)
subplot(2,1,1)
plot(td,ud(3,:))
legend('u_{d3}')
line([td(1), td(end)], [1180,1180],'LineStyle','--','LineWidth',1.5,'Color',[0,0,0]);xlim([td(1),td(end)]);
line([td(1), td(end)], [1000,1000],'LineStyle','--','LineWidth',1.5,'Color',[0,0,0]);xlim([td(1),td(end)]);
axis([1 601 900 1300])
subplot(2,1,2)
plot(td,ud(4,:))
legend('u_{d4}')
line([td(1), td(end)], [33,33],'LineStyle','--','LineWidth',1.5,'Color',[0,0,0]);xlim([td(1),td(end)]);
line([td(1), td(end)], [15,15],'LineStyle','--','LineWidth',1.5,'Color',[0,0,0]);xlim([td(1),td(end)]);
axis([1 601 10 40])

figure(4)
subplot(2,1,1)
plot(td,ud(5,:))
legend('u_{d5}')
line([td(1), td(end)], [3.4,3.4],'LineStyle','--','LineWidth',1.5,'Color',[0,0,0]);xlim([td(1),td(end)]);
line([td(1), td(end)], [0,0],'LineStyle','--','LineWidth',1.5,'Color',[0,0,0]);xlim([td(1),td(end)]);
axis([1 601 -1 5])
subplot(2,1,2)
plot(td,ud(6,:))
legend('u_{d6}')
line([td(1), td(end)], [180,180],'LineStyle','--','LineWidth',1.5,'Color',[0,0,0]);xlim([td(1),td(end)]);
line([td(1), td(end)], [80,80],'LineStyle','--','LineWidth',1.5,'Color',[0,0,0]);xlim([td(1),td(end)]);
axis([1 601 50 250])

figure(5)
subplot(2,1,1)
plot(td,ud(7,:))
legend('u_{d7}')
line([td(1), td(end)], [5,5],'LineStyle','--','LineWidth',1.5,'Color',[0,0,0]);xlim([td(1),td(end)]);
line([td(1), td(end)], [3.6,3.6],'LineStyle','--','LineWidth',1.5,'Color',[0,0,0]);xlim([td(1),td(end)]);
axis([1 601 3 6])
subplot(2,1,2)
plot(td,ud(8,:))
legend('u_{d8}')
line([td(1), td(end)], [430,430],'LineStyle','--','LineWidth',1.5,'Color',[0,0,0]);xlim([td(1),td(end)]);
line([td(1), td(end)], [330,330],'LineStyle','--','LineWidth',1.5,'Color',[0,0,0]);xlim([td(1),td(end)]);
axis([1 601 250 500])

figure(6)
plot(td(:,2:601),yd(:,2:601),'.-')
legend('y_{d}')
xlim([1 600])

