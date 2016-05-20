clc;clear all; close all;
s = tf('s');
xi = [0.2 0.1 0.1 0.3 0.3];
wc = [0.8 0.5 1.2 0.5 1.2];

sys_idx = 3;
G = 1/(s^2+2*xi(sys_idx)*wc(sys_idx)*s+wc(sys_idx)^2);
% [A,B,C,D] = tf2ss(1,[1 2*xi(sys_idx)*wc(sys_idx) wc(sys_idx)^2]);
% k = acker(A,B,[pole-3])
T=feedback(G,1);
p=1/bode(T,0);
T = p*T;
[y,t] = step(T);
info = stepinfo(T);
fprintf('Kp:%f\tSettlingTime: %f\tsteady state error:%f\n',p,info.SettlingTime,y(end)-1)

