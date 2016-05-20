clc; clear all; close all;
A = [0 1 ; -2 -1];
B = [0 ; 1];
C = [1 0];
D = 0;
sys = ss(A,B,C,D);
xini = [1;1];
[y,t,x] = initial(sys, xini);
[V , D] = eig(A)
plot(x(:,1) , x(:,2))
hold on

xini = [1;-1];
[y,t,x] = initial(sys, xini);
plot(x(:,1) , x(:,2))

xini = [-1;1];
[y,t,x] = initial(sys, xini);
plot(x(:,1) , x(:,2))

xini = [-1;-1];
[y,t,x] = initial(sys, xini);
plot(x(:,1) , x(:,2))

line([0 0],V(:,1)')




hold off
figure
plot(t,x(:,1) ,t, x(:,2))

