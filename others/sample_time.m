clc,clear;close all;
A=[0 1;-2 -3];
B=[0;1];
C=[1 1];
D=[];
[G,H]=c2d(A,B,0.01);
sys1=ss(A,B,C,D)
sys2=ss(G,H,C,D,0.01)
figure
step(sys1)
hold on
step(sys2)