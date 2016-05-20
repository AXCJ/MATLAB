clc; clear; close all;

%% original system

A = [-0.0558 -0.9968  0.0802 0.0415;
      0.5980 -0.1150 -0.0318      0;
     -3.0500  0.3880 -0.4650      0;
           0  0.0805       1      0];
       
B=[ 0.00729  0.0583;
   -0.4750  -2.0100;
    0.1530   0.0241;
         0        0];
     
C=[0 1 0 0;
   0 0 1 0];

D=[0 0;
   0 0];

syss=ss(A,B,C,D);
%% pole placement

[T,D1]=eig(A);
invT=inv(T);

% desire pole
POLE=[-0.5-1i -0.5+1i -0.5 -2];

L_polePlace=place(A,B,POLE) % L=acker(A,B,POLE) for single input
% Complete Modal Control
% [L,lambdaCL]=Roppenecker(A,B,D1,T,POLE)

[p,z]=pzmap(A-B*L,B,C,D)

L=[-20.1334 3.6117 13.7112 6.7387;
    5.0231 -1.2026 -3.3337 -1.6468];
sys=ss(A-B*L,B,C,D);
x0=[0.15 0 0 0];
x1=[0 0 0 0.5];

%% figure1 (a)

figure
[y,t,x]=initial(sys,x0);

plot(t,x(:,1),'r','LineWidth',2)
xlabel('t(s)','fontsize',16,'fontweight','b');
title('(a)','fontsize',16,'fontweight','b');
text(t(20,1),x(20,1),'\color{red} \leftarrow \beta',...
     'HorizontalAlignment','left',...
     'FontSize',18)
set(gca,'fontsize',14);
hold on
plot(t,x(:,4),'LineWidth',2);
text(t(51,1),x(51,4),'\color{blue} \downarrow \phi',...
     'HorizontalAlignment','left',...
     'FontSize',18)
hold off
grid

axis([0 10 -0.05 0.2])


%% figure2 (b)

[y1,t,x1]=initial(sys,x1);
figure
% beta
plot(t,x1(:,1),'r','LineWidth',2);
xlabel('t(s)','fontsize',16,'fontweight','b');
title('(b)','fontsize',16,'fontweight','b');
text(t(29),x1(29,1),'\color{red} \downarrow \beta',...
     'HorizontalAlignment','left',...
     'FontSize',18)
set(gca,'fontsize',14);
hold on
% phi
plot(t,x1(:,4),'LineWidth',2);
text(t(54),x1(54,4),'\color{blue} \leftarrow \phi',...
     'HorizontalAlignment','left',...
     'FontSize',18)
hold off
grid

axis([0 10 -0.2 0.7])

% close all