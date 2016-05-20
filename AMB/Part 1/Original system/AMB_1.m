clear all;

%---Data from Table 8.1---%
ksA = -1e7, ksB = ksA;
kiA = 250, kiB = ksA;
m = 100;
Ix = 8.3333, Iy = Ix;
Iz = 0.75;
PA = -2*ksA/kiA, PB = PA;
DA = sqrt(-m*ksA/2)/kiA, DB = DA;
b = 0.4, a = -b;
d = 0.45, c = -d;
omega = 500*2*pi;
%-------------------------%
M = diag([Ix m Ix m]);
G = [zeros(2) [Iz*omega 0; 0 0]; [-Iz*omega 0; 0 0] zeros(2)];
B = [[a b;1 1] zeros(2); zeros(2) [a b; 1 1]];
C = [[c 1; d 1] zeros(2); zeros(2) [c 1; d 1]];
Ki = diag([kiA kiB kiA kiB]);

tempKss = [ksA*a^2 + ksB*b^2 ksA*a + ksB*b;
           ksA*a^2 + ksB*b   ksA + ksB];
Kss = [tempKss zeros(2,2);
       zeros(2,2) tempKss];
P = diag([PA PB PA PB]);
D = diag([DA DB DA DB]);

%from PD feedback control
Kc = B*Ki*P*C;  %stiffness matrices
Dc = B*Ki*D*C;  %damping matrices


A_final = [zeros(4,4) eye(4); -B'*inv(M)*(Kss)*inv(B') -B'*inv(M)*(G)*inv(B')];
B_final = [zeros(4); B'*inv(M)*B*Ki]
C_final = [C*inv(B') zeros(4)];

% i = [ixA ixB iyA iyB]'
%x = [beta x -alpha y]
%y = [XseA XseB YseA YseB]
sys = ss(A_final , B_final, C_final, [0]);
t = 0:0.0002:5;
[y,t] = step(sys);
figure(1);
step(sys);
figure(2);
yy = sum(y,3);
subplot(221);
plot(t,yy(:,1));
subplot(222);
plot(t,yy(:,2));
subplot(223);
plot(t,yy(:,3));
subplot(224);
plot(t,yy(:,4));
