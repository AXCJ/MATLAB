clear; close all; clc;

Tf = 1e-2;
Tend = 10;
[A,B,C,D]=zp2ss([],[-3 -5],1);
Tseq = 0:Tf:Tend;


% n : the number of states of plant
% r : the number of inputs of plant
% p : the number of outputs of plant
[n,r] = size(B);
p = size(C,1);


% ~~~ Controller design
Qc = 1e6*eye(p);
Rc = 1e0*eye(r);
Kc = lqr(A, B, C'*Qc*C, Rc); % Tracker --> C'*Q*C
Ec = -inv(Rc)*B'*inv(A - B*Kc)'*C'*Qc ;
% Ec = -pinv(C*inv(A - B*Kc)*B)
% ~~~ Controller design ~~~ End


r = sin(Tseq);
simin = [Tseq; r]';
[tt , xx , y1, u] = sim('control_input', Tseq);
figure(1)
plot(Tseq, r, Tseq, y1)
figure(2)
plot(Tseq, u)


%%
Kc = Kc*0;
Ec = Ec / Ec;
simin = [Tseq; u']';
[tt , xx , y2, uu] = sim('control_input', Tseq);
figure(3)
plot(Tseq, r, Tseq, y2)
figure(4)
plot(Tseq, y1, Tseq, y2)


% 所以只要求出desired trajetory對應的control input，即可得到理想的輸出軌跡





