function [ y_c,y_o,tf ] = solver(xi, wc, sys_idx, Ts, Tend, ref, Tf)


num = [1];
den = [1 2*xi(sys_idx)*wc(sys_idx) wc(sys_idx)^2];
[A,B,C,D] = tf2ss(num,den);
sys = ss(A,B,C,D);

x_ini = [ 0  ; 0 ];
% ~~~ System setting ~~~ End  

[n,r] = size(B);
p = size(C,1); % the size of the first dimension of C 
% n : the number of states  of plant
% r : the number of inputs  of plant
% p : the number of outputs of plant

% ~~~ Observer setting
[Ao,Bo,Co,Do] = tf2ss(num,den);

xo_ini = [ 0.5  ; 0.5 ];
% ~~~ Observer setting ~~~ End  

[no,ro] = size(Bo);
po = size(Co,1);
% no : the number of states  of observer
% ro : the number of inputs  of observer
% po : the number of outputs of observer

% ~~~ Observer design
Ro = 1e0*eye(ro);
Qo = 5e3*eye(no); % Regulator --> the size of Qo is no by no

[Ko,Po] = lqr(A', C', Qo, Ro); % Regulator --> Qo
Ko = Ko'; % (A, B ,Kc) <---> (A', C', Ko') for continuous-time system
% Design Ko such that ( A - Ko*C ) is asymptotically stable
% ~~~ Observer design ~~~ End

% --- Controller design
R = 1e0*eye(ro);
Q = 5e6*eye(po); % Tracker --> the size of Q is po by po

[Kc,Pc] = lqr(Ao, Bo, Co'*Q*Co, Ro); % Tracker --> C'*Q*C
% Ec = -inv(R)*Bo'*inv(Ao - Bo*Kc)'*Co'*Q;
Ec = -inv(Co*inv(Ao-Bo*Kc)*Bo);
% --- Controller design --- end
sys_ob = ss(A-B*Kc-Ko*C,Ko,-Kc,[]);
[p,z] = pzmap(sys_ob);
x_c(:,1) = x_ini;
x_o(:,1) = xo_ini;
ii = 0;
kk = 0;
% ~~~ Perform simulation
for t_ = 0 : Tf : Tend
  ii = ii + 1;
  tf(:,ii) = t_; % continuous timespan
  if ~mod(ii-1,Ts/Tf)
        kk = kk + 1;
  end
  ref_c(:,ii) =ref(:,kk);
  u_c(:, ii) = Ec*ref_c(:,ii) - Kc*x_o(:,ii); % Tracker
  options = odeset('RelTol',1e-7);
  tspan = [ t_ , t_+Tf ];

  [Tc,Xc] = ode45(@(tt, xx) plant1(tt, xx, u_c(:, ii), A, B), tspan, x_c(:,ii), options);
  y_c(:, ii) = C*x_c(:, ii);

  [To,Xo] = ode45(@(tt, xx) observer1(tt, xx, u_c(:, ii), Ao, Bo, Co, Ko, y_c(:, ii)), tspan, x_o(:,ii), options);
  y_o(:, ii) = Co*x_o(:, ii);

  if(t_ >= Tend)
      break;
  end

  x_c(:,ii+1) = Xc(end,:)';
  x_o(:,ii+1) = Xo(end,:)';

end


end

