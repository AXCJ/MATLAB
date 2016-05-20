clc; clear; close all;
Ts = 0.1;
Tend = 10;
Sam_t = 0:Ts:Tend-Ts;
T_length = length(Sam_t);
%% Input and Output Data
% u = column vector,u = [ u1 u2...up...ul(uL)]
% y = column vector,y = [ y1 y2...yp...yl(yL)]
load OKIDData


%ªì©l­È
x_initial = [-351.2222;11.2992;8.3836;8.0217;-4.4119;7.6555;37.6572;-23.7958];
%200µ§«áª¬ºA
%x_initial = [-179.849787951128;-123.034404474497;-26.4516361487866;10.1089732218207;23.0374552349171;-8.44590222403733;26.2194539226039;46.1741280143532];
u_star = [5000;3.4321;1090;24;2.1837;130;4.4834;380];
Kss = [-0.0722;0.0896;-0.0220;-0.0220;-0.0234;-0.0257;-0.0360;0.0164];

i = 0;
for t = 1:1:100
    i = i + 1;
    u_randn = [    5000*0.03*randn(1,1);
                  3.4321*0.03*randn(1,1); 
                   1090*0.03*randn(1,1);
                   24*0.03*randn(1,1);
                  2.1837*0.03*randn(1,1);
                   130*0.03*randn(1,1);
                  4.4834*0.03*randn(1,1);
                   380*0.03*randn(1,1)];
    esys(:,i) = -0.0073 + 13.3067*randn(1,1);

    u(:,i) = u_star + u_randn;
end
i = 0;
xd(:,1) = x_initial;
for t = 1:1:100
    i = i + 1;
    xd(:,i+1) = G*xd(:,i) + H*u(:,i) + Kss*esys(:,i);
    y(:,i) = C*xd(:,i);
end

%% Initial parameters
%m and "l" => small L
%* In textbook u = 0~l-1, in there u = 1~l
%* In textbook u(p)     , in there u(p+1)
%* to normalize formula, let (textbook p) = (index_p) => p+1
Alpha = 20;
Beta = 20;
n = 8;
[m, l] = size(y);
r = size(u,1);
v = [u; y];
p = 8;
index_p = p+1;


V_bar = u(:, index_p:l);
V_bar = [];%?????????????
for i = 1 : p
    V_bar = [V_bar; v(:, index_p-i:(l-i))];
end

%% Pseudoinverse
% Y_bar formal: Y_bar = [D C(B_bar) C(A_bar)B_bar ... C(A^(p-1)_bar)B_bar]
y_bar = y(:,index_p:l);
%Y_bar = y_bar*V_bar'*inv(V_bar*V_bar');
Y_bar = y_bar*pinv(V_bar);

%% System Markov parameters
%Yk_bar belong m x [r+(r+m)p]
for i = 1:p
    Yk_bar(:, :, i) = Y_bar(:, ((i-1)*(r+m)+1):(i*(r+m))); 
end
for i = 1:p
    Yk1_bar(:, :, i) = Yk_bar(:, 1:r, i  );
    Yk2_bar(:, :, i) = -Yk_bar(:, ((r+1):(r+m)), i);
end

Pk(:, :, 1) = [ (Yk1_bar(:, :, 1)-Yk2_bar(:, :, 1)*D),  Yk2_bar(:, :, 1)];
Yk(:, :, 1) = Pk(:, 1:r, 1);
Yok(:, :, 1) = Pk(:, (r+1):(r+m), 1);

for i = 2:Alpha + Beta
    if (i <= p)
        sumTempPk = zeros(size(Pk(:, :, 1) ));
        for j = 1:(i-1)
            sumTempPk = sumTempPk + Yk2_bar(:, :, j)*Pk(:, :, i-j);
        end
        Pk(:, :, i) = [ (Yk1_bar(:, :, i)-Yk2_bar(:, :, i)*D),  Yk2_bar(:, :, i)] - sumTempPk;
    else
        sumTempPk = zeros(size(Pk(:, :, 1) ));
        for j = 1:p
            sumTempPk = sumTempPk + Yk2_bar(:, :, j)*Pk(:, :, i-j);
        end
        Pk(:, :, i) = -sumTempPk;
    end
    
    Yk(:, :, i) = Pk(:, 1:r, i);
    Yok(:, :, i) = Pk(:, (r+1):(r+m), i);
end

%% Hankel matrix
H =  [];
for i = 1:Alpha
    Htemp = [];
    for j = 1:(Beta+1)
        Htemp = [Htemp Pk(:, :, i+j-1)];
    end
    H = [H; Htemp];
end

H0 = H;
H0(:, (end-(r+m)+1):end) = [];
H1 = H;
H1(:, 1:(r+m) ) = [];

[A_ok, Btemp, C_ok, er, U, Sigma, V] = era_WT( H1, H0, (r+m), m, n);

B_ok = Btemp(:, 1:r); % B : ID systeim => B
G_ok = Btemp(:, (r+1):(r+m) ); % G : ID systeim => G


%% Test
x_o(:,1) = pinv(C_ok)*(C_ok*x_initial);
for i = 1 : T_length
    y_o(:,i) = C_ok*x_o(:,i);
    x_o(:,i+1) = A_ok*x_o(:,i)  + B_ok*u(:,i)-G_ok*( y(:,i) -  y_o(:,i));
end

i = 0;
for t = 0 : 0.1 : Tend-0.1
    i = i + 1;
    td(:,i) = i;
end

%plot(td,y,'.-');
plot(td,y,td,y_o,'.-');
legend('y','y\_ok');