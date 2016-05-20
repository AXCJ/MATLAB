clc; clear; close all;
%%
load A_ok
load B_ok
load C_ok


Ts = 0.1;
Tend = 10;
Sam_t = 0:Ts:Tend-Ts;
T_length = size(Sam_t,2);

u_star = [5000;3.4321;1090;24;2.1837;130;4.4834;380];
Kss = [-0.0722;0.0896;-0.0220;-0.0220;-0.0234;-0.0257;-0.0360;0.0164];

Tend = 10;
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
%% 系統創建

%初始狀態
%x_initial = [-351.2222;11.2992;8.3836;8.0217;-4.4119;7.6555;37.6572;-23.7958];
%200筆後狀態
x_initial = [-179.849787951128;-123.034404474497;-26.4516361487866;10.1089732218207;23.0374552349171;-8.44590222403733;26.2194539226039;46.1741280143532];
x_o(:,1) = pinv(C_ok)*(C_ok*x_initial);
xd(:,1) = x_initial;
for i = 1 : T_length

    y_o(:,i) = C_ok*x_o(:,i);
    x_o(:,i+1) = A_ok*x_o(:,i)  + B_ok*u(:,i);

    xd(:,i+1) = A_ok*xd(:,i) + B_ok*u(:,i) + Kss*esys(:,i);
    y(:,i) = C_ok*xd(:,i);


end



%% 繪圖
i = 0;
for t = 0 : 0.1 : Tend-0.1
    i = i + 1;
    td(:,i) = i;
end
%plot(td,y,'.-');
plot(td,y,'b',td,y_o,'r');
legend('y','y\_ok');