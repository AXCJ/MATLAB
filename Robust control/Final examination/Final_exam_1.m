clear;clc;close all;
% scrsz = get(groot,'ScreenSize');
 load('angle.mat')
 theta1 = angle1;
 theta2 = angle2;
% %  method1
%  position=find(theta1<theta2);
%  compensation=max(abs(theta1(position)-theta2(position)));
%  theta1(position)= theta1(position)+compensation+0.001;
 %  method2
 %角度補償計算部分
 position=find(theta1<theta2);
 compensation=max(abs(theta1(position)-theta2(position)))+0.01;
 theta1= theta1+compensation;
 %補償結束
 %method3
%  position=find(theta1<theta2);
%  compensation=abs(theta1(position)-theta2(position));
%  theta1(position)= theta1(position)+compensation+0.01;

 
for ii=0:1:850
i=ii+1;
x(i)=(50*tan(theta1(i))+20*tan(theta2(i)))/(tan(theta1(i))-tan(theta2(i)));
 y(i)=(70*tan(theta1(i))*tan(theta2(i)))/(tan(theta1(i))-tan(theta2(i)));
end
 Z=[x;y];



%%
%   Description : kalmanFiltering
a=[1;2];%加速度值
X=[0 ; 0 ; 0 ; 0]; %初始狀態
P=eye(4); %狀態協方差矩陣
F=[1 0 1 0;0 1 0 1;0 0 1 0;0 0 0 1]; %狀態轉移矩陣
B=[0.5 0;0 0.5;1 0;0 1];%輸入矩陣
Q=1e-4*eye(4); %狀態轉移協方差矩陣
H=[1 0 0 0;0 1 0 0]; %觀測矩陣
R=4*eye(2); %觀測雜訊協方差矩陣

for i = 1:length(Z)
%基於上一狀態預測當前狀態  
 X_ = F*X;%+B*a;
% 更新協方差  Q系統過程的協方差  這兩個公式是對系統的預測
P_ = F*P*F'+Q;
% 計算卡爾曼增益
K = P_*H'/(H*P_*H'+R);
% 得到當前狀態的最優化估算值 +增益乘以殘差
X = X_+K*(Z(:,i)-H*X_);
%更新K狀態的協方差
P = (eye(4)-K*H)*P_;
Y(:,i)=H*X;
end
% plot(Y(1,:),'r')
% figure
% plot(Z(1,:),'b')

% figure
% % figure('Name', 're' ,'Position',[9 scrsz(4)/1.7 scrsz(3)/3 scrsz(4)/3]) % Left_Top of screen
% % plot( 1:length(Y(2,:)) ,Y(2,:),'-',1:length(Z(2,:)), Z(2,:),'-')
% plot( 1:length(Y(2,:)) ,Y(2,:))
% figure
% % figure('Name', 'Op' ,'Position',[scrsz(3)/1.5 scrsz(4)/1.7 scrsz(3)/3 scrsz(4)/3]) % Right_Top of screen
% plot(1:length(Z(2,:)), Z(2,:),'-')


x = Y(1,:);
y = Y(2,:);


% for ii = 0:length(x)-1
%     i=ii+1;
%     % RIGHT
%   eq1 = [int2str(x(i)) '*(tan(K1)-tan(K2))-50*tan(K1)-20*tan(K2)=0'];
%   eq2 =[ int2str(y(i)) '*(tan(K1)-tan(K2))-70*tan(K1)*tan(K2)=0'];
%   [K1(:,i),K2(:,i)]=solve(eq1,eq2);
%   i
% end
%     theta1hat=eval(sum(K1));
%     theta2hat=eval(sum(K2));
for ii = 0:length(Z)-1
    i=ii+1;
    % RIGHT
    if x(i)>=50
        theta1hat(i) = atan(y(i)/(x(i)-50));
        theta2hat(i) = atan(y(i)/(x(i)+20)); 
    elseif x(i)<=-20
    % left
        theta1hat(i) =  pi-atan(y(i)/(50-x(i)));
        theta2hat(i) =  pi-atan(y(i)/(-20-x(i)));
    else
    % middle
        theta1hat(i) = pi- atan(y(i)/(50-x(i)));
        theta2hat(i) = atan(y(i)/(x(i)+20));
    end
end

for ii = 0:length(Z)-1
    i=ii+1;
    if x(i)>=0
        thetak(i)=atan(y(i)/(x(i)));
    else
         thetak(i)=pi-atan(y(i)/(x(i)));
    end
end


figure
plot(theta1,'b'); 
hold on;
plot(theta1hat,'r');
hold off
figure
plot(theta2,'b');
hold on; 
plot(theta2hat,'r'); 
hold off
save thetak.mat 'thetak'
save xy.mat 'x' 'y'


