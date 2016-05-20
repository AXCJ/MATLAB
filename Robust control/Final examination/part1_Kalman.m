%%
%
%   Description : kalmanFiltering

clc;clear all;close all;

scrsz = get(groot,'ScreenSize');
load('angle.mat')
pos = find(angle2>angle1);
offset = max(angle2(pos)-angle1(pos))+0.001;
angle1 = angle1+offset;

for ii=0:1:850
i=ii+1;
    x(i)=(50*tan(angle1(i))+20*tan(angle2(i)))/(tan(angle1(i))-tan(angle2(i)));
    y(i)=70*tan(angle1(i))*tan(angle2(i))/(tan(angle1(i))-tan(angle2(i)));
%     if(angle1(i) < angle2(i))
%         fprintf('i: %03d    angle1: %5f   angle2: %5f   x: %5f   y: %5f\n',i,angle1(i)*180/pi,angle2(i)*180/pi,x(i+1),y(i+1));
%     end
end



%%
X=[x(1) ; y(1) ; 0 ; 0]; %初始狀態
P=eye(4); %狀態協方差矩陣
F=[1 0 1 0;0 1 0 1;0 0 1 0;0 0 0 1]; %狀態轉移矩陣
B=[0.5 0;0 0.5;1 0;0 1];
Q=1e-3*eye(4); %狀態轉移協方差矩陣
H=[1 0 0 0;0 1 0 0]; %觀測矩陣
R=4*eye(2); %觀測雜訊協方差矩陣
Z=[x;y];
for i = 1:length(x)
%基於上一狀態預測當前狀態  
X_ = F*X;
% 更新協方差  Q系統過程的協方差  這兩個公式是對系統的預測
P_ = F*P*F'+Q;
% 計算卡爾曼增益
K = P_*H'*inv(H*P_*H'+R);
% 得到當前狀態的最優化估算值  增益乘以殘差
X = X_+K*(Z(:,i)-H*X_);
%更新K狀態的協方差
P = (eye(4)-K*H)*P_;
Kal_Output(:,i)=H*X;
end
figure('Name', 'Y' ,'Position',[9 scrsz(4)/1.9 scrsz(3)/3 scrsz(4)/2.55]) % Left_Top of screen
plot( 1:length(Kal_Output(2,:)) ,Kal_Output(2,:),'-',1:length(Z(2,:)), Z(2,:),'-')

figure('Name', 'X' ,'Position',[9 scrsz(4)/1.9/10 scrsz(3)/3 scrsz(4)/2.55]) % Left_Bottom of screen
plot( 1:length(Kal_Output(1,:)) ,Kal_Output(1,:),'-',1:length(Z(1,:)), Z(1,:),'-')

cal_x = Kal_Output(1,:);
cal_y = Kal_Output(2,:);
% cal_x = x;
% cal_y = y;
for ii = 0:length(cal_x)-1
    i=ii+1;
    if(cal_y(i)>0)
        % right
        if cal_x(i)>50
            angle1hat(i) = atan(cal_y(i)/(cal_x(i)-50));
            angle2hat(i) = atan(cal_y(i)/(cal_x(i)+20));
        elseif -20<cal_x(i) && cal_x(i)<50
        % middle
            angle1hat(i) = pi - atan(cal_y(i)/(50-cal_x(i)));
            angle2hat(i) = atan(cal_y(i)/(cal_x(i)+20));
        elseif cal_x(i)<-20
        % left
            angle1hat(i) = pi - atan(cal_y(i)/(50-cal_x(i)));
            angle2hat(i) = pi - atan(cal_y(i)/(-20-cal_x(i)));

        end
    else
        % right
        if cal_x(i)>50
            angle1hat(i) = pi + atan(cal_y(i)/(cal_x(i)-50));
            angle2hat(i) = pi + atan(cal_y(i)/(cal_x(i)+20));
         elseif -20<cal_x(i) && cal_x(i)<50
        % middle
            angle1hat(i) = atan(cal_y(i)/(cal_x(i)-50));
            angle2hat(i) = atan(cal_y(i)/(cal_x(i)+20));
        elseif cal_x(i)<-20
        % left
            angle1hat(i) = atan(cal_y(i)/(cal_x(i)-50));
            angle2hat(i) = atan(cal_y(i)/(20+cal_x(i)));
        end
    end
    if(i<40 && angle2hat(i)<1)
        fprintf('kal_x: %5f   kal_y: %5f\nori_x: %5f   ori_y: %5f i=%d\n',cal_x(i),cal_y(i),x(i),y(i),i);
    end
end
figure('Name', 'angle1' ,'Position',[scrsz(3)/3 scrsz(4)/1.9 scrsz(3)/3 scrsz(4)/2.55]) % Middle_Top of screen
plot(angle1,'.');hold on; plot(angle1hat,'linewidth',2); hold off;legend({'$\hat{\theta_1}$','$\theta_1$'},'Interpreter','latex', 'fontsize',20)
title({'$\hat{\theta_1} vs. \theta_1$'},'Interpreter','latex', 'fontsize',30)
set(gca,'FontSize',20)
figure('Name', 'angle2' ,'Position',[scrsz(3)/3 scrsz(4)/1.9/10 scrsz(3)/3 scrsz(4)/2.55]) % Middle_Bottom of screen
plot(angle2,'.'); hold on; plot(angle2hat,'linewidth',2); hold off;legend({'$\hat{\theta_2}$','$\theta_2$'},'Interpreter','latex', 'fontsize',20)
title({'$\hat{\theta_2} vs. \theta_2$'},'Interpreter','latex', 'fontsize',30)
set(gca,'FontSize',20)
for i = 1:length(Kal_Output(1,:))
    if(Kal_Output(1,i)>0)
        laserAngle(i) = atan(Kal_Output(2,i)/Kal_Output(1,i));
    else
        laserAngle(i) = pi - atan(-Kal_Output(2,i)/Kal_Output(1,i));
    end
end

figure('Name', 'Laser' ,'Position',[scrsz(3)/1.5 scrsz(4)/1.9 scrsz(3)/3 scrsz(4)/2.55]) % Right_Top of screen
plot(laserAngle)
title('\bf\fontsize{30}\fontname{Cambria} Laser')
title({'$\theta_{Laser}$'},'Interpreter','latex', 'fontsize',30)
set(gca,'FontSize',20)


figure('Name', 'x-y' ,'Position',[scrsz(3)/1.5 scrsz(4)/1.9/10 scrsz(3)/3 scrsz(4)/2.55]) % Right_Bottom of screen

% h = plot(0,0,'XDataSource','x','YDataSource','y');
% ylim([0 1000]);
% xlim([-1000 1000]);
% for k=1:length(Kal_Output(2,:))
%     x=Kal_Output(1,1:k);
%     y=Kal_Output(2,1:k);
%     refreshdata(h)
%     pause(1e-9)
%     drawnow
% end

hold on
title('\bf\fontsize{30}\fontname{Cambria} X vs. Y')
ylim([0 1000]);
xlim([-1000 1000]);
for i=1:length(Kal_Output(2,:))
    h = plot(0,0,'linewidth',1);
    plot(Kal_Output(1,1:i),Kal_Output(2,1:i),'r','linewidth',1)
%     ylabel(['i = ' num2str(i)],'rotation',0,'position',get(get(gca,'YLabel'),'Position') - [0.090 0.01 0])
    pause(1e-9)
end
% set(gca,'FontSize',20)

x = Kal_Output(1,:);
y = Kal_Output(2,:);
save laserAngle 'laserAngle' 'x' 'y' 'angle1hat' 'angle2hat'

figure('Name', 'Laser' ,'Position',[scrsz(3)/1.5 scrsz(4)/1.9 scrsz(3)/3 scrsz(4)/2.55]) % Right_Top of screen
hold on;plot(angle1hat,':','linewidth',2);plot(angle2hat,'--','linewidth',2);
% title('\bf\fontsize{30}\fontname{Cambria} Laser vs. \')
set(gca,'FontSize',20)
legend({'$\hat{\theta_1}$','$\hat{\theta_2}$'},'Interpreter','latex', 'fontsize',20)
