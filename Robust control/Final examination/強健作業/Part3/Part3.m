clc; clear all; close all;
load('data.mat');
t = track1.time;
dt = 1;

trackerLocation = [50,0; -20,0];
laserLocation = [0 0];

figure;
plot(t,track1.angle);
figure;
plot(t,track2.angle);
figure;
plot(t,track1.angle,t,track2.angle);

a1 = smooth(track1.angle,15)';
a2 = smooth(track2.angle,15)';
% a1 = track1.angle;
% a2 = track2.angle;

figure;
plot(t,a1);
figure;
plot(t,a2);
figure;
plot(t,a1,t,a2);

for i = 1:size(t,2)
    s1 = a1(1,i);
    s2 = a2(1,i);
    xt(:,i) = (50*tan(s1)+20*tan(s2))/(tan(s1)-tan(s2));
    yt(:,i) = tan(s1)*(xt(:,i)-50);
end

figure;
plot(t,xt);
figure;
plot(t,yt);

figure;
plot(xt,yt);

X = [xt;yt];
%X = [x, y, dx, dy]';
A = [1 0 dt 0; 0 1 0 dt; 0 0 1 0; 0 0 0 1]; %state update matrice
B = [(dt^2/2); (dt^2/2); dt; dt];
C = [1 0 0 0; 0 1 0 0];

detectedLocations = X;

figure;
hold on;
ylabel('Location');
ylim([0,800]);
xlabel('Time');
xlim([-1200,1200]);
kalman = [];

for idx = 1: length(detectedLocations)
    location = detectedLocations(:,idx);
    if isempty(kalman)
        if ~isnan(location)
            
            stateModel = A;
            %ControlModel = B
            measurementModel = C;
            kalman = vision.KalmanFilter(stateModel,measurementModel,'ProcessNoise',1e-2,'MeasurementNoise',4);
            kalman.State = [location; 0; 0];
        end
    else
        trackedLocation = predict(kalman);
        xk(idx) = trackedLocation(1);
        yk(idx) = trackedLocation(2);
        if ~isnan(location)
            plot(location(1), location(2),'k--.');
            d = distance(kalman,location);
            title(sprintf('Distance:%f', d));
            trackedLocation = correct(kalman,location);
        else
            title('Missing detection');
        end
%         pause(1e-14);
%         plot(xk(1:idx),yk(1:idx),'r');
    end
end
legend('Detected locations','Predicted/corrected locations');

figure;
plot(t,xk);

figure;
plot(t,yk);

[angle1hat, angle2hat] = returnSita(xk,yk,trackerLocation(1,1),trackerLocation(2,1));

figure('Name', 'angle1' ) % Middle_Top of screen
plot(angle1hat); hold on;plot(track1.angle); hold off;legend({'$\hat{\theta_1}$','$\theta_1$'},'Interpreter','latex')
figure('Name', 'angle2' ) % Middle_Bottom of screen
plot(angle2hat); hold on;plot(track2.angle); hold off;legend({'$\hat{\theta_2}$','$\theta_2$'},'Interpreter','latex')

for i = 1:size(xk,2)
    if(xk(i) < 0)
        laserSai(i) = pi + atan(yk(i)/xk(i));
    else
        laserSai(i) = atan(yk(i)/xk(i));
    end
end
laserSai(1) = atan(yt(1)/xt(1));
figure;
plot(laserSai);

%% Control 
r = laserSai(2:end);
t2 = t;
t2(end) = [];

x = 0;
x0_o = [0.5 0.5]';
xi = 0.2;
omegaC = 0.8;

sys = tf([1],[1 2*xi*omegaC omegaC^2]);
[A,B,C,D] = tf2ss([1],[1 2*xi*omegaC omegaC^2]);

nr = size(A,1);     % No. of states
mr = size(B,2);     % No. of inputs
pr = size(C,1);     % No. of outputs

Qc = 1e2*eye(pr);
Rc = eye(mr);

Rc_bar = Rc + D'*Qc*D;
[F,Pc] = lqr(A, B, C'*Qc*C, Rc_bar);
F = -F;
% F = -place(A,B,[-30+5i -30-5i]);
M = -inv(C*inv(A+B*F)*B);

Ao=A';
Bo=C';
Co = B';
Qo = 1e4*eye(size(Ao,1));
[H,Po] = lqr(Ao, Bo, Qo, Rc_bar);
H = -H';

opts_sim = simset('Solver', 'ode45', 'FixedStep', dt, 'InitialStep', dt, 'MaxStep', dt, 'RelTol', 1e-7);
[T, XX, U, Y] = sim('trackerSys', t2, opts_sim, [ t2', [r'] ]);
T = T';
Y = Y';

figure();
plot(T, Y,t2, r,'--.'); 
grid on;

Error = Y-r;
figure;
plot(Error);
% (x1y2+x2y3+x3y1-x2y1-x3y2-x1y3) /2
x23 = xk(1:end-1);
y23 = yk(1:end-1);
for i = 1:size(t2,2)
    sL = Y(1,i);
    s1 = angle1hat(1,i);
    s2 = angle2hat(1,i);
    x12(i) = 20*tan(s2)/(tan(sL)-tan(s2));
    y12(i) = tan(sL)*x12(i);
    x13(i) = -50*tan(s1)/(tan(sL)-tan(s1));
    y13(i) = tan(sL)*x12(i);
    area(i) = 0.5*abs((x13(i)*y12(i)+x12(i)*y23(i)+x23(i)*y13(i))-(x12(i)*y13(i)+x23(i)*y12(i)+x13(i)*y23(i)));
end

figure;
plot(area);

