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

a1 = smooth(track1.angle,30)';
a2 = smooth(track2.angle,30)';
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
       plot(location(1), location(2),'k+');
      d = distance(kalman,location);
       title(sprintf('Distance:%f', d));
       trackedLocation = correct(kalman,location);
     else
       title('Missing detection');
     end
     pause(1e-9);
     plot(trackedLocation(1),trackedLocation(2),'ro');
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
        laserSai(i) = pi-atan(yk(i)/xk(i));
    else
        laserSai(i) = atan(yk(i)/xk(i));
    end
end
figure;
plot(laserSai);