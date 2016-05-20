% hyyly
% Reference
% http://www.mathworks.cn/matlabcentral/fileexchange/13634-axescoord2figurecoord
% http://www.ilovematlab.cn/thread-245377-1-1.html
% http://www.ilovematlab.cn/thread-246998-1-1.html

x = linspace(-2*pi,8*pi,200);
y = sin(x);
f=figure(1);
h = plot(0,0,'linewidth',1);
xMin=-2*pi;  %X-axis minimum
xMax=2*pi;   %X-axis maximum
xlim([xMin,xMax]);  %X-axis initial range
ylim([-1.5 1.5]);
title('Power by hyyly');
set(f,'Units','normalized');
% % % % % % % % % % % % 
m=2;
xa = [x(m-1) x(m)];
ya = [y(m-1) y(m)];
[xaf,yaf] = axescoord2figurecoord(xa,ya);  %axes coordinate to figure coordinate
ar=annotation(f,'arrow',xaf,yaf,'Color',[1 0 0],'LineStyle','none','HeadStyle','vback1','HeadLength',12,'HeadWidth',10);
    

for m = 2:numel(x)
    if(x(m)+pi>xMax)
        xdif=x(m)+pi-xMax;
        xlim([xMin+xdif,xMax+xdif]);
    end
    set(h,'xData',x(1:m),'yData',y(1:m));
    xa = [x(m-1) x(m)];
    ya = [y(m-1) y(m)];
    [xaf,yaf] = axescoord2figurecoord(xa,ya);
    %annotation(f,'arrow',xaf,yaf,'Color',[1 0 0],'LineStyle','none');
    set(ar,'X',xaf,'Y',yaf);
    drawnow;
    pause(0.1);
end

