function [angle1hat, angle2hat] = returnSita(cal_x,cal_y,x1,x2)

for ii = 0:length(cal_x)-1
    i=ii+1;
    if(cal_y(i)>0)
        % right
        if cal_x(i)>x1
            angle1hat(i) = atan(cal_y(i)/(cal_x(i)-x1));
            angle2hat(i) = atan(cal_y(i)/(cal_x(i)-x2));
        elseif x2<cal_x(i) && cal_x(i)<x1
        % middle
            angle1hat(i) = pi - atan(cal_y(i)/(x1-cal_x(i)));
            angle2hat(i) = atan(cal_y(i)/(cal_x(i)-x2));
        elseif cal_x(i)<x2
        % left
            angle1hat(i) = pi - atan(cal_y(i)/(x1-cal_x(i)));
            angle2hat(i) = pi - atan(cal_y(i)/(x2-cal_x(i)));

        end
    else
        % right
        if cal_x(i)>x1
            angle1hat(i) = pi + atan(cal_y(i)/(cal_x(i)-x1));
            angle2hat(i) = pi + atan(cal_y(i)/(cal_x(i)-x2));
         elseif x2<cal_x(i) && cal_x(i)<x1
        % middle
            angle1hat(i) = atan(cal_y(i)/(cal_x(i)-x1));
            angle2hat(i) = atan(cal_y(i)/(cal_x(i)-x2));
        elseif cal_x(i)<x2
        % left
            angle1hat(i) = atan(cal_y(i)/(cal_x(i)-x1));
            angle2hat(i) = atan(cal_y(i)/(-x2+cal_x(i)));
        end
    end
end

