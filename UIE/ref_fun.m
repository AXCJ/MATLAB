% t : time
% r : reference input
function r = ref_fun(t, mode)
if(mode == 1)
        r(:,1) = 3*sin(t);
  %      r(:,2) = cos(t);
else
    if (t >= 0 && t < 1)
        r(:,1) = cos(2*pi*t);
        r(:,2) = 1.2*(t^2)*(1-t);
    elseif (t >= 1 && t < 2)
        r(:,1) = 0.5*(t^2)*(1-t);
        r(:,2) = cos(2*pi*t);
    else
        r(:,1) = 0.5*cos(4*pi*t+1);
        r(:,2) = 0.2*sin(4*pi*t)-0.5;
    end
end
end