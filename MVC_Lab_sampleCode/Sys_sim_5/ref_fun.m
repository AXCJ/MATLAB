% <<< Generate a reference input >>> %
function r = ref_fun( t )
  
  if t > 5
    r = -0.3;
  elseif t > 3
    r = 2;
  else
    r = 1; 
  end

end
% t : time
% r : reference input
