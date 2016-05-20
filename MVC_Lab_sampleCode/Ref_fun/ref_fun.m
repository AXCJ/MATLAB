function [u1,u2] = ref_fun( t )

  if t < 5
     u1  = cos(pi/2*t);
     u2  = sin(pi/2*t);
  else
     u1 = -1;
     u2 =  1;
  end

end

