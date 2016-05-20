function [polrac] = Get_Eigenvalue(nval,Ts,maxident)
%
% Acquisition of eigenvalues ("poles"), real or complex conjugated, continuous
%
%  [polrac] = Get_Eigenvalue(nval,Ts,maxident)  prompts the user to enter a number nval
%  of continuous-time eigenvalues, which are then converted to discrete-time
%  eigenvalues if a strictly positive sampling time is transmitted by the input
%  argument "Ts".
%  The entered eigenvalues or their discrete equivalents are returned by the vector
%  "polrac".
%  The input parameter "maxident" stops the acquisition process if a number maxident
%  of identical values entered consecutively is exceeded, and obliges choosing another
%  value
%
%  Author: E. Ostertag, 10 February 2010
%  Last update: 19 February 2011
%

ni = nargin;
if ni < 3,
  maxident = nval;
end
if ni > 1 && Ts > 0,
  disp('*** Acquisition of desired continuous-time eigenvalues (program will convert them to discrete values): ***');
else
  disp('*** Acquisition of desired continuous-time eigenvalues: ***');
end
ieig = 1;
if nval == 1
  lambdades(ieig) = input(' Desired eigenvalue (unique, thus necessarily real) = ');
else
  while ieig <= nval,
    if ieig <= nval-1,
      lambdades(ieig) = input([' Desired eigenvalue #' num2str(ieig) ' (real or complex) = ']);
    else
      lambdades(ieig) = input([' Desired eigenvalue #' num2str(ieig) ' (only real) = ']);
      if ~isreal(lambdades(ieig)),
        lambdades(ieig) = real(lambdades(ieig));
        disp(['      Converted to real value: ' num2str(lambdades(ieig)) ]);
      end
    end
    if ieig == 1,
      nident = 1;
    else
      if lambdades(ieig) == lambdades(ieig-1),
        nident = nident+1;
      else
        nident = 1;
      end
    end
    if nident > maxident,
      while lambdades(ieig) == lambdades(ieig-1),
        disp('Number of identical values > allowed number: enter again --->'); beep;
        lambdades(ieig) = input([' Desired eigenvalue #' num2str(ieig) ' = ']);
      end
      nident = 1;
    end
    if isreal(lambdades(ieig)),
      ieig = ieig+1;
    else
      lambdades(ieig+1) = conj(lambdades(ieig));
      disp(['Desired eigenvalue #' num2str(ieig+1) ' (complex conjugated) = ' ....
             num2str(lambdades(ieig+1))]);
      ieig = ieig+2;
    end
  end
end
if ni > 1 && Ts > 0,
  lambdades = exp(lambdades*Ts);
end
polrac = lambdades;
