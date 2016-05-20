%*****************************************************************************
%    MMCE.m :
%    Design and simulation of monovariable (SISO) and multivariable (MIMO)    
%    controllers and observers.
%
%    This software allows solving all the exercises of:
%
%              MONO- AND MULTIVARIABLE CONTROL AND ESTIMATION
%                Linear, quadratic and LMI methods
%                           Eric Ostertag
%
%    Ed. Springer, 2011
%
%*****************************************************************************
%  MMCE allows to compute state-feedback control laws or observers
%  (full-state or reduced-order, disturbance observers, Kalman filters), 
%  for five plants and an academic example used in the solved exercises of the book,
%  represented by continuous-time or discrete-time models. Three of them are
%  multivariable (MIMO) plants.
%
%  The proposed methods include pole placement, complete modal design,
%  decoupling methods, quadratic optimisation methods (LQC), a general algorithm
%  which is introduced in this book, and LMI methods.
%
%  The programm allows plotting various response curves, then switching
%  to simulation, in order to test the algorithms, without leaving it, so that
%  the user can go on checking other methods as long as he works with the same
%  plant.
%
%  Objects used by the software: 
%    called subroutines (function m-files):
%        Roppenecker, Falb_Wolovich_or_Roppenecker, Becker_Ostertag,
%        Get_Eigenvalue, Choice_Eigenvec_Paramvec, Rounding_Matrix,
%        yalmip_lmi_reg;
%    independent modules (script m-files):
%        Responses, Simulation_Scopes;
%    simulation models for Matlab R2009b (Matlab 7.9 / Simulink 7.4):
%        Universal_continuous_simulation.mdl, Universal_discrete_simulation.mdl,
%        Pendulum_and_frictions.mdl.
%
%    N.B.: due to the upwards (but not backwards) compatibility of Simulink
%    models, it is recommended to users of previous Simulink versions to use
%    instead, by renaming them as above, the following simulation files:
%       simulation models for Matlab 7.5 / Simulink 7.0.1 (R2007b):
%          Universal_continuous_simulation_75.mdl,
%          Universal_discrete_simulation_75.mdl,
%          Pendulum_and_frictions_75.mdl,
%
%  Author: Eric Ostertag
%  Last update: 18 March 2011
%  Note for the last update: the symbol used here for the time unit "minute"
%    is "min", which is accepted by the International Committee for Weights
%    and Measures for use with the International System of Units (SI); the
%    abbreviation "mn", which was used formerly to avoid confusion with the
%    mathematical abbreviation for "minimum" and occurs still in the book,
%    will be abandonned and corrected in future printings of the book.
%
%*****************************************************************************

echo on;
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%  Monovariable and Multivariable Controllers and Estimators(MMCE)  %%%%%%%%%
       %                             (MMCE.m)                              %
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
echo off;
clear; close all; format compact; format short;
disp(['Date and time: ' num2str(datestr(now)) ])
echo on;
% ***
% Choice of plant
% ***
echo off;
disp('*');
disp('Plant under study:');
disp('     1) inverted pendulum LIP 100 (Amira)');
disp('     2) three-tank setup');
disp('     3) air-flow heater: process trainer PT326 (Feedback)');
disp('     4) lateral motion of a Boeing 747');
disp('     5) magnetic tape drive');
disp('     6) plant of Exercise 3.7');
disp('     7) academic example: uncontrollable, unobservable, with');
disp('                selectable stabilizability and detectability');
proc = 0;
while proc < 1 || proc > 7,
  proc = input('Your choice (1 to 7) --> ');
  if isempty(proc), proc = 0; beep; end
end
name_y = cell(3,1);
if proc == 1,      % Inverted Pendulum LIP 100 (Amira)
  % State variables as given in electric units by the sensors
  %  x1 = cart position (xc) ; x2 = pendulum angle (theta)
  %  x3 = cart velocity_dot) ; x4 = pendulym angular velocity (theta_dot)
  %  time unit: s
  tscale = 't (s)';
  Ac = [0 0 -1.9503 0; 0 0 0 1; 0 -0.1289 -1.9148 0.0008185; 0 21.473 26.3389 -0.1362];
  Bc = [0; 0; -6.1344; 84.3027];
  Em = [1 0 0 0]';       % E = Em = input matrix for the equation noise or load disturbance: xdot = Ax + Bu + Ev
  Cmm = [1 0 0 0; 0 1 0 0];      % cart position and pendulum angle are the only measurements
  name_y{1} = 'y(1) = cart position (sensor unit, V)';
  name_y{2} = 'y(2) = pendulum angle (sensor unit, V)';
  three_measurements = input('Two measurements: position and angle (default); change to three measurements (Y/N)? [N] --> ','s');
  if isempty(three_measurements), three_measurements = 'N'; disp('N'); end
  if upper(three_measurements) == 'Y'
    Cmm(3,:) = [0 0 1 0];    % Cart velocity measured in addition
    name_y{3} = 'y(3) = cart velocity (sensor unit, V)';
  end
  Ts = 0.03;      % Sampling time of discrete-time model: Ts = 0.03 s
elseif proc == 2,    % Three-tank setup
  %  x1, x2, x3 = errors for the 3 liquid heights (m)
  %  u1 ,u2 = feeding flows (m3/min)
  %  time unit: min
  tscale = 't (min)';
  Ac = [-0.332 0.332 0; 0.332 -0.664 0.332; 0 0.332 -0.524];
  Bc = [0.764  0; 0  0; 0 0.764];
  Cmm = [1 0 0; 0 0 1]; % Only x1 and x3 are measured
  name_y{1} = 'y(1) = liquid height error of 1st tank (m)';
  name_y{2} = 'y(2) = liquid height error of 3rd tank (m)';
  Em = [0 1 0]';
  Ts = 0.1;       % Discrete model sampling time: Ts = 0.1 min
elseif proc == 3,    % Process Trainer PT326 (Feedback)
  s = tf('s');
  Gs = 1/(1+0.25*s)^2;
  Gs.inputdelay = 0.2;
  [num den] = tfdata(Gs, 'v');
  [Ar Br Cr Dr] = tf2ss(num, den);
  Ac = Ar'; Bc = Cr'; Cmm = Br'; D = Dr';  % Transformation to observability canonical form
  % x1 = output temperature (sensor units)
  % x2 = non physical quantity
  % time unit: s
  tscale = 't (s)';
  name_y{1} = 'y(1) = temperature at Process Trainer output (sensor units)';
  Em = [1 1]';
  Ts = 0.2;       % Discrete model sampling time: Ts = 0.2 s
elseif proc == 4,    % Lateral motion of a Boeing 747
  %  x1 = side-slip angle, bêta (rad); x2 = yaw rate, r (rad/s)
  %  x3 = roll rate, p(rad/s); x4 = roll angle, phi (rad)
  %  u1 = rudder angle, delta_g (rad)
  %  u2 = ailerons angle, delta_a (rad)
  %  time unit: s
	tscale = 't (s)';
  Ac = [-0.0558 -0.9968 0.0802 0.0415; 0.598 -0.115 -0.0318 0; -3.05 0.388 -0.4650 0; 0 0.0805 1 0];
  Bc = [0.00729  0.0583; -0.475  -2.009; 0.153 0.0241; 0  0];
  Cmm = [0 1 0 0; 0 0 1 0]; % Only x2 and x3 are measured
  name_y{1} = 'y(1) = yaw rate, r (rad/s)';
  name_y{2} = 'y(2) = roll rate, p(rad/s)';
  Em = [1 0 0 1]';
  Ts = 0.1;            % Discrete model sampling time: Ts = 0.1 s
elseif proc == 5,    % Magnetic tape drive
  %  x1, x2 = horizontal position of tape at capstan (mm)
  %  x3, x4 = angular rates of motor/capstan assembly (rad/s);
  %  u1, u2 = current into drive motors 1 and 2, respectively (A)
  %  time unit: s
  tscale = 't (s)';
  Ac = [0 0 -100 0; 0 0 0 100; 1.25 -1.25 -0.2 -0.2; 1.25 -1.25 -0.2 -0.2];
  Bc = [0  0; 0  0; 0.4 0; 0  0.4];
  Cmm = [0.5 0.5 0 0; -2 2 0.32 0.32]; 
  % tape position over read head, (x1+x2)/2, and tape tension Te are the outputs
  name_y{1} = 'y(1) = position of tape over read head (mm)';
  name_y{2} = 'y(2) = tape tension (N)';
  Em = [0 0; 0 0; 1 0; 0 1];
  Ts = 0.01;            % Discrete model sampling time: Ts = 0.01 s
elseif proc == 6,    % Plant of Exercise 3.7
  tscale = 't (arbitrary units)';
  s = tf('s');
  Gs = 1/(s*(1-s));
% % %   diag_37 = canon(Gs,'modal');
% % %   [Ac Bc Cmm Dc] = ssdata(diag_37);
  Ac = [0 0; 0 1];
  Bc = [1 1]';
  Cmm = [1 -1];
  name_y{1} = 'y(1) = unique output';
  Em = [1 1]';
  Ts = 1;
elseif proc == 7,    % Academic example; uncontrollable, unobservable, with
  %                             selectable stabilizability and detectability;
	tscale = 't (arbitrary units)';
  disp('Select the type of plant: ');
  disp('      0) stabilizable and detectable (default)');
  disp('      1) unstabilizable, but still detectable');
  disp('      2) undetectable, but still stabilizable');
  disp('      3) both unstabilizable and undetectable');
  NCNO = input('Your choice --> ');
  if isempty(NCNO),
    NCNO = 0;
  end
  if NCNO == 0 || NCNO < 0 || NCNO > 3
    Ac = [1 0 0; 0 -2 0; 0 0 -4];
  elseif NCNO == 1
    Ac = [1 0 0; 0  2 0; 0 0 -4];
  elseif NCNO == 2
    Ac = [1 0 0; 0 -2 0; 0 0  4];
  elseif NCNO == 3
    Ac = [1 0 0; 0  2 0; 0 0  4];
  end
  Bc = [1 0 1]';
  Cmm = [0.5 0.5 0];
  name_y{1} = 'y(1) = unique output';
  Em = [1 1 1]';
  Ts = 1;            % Discrete model sampling time: Ts = 1 s
  proc7 = ss(Ac,Bc,Cmm,0);
  [Num7 Den7] = tfdata(proc7,'v');
  zpk(proc7)
end
n = size(Ac,1); p = size(Bc,2); q = size(Cmm,1);
D = zeros(q,p); re = size(Em,2);
disp('Checking (complete) controllability of plant:')
Qc = ctrb(Ac,Bc); disp(['    Rank of Qc = '  num2str(rank(Qc))])
disp('svd(Qc) = '); svd(Qc)'
disp('Checking (complete) observability of plant:')
Qo = obsv(Ac,Cmm); disp(['    Rank of Qo = '  num2str(rank(Qo))])
disp('svd(Qo) = '); svd(Qo)'
disp('Poles and zeros of continuous-time plant, in open loop: ')
[Poles Zeros] = pzmap(Ac, Bc, Cmm(1:min(p,q),:), D(1:min(p,q),:)),
[Phi,Gamma,Cmm,D] = c2dm(Ac,Bc,Cmm,D,Ts,'zoh');
% Phi and Gamma have same dimensions as A and B
algo = 0; isfb = 0; iobs = 0; integ = 0; old_plant = ' '; isimul = 0; memsfb = 0;
Q_ipoc = []; R_ipoc = [];
while algo ~= 100,
  echo off;
  disp('*');
  disp('Choice of model:');
  disp('      C) continuous (défault)');
  if proc == 2
    disp('      D) discrete (sampling time Ts = 0.1 min)');
  else
    disp(['      D) discrete (sampling time Ts = ' num2str(Ts) ' s)']);
  end
  disp('      Q) quit program');
  plant = ' ';
  while ~(plant == 'C' || plant == 'D' || plant == 'Q')    
    plant = input('Your choice --> ','s');
    if isempty(plant)
      plant = 'C'; disp('C');
    else
      plant = upper(plant);
    end
  end
  if plant == 'Q',
    algo = 100;
  elseif plant == 'D',
    Phi, A = Phi; Gamma, B = Gamma; E = Em; T = Ts; VpVp = ' discrete-time ';
    disp('Poles and zeros of discrete-time model, in open loop: ');
    [Poles_z Zeros_z] = pzmap(A, B, Cmm(1:min(p,q),:), D(1:min(p,q),:)),
  else
    A = Ac; B = Bc; E = Em; T = 0; VpVp = ' ';
  end
  if algo ~= 100,
    displaydetails = 'N';
    rep = input('Do you want all details to be displayed (Y/N)? [N] --> ','s');
    if isempty(rep), rep = 'N'; disp('N'); end
    if upper(rep) == 'Y', displaydetails = 'Y'; end
    n = size(A,1); Cm = Cmm; integ = 0; qw = 0;
    if plant ~= old_plant,
      isfb = 0; iobs = 0; isimul = 0; old_plant = plant; KEST = [];
    end
    disp('*');
    disp('Main menu: synthesis of ...');
    disp('     A) state feedback (default)');
    disp('     B) full-state (Luenberger) observer');
    disp('     C) reduced-order observer');
    if proc == 1,
      disp('     D) disturbance observer (inverted pendulum)');
    end
    disp('     E) Kalman filter');
    disp('     F) LMI methods');
    syn = input('Your choice --> ','s');
    if isempty(syn), syn = 'A'; disp('A'); else syn = upper(syn); end
    if syn == 'A' || syn == 'F';
      disp(' Measurement-vector (output vector) components:');
      for i=1:q
        disp(['     ' name_y{i}]);
      end
      if p == q,
        if p == 1,
          disp('SISO system; the single output will be set to its unique reference value');
        else
          disp(['Number of inputs = number of outputs = ' num2str(q) ' ==>']);
          disp(['      the ' num2str(q) ' components of y will be set to their reference value']);
        end
        Sc = eye(p);
      elseif p < q,
        if p == 1,
          l_c = 0;
          while l_c < 1 || l_c > q,
            l_c = input(['Scalar control signal u ==> \n' ...
          '      index of unique component of y to be set to its reference value --> ']);
            if isempty(l_c), l_c = 0; end
          end
          qc = 1;
        else
          qc = p+1;
          while qc > p,
            l_c = input(['Row-vector with the indices of y-components to set \n' ...
            '      to a reference value (' num2str(p) ' values at most) --> ']);
            if isempty(l_c),
              l_c = 1;
            end
            qc = length(l_c); l_c = sort(l_c);
          end
        end
        Sc = zeros(qc, q);
        for i=1:qc,
          Sc(i,l_c(i)) = 1;
        end
      else
        disp(['System has p > q: suppress ' num2str(p-q) ' control input(s)']);
        beep
        break
      end
      integ = input('Integral action: 0 = without (default), 1 = with; Your choice --> ');
      if isempty(integ) || integ ~= 1,
        integ = 0; qw = 0; disp('0');
      end
      if integ == 1,
        qi = q+1;
        while qi > q,
          l_int = input(['Row-vector with indices of y-components to feed back ' ...
          'through integrator(s) \n   (' num2str(q) ' values at most) --> ']);
          qi = length(l_int); l_int = sort(l_int);
          for i=1:qi,
            while (l_int(i) < 1 || l_int(i) > q)
              l_int(i)=input(['Value out of allowed range, [1,' num2str(q) '] : reenter it ==> ']);
            end
          end
        end
        Si = zeros(qi, q);
        for i=1:qi,
          Si(i,l_int(i)) = 1;
        end
        Ci = Si*Cm;
        B = [B; zeros(qi,p)]; Cm = [Cm zeros(q,qi)]; B_ref = [zeros(n,q); Si];
        if plant == 'C',
          a22 = zeros(qi);
        else
          a22 = eye(qi);
        end
        disp('*');
        disp('*** State space representation of the augmented system (with integrator(s)):');
        A = [A zeros(n,qi); -Ci a22], B, Cm, n = n+qi; qw = qi;
      else
        Si = zeros(q,q); qi = 0;
      end
    end
    if syn == 'A', algo = 0;
      Qc = ctrb(A,B);
      disp(['Rank of the controllabiblity matrix Qc of the present system: ' num2str(rank(Qc))])
      if rank(Qc) ~= n,
        disp(' *** Warning: plant is not completely controllable: some design methods will not be available ***'); beep;
      end
      disp('Algorithm to be used:');
      if rank(Qc) == n,
        disp('      1) pole placement (MATLAB modules "acker.m" or "place.m")');
      end
      if integ == 0,
        disp('      2) simple modal control');
      end
      disp('      3) complete modal synthesis (Roppenecker''s formula)');
      if integ == 0,
        disp('      4) decoupling method (Falb-Wolovich or Roppenecker)');
      end
      if rank(Qc) == n,
        disp('      5) general formula (Becker-Ostertag)');
      end
      disp('      6) quadratic criterion: LQC (solution of the ARE or DARE), with state weighting: xT*Q*x');
      disp('      7) quadratic criterion: LQC (solution of the ARE or DARE), with output weighting: xT*CT*Q*C*x');
    elseif syn == 'B' || syn == 'C'
      algo = 0;
      if syn == 'C',
        c_shape = input('Shape of C: 1 = (0|Iq), 2 = (Iq|0); Your choice [1] --> ');
        if isempty(c_shape), c_shape = 1; disp('1'); end
        if c_shape == 1,          % Case C=(0|Iq)
          A11 = A(1:n-q, 1:n-q);     A12 = A(1:n-q, n-q+1:end);
          A21 = A(n-q+1:end, 1:n-q); A22 = A(n-q+1:end, n-q+1:end);
          B1 = B(1:n-q, :); B2 = B(n-q+1:end, :);
        else                      % Case C=(Iq|0): permutation of subscripts 1 and 2
          A22 = A(1:q, 1:q);     A21 = A(1:q, q+1:end);
          A12 = A(q+1:end, 1:q); A11 = A(q+1:end, q+1:end);
          B2 = B(1:q, :); B1 = B(q+1:end, :);
        end
      end
      if syn == 'B', Qo = obsv(A,Cm); else Qo = obsv(A11,A21); end
      if (syn == 'B' && rank(Qo) ~= n) || (syn == 'C' && rank(Qo) ~= n-q),
        disp(' *** Warning: plant or its reduced part is not completely observable: ***');
        disp(' ***                            some algorithms will not be available ***'); beep;
        disp('Algorithm to be used:');
        disp('      3) Roppenecker''s formula');
        disp('      6) quadratic criterion "LQC"');
      else
        disp('Algorithm to be used:');
        disp('      1) pole placement (MATLAB modules "acker.m" or "place.m")');
        disp('      3) Roppenecker''s formula');
        disp('      5) general formula (Becker-Ostertag)');
        disp('      6) quadratic criterion "LQC"');
      end
    elseif syn == 'D',
      algo = 10;
    elseif syn == 'E',
      algo = 11;
    elseif syn == 'F',
      Qc = ctrb(A,B);
      disp(['Rank of the controllabiblity matrix Qc of the present system: ' num2str(rank(Qc))])
      if rank(Qc) ~= n,
        disp(' *** Warning: plant is not completely controllable ***'); beep;
      end
      Qo = obsv(A,Cm);
      disp(['Rank of the observabiblity matrix Qo of the present system: ' num2str(rank(Qo))])
      if rank(Qo) ~= n
        disp(' *** Warning: plant is not completely observable ***'); beep;
      end
      algo = 12;
    end
    if algo ~= 10 && algo ~= 11 && algo ~= 12,
      while algo<1 || algo>7 && algo~=99,
        algo = input('Your choice --> ');
        if isempty(algo),
          algo = 99;
        end
      end
    end
    disp (['Plant eigenvalues and eigenvectors' VpVp ':']);
    [Vp, Lambda] = eig(A,'nobalance')
  end
  switch algo

  case 1,
    echo on;
%
% *** 1.- SISO or MIMO pole placement ("acker.m" or "place.m")
% 
    echo off;
    if syn == 'A' && rank(Qc) ~= n,
      disp('*** Plant uncontrollable: design method unavailable ***'); isfb = 0; beep;
      return
    end
    if (syn == 'B' && rank(Qo) ~= n) || (syn == 'C' && rank(Qo) ~= n-q),
      disp('*** Plant unobservable: design method unavailable ***'); iobs = 0; beep;
      return
    end
    if syn == 'A',
      disp(['N.B.: To be able to use "place.m", eigenvalues order of multiplicity must be'...
          ' <= ' num2str(p) ' (number of inputs)']);
    else
          disp(['N.B.: To be able to use "place.m", eigenvalues order of multiplicity must be'...
          ' <= ' num2str(q) ' (number of outputs)']);
    end
    if syn == 'A' || syn == 'B',
      lambda = Get_Eigenvalue(n,T); nb = n;
      if syn == 'A', rB = rank(B); else rB = rank(Cm); end
    end
    if syn == 'C',
      lambda = Get_Eigenvalue(n-q,T); rB = rank(A21); nb = n-q;
    end
    ps = sort(lambda); pmult = 0;
    for i=1:nb-rB
      imax = min(nb,i+rB);
      if all(ps(i:imax) == ps(i)),
        pmult = 1;
      end
    end
    if (syn == 'A' && p == 1) || (syn == 'B' && q == 1) || (syn == 'C' && n-q == 1),
      if pmult == 0,
        method = input('Choice of MATLAB module: 1 = acker.m, 2 = place.m; Your choice [1] --> ');
        if isempty(method), method = 1; disp('1'); end
      else
        method = 1;
      end
    else
      if pmult == 0,
        method = 2;
      else
        method = 0;
      end
    end
    if method == 0,
      if syn == 'A',
        disp('There are poles with multiplicity > p ==> make another choice');
      else
        disp('There are poles with multiplicity > q (resp., n-q) ==> make another choice');
      end
      iobs = 0;
      algo = 99;
      beep;
    else
      if method == 1,
        if syn == 'A',
          Lpp = acker(A,B,lambda); L = Lpp, isfb = 1;
        end
        if syn == 'B',
          Gpp = acker(A',Cm',lambda); G = Gpp', iobs = 1;
        end
        if syn == 'C',
          Gpp = acker(A11',A21',lambda); Gr = Gpp', iobs = 1;
        end
      elseif method == 2,
        if syn == 'A',
          Lpp = place(A,B,lambda); L = Lpp, isfb = 1;
        end
        if syn == 'B',
          Gpp = place(A',Cm',lambda); G = Gpp', iobs = 1;
        end
        if syn == 'C',
          Gpp = place(A11',A21',lambda); Gr = Gpp', iobs = 1;
        end
      end
    end
    
  case 2,
    echo on;
%
% *** 2.- Design of a simple modal control
%
    echo off;
    if integ == 1,
      disp('*** Design method available only without integrator ***'); isfb = 0; beep;
    else
      Wp = inv(Vp); % the rows of Wp are the left eigenvectors of A
      disp(['Only ' num2str(p) ' eigenvalue(s) will apply for shifting']);
      j = 0;
      for i=1:n
        shift = input(['Shift eigenvalue lambda' num2str(i) ' = ' num2str(Lambda(i,i)) ' (Y/N)? [N] --> '],'s');
        if upper(shift) == 'Y'
          j = j+1; lambda = input('      New value (continuous) --> ');
          if plant == 'D',
            lambda = exp(lambda*Ts);
          end
          Ti_p(j,:) = Wp(i,:); lambdadiff(j) = Lambda(i,i)-lambda;
        end
        if j == p,
          break;
        end
      end
      if j>0,
        Bp_hat = Ti_p*B; L_modal = inv(Bp_hat)*diag(lambdadiff)*Ti_p;
        L = L_modal, isfb = 1;
      end
    end

  case 3,
    echo on;
%
% *** 3.- Design of a complete modal control (Roppenecker's formula)
%
    echo off;
    if syn == 'A',
      [Lrop,lambda_cre] = Roppenecker(A,B,Lambda,Vp,[],T); L = Lrop, isfb = 1;
    end
    if syn == 'B',
      [Grop,lambda_obs] = Roppenecker(A',Cm',Lambda,Vp,[],T); G = Grop',iobs = 1;
    end
    if syn == 'C',
      [Grop,lambda_obs] = Roppenecker(A11',A21',Lambda,Vp,[],T);Gr=Grop', iobs = 1;
    end
    
  case 4,
    echo on;
%
% *** 4.- Decoupling method, according to Falb-Wolovich or Roppenecker
%
    echo off;
    if integ == 1,
      disp('*** Design only without integrator ***'); isfb = 0; beep;
    else
      [L, M] = Falb_Wolovich_or_Roppenecker(A,B,Cm,Lambda,Vp,T);
      if L == 0,
        isfb = 0;
      else
        L, M,
% Display decoupled transfer matrix and
% check the zero's position(s) (RHP or not ?)
        A_CL = Rounding_Matrix(A-B*L); G_CL = Rounding_Matrix(Cm(1:p,:)*inv(-A_CL)*B);
        B_CL = B*M*Sc; sys_CL = ss(A_CL, B_CL, Cm, D, T);
        disp('Closed-loop transfer matrix:');
        Gw = tf(sys_CL)
        disp('Closed-loop zero(s): ');
        Z = zero(sys_CL)
        disp(' ');
        disp('Minimal realization (suppress uncontrollable and unobservable modes');
        Msys_CL = minreal(sys_CL);  % suppress uncontrollable and unobservable modes
        disp(' ');
        isfb = 2;
      end
    end

  case 5,
    echo on;
%
% *** 5.- General formula (Becker-Ostertag)
%
    echo off;
    if syn == 'A' && rank(Qc) ~= n,
      disp('*** Plant uncontrollable: design method unavailable ***'); isfb = 0; beep;
      return
    end
    if (syn == 'B' && rank(Qo) ~= n) || (syn == 'C' && rank(Qo) ~= n-q),
      disp('*** Plant unobservable: design method unavailable ***'); iobs = 0; beep;
      return
    end
    if syn == 'A',
      [Lbo,lambda_cre] = Becker_Ostertag(A,B,Lambda,[],T); L = Lbo, isfb = 1;
    end
    if syn == 'B',
      [Gbo,lambda_obs] = Becker_Ostertag(A',Cm',Lambda,[],T); G = Gbo', iobs = 1;
    end
    if syn == 'C',
      [Gbo,lambda_obs] = Becker_Ostertag(A11',A21',Lambda,[],T); Gr = Gbo', iobs = 1;
    end
    
  case {6,7},
    echo on;
%
% *** 6-7.- Design by quadratic criterion (LQC), with state or output weighting
%
    echo off;
    Q = []; diag_Q = []; R = []; diag_R = [];
    if algo == 6,
      disp('*** State weighting: xT*Q*x ***');
      if syn == 'C'
        Cw = eye(n-q);
      else
        Cw = eye(n);
      end
    else
      disp('*** Output weighting: yT*Q*y ***');
      Cw = [Cm; zeros(qw,n-qw) eye(qw)];
    end
    nw = size(Cw,1);
    if isempty(Q_ipoc) && nw == 1, 
      Q = input('Q is scalar: enter single (positive) value --> ');
    else
      diag_yes = input('Choose matrix Q in diagonal form (Y/N)? [Y] --> ','s');
      if isempty(diag_yes)
        diag_yes = 'Y';
      else
        diag_yes = upper(diag_yes(1));
      end
      while size(Q,1) ~= nw, 
        if diag_yes == 'Y'
          diag_Q = input(['Row-vector containing the ' num2str(nw) ' diagonal elements qii of Q --> ']);
          Q = diag(diag_Q);
        else
          Q = input(['Matrix Q (' num2str(nw) 'x' num2str(nw) ') = ']);
        end
      end
    end
    Qw = Cw'*Q*Cw;
    if syn == 'A',
      nR = p;
    elseif syn == 'B',
      nR = q;
    elseif syn == 'C',
      nR = n-q;
    end
    if isempty(R_ipoc) && nR == 1
      R = input('R is scalar: enter single (positive) value --> ');
    else
      diag_yes = input('Choose matrix R in diagonal form (Y/N)? [Y] --> ','s');
      if isempty(diag_yes)
        diag_yes = 'Y';
      else
        diag_yes = upper(diag_yes(1));
      end
      while size(R,1) ~= nR,
        if diag_yes == 'Y'
          diag_R = input(['Row-vector containing the ' num2str(nR) ' diagonal elements rii of R --> ']);
          R = diag(diag_R);
        else
          R = input(['Matrix R (' num2str(nR) 'x' num2str(nR) ') = ']);
        end
      end
    end
    if plant == 'C',
      if syn == 'A',
        [Lopt, P, Lambda_CL] = lqr(A,B,Qw,R); L = Lopt, isfb = 1;
        if displaydetails == 'Y'
          disp('Solution of the associated ARE: '); P
          disp('Closed-loop eigenvalues: '); Lambda_CL
        end
      end
      if syn == 'B',
        [Gopt, Sigmac, Lambda_obs] = lqr(A',Cm',Qw,R); G = Gopt', iobs = 1;
        if displaydetails == 'Y'
          disp('Solution of the associated dual ARE: '); Sigmac
          disp('Observer eigenvalues: '); Lambda_obs
        end
      end
      if syn == 'C',
        [Gopt, Sigmac, Lambda_obs] = lqr(A11',A21',Qw,R); Gr = Gopt', iobs = 1;
        if displaydetails == 'Y'
          disp('Solution of the associated dual ARE: '); Sigmac
          disp('Observer eigenvalues: '); Lambda_obs
        end
      end
    elseif plant == 'D',
      if syn == 'A',
        [Lopt, P, Lambda_CL] = dlqr(A,B,Qw,R); L = Lopt, isfb = 1;
        if displaydetails == 'Y'
          disp('Solution of the associated DARE: '); P
          disp('Closed-loop eigenvalues: '); Lambda_CL
        end
      end
      if syn == 'B',
        [Gopt, Sigmad, Lambda_obs] = dlqr(A',Cm',Qw,R); G = Gopt', iobs = 1;
        if displaydetails == 'Y'
          disp('Solution of the associated dual DARE: '); Sigmad
          disp('Observer eigenvalues: '); Lambda_obs
        end
      end
      if syn == 'C',
        [Gopt, Sigmad, Lambda_obs] = dlqr(A11',A21',Qw,R); Gr = Gopt', iobs = 1;
        if displaydetails == 'Y'
          disp('Solution of the associated dual DARE: '); Sigmad
          disp('Observer eigenvalues: '); Lambda_obs
        end
      end
    end

  case 10,
    echo on;
%
% *** 10.- Design of an observer for the discrete-time plant augmented with the constant
%          disturbance (dry friction), for the inverted pendulum, only.
%
    echo off;
    if plant == 'C',
      disp('*** Design in discrete time, exclusively ***'); algo = 0; beep;
    elseif plant == 'D',
      if proc ~= 1,
        disp('*** Inverted pendulum only ***'); algo = 0; beep;
      else
% Continuous-time model augmented with the constant disturbance (dry friction)
        A_tild = [Ac Bc; zeros(1,n) zeros(1,p)], B_tild = [Bc; zeros(1,p)],
        C_tild = [Cm zeros(q,p)], D_tild = D,
% Corresponding augmented discrete-time model
        [A,B,Cm,D] = c2dm(A_tild,B_tild,C_tild,D_tild,Ts,'zoh'); Phi_tild = A, Gamma_tild = B,
        n = size(A,1); q = size(Cm,1);
        [Vp_aug,Lambda_aug] = eig(A,'nobalance');
        typeobs = input('Observer type: 1 = identity (full-state), 2 = reduced-order; Your choice [1] --> ');
        if isempty(typeobs), typeobs = 1; disp('1'); end
        disp('Choice of design method:');
        disp('                 1 = pole placement (module "place.m") (default)');
        disp('                 3 = Roppenecker''s formula');
        disp('                 5 = Becker-Ostertag''s formula');
        method = input('Your choice --> ');
        if isempty(method)
          method = 1; disp('1');
        end
        if typeobs == 1
          switch method;
          case 3
            [Lobs1,lambda_obs] = Roppenecker(A',Cm',Lambda_aug,Vp_aug,[],Ts);
            G = Lobs1'
          case 5
            [Lobs2,lambda_obs] = Becker_Ostertag(A',Cm',Lambda_aug,[],Ts);
            G = Lobs2'
          otherwise
            disp(['Enter n = ' num2str(n) ' (continuous) eigenvalues ' ...
                  ' for the identity (full-state) observer ===> ...']);
            disp(['N.B.: order of multiplicity <= ' num2str(q) ' (number of outputs)']);
            lambda_obs = Get_Eigenvalue(n,Ts);
            G = (place(A',Cm',lambda_obs))'
          end
        elseif typeobs == 2
  % Partitioning of the matrices between measured parts: x1 to xq,
  % and estimated parts: x(q+1) to xn, plus perturbation
          c_shape = 2;
          A22 = A(1:q, 1:q); A21 = A(1:q, q+1:end);
          A12 = A(q+1:end, 1:q); A11 = A(q+1:end, q+1:end);
          B2 = B(1:q, :); B1 = B(q+1:end, :);
          switch method;
          case 3
            [Grop,lambda_obs] = Roppenecker(A11',A21',Lambda_aug,Vp_aug,[],Ts); Gr = Grop'
          case 5
            [Gbo,lambda_obs] = Becker_Ostertag(A11',A21',Lambda_aug,[],Ts); Gr = Gbo'
          otherwise
            disp(['Enter r=' num2str(n-q) ' (continuous) eigenvalues' ...
                  ' for the reduced-order observer ===> ...']);
            disp(['N.B.: order of multiplicity <= ' num2str(q) ' (number of outputs)']);
            lambda_obs = Get_Eigenvalue(n-q,Ts);
            Gr = (place(A11',A21',lambda_obs))'
          end
        end
        iobs = 2;
      end
    end

  case 11,
    echo on;
%
% *** 11.- Design of a Kalman filter ("kalman.m")
%
    echo off;
    diag_Qbar = []; diag_Rbar = [];
% a) Adaptation of the LTI model of the plant to the structure of "kalman.m" (see on line help)
    ext_mod = ss(Phi,[Gamma Gamma],Cm,[D zeros(q,p)],Ts);
% b) Calculation of the Kalman filter
    disp('Variance matrix of equation noise (process noise):');
    if p == 1, 
      diag_Qbar = input('Q_bar is scalar: enter single (positive) value --> ');
    else
      while length(diag_Qbar) ~= p, 
        diag_Qbar = input(['Row-vector containing the ' num2str(p) ' diagonal elements qii of Q_bar --> ']);
      end
    end
    disp('Variance matrix of measurement noise:');
    if q == 1
      diag_Rbar = input('R_bar is scalar: enter single (positive) value --> ');
    else
      while length(diag_Rbar) ~= q, 
        diag_Rbar = input(['Row-vector containing the ' num2str(q) ' diagonal elements rii of R_bar --> ']);
      end
    end
    Qbar = diag(diag_Qbar); Rbar = diag(diag_Rbar);
    [KEST,PhiK,Sigma_moins,K_] = kalman(ext_mod,Qbar,Rbar);
    disp('Steady-state Kalman gain:'); K_,
    disp('Matrices of the Kalman filter state equation (full-state observer form):');
    E = Gamma; G = PhiK, iobs = 1;
    
  case 12,
    echo on;
%
% *** 12.- Use of LMI methods to solve some control problems
%
    echo off;
    lmi_method = 0;
    options = sdpsettings();
    opts = sdpsettings(options,'verbose',0,'warning',1);
% % %     opts = sdpsettings(opts,'solver','sedumi','sedumi.eps',1e-12);
    disp('LMI method to be used:');
    disp('      1) stabilizability and stabilization by state-feedback');
    disp('      2) detectability and detection by observer');
    disp('      3) alpha-stabilization by state-feedback')
    disp('      4) LQ regulator design');
    disp('      5) state-feedback design with regional pole placement');
    disp('      6) LQ regulator design with regional pole placement');
    disp('      7) inverse problem of optimal control (ipoc)');
    while lmi_method < 1 || lmi_method > 7 && lmi_method~=99,
      lmi_method = input('Your choice --> ');
      if isempty(lmi_method),
        lmi_method = 99;
        algo = 99;
      end
    end
    if lmi_method == 1      % stabilizability and stabilization by state-feedback
      disp('*');
      clear Y W; 
      Y = sdpvar(n,n);   % Y = P^-1
      W = sdpvar(p,n);   % W = L*P^-1 = L*Y ===> L = W*Y^-1
      F = set(Y > 0);
      if plant == 'C'
        F = F + set(A*Y + Y*A' - B*W - W'*B' < 0);
      elseif plant == 'D'
        F = F + set([-Y       Y*A'-W'*B'
                     A*Y-B*W  -Y        ] < 0);
      end
      F = F + set(-100 < recover(depends(F))< 100);
      sol = solvesdp(F,[],opts);
      disp(['Solving the stabilizability LMI: ' num2str(sol.info) ])
      P = inv(double(Y));
      Lstab_lmi = double(W)*P;
      if displaydetails == 'Y'
        disp('Lyapunov matrix solution of the LMI: '); P
      end
      trace_P = trace(P);
      if abs(trace_P) > 1e8
        disp('*** Plant is not stabilizable by state feedback ***'); beep
        isfb = 0;
      else
        L = Lstab_lmi, isfb = 1;
      end
    elseif lmi_method == 2    % detectability and detection by observer
      disp('*');
      clear Y V; 
      Y = sdpvar(n,n);   % Y = P^-1
      V = sdpvar(n,q);   % V = P^-1*G = Y*G ===> G = Y^-1*V
      F = set(Y > 0);
      if plant == 'C'
        F = F + set(A'*Y + Y*A - V*Cm - Cm'*V' < 0);
      elseif plant == 'D'
        F = F + set([-Y           Y*A-V*Cm
                     A'*Y-Cm'*V'  -Y      ] < 0);
      end
      F = F + set(-100 < recover(depends(F))< 100);
      sol = solvesdp(F,[],opts);
      disp(['Solving the detectability LMI: ' num2str(sol.info) ])
      P = inv(double(Y));
      Gdet_lmi = P*double(V);
      if displaydetails == 'Y'
        disp('Lyapunov matrix solution of the LMI: '); P
      end
      trace_P = trace(P);
      if abs(trace_P) > 1e8
        disp('*** Plant is not detectable by an observer ***'); beep
        iobs = 0;
      else
        G = Gdet_lmi, iobs = 1;
      end
      
    elseif lmi_method == 3    % alpha-stabilization by state-feedback
      disp('*');
      clear Y W; 
      Y = sdpvar(n,n);   % Y = P^-1
      W = sdpvar(p,n);   % W = L*P^-1 = L*Y ===> L = W*Y^-1
      F = set(Y > 0);
      if plant == 'C'
        alpha = input('  Value of alpha > 0:  ');
        if isempty(alpha), alpha = -1; end
        while length(alpha)~=1 | ischar(alpha) | alpha <= 0 | ~isreal(alpha),
          beep; alpha = input('  Enter a positive real number:  ');
          if isempty(alpha), alpha=-1; end
        end
        F = F + set(A*Y + Y*A' + 2*alpha*Y - B*W - W'*B' <= 0);
      elseif plant == 'D'
        alpha = input('  Value of alpha > 1:  ');
        if isempty(alpha), alpha = -1; end
        while length(alpha)~=1 | ischar(alpha) | alpha <= 1 | ~isreal(alpha),
          beep; alpha = input('  Enter a real number greater than 1:  ');
          if isempty(alpha), alpha=-1; end
        end
        F = F + set([-Y             alpha*Y*A'-W'*B'
                     alpha*A*Y-B*W        -Y        ] < 0);
      end
      sol = solvesdp(F,-alpha,opts);
      disp(['Solving the LMI for alpha-stabilization: ' num2str(sol.info) ])
      P = inv(double(Y));
      if displaydetails == 'Y'
        disp('Lyapunov matrix solution of the LMI: '); P
      end
      Lalpha_stab_lmi = double(W)*P;
      L = Lalpha_stab_lmi, isfb = 1;
      
    elseif lmi_method == 4 || lmi_method == 6         % LQ regulator design
      Q = []; diag_Q = []; R = []; diag_R = [];
      disp('*** State weighting: xT*Q*x ***');
      Cw = eye(n);
      nw = size(Cw,1);
      disp('Choose matrix Q in diagonal form');
      while size(Q,1) ~= nw, 
        diag_Q = input(['Row-vector containing the ' num2str(nw) ' diagonal elements qii of Q --> ']);
        Q = diag(diag_Q);
        diag_H = sqrt(diag_Q);
        H = diag(diag_H);
      end
      nR = p;
      disp('Choose matrix R in diagonal form');
      while size(R,1) ~= nR,
        diag_R = input(['Row-vector containing the ' num2str(nR) ' diagonal elements rii of R --> ']);
        R = diag(diag_R);
      end
      clear Y W;
      Y = sdpvar(n,n);   % Y = P^-1
      W = sdpvar(p,n);   % W = L*P^-1 = L*Y ===> L = W*Y^-1
      gamma = sdpvar(1,1);
      F = set(Y > 0);
      if plant == 'C'
        F = F + set([A*Y+Y*A'-B*W-W'*B'  Y*H'        W'
                      H*Y              -eye(n)      zeros(n,p)
                      W                zeros(p,n) -inv(R)    ] < 0);
      elseif plant =='D'
        F = F + set([-Y        Y*A'-W'*B'  Y*H'        W'
                      A*Y-B*W      -Y      zeros(n,n) zeros(n,p)
                      H*Y      zeros(n,n) -eye(n)     zeros(n,p)
                      W        zeros(p,n) zeros(p,n)  -inv(R)    ] < 0);
      end
      if lmi_method == 6          % LQ regulator design with regional pole placement
        display = 'Off';
        [Fout,Lreg,Mreg] = yalmip_lmi_reg(F,A,B,Y,W,display);
        F = Fout;
      end
      F = F + set([gamma*eye(n) eye(n)
                   eye(n)       Y     ] > 0);
      sol = solvesdp(F,gamma,opts);
      disp(['Solving the LQC design LMI: ' num2str(sol.info) ])
      P = inv(double(Y));
      if displaydetails == 'Y'
        disp('Lyapunov matrix solution of the LMI: '); P
        disp(['Minimum value obtained for gamma: ' num2str(double(gamma)) ])
      end
      Lopt_lmi = double(W)*P;
      L = Lopt_lmi, isfb = 1;
      
    elseif lmi_method == 5    % state-feedback design with regional pole placement
      disp('*');
      clear Y W; 
      Y = sdpvar(n,n);   % Y = P^-1
      W = sdpvar(p,n);   % W = L*P^-1 = L*Y ===> L = W*Y^-1
      F = set(Y > 0);
      display = 'Off';
      [Fout,Lreg,Mreg] = yalmip_lmi_reg(F,A,B,Y,W,display);
      F = Fout;
      sol = solvesdp(F,[],opts);
      disp(['Solving the state-feedback with regional pole plecement: ' num2str(sol.info) ])
      P = inv(double(Y));
      if displaydetails == 'Y'
        disp('Lyapunov matrix solution of the LMI: '); P
      end
      Lpconf_lmi = double(W)*P;
      L = Lpconf_lmi, isfb = 1;
      
    elseif lmi_method == 7    % inverse problem of optimal control
      disp('*');
      if memsfb == 1
        old_L = input('Do you want to use the last calculated L matrix (Y/N)? [Y] --> ','s');
      else
        old_L = 'N';
      end
      if upper(old_L) ~= 'Y'
        L = input('Matrix L = ');
      end
      clear P P1 Q R; 
      P = sdpvar(n,n);
      P1 = sdpvar(n,n);
      Q = sdpvar(n,n);
      R = sdpvar(p,p);
      F = set(P >= 0) + set(P1 > 0);
      F = F + set(Q >= 0) + set(R > 0);
      if plant == 'C'
        F = F + set((A-B*L)'*P + P*(A-B*L) + L'*R*L + Q == 0);
        F = F + set(B'*P - R*L == 0);
        F = F + set(A'*P1 + P1*A < Q);
      elseif plant == 'D'
        F = F + set((A-B*L)'*P*(A-B*L) -P + L'*R*L + Q == 0);
        F = F + set(B'*P*A - (R+B'*P*B)*L == 0);
        F = F + set(A'*P1*A - P1 < Q);
      end
      sol = solvesdp(F,[],opts);
      disp(['Solving the inverse problem of optimal control: ' num2str(sol.info) ])
      P = double(P);
      if displaydetails == 'Y'
        disp('Lyapunov matrix solution of the LMI: '); P
      end
      Q_ipoc = double(Q)
      R_ipoc = double(R)
      algo = 99;
    end
    
  case 100,
    break
  otherwise
    disp('*** Make a new choice ***'); beep;
  end

%
% Completing the designs (observer and closed-loop system)
%
  if algo ~= 99,
    if iobs == 1 || iobs == 2,
      if syn == 'B' || (syn == 'D' && typeobs == 1) || syn == 'E' || syn == 'F'
        F = A-G*Cm, J = B; H = eye(n); K = zeros(n,q);  % Completing observer design
      elseif syn == 'C' || (syn == 'D' && typeobs == 2),
        F = A11-Gr*A21, G = F*Gr+A12-Gr*A22, J = B1-Gr*B2
        if c_shape == 1
          H = [eye(n-q); zeros(q,n-q)], K = [Gr; eye(q)]
        else
          H = [zeros(q,n-q); eye(n-q)], K = [eye(q); Gr]
        end
      end
    end
    if isfb == 1 || isfb == 2,    % Completing state-feedback and M matrix
      memsfb = 1;     % Memorize that a state feedback has already been designed
      if integ == 1
        L1 = L(:,1:n-qi), L2 = L(:,n-qi+1:n),
      else
        L1 = L; L2 = zeros(p,q);
      end
      if plant == 'C',
        A_CL1 = Ac - Bc*L1; G_CL = Sc*Cmm*inv(-A_CL1)*Bc;
      else
        A_CL1 = Phi - Gamma*L1; G_CL = Sc*Cmm*inv(eye(size(A_CL1))-A_CL1)*Gamma;
      end
      if isfb ~= 2,
        if rank(G_CL) < p,
          M = eye(p),
        else
          M = inv(G_CL),
        end
      end
      disp('Closed-loop poles and zeros: ')
      if integ == 1,
        A_CL = A-B*L;
        B_CL = B*M*Sc + B_ref;
        C_CL = Cm;
        D_CL = zeros(q,q);
      else
        A_CL = A_CL1;
        B_CL = B*M*Sc;
        C_CL = Cmm;
        D_CL = zeros(q,q);
      end
      sys_CL = ss(A_CL, B_CL, C_CL, D_CL, T);
      [Poles Zeros] = pzmap(A_CL,B_CL,C_CL(1:min(p,q),:),D_CL(1:min(p,q),:))
      CL_eigv = eig(A_CL);
      badEVs = [];
      if plant == 'C'
        for i = 1:length(CL_eigv),
          if real(CL_eigv(i)) > 0, badEVs = [badEVs CL_eigv(i)]; end
        end
        if ~isempty(badEVs),
          disp(['!!! Closed-loop unstable eigenvalues (Re > 0) : ' num2str(badEVs)]); beep
        end
      elseif plant =='D'
        for i = 1:length(CL_eigv),
          if abs(CL_eigv(i)) > 1, badEVs = [badEVs CL_eigv(i)]; end
        end
        if ~isempty(badEVs),
          disp(['!!! Closed-loop unstable eigenvalues (modulus > 1) : ' num2str(badEVs)]); beep
        end
      end
%
% Parameters of the various plants used in the plots and the simulations
%
      if proc == 1      % Inverted pendulum LIP 100 (Amira)
        x0 = [0; -2; 0.5; 0]; yr = [1; 0]; start_reference = [0; 3]; tfinal = 10; max_samples = 300;
        p_load = 0.5; start_p_load = 5;
        p_measurements = [0.05; 0.1]; start_p_measurements = 7;
        if plant == 'D', simulation_horizon = 300*Ts; else simulation_horizon = 10; end
        eqn_pwr = 0.01*Ts; eqn_seed = 23341;
        meas_pwr = [0.02; 0.01]*Ts ; meas_seed = [14223; 42313];
        if upper(three_measurements) == 'Y'
          yr = [yr; 0]; start_reference = [start_reference; 0];
          p_measurements = [p_measurements; 0];
          meas_pwr = [meas_pwr; 0.01*Ts] ; meas_seed = [meas_seed; 31424];
        end
      end
      if proc == 2,    % Three-tank setup
        x0 = [-0.8; -0.5; -0.2]; yr = [0.4; 0.8]; start_reference = [0; 10]; tfinal = 30; max_samples = 300;
        p_load = 0.5; start_p_load = 15;
        p_measurements = [1; 1]; start_p_measurements = [20; 25]; simulation_horizon = 300*Ts;
        eqn_pwr = 0.01*Ts; eqn_seed = 23341;
        meas_pwr = [0.02; 0.01]*Ts ; meas_seed = [14223; 42313];
      end
      if proc == 3,    % Process Trainer PT326 (Feedback)
        x0 = [10; 0]; yr = 1; start_reference = 0; tfinal = 20; max_samples = 100;
        p_load = 0.5; start_p_load = 15;
        p_measurements = 1; start_p_measurements = 20; simulation_horizon = 100*Ts;
        eqn_pwr = 0.01*Ts; eqn_seed = 23341; meas_pwr = 0.01*Ts ; meas_seed = 14223;
      end
      if proc == 4,    % Lateral motion of a Boeing 747
        x0 = [0.15; 0; 0; 0.5]; yr = [1; 1]; start_reference = [0; 10]; tfinal = 30; max_samples = 100;
        p_load = 0.5; start_p_load = 15;
        p_measurements = 1; start_p_measurements = 20; simulation_horizon = 100*Ts;
        eqn_pwr = 0.01*Ts; eqn_seed = 23341;
        meas_pwr = [0.02; 0.01]*Ts ; meas_seed = [14223; 42313];
      end
      if proc == 5,    % Magnetic tape drive
        x0 = [-1; 2; 0; 0]; yr = [1; 2]; start_reference = [0; 0]; tfinal = 2; max_samples = 100;
% % %         p_load = [1 2]; start_p_load = [50*Ts 60*Ts];
        p_load = [1; 0]; start_p_load = [50*Ts; 60*Ts];
        p_measurements = 0.2; start_p_measurements = 75*Ts; simulation_horizon = 100*Ts;
        eqn_pwr = 0.01*Ts; eqn_seed = [23341; 34412];
        meas_pwr = [0.02^2; 0.01^2]*Ts ; meas_seed = [14223; 42313];
      end
      if proc == 6,    % Plant of Exercise 3.7
        x0 = [1; 0]; yr = 1; start_reference = 0; tfinal = 30; max_samples = 30;
        p_load = 0.5; start_p_load = 15;
        p_measurements = 1; start_p_measurements = 20; simulation_horizon = 10*Ts;
        eqn_pwr = 0.01*Ts; eqn_seed = 23341; meas_pwr = 0.01*Ts ; meas_seed = 14223;
      end
      if proc == 7,    % Academic example (non controllable, non observable)
        x0 = [-1; 2; 0]; yr = 1; start_reference = 0; tfinal = 10; max_samples = 15;
        p_load = 1; start_p_load = 5*Ts;
        p_measurements = -0.5; start_p_measurements = 10*Ts; simulation_horizon = 15*Ts;
        eqn_pwr = 0.01*Ts; eqn_seed = 23341;
        meas_pwr = [0.02^2; 0.01^2]*Ts ; meas_seed = [14223; 42313];
      end
      xx0 = x0;
      if T == 0,
        t = 0:tfinal/200:tfinal;
      else
        t = 0:T:max_samples*T;
      end
%
%  Closed-loop step-response plots
%
      isfb = 0;    % ready to design another controller and/or an observer
      isimul = 1;  % ready for simulation
      Responses
    end
%
% Simulations, if desired
%
    if isimul == 0
      disp('**************** WARNING! ************************')
      disp('A controller is required to perform simulations.')
      disp('If you want to do so, please type ''C'' or ''D'' at the next prompt,')
      disp('according to your model, and design a controller.')
      disp('Otherwise, type ''Q'' to quit')
    else
      sim = input('Do you want to simulate the control loop (Y/N)? [Y] --> ','s');
      if isempty(sim)
        sim = 'Y'; disp('Y');
      else
        sim = sim(1);
      end
      if upper(sim) == 'Y',
        if integ == 1 || iobs == 2
          if plant == 'C'
            A = Ac; B = Bc;
          else
            A = Phi; B = Gamma;
          end
        end
        n = size(A,1); Cm = Cmm; q = size(Cm,1); z0 = 0;
        if iobs == 0,
          F = zeros(n); G=zeros(n,q); J=zeros(n,p); H=zeros(n); K=zeros(n,q);
        end
        x0_nul = input('Zero initial state of plant (Y/N)? [Y] --> ','s');
        if isempty(x0_nul)
          x0_nul = 'Y'; disp('Y');
        else
          x0_nul = x0_nul(1);
        end
        if upper(x0_nul) == 'Y',
          x0 = zeros(n,1);
        else
          x0 = xx0
        end
        if plant == 'C',
          disp('Simulate with "Universal_continuous_simulation.mdl" ===>');
        elseif plant == 'D',
          if iobs == 2,
            disp('Simulate with "Pendulum_and_frictions.mdl" ===>');
          else
            disp('Simulate with "Universal_discrete_simulation.mdl" ===>');
          end
        end
        disp('(After simulation, type "return" to return to program)')
        keyboard
        Simulation_Scopes      % Plotting oscillograms obtained during simulation
      end
    end
  end
end