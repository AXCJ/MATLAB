function [L,lambdaCL] = Becker_Ostertag(A,B,Lambda,lambdades,Ts)
%
% Design of a multivariable controller or observer by the formula of
% Becker-Ostertag
%
%  [L,lambdaCL] = Becker_Ostertag(A,B,Lambda,[],Ts)  begins acquiring the desired
%  closed-loop eigenvalues. Two cases are possible, according to the choice
%  made for the canonical form of the closed-loop system matrix:
%
%    * case 1: controllability canonical form 
%    * case 2: diagonal canonical form
%
%  In case 1, an arbitrary number of identical eigenvalues can be chosen.
%  In case 2, a given eigenvalue must not have an order of multiplicity greater than 
%  the number of plant inputs. The entered eigenvalues are returned to the calling program
%  by means of the vector lambdaCL:
%    if Ts = 0, the entered eigenvalues are returned without modification;
%    if Ts > 0, the entered eigenvalues are converted to their discrete equivalent
%              (at sampling period Ts) by exp(lambda*Ts), before being returned.
%  The input matrix argument Lambda is the diagonal matrix of the eigenvalues of A, thus of the
%  open-loop system, for comparison purposes with the eigenvalues prescribed for the 
%  closed loop. The subroutine prompts then the user to enter the invariant parameter vectors,
%  correlated or not with the choice of the closed-loop eigenvectors, and computes then, 
%  with the formula of Becker-Ostertag, the state-feedback matrix L which assigns the 
%  entered eigenvalues to the matrix A-B*L.
%
%  [L,lambdaCL] = Becker_Ostertag(A,B,Lambda,lambdades,Ts) does the same thing, except that
%  the eigenvalues desired for the closed loop are given as input argument by means of the
%  vector lambdades. They are thus not entered anymore by the user.
%
%  N.B.: by duality, this subroutine is also used for the computation of observers
%  (full or reduced-order).
%
%  Called subroutines: Get_Eigenvalue, Choice_Eigenvec_Paramvec, Rounding_Matrix.
%
%  Author: E. Ostertag  10 February 2010
%  Last update: 12 August 2010
%

n = size(A,1); p = size(B,2);
ni = nargin;
if ni<3,
  disp('***Not enough arguments***'); beep; L=0; lambdaCL=0;
  return;
end;
%
% 1. Choice of the canonical model of the system matrix having the desired eigenvalues
%
disp('Canonical form chosen for the closed-loop system matrix:');
disp('                 1 = controllability form');
disp('                 2 = diagonal form');
form = input('Your choice --> ');
if isempty(form) || form ~= 2,
  form = 1;
end
D = zeros(n);
while rank(D) ~= n,
%
% 2. Choice of the eigenvalues desired for the closed-loop system, except if they belong to the
%    input arguments
%
  if isempty(lambdades),
    if form == 2,
      disp('Enter arbitrary eigenvalues, with an order of multiplicity <= number of inputs')
      lambdades = Get_Eigenvalue(n,0,p);
    else
      disp('Enter arbitrary eigenvalues, distinct or not, real or complex)');
      lambdades = Get_Eigenvalue(n);
    end
  end
  if ni == 5 && Ts > 0,
    lambdades = exp(lambdades*Ts);
  end
  polchar = real(poly(lambdades)); polchar = polchar/polchar(1);
  if form == 2,
    Acan = diag(lambdades);
  else
    [Acan, Bcan, Ccan, Dcan] = tf2ss(1,polchar);
  end
  Qc = ctrb(A,B);
%
% 3. Choice of the eigenvectors and of the parameter vectors
%
  [V, P] = Choice_Eigenvec_Paramvec(A,B,n,p,Lambda,lambdades);
%
% 4. Application of the Becker-Ostertag formula
%
  q0 = polyvalm(polchar,A);
  for i=1:n,
    polmatq{i} = polyvalm(polchar(1:n+1-i),Acan);
    P_star((i-1)*p+1:i*p,:) = P*polmatq{i};
  end
  D = Qc*P_star;
  if rank(D) ~= n,;    % Checking the rank of D
    lambdades = [];
    disp('The matrix D = Qc * [P*q1(Acan); ...; P*qn(Acan)] must be invertible');
    disp('Start again choosing different eigenvalues and parameter vectors'); beep;
  end
end
V = inv(q0)*D, L_Bec_Ost = P*inv(V); P
L = Rounding_Matrix(real(L_Bec_Ost));
lambdaCL = lambdades;