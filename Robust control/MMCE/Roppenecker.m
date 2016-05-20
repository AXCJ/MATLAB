function [L,lambdaCL] = Roppenecker(A,B,Lambda,Vp,lambdades,Ts)
%
% Design of a multivariable controller or observer by the Roppenecker's formula,
% or "Complete modal control" design.
%
%  [L,lambdaCL] = Roppenecker(A,B,Lambda,Vp,[],Ts) begins acquiring the continuous-time
%  eigenvalues desired for the closed-loop system, and returns them as an output argument
%  in the vector lambdaCL. No eigenvalue can have an order of multiplicity greater than
%  the number p of plant inputs.
%     If Ts = 0, the entered eigenvalues are returned without modification;
%     if Ts > 0, the entered eigenvalues are converted to their discrete
%     equivalent (sampling period Ts) by exp(lambda*Ts), before being returned.
%  The input matrix argument Lambda is the diagonal matrix of the eigenvalues of A, thus of the
%  open-loop system, for comparison purposes with the eigenvalues prescribed for the 
%  closed loop. The subroutine prompts then the user to enter the invariant parameter vectors,
%  correlated or not with the choice of the closed-loop eigenvectors, and computes with
%  the Roppenecker's formula the state-feedback matrix L which assigns the entered 
%  eigenvalues to the matrix A-B*L.
%
%  [L,lambdaCL] = Roppenecker(A,B,Lambda,Vp,lambdades,Ts) does the same thing, except that
%  the eigenvalues desired for the closed loop are given as input argument by means of the
%  vector lambdades. They are thus not entered anymore by the user.
%
%  Vp is the matrix whose columns are the plant eigenvectors, for the case one of the
%  desired eigenvalues is identical to one of the plant eigenvalues (non shifted eigenvalue).
%
%  N.B.: this subroutine is also used, by duality, for the computation of observers 
%  (full or reduced-order).
%
%  Called subroutines: Get_Eigenvalue, Choice_Eigenvec_Paramvec, Rounding_Matrix.
%
%  Author: E. Ostertag, 10 February 2010 
%  Last update: 11 August 2010
%

n = size(A,1); p = size(B,2);
%
% 1. Choice of the eigenvalues desired for the closed-loop system, except if they belong to the
%    input arguments
%
ni = nargin;
entry = 0;
if ni<5,
  entry = 1;
elseif isempty(lambdades),
  entry = 1;
end
if entry == 1,
  disp('N.B.: Order of multiplicity of the chosen eigenvakues <= number of inputs')
  disp(['Enter n=' num2str(n) ' eigenvalues for the closed-loop system ===>']);
  lambdades = Get_Eigenvalue(n,0,p);
end
if ni == 6 && Ts > 0,
  lambdades = exp(lambdades*Ts);
end
%
% 2. Choice of the parameter vectors and application of the formula
%
V = zeros(n);
while det(V) == 0,
  [V, P, choice, kident] = Choice_Eigenvec_Paramvec(A,B,n,p,Lambda,lambdades);
  if choice == 1,
% In case a desired closed-loop eigenvalues is identical to one open-loop
% eigenvalue, parameter vector = open-loop eigenvector
    for i = 1:n,
      if kident(i) == 0;
        V(:,i) = inv(A-lambdades(i)*eye(n))*B*P(:,i);
      else
        P(:,i) = 0;  V(:,i) = Vp(:,kident(i));
      end
    end
  end
%
% 3. Checking the linear independence of the new eigenvectors, vRi:
%
  if det(V) == 0;
    disp('The vectors vLi = inv(A-lambdades(i)*I)*B*pi must be linearly independent');
    disp('Start again choosing different pi vectors'); beep;
  else
    V, P
  end
end
L_Roppenecker = P*inv(V);
L = Rounding_Matrix(real(L_Roppenecker));
lambdaCL = lambdades;