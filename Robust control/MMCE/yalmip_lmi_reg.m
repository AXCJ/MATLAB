%YALMIP_LMI_REG   Specify LMI regions for pole confinement
%
%    If < 7 input arguments, several partial LMI regions can be specified 
%    successively for state-feedback design with pole confinement constraints.
%    At each step, the matrices Lreg and Mreg of the total LMI region for the 
%    pole confinement are augmented with the corresponding matrices of the 
%    new partial LMI region. At the end of these acquisitions, the resulting
%    total LMI region is added to the input LMI constraint "F" as a set 
%    representing the pole confinement constraint, and returned to the 
%    calling program as output argument "F_polreg" (mode = 1).
%
%    If display = "On" (default is "Off"), the partial LMI region matrices 
%    are displayed as they are created, and the total LMI region matrices are
%    displayed at the end of the program.
%    If < 6 arguments, display is automatically set to "Off".
%
%    If there are 8 arguments, no interactive region creation takes place; 
%    the Lreg and Mreg matrices transmitted by input arguments and defining
%    the chosen LMI region in the calling program are used directly to build
%    the corresponding LMI constraint (mode = 0).
%
% Syntax:
%    [F_polreg,Lreg,Mreg] = yalmip_lmi_reg(F,A,B,Y,W,display,Lreg,Mreg)
%
% Author : Eric Ostertag, 17 April 2010
% Last update: 11 August 2010
%
% N.B. : notations: same as Scherer, Gahinet and Chilali (1997);
%        region: R = {z in C: Lreg + z*Mreg + conj(z)*Mreg' < 0}
%        lmi: [lij*P + mij*A'P + mji*PA]i,j < 0,  P>0
%
% This program contains statements from lmireg.m in Robust Control Toolbox
% (MATLAB)

function [F_polreg,Lreg,Mreg] = yalmip_lmi_reg(F,A,B,Y,W,display,Lreg,Mreg)

if nargin < 6
  display = 'Off';
  mode = 1;
  Lreg = []; Mreg = [];
elseif nargin < 7
  mode = 1;
  Lreg = []; Mreg = [];
else
  mode = 0;
end

while mode,
 disp(sprintf('Select a region among the following:'));
 disp('     h)   Half-plane ');
 disp('     d)   Disk ');
 disp('     c)   Conic sector ');
 disp('     e)   Ellipse ');
 disp('     p)   Parabola ');
 disp('     s)   Strip (horizontal) ');
 disp('     q)   Quit');
 choice = input('Your choice: ','s');
 choice = lower(choice(1:min(1,length(choice))));

 if strcmp(choice,'q'),          %%%%%%%% quit
  mode=0;

 elseif strcmp(choice,'h'),  %%%%%%%% Vertical half-plane
  x0 = [];
  e = input('  Orientation (x < x0 -> l (left), x > x0 -> r (right)):  ','s');
  e = lower(e);
  while ~strcmp(e,'l') & ~strcmp(e,'r'),
     e = input('  Enter l or r:  ','s');
     e = lower(e);    
  end
  if e=='l', e=1; else e=-1; end
  x0 = input('  Specify x0:  ');
  while length(x0)~=1 | ischar(x0) | ~isreal(x0),
     x0 = input('  Specify a real number:  ');
  end
  Lreg = mdiag(Lreg,-2*e*x0); 
  Mreg = mdiag(Mreg,e);

 elseif choice == 'c'          %%%%%%%% Conic sector
  x0=[]; t=[];
  x0 = input('  Abscissa x0 of the tip of the sector:  ');
  while length(x0)~=1 | ischar(x0) | ~isreal(x0),
     x0 = input('  Enter a real number:  ');
  end
  t = input('  Inner half-angle (half-angle < pi/2 ==> sector contains x = -Inf):  ');
  if isempty(t), t=-1; end
  while length(t)~=1 | ischar(t) | t<=0 | t>= pi | imag(t),
  t = input('  Enter a half-angle in [0, pi]:  ');
  if isempty(t), t=-1; end
  end
  s=sin(t); c=cos(t);
  Lreg = mdiag(Lreg,2*[-s*x0 0;0 -s*x0]);
  Mreg = mdiag(Mreg,[s -c;c s]);

 elseif choice == 'd'         %%%%%%%% Disk
  q=[]; r=[];
  q = input('  Abscissa q of the center:  ');
  while length(q)~=1 | ischar(q) | ~isreal(q),
     q = input('  Enter a real number:  ');
  end
  r = input('  Radius r:  ');
  if isempty(r), r=-1; end
  while length(r)~=1 | ischar(r) | r<=0 | ~isreal(r),
    r = input('  Enter a positive real number:  ');
    if isempty(r), r=-1; end
  end
  Lreg = mdiag(Lreg,[-r -q;-q -r]);
  Mreg = mdiag(Mreg,[0 0;1 0]);
  
 elseif choice == 'e'  %%%%%%%% Ellipse
  q=[]; a=[]; b=[];
  disp(' ');
  disp('  Ellipse of equation  ((x-q)/a)^2 + (y/b)^2 < 1 ');
  q = input('  Abscissa q of the center:  ');
  while length(q)~=1 | ischar(q) | ~isreal(q),
     q = input('  Enter a real number:  ');
  end
  a = input('  Half-length a of the horizontal axis:  ');
  if isempty(a), a=-1; end
  while length(a)~=1 | ischar(a) | a <=0 | ~isreal(a),
    a = input('  Enter a positive real number:  ');
    if isempty(a), a=-1; end
  end
  b = input('  Half-length b of the vertical axis:  ');
  if isempty(b), b=-1; end
  while length(b)~=1 | ischar(b) | b <=0 | ~isreal(b),
    b = input('  Enter a positive real number:  ');
    if isempty(b), b=-1; end
  end
  a=1/abs(a);b=1/abs(b);
  Lreg = mdiag(Lreg,[-1 -q*a;-q*a -1]);
  Mreg = mdiag(Mreg,[0 (a-b)/2;(a+b)/2 0]);

 elseif choice == 'p'  %%%%%%%%% Parabola
  x0=[]; p=[];
  disp(' ');
  disp('  Parabola of equation  y^2 + p*(x-x0) < 0');
  x0 = input('  Enter x0:  ');
  while length(x0)~=1 | ischar(x0) | ~isreal(x0),
     x0 = input('  Enter a real number:  ');
  end
  p = input('  Enter the parameter p:  ');
  while length(p)~=1 | ischar(p) | ~isreal(p),
     p = input('  Enter a real number:  ');
  end
  q=-p*x0;
  Lreg = mdiag(Lreg,2*[q-1 q+1;q+1 q-1]);
  Mreg = mdiag(Mreg,[p p-2;p+2 p]);

 elseif choice == 's'  %%%%%%%% Strip (horizontal)
  r=[];
  disp(' ');
  disp('  Horizontal strip of equation  -r < y < r');
  r = input('  Enter r:  ');
  if isempty(r), r=-1; end
  while length(r)~=1 | ischar(r) | r<=0 | ~isreal(r),
    r = input('  Enter a positive real number:  ');
    if isempty(r), r=-1; end
  end
  Lreg = mdiag(Lreg,[-r 0;0 -r]);
  Mreg = mdiag(Mreg,[0 -1;1 0]); 
  
 end
 if mode == 1 && strcmp(display,'On')
   disp('LMI region obtained at this step:')
   Lreg,Mreg
 end
end
F = F + set(kron(Lreg,Y) + kron(Mreg,(A*Y-B*W)') + kron(Mreg',(A*Y-B*W)) < 0);
if nargin < 8 && strcmp(display,'On')
    disp('Total LMI region:')
    Lreg, Mreg
end
F_polreg = F;
