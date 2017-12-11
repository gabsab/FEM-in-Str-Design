function [ qg, Kg, Ks ] = corotbeam ( EA, EI, x, d )
%corotbeam returns the internal force vector and the tangent stiffness
%matrix for a corotational beam element with linear formulation
%   EA - axial stiffness
%   EI - bending stiffness 
%   x - vector with coordinates of points [x1 y1 x2 y2]
%   d - displacement vector [u1 v1 theta1 u2 v2 theta2 ]

x21=x(3)-x(1);
y21=x(4)-x(2);     
u21=d(4)-d(1);
v21=d(5)-d(2);
theta1=d(3);
theta2=d(6);
%define initial and current lengths
Lo=sqrt(x21^2+y21^2);
Ln=sqrt((x21+u21)^2+(y21+v21)^2);
%rigid rotation
cos_betha0=x21/Lo;
sin_betha0=y21/Lo;
cos_betha=(x21+u21)/Ln;
c=cos_betha;
sin_betha=(y21+v21)/Ln;
s=sin_betha;
sin_alpha= sin_betha*cos_betha0-cos_betha*sin_betha0;
cos_alpha= cos_betha*cos_betha0+sin_betha*sin_betha0;
if cos_alpha >=0
    alpha=asin(sin_alpha);
elseif cos_alpha < 0 && sin_alpha >= 0
    alpha=acos(cos_alpha);
elseif  cos_alpha < 0 && sin_alpha < 0
    alpha=-acos(cos_alpha);
end
%alpha=asin(sin_betha)-asin(sin_betha0);
%calculate components of the local displacement vector
u_bar=Ln-Lo;
theta1_bar=theta1-alpha;
theta2_bar=theta2-alpha;
%internal force vector
%qe=[N M1 M2]';

% transformation matrix
B=[  -c        -s        0       c        s        0
   -s/Ln      c/Ln       1      s/Ln    -c/Ln      0
   -s/Ln      c/Ln       0      s/Ln    -c/Ln      1 ];
%notations
r=[-c  -s  0  c  s  0]';
z=[ s  -c  0  -s c  0]';
%b1=r;
%b2=[0  0  1  0  0  0]'*(-z/Ln);
%b3=[0  0  0  0  0  1]'*(-z/Ln);
%B=[b1';b2';b3'];
BN=z*z'/Ln;
BM=1/Ln^2*(r*z'+z*r');



%% Linear strain definition
%local displacement vector
pe=[u_bar theta1_bar theta2_bar]';
%local element stiffness matrix
Ke=[EA/Lo   0        0
    0     4*EI/Lo   2*EI/Lo 
    0     2*EI/Lo   4*EI/Lo ];
%local internal force vector
qe=Ke*pe;
N=qe(1);
M1=qe(2);
M2=qe(3);
%calculate global stiffness matrix
Kg=B'*Ke*B+BN*N+BM*(M1+M2);
%global internal force vector
qg=B'*qe;

end

