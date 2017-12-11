function [ qg, Kg, sn, epsn ] = corotbeamplastic ( b, h, x, d, so, epso, do )
%corotbeam returns the internal force vector and the tangent stiffness
%matrix for a corotational beam element with linear formulation
%   EA - axial stiffness
%   EI - bending stiffness
%   x - vector with coordinates of points [x1 y1 x2 y2]
%   d - displacement vector [u1 v1 theta1 u2 v2 theta2 ]
%   do - last converged displacement [u1 v1 theta1 u2 v2 theta2 ]
global wx wy E Et yield H yg xg
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

%local displacement vector of new state
pen=[u_bar theta1_bar theta2_bar]';

%% calculate displacement vector for previous state
d=do;

u21=d(4)-d(1);
v21=d(5)-d(2);
theta1=d(3);
theta2=d(6);
%define initial and current lengths
Ln=sqrt((x21+u21)^2+(y21+v21)^2);
%rigid rotation
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

peo=[u_bar theta1_bar theta2_bar]';

%% Linear strain definition with elastoplastic material
%%define Gauss points
Ket=zeros(3,3);
qet=zeros(3,1);
sn=zeros(7,2);
epsn=zeros(7,2);

for j=1:2
    for i=1:7
        
        %calculate strain variation
        St=[1/Lo -yg(i)*h/2*(6*xg(j)*Ln/Lo^2-4/Lo) -yg(i)*h/2*(6*xg(j)*Ln/Lo^2-2/Lo)];
        
        de=St*(pen-peo);
     
        %calculate stress and plastic strain from incremental strain
        [sn(i,j), epsn(i,j), En ] = pstress1d( so(i,j),epso(i,j),de, E, Et, H, yield );
        %local element stiffness matrix (not valid in elasto-plastic domain)
%         Ke=[En*A/Lo   0        0
%             0     4*En*I/Lo   2*En*I/Lo
%             0     2*En*I/Lo   4*En*I/Lo ];
        %local internal force vector
        
        qe=St'*sn(i,j);
        N=qe(1);
        M1=qe(2);
        M2=qe(3);
        Ke=St'*En*St;
        Ket=Ket+b*Ke*wx(j)*Lo*wy(i)*h/2;
        qet=qet+b*qe*wx(j)*Lo*wy(i)*h/2;
    end
end
%calculate global stiffness matrix
Kg=B'*Ket*B+BN*N+BM*(M1+M2);        
%global internal force vector
qg=B'*qet;

end

