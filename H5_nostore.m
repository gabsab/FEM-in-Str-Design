% Lee's frame 10 elements
% non-linear elastoplastic analysis (corotational formulation)

%% Initialise
clear all;
global wx wy E Et yield H xg yg  

%% weight factors
wx=[0.5 0.5];       xg=[0.21132 0.78868];
wy=[0.12948 0.27971 0.38183 0.41796 0.38183 0.27971 0.12948];    yg=[-0.94911 -0.74153 -0.40585 0 0.40585 0.74153 0.94911];

%% monitors
forc=[];
dispx=[];
dispy=[];

%% tolerance
tol=1E-4; 

%% geometry
L=120;
a=3;
b=2;
I=a*b^3/12;
A=a*b;
%% material properties
E=720;
Et=E/10;
yield=10.44;
H=Et/(1-Et/E);

%% nodes coordinates
coor=zeros(11,8);
elem=zeros(10,2);

Ne=size(elem,1);
Np=size(coor,1);

%     x   y  BC_x  BC_y   BC_xy  Load_x Load_y M_xy
coor=[0   0   1     1       0      0      0     0    
      0   24  0     0       0      0      0     0
      0   48  0     0       0      0      0     0
      0   72  0     0       0      0      0     0
      0    96 0     0       0      0      0     0
      0   120 0     0       0      0      0     0
      24  120 0     0       0      0     -1     0
      48  120 0     0       0      0      0     0
      72  120 0     0       0      0      0     0
      96  120 0     0       0      0      0     0
      120 120 1     1       0      0      0     0];
  
  elem=[1 2
        2 3
        3 4
        4 5
        5 6
        6 7
        7 8 
        8 9
        9 10
        10 11];

afg=[];
P=[];

for i=1:Np
    if coor(i,3)==0 
       afg=[afg;3*(i-1)+1];
       P=[P;coor(i,6)];       
    end
    if coor(i,4)==0 
       afg=[afg;3*(i-1)+2];
       P=[P;coor(i,7)];       
    end    
    if coor(i,5)==0 
       afg=[afg;3*(i-1)+3];
       P=[P;coor(i,8)];       
    end    
end  

nfg=size(afg,1);

%% assembling 
%construct stiffness matrix
Kg=zeros(3*Np,3*Np); % size of matrix determined by number of degrees of freedom per node
qg=zeros(3*Np,1); % internal force vector
d=zeros(3*Np,1); %displacement matrix
ddo=d;
do=ddo;
Pr=d; %full matrix of external forces
residualg=d; %residual forces global matrix
sgo=zeros(7*Ne,2); %stress matrix for gausspoints
epsgo=zeros(7*Ne,2);
ego=zeros(7*Ne,2);
sn=zeros(7,2);
epsn=zeros(7,2);
%initialisation
for i=1:Ne  
  m1=elem(i,1);
  m2=elem(i,2);
  %coordinates matrix
  x=[coor(m1,1:2) coor(m2,1:2)]';
  v=[3*m1-2:3*m1 3*m2-2:3*m2]';
  gp=[7*m1-6:7*m2-7]';
  [q,K, sn, epsn]=corotbeamplastic ( a, b, x, d(v), sgo(gp,1:2) ,epsgo(gp,1:2), ddo(v) );
  Kg(v,v)=Kg(v,v)+K;
  
end
Kr=Kg(afg,afg);  %initial tangent stiffness matrix

%displacement increment
d_inc=-0.1;
%d(afg,1)=d_inc*P;
aux=[zeros(17,1); 1 ; zeros(11,1)];
dd=d_inc*aux;

%% apply displacement control 
%define tangent stiffnes matrix
Ka=[Kr,-P;aux',0];
%lambda
la=0;
max_idx=18;
% %previous displacement
% d1=d(afg,1);
% %current displacement
% d2=d1;
%direction vector
% s=1;
%perform analysis
step3=0;
%initialize counter 
%c=1;
while d(19) < 100
    
    %assemble tangent vector 
    t=[Kr\P;1];
    %predictor
    %identify maximum displacement: node number and direction
    if d(19) < 90 && abs(d(20)) > 48
        max_idx=17;
        d_inc=0.1;
    end 
%     
%     if abs(d(19)) >= 90    
%         max_idx=18; 
%         d_inc=-0.02;
%         step3=1;
    % end
    
    
    
    %change displacement control to the node with the highest displacement,
    aux=[zeros(max_idx-1,1); 1 ; zeros(length(t)-max_idx-1,1) ]; 
    dd=aux*d_inc;
    predictor=dd(max_idx)/t(max_idx)*t;
    %write dlambda and ddisplacement
    dd=predictor(1:end-1);
    dF=predictor(end);
    %update displacement matrix
    d(afg,1)=d(afg,1)+dd;
    
    la=la+dF;
    %create the big matrix
    resid=100;
    
        
  
    while resid>tol
        %assemble        
        Kg=zeros(3*Np,3*Np);
        qg=zeros(3*Np,1);
        sgn=zeros(7*Ne,2);
        epsgn=zeros(7*Ne,2);
        for i=1:Ne
            m1=elem(i,1);
            m2=elem(i,2);
            %coordinates matrix
            x=[coor(m1,1:2) coor(m2,1:2)]';
            v=[3*m1-2:3*m1 3*m2-2:3*m2]';
            gp=[7*m1-6:7*m2-7]';
            [q,K,sn,epsn]=corotbeamplastic( a, b, x, d(v),sgo(gp,1:2), epsgo(gp,1:2),do(v) );
            Kg(v,v)=Kg(v,v)+K;
            qg(v,1)=qg(v,1)+q;
            sgn(gp,1:2)=sgn(gp,1:2)+sn;
            epsgn(gp,1:2)=epsgn(gp,1:2)+epsn;
            %egn(gp,1:2)=egn(gp,1:2)+en;
        end
        qr=qg(afg,1);
        Kr=Kg(afg,afg);  %converged stiffness matrix
        Ka=[Kr,-P;aux',0];
        residualr=qr-la*P;
        resid=norm(residualr);
        if resid>tol
            um = -Ka\[residualr;0];
            dd=um(1:end-1);   
            dF=um(end);
            d(afg,1)=d(afg,1)+dd;
            la=la+dF;            
        end          
    end
    
    %store incremental displacements at every converged step
    %ddo(afg,1)=d(afg,1)-do(afg,1);
    %do(afg,s)=d(afg,1);
    
    %converged step counter 
    %c=c+1;
    do= d;
    sgo= sgn; 
    %stored the displacement matrix at converged step
    
    epsgo= epsgn;
    %ego=egn;  
        
    forc=[forc,la];
    dispx=[dispx,d(19)];
    dispy=[dispy,d(20)];
    
%     defbeam(coor,elem,d,1)
%     hold on
    
    
end

% d=zeros(3*Np);
% d(afg)=dr;

P=1
u=d(19)
v=-d(20)

defbeam(coor,elem,d,0.5)

 