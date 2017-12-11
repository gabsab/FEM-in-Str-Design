function [ de ] = strain_increment( do,dn )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
for j=1:2
    for i=1:7
        
        %calculate strain variation
        St=[1/Lo -yg(i)/2*(6/Lo^2-12*xg(j)/Lo^3) -yg(i)/2*(6*xg(j)/Lo^2-2/Lo)];
        de=St*pe;
        %calculate stress and plastic strain
        [sn(i,j), epsn(i,j), En ] = pstress1d( so(i,j),epso(i,j),de, E, Et, H, yield );
        %local element stiffness matrix
        Ke=[En*A/Lo   0        0
            0     4*En*I/Lo   2*En*I/Lo
            0     2*En*I/Lo   4*En*I/Lo ];
        %local internal force vector
        
        qe=Ke*pe;
        N=qe(1);
        M1=qe(2);
        M2=qe(3);
        
        Ket=Ket+Ke*wx(j)*wy(i)/2;
        qet=qet+qe*wx(j)*wy(i)/2;
    end
end

end

