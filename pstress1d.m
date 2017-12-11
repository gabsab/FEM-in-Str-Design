function [sn, epsn, En ] = pstress1d( so,epso,de, E, Et, H, yield )
%calculates, under one dimension assumption, the new stress sn, the 
%new equivalent plastic strain epsn and the tangent operator En as function
%of the old stress so, the old equivalent plastic strain epso, the strain increment de
%and the material parameters.
%   sn - new stress
%   epsn - new plastic strain
%   En - tangent operator
%   s0- previous stress
%   espo- previoues strain
%   de - strain increment
%   E - elastic modulus
%   Et - plastic modulus
%   H - hardening coefficient
%   yield - current yield point

%the yield stress at state 1
sy=yield+H*epso;
% elastic predictor
sb=so+E*de;
ftol=yield*10e-6;
% fyld=abs(sb)-abs(sy);
if abs(sb) - abs(sy) <= ftol
    sn=sb;
    epsn=epso;
    En=E;
else
    
   epsn = epso +(abs(sb)-sy)/(E+H);
   sn = sb-E/(E+H)*(sb-sign(sb)*sy);%+E*de;
   %sn = sb-sign(sb)*E*(abs(sb)-sy)/(E+H);
   En=Et;
%     sn=sb;
%     epsn=epso;
%     En=E;
   
   
end

