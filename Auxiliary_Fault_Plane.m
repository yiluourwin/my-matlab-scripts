% -------------------------------------------
% Calculate Auxiliary Fault Plane
% 
% Author: Yi LUO
% Built: 18/11/21 11:00
% -------------------------------------------
% First 3 input arguments are required. Default isdegree=1. 
% -------------------------------------------

function [phi2,sigma2,lambda2]=Auxiliary_Fault_Plane(phi,sigma,lambda,isdegree)

if nargin<4
    isdegree=1;
end

if isdegree
    phi=phi./180.*pi;
    sigma=sigma./180.*pi;
    lambda=lambda./180.*pi;
end

sigma2=acos(sin(lambda).*sin(sigma));
% if sigma2>pi/2
%     sigma2=pi-sigma2;
% end

sinphi1_2=real(cos(lambda)./sin(sigma2));
cosphi1_2=real(-1./tan(sigma)./tan(sigma2));

if sinphi1_2 > 0 && cosphi1_2 > 0
    phi1_2=real(asin(complex(sinphi1_2)));
elseif sinphi1_2 > 0 && cosphi1_2 < 0
    phi1_2=real(pi-asin(complex(sinphi1_2)));
elseif sinphi1_2 < 0 && cosphi1_2 < 0
    phi1_2=real(pi-asin(complex(sinphi1_2)));
else
    phi1_2=real(2.*pi+asin(complex(sinphi1_2)));
end

phi2=phi-phi1_2;

sinlam2=cos(sigma)/sin(sigma2);
coslam2=sin(sigma)*sin(phi2-phi);

if sinlam2 > 0 && coslam2 > 0
    lambda2=real(asin(complex(sinlam2)));
elseif sinlam2 > 0 && coslam2 < 0
    lambda2=real(pi-asin(complex(sinlam2)));
elseif sinlam2 < 0 && coslam2 < 0
    lambda2=real(pi-asin(complex(sinlam2)));
else
    lambda2=real(2.*pi+asin(complex(sinlam2)));
end


if sigma2 > pi/2
    phi2=pi+phi2;
    sigma2=pi-sigma2;
    lambda2=2*pi-lambda2;
end

if phi2 < 0
    phi2=phi2+2.*pi;
elseif phi2 > 2*pi
    phi2=phi2-2.*pi;
end


if isdegree
    phi2=phi2.*180./pi;
    sigma2=sigma2.*180./pi;
    lambda2=lambda2.*180./pi;
end

phi2=real(phi2);
sigma2=real(sigma2);
lambda2=real(lambda2);
