% --------------------------------
% Qkappa versus Qp
%
% Author: Yi LUO
% Built: 2019/05/28 14:30
% --------------------------------

function [Qkappa, Qmu] = Qp2Qkappa( Qp, Qs, cp, cs)
Qmu=Qs;

Qkappa=(1./(1-(cs./cp).^2).*(Qp.^-1-(cs./cp).^2.*Qs.^-1)).^-1;

end
