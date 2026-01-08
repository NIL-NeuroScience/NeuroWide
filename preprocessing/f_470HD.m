%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% f_470HD %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Hemodynamic correction code for green fluorescence channel
%
% INPUTS:
%   gfp: intensity of green fluorescence channel
%   HbO: estimated change in [HbO]
%   HbR: estimated change in [HbR]
%
% OUTPUTS:
%   gfp: Delta F/F of green fluorescence channel
%   gfp_HD: hemodynamic corrected Delta F/F of green fluorescence channel
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [gfp,gfp_HD] = f_470HD(gfp,HbO,HbR)

if nargin > 1
    tmpWL = [470 515];
    tmpExtinction = f_GetExtinctions(tmpWL);     % in cm
    tmpPathEx = f_pathlengths(tmpWL(1),0.4)/2;       % pathlengths returns in cm
    tmpPathEm = f_pathlengths(tmpWL(2),0.4)/2;
    tmpMuaEx = (tmpExtinction(1,1).*HbO) + (tmpExtinction(1,2).*HbR);
    tmpMuaEm = (tmpExtinction(2,1).*HbO) + (tmpExtinction(2,2).*HbR);
    gfp = gfp./mean(gfp,3);
    gfp_HD = (gfp)./exp(-(tmpMuaEx.*tmpPathEx + tmpMuaEm.*tmpPathEm))-1;
    gfp = gfp-1;
else
    gfp = gfp./mean(gfp,3)-1;
    gfp_HD = [];
end