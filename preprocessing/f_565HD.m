%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% f_565HD %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Hemodynamic correction code for red fluorescence channel. Algorithm was
% taken from Shahsavarani et al., 2023, Cell Reports.
% 
% INPUTS:
%   rfp: intensity of red fluorescence channel
%   HbO: intensity of 525 nm channel
%   HbR: intensity of 625 nm channel
% 
% OUTPUTS:
%   rfp: Delta F/F of red fluorescence channel
%   rfp_HD: hemodynamic corrected Delta F/F of red fluorescence channel
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [rfp,rfp_HD] = f_565HD(rfp,HbO,HbR)

rfp = rfp./mean(rfp,3);

if nargin > 1
    tmpHbRed = HbR./mean(HbR,3);
    tmpHbGreen = HbO./mean(HbO,3);
    rfp_HD = rfp./(tmpHbRed.^0.8.*tmpHbGreen.^0.4);
    rfp_HD = rfp_HD-1;
else
    rfp_HD = [];
end

rfp = rfp-1;
